use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Write, BufWriter};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};
use std::thread; // Added for available_parallelism

use clap::Parser;
use rayon::prelude::*;

// --- Random Number Generator (LCG) ---
struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        LcgRng { state: if seed == 0 { 1 } else { seed } }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.state
    }

    fn gen_range_inclusive(&mut self, min: usize, max: usize) -> usize {
        if min > max { return min; }
        let range_size = max - min + 1;
        min + (self.next_u64() % range_size as u64) as usize
    }
    
    fn gen_range_exclusive(&mut self, upper_bound: usize) -> usize {
        if upper_bound == 0 { return 0; }
        (self.next_u64() % upper_bound as u64) as usize
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = self.gen_range_exclusive(i + 1);
            slice.swap(i, j);
        }
    }
}

// --- FASTA I/O ---
#[derive(Debug, Clone)]
pub struct FastaRecord {
    header: String,
    sequence: Vec<u8>,
}

mod fasta_io {
    use super::FastaRecord;
    use std::io::{self, BufRead, Write};

    pub struct FastaParser<R: BufRead> {
        reader: R,
        buffer: String,
        current_header: Option<String>,
        current_sequence_parts: Vec<Vec<u8>>,
    }

    impl<R: BufRead> FastaParser<R> {
        pub fn new(reader: R) -> Self {
            FastaParser {
                reader,
                buffer: String::new(),
                current_header: None,
                current_sequence_parts: Vec::new(),
            }
        }

        fn process_buffered_record(&mut self) -> Option<io::Result<FastaRecord>> {
            if let Some(header) = self.current_header.take() {
                if !self.current_sequence_parts.is_empty() {
                    let mut sequence_parts_len = 0;
                    for part in &self.current_sequence_parts {
                        sequence_parts_len += part.len();
                    }
                    let mut sequence = Vec::with_capacity(sequence_parts_len);
                    for part in self.current_sequence_parts.drain(..) {
                        sequence.extend(part);
                    }
                    
                    for base in &mut sequence {
                        *base = base.to_ascii_uppercase();
                        if *base == b'U' {
                            *base = b'T';
                        }
                    }
                    return Some(Ok(FastaRecord { header, sequence }));
                }
            }
            None
        }
    }
    
    impl<R: BufRead> Iterator for FastaParser<R> {
        type Item = io::Result<FastaRecord>;

        fn next(&mut self) -> Option<Self::Item> {
            loop {
                self.buffer.clear();
                match self.reader.read_line(&mut self.buffer) {
                    Ok(0) => { 
                        return self.process_buffered_record();
                    }
                    Ok(_) => { 
                        let is_header_line = {
                            let trimmed_line_ref = self.buffer.trim_end();
                            trimmed_line_ref.starts_with('>')
                        };

                        if is_header_line {
                            let record_to_return = self.process_buffered_record(); 
                            let new_header_str = self.buffer.trim_end()[1..].to_string();
                            self.current_header = Some(new_header_str);
                            if record_to_return.is_some() {
                                return record_to_return; 
                            }
                        } else { 
                            if self.current_header.is_some() { 
                                let sequence_line_trimmed = self.buffer.trim_end();
                                if !sequence_line_trimmed.is_empty() {
                                    self.current_sequence_parts.push(sequence_line_trimmed.as_bytes().to_vec());
                                }
                            }
                        }
                    }
                    Err(e) => return Some(Err(e)),
                }
            }
        }
    }

    pub fn write_fasta_record<W: Write>(
        writer: &mut W,
        record: &FastaRecord,
        line_width: usize,
    ) -> io::Result<()> {
        writeln!(writer, ">{}", record.header)?;
        for chunk in record.sequence.chunks(line_width) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
        Ok(())
    }
}

// --- Genetic Code ---
#[derive(Clone, Debug)]
pub struct GeneticCode {
    codon_to_aa: HashMap<[u8; 3], u8>,
    aa_to_codons: HashMap<u8, Vec<[u8; 3]>>,
    start_codons: Vec<[u8; 3]>,
    stop_codons: Vec<[u8; 3]>,
}

mod genetic_codes {
    use super::{GeneticCode, HashMap};
    use std::io::{self, BufRead, ErrorKind};

    impl GeneticCode {
        pub fn get_aa(&self, codon: &[u8; 3]) -> Option<u8> {
            self.codon_to_aa.get(codon).cloned()
        }

        pub fn get_synonymous_codons(&self, aa: u8) -> Option<&Vec<[u8; 3]>> {
            self.aa_to_codons.get(&aa)
        }

        #[allow(dead_code)]
        pub fn is_start_codon(&self, codon: &[u8; 3]) -> bool {
            self.start_codons.contains(codon)
        }
        
        #[allow(dead_code)]
        pub fn is_stop_codon(&self, codon: &[u8; 3]) -> bool {
            self.stop_codons.contains(codon)
        }
    }

    pub fn canonical() -> GeneticCode {
        let str_definitions: Vec<(&str, u8, bool)> = vec![
            ("TTT", b'F', false), ("TTC", b'F', false), ("TTA", b'L', false), ("TTG", b'L', false),
            ("TCT", b'S', false), ("TCC", b'S', false), ("TCA", b'S', false), ("TCG", b'S', false),
            ("TAT", b'Y', false), ("TAC", b'Y', false), ("TAA", b'*', false), ("TAG", b'*', false),
            ("TGT", b'C', false), ("TGC", b'C', false), ("TGA", b'*', false), ("TGG", b'W', false),
            ("CTT", b'L', false), ("CTC", b'L', false), ("CTA", b'L', false), ("CTG", b'L', true),
            ("CCT", b'P', false), ("CCC", b'P', false), ("CCA", b'P', false), ("CCG", b'P', false),
            ("CAT", b'H', false), ("CAC", b'H', false), ("CAA", b'Q', false), ("CAG", b'Q', false),
            ("CGT", b'R', false), ("CGC", b'R', false), ("CGA", b'R', false), ("CGG", b'R', false),
            ("ATT", b'I', false), ("ATC", b'I', false), ("ATA", b'I', false), ("ATG", b'M', true),
            ("ACT", b'T', false), ("ACC", b'T', false), ("ACA", b'T', false), ("ACG", b'T', false),
            ("AAT", b'N', false), ("AAC", b'N', false), ("AAA", b'K', false), ("AAG", b'K', false),
            ("AGT", b'S', false), ("AGC", b'S', false), ("AGA", b'R', false), ("AGG", b'R', false),
            ("GTT", b'V', false), ("GTC", b'V', false), ("GTA", b'V', false), ("GTG", b'V', true),
            ("GCT", b'A', false), ("GCC", b'A', false), ("GCA", b'A', false), ("GCG", b'A', false),
            ("GAT", b'D', false), ("GAC", b'D', false), ("GAA", b'E', false), ("GAG", b'E', false),
            ("GGT", b'G', false), ("GGC", b'G', false), ("GGA", b'G', false), ("GGG", b'G', false),
        ];
        let mut definitions = Vec::with_capacity(str_definitions.len());
        for (codon_str, aa, is_start) in str_definitions {
            let bytes = codon_str.as_bytes();
            definitions.push(([bytes[0], bytes[1], bytes[2]], aa, is_start));
        }
        build_from_definitions(definitions)
    }

    pub fn from_file<R: BufRead>(reader: R) -> io::Result<GeneticCode> {
        let mut definitions = Vec::new();
        for line_res in reader.lines() {
            let line = line_res?;
            let parts: Vec<&str> = line.trim().split_whitespace().collect();
            if !(2..=3).contains(&parts.len()) {
                return Err(io::Error::new(ErrorKind::InvalidData, format!("Invalid genetic code line format: {}", line)));
            }
            
            let codon_str = parts[0].to_uppercase();
            if codon_str.len() != 3 {
                 return Err(io::Error::new(ErrorKind::InvalidData, format!("Invalid codon length: {}. Expected 3, got {}", codon_str, codon_str.len())));
            }
            let aa_char = parts[1].chars().next().ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Missing amino acid character"))?;
            
            let mut codon_arr = [0u8; 3];
            for (idx, char_base) in codon_str.chars().enumerate() {
                let base_u8 = char_base.to_ascii_uppercase() as u8;
                codon_arr[idx] = if base_u8 == b'U' { b'T' } else { base_u8 };
                if !matches!(codon_arr[idx], b'A' | b'T' | b'C' | b'G') { 
                     return Err(io::Error::new(ErrorKind::InvalidData, format!("Invalid base '{}' in codon {}", char_base, codon_str)));
                }
            }

            let is_start = parts.get(2).map_or(false, |&s| s.to_uppercase() == "S");
            definitions.push((codon_arr, aa_char as u8, is_start));
        }
        Ok(build_from_definitions(definitions))
    }

    fn build_from_definitions(definitions: Vec<([u8; 3], u8, bool)>) -> GeneticCode {
        let mut codon_to_aa = HashMap::new();
        let mut aa_to_codons: HashMap<u8, Vec<[u8; 3]>> = HashMap::new();
        let mut start_codons = Vec::new();
        let mut stop_codons = Vec::new();

        for (codon, aa, is_start) in definitions {
            codon_to_aa.insert(codon, aa);
            aa_to_codons.entry(aa).or_default().push(codon);
            if is_start {
                start_codons.push(codon);
            }
            if aa == b'*' { 
                stop_codons.push(codon);
            }
        }
        GeneticCode { codon_to_aa, aa_to_codons, start_codons, stop_codons }
    }
}

// --- Augmentation Strategies ---
mod augmenters {
    use super::{FastaRecord, GeneticCode, LcgRng};

    fn new_header(original_header: &str, aug_idx: usize, strategy_tag: &str) -> String {
        format!("{}_aug{}_{}", original_header, aug_idx, strategy_tag)
    }

    pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        let mut rc_seq = Vec::with_capacity(seq.len());
        for &base in seq.iter().rev() {
            rc_seq.push(match base {
                b'A' => b'T', b'T' => b'A',
                b'C' => b'G', b'G' => b'C',
                b'N' => b'N', 
                other => other, 
            });
        }
        rc_seq
    }

    pub fn sliding_windows(
        original_record: &FastaRecord, count: usize, window_size: usize, step_size: usize,
        base_aug_idx: &mut usize, rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        if window_size == 0 || step_size == 0 || window_size > original_record.sequence.len() {
            return augmented_records;
        }

        let mut potential_windows = Vec::new();
        for i in (0..).step_by(step_size) {
            if i + window_size > original_record.sequence.len() { break; }
            potential_windows.push(original_record.sequence[i..i + window_size].to_vec());
        }

        if potential_windows.is_empty() { return augmented_records; }
        rng.shuffle(&mut potential_windows);

        for seq_slice in potential_windows.into_iter().take(count) {
            augmented_records.push(FastaRecord {
                header: new_header(&original_record.header, *base_aug_idx, &format!("slide_w{}_s{}", window_size, step_size)),
                sequence: seq_slice,
            });
            *base_aug_idx += 1;
        }
        augmented_records
    }

    pub fn chunks(
        original_record: &FastaRecord, count: usize, num_chunks: usize,
        base_aug_idx: &mut usize, _rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        if num_chunks == 0 || original_record.sequence.is_empty() { return augmented_records; }

        for _ in 0..count {
            let seq_len = original_record.sequence.len();
            let chunk_len = seq_len / num_chunks;
            let remainder = seq_len % num_chunks;
            let mut current_pos = 0;
            for i in 0..num_chunks {
                let current_chunk_len = chunk_len + if i < remainder { 1 } else { 0 };
                if current_chunk_len == 0 { continue; }
                let chunk_seq = original_record.sequence[current_pos..current_pos + current_chunk_len].to_vec();
                augmented_records.push(FastaRecord {
                    header: new_header(&original_record.header, *base_aug_idx, &format!("chunk{}_of{}", i + 1, num_chunks)),
                    sequence: chunk_seq,
                });
                *base_aug_idx += 1;
                current_pos += current_chunk_len;
            }
        }
        augmented_records
    }
    
    pub fn centered_window(
        original_record: &FastaRecord, count: usize, window_size: usize,
        base_aug_idx: &mut usize, _rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        if window_size == 0 || window_size > original_record.sequence.len() { return augmented_records; }
        for _ in 0..count {
            let seq_len = original_record.sequence.len();
            let start = (seq_len.saturating_sub(window_size)) / 2; 
            let end = (start + window_size).min(seq_len);
            let sequence = original_record.sequence[start..end].to_vec();
            if sequence.is_empty() && window_size > 0 { continue; } 
            augmented_records.push(FastaRecord {
                header: new_header(&original_record.header, *base_aug_idx, &format!("center_w{}", window_size)),
                sequence,
            });
            *base_aug_idx += 1;
        }
        augmented_records
    }

    pub fn beginning_truncations(
        original_record: &FastaRecord, count: usize, min_keep_len: usize, max_keep_len: usize,
        base_aug_idx: &mut usize, rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        let seq_len = original_record.sequence.len();
        if seq_len == 0 { return augmented_records; }
        for _ in 0..count {
            let keep_len = rng.gen_range_inclusive(min_keep_len.min(seq_len), max_keep_len.min(seq_len));
            if keep_len == 0 { continue; }
            let sequence = original_record.sequence[..keep_len].to_vec();
            augmented_records.push(FastaRecord {
                header: new_header(&original_record.header, *base_aug_idx, &format!("begtrunc_k{}", keep_len)),
                sequence,
            });
            *base_aug_idx += 1;
        }
        augmented_records
    }

    pub fn end_truncations(
        original_record: &FastaRecord, count: usize, min_keep_len: usize, max_keep_len: usize,
        base_aug_idx: &mut usize, rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        let seq_len = original_record.sequence.len();
        if seq_len == 0 { return augmented_records; }
        for _ in 0..count {
            let keep_len = rng.gen_range_inclusive(min_keep_len.min(seq_len), max_keep_len.min(seq_len));
            if keep_len == 0 { continue; }
            let sequence = original_record.sequence[seq_len - keep_len..].to_vec();
            augmented_records.push(FastaRecord {
                header: new_header(&original_record.header, *base_aug_idx, &format!("endtrunc_k{}", keep_len)),
                sequence,
            });
            *base_aug_idx += 1;
        }
        augmented_records
    }

    pub fn random_truncations(
        original_record: &FastaRecord, count: usize, min_output_len: usize, max_output_len: usize,
        base_aug_idx: &mut usize, rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        let seq_len = original_record.sequence.len();
        if seq_len == 0 { return augmented_records; }
        for _ in 0..count {
            let out_len = rng.gen_range_inclusive(min_output_len.min(seq_len), max_output_len.min(seq_len));
            if out_len == 0 { continue; }
            if out_len > seq_len { continue; }
            let start_max = seq_len - out_len;
            let start = rng.gen_range_inclusive(0, start_max);
            let sequence = original_record.sequence[start..start + out_len].to_vec();
            augmented_records.push(FastaRecord {
                header: new_header(&original_record.header, *base_aug_idx, &format!("randtrunc_o{}", out_len)),
                sequence,
            });
            *base_aug_idx += 1;
        }
        augmented_records
    }

    pub fn synonymous_variants(
        original_record: &FastaRecord, count: usize, genetic_code: &GeneticCode,
        base_aug_idx: &mut usize, rng: &mut LcgRng,
    ) -> Vec<FastaRecord> {
        let mut augmented_records = Vec::new();
        let seq = &original_record.sequence;
        if seq.len() < 6 || seq.len() % 3 != 0 { return augmented_records; } 

        for _i_variant in 0..count {
            let mut new_sequence = seq.clone();
            let mut changed = false;
            
            for c_idx in (3..(seq.len() - 3)).step_by(3) { 
                if c_idx + 3 > new_sequence.len() { break; }
                let codon_slice = &new_sequence[c_idx..c_idx+3];
                let codon: [u8; 3] = [codon_slice[0], codon_slice[1], codon_slice[2]];

                if let Some(aa) = genetic_code.get_aa(&codon) {
                    if aa == b'*' { continue; } 
                    if let Some(syn_codons) = genetic_code.get_synonymous_codons(aa) {
                        if syn_codons.len() > 1 {
                            let mut available_syn_codons: Vec<&[u8;3]> = syn_codons.iter().filter(|sc| **sc != codon).collect();
                            if !available_syn_codons.is_empty() {
                                rng.shuffle(&mut available_syn_codons);
                                let new_codon = available_syn_codons[0];
                                new_sequence[c_idx..c_idx+3].copy_from_slice(new_codon);
                                changed = true;
                            }
                        }
                    }
                }
            }
            if changed {
                 augmented_records.push(FastaRecord {
                    header: new_header(&original_record.header, *base_aug_idx, "synonymous"),
                    sequence: new_sequence,
                });
                *base_aug_idx += 1;
            }
        }
        augmented_records
    }
}

// --- CLI Parsing with Clap ---
#[derive(Debug, Clone, Copy)]
enum StrategyType {
    SlidingWindow, Chunks, CenteredWindow, BeginningTruncation, EndTruncation, RandomTruncation, Synonymous,
}

#[derive(Debug, Clone)]
struct StrategyCliConfig {
    stype: StrategyType,
    count: usize,
    params: Vec<usize>,
}

fn parse_strategy_clap(s: &str) -> Result<StrategyCliConfig, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() < 2 {
        return Err(format!("Invalid strategy format: {}. Expected type:count[:param1:...]", s));
    }

    let strategy_name = parts[0];
    let count: usize = parts[1].parse().map_err(|e| format!("Invalid count in strategy '{}': {}", parts[1], e))?;
    let params_str = parts.get(2..).unwrap_or(&[]);
    let mut params: Vec<usize> = Vec::with_capacity(params_str.len());
    for p_str in params_str {
        params.push(p_str.parse().map_err(|e| format!("Invalid parameter in strategy '{}': {}", p_str, e))?);
    }

    let stype = match strategy_name {
        "sliding" => { if params.len() != 2 { return Err("Sliding strategy needs 2 params: window_size:step_size".to_string()); } StrategyType::SlidingWindow },
        "chunks" => { if params.len() != 1 { return Err("Chunks strategy needs 1 param: num_chunks".to_string()); } StrategyType::Chunks },
        "centered" => { if params.len() != 1 { return Err("Centered strategy needs 1 param: window_size".to_string()); } StrategyType::CenteredWindow },
        "begin_trunc" => { if params.len() != 2 { return Err("Beginning truncation needs 2 params: min_keep_len:max_keep_len".to_string()); } StrategyType::BeginningTruncation },
        "end_trunc" => { if params.len() != 2 { return Err("End truncation needs 2 params: min_keep_len:max_keep_len".to_string()); } StrategyType::EndTruncation },
        "random_trunc" => { if params.len() != 2 { return Err("Random truncation needs 2 params: min_output_len:max_output_len".to_string()); } StrategyType::RandomTruncation },
        "synonymous" => { if !params.is_empty() { return Err("Synonymous strategy takes no parameters beyond count".to_string()); } StrategyType::Synonymous },
        _ => return Err(format!("Unknown strategy type: {}", strategy_name)),
    };
    Ok(StrategyCliConfig { stype, count, params })
}

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
struct CliConfig {
    #[arg(short, long, value_name = "FILE")]
    input_path: PathBuf,

    #[arg(short, long, value_name = "FILE")]
    output_path: PathBuf,

    #[arg(short = 'S', long = "strategy", value_parser = parse_strategy_clap, action = clap::ArgAction::Append)]
    strategies: Vec<StrategyCliConfig>,

    #[arg(long)]
    revcomp: bool,

    #[arg(long, value_name = "FILE")]
    genetic_code_file: Option<PathBuf>,

    #[arg(short, long, value_name = "NUM")] // Removed default_value_t
    threads: Option<usize>, // Changed to Option<usize>
}

// --- Main Application Logic ---
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli_config = CliConfig::parse();

    // Determine the number of threads to use and initialize Rayon's global pool.
    // This must be done *before* any other Rayon operations that might auto-initialize the pool.
    let num_threads = cli_config.threads.unwrap_or_else(|| {
        thread::available_parallelism().map_or(1, |n| n.get())
    });

    if num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .map_err(|e| format!("Failed to build Rayon thread pool: {}. This can happen if the pool is already initialized by another part of the program or a dependency before this configuration attempt.", e))?;
    } else {
        // If num_threads is 0, Rayon might default to a sensible number or this could be an error condition
        // depending on desired behavior. For now, we let Rayon decide if build_global isn't called.
        // Or, one might choose to error out if 0 threads are specified.
        // For simplicity, if 0 is somehow passed, we skip explicit configuration. Rayon will then auto-init.
        eprintln!("Warning: Number of threads specified or defaulted to 0. Rayon will use its default configuration.");
    }
    
    let genetic_code_arc = Arc::new(if let Some(ref gcf_path) = cli_config.genetic_code_file {
        let file = File::open(gcf_path).map_err(|e| format!("Failed to open genetic code file {:?}: {}", gcf_path, e))?;
        genetic_codes::from_file(BufReader::new(file)).map_err(|e| format!("Failed to parse genetic code file {:?}: {}", gcf_path, e))?
    } else {
        genetic_codes::canonical()
    });

    let input_file = File::open(&cli_config.input_path)
        .map_err(|e| format!("Failed to open input file {:?}: {}", cli_config.input_path, e))?;
    
    let all_input_records: Vec<FastaRecord> = fasta_io::FastaParser::new(BufReader::new(input_file))
        .collect::<Result<Vec<_>,_>>() 
        .map_err(|e: std::io::Error| format!("Error parsing FASTA input: {}", e))?; 

    let cli_config_arc = Arc::new(cli_config.clone()); // cli_config is still needed for output_path etc.

    let augmented_batches: Vec<Vec<FastaRecord>> = all_input_records
        .into_par_iter()
        .map_init(
            || { 
                let time_seed = SystemTime::now().duration_since(UNIX_EPOCH)
                    .map_or(0, |d| d.as_nanos() as u64);
                let thread_idx_seed = rayon::current_thread_index().unwrap_or(0) as u64;
                LcgRng::new(time_seed.wrapping_add(thread_idx_seed).wrapping_add(1)) 
            },
            |rng, record| { 
                let mut base_aug_idx_for_record: usize = 1;
                let mut all_augmented_for_this_record = Vec::new();

                for strategy_conf in &cli_config_arc.strategies {
                    let mut current_augs = Vec::new();
                    match strategy_conf.stype {
                        StrategyType::SlidingWindow => current_augs.extend(augmenters::sliding_windows(&record, strategy_conf.count, strategy_conf.params[0], strategy_conf.params[1], &mut base_aug_idx_for_record, rng)),
                        StrategyType::Chunks =>  current_augs.extend(augmenters::chunks(&record, strategy_conf.count, strategy_conf.params[0], &mut base_aug_idx_for_record, rng)),
                        StrategyType::CenteredWindow => current_augs.extend(augmenters::centered_window(&record, strategy_conf.count, strategy_conf.params[0], &mut base_aug_idx_for_record, rng)),
                        StrategyType::BeginningTruncation => current_augs.extend(augmenters::beginning_truncations(&record, strategy_conf.count, strategy_conf.params[0], strategy_conf.params[1], &mut base_aug_idx_for_record, rng)),
                        StrategyType::EndTruncation => current_augs.extend(augmenters::end_truncations(&record, strategy_conf.count, strategy_conf.params[0], strategy_conf.params[1], &mut base_aug_idx_for_record, rng)),
                        StrategyType::RandomTruncation => current_augs.extend(augmenters::random_truncations(&record, strategy_conf.count, strategy_conf.params[0], strategy_conf.params[1], &mut base_aug_idx_for_record, rng)),
                        StrategyType::Synonymous => {
                            all_augmented_for_this_record.extend(augmenters::synonymous_variants(&record, strategy_conf.count, &genetic_code_arc, &mut base_aug_idx_for_record, rng));
                            continue; 
                        }
                    }
                    
                    if cli_config_arc.revcomp {
                        let mut rev_comp_augs = Vec::with_capacity(current_augs.len());
                        for aug_rec in &current_augs {
                            rev_comp_augs.push(FastaRecord {
                                header: format!("{}_rc", aug_rec.header),
                                sequence: augmenters::reverse_complement(&aug_rec.sequence),
                            });
                        }
                        current_augs.extend(rev_comp_augs);
                    }
                    all_augmented_for_this_record.extend(current_augs);
                }
                all_augmented_for_this_record
            }
        )
        .collect();

    // Use cli_config.output_path from the original (not Arc-wrapped) config
    let output_file = File::create(&cli_config.output_path) 
        .map_err(|e| format!("Failed to create output file {:?}: {}", cli_config.output_path, e))?;
    let mut writer = BufWriter::new(output_file);
    let output_line_width = 70;

    for batch in augmented_batches {
        for aug_rec in batch {
            fasta_io::write_fasta_record(&mut writer, &aug_rec, output_line_width)?;
        }
    }
    
    writer.flush()?;
    Ok(())
}