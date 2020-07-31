use bio_types::sequence::SequenceRead;
use clap::Clap;
use std::io::BufRead;

lazy_static::lazy_static! {
    /// The expression used to validate non-regex patterns
static ref UMI_PATTERN: regex::Regex = regex::Regex::new("^(N{2,})([ATCG]*)$").unwrap();
}

/// The metrics data structure written to provide statistics on processing
#[derive(serde::Serialize)]
struct ExtractionMetrics<'a> {
    #[serde(rename = "total reads/pairs")]
    total: usize,
    #[serde(rename = "reads/pairs with matching pattern")]
    matching: usize,
    #[serde(rename = "discarded reads/pairs")]
    discards: usize,
    #[serde(rename = "discarded reads/pairs due to unknown UMI")]
    unknown_umi: usize,
    #[serde(rename = "umi-list file")]
    acceptable_umi_list: &'a str,
    pattern1: Option<&'a str>,
    pattern2: Option<&'a str>,
}
/// The mechanism to process inline UMIs
enum InlineHandler {
    /// Contributes sequence and UMI and uses a leading UMI followed by a fixed spacer sequence
    Nucleotide(usize, String),
    /// Contributes sequence and UMI information and extracts the UMI using a regex
    Regex {
        /// The regular expression we are capturing
        regex: regex::bytes::Regex,
        /// Whether the user requested a full match or not
        full_match: bool,
        /// The capture groups corresponding to UMI regions
        umi_targets: Vec<usize>,
        /// The capture groups corresponding to UMI or discard regions (since both get written to the extracted file)
        extracted_targets: Vec<usize>,
    },
}
// The kind of input files we can handle
enum InputReadHandler<'a, R>
where
    R: std::io::Read,
{
    // This file contributes sequence and UMI
    Inline(bio::io::fastq::Records<R>, &'a InlineHandler),
    // This file contributes only sequence information
    SequenceOnly(bio::io::fastq::Records<R>),
    // This file contributes only UMI information
    UmiOnly(bio::io::fastq::Records<R>),
}
/// A pair of FASTQ nucleotides and their corresponding scores
type ScoredSequence = (Vec<u8>, Vec<u8>);

/// The result of processing a single record in a FASTQ file
#[derive(PartialEq)]
enum Output {
    /// There is no more data
    Empty,
    /// There is data, but it fails to match and should be junked
    Discard(bio::io::fastq::Record),
    /// Actual good data!
    Valid {
        input: bio::io::fastq::Record,
        umi: Vec<u8>,
        main: ScoredSequence,
        extracted: ScoredSequence,
    },
}
/// Type alias for the output file handler mess
type OutputFileWriter = flate2::write::GzEncoder<std::fs::File>;
/// Holder for FASTQ files that are being written to per read
struct OutputReadHandler {
    main_file: Option<bio::io::fastq::Writer<OutputFileWriter>>,
    extracted_file: bio::io::fastq::Writer<OutputFileWriter>,
    discard_file: bio::io::fastq::Writer<OutputFileWriter>,
}
impl InlineHandler {
    /// Generates output from an input sequence that will be processed in an appropriate way for the inline extraction method
    fn extract(&self, read: bio::io::fastq::Record) -> Output {
        match self {
            // If the input read is long enough, check that the spacer matches and put everything before the spacer into UMI and everything after into extracted sequence
            InlineHandler::Nucleotide(offset, spacer) => {
                let end = *offset + spacer.len();
                if end <= read.seq().len() && read.seq()[*offset..end] == *spacer.as_bytes() {
                    Output::Valid {
                        main: (
                            read.seq()[end..read.seq().len()].into(),
                            read.seq()[end..read.qual().len()].into(),
                        ),
                        umi: read.seq()[0..*offset].into(),
                        extracted: (read.seq()[0..end].into(), read.seq()[0..end].into()),
                        input: read,
                    }
                } else {
                    Output::Discard(read)
                }
            }
            // Process the read through the regular expression, anything in the `umi` capture groups becomes UMI, everything in `discard` is written to the discard file, and everything else is extracted sequence
            InlineHandler::Regex {
                regex,
                full_match,
                umi_targets,
                extracted_targets,
            } => match regex.captures(read.seq()) {
                Some(result) => {
                    // Check that some of the capture groups and something and, if we were asked for a full match, the regular expression makes it to the end of the sequece
                    if result
                        .iter()
                        .skip(1)
                        .all(|m| m.map(|mat| mat.as_bytes().is_empty()).unwrap_or(true))
                        || *full_match && result.get(0).unwrap().end() < read.seq().len()
                    {
                        Output::Discard(read)
                    } else {
                        let umi_indices = indices_from_regex(&result, umi_targets);
                        let extracted_indices = indices_from_regex(&result, extracted_targets);
                        Output::Valid {
                            main: (0..read.seq().len())
                                .filter(|i| extracted_indices.binary_search(i).is_err()) // Since the extracted region is a superset of the UMI, we don't need to search both
                                .map(|i| {
                                    (*read.seq().get(i).unwrap(), *read.qual().get(i).unwrap())
                                })
                                .unzip(),
                            extracted: (
                                extracted_indices
                                    .iter()
                                    .flat_map(|&i| read.seq().get(i))
                                    .copied()
                                    .collect(),
                                extracted_indices
                                    .iter()
                                    .flat_map(|&i| read.qual().get(i))
                                    .copied()
                                    .collect(),
                            ),
                            umi: umi_indices
                                .iter()
                                .flat_map(|&i| read.seq().get(i))
                                .copied()
                                .collect(),
                            input: read,
                        }
                    }
                }
                None => Output::Discard(read),
            },
        }
    }
    /// Try to find a fixed UMI pattern in an input string
    fn parse(pattern: &str, full_match: bool) -> Option<InlineHandler> {
        if let Some(captures) = UMI_PATTERN.captures(pattern) {
            Some(InlineHandler::Nucleotide(
                captures.get(1)?.end(),
                captures.get(2)?.as_str().into(),
            ))
        } else {
            let regex = regex::bytes::Regex::new(pattern).unwrap();
            Some(InlineHandler::Regex {
                umi_targets: regex
                    .capture_names()
                    .enumerate()
                    .filter(|(_, capture)| capture.map(|c| c.starts_with("umi")).unwrap_or(false))
                    .map(|(index, _)| index)
                    .collect(),
                extracted_targets: regex
                    .capture_names()
                    .enumerate()
                    .filter(|(_, capture)| {
                        capture
                            .map(|c| c.starts_with("umi") || c.starts_with("discard"))
                            .unwrap_or(false)
                    })
                    .map(|(index, _)| index)
                    .collect(),
                regex,
                full_match,
            })
        }
    }
}

impl<'a, R> InputReadHandler<'a, R>
where
    R: std::io::Read,
{
    /// Process a read for UMI content in a way that is appropriate for this file type
    fn next(self: &mut Self) -> Output {
        match self {
            // Extract both UMI and sequence
            InputReadHandler::Inline(reader, handler) => match reader.next() {
                None => Output::Empty,
                Some(record) => handler.extract(record.unwrap()),
            },
            // Dump everything into extracted sequence
            InputReadHandler::SequenceOnly(reader) => match reader.next() {
                None => Output::Empty,
                Some(record) => {
                    let read = record.unwrap();
                    Output::Valid {
                        umi: vec![],
                        main: (read.seq().into(), read.qual().into()),
                        extracted: (vec![], vec![]),
                        input: read,
                    }
                }
            },
            // Dump everything into UMI
            InputReadHandler::UmiOnly(reader) => match reader.next() {
                None => Output::Empty,
                Some(record) => {
                    let read = record.unwrap();
                    Output::Valid {
                        umi: read.seq().into(),
                        main: (vec![], vec![]),
                        extracted: (read.seq().into(), read.qual().into()),
                        input: read,
                    }
                }
            },
        }
    }
}
impl OutputReadHandler {
    /// Write this whole read to the discard file
    fn discard(self: &mut Self, read: &bio::io::fastq::Record) {
        self.discard_file.write_record(read).unwrap()
    }
    /// Create a new output handler that writes to the main, discard, and extracted FASTQs
    fn new(prefix: &str, read: usize, has_main: bool) -> OutputReadHandler {
        OutputReadHandler {
            main_file: if has_main {
                Some(write_fastq(prefix, read, ""))
            } else {
                None
            },
            extracted_file: write_fastq(prefix, read, ".extracted"),
            discard_file: write_fastq(prefix, read, ".discarded"),
        }
    }
    /// Write this processed read to the main and extracted output files
    fn write(
        self: &mut Self,
        header: &[u8],
        desc: Option<&str>,
        main: &ScoredSequence,
        extracted: &ScoredSequence,
    ) {
        if let Some(main_file) = self.main_file.as_mut() {
            main_file
                .write(std::str::from_utf8(header).unwrap(), desc, &main.0, &main.1)
                .unwrap();
        }
        self.extracted_file
            .write(
                std::str::from_utf8(header).unwrap(),
                desc,
                &extracted.0,
                &extracted.1,
            )
            .unwrap();
    }
}

impl Output {
    /// Get the name of the input sequence associated which generated this output.
    ///
    /// Panics if there was no input
    fn name(self: &Self) -> &[u8] {
        self.read().name()
    }
    /// Get the input sequence associated which generated this output.
    ///
    /// Panics if there was no input
    fn read(self: &Self) -> &bio::io::fastq::Record {
        match self {
            Output::Empty => panic!("Trying to get seqeuence of empty record."),
            Output::Discard(record) => record,
            Output::Valid {
                input,
                umi: _,
                main: _,
                extracted: _,
            } => input,
        }
    }
}
/// Take a bunch of input files and process all the reads
fn extract_barcodes<R>(
    inputs: &mut [InputReadHandler<R>],
    outputs: &mut [OutputReadHandler],
    separator: &str,
    metrics: &mut ExtractionMetrics,
    counts: &mut std::collections::HashMap<String, usize>,
) where
    R: std::io::Read,
{
    loop {
        // For every input file, try to extract a read and process it
        let reads: Vec<Output> = inputs.iter_mut().map(|handler| handler.next()).collect();
        // See how many files are out of reads
        let empty_files = reads
            .iter()
            .filter(|&output| *output == Output::Empty)
            .count();
        // If all of them are done, we're done
        if empty_files == inputs.len() {
            return;
        }
        // If only some, then we have a problem
        if empty_files > 0 {
            panic!("Reads are unavailable from one or more input files")
        }
        // Check that the headers match up
        if reads.iter().skip(1).any(|o| reads[0].name() != o.name()) {
            panic!("Reads are desynchronized")
        }
        metrics.total += 1;
        // Every file gets to vote what to do next; if anyone thinks their read is garbage, we reject this read
        if reads.iter().any(|output| match output {
            Output::Discard(_) => true,
            _ => false,
        }) {
            metrics.discards += 1;
            for (index, output) in outputs.iter_mut().enumerate() {
                if let Some(o) = reads.get(index) {
                    output.discard(o.read());
                }
            }
        } else {
            // Combine all the bits of UMI we have found across all the reads and combine them
            let umi: Vec<u8> = reads
                .iter()
                .map(|o| match o {
                    Output::Valid {
                        umi,
                        extracted: _,
                        main: _,
                        input: _,
                    } => umi.as_slice(),
                    _ => panic!("Unexpected invalid record"),
                })
                .filter(|u| !u.is_empty())
                .collect::<Vec<&[u8]>>()
                .join(&b'-');
            let umi_str: String = std::str::from_utf8(&umi).unwrap().into();
            // If our counts hashmap has this UMI, then increment it and output the modified read
            match counts.get_mut(&umi_str) {
                Some(count) => {
                    *count += 1;
                    metrics.matching += 1;
                    let header = &[&reads[0].name(), separator.as_bytes(), &umi].concat();
                    for (index, output) in outputs.iter_mut().enumerate() {
                        if let Some(Output::Valid {
                            umi: _,
                            extracted,
                            main,
                            input,
                        }) = reads.get(index)
                        {
                            output.write(&header, input.desc(), main, extracted);
                        }
                    }
                }
                // Otherwise, this UMI is unknown to us and we will discard these reads
                None => {
                    metrics.unknown_umi += 1;
                    for (index, output) in outputs.iter_mut().enumerate() {
                        if let Some(o) = reads.get(index) {
                            output.discard(o.read());
                        }
                    }
                }
            }
        }
    }
}
/// Create a list of indices of all the positions that were part of a regular expression capture group
fn indices_from_regex(result: &regex::bytes::Captures, target_captures: &[usize]) -> Vec<usize> {
    let mut indices: Vec<usize> = target_captures
        .iter()
        .flat_map(|target| result.get(*target))
        .flat_map(|r| r.range())
        .collect();
    indices.sort();
    indices
}

/// Read a gzipped FASTQ
fn read_fastq(
    path: &str,
) -> bio::io::fastq::Reader<flate2::bufread::MultiGzDecoder<std::io::BufReader<std::fs::File>>> {
    std::fs::File::open(path)
        .map(std::io::BufReader::new)
        .map(flate2::bufread::MultiGzDecoder::new)
        .map(bio::io::fastq::Reader::new)
        .unwrap()
}

/// Create a new gzipped output FASTQ
fn write_fastq(
    prefix: &str,
    read: usize,
    suffix: &str,
) -> bio::io::fastq::Writer<OutputFileWriter> {
    std::fs::File::create(format!("{}_R{}{}.fastq.gz", prefix, read, suffix))
        .map(|w| flate2::write::GzEncoder::new(w, flate2::Compression::best()))
        .map(bio::io::fastq::Writer::new)
        .unwrap()
}

/// Write out metrics and counts to JSON files
fn write_stats(
    prefix: &str,
    metrics: &ExtractionMetrics,
    counts: &mut std::collections::HashMap<String, usize>,
) {
    // We populate the hash map with all allowed UMIs; drop any from the output that weren't found
    counts.retain(|_, &mut v| v > 0);
    serde_json::to_writer(
        std::fs::File::create(format!("{}_UMI_counts.json", prefix)).unwrap(),
        counts,
    )
    .unwrap();
    serde_json::to_writer(
        std::fs::File::create(format!("{}_extraction_metrics.json", prefix)).unwrap(),
        metrics,
    )
    .unwrap();
}

#[derive(clap::Clap)]
#[clap(
    version = "0.1",
    author = "Andre Masella <andre.masella@oicr.on.ca>",
    about = "A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data"
)]
struct Opts {
    #[clap(long, about = "Path to file with valid UMIs (1st column)")]
    umilist: String,
    #[clap(long, about = "The prefix for output data files.")]
    prefix: String,
    #[clap(
        long,
        default_value = "_",
        about = "String separating the UMI sequence in the read name"
    )]
    separator: String,
    #[clap(subcommand)]
    subcmd: SubCommand,
}
#[derive(Clap)]
enum SubCommand {
    #[clap(name = "inline")]
    Inline {
        #[clap(long, required = true, about = "Path to input FASTQ 1")]
        r1_in: Vec<String>,
        #[clap(long, about = "Path to input FASTQ 2")]
        r2_in: Vec<String>,
        #[clap(long, about = "Barcode string or regex for extracting UMIs in read 1")]
        pattern1: String,
        #[clap(long, about = "Barcode string or regex for extracting UMIs in read 2")]
        pattern2: Option<String>,
        #[clap(
            long,
            about = "Require the regex to inclreach to the enf the input read."
        )]
        full_match: bool,
    },
    #[clap(name = "separate")]
    Separate {
        #[clap(
            long,
            required = true,
            about = "Path to input FASTQ containing the UMI"
        )]
        ru_in: Vec<String>,
        #[clap(long, required = true, about = "Path to input FASTQ 1")]
        r1_in: Vec<String>,
        #[clap(long, about = "Path to input FASTQ 2")]
        r2_in: Vec<String>,
    },
}

fn main() {
    let args = Opts::parse();
    // Prepopulate the counts dictionary with all the UMIs in the supplied file; we will use these for validation later
    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let file = std::fs::File::open(&args.umilist).unwrap();
    for line in std::io::BufReader::new(file).lines() {
        counts.insert(line.unwrap(), 0);
    }
    // Figure out which case we are in, because there are so many
    std::process::exit(match args.subcmd {
        SubCommand::Inline {
            r1_in,
            r2_in,
            pattern1,
            pattern2,
            full_match,
        } => {
            let p1 = InlineHandler::parse(&pattern1, full_match).expect("Cannot parse pattern 1");
            if r2_in.is_empty() {
                let mut metrics = ExtractionMetrics {
                    total: 0,
                    matching: 0,
                    discards: 0,
                    unknown_umi: 0,
                    acceptable_umi_list: &args.umilist,
                    pattern1: Some(&pattern1),
                    pattern2: None,
                };
                for r1 in r1_in {
                    extract_barcodes(
                        &mut [InputReadHandler::Inline(read_fastq(&r1).records(), &p1)],
                        &mut [OutputReadHandler::new(&args.prefix, 1, true)],
                        &args.separator,
                        &mut metrics,
                        &mut counts,
                    );
                }
                write_stats(&args.prefix, &metrics, &mut counts);
                0
            } else if r2_in.len() != r1_in.len() {
                eprintln!("Number of R1 files does not match R2 files.");
                1
            } else {
                let pattern2_str = pattern2.expect("A pattern must be provided.");
                let p2 = InlineHandler::parse(&pattern2_str, full_match)
                    .expect("Valid pattern2 is required");
                let mut metrics = ExtractionMetrics {
                    total: 0,
                    matching: 0,
                    discards: 0,
                    unknown_umi: 0,
                    acceptable_umi_list: &args.umilist,
                    pattern1: Some(&pattern1),
                    pattern2: Some(&pattern2_str),
                };
                for (r1, r2) in r1_in.iter().zip(r2_in.iter()) {
                    extract_barcodes(
                        &mut [
                            InputReadHandler::Inline(read_fastq(r1).records(), &p1),
                            InputReadHandler::Inline(read_fastq(r2).records(), &p2),
                        ],
                        &mut [
                            OutputReadHandler::new(&args.prefix, 1, true),
                            OutputReadHandler::new(&args.prefix, 2, true),
                        ],
                        &args.separator,
                        &mut metrics,
                        &mut counts,
                    );
                }
                write_stats(&args.prefix, &metrics, &mut counts);
                0
            }
        }
        SubCommand::Separate {
            ru_in,
            r1_in,
            r2_in,
        } => {
            let mut metrics = ExtractionMetrics {
                total: 0,
                matching: 0,
                discards: 0,
                unknown_umi: 0,
                acceptable_umi_list: &args.umilist,
                pattern1: None,
                pattern2: None,
            };
            if ru_in.len() != r1_in.len() {
                eprintln!("Mismatched file counts for R1 and UMIs.");
                1
            } else if r2_in.is_empty() {
                for (r1, ru) in r1_in.iter().zip(ru_in.iter()) {
                    extract_barcodes(
                        &mut [
                            InputReadHandler::SequenceOnly(read_fastq(r1).records()),
                            InputReadHandler::UmiOnly(read_fastq(ru).records()),
                        ],
                        &mut [
                            OutputReadHandler::new(&args.prefix, 1, true),
                            OutputReadHandler::new(&args.prefix, 2, false),
                        ],
                        &args.separator,
                        &mut metrics,
                        &mut counts,
                    );
                }
                write_stats(&args.prefix, &metrics, &mut counts);
                0
            } else if r2_in.len() != ru_in.len() {
                eprintln!("Mismatched file counts for R2 and UMIs.");
                1
            } else {
                for ((r1, r2), ru) in r1_in.iter().zip(r2_in.iter()).zip(ru_in.iter()) {
                    extract_barcodes(
                        &mut [
                            InputReadHandler::SequenceOnly(read_fastq(r1).records()),
                            InputReadHandler::SequenceOnly(read_fastq(r2).records()),
                            InputReadHandler::UmiOnly(read_fastq(ru).records()),
                        ],
                        &mut [
                            OutputReadHandler::new(&args.prefix, 1, true),
                            OutputReadHandler::new(&args.prefix, 2, true),
                            OutputReadHandler::new(&args.prefix, 3, false),
                        ],
                        &args.separator,
                        &mut metrics,
                        &mut counts,
                    );
                }
                write_stats(&args.prefix, &metrics, &mut counts);
                0
            }
        }
    })
}
