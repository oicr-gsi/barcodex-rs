use clap::Parser;
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
struct OutputReadHandler<W: std::io::Write> {
    main_file: Option<bio::io::fastq::Writer<W>>,
    extracted_file: bio::io::fastq::Writer<W>,
    discard_file: bio::io::fastq::Writer<W>,
}
impl<'a> ExtractionMetrics<'a> {
    fn new(
        umi_list: &'a str,
        pattern1: Option<&'a str>,
        pattern2: Option<&'a str>,
    ) -> ExtractionMetrics<'a> {
        ExtractionMetrics {
            total: 0,
            matching: 0,
            discards: 0,
            unknown_umi: 0,
            acceptable_umi_list: umi_list,
            pattern1,
            pattern2,
        }
    }
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
                            read.qual()[end..read.qual().len()].into(),
                        ),
                        umi: read.seq()[0..*offset].into(),
                        extracted: (read.seq()[0..end].into(), read.qual()[0..end].into()),
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
impl<W: std::io::Write> OutputReadHandler<W> {
    /// Write this whole read to the discard file
    fn discard(self: &mut Self, read: &bio::io::fastq::Record) {
        self.discard_file.write_record(read).unwrap()
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
        self.read().id().as_bytes()
    }
    /// Get the input sequence associated which generated this output.
    ///
    /// Panics if there was no input
    fn read(self: &Self) -> &bio::io::fastq::Record {
        match self {
            Output::Empty => panic!("Trying to get sequence of empty record."),
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
fn extract_barcodes<R, W: std::io::Write>(
    inputs: &mut [InputReadHandler<R>],
    outputs: &mut [OutputReadHandler<W>],
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
                .join(&b'.');
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
/// Create a new output handler that writes to the main, discard, and extracted FASTQs
fn new_output(prefix: &str, read: usize, has_main: bool) -> OutputReadHandler<OutputFileWriter> {
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

#[derive(Parser)]
#[clap(
    version = "0.1",
    author = "Andre Masella <andre.masella@oicr.on.ca>",
    //help = "A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data"
)]
struct Opts {
    #[clap(long, help = "Path to file with valid UMIs (1st column)")]
    umilist: String,
    #[clap(long, help = "The prefix for output data files.")]
    prefix: String,
    #[clap(
        long,
        default_value = ":",
        help = "String separating the UMI sequence in the read name"
    )]
    separator: String,
    #[clap(subcommand)]
    subcmd: SubCommand,
}
#[derive(Parser)]
enum SubCommand {
    #[clap(
        name = "inline",
        //help = "Extract UMIs where the UMI is mixed with the sequence."
    )]
    Inline {
        #[clap(long, required = true, help = "Path to input FASTQ 1")]
        r1_in: Vec<String>,
        #[clap(long, help = "Path to input FASTQ 2")]
        r2_in: Vec<String>,
        #[clap(long, help = "Barcode string or regex for extracting UMIs in read 1")]
        pattern1: String,
        #[clap(long, help = "Barcode string or regex for extracting UMIs in read 2")]
        pattern2: Option<String>,
        #[clap(
            long,
            help = "Require the regex to inclreach to the enf the input read."
        )]
        full_match: bool,
    },
    #[clap(
        name = "separate",
        //help = "Extract UMIs where the UMI is in a separate read."
    )]
    Separate {
        #[clap(long, required = true, help = "Path to input FASTQ containing the UMI")]
        ru_in: Vec<String>,
        #[clap(long, required = true, help = "Path to input FASTQ 1")]
        r1_in: Vec<String>,
        #[clap(long, help = "Path to input FASTQ 2")]
        r2_in: Vec<String>,
    },
}

/// Extract allowed UMIs from a file for single-UMI situations
///
/// We're only expecting one UMI, so take lines that looke like `NNNN` or `NNNN 1` or `NNNN 1 2`
fn read_umi_list_single(umilist: &str) -> std::collections::HashMap<String, usize> {
    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let file = std::fs::File::open(umilist).unwrap();
    for line in std::io::BufReader::new(file).lines() {
        let l = line.unwrap();
        if l.split_whitespace().count() == 1 || l.split_whitespace().skip(1).any(|c| c == "1") {
            counts.insert(l.split_whitespace().next().unwrap().into(), 0);
        }
    }
    counts
}

/// Extract allowed UMIs for dual-UMI situations
///
/// We're expecting two UMIs, so treat `NNNN 1` as UMIs for read 1, `NNNN 2` as UMIs for read 2, and `NNNN` or `NNNN 1 2` as UMIs for both reads
fn read_umi_list_paired(umilist: &str) -> std::collections::HashMap<String, usize> {
    let mut umi1 = Vec::new();
    let mut umi2 = Vec::new();
    let file = std::fs::File::open(umilist).unwrap();
    for line in std::io::BufReader::new(file).lines() {
        let l = line.unwrap();
        if l.split_whitespace().count() == 1 {
            let sequence = l.split_whitespace().next().unwrap();
            umi1.push(sequence.to_string());
            umi2.push(sequence.to_string());
        }
        if l.split_whitespace().skip(1).any(|c| c == "1") {
            umi1.push(l.split_whitespace().next().unwrap().to_string());
        }
        if l.split_whitespace().skip(1).any(|c| c == "2") {
            umi2.push(l.split_whitespace().next().unwrap().to_string());
        }
    }
    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for u1 in &umi1 {
        for u2 in &umi2 {
            counts.insert(format!("{}.{}", u1, u2), 0);
        }
    }
    counts
}

fn main() {
    let args = Opts::parse();

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
                // Prepopulate the counts dictionary with all the UMIs in the supplied file; we will use these for validation later
                let mut counts = read_umi_list_single(&args.umilist);
                let mut metrics = ExtractionMetrics::new(&args.umilist, Some(&pattern1), None);
                for r1 in r1_in {
                    extract_barcodes(
                        &mut [InputReadHandler::Inline(read_fastq(&r1).records(), &p1)],
                        &mut [new_output(&args.prefix, 1, true)],
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
                // Prepopulate the counts dictionary with all the UMIs in the supplied file; we will use these for validation later
                let mut counts = read_umi_list_paired(&args.umilist);
                let pattern2_str = pattern2.expect("A pattern must be provided.");
                let p2 = InlineHandler::parse(&pattern2_str, full_match)
                    .expect("Valid pattern2 is required");
                let mut metrics =
                    ExtractionMetrics::new(&args.umilist, Some(&pattern1), Some(&pattern2_str));
                for (r1, r2) in r1_in.iter().zip(r2_in.iter()) {
                    extract_barcodes(
                        &mut [
                            InputReadHandler::Inline(read_fastq(r1).records(), &p1),
                            InputReadHandler::Inline(read_fastq(r2).records(), &p2),
                        ],
                        &mut [
                            new_output(&args.prefix, 1, true),
                            new_output(&args.prefix, 2, true),
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
            // Prepopulate the counts dictionary with all the UMIs in the supplied file; we will use these for validation later
            let mut counts = read_umi_list_single(&args.umilist);
            let mut metrics = ExtractionMetrics::new(&args.umilist, None, None);
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
                            new_output(&args.prefix, 1, true),
                            new_output(&args.prefix, 2, false),
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
                            new_output(&args.prefix, 1, true),
                            new_output(&args.prefix, 2, true),
                            new_output(&args.prefix, 3, false),
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
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    fn extract_separate<W: std::io::Write>(
        read1: &[u8],
        read2: &[u8],
        output1: &mut W,
        output1_ex: &mut W,
        output1_discard: &mut W,
        output2_ex: &mut W,
        output2_discard: &mut W,
        counts: &mut std::collections::HashMap<String, usize>,
    ) -> ExtractionMetrics<'static> {
        let mut metrics = ExtractionMetrics::new("/tmp/test", None, None);
        extract_barcodes(
            &mut [
                InputReadHandler::SequenceOnly(
                    bio::io::fastq::Reader::new(std::io::Cursor::new(read1)).records(),
                ),
                InputReadHandler::UmiOnly(
                    bio::io::fastq::Reader::new(std::io::Cursor::new(read2)).records(),
                ),
            ],
            &mut [
                OutputReadHandler {
                    main_file: Some(bio::io::fastq::Writer::new(output1)),
                    extracted_file: bio::io::fastq::Writer::new(output1_ex),
                    discard_file: bio::io::fastq::Writer::new(output1_discard),
                },
                OutputReadHandler {
                    main_file: None,
                    extracted_file: bio::io::fastq::Writer::new(output2_ex),
                    discard_file: bio::io::fastq::Writer::new(output2_discard),
                },
            ],
            ":",
            &mut metrics,
            counts,
        );
        metrics
    }
    fn extract_inline<W: std::io::Write>(
        read1: &[u8],
        handler1: &InlineHandler,
        read2: &[u8],
        handler2: &InlineHandler,
        output1: &mut W,
        output1_ex: &mut W,
        output1_discard: &mut W,
        output2: &mut W,
        output2_ex: &mut W,
        output2_discard: &mut W,
        counts: &mut std::collections::HashMap<String, usize>,
    ) -> ExtractionMetrics<'static> {
        let mut metrics = ExtractionMetrics::new("/tmp/test", None, None);
        extract_barcodes(
            &mut [
                InputReadHandler::Inline(
                    bio::io::fastq::Reader::new(std::io::Cursor::new(read1)).records(),
                    handler1,
                ),
                InputReadHandler::Inline(
                    bio::io::fastq::Reader::new(std::io::Cursor::new(read2)).records(),
                    handler2,
                ),
            ],
            &mut [
                OutputReadHandler {
                    main_file: Some(bio::io::fastq::Writer::new(output1)),
                    extracted_file: bio::io::fastq::Writer::new(output1_ex),
                    discard_file: bio::io::fastq::Writer::new(output1_discard),
                },
                OutputReadHandler {
                    main_file: Some(bio::io::fastq::Writer::new(output2)),
                    extracted_file: bio::io::fastq::Writer::new(output2_ex),
                    discard_file: bio::io::fastq::Writer::new(output2_discard),
                },
            ],
            ":",
            &mut metrics,
            counts,
        );
        metrics
    }

    #[test]
    fn test_s0() {
        let mut output1 = Vec::new();
        let mut output1_ex = Vec::new();
        let mut output1_discard = Vec::new();
        let mut output2_ex = Vec::new();
        let mut output2_discard = Vec::new();
        let mut counts = std::collections::HashMap::new();
        counts.insert("AAA".into(), 0);
        let results = extract_separate(
            b"",
            b"",
            &mut output1,
            &mut output1_ex,
            &mut output1_discard,
            &mut output2_ex,
            &mut output2_discard,
            &mut counts,
        );
        assert_eq!(results.total, 0);
        assert_eq!(results.discards, 0);
        assert_eq!(results.matching, 0);
        assert_eq!(results.unknown_umi, 0);
        assert_eq!(&output1, b"");
        assert_eq!(&output1_ex, b"");
        assert_eq!(&output1_discard, b"");
        assert_eq!(&output2_ex, b"");
        assert_eq!(&output2_discard, b"");
    }

    #[test]
    fn test_s1() {
        let mut output1 = Vec::new();
        let mut output1_ex = Vec::new();
        let mut output1_discard = Vec::new();
        let mut output2_ex = Vec::new();
        let mut output2_discard = Vec::new();
        let mut counts = std::collections::HashMap::new();
        counts.insert("AAA".into(), 0);
        let results = extract_separate(
            b"@M00000:0:000000000-00000:1:1:1:1 1:N:0\nACGTACGT\n+\nBBBBBBBB\n",
            b"@M00000:0:000000000-00000:1:1:1:1 2:N:0\nAAA\n+\nBBB\n",
            &mut output1,
            &mut output1_ex,
            &mut output1_discard,
            &mut output2_ex,
            &mut output2_discard,
            &mut counts,
        );
        assert_eq!(results.total, 1);
        assert_eq!(results.discards, 0);
        assert_eq!(results.matching, 1);
        assert_eq!(results.unknown_umi, 0);
        println!("{}", std::str::from_utf8(&output1).unwrap());
        println!("{}", std::str::from_utf8(&output2_ex).unwrap());
        assert!(output1
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AAA 1:N:0\nACGTACGT\n+\nBBBBBBBB\n".iter()));
        assert!(output1_ex
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AAA 1:N:0\n\n+\n\n".iter()));
        assert_eq!(&output1_discard, b"");
        assert!(output2_ex
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AAA 2:N:0\nAAA\n+\nBBB\n".iter()));
        assert_eq!(&output2_discard, b"");
    }
    #[test]
    fn test_s2() {
        let mut output1 = Vec::new();
        let mut output1_ex = Vec::new();
        let mut output1_discard = Vec::new();
        let mut output2_ex = Vec::new();
        let mut output2_discard = Vec::new();
        let mut counts = std::collections::HashMap::new();
        counts.insert("TTT".into(), 0);
        let results = extract_separate(
            b"@M00000:0:000000000-00000:1:1:1:1 1:N:0\nACGTACGT\n+\nBBBBBBBB\n",
            b"@M00000:0:000000000-00000:1:1:1:1 2:N:0\nAAA\n+\nBBB\n",
            &mut output1,
            &mut output1_ex,
            &mut output1_discard,
            &mut output2_ex,
            &mut output2_discard,
            &mut counts,
        );
        assert_eq!(results.total, 1);
        assert_eq!(results.discards, 0);
        assert_eq!(results.matching, 0);
        assert_eq!(results.unknown_umi, 1);
        assert_eq!(&output1, b"");
        assert_eq!(&output1_ex, b"");
        assert!(output1_discard
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1 1:N:0\nACGTACGT\n+\nBBBBBBBB\n".iter()));
        assert_eq!(&output2_ex, b"");
        assert!(output2_discard
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1 2:N:0\nAAA\n+\nBBB\n".iter()));
    }
    #[test]
    fn test_i0() {
        let mut output1 = Vec::new();
        let mut output1_ex = Vec::new();
        let mut output1_discard = Vec::new();
        let mut output2 = Vec::new();
        let mut output2_ex = Vec::new();
        let mut output2_discard = Vec::new();
        let mut counts = std::collections::HashMap::new();
        counts.insert("AAA".into(), 0);
        let results = extract_inline(
            b"",
            &InlineHandler::Nucleotide(2, "A".into()),
            b"",
            &InlineHandler::Nucleotide(2, "T".into()),
            &mut output1,
            &mut output1_ex,
            &mut output1_discard,
            &mut output2,
            &mut output2_ex,
            &mut output2_discard,
            &mut counts,
        );
        assert_eq!(results.total, 0);
        assert_eq!(results.discards, 0);
        assert_eq!(results.matching, 0);
        assert_eq!(results.unknown_umi, 0);
        assert_eq!(&output1, b"");
        assert_eq!(&output1_ex, b"");
        assert_eq!(&output1_discard, b"");
        assert_eq!(&output2, b"");
        assert_eq!(&output2_ex, b"");
        assert_eq!(&output2_discard, b"");
    }

    #[test]
    fn test_i1() {
        let mut output1 = Vec::new();
        let mut output1_ex = Vec::new();
        let mut output1_discard = Vec::new();
        let mut output2 = Vec::new();
        let mut output2_ex = Vec::new();
        let mut output2_discard = Vec::new();
        let mut counts = std::collections::HashMap::new();
        counts.insert("AC.AA".into(), 0);
        let results = extract_inline(
            b"@M00000:0:000000000-00000:1:1:1:1 1:N:0\nACGTACGT\n+\nBBBBBBBB\n",
            &InlineHandler::parse("(?P<umi>.{2})(?P<discard>G)", false).unwrap(),
            b"@M00000:0:000000000-00000:1:1:1:1 2:N:0\nAATGG\n+\nBBBBB\n",
            &InlineHandler::parse("NNT", false).unwrap(),
            &mut output1,
            &mut output1_ex,
            &mut output1_discard,
            &mut output2,
            &mut output2_ex,
            &mut output2_discard,
            &mut counts,
        );
        assert_eq!(results.total, 1);
        assert_eq!(results.discards, 0);
        assert_eq!(results.matching, 1);
        assert_eq!(results.unknown_umi, 0);
        println!("{}", std::str::from_utf8(&output1).unwrap());
        println!("{}", std::str::from_utf8(&output2_ex).unwrap());
        assert!(output1
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AC.AA 1:N:0\nTACGT\n+\nBBBBB\n".iter()));
        assert!(output1_ex
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AC.AA 1:N:0\nACG\n+\nBBB\n".iter()));
        assert_eq!(&output1_discard, b"");
        assert!(output2
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AC.AA 2:N:0\nGG\n+\nBB\n".iter()));
        assert!(output2_ex
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1:AC.AA 2:N:0\nAAT\n+\nBBB\n".iter()));
        assert_eq!(&output2_discard, b"");
    }
    #[test]
    fn test_i2() {
        let mut output1 = Vec::new();
        let mut output1_ex = Vec::new();
        let mut output1_discard = Vec::new();
        let mut output2 = Vec::new();
        let mut output2_ex = Vec::new();
        let mut output2_discard = Vec::new();
        let mut counts = std::collections::HashMap::new();
        counts.insert("TT".into(), 0);
        let results = extract_inline(
            b"@M00000:0:000000000-00000:1:1:1:1 1:N:0\nACGTACGT\n+\nBBBBBBBB\n",
            &InlineHandler::Nucleotide(2, "A".into()),
            b"@M00000:0:000000000-00000:1:1:1:1 2:N:0\nAATGG\n+\nBBBBB\n",
            &InlineHandler::Nucleotide(2, "T".into()),
            &mut output1,
            &mut output1_ex,
            &mut output1_discard,
            &mut output2,
            &mut output2_ex,
            &mut output2_discard,
            &mut counts,
        );
        assert_eq!(results.total, 1);
        assert_eq!(results.discards, 1);
        assert_eq!(results.matching, 0);
        assert_eq!(results.unknown_umi, 0);
        assert_eq!(&output1, b"");
        assert_eq!(&output1_ex, b"");
        assert!(output1_discard
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1 1:N:0\nACGTACGT\n+\nBBBBBBBB\n".iter()));
        assert_eq!(&output2_ex, b"");
        assert!(output2_discard
            .iter()
            .eq(b"@M00000:0:000000000-00000:1:1:1:1 2:N:0\nAATGG\n+\nBBBBB\n".iter()));
    }
}
