use bio_types::sequence::SequenceRead;
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
    pattern1: &'a Option<String>,
    pattern2: &'a Option<String>,
}
/// The mechanism to process inline UMIs
enum InlineHandler {
    // Contributes sequence and UMI and uses a leading UMI followed by a fixed spacer sequence
    Nucleotide(usize, String),
    // Contributes sequence and UMI information and extracts the UMI using a regex
    Regex(regex::bytes::Regex, bool),
}
// The kind of input files we can handle
enum InputReadHandler<'a> {
    // This file contributes sequence and UMI
    Inline(&'a InlineHandler),
    // This file contributes only sequence information
    SequenceOnly,
    // This file contributes only UMI information
    UmiOnly,
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
struct OutputReadHandler<W>
where
    W: std::io::Write,
{
    main_file: bio::io::fastq::Writer<W>,
    extracted_file: bio::io::fastq::Writer<W>,
    discard_file: bio::io::fastq::Writer<W>,
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
            InlineHandler::Regex(regex, full_match) => match regex.captures(read.seq()) {
                Some(result) => {
                    if *full_match && result.get(0).unwrap().end() < read.seq().len() {
                        Output::Discard(read)
                    } else {
                        let umi_indices = indices_from_regex(regex, &result, "umi");
                        let mut discard_indices = indices_from_regex(regex, &result, "discard");
                        discard_indices.extend_from_slice(&umi_indices);
                        discard_indices.sort();
                        Output::Valid {
                            main: (0..read.seq().len())
                                .filter(|i| {
                                    !umi_indices.contains(i) && !discard_indices.contains(i)
                                })
                                .map(|i| {
                                    (*read.seq().get(i).unwrap(), *read.qual().get(i).unwrap())
                                })
                                .unzip(),
                            extracted: (
                                discard_indices
                                    .iter()
                                    .flat_map(|&i| read.seq().get(i))
                                    .map(|&v| v)
                                    .collect(),
                                discard_indices
                                    .iter()
                                    .flat_map(|&i| read.qual().get(i))
                                    .map(|&v| v)
                                    .collect(),
                            ),
                            umi: umi_indices
                                .iter()
                                .flat_map(|&i| read.seq().get(i))
                                .map(|&v| v)
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
            Some(InlineHandler::Regex(
                regex::bytes::Regex::new(pattern).unwrap(),
                full_match,
            ))
        }
    }
}

impl<'a> InputReadHandler<'a> {
    /// Process a read for UMI content in a way that is appropriate for this file type
    fn process(self: &mut Self, read: bio::io::fastq::Record) -> Output {
        match self {
            // Extract both UMI and sequence
            InputReadHandler::Inline(handler) => handler.extract(read),
            // Dump everything into extracted sequence
            InputReadHandler::SequenceOnly => Output::Valid {
                umi: vec![],
                main: (read.seq().into(), read.qual().into()),
                extracted: (vec![], vec![]),
                input: read,
            },
            // Dump everything into UMI
            InputReadHandler::UmiOnly => Output::Valid {
                umi: read.seq().into(),
                main: (vec![], vec![]),
                extracted: (vec![], vec![]),
                input: read,
            },
        }
    }
}
impl<W> OutputReadHandler<W>
where
    W: std::io::Write,
{
    /// Write this whole read to the discard file
    fn discard(self: &mut Self, read: &bio::io::fastq::Record) -> () {
        self.discard_file.write_record(read).unwrap()
    }
    /// Write this processed read to the main and extracted output files
    fn write(
        self: &mut Self,
        header: &[u8],
        desc: Option<&str>,
        main: &ScoredSequence,
        extracted: &ScoredSequence,
    ) -> () {
        self.main_file
            .write(std::str::from_utf8(header).unwrap(), desc, &main.0, &main.1)
            .unwrap();
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
    fn name<'a>(self: &'a Self) -> &'a [u8] {
        self.read().name()
    }
    /// Get the input sequence associated which generated this output.
    ///
    /// Panics if there was no input
    fn read<'a>(self: &'a Self) -> &'a bio::io::fastq::Record {
        match self {
            Output::Empty => panic!("Trying to get c of empty record."),
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
fn extract_barcodes<R, W>(
    inputs: &mut [(bio::io::fastq::Records<R>, InputReadHandler)],
    outputs: &mut [OutputReadHandler<W>],
    separator: &str,
    metrics: &mut ExtractionMetrics,
    counts: &mut std::collections::HashMap<String, usize>,
) -> ()
where
    R: std::io::Read,
    W: std::io::Write,
{
    loop {
        // For every input file, try to extract a read and process it
        let reads: Vec<Output> = inputs
            .iter_mut()
            .map(|(reader, handler)| match reader.next() {
                None => Output::Empty,
                Some(record) => handler.process(record.unwrap()),
            })
            .collect();
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
                match reads.get(index) {
                    Some(o) => output.discard(o.read()),
                    None => (),
                }
            }
        } else {
            // Combine all the bits of UMI we have found across all the reads and combine them
            let umi: Vec<u8> = reads
                .iter()
                .flat_map(|o| match o {
                    Output::Valid {
                        umi,
                        extracted: _,
                        main: _,
                        input: _,
                    } => umi.into_iter(),
                    _ => panic!("Unexpected invalid record"),
                })
                .map(|v| *v)
                .collect();
            let umi_str: String = std::str::from_utf8(&umi).unwrap().into();
            // If our counts hashmap has this UMI, then increment it and output the modified read
            match counts.get_mut(&umi_str) {
                Some(count) => {
                    *count += 1;
                    metrics.matching += 1;
                    let header = &[&reads[0].name(), separator.as_bytes(), &umi].concat();
                    for (index, output) in outputs.iter_mut().enumerate() {
                        match reads.get(index) {
                            Some(Output::Valid {
                                umi: _,
                                extracted,
                                main,
                                input,
                            }) => output.write(&header, input.desc(), main, extracted),
                            _ => (),
                        }
                    }
                }
                // Otherwise, this UMI is unknown to us and we will discard these reads
                None => {
                    metrics.unknown_umi += 1;
                    for (index, output) in outputs.iter_mut().enumerate() {
                        match reads.get(index) {
                            Some(o) => output.discard(o.read()),
                            None => (),
                        }
                    }
                }
            }
        }
    }
}
/// Create a list of indices of all the positions that were part of a regular expression capture group with a name matching the prefix provided
fn indices_from_regex(
    regex: &regex::bytes::Regex,
    result: &regex::bytes::Captures,
    prefix: &str,
) -> Vec<usize> {
    let mut indices: Vec<usize> = regex
        .capture_names()
        .filter(|name| {
            name.map(|name_str| name_str.starts_with(prefix))
                .unwrap_or(false)
        })
        .flatten()
        .flat_map(|name| result.name(name))
        .flat_map(|r| r.range().into_iter())
        .collect();
    indices.sort();
    indices
}

/// Create a new output handler that writes to the main, discard, and extracted FASTQs
fn new_output(prefix: &Option<String>, read: usize) -> OutputReadHandler<OutputFileWriter> {
    OutputReadHandler {
        main_file: write_fastq(prefix, read, ""),
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
    prefix: &Option<String>,
    read: usize,
    suffix: &str,
) -> bio::io::fastq::Writer<OutputFileWriter> {
    std::fs::File::create(format!(
        "{}_R{}{}.fastq.gz",
        prefix.as_ref().expect("--prefix is required"),
        read,
        suffix
    ))
    .map(|w| flate2::write::GzEncoder::new(w, flate2::Compression::best()))
    .map(bio::io::fastq::Writer::new)
    .unwrap()
}

/// Write out metrics and counts to JSON files
fn write_stats(
    prefix: &str,
    metrics: &ExtractionMetrics,
    counts: &mut std::collections::HashMap<String, usize>,
) -> () {
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

fn main() {
    let mut data: String = "single".into();
    let mut full_match = false;
    let mut inline = false;
    let mut pattern1: Option<String> = None;
    let mut pattern2: Option<String> = None;
    let mut prefix: Option<String> = None;
    let mut read1: Vec<String> = vec![];
    let mut read2: Vec<String> = vec![];
    let mut read3: Vec<String> = vec![];
    let mut separator: String = "_".into();
    let mut umi_list: Option<String> = None;
    {
        let mut parser = argparse::ArgumentParser::new();
        parser.set_description("A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data");
        parser.refer(&mut read1).add_option(
            &["--r1_in"],
            argparse::ParseList,
            "Path to input FASTQ 1",
        );
        parser.refer(&mut read2).add_option(
            &["--r2_in"],
            argparse::ParseList,
            "Path to input FASTQ 2",
        );
        parser.refer(&mut read3).add_option(
            &["--r3_in"],
            argparse::ParseList,
            "Path to input FASTQ 3",
        );
        parser.refer(&mut pattern1).add_option(
            &["--pattern1"],
            argparse::StoreOption,
            "Barcode string of regex for extracting UMIs in read 1",
        );
        parser.refer(&mut pattern2).add_option(
            &["--pattern2"],
            argparse::StoreOption,
            "Barcode string of regex for extracting UMIs in read 2",
        );
        parser.refer(&mut inline).add_option(
            &["--inline"],
            argparse::StoreTrue,
            "UMIs inline with reads or not. True if activated",
        );
        parser.refer(&mut prefix).add_option(
            &["--prefix"],
            argparse::StoreOption,
            "The prefix for output data files.",
        );
        parser.refer(&mut data).add_option(
            &["--data"],
            argparse::Store,
            "Paired or single end sequencing",
        );
        parser.refer(&mut separator).add_option(
            &["--separator"],
            argparse::Store,
            "String separating the UMI sequence in the read name",
        );
        parser.refer(&mut full_match).add_option(
            &["--full_match"],
            argparse::StoreTrue,
            "Requires the regex pattern to match the entire read sequence. True if activated",
        );
        parser.refer(&mut umi_list).add_option(
            &["--umilist"],
            argparse::ParseOption,
            "Path to file with valid UMIs (1st column)",
        );
        parser.parse_args_or_exit();
    }

    // Prepopulate the counts dictionary with all the UMIs in the supplied file; we will use these for validation later
    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let file = std::fs::File::open(
        umi_list
            .as_ref()
            .expect("A list of allowed UMIs must be provided."),
    )
    .unwrap();
    for line in std::io::BufReader::new(file).lines() {
        counts.insert(line.unwrap(), 0);
    }
    let mut metrics = ExtractionMetrics {
        total: 0,
        matching: 0,
        discards: 0,
        unknown_umi: 0,
        acceptable_umi_list: &umi_list.as_ref().unwrap(),
        pattern1: &pattern1,
        pattern2: &pattern2,
    };

    // Figure out which case we are in, because there are so many
    std::process::exit(match data.as_str() {
        "single" => {
            if !read3.is_empty() {
                eprint!("Single ended data requested but R3 files provided.");
                1
            } else if inline {
                if !read2.is_empty() {
                    eprint!("Single ended data requested but R2 files provided.");
                    1
                } else if pattern2 != None {
                    eprint!("Single ended data requested but pattern2 provided.");
                    1
                } else {
                    let pattern = InlineHandler::parse(
                        &pattern1.as_ref().expect("--pattern1 is required"),
                        full_match,
                    )
                    .expect("Valid pattern1 is required");
                    for r1 in read1 {
                        extract_barcodes(
                            &mut [(
                                read_fastq(&r1).records(),
                                InputReadHandler::Inline(&pattern),
                            )],
                            &mut [new_output(&prefix, 1)],
                            &separator,
                            &mut metrics,
                            &mut counts,
                        );
                    }
                    write_stats(&prefix.unwrap(), &metrics, &mut counts);
                    0
                }
            } else if read1.len() != read2.len() {
                eprint!("Mismatched file counts for R1 and R2.");
                1
            } else if pattern1 != None || pattern2 != None {
                eprint!("Patterns are provided for non-inline UMI.");
                1
            } else {
                for (r1, r2) in read1.iter().zip(read2.iter()) {
                    extract_barcodes(
                        &mut [
                            (read_fastq(&r1).records(), InputReadHandler::SequenceOnly),
                            (read_fastq(&r2).records(), InputReadHandler::UmiOnly),
                        ],
                        &mut [new_output(&prefix, 1)],
                        &separator,
                        &mut metrics,
                        &mut counts,
                    );
                }
                write_stats(&prefix.unwrap(), &metrics, &mut counts);
                0
            }
        }
        "paired" => {
            if inline {
                if !read3.is_empty() {
                    eprint!("Paired end inline data requested but R3 files provided.");
                    1
                } else if read1.len() != read2.len() {
                    eprint!("Mismatched file counts for R1 and R2.");
                    1
                } else {
                    let p1 = InlineHandler::parse(
                        &pattern1.as_ref().expect("--pattern1 is required"),
                        full_match,
                    )
                    .expect("Valid pattern1 is required");
                    let p2 = InlineHandler::parse(
                        &pattern2.as_ref().expect("--pattern2 is required"),
                        full_match,
                    )
                    .expect("Valid pattern2 is required");
                    for (r1, r2) in read1.iter().zip(read2.iter()) {
                        extract_barcodes(
                            &mut [
                                (read_fastq(&r1).records(), InputReadHandler::Inline(&p1)),
                                (read_fastq(&r2).records(), InputReadHandler::Inline(&p2)),
                            ],
                            &mut [new_output(&prefix, 1), new_output(&prefix, 2)],
                            &separator,
                            &mut metrics,
                            &mut counts,
                        );
                    }
                    write_stats(&prefix.unwrap(), &metrics, &mut counts);
                    0
                }
            } else if read1.len() != read2.len() && read1.len() != read3.len() {
                eprint!("Mismatched file counts for R1, R2, and R3.");
                1
            } else if pattern1 != None || pattern2 != None {
                eprint!("Patterns are provided for non-inline UMI.");
                1
            } else {
                for ((r1, r2), r3) in read1.iter().zip(read2.iter()).zip(read3.iter()) {
                    extract_barcodes(
                        &mut [
                            (read_fastq(&r1).records(), InputReadHandler::SequenceOnly),
                            (read_fastq(&r2).records(), InputReadHandler::SequenceOnly),
                            (read_fastq(&r3).records(), InputReadHandler::UmiOnly),
                        ],
                        &mut [new_output(&prefix, 1), new_output(&prefix, 2)],
                        &separator,
                        &mut metrics,
                        &mut counts,
                    );
                }
                write_stats(&prefix.unwrap(), &metrics, &mut counts);
                0
            }
        }
        _ => {
            eprint!(
                "Unknown option for --data: {} should be single or pairied",
                data
            );
            1
        }
    })
}
