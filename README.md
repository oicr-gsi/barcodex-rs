# BarcodEX

[BarcodEX](https://github.com/oicr-gsi/barcodex) is tool for extracting Unique
Molecular Identifiers (UMIs) from single or paired-end read sequences.  It can
handle UMIs inline with reads or located in separate fastqs. This is a
reimplementation of the main package in Rust for improved performance when
handling large datasets.

## Installation
Obtain [Rust](https://www.rust-lang.org/tools/install) for your system. Then invoke:

    cargo build --release
    sudo cp target/release/barcodex-rs /usr/local/bin

You can copy the binary to any installation directory you choose.

## Running
BarcodEX requires input gzipped FASTQs. It can operate in two modes: the UMIs
and sequence are mixed in the same file (`--inline`) or the UMIs and sequence
are in separate files.

|Scenario                           | Arguments                                                                                                |
|-----------------------------------|----------------------------------------------------------------------------------------------------------|
| 1 read with an embedded UMI       | `--data single --inline --pattern1` _pattern_ `--r1_in` _fastq_                                          |
| 1 data read and 1 UMI read        | `--data single --r1_in` _fastq_ `--r2_in` _fastq_                                                        |
| 2 reads each with an embedded UMI | `--data paired --inline --pattern1` _pattern_ `--r1_in` _fastq_ `--pattern2` _pattern_ `--r2_in` _fastq_ |
| 2 data reads and 1 UMI read       | `--data paired --r1_in` _fastq_ `--r2_in` _fastq_ `--r3_in` _fastq_                                      |

In every case `--prefix` must be specified which specifies the start of the
output files. (_e.g._ `--prefix x` will result in `x_R1.fastq.gz`,
`x_R1.discarded.fastq.gz`, `x_R1.extracted.fastq.gz`, `x_UMI_counts.json`,
`x_extraction_metrics.json`).

The UMIs will only be accepted if they match an allow list provided with
`--umilist`. The list is a text file with one UMI per line. In the case of 2
reads with embedded UMIs, the two parts of the UMI must be separated by a `-`
(_e.g._, if the UMI on read 1 is `AAA` and the UMI on read 2 is `CCC`, then
`AAA-CCC` should be in the allow list)

The UMI will be placed in the header of the output file, separated by
`--separator` or `_` if unspecified.

## UMI Patterns
There are two ways to specify a pattern for the UMI: a nucleotide sequence or a
regular expression.

A nucleotide sequence is two or more `N`s followed by a spacer sequence.  For
instance the pattern `NNNNN` extracts the first 5 nucleotides from the read
whereas pattern `NNNNNATCG` extracts the first 9 nucleotides, using the first 5
nucleotides as the UMI and checks that the next 4 nucleotides match the spacer
`ATCG`.

Regular expressions allow more flexibility for extracting UMIs, in particular
UMIs with complex design and UMIs not starting at the beginning of the read.  A
good introduction to regular expression can be found in this [Regular
Expression HOWTO](https://docs.python.org/3/howto/regex.html). 

Sequences are extracted from the read using named groups within the regex.
Groups that have names that start with `umi` will be used for the UMI and
groups with names that start with `discard` will be matched but not included in
the output. Group names must be unique, so suffixing with `_` _number_ is
recommended.

For example, this expression extracts a 3bp UMI followed by TT spacer that is
removed from read and discarded:

     (?<umi>.{3})(?<discard>T{2})

Any sequence not contained in `umi` and `discard` groups will remain in the
read. Thus, it is important to construct the regular expression such that the
beginning of the read is captured in groups.

Normally, the regular expression is matched at the beginning of the read and
any unmatched bases at the end are assumed to be sequence. If the whole read
must be matched, use `--full_match`.
