# SHI7 (Short-read Iterative Trimmer) Documentation

## Overview

SHI7 (Short-read Iterative Trimmer) is a powerful and flexible tool designed for quality control and trimming of short-read DNA sequences, with a focus on microbiome data. This C-based program offers a wide range of features for processing both single-end and paired-end reads, making it suitable for various sequencing platforms and experimental designs.

## Version

This documentation covers SHI7 version 0.92f.

## Features

1. Native gzip support
2. Adapter autodetection and trimming
3. Support for newer Illumina 2-color instruments (poly-G repeat removal)
4. PALINCUT mode for detection and removal of novel adapters or suspected technical sequences (paired-end only)
5. Header stripping/anonymizing
6. Maximum depth output (rarefaction)
7. FASTA and FASTQ output format support
8. Illumina machine discard signal detection
9. Various quality trimming modes
10. Different trimming quality thresholds for read start and end
11. High-speed processing

## Compilation

The provided code is a standalone C program. To compile it, use a C compiler such as GCC:

```
gcc -O3 -o shi7_trimmer shi7_trimmer.c -lz
```

Note: The `-lz` flag is necessary to link against the zlib library for gzip support.

## Usage

```
shi7_trimmer in_seqs.fastq out_prefix MINLEN QUAL <QUAL_RIGHT> [options]
```

### Required Arguments

1. `in_seqs.fastq`: Input FASTQ file (can be gzipped)
2. `out_prefix`: Prefix for output files
3. `MINLEN`: Minimum allowed read length after trimming
4. `QUAL`: Quality threshold for trimming (left side of read)
5. `<QUAL_RIGHT>`: Quality threshold for trimming (right side of read, optional)

### Optional Arguments

- `[ROLLING X [PREROLL/POSTROLL]]`: Use rolling average over X bases
  - `PREROLL`: Allow initial windows < X
  - `POSTROLL`: Add back X bases to both sides after trimming
- `[FLOOR X]`: Cut until X consecutive bases have minimum quality
- `[CUT]`: Treat 'QUAL' and 'QUAL_RIGHT' as number of bases to blindly cut
- `[ASS_QUALITY X]`: Discard reads with average quality < X (after trimming)
- `[N-IFY X]`: Convert bases with quality <= X to N
- `[CASTN X]`: Convert N bases to quality score X
- `[DITCH X]`: Discard reads with any quality score below X
- `[STRIP]`: Replace header names with numerical indices
- `[HEAD X]`: Only output the first X good sequences
- `[AUTOADAP]`: Automatically detect and remove common adapter sequences and tails
- `[ADAP ...]`: Manually specify adapter sequence to remove (up to 5 times)
- `[ADAP2 ...]`: Manually specify paired partial adapter sequence (up to 5 times)
- `[PALINCUT X]`: Auto-detect and remove paired-end adapters (X total - 2 end)
- `[R2 in_seqsR2.fastq]`: Specify R2 file for paired-end mode
- `[OUTFASTA]`: Output in FASTA format instead of FASTQ

Note: Add 31 to all quality arguments if using PHRED64 encoding.

## Trimming Methods

SHI7 offers several trimming methods:

1. Default: Trim until a base with quality > threshold is encountered
2. Rolling Average: Use a rolling average of quality scores
3. Floor: Trim until X consecutive bases meet the quality threshold
4. Cut: Blindly cut a specified number of bases from each end

## Adapter Removal

SHI7 provides multiple options for adapter removal:

1. Automatic adapter detection (AUTOADAP)
2. Manual specification of up to 5 adapter sequences (ADAP)
3. Paired-end partial adapter removal (ADAP2)
4. Auto-detection of paired-end adapters (PALINCUT)

## Output

SHI7 generates trimmed and quality-controlled sequences in either FASTQ (default) or FASTA format. For paired-end data, it produces separate output files for R1 and R2 reads.

## Performance Considerations

1. The program uses buffered I/O for improved performance.
2. It allocates a maximum read length of 16MB, which should be sufficient for most sequencing technologies.
3. The code includes optimizations for speed, such as using rolling calculations for quality assessments.

## Error Handling

SHI7 includes basic error handling:

1. It checks for correct number of arguments (5-25).
2. It verifies that input and output files can be opened successfully.
3. It provides warnings for incompatible option combinations (e.g., PALINCUT with single-end data).

## Limitations

1. The maximum read length is hard-coded to 16MB.
2. The program assumes PHRED33 quality encoding by default.
3. There's no built-in multi-threading support.

## Example Usage

1. Basic usage (single-end):
   ```
   shi7_trimmer in_seqs.fastq out_prefix 50 20
   ```
   This will trim reads to a minimum length of 50bp, cutting bases with quality < 20.

2. Paired-end with automatic adapter removal:
   ```
   shi7_trimmer in_seqs_R1.fastq out_prefix 100 30 30 AUTOADAP R2 in_seqs_R2.fastq
   ```
   This will process paired-end data, trimming to a minimum length of 100bp, using a quality threshold of 30 on both ends, and automatically removing adapters.

3. Using rolling average quality:
   ```
   shi7_trimmer in_seqs.fastq out_prefix 75 25 ROLLING 10 PREROLL
   ```
   This will use a rolling average of 10 bases for quality assessment, allowing initial windows smaller than 10 bases.

## Output Statistics

After processing, SHI7 provides summary statistics including:

1. Number of low-quality sequences discarded
2. Number of sequences that passed QC
3. Average trimmed read length
4. Average number of bases trimmed from left and right ends
5. Average read quality after trimming

These statistics help in assessing the overall quality of the sequencing run and the effect of the trimming process.

## Advanced Features

1. N-ification: Convert low-quality bases to N instead of trimming them.
2. Quality score reassignment: Assign a specific quality score to N bases.
3. Strict quality filtering: Discard reads with any base below a specified quality threshold.
4. Read subsampling: Output only a specified number of high-quality reads.
5. PALINCUT: A specialized mode for detecting and removing adapters in paired-end data by identifying palindromic sequences.

## Notes for Developers

1. The code uses a mix of C89 and C99 features. Ensure your compiler supports these standards.
2. The program relies on zlib for gzip support. Make sure this library is available on your system.
3. The code includes several macros and inline functions for performance optimization.
4. Error messages and usage information are printed to stdout, not stderr.
5. The code uses a single-pass algorithm for most operations, optimizing for speed and memory usage.

## Conclusion

SHI7 is a versatile and efficient tool for short-read quality control and trimming, particularly suited for microbiome studies. Its range of features and optimizations make it a powerful choice for processing high-throughput sequencing data. Users should carefully consider their experimental design and sequencing technology when choosing trimming parameters to achieve optimal results.
