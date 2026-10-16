# Bioinformatics utilities

A collection of functions for working with biological sequences and BLAST results.

## Project Structure

- `bio_utils.py` - main functions.
- `helpers/fastq.py` - helper functions for FASTQ processing.
- `helpers/bio_files.py` - helper functions for working with FASTA files and BLAST results.

## FASTQ Filtering

The `filter_fastq` function filters reads by GC content, sequence length, and average quality using the Phred+33 scale.

A read is saved only if it satisfies all specified conditions. The output file is written to the `filtered/` directory.

Example:

```python
from bio_utils import filter_fastq

filter_fastq(
    input_fastq="reads.fastq",
    output_fastq="filtered_reads.fastq",
    gc_bounds=(30, 70),
    length_bounds=(50, 300),
    quality_threshold=20,
)
```
## FASTA Conversion

The `convert_multiline_fasta_to_oneline` function combines multiline FASTA sequences into single lines.

If the output file is not specified, a new file with the `_oneline` suffix is created.

Example:

```python
from bio_utils import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline(
    "sequences.fasta",
    "sequences_oneline.fasta",
)
```
## BLAST Processing

The `parse_blast_output` function extracts the first hit description from each BLAST results block, sorts the descriptions alphabetically, and writes them to an output file.

Example:

```python
from bio_utils import parse_blast_output

parse_blast_output(
    "blast_report.txt",
    "blast_descriptions.txt",
)
```
