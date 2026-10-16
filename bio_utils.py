from pathlib import Path

from helpers.fastq import (
    read_fastq_record,
    check_gc,
    check_length,
    check_quality,
    write_fastq_record,
)
from helpers.bio_files import (
    read_fasta_records,
    write_fasta_records,
    get_blast_descriptions,
)

def normalize_bounds(
    bounds: tuple[float, float] | tuple[int, int] | float | int | None,
    default: tuple[float, float] | tuple[int, int],
) -> tuple[float, float] | tuple[int, int]:
    """
    Converts filter bounds into a pair of values.

    Args:
        bounds: A single upper bound or a pair of lower and upper bounds.
            If None, the default bounds are used.
        default: Bounds to use when no value is provided.

    Returns:
        A pair of lower and upper filter bounds.
    """
    if bounds is None:
        return default

    if isinstance(bounds, (int, float)):
        return default[0], bounds

    return bounds

def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters a FASTQ file by GC content, sequence length, and average quality.

    Args:
        input_fastq: Path to the input FASTQ file.
        output_fastq: Name of the output file in the filtered directory.
        gc_bounds: Allowed GC content range in percent, or only its upper bound.
        length_bounds: Allowed read length range, or only its upper bound.
        quality_threshold: Minimum average read quality using Phred+33.

    Returns:
        None. Reads that pass all filters are written to the filtered directory.
    """
    gc_bounds = normalize_bounds(gc_bounds, (0, 100))
    length_bounds = normalize_bounds(length_bounds, (0, 2**32))

    output_dir = Path("filtered")
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / output_fastq

    with open(input_fastq, "r") as input_file, open(
        output_path, "w"
    ) as output_file:
        while True:
            record = read_fastq_record(input_file)

            if record is None:
                break

            header, sequence, separator, quality = record
            sequence = sequence.strip()
            quality = quality.strip()

            if (
                check_gc(sequence, gc_bounds)
                and check_length(sequence, length_bounds)
                and check_quality(quality, quality_threshold)
            ):
                write_fastq_record(
                    output_file,
                    (header, sequence + "\n", separator, quality + "\n"),
                )

def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str | None = None,
) -> None:
    """
    Converts a multiline FASTA file to a format with one line per sequence.

    Args:
        input_fasta: Path to the input FASTA file.
        output_fasta: Path to the output file. If not specified, a file with
            the _oneline suffix is created before the file extension.

    Returns:
        None. The result is written to the output FASTA file.
    """
    if output_fasta is None:
        input_path = Path(input_fasta)
        output_fasta = str(
            input_path.with_name(input_path.stem + "_oneline" + input_path.suffix)
        )

    records = read_fasta_records(input_fasta)
    write_fasta_records(records, output_fasta)

def parse_blast_output(input_file: str, output_file: str) -> None:
    """
     Extracts the first hit description from each BLAST results block and sorts
    the descriptions alphabetically.

    Args:
        input_file: Path to the BLAST text report.
        output_file: Path to the file for sorted descriptions.

    Returns:
        None. Each description is written on a separate line.
    """
    descriptions = get_blast_descriptions(input_file)
    descriptions.sort()

    with open(output_file, "w") as file:
        for description in descriptions:
            file.write(description + "\n")

if __name__ == "__main__":
    print("FASTQ utilities are ready.")