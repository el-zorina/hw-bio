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

def normalize_bounds(bounds, default) -> tuple[float, float] | tuple[int, int]:
    """
    Приводит границы фильтра к паре значений.

    Одно число считается верхней границей интервала.
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
    Фильтрует FASTQ-файл по GC-составу, длине и среднему качеству.

    Подходящие риды записываются в папку filtered.
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


if __name__ == "__main__":
    print("FASTQ utilities are ready.")

def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str | None = None,
) -> None:
    """
    Преобразует многострочный FASTA в формат с одной строкой на последовательность.

    Если выходной файл не указан, создаёт файл с суффиксом _oneline.
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
    Сохраняет первые описания совпадений из BLAST-отчёта по алфавиту.
    """
    descriptions = get_blast_descriptions(input_file)
    descriptions.sort()

    with open(output_file, "w") as file:
        for description in descriptions:
            file.write(description + "\n")
