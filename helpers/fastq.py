def read_fastq_record(file) -> tuple[str, str, str, str] | None:
    """
    Reads one record from a FASTQ file.

    Returns four lines of a record, or None if the end of the file is reached.
    """
    header = file.readline()

    if not header:
        return None

    sequence = file.readline()
    separator = file.readline()
    quality = file.readline()

    return header, sequence, separator, quality


def calculate_gc(sequence: str) -> float:
    """
    Calculates the percentage of G and C nucleotides in a sequence.
    """
    sequence = sequence.upper()
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) * 100


def calculate_average_quality(quality: str) -> float:
    """
    Calculates the average read quality using the Phred+33 scale.
    """
    if not quality:
        return 0.0

    scores = [ord(symbol) - 33 for symbol in quality.rstrip("\n")]
    if not scores:
        return 0.0

    return sum(scores) / len(scores)


def check_gc(sequence: str, gc_bounds: tuple[float, float]) -> bool:
    """
    Checks whether the GC content falls within the specified range.
    """
    gc_content = calculate_gc(sequence)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def check_length(sequence: str, length_bounds: tuple[int, int]) -> bool:
    """
    Checks whether the sequence length falls within the specified range.
    """
    length = len(sequence)
    return length_bounds[0] <= length <= length_bounds[1]


def check_quality(quality: str, quality_threshold: float) -> bool:
    """
     Checks whether the average quality meets the specified threshold.
    """
    return calculate_average_quality(quality) >= quality_threshold


def write_fastq_record(file, record: tuple[str, str, str, str]) -> None:
    """
     Writes one FASTQ record to a file.
    """
    file.writelines(record)
