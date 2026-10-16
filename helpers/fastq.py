def read_fastq_record(file) -> tuple[str, str, str, str] | None:
    """
    Читает одну запись из FASTQ-файла.

    Возвращает четыре строки записи или None, если файл закончился.
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
    Вычисляет процент G и C в последовательности.
    """
    sequence = sequence.upper()
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) * 100


def calculate_average_quality(quality: str) -> float:
    """
    Вычисляет среднее качество рида по шкале Phred+33.
    """
    if not quality:
        return 0.0

    scores = [ord(symbol) - 33 for symbol in quality.rstrip("\n")]
    if not scores:
        return 0.0

    return sum(scores) / len(scores)


def check_gc(sequence: str, gc_bounds: tuple[float, float]) -> bool:
    """
    Проверяет, попадает ли GC-состав в заданный интервал.
    """
    gc_content = calculate_gc(sequence)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def check_length(sequence: str, length_bounds: tuple[int, int]) -> bool:
    """
    Проверяет, попадает ли длина последовательности в заданный интервал.
    """
    length = len(sequence)
    return length_bounds[0] <= length <= length_bounds[1]


def check_quality(quality: str, quality_threshold: float) -> bool:
    """
    Проверяет, не ниже ли среднее качество заданного порога.
    """
    return calculate_average_quality(quality) >= quality_threshold


def write_fastq_record(file, record: tuple[str, str, str, str]) -> None:
    """
    Записывает одну FASTQ-запись в файл.
    """
    file.writelines(record)
