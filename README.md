# Bioinformatics utilities

Набор функций для работы с биологическими последовательностями и результатами BLAST.

## Структура проекта

- `bio_utils.py` - основные функции.
- `helpers/fastq.py` - вспомогательные функции для обработки FASTQ.
- `helpers/bio_files.py` - вспомогательные функции для работы с FASTA и BLAST.

## Фильтрация FASTQ

Функция `filter_fastq` отбирает риды по GC-составу, длине последовательности и среднему качеству по шкале Phred+33.

Рид сохраняется, если удовлетворяет всем заданным условиям. Результат записывается в каталог `filtered/`.

Пример использования:

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
## Преобразование FASTA

Функция `convert_multiline_fasta_to_oneline` объединяет строки каждой последовательности FASTA в одну строку.

Если выходной файл не указан, создаётся файл с суффиксом `_oneline`.

Пример использования:

```python
from bio_utils import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline(
    "sequences.fasta",
    "sequences_oneline.fasta",
)
```
## Обработка BLAST

Функция `parse_blast_output` извлекает первое описание совпадения из каждого блока результатов BLAST, сортирует описания по алфавиту и записывает их в выходной файл.

Пример использования:

```python
from bio_utils import parse_blast_output

parse_blast_output(
    "blast_report.txt",
    "blast_descriptions.txt",
)
```
