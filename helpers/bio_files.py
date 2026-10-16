def read_fasta_records(input_fasta: str) -> list[tuple[str, str]]:
    """
    Reads a FASTA file and combines the lines of each sequence.
    """
    records = []
    header = None
    sequence_parts = []

    with open(input_fasta, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_parts)))

                header = line
                sequence_parts = []
            elif header is not None:
                sequence_parts.append(line)

    if header is not None:
        records.append((header, "".join(sequence_parts)))

    return records


def write_fasta_records(
    records: list[tuple[str, str]], output_fasta: str
) -> None:
    """
    Writes FASTA records with each sequence on a single line.
    """
    with open(output_fasta, "w") as file:
        for header, sequence in records:
            file.write(header + "\n")
            file.write(sequence + "\n")


def get_blast_descriptions(input_file: str) -> list[str]:
    """
     Extracts the first hit description from each BLAST results block.
    """
    descriptions = []
    in_results = False
    description_end = None

    with open(input_file, "r") as file:
        for line in file:
            if "Sequences producing significant alignments:" in line:
                in_results = True
                description_end = None
                continue

            if not in_results:
                continue

            if not line.strip():
                in_results = False
                continue

            if description_end is None:
                if "Scientific" in line:
                    description_end = line.find("Scientific")
                continue

            if line.lstrip().startswith("Description"):
                continue

            description = line[:description_end].strip()

            if description:
                descriptions.append(description)
                in_results = False

    return descriptions
