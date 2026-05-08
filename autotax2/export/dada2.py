"""DADA2 export helpers."""

from __future__ import annotations


DADA2_TOGENUS_FORMAT = "dada2-togenus-trainset"
DADA2_SPECIES_FORMAT = "dada2-assignspecies"


def format_dada2_species_header(sequence_id: str, genus_species: str) -> str:
    """Format a DADA2 assignSpecies FASTA header."""
    return f"{sequence_id} {genus_species}"


def format_dada2_species_row(sequence_id: str, genus: str, species: str) -> str:
    """Format a DADA2 assignSpecies row."""
    return f"{sequence_id}\t{genus} {species}"
