"""Default rank identity thresholds for autotax2."""

from __future__ import annotations

DEFAULT_SPECIES_ID = 0.972
DEFAULT_GENUS_ID = 0.901
DEFAULT_FAMILY_ID = 0.801
DEFAULT_ORDER_ID = 0.729
DEFAULT_CLASS_ID = 0.722
DEFAULT_PHYLUM_ID = 0.696

RANK_IDENTITY_DEFAULTS = {
    "species": DEFAULT_SPECIES_ID,
    "genus": DEFAULT_GENUS_ID,
    "family": DEFAULT_FAMILY_ID,
    "order": DEFAULT_ORDER_ID,
    "class": DEFAULT_CLASS_ID,
    "phylum": DEFAULT_PHYLUM_ID,
}
