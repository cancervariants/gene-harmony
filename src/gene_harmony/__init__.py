"""Gene Harmony: a tool for harmonizing ambiguous gene symbols."""

from gene_harmony.resolver import GeneHarmony

from importlib.metadata import version as _version

__all__ = ["GeneHarmony"]
__version__ = _version("gene-harmony")