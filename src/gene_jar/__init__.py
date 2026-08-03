"""Gene Harmony: a tool for harmonizing ambiguous gene symbols."""

from gene_jar.resolver import GeneJar

from importlib.metadata import version as _version

__all__ = ["GeneJar"]
__version__ = _version("gene-jar")