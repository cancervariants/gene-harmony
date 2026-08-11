"""Utilities for resolving gene symbols across multiple harmonized reference datasets."""

from enum import StrEnum
from dataclasses import dataclass
from typing import Any

import pandas as pd


class MatchType(StrEnum):
    """Match options"""

    IDENTICAL = "identical"
    PARTIAL = "partial"

@dataclass
class AmbiguityResult:
    """Result of a gene symbol ambiguity check."""

    ambiguous: bool
    primary_gene_symbol: list[str]
    hgnc_id: list[Any]
    ncbi_id: list[Any]
    ensg_id: list[Any]



class GeneJar:
    """Resolve gene symbols across multiple harmonized gene reference datasets.

    :param primary_df: DataFrame containing primary gene symbols. Must include
        ``gene_symbol`` and ``primary_gene_symbol``.
    :param ortholog_df: DataFrame containing ortholog symbols.
    :param flj_clone_df: DataFrame containing FLJ clone symbols.
    :param disease_df: DataFrame containing disease-related prefix symbols.
    :param hgnc_gene_group_df: DataFrame containing HGNC gene group symbols.
    :param gene_id_df: DataFrame containing gene identifier symbols.
    :param mgi_withdrawn_df: DataFrame containing withdrawn MGI symbols.
    :param ncbi_gene_group_df: DataFrame containing NCBI gene group symbols.
    :param ncbi_gene_interaction_df: DataFrame containing NCBI gene interaction symbols.
    :param ncbi_gene_neighbor_df: DataFrame containing NCBI gene neighbor symbols.
    :param placeholder_df: DataFrame containing placeholder symbols.
    :param previous_df: DataFrame containing previous gene symbols.
    :param protein_mass_df: DataFrame containing protein mass symbols.
    """

    def __init__(
        self,
        *,
        primary_df: pd.DataFrame,
        ortholog_df: pd.DataFrame,
        flj_clone_df: pd.DataFrame,
        disease_df: pd.DataFrame,
        hgnc_gene_group_df: pd.DataFrame,
        gene_id_df: pd.DataFrame,
        mgi_withdrawn_df: pd.DataFrame,
        ncbi_gene_group_df: pd.DataFrame,
        ncbi_gene_interaction_df: pd.DataFrame,
        ncbi_gene_neighbor_df: pd.DataFrame,
        placeholder_df: pd.DataFrame,
        previous_df: pd.DataFrame,
        protein_mass_df: pd.DataFrame,
    ):
        """Initialize with dataframes.
        primary_df and ortholog_df must have columns: 'gene_symbol', 'primary_gene_symbol'
        """
        self.dfs = {
            "Primary": primary_df,
            "Ortholog Symbol": ortholog_df,
            "Clone Symbol": flj_clone_df,
            "Prefix Condition Symbol": disease_df,
            "Prefix Gene Group Symbol": hgnc_gene_group_df,
            "Gene Identifier Symbol": gene_id_df,
            "Withdrawn MGI Symbol": mgi_withdrawn_df,
            "Gene Group Symbol": ncbi_gene_group_df,
            "Gene Interaction Symbol": ncbi_gene_interaction_df,
            "Gene Neighbor Symbol": ncbi_gene_neighbor_df,
            "Placeholder Symbol": placeholder_df,
            "Previous Symbol": previous_df,
            "Protein Mass Symbol": protein_mass_df,
        }
        default_cols = ("alias_symbol", "primary_gene_symbol")
        self.column_map = {
            "Primary": ("gene_symbol", "primary_gene_symbol"),
        }

        for key in self.dfs:
            self.column_map.setdefault(key, default_cols)


    def symbol_categories(self) -> list[str]:
        """Return the valid category names accepted by resolve()."""
        return sorted(self.dfs)

    def resolve(
        self,
        symbol: str,
        symbol_category: str = "Primary",
        match_type: MatchType = MatchType.IDENTICAL,
    ) -> pd.DataFrame:
        """Resolve a gene symbol from a specific symbol category.

        :param symbol: Gene symbol to search for.
        :param symbol_category: Source dataframe key.
        :param match_type: Type of matching to perform.
        :return: Filtered DataFrame.
        """
        if symbol_category not in self.dfs:
            message = (
                f"Unknown symbol category '{symbol_category}'. "
                f"Available: {sorted(self.dfs)}"
            )
            raise ValueError(message)

        df = self.dfs[symbol_category]
        col_main, col_primary = self.column_map[symbol_category]
        target = symbol.upper()

        if match_type is MatchType.IDENTICAL:
            mask = df[col_main].astype(str).str.upper().eq(target) | df[
                col_primary
            ].astype(str).str.upper().eq(target)
        elif match_type is MatchType.PARTIAL:
            mask = df[col_main].astype(str).str.upper().str.contains(
                target, na=False
            ) | df[col_primary].astype(str).str.upper().str.contains(target, na=False)

        return df.loc[mask].copy()
    def _flatten_unique(self, values: pd.Series) -> list[Any]:
        """Flatten nested identifier values and return unique non-null entries."""
        flattened = []

        for value in values.dropna():
            if isinstance(value, (set, list, tuple)):
                flattened.extend(value)
            else:
                flattened.append(value)

        return list(dict.fromkeys(flattened))

    def ambiguity_check(self, symbol: str) -> AmbiguityResult:
        """Check whether a symbol is associated with multiple primary genes."""
        result = self.resolve(
            symbol=symbol,
            match_type=MatchType.IDENTICAL,
        )

        primary_symbols = self._flatten_unique(
            result["primary_gene_symbol"]
        )

        return AmbiguityResult(
            ambiguous=len(primary_symbols) > 1,
            primary_gene_symbol=primary_symbols,
            hgnc_id=self._flatten_unique(result["HGNC_ID"]),
            ncbi_id=self._flatten_unique(result["NCBI_ID"]),
            ensg_id=self._flatten_unique(result["ENSG_ID"]),
        )