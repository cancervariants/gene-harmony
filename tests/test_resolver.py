"""Tests for the GeneJar resolver."""

import pandas as pd

from gene_jar import GeneJar


def test_resolve_finds_a1b_in_gene_symbol():
    """Test that resolving A1BG returns a gene_symbol containing A1B."""
    primary_df = pd.DataFrame(
        {
            "gene_symbol": ["A1B", "TP53", "KRAS"],
            "primary_gene_symbol": ["A1BG", "TP53", "KRAS"],
        }
    )

    alias_df = pd.DataFrame(
        {
            "alias_symbol": [],
            "primary_gene_symbol": [],
        }
    )

    gj = GeneJar(
        primary_df=primary_df,
        ortholog_df=alias_df.copy(),
        flj_clone_df=alias_df.copy(),
        phenotype_df=alias_df.copy(),
        hgnc_gene_group_df=alias_df.copy(),
        gene_id_df=alias_df.copy(),
        mgi_withdrawn_df=alias_df.copy(),
        ncbi_related_gene_df=alias_df.copy(),
        ncbi_gene_interaction_df=alias_df.copy(),
        ncbi_gene_neighbor_df=alias_df.copy(),
        placeholder_df=alias_df.copy(),
        previous_df=alias_df.copy(),
        protein_mass_df=alias_df.copy(),
        alternate_abbreviation_df=alias_df.copy(),
    )

    result = gj.resolve(symbol="A1BG")

    assert "A1B" in result["gene_symbol"].values
