import pandas as pd

class GeneHarmony:
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
        """
        Initialize with dataframes.
        primary_df and ortholog_df must have columns: 'gene_symbol', 'primary_gene_symbol'
        """
        self.dfs = {
            "Primary": primary_df,
            "Ortholog Symbol": ortholog_df,
            "Clone Symbol":flj_clone_df, 
            "Prefix Condition Symbol":disease_df, 
            "Prefix Gene Group Symbol":hgnc_gene_group_df, 
            "Gene Identifier Symbol":gene_id_df,
            "Withdrawn MGI Symbol":mgi_withdrawn_df, 
            "Gene Group Symbol":ncbi_gene_group_df, 
            "Gene Interaction Symbol":ncbi_gene_interaction_df, 
            "Gene Neighbor Symbol":ncbi_gene_neighbor_df, 
            "Placeholder Symbol":placeholder_df, 
            "Previous Symbol":previous_df, 
            "Protein Mass Symbol":protein_mass_df
        }
        default_cols = ("alias_symbol", "primary_gene_symbol")
        self.column_map = {
            "Primary": ("gene_symbol", "primary_gene_symbol"),
        }

        for key in self.dfs:
            self.column_map.setdefault(key, default_cols)

    def resolve(self, symbol: str, source: str = "Primary", match_type: str = "identical") -> pd.DataFrame:
        """
        Resolve a gene symbol from a specific source.

        :param symbol: Gene symbol to search for.
        :param source: Source dataframe key.
        :param match_type: "identical" for exact matches, "partial" for substring matches.
        :return: Filtered DataFrame.
        """
        if source not in self.dfs:
            raise ValueError(f"Unknown source '{source}'. Available: {sorted(self.dfs)}")
        
        df = self.dfs[source]
        col_main, col_primary = self.column_map[source]
        target = symbol.upper()
        
        if match_type == "identical":
            mask = (
                df[col_main].astype(str).str.upper().eq(target)
                | df[col_primary].astype(str).str.upper().eq(target)
            )
        elif match_type == "partial":
            mask = (
                df[col_main].astype(str).str.upper().str.contains(target, na=False)
                | df[col_primary].astype(str).str.upper().str.contains(target, na=False)
            )
        else:
            raise ValueError("match_type must be 'identical' or 'partial'")
        
        return df.loc[mask].copy()