import json
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict


class MedianPolish:
    def fit(self, df: pd.DataFrame, maxiter=3):
        self.overall = np.nanmedian(df)
        self.row = np.zeros(df.shape[0])
        self.col = np.zeros(df.shape[1])

        df -= self.overall

        for _ in range(maxiter):
            r = df.median(axis=1)
            self.row += r
            df -= r.values[:, np.newaxis]
            c = df.median()
            self.col += c
            df -= c
        self.df = df
        return self

class PeptideDataset:
    @classmethod
    def load_sage(cls, folder: str) -> "PeptideDataset":
        self = cls.__new__(cls)

        """Load Sage search and quant results from an output folder"""
        df = (
            pd.read_csv(f"{folder}/results.sage.tsv", sep="\t")
            .loc[
                lambda x: (x.spectrum_q <= 0.01)
                & (x.peptide_q <= 0.01)
                & (x.protein_q <= 0.01)
                & (x.label == 1)  # Remove decoy hits
                & (x.num_proteins == 1)  # Unique only
            ]
            .drop_duplicates(subset="peptide")
            .assign(organism=lambda x: x.proteins.str.extract("_(\w+)$"))
        )

        lfq = (
            pd.read_csv(f"{folder}/lfq.tsv", sep="\t")
            .set_index("peptide")
            .loc[lambda df: df.q_value <= 0.05]
            .drop(["spectral_angle", "score", "q_value", "proteins"], axis=1)
        )
        self.df = pd.merge(
            df[["peptide", "proteins", "organism"]],
            lfq,
            left_on="peptide",
            right_on="peptide",
            how="left",
        )
        return self

    @classmethod
    def load_fragpipe(cls, folder: str) -> "PeptideDataset":
        self = cls.__new__(cls)
        self.df = (
            pd.read_csv(f"{folder}/MSstats.csv")
            .rename({"ProteinName": "proteins", "PeptideSequence": "peptide"}, axis=1)
            .assign(organism=lambda x: x.proteins.str.extract("_(\w+)$"))
            .assign(variable=lambda df: df.Run + ".mzML")
            .pivot_table(
                index=["proteins", "peptide", "organism"],
                columns="variable",
                values="Intensity",
            )
            .reset_index()
        )
        return self

    @classmethod
    def load_maxquant(cls, folder: str) -> "PeptideDataset":
        self = cls.__new__(cls)
        self.df = (
            pd.read_csv(f"{folder}/evidence.txt", sep="\t")
            .loc[
                lambda df: (df.Proteins.str.count(";") == 0)
                & (df.Reverse.astype(str) != "+")
            ]
            .rename({"Sequence": "peptide", "Proteins": "proteins"}, axis=1)
            .assign(organism=lambda x: x.proteins.str.extract("_(\w+)$"))
            .assign(variable=lambda df: df["Raw file"] + ".mzML")
            .pivot_table(
                index=["peptide", "proteins", "organism"],
                columns="variable",
                values="Intensity",
            )
            .reset_index()
        )
        return self

    def prep_msstats(self, groups) -> pd.DataFrame:
        return (
            self.df.rename(
                {"peptide": "PeptideSequence", "proteins": "ProteinName"}, axis=1
            )
            .drop(["organism","charge"], axis=1)
            .melt(id_vars=["PeptideSequence", "ProteinName"])
            .assign(
                PrecursorCharge=2,
                FragmentIon=np.nan,
                ProductCharge=np.nan,
                IsotopeLabelType="L",
                BioReplicate=1,
            )
            .assign(Condition=lambda df: df.variable.apply(lambda x: groups[x]))
            .rename({"variable": "Run", "value": "Intensity"}, axis=1)
        )

    def peptide_median_polish(self) -> pd.DataFrame:
        def summarize(df):
            m = MedianPolish()
            m.fit(df)
            return m.col + m.overall

        return (
            self.df.set_index(["proteins", "peptide", "organism", "charge"])
            .dropna(how="all")
            .applymap(lambda x: np.log2(x + 1))
            .groupby(["proteins", "organism"])
            .apply(summarize)
        )




def groups() -> Dict[str, str]:
    r = re.compile(r"(B|T4)\.naive")
    metadata = json.loads(open(f"results.json").read())
    files = [c.split("/")[-1] for c in metadata["mzml_paths"]]
    groups = {c: r.search(c).groups()[0] for c in files}
    return groups


s = PeptideDataset.load_sage(".")
s.prep_msstats(groups()).to_csv("MSstats-Sage.csv")
s.peptide_median_polish().to_csv("lfq.proteins.tsv", sep='\t')

# PeptideDataset.load_maxquant(".").prep_msstats(groups()).to_csv("MSstats-MaxQuant.csv")
# PeptideDataset.load_fragpipe(".").prep_msstats(groups()).to_csv("MSstats-FragPipe.csv")
