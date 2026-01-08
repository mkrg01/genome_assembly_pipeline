from snakemake.script import snakemake
import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t", header=None, names=["contig", "start", "end", "repeat_class", "score", "strand"])

if snakemake.wildcards.repeat_class == "LTR":
    df = df[df["repeat_class"].str.startswith("LTR")]
    df.loc[:,"repeat_class"] = "LTR"
elif snakemake.wildcards.repeat_class == "Copia":
    df = df[df["repeat_class"].str.startswith("LTR/Copia")]
    df.loc[:,"repeat_class"] = "Copia"
elif snakemake.wildcards.repeat_class == "Gypsy":
    df = df[df["repeat_class"].str.startswith("LTR/Gypsy")]
    df.loc[:,"repeat_class"] = "Gypsy"
elif snakemake.wildcards.repeat_class == "LINE":
    df = df[df["repeat_class"].str.startswith("LINE")]
    df.loc[:,"repeat_class"] = "LINE"
elif snakemake.wildcards.repeat_class == "SINE":
    df = df[df["repeat_class"].str.startswith("SINE")]
    df.loc[:,"repeat_class"] = "SINE"
elif snakemake.wildcards.repeat_class == "DNA_transposon":
    df = df[df["repeat_class"].str.startswith("DNA")]
    df.loc[:,"repeat_class"] = "DNA_transposon"
elif snakemake.wildcards.repeat_class == "satellite":
    df = df[df["repeat_class"].str.startswith("Satellite")]
    df.loc[:,"repeat_class"] = "satellite"
else:
    raise ValueError(f"Unsupported repeat class: {snakemake.wildcards.repeat_class}")

df.to_csv(snakemake.output[0], sep="\t", header=False, index=False)
