import pandas as pd
import pyranges as pr

df1 = pr.read_gff3("z_tritici/helixer_remapped.gff3").df
df1 = df1[df1["Feature"] == "exon"]
c1 = pr.PyRanges(df1).cluster(slack=0).df
print("Helixer clusters:", c1["Cluster"].nunique())

df2 = pr.read_gtf("z_tritici/stringtie_geneset.gtf").df
df2 = df2[df2["Feature"] == "exon"]
c2 = pr.PyRanges(df2).cluster(slack=0).df
print("StringTie clusters:", c2["Cluster"].nunique())

df3 = pd.concat([df1, df2], ignore_index=True)
c3 = pr.PyRanges(df3).cluster(slack=0).df
print("Combined clusters:", c3["Cluster"].nunique())
