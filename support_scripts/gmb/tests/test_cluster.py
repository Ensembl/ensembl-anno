import pandas as pd
import pyranges as pr

df = pd.DataFrame({
    'Chromosome': ['1', '1'],
    'Start': [100, 100],
    'End': [200, 200],
    'Strand': ['+', '-']
})

gr = pr.PyRanges(df)
print("Original:")
print(gr)

print("\ncluster(strand=None) [Default]:")
try:
    c = gr.cluster(strand=None)
    print(c)
    print(f"Clusters: {c.Cluster.unique()}")
except Exception as e:
    print(e)

print("\ncluster(strand=True):")
try:
    c = gr.cluster(strand=True)
    print(c)
    print(f"Clusters: {c.Cluster.unique()}")
except Exception as e:
        print(e)

print("\ncluster(strand=False):")
try:
    c = gr.cluster(strand=False)
    print(c)
    print(f"Clusters: {c.Cluster.unique()}")
except Exception as e:
    print(e)
