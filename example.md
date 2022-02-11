# Imports and Data Loading
Import pandas for data manipulation, plotly for plotting, and molplot for visualising structures!

```python
import pandas as pd
import plotly.express as px

import molplotly
```

Let's load the ESOL dataset from [ESOL: Estimating Aqueous Solubility Directly from Molecular Structure](https://doi.org/10.1021/ci034243x) (included as `example.csv` in the repo)

```python
# df_esol = pd.read_csv('example.csv')
df_esol = pd.read_csv(
    'https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv')
df_esol['y_pred'] = df_esol['ESOL predicted log solubility in mols per litre']
df_esol['y_true'] = df_esol['measured log solubility in mols per litre']
```

# Simple Examples
Let's make a scatter plot comparing the measured vs predicted solubilities using [`plotly`](https://plotly.com/python/)

```python
df_esol['delY'] = df_esol["y_pred"] - df_esol["y_true"]
fig_scatter = px.scatter(df_esol,
                         x="y_true",
                         y="y_pred",
                         color='delY',
                         title='ESOL Regression (default plotly)',
                         labels={'y_pred': 'Predicted Solubility',
                                 'y_true': 'Measured Solubility',
                                 'delY': 'ΔY'},
                         width=1200,
                         height=800)

# This adds a dashed line for what a perfect model _should_ predict
y = df_esol["y_true"].values
fig_scatter.add_shape(
    type="line", line=dict(dash='dash'),
    x0=y.min(), y0=y.min(),
    x1=y.max(), y1=y.max()
)

fig_scatter.show()
```
![default plotly](https://raw.githubusercontent.com/wjm41/molplotly/main/images/default.png)

now all we have to do is `add_molecules`!

```python
fig_scatter.update_layout(title='ESOL Regression (with add_molecules!)')

app_scatter = molplotly.add_molecules(fig=fig_scatter,
                                      df=df_esol,
                                      smiles_col='smiles',
                                      title_col='Compound ID'
                                      )

# change the arguments here to run the dash app on an external server and/or change the size of the app!
app_scatter.run_server(mode='inline', port=8001, height=1000)
```
![simple molplotly](https://raw.githubusercontent.com/wjm41/molplotly/main/images/example.gif)

Cool right? Let's explore some more options:

Apart from showing the $(x,y)$ coordinates (you can turn them off using `show_coords=False`), we can add extra values to show up in the mouse tooltip by specifying `caption_cols` - the values in these columns of `df_esol` are also shown in the hover box.

We can also apply some function transformations to the captions via `caption_transform` - in this example, rounding all our numbers to 2 decimal places.

```python
fig_scatter.update_layout(
    title='ESOL Regression (with add_molecules & extra captions)')

app_scatter_with_captions = molplotly.add_molecules(fig=fig_scatter,
                                                    df=df_esol,
                                                    smiles_col='smiles',
                                                    title_col='Compound ID',
                                                    caption_cols=['Molecular Weight', 'Number of Rings'],
                                                    caption_transform={'Predicted Solubility': lambda x: f"{x:.2f}",
                                                                       'Measured Solubility': lambda x: f"{x:.2f}",
                                                                       'Molecular Weight': lambda x: f"{x:.2f}"
                                                                       },
                                                    show_coords=True)

app_scatter_with_captions.run_server(mode='inline', port=8002, height=1000)
```
![with caption](https://raw.githubusercontent.com/wjm41/molplotly/main/images/caption.gif)

What about adding colors? Here I've made an arbitrary random split of the dataset into `train` and `test`. When plotting, this leads to two separate plotly "curves" so the condition determining the color of the points needs to be passed in to the `add_molecules` function in order for the correct SMILES to be selected for visualisation - this is done via `color_col`. Notice that the `title` for the molecules in the hover box have the same color as the data point! 

For fun I also used the `size` argument in the scatter plot to change the size of the markers in proportion to the molecular weight.

(notice I've been choosing different `port` numbers in all my plots, this is so that they don't interfere with each other!)

```python
from sklearn.model_selection import train_test_split

train_inds, test_inds = train_test_split(df_esol.index)
df_esol['dataset'] = [
    'Train' if x in train_inds else 'Test' for x in df_esol.index]

fig_train_test = px.scatter(df_esol,
                            x="y_true",
                            y="y_pred",
                            size='Molecular Weight',
                            color='dataset',
                            title='ESOL Regression (colored by random train/test split)',
                            labels={'y_pred': 'Predicted Solubility',
                                    'y_true': 'Measured Solubility'},
                            width=1200,
                            height=800)
# fig.show()
app_train_test = molplotly.add_molecules(fig=fig_train_test,
                                         df=df_esol,
                                         smiles_col='smiles',
                                         title_col='Compound ID',
                                         color_col='dataset')

app_train_test.run_server(mode='inline', port=8003, height=1000)
```
![with color](https://raw.githubusercontent.com/wjm41/molplotly/main/images/color.gif)

# More complex examples

Let's go beyond scatter plots and explore a few other graphs that might be relevant for cheminformatics, hopefully letting you see how `molplotly` could be useful for you when looking through (messy) data!

### Strip plots

Strip plots are useful for visualising how the same property is distributed between data from different groups. Here I plot how the measured solubility changes with the number of rings on a molecule (it goes down, surprising I know).

Violin plots can also useful for this purpose but it's not compatible with `plotly` (see section ["violin plots"](#violin)) 

```python
fig_strip = px.strip(df_esol.sort_values('Number of Rings'), # sorting so that the colorbar is sorted!
                     x='Number of Rings',
                     y='y_true',
                     color='Number of Rings',
                     labels={'y_true': 'Measured Solubility'},
                     width=1000,
                     height=800)

app_strip = molplotly.add_molecules(fig=fig_strip,
                          df=df_esol,
                          smiles_col='smiles',
                          title_col='Compound ID',
                          color_col='Number of Rings',
                          caption_transform={'Measured Solubility': lambda x: f"{x:.2f}"},
                          wrap=True,
                          wraplen=25,
                          width=150,
                          show_coords=True)

app_strip.run_server(mode='inline', port=8004, height=850)
```
![strip plot](https://raw.githubusercontent.com/wjm41/molplotly/main/images/strip.gif)

### Scatter Matrices

For visualising the relationship between multiple variables at once, use a matrix of scatter plots!

Here I've increased the width of the hover box using the `width` parameter because the caption titles were getting long; also I've used `show_coords=False` because $(x, y)$ coordinates for non-trivial scatter plots become messy.

```python
features = ['Number of H-Bond Donors',
            'Number of Rings',
            'Number of Rotatable Bonds',
            'Polar Surface Area']
fig_matrix = px.scatter_matrix(df_esol,
                               dimensions=features,
                               width=1200,
                               height=800,
                               title='Scatter matrix of molecular properties')

app_matrix = molplotly.add_molecules(fig=fig_matrix,
                                     df=df_esol,
                                     smiles_col='smiles',
                                     title_col='Compound ID',
                                     caption_cols=features,
                                     width=200,
                                     show_coords=False)

# Only show informative lower triangle
fig_matrix.update_traces(diagonal_visible=False, showupperhalf=False)
app_matrix.run_server(mode='inline', port=8005, height=1000)
```
![scatter matrix](https://raw.githubusercontent.com/wjm41/molplotly/main/images/matrix.png)

### Visualising MorganFP PCA components

A common way to visualise a molecular dataset is to calculate the morgan fingerprints of the molecules and visualise them in a 2D embedding (eg PCA/t-SNE). In this example I'm going to plot the 2 largest PCA components for ESOL and inspect the data.

Let's calculate the PCA components first!

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.decomposition import PCA


def smi_to_fp(smi):
    fp = AllChem.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles(smi), 2, nBits=1024)
    arr = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

esol_fps = np.array([smi_to_fp(smi) for smi in df_esol['smiles']])
pca = PCA(n_components=2)
components = pca.fit_transform(esol_fps.reshape(-1, 1024))
df_esol['PCA-1'] = components[:, 0]
df_esol['PCA-2'] = components[:, 1]
```

and now let's look at them!

with `molplotly`, it's super easy to see which molecules are where - steroid molecules at the top, alcohols in the bottom left, chlorinated aromatic compounds in the bottom right.

```python
fig_pca = px.scatter(df_esol,
                     x="PCA-1",
                     y="PCA-2",
                     color='y_true',
                     title='ESOL PCA of morgan fingerprints',
                     labels={'y_true': 'Measured Solubility'},
                     width=1200,
                     height=800)

app_pca = molplotly.add_molecules(fig=fig_pca,
                                  df=df_esol,
                                  smiles_col='smiles',
                                  title_col='Compound ID',
                                  caption_cols=['y_true'],
                                  color_col='y_true',
                                  show_coords=False)

app_pca.run_server(mode='inline', port=8006, height=850)
```
![PCA plot](https://raw.githubusercontent.com/wjm41/molplotly/main/images/pca.gif)

### Clustering

Let's do some clustering of the ESOL molecules, borrowing useful functions from Pat Walters' excellent blog post on [clustering](http://practicalcheminformatics.blogspot.com/2021/11/picking-highest-scoring-molecules-from.html).

```python
from rdkit.ML.Cluster import Butina

def smi2fp(smi):
    fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2)
    return fp


def taylor_butina_clustering(fp_list, cutoff=0.35):
    dists = []
    nfps = len(fp_list)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dists.extend([1-x for x in sims])
    mol_clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return mol_clusters


cluster_res = taylor_butina_clustering(
    [smi2fp(smi) for smi in df_esol['smiles']])
cluster_id_list = np.zeros(len(df_esol), dtype=int)
for cluster_num, cluster in enumerate(cluster_res):
    for member in cluster:
        cluster_id_list[member] = cluster_num
df_esol['cluster'] = cluster_id_list
```

Now let's make a strip plot of the top-10 clusters, see what they look like and how soluable they are!

```python
df_cluster = df_esol.query('cluster < 10').copy().reset_index()
# sorting is needed to make the legend appear in order!
df_cluster = df_cluster.sort_values('cluster')

fig_cluster = px.strip(df_cluster,
                      y='y_true',
                      color='cluster',
                      labels={'y_true': 'Measured Solubility'},
                      width=1000,
                      height=800)

app_cluster = molplotly.add_molecules(fig=fig_cluster,
                           df=df_cluster,
                           smiles_col='smiles',
                           title_col='Compound ID',
                           color_col='cluster'
                           )

app_cluster.run_server(mode='inline', port=8007, height=850)
```
![Clustered strip plot](https://raw.githubusercontent.com/wjm41/molplotly/main/images/cluster.gif)

# Incompatible `plotly` functionality with molplotly

`Plotly` is a graphing library that does far more than just scatter plots - it has lots of cool functionalities that unfortunately clash with how `molplotly` implements the hover box (for now at least). Here are some examples of known incompatibilities, which are still very useful data visualisations in vanilla `plotly`!

### Marginals on scatter plots 

I like having marginals on the sides by default because the data density in a dataset can often vary a lot. Anything to do with histogram/violin plots don't work yet with `molplotly`.

```python
fig_marginal = px.scatter(df_esol,
                 x="y_true",
                 y="y_pred",
                 title='ESOL Regression (with histogram marginals)',
                 labels={'y_pred': 'Predicted Solubility',
                         'y_true': 'Measured Solubility'},
                 marginal_x='violin',
                 marginal_y='histogram',
                 width=1200,
                 height=800)
fig_marginal.show()
```
![Marginal plot](https://raw.githubusercontent.com/wjm41/molplotly/main/images/marginal.png)

### Violin plots

The aesthetic of violin plots are nice, especially when there's a lot of datapoints but if there's not much data (often the case in drug discovery!) then those nice smooth KDE curves can be misleading so I usually prefer strip plots. `plotly` has cool mouseover data on violin plots which are incompatible with `molplotly` but at least if there's enough data that I prefer using a violin plot, it's probably too memory consuming to run a strip plot with `molplotly` anyway!

<a name="violin"></a>

```python
fig_violin = px.violin(df_esol,
                       y="y_true",
                       title='ESOL violin plot of measured solubility',
                       labels={'y_true': 'Measured Solubility'},
                       box=True,
                       points='all',
                       width=1200,
                       height=800)
fig_violin.show()
```
![Violin plot](https://raw.githubusercontent.com/wjm41/molplotly/main/images/violin.png)
