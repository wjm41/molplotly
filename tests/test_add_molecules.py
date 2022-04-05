import molplotly
import multiprocessing
import time
import pandas as pd
import plotly.express as px
from jupyter_dash import JupyterDash

from . import ROOT

df_esol = pd.read_csv(f"{ROOT}/examples/example.csv")
df_esol["y_pred"] = df_esol["ESOL predicted log solubility in mols per litre"]
df_esol["y_true"] = df_esol["measured log solubility in mols per litre"]


df_esol["delY"] = df_esol["y_pred"] - df_esol["y_true"]
fig_scatter = px.scatter(
    df_esol,
    x="y_true",
    y="y_pred",
    color="delY",
    title="ESOL Regression (default plotly)",
    labels={
        "y_pred": "Predicted Solubility",
        "y_true": "Measured Solubility",
        "delY": "ΔY",
    },
    width=1200,
    height=800,
)


def test_add_molecules():
    app_scatter = molplotly.add_molecules(
        fig=fig_scatter, df=df_esol, smiles_col="smiles", title_col="Compound ID"
    )

    assert isinstance(app_scatter, JupyterDash)
