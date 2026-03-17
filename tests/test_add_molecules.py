import molplotly
import pandas as pd
import plotly.express as px
from dash import Dash

from . import ROOT

df_esol = pd.read_csv(f"{ROOT}/examples/example.csv")
df_esol["y_pred"] = df_esol["ESOL predicted log solubility in mols per litre"]
df_esol["y_true"] = df_esol["measured log solubility in mols per litre"]
df_esol["delY"] = df_esol["y_pred"] - df_esol["y_true"]


def _make_scatter(**kwargs):
    return px.scatter(
        df_esol,
        x="y_true",
        y="y_pred",
        labels={
            "y_pred": "Predicted Solubility",
            "y_true": "Measured Solubility",
        },
        **kwargs,
    )


def test_add_molecules_basic():
    """Basic scatter plot returns a Dash app."""
    fig = _make_scatter()
    app = molplotly.add_molecules(
        fig=fig, df=df_esol, smiles_col="smiles", title_col="Compound ID"
    )
    assert isinstance(app, Dash)


def test_add_molecules_with_color():
    """Scatter with discrete color column."""
    df_esol["solubility_class"] = pd.cut(
        df_esol["y_true"], bins=3, labels=["low", "mid", "high"]
    )
    fig = px.scatter(
        df_esol,
        x="y_true",
        y="y_pred",
        color="solubility_class",
        labels={
            "y_pred": "Predicted Solubility",
            "y_true": "Measured Solubility",
        },
    )
    app = molplotly.add_molecules(
        fig=fig,
        df=df_esol,
        smiles_col="smiles",
        color_col="solubility_class",
    )
    assert isinstance(app, Dash)


def test_add_molecules_with_symbol():
    """Scatter with symbol column."""
    df_esol["solubility_class"] = pd.cut(
        df_esol["y_true"], bins=2, labels=["low", "high"]
    )
    fig = px.scatter(
        df_esol,
        x="y_true",
        y="y_pred",
        symbol="solubility_class",
        labels={
            "y_pred": "Predicted Solubility",
            "y_true": "Measured Solubility",
        },
    )
    app = molplotly.add_molecules(
        fig=fig,
        df=df_esol,
        smiles_col="smiles",
        symbol_col="solubility_class",
    )
    assert isinstance(app, Dash)


def test_add_molecules_with_captions():
    """Captions and caption transforms."""
    fig = _make_scatter()
    app = molplotly.add_molecules(
        fig=fig,
        df=df_esol,
        smiles_col="smiles",
        caption_cols=["Compound ID", "y_pred"],
        caption_transform={"y_pred": lambda x: f"{x:.2f}"},
    )
    assert isinstance(app, Dash)


def test_add_molecules_no_img():
    """Show img disabled."""
    fig = _make_scatter()
    app = molplotly.add_molecules(
        fig=fig, df=df_esol, smiles_col="smiles", show_img=False
    )
    assert isinstance(app, Dash)


def test_add_molecules_no_coords():
    """Show coords disabled."""
    fig = _make_scatter()
    app = molplotly.add_molecules(
        fig=fig, df=df_esol, smiles_col="smiles", show_coords=False
    )
    assert isinstance(app, Dash)


def test_add_molecules_multiple_smiles():
    """Multiple SMILES columns produce a dropdown."""
    df_esol["smiles2"] = df_esol["smiles"]
    fig = _make_scatter()
    app = molplotly.add_molecules(
        fig=fig, df=df_esol, smiles_col=["smiles", "smiles2"]
    )
    assert isinstance(app, Dash)


def test_add_molecules_custom_svg_size():
    """Custom SVG dimensions."""
    fig = _make_scatter()
    app = molplotly.add_molecules(
        fig=fig,
        df=df_esol,
        smiles_col="smiles",
        svg_height=300,
        svg_width=400,
    )
    assert isinstance(app, Dash)


def test_add_molecules_version():
    """Package exposes a version string."""
    assert hasattr(molplotly, "__version__")
    assert isinstance(molplotly.__version__, str)
