from __future__ import annotations

import base64
import textwrap
from io import BytesIO
from typing import Callable
import itertools
import re

import pandas as pd
import numpy as np
from dash import Input, Output, dcc, html, no_update
from jupyter_dash import JupyterDash
from pandas.core.groupby import DataFrameGroupBy
from plotly.graph_objects import Figure
import plotly.graph_objects as go

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D


def str2bool(v: str) -> bool:
    return v.lower() in ("yes", "true", "t", "1")


def test_groups(fig: Figure, df_grouped: DataFrameGroupBy):
    """Test if plotly figure curve names match up with pandas dataframe groups

    Args:
        fig (plotly figure): _description_
        groups (pandas groupby object): _description_

    Returns:
        _type_: Bool describing whether or not groups is the correct dataframe grouping descbrining the data in fig
    """
    str_groups = {}
    for name, group in df_grouped:
        if isinstance(name, tuple):
            str_groups[", ".join(str(x) for x in name)] = group
        else:
            str_groups[name] = group

    for data in fig.data:
        if data.name in str_groups:
            if len(data.y) == len(str_groups[data.name]):
                continue
        else:
            return False
    return True


def find_correct_column_order(cols, col_names):

    correctly_ordered_cols = []
    for col in col_names:
        if col in cols:
            correctly_ordered_cols.append(col)
    return correctly_ordered_cols


def find_grouping(
    fig: Figure, df_data: pd.DataFrame, cols: list[str]
) -> tuple[DataFrameGroupBy, dict]:

    if fig.data[0].hovertemplate is not None:
        col_names = re.findall(r"(.*?)=.*?<.*?>", fig.data[0].hovertemplate)
        col_names = [re.sub(r"(.*)>", "", col_name) for col_name in col_names]
        if set(cols).issubset(set(col_names)) is False:
            raise ValueError(
                f"symbol_col/color_col/facet_col is misspecified because the specified dataframe grouping names {cols} don't match the names in the plotly figure {col_names}.",
            )

        cols = find_correct_column_order(cols, col_names)
        df_grouped = df_data.groupby(cols)

        str_groups = {}
        for name, group in df_grouped:
            if isinstance(name, tuple):
                str_groups[", ".join(str(x) for x in name)] = group
            else:
                str_groups[name] = group

        curve_dict = {}
        for index, data in enumerate(fig.data):
            curve_name = re.findall(r".*?=(.*?)<.*?>", data.hovertemplate)
            curve_name = ", ".join(str(x) for x in curve_name)
            if "%{x}" in curve_name:
                unique_x_values = np.unique(data.x)
                if len(unique_x_values) == 1 and len(data.x) != 1:
                    curve_name = curve_name.replace("%{x}", str(unique_x_values[0]))
                else:
                    curve_name = curve_name.replace(", %{x}", "")
            if "%{y}" in curve_name:
                unique_y_values = np.unique(data.y)
                if len(unique_y_values) == 1 and len(data.y) != 1:
                    curve_name = curve_name.replace("%{y}", str(unique_y_values[0]))
                else:
                    curve_name = curve_name.replace(", %{y}", "")
            if "%{marker.size}" in curve_name:
                unique_size_values = np.unique(data.marker.size)
                if len(unique_size_values) == 1 and len(data.marker.size) != 1:
                    curve_name = curve_name.replace(
                        "%{marker.size}", str(unique_size_values[0])
                    )
                else:
                    curve_name = curve_name.replace(", %{marker.size}", "")
            curve_dict[index] = str_groups[curve_name]

        return df_grouped, curve_dict

    else:
        grouping_found = False
        for combo in itertools.permutations(cols):
            df_grouped_tmp = df_data.groupby(list(combo))
            if test_groups(fig, df_grouped_tmp):
                df_grouped = df_grouped_tmp
                str_groups = {}
                for name, group in df_grouped:
                    if isinstance(name, tuple):
                        str_groups[", ".join(str(x) for x in name)] = group
                    else:
                        str_groups[name] = group
                curve_dict = {
                    index: str_groups[x["name"]] for index, x in enumerate(fig.data)
                }
                grouping_found = True

                return df_grouped, curve_dict

            else:
                continue

        if not grouping_found:
            raise ValueError(
                "symbol_col/color_col/facet_col is misspecified because the dataframe grouping names don't match the names in the plotly figure."
            )


def add_molecules(
    fig: Figure,
    df: pd.DataFrame,
    smiles_col: str | list[str] = "SMILES",
    mol_col: Mol | list[Mol] = None,
    show_img: bool = True,
    svg_size: int = 200,
    svg_height: int | None = None,
    svg_width: int | None = None,
    alpha: float = 0.75,
    mol_alpha: float = 0.7,
    title_col: str = None,
    show_coords: bool = True,
    caption_cols: list[str] = None,
    caption_transform: dict[str, Callable] = {},
    color_col: str = None,
    symbol_col: str = None,
    facet_col: str = None,
    wrap: bool = True,
    wraplen: int = 20,
    width: int = 150,
    fontfamily: str = "Arial",
    fontsize: int = 12,
    reaction: bool = False,
) -> JupyterDash:
    """
    A function that takes a plotly figure and a dataframe with molecular SMILES
    and returns a dash app that dynamically generates an image of molecules in the hover box
    when hovering the mouse over datapoints.
    ...

    Attributes
    ----------
    fig : plotly.graph_objects.Figure object
        a plotly figure object containing datapoints plotted from df.
    df : pandas.DataFrame object
        a pandas dataframe that contains the data plotted in fig.
    smiles_col : str | list[str], optional
        name of the column in df containing the smiles plotted in fig (default 'SMILES').
        If provided as a list, will add a slider to choose which column is used for rendering the structures.
    mol_col : Mol | list[Mol], optional
        name of the column in df containing RDKit Mol objects of the molecules plotted in fig (default None).
        If not None, the structures will be drawn using the coordinates of the Mol objects.
        If provided as a list, will add a slider to choose which column is used for rendering the structures.
    show_img : bool, optional
        whether or not to generate the molecule image in the dash app (default True).
    svg_size: int, optional
        the size in pixels of the height and width of the drawing. Is overridden by svg_height or svg_width (default 200)
    svg_height : int, optional
        the svg_height in pixels of the molecule drawing (default None).
    svg_width : int, optional
        the svg_width in pixels of the molecule drawing (default None).
    alpha : float, optional
        the transparency of the hoverbox, 0 for full transparency 1 for full opaqueness (default 0.7).
    mol_alpha : float, optional
        the transparency of the SVG molecule image, 0 for full transparency 1 for full opaqueness (default 0.7).
    title_col : str, optional
        name of the column in df to be used as the title entry in the hover box (default None).
    show_coords : bool, optional
        whether or not to show the coordinates of the data point in the hover box (default True).
    caption_cols : list[str], optional
        list of column names in df to be included in the hover box (default None).
    caption_transform : dict[str, callable], optional
        Functions applied to specific items in all cells. The dict must follow a key: function structure where
        the key must correspond to one of the columns in subset or tooltip (default {}).
    color_col : str, optional
        name of the column in df that is used to color the datapoints in df - necessary when there is discrete conditional coloring (default None).
    symbol_col : str, optional
        name of the column in df that is used to determine the symbols of the datapoints in df (default None).
    facet_col : str, optional
        name of the column in df that is used to facet the data to multiple plots (default None).
    wrap : bool, optional
        whether or not to wrap the title text to multiple lines if the length of the text is too long (default True).
    wraplen : int, optional
        the threshold length of the title text before wrapping begins - adjust when changing the width of the hover box (default 20).
    width : int, optional
        the width in pixels of the hover box (default 150).
    fontfamily : str, optional
        the font family used in the hover box (default 'Arial').
    fontsize : int, optional
        the font size used in the hover box - the font of the title line is fontsize+2 (default 12).
    reaction: bool, optional
        toggles rdkit to process the image like a reaction
    """
    df_data = df.copy()
    if color_col is not None:
        df_data[color_col] = df_data[color_col].astype(str)
    if symbol_col is not None:
        df_data[symbol_col] = df_data[symbol_col].astype(str)
    if facet_col is not None:
        df_data[facet_col] = df_data[facet_col].astype(str)

    if len(fig.data) != 1:
        colors = {index: x.marker["color"] for index, x in enumerate(fig.data)}

        cols = []

        if color_col is not None:
            cols.append(color_col)
        if symbol_col is not None:
            cols.append(symbol_col)
        if facet_col is not None:
            cols.append(facet_col)
        cols = list(set(cols))
        _, curve_dict = find_grouping(fig, df_data, cols)
    else:
        colors = {0: "black"}

    if not svg_height:
        svg_height = svg_size
    if not svg_width:
        svg_width = svg_size

    app = JupyterDash(__name__)
    if smiles_col is None and mol_col is None:
        raise ValueError("Either smiles_col or mol_col has to be specified!")

    if isinstance(smiles_col, str):
        smiles_col = [smiles_col]
    if isinstance(mol_col, str):
        mol_col = [mol_col]

    if mol_col is not None:
        if len(mol_col) > 1:
            menu = dcc.Dropdown(
                options=[{"label": x, "value": x} for x in mol_col],
                value=mol_col[0],
                multi=True,
                id="smiles-menu",
                placeholder="Select a mol column to display",
            )
        else:
            menu = dcc.Dropdown(
                options=[{"label": mol_col, "value": mol_col}],
                value=mol_col,
                multi=False,
                id="smiles-menu",
                disabled=True,
            )
    elif smiles_col is not None:
        if len(smiles_col) > 1:
            menu = dcc.Dropdown(
                options=[{"label": x, "value": x} for x in smiles_col],
                value=smiles_col[0],
                multi=True,
                id="smiles-menu",
                placeholder="Select a SMILES column to display",
                searchable=True,
            )
        else:
            menu = dcc.Dropdown(
                options=[{"label": smiles_col, "value": smiles_col}],
                value=smiles_col,
                multi=False,
                id="smiles-menu",
                disabled=True,
            )
    else:
        menu = dcc.Dropdown(
            options=None,
            value=None,
            id="smiles-menu",
            placeholder="Select a mol column to display",
            disabled=True,
        )

    fig_copy = go.Figure(fig)
    fig_copy.update_traces(hoverinfo="none", hovertemplate=None)

    app.layout = html.Div(
        [
            menu,
            dcc.Graph(id="graph-basic-2", figure=fig_copy, clear_on_unhover=True),
            dcc.Tooltip(
                id="graph-tooltip", background_color=f"rgba(255,255,255,{alpha})"
            ),
        ]
    )

    @app.callback(
        output=[
            Output("graph-tooltip", "show"),
            Output("graph-tooltip", "bbox"),
            Output("graph-tooltip", "children"),
        ],
        inputs=[
            Input("graph-basic-2", "hoverData"),
            Input("smiles-menu", "value"),
        ],
    )
    def display_hover(hoverData, value):
        if hoverData is None:
            return False, no_update, no_update

        if value is None:
            if mol_col is not None:
                value = mol_col
            elif smiles_col is not None:
                value = smiles_col
        if isinstance(value, str):
            chosen_smiles = [value]
        else:
            chosen_smiles = value

        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        num = pt["pointNumber"]
        curve_num = pt["curveNumber"]

        if len(fig.data) != 1:
            df_curve = curve_dict[curve_num].reset_index(drop=True)
            df_row = df_curve.iloc[num]
        else:
            df_row = df.iloc[num]

        hoverbox_elements = []

        if show_img:
            for col in chosen_smiles:
                smiles = df_row[col]
                buffered = BytesIO()
                if isinstance(smiles, str):
                    # Generate 2D SVG if smiles column is a string

                    d2d = rdMolDraw2D.MolDraw2DSVG(svg_width, svg_height)
                    opts = d2d.drawOptions()
                    opts.clearBackground = False
                    if reaction:
                        try:
                            d2d.DrawReaction(ReactionFromSmarts(smiles, useSmiles=True))
                        except:
                            d2d.DrawMolecule(Chem.MolFromSmiles(smiles))
                    else:
                        d2d.DrawMolecule(Chem.MolFromSmiles(smiles))
                    d2d.FinishDrawing()
                    img_str = d2d.GetDrawingText()
                    buffered.write(str.encode(img_str))
                    img_str = base64.b64encode(buffered.getvalue())
                    img_str = f"data:image/svg+xml;base64,{repr(img_str)[2:-1]}"

                elif isinstance(smiles, Mol):
                    # if smiles column is a Mol object, use the 3D coordinates of the mol object
                    img = Chem.Draw.MolToImage(smiles)
                    img.save(buffered, format="PNG")
                    img_str = base64.b64encode(buffered.getvalue())
                    img_str = "data:image/png;base64,{}".format(repr(img_str)[2:-1])

                else:
                    raise TypeError(
                        "smiles_col or mol_col not specified with the correct type."
                    )
                if len(smiles_col) > 1:
                    hoverbox_elements.append(
                        html.H2(
                            f"{col}",
                            style={
                                "color": colors[curve_num],
                                "font-family": fontfamily,
                                "fontSize": fontsize + 2,
                            },
                        )
                    )
                hoverbox_elements.append(
                    html.Img(
                        src=img_str,
                        style={
                            "width": "100%",
                            "background-color": f"rgba(255,255,255,{mol_alpha})",
                        },
                    )
                )

        if title_col is not None:
            title = str(df_row[title_col])
            if title_col in caption_transform:
                title = caption_transform[title_col](title)

            if len(title) > wraplen:
                if wrap:
                    title = textwrap.fill(title, width=wraplen)
                else:
                    title = title[:wraplen] + "..."

            # TODO colorbar color titles
            hoverbox_elements.append(
                html.H4(
                    f"{title}",
                    style={
                        "color": colors[curve_num],
                        "font-family": fontfamily,
                        "fontSize": fontsize,
                    },
                )
            )
        if show_coords:
            x_label = fig.layout.xaxis.title.text
            y_label = fig.layout.yaxis.title.text
            if x_label in caption_transform:
                style_str = caption_transform[x_label](pt["x"])
                hoverbox_elements.append(
                    html.P(
                        f"{x_label} : {style_str}",
                        style={
                            "color": "black",
                            "font-family": fontfamily,
                            "fontSize": fontsize,
                        },
                    )
                )
            else:
                hoverbox_elements.append(
                    html.P(
                        f"{x_label}: {pt['x']}",
                        style={
                            "color": "black",
                            "font-family": fontfamily,
                            "fontSize": fontsize,
                        },
                    )
                )
            if y_label in caption_transform:
                style_str = caption_transform[y_label](pt["y"])
                hoverbox_elements.append(
                    html.P(
                        f"{y_label} : {style_str}",
                        style={
                            "color": "black",
                            "font-family": fontfamily,
                            "fontSize": fontsize,
                        },
                    )
                )
            else:
                hoverbox_elements.append(
                    html.P(
                        f"{y_label} : {pt['y']}",
                        style={
                            "color": "black",
                            "font-family": fontfamily,
                            "fontSize": fontsize,
                        },
                    )
                )
        if caption_cols is not None:
            for caption in caption_cols:
                caption_val = df_row[caption]
                if caption in caption_transform:
                    style_str = caption_transform[caption](caption_val)
                    hoverbox_elements.append(
                        html.P(
                            f"{caption} : {style_str}",
                            style={
                                "color": "black",
                                "font-family": fontfamily,
                                "fontSize": fontsize,
                            },
                        )
                    )
                else:
                    hoverbox_elements.append(
                        html.P(
                            f"{caption} : {caption_val}",
                            style={
                                "color": "black",
                                "font-family": fontfamily,
                                "fontSize": fontsize,
                            },
                        )
                    )
        children = [
            html.Div(
                hoverbox_elements,
                style={
                    "width": f"{width}px",
                    "white-space": "normal",
                },
            )
        ]

        return True, bbox, children

    return app
