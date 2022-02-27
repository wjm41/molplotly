from io import BytesIO
import os
import base64
import textwrap
from PIL import ImageFont

import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from jupyter_dash import JupyterDash

import plotly.express as px
from dash import dcc, html, Input, Output, no_update


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


def test_groups(fig, df_grouped):
    """Test if plotly figure curve names match up with pandas dataframe groups

    Args:
        fig (plotly figure): _description_
        groups (pandas groupby object): _description_

    Returns:
        _type_: Bool describing whether or not groups is the correct dataframe grouping descbrining the data in fig
    """
    str_groups = {}
    for name, group in df_grouped:
        # if isinstance(name, bool) or isinstance(name, int):
        #     str_groups[str(name)] = group
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


def find_grouping(fig, df_data, cols):

    if len(cols) == 1:
        df_grouped = df_data.groupby(cols)
        if not test_groups(fig, df_grouped):
            raise ValueError(
                "marker_col is mispecified because the dataframe grouping names don't match the names in the plotly figure."
            )

    elif len(cols) == 2:  # color_col and marker_col

        df_grouped_x = df_data.groupby(cols)
        df_grouped_y = df_data.groupby([cols[1], cols[0]])

        if test_groups(fig, df_grouped_x):
            df_grouped = df_grouped_x

        elif test_groups(fig, df_grouped_y):
            df_grouped = df_grouped_y
        else:
            raise ValueError(
                "color_col and marker_col are mispecified because their dataframe grouping names don't match the names in the plotly figure."
            )
    else:
        raise ValueError("Too many columns specified for grouping.")
    
    str_groups = {}
    for name, group in df_grouped:
        if isinstance(name, tuple):
            str_groups[", ".join(str(x) for x in name)] = group
        else:
            str_groups[name] = group

    curve_dict = {index: str_groups[x["name"]] for index, x in enumerate(fig.data)}
    return df_grouped, curve_dict


def add_molecules(
    fig,
    df,
    smiles_col="SMILES",
    show_img=True,
    svg_size=200,
    alpha=0.75,
    mol_alpha=0.7,
    title_col=None,
    show_coords=True,
    caption_cols=None,
    caption_transform={},
    color_col=None,
    marker_col=None,
    wrap=True,
    wraplen=20,
    width=150,
    fontfamily="Arial",
    fontsize=12,
):
    """
    A function that takes a plotly figure and a dataframe with molecular SMILES
    and returns a dash app that dynamically generates an image of molecules in the hover box
    when hovering the mouse over datapoints.
    ...

    Attributes
    ----------
    fig : plotly.graph_objects.Figure object
        a plotly figure object containing datapoints plotted from df
    df : pandas.DataFrame object
        a pandas dataframe that contains the data plotted in fig
    smiles_col : str|list[str], optional
        name of the column in df containing the smiles plotted in fig (default 'SMILES').
        If provided as a list, will add a slider to choose which column is used for rendering the structures.
    show_img : bool, optional
        whether or not to generate the molecule image in the dash app (default True)
    title_col : str, optional
        name of the column in df to be used as the title entry in the hover box (default None)
    show_coords : bool, optional
        whether or not to show the coordinates of the data point in the hover box (default True)
    caption_cols : list, optional
        list of column names in df to be included in the hover box (default None)
    caption_transform : dict, optional
        Functions applied to specific items in all cells. The dict must follow a key: function structure where
        the key must correspond to one of the columns in subset or tooltip. (default {})
    color_col : str, optional
        name of the column in df that is used to color the datapoints in df - necessary when there is discrete conditional coloring (default None)
    wrap : bool, optional
        whether or not to wrap the title text to multiple lines if the length of the text is too long (default True)
    wraplen : int, optional
        the threshold length of the title text before wrapping begins - adjust when changing the width of the hover box (default 20)
    width : int, optional
        the width in pixels of the hover box (default 150)
    fontfamily : str, optional
        the font family used in the hover box (default 'Arial')
    fontsize : int, optional
        the font size used in the hover box - the font of the title line is fontsize+2 (default 12)
    """
    fig.update_traces(hoverinfo="none", hovertemplate=None)
    df_data = df.copy()
    if color_col is not None:
        df_data[color_col] = df_data[color_col].astype(str)
    if marker_col is not None:
        df_data[marker_col] = df_data[marker_col].astype(str)
    colors = {0: "black"}

    if len(fig.data) != 1:
        if color_col is None and marker_col is None:
            raise ValueError(
                "More than one plotly curve in figure - color_col and/or marker_col needs to be specified."
            )
        if color_col is None:
            df_grouped, curve_dict = find_grouping(fig, df_data, [marker_col])
        elif marker_col is None:
            df_grouped, curve_dict = find_grouping(fig, df_data, [color_col])
        else:
            df_grouped, curve_dict = find_grouping(
                fig, df_data, [color_col, marker_col]
            )

        colors = {index: x.marker["color"] for index, x in enumerate(fig.data)}

    app = JupyterDash(__name__)
    if isinstance(smiles_col, str):
        smiles_col = [smiles_col]

    if len(smiles_col) > 1:
        slider = dcc.Slider(
            min=0,
            max=len(smiles_col) - 1,
            step=1,
            marks={i: smiles_col[i] for i in range(len(smiles_col))},
            value=0,
            id="smiles-slider",
        )
    else:
        slider = dcc.Store(id="smiles-slider", data=0)

    app.layout = html.Div(
        [
            dcc.Graph(id="graph-basic-2", figure=fig, clear_on_unhover=True),
            dcc.Tooltip(
                id="graph-tooltip", background_color=f"rgba(255,255,255,{alpha})"
            ),
            slider,
        ]
    )

    @app.callback(
        output=[
            Output("graph-tooltip", "show"),
            Output("graph-tooltip", "bbox"),
            Output("graph-tooltip", "children"),
        ],
        inputs=[Input("graph-basic-2", "hoverData"), Input("smiles-slider", "value")],
    )
    def display_hover(hoverData, value):
        if hoverData is None:
            return False, no_update, no_update
        if value is None:
            value = 0
        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        num = pt["pointNumber"]
        curve_num = pt["curveNumber"]

        # print(hoverData)
        # print(pt)

        if len(fig.data) != 1:
            # TODO replace with query
            # df_curve = df_grouped.get_group(curve_dict[curve_num]).reset_index(
            #     drop=True
            # )
            df_curve = curve_dict[curve_num].reset_index(drop=True)
            # df_curve = df[df[color_col] == curve_dict[curve_num]]
            df_row = df_curve.iloc[num]
        else:
            df_row = df.iloc[num]

        hoverbox_elements = []

        if show_img:
            # The 2D image of the molecule is generated here
            smiles = df_row[smiles_col[value]]
            buffered = BytesIO()
            img = Chem.Draw.MolToImage(Chem.MolFromSmiles(smiles))
            d2d = rdMolDraw2D.MolDraw2DSVG(svg_size, svg_size)
            opts = d2d.drawOptions()
            opts.clearBackground = False
            d2d.DrawMolecule(Chem.MolFromSmiles(smiles))
            d2d.FinishDrawing()
            img_str = d2d.GetDrawingText()
            buffered.write(str.encode(img_str))
            img_str = base64.b64encode(buffered.getvalue())

            img_str = "data:image/svg+xml;base64,{}".format(repr(img_str)[2:-1])

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
            title = df_row[title_col]
            if len(title) > wraplen:
                if wrap:
                    title = textwrap.fill(title, width=wraplen)
                else:
                    title = title[:wraplen] + "..."

            # TODO colorbar color titles
            hoverbox_elements.append(
                html.H2(
                    f"{title}",
                    style={
                        "color": colors[curve_num],
                        "font-family": fontfamily,
                        "fontSize": fontsize + 2,
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
