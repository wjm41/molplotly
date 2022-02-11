from io import BytesIO
import os
import base64
import textwrap
from PIL import ImageFont

from rdkit import Chem

from jupyter_dash import JupyterDash

import plotly.express as px
from dash import dcc, html, Input, Output, no_update


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


def add_molecules(fig,
                  df,
                  smiles_col='SMILES',
                  show_img=True,
                  title_col=None,
                  show_coords=True,
                  caption_cols=None,
                  caption_transform={},
                  color_col=None,
                  wrap=True,
                  wraplen=20,
                  width=150,
                  fontfamily='Arial',
                  fontsize=12):
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
    smiles_col : str, optional
        name of the column in df containing the smiles plotted in fig (default 'SMILES')
    show_img : bool, optional
        whether or not to generate the molecule image in the dash app (default True)
    title_col : str, optional
        name of the column in df to be used as the title entry in the hover box (default None)
    show_coords : bool, optional
        whether or not to show the coordinates of the data point in the hover box (default True)
    caption_cols : list, optional
        list of column names in df to be included in the hover box (default None)
    caption_transform : dict, optional
        Functions applied to specific items in all cells. The dict must follow a key: function structure where the key must correspond to one of the columns in subset or tooltip. (default {})
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

    colors = {0: 'black'}
    if len(fig.data) != 1:
        if color_col is not None:
            colors = {index: x.marker['color']
                      for index, x in enumerate(fig.data)}
            if df[color_col].dtype == bool:
                curve_dict = {index: str2bool(x['name'])
                              for index, x in enumerate(fig.data)}
            elif df[color_col].dtype == int:
                curve_dict = {index: int(x['name'])
                              for index, x in enumerate(fig.data)}
            else:
                curve_dict = {index: x['name']
                              for index, x in enumerate(fig.data)}
        else:
            raise ValueError(
                'color_col needs to be specified if there is more than one plotly curve in the figure!')

    app = JupyterDash(__name__)
    app.layout = html.Div([
        dcc.Graph(id="graph-basic-2", figure=fig, clear_on_unhover=True),
        dcc.Tooltip(id="graph-tooltip"),
    ])

    @app.callback(
        output=[Output("graph-tooltip", "show"), Output("graph-tooltip",
                                                        "bbox"), Output("graph-tooltip", "children")],
        inputs=[Input("graph-basic-2", "hoverData")]
    )
    def display_hover(hoverData):
        if hoverData is None:
            return False, no_update, no_update

        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        num = pt["pointNumber"]
        curve_num = pt['curveNumber']

        if len(fig.data) != 1:
            df_curve = df[df[color_col] ==
                          curve_dict[curve_num]].reset_index(drop=True)
            df_row = df_curve.iloc[num]
        else:
            df_row = df.iloc[num]

        hoverbox_elements = []

        if show_img:
            # The 2D image of the molecule is generated here
            smiles = df_row[smiles_col]
            buffered = BytesIO()
            img = Chem.Draw.MolToImage(Chem.MolFromSmiles(smiles))
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue())
            img_str = "data:image/png;base64,{}".format(
                repr(img_str)[2:-1])
            hoverbox_elements.append(
                html.Img(src=img_str, style={"width": "100%"}))

        if title_col is not None:
            title = df_row[title_col]
            if len(title) > wraplen:
                if wrap:
                    title = textwrap.fill(title, width=wraplen)
                else:
                    title = title[:wraplen] + '...'
            hoverbox_elements.append(
                html.H2(f"{title}", style={"color": colors[curve_num],
                        "font-family": fontfamily, "fontSize": fontsize+2}))
        if show_coords:
            x_label = fig.layout.xaxis.title.text
            y_label = fig.layout.yaxis.title.text
            if x_label in caption_transform:
                style_str = caption_transform[x_label](pt['x'])
                hoverbox_elements.append(html.P(f"{x_label} : {style_str}",
                                                style={"color": "black", "font-family": fontfamily, "fontSize": fontsize}))
            else:
                hoverbox_elements.append(html.P(f"{x_label}: {pt['x']}",
                                                style={"color": "black", "font-family": fontfamily, "fontSize": fontsize}))
            if y_label in caption_transform:
                style_str = caption_transform[y_label](pt['y'])
                hoverbox_elements.append(html.P(f"{y_label} : {style_str}",
                                                style={"color": "black", "font-family": fontfamily, "fontSize": fontsize}))
            else:
                hoverbox_elements.append(html.P(f"{y_label} : {pt['y']}",
                                                style={"color": "black", "font-family": fontfamily, "fontSize": fontsize}))
        if caption_cols is not None:
            for caption in caption_cols:
                caption_val = df_row[caption]
                if caption in caption_transform:
                    style_str = caption_transform[caption](caption_val)
                    hoverbox_elements.append(html.P(f"{caption} : {style_str}",
                                                    style={"color": "black", "font-family": fontfamily, "fontSize": fontsize}))
                else:
                    hoverbox_elements.append(html.P(f"{caption} : {caption_val}",
                                                    style={"color": "black", "font-family": fontfamily, "fontSize": fontsize}))
        children = [html.Div(hoverbox_elements, style={
            'width': f'{width}px', 'white-space': 'normal'})]

        return True, bbox, children
    return app
