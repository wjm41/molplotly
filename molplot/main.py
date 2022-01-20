from io import BytesIO
import os
import base64

from rdkit import Chem

from jupyter_dash import JupyterDash

import plotly.express as px
from dash import dcc, html, Input, Output, no_update


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


def add_molecules(fig, df, smiles_col='SMILES', show_img=True, title_col=None,  show_coords=True, caption_cols=None, condition_col=None):
    fig.update_traces(hoverinfo="none", hovertemplate=None)

    colors = {0: 'black'}
    if len(fig.data) != 1:
        colors = px.colors.qualitative.Plotly
        if condition_col is not None:
            if df[condition_col].dtype == bool:
                curve_dict = {index: str2bool(x['name'])
                              for index, x in enumerate(fig.data)}
            else:
                curve_dict = {index: x['name']
                              for index, x in enumerate(fig.data)}
        else:
            raise ValueError(
                'condition_col needs to be specified if there is more than one plotly curve in the figure!')

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
            df_curve = df[df[condition_col] ==
                          curve_dict[curve_num]].reset_index(drop=True)
            df_row = df_curve.iloc[num]
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

            hoverbox_elements.append(
                html.H2(f"{df_row[title_col]}", style={"color": colors[curve_num],
                        "font-family": 'Arial'}))
        if show_coords:
            hoverbox_elements.append(html.P(f"{fig.layout.xaxis.title.text}: {pt['x']}",
                                            style={"color": "black", "font-family": 'Arial'}))
            hoverbox_elements.append(html.P(f"{fig.layout.yaxis.title.text} : {pt['y']}",
                                            style={"color": "black", "font-family": 'Arial'}))
        if caption_cols is not None:
            for caption in caption_cols:
                caption_val = df_row[caption]
                hoverbox_elements.append(html.P(f"{caption} : {caption_val}",
                                                style={"color": "black", "font-family": 'Arial'}))
        children = [html.Div(hoverbox_elements, style={
            'width': '200px', 'white-space': 'normal'})]

        return True, bbox, children
    return app


# def add_molecules(df, fig, condition):
#     colors = px.colors.qualitative.Plotly
#     fig.update_traces(hoverinfo="none", hovertemplate=None)
#     if df[condition].dtype == bool:
#         curve_dict = {index: str2bool(x['name'])
#                       for index, x in enumerate(fig.data)}
#     # elif df[condition].dtype == str:

#     else:
#         curve_dict = {index: x['name'] for index, x in enumerate(fig.data)}
#         # raise Exception('Incompatible Datatype!')
#     app = JupyterDash(__name__)
#     app.layout = html.Div([
#         dcc.Graph(id="graph-basic-2", figure=fig, clear_on_unhover=True),
#         dcc.Tooltip(id="graph-tooltip"),
#     ])

#     @app.callback(
#         output=[Output("graph-tooltip", "show"), Output("graph-tooltip",
#                                                         "bbox"), Output("graph-tooltip", "children")],
#         inputs=[Input("graph-basic-2", "hoverData")]
#     )
#     def display_hover(hoverData):
#         if hoverData is None:
#             return False, no_update, no_update

#         # demo only shows the first point, but other points may also be available
#         pt = hoverData["points"][0]
#         # print(hoverData)
#         bbox = pt["bbox"]
#         num = pt["pointNumber"]

#         # print(curve_dict)
#         # print(curve_dict[pt['curveNumber']])
#         # print(df[curve_dict[pt['curveNumber']]])
#         # print(curve_dict[pt['curveNumber']])
#         # print(df[df[condition]==curve_dict[pt['curveNumber']]])
#         curve_num = pt['curveNumber']
#         # print(curve_num)
#         df_curve = df[df[condition] ==
#                       curve_dict[curve_num]].reset_index(drop=True)
#         # print(df_curve)
#         # print(num)
#         df_row = df_curve.iloc[num]
#         # print(df_row.to_string())
#         name = df_row['Molecule Name']
#         smiles = df_row['SMILES']
#         pred = df_row['50th quantile']
#         measurement = df_row['inhibition']

#         buffered = BytesIO()
#         img = Draw.MolToImage(Chem.MolFromSmiles(smiles))
#         img.save(buffered, format="PNG")
#         img_str = base64.b64encode(buffered.getvalue())
#         img_str = "data:image/png;base64,{}".format(
#             repr(img_str)[2:-1])

#         children = [
#             html.Div([
#                 html.Img(src=img_str, style={"width": "100%"}),
#                 html.H2(f"{name}", style={"color": colors[curve_num],
#                         "font-family": 'Arial'}),
#                 html.P(f"Prediction (%): {pred:.1f}", style={"color": "black",
#                                                              "font-family": 'Arial'}),
#                 html.P(f"Measurement (%): {measurement:.1f}", style={"color": "black",
#                                                                      "font-family": 'Arial'}),
#             ], style={'width': '200px', 'white-space': 'normal'})
#         ]

#         return True, bbox, children
#     return app
