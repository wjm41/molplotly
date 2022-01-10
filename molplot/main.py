from io import BytesIO
import os
import base64

from rdkit import Chem

from jupyter_dash import JupyterDash
from dash import dcc, html, Input, Output, no_update


def add_molecules(df, fig):
    fig.update_traces(hoverinfo="none", hovertemplate=None)
    
    app = JupyterDash(__name__)
    app.layout = html.Div([
        dcc.Graph(id="graph-basic-2", figure=fig, clear_on_unhover=True),
        dcc.Tooltip(id="graph-tooltip"),
    ])

    @app.callback(
        output=[Output("graph-tooltip", "show"), Output("graph-tooltip", "bbox"), Output("graph-tooltip", "children")],
        inputs=[Input("graph-basic-2", "hoverData")]
        )
    def display_hover(hoverData):
        if hoverData is None:
            return False, no_update, no_update

        # demo only shows the first point, but other points may also be available
        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        num = pt["pointNumber"]

        df_row = df.iloc[num]
        name = df_row['name']
        smiles = df_row['SMILES']
        y = df_row['y']
        x = df_row['x']

        buffered = BytesIO()
        img = Chem.Draw.MolToImage(Chem.MolFromSmiles(smiles))
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue())
        img_str = "data:image/png;base64,{}".format(
            repr(img_str)[2:-1])

        children = [
            html.Div([
                html.Img(src=img_str, style={"width": "100%"}),
                html.H2(f"{name}", style={"color": "darkblue",
                        "font-family": 'Arial'}),
                html.P(f"y : {y}", style={"color": "black", "font-family": 'Arial'}),
                html.P(f"x: {x}", style={"color": "black", "font-family": 'Arial'}),
            ], style={'width': '200px', 'white-space': 'normal'})
        ]

        return True, bbox, children
    return app
