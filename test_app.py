import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash(__name__)

app.layout = html.Div(
    dcc.Dropdown(
        options=[
            {"label": "New York City", "value": "NYC"},
            {"label": "Montréal", "value": "MTL"},
            {"label": "San Francisco", "value": "SF"},
        ],
        value="MTL",
    )
)

if __name__ == "__main__":
    app.run_server(port=8003)
