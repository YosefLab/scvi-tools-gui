import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import base64
import datetime
import io
import plotly.express as px



import json

import scvi
import scanpy as sc

sc.set_figure_params(figsize=(4, 4))


app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}
LOGO_STYLE = {
    "height" : "auto",
    "width" : "100%"
}

sidebar = html.Div(
    [
        html.Img(src="https://docs.scvi-tools.org/en/stable/_static/logo.png", style=LOGO_STYLE),
        html.Hr(),
        html.P(
            "A gui for scvi.", className="lead"
        ),
        dbc.Nav(
            [
                dbc.NavLink("Upload Data", href="/", active="exact"),
                dbc.NavLink("Preprocess", href="/preprocess", active="exact"),
                dbc.NavLink("Setup anndata", href="/setup-anndata", active="exact"),
                dbc.NavLink("Train model", href="/train-model", active="exact"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)



content = html.Div(id="page-content", style=CONTENT_STYLE)

app.layout = html.Div([dcc.Location(id="url"), sidebar, content])

def upload_page():
    return( 
        html.Div(
            [
                html.H2("Loading your .h5ad datafile"),
                html.Hr(),

                html.H5("Local files"),
                
                dcc.Upload(
                    id='upload-data',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select Files')
                    ]),
                    style={
                        'width': '100%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'margin': '10px'
                    },
                    # Allow multiple files to be uploaded
                    multiple=True
                ),
                html.H5("Google drive file link (faster): "),
                dcc.Input(
                    placeholder='Enter a google drive link.',
                    type='text',
                    value='',
                    style = {
                        "width" : "100%"
                    }
                ),
                html.H5("Load demo dataset"),
                dcc.Dropdown(
                    options=[
                        {'label' : 'dog', 'value' : "chicken"},
                    ]
                ),
                html.Div(id='output-data-upload',
                ),

            ]
        )
    )

def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    open("data_original.h5ad", 'wb').write(decoded)
    try:
        if 'h5ad' in filename:
            # Assume that the user uploaded a h5ad file
            
            adata = scvi.data.read_h5ad("data.h5ad")
            print (adata)
            
        else:
            raise Exception()
            
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file. Make sure you are loading a .h5ad anndata file.'
        ])
    

    return html.Div([
        html.H5(filename),
        
    ])
@app.callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children

def preprocess_page():
    return (
        html.Div(
            children = [
                html.H2("Preprocessing your data with scanpy"),
                html.Hr(),
                dbc.Row(
                    children = [
                        dbc.Col(dcc.Graph(id='n_genes_by_counts')),
                        dbc.Col(dcc.Graph(id='total_counts'))
                    ]
                ),
                dbc.Row(children = [
                    html.Div(id="status_dialog", className="mt-3"),
                    dbc.Col([
                        html.H5('Filter genes options'),
                        html.P('min_count'),
                        dcc.Input(
                            placeholder='3',
                            type='number',
                            value=3,
                            id="filter_genes_min_count",
                            className = "mb-3",
                        ),


                    ]), 
                    dbc.Col([
                            html.H5('Highly Variable Genes Options'),
                            html.P('Features selected'),
                            dcc.Input(
                                placeholder='2000',
                                type='number',
                                value=2000,
                                className = "mb-3",
                                id="num_features_selected"
                            ),

                            html.P('Inplace subset to highly variable genes'),
                            dcc.RadioItems(
                                options=[
                                    {'label': ' True ', 'value': "True"},
                                    {'label': ' False ', 'value': "False"},
                                ],
                                value='True',
                                className = "mb-3",
                                id="subset"
                            ),  

                            html.P('Flavor'),
                            dcc.RadioItems(
                                options=[
                                    {'label': ' seurat ', 'value': "seurat"},
                                    {'label': ' cell_ranger ', 'value': "cell_ranger"},
                                    {'label': ' seurat_v3 ', 'value': "seurat_v3"},
                                ],
                                value='seurat_v3',
                                className = "mb-3",
                                id="flavor"
                            ), 
                        ])
                    ]), dbc.Button("Submit", id="submit-button", n_clicks=0)
                
            ], style={})
    )

@app.callback([
    Output("status_dialog","children"),
    Output('n_genes_by_counts', 'figure'),
    Output('total_counts', 'figure'),
    Input('submit-button', 'n_clicks'),
    State("filter_genes_min_count", "value"),
    State("num_features_selected", "value"),
    State("subset", "value"),
    State("flavor", "value"),
])
def preprocess_callback(n_clicks, min_counts, num_features, subset, flavor):
    
    adata = scvi.data.read_h5ad("data_original.h5ad")
    sc.pp.filter_genes(adata, min_counts=min_counts)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=int(num_features),
        subset=subset,
        layer="counts",
        flavor=flavor
    )
    adata.write_h5ad("data.h5ad")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    # return [dbc.Alert(
    #     [
    #         html.H4("Success!", className="alert-heading"),
    #         html.Hr(),
    #         html.P(
    #             str(adata),
    #             className="mb-0",
    #         ),

    #     ], dismissable=True, is_open=True,
    # )], 
    return [], px.violin(adata.obs, y="n_genes_by_counts"), px.violin(adata.obs, y="total_counts")
    

def setup_anndata():
    adata = scvi.data.read_h5ad("data.h5ad")
    batch_key_names = adata.obs.keys()
    obsm_key_names = adata.obsm.keys()
    uns_key_names = adata.uns.overloaded.keys()
    create_options = lambda lis : [{'label' : i, 'value' : i} for i in lis]
    batch_key_options = create_options(batch_key_names)
    obsm_key_options = create_options(obsm_key_names)
    uns_key_options = create_options(uns_key_names)
    return (
        html.Div([
            html.H2("Setup anndata for training"),
            html.Hr(),
            html.P('batch_key'),
            dcc.Dropdown(
                id = "batch_key",
                options=batch_key_options,
                value='None'
            ),
            html.P('protein_expression_obsm_key'),
            dcc.Dropdown(
                id = "obsm_key",
                options=obsm_key_options,
                value='None'
            ),
            html.P('protein_names_uns_key'),
            dcc.Dropdown(
                id = "uns_key",
                options=uns_key_options,
                value='None'
            ),
            dbc.Button("Submit", id="setup-anndata-submit", n_clicks=0),
            html.Div(id="setup-anndata-status", children= [
                
            ]), 

        ], style={})

    )

@app.callback(
    Output("setup-anndata-status",  "children"),
    Input("setup-anndata-submit", "n_clicks"),
    State("batch_key", "value"),
    State("obsm_key", "value"),
    State("uns_key", "value"),
)
def setup_anndata_callback(n_clicks, batch_key, obsm_key, uns_key):
    adata = scvi.data.read_h5ad("data.h5ad")
    anndata_config = {

        'batch_key' : batch_key,
        'protein_expression_obsm_key' : obsm_key,
        'protein_names_uns_key' : uns_key,
    }
    anndata_config = {key : anndata_config[key] for key in anndata_config  if anndata_config[key] and anndata_config[key] != 'None'}
    json.dump(anndata_config, open('anndata_config.json','w'))


    return [dbc.Alert(
        [
            html.H4("Success!", className="alert-heading"),
            html.Hr(),
            html.P(
                str(adata),
                className="mb-0",
            ),
        ], dismissable=True, is_open=True,
    )]

def train_model():
    return (
        html.Div([


            dbc.Row([
                html.H2("Select model parameters and train model"),
                
            ]),
            dbc.Row([
                html.Hr(),
            ]),
            dbc.Row(
                [
                    dcc.Dropdown(
                        id='model-type',
                        options=[
                            {'label': 'SCVI', 'value': 'scvi'},
                            {'label': 'SCANVI', 'value': 'scanvi'},
                            {'label': 'TOTALVI', 'value': 'totalvi'},
                            {'label': 'LinearSCVI', 'value': 'linearscvi'},
                            {'label': 'AUTOZI', 'value': 'autozi'},
                            {'label': 'GIMVI', 'value': 'gimvi'},
                        ],
                        value='scvi',
                        className="mb-3",
                        style={
                            "width" : "200px"
                        }
                    )
                ]
            ),
            
            dbc.Row([
                dbc.Button(
                    "Advanced options",
                    id="collapse-button",
                    className="mb-3",
                    color="primary",
                ),
            ]),
            dbc.Collapse([
                html.Div(
                    [
                        html.H5('Model params'),
                        html.P('n_hidden'),
                        dcc.Input(
                            placeholder='128',
                            type='number',
                            value=128,
                            className = "mb-3",
                            id="n_hidden"
                        ),
                    ]
                ),
                html.Div(
                    [
                        html.H5('Training params'),
                        html.P('lr'),
                        dcc.Input(
                            placeholder='0.001',
                            type='number',
                            value=0.001,
                            className = "mb-3",
                            id="lr"
                        ),
                    ]
                ),
                ],
                id="collapse",
            ),
            dbc.Row([
                dbc.Button(
                    "Train model",
                    color="primary",
                    id="train-button",
                ),
            ]),
            html.Div(id="train-model-status")
        ])
    )

@app.callback(
    Output("train-model-status", "children"),
    Input("train-button", "n_clicks"),
    State("n_hidden", "value"),
    State("lr", "value"),
)
def train_callback(n_clicks, n_hidden, lr):
    adata = scvi.data.read_h5ad("data.h5ad")
    anndata_config = json.loads(open("anndata_config.json", "r").read())
    scvi.data.setup_anndata(
        adata, 
        **anndata_config
    )
    vae = scvi.model.SCVI(adata, n_hidden=n_hidden)
    vae.train(lr=lr)
    vae.save("my_model/")
    return [dbc.Alert(
        [
            html.H4("Success!", className="alert-heading"),
            html.Hr(),
            html.P(
                "Trained model and saved it at 'my_model'."
            ),

        ], dismissable=True, is_open=True,
    )]

@app.callback(
    Output("collapse", "is_open"),
    [Input("collapse-button", "n_clicks")],
    [State("collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname == "/":
        return upload_page()
    elif pathname == "/preprocess":
        return preprocess_page()
    elif pathname == "/setup-anndata":
        return setup_anndata()
    elif pathname == '/train-model':
        return train_model()
    # If the user tries to reach a different page, return a 404 message
    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    )


if __name__ == "__main__":
    app.run_server()
