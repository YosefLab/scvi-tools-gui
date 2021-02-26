import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_uploader as du
from dash.dependencies import Input, Output, State
import base64
import datetime
import io
import sys

import subprocess
import os
import plotly.express as px

from utils import *

import uuid

import gdown


from styles import SIDEBAR_STYLE, CONTENT_STYLE, LOGO_STYLE

import json

import scvi
import scanpy as sc

import time

sc.set_figure_params(figsize=(4, 4))


app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])


du.configure_upload(app, r"./data")

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
                dbc.NavLink("cellxgene visualization", href="/visualize", active="exact"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)

progress = html.Div(
    [
        dcc.Interval(id="progress-interval", n_intervals=0, interval=500),
        dbc.Progress(id="progress"),
    ],
    className = "mt-3"
)


builtins = {
    "cortex" : scvi.data.cortex
}

models = {
    "scvi" : scvi.model.SCVI,
    "totalvi" : scvi.model.TOTALVI
}

model_options = [
    {"label" : "SCVI", "value" : "scvi"},
    {"label" : "TOTALVI", "value" : "totalvi"}
]

content = html.Div(id="page-content", style=CONTENT_STYLE)

app.layout = html.Div([dcc.Location(id="url"), sidebar, content])



def upload_page():
    return( 
        html.Div(
            [
                html.H2("Loading anndata"),
                html.Hr(),
                
                
                html.Div(
                    id="upload-message"
                ),

                html.H5("Choose dataset"),
            
                dcc.Dropdown(
                    options=get_datasets(),
                    id = "dataset-dropdown",
                    className="mb-3"
                ),
                dbc.Button(
                    "Submit",
                    id="submit-datafile",
                    className="mb-3 mr-3 mt-3",
                    color="success",
                ),
                dbc.Button("Upload your own .h5ad data", id="open", color="primary"),

                html.Div(
                    [
                        dbc.Modal(
                            [
                                dbc.ModalHeader("Header"),
                                dbc.ModalBody([

                                    html.H5("Local files"),
                                    html.Div(
                                        [
                                            du.Upload(
                                                id='dash-uploader',
                                                upload_id=uuid.uuid1(),  # Unique session id
                                            ),                        
                                            html.Div(id='callback-output'),
                                        ],
                                        className="mb-3",
                                    ),
                                    html.H5("Google drive file link (faster): "),
                                    dbc.Row([
                                        dbc.Col(
                                            [
                                                dcc.Input(
                                                    placeholder='Dataset name',
                                                    type='text',
                                                    value='',
                                                    style = {
                                                        "width" : "100%"
                                                    },
                                                    id="name"
                                                ),
                                            ],
                                            width = 4
                                        ),
                                        dbc.Col(
                                            [
                                                dcc.Input(
                                                    placeholder='Public google drive link',
                                                    type='text',
                                                    value='',
                                                    style = {
                                                        "width" : "100%"
                                                    },
                                                    id="link"
                                                ),
                                            ]
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Button("Download", id="download", )
                                            ],
                                            width = 2
                                        )
                                    ]),
                                    html.Div(
                                        [
                                            
                                            dbc.Spinner(html.Div(id="loading-output")),
                                        ]
                                    )




                                ]),
                                dbc.ModalFooter(
                                    dbc.Button("Done", id="close", className="ml-auto")
                                ),
                            ],
                            id="modal",
                            size="xl",
                        ),
                    ]
                ),

            ]
        )
    )


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
                    ]),
                    dbc.Spinner(html.Div(id="loading-output-preprocess")),
                    dbc.Button("Submit", id="submit-button", n_clicks=0),
                    
                
            ], style={})
    )

def setup_anndata_page():
    adata = scvi.data.read_h5ad("data/preprocessed_data.h5ad")
    batch_key_names = ["None"] + list(adata.obs.keys())
    obsm_key_names = ["None"] + list(adata.obsm.keys())
    uns_key_names = ["None"] + list(adata.uns.overloaded.keys())
    
    batch_key_options = make_options(batch_key_names)
    obsm_key_options = make_options(obsm_key_names)
    uns_key_options = make_options(uns_key_names)
    return (
        html.Div([
            html.H2("Setup anndata for training"),
            html.Hr(),
            html.Div(id="setup-anndata-status", children= [], className="mb-2 mt-2"), 
            html.Label('batch_key'),
            dcc.Dropdown(
                id = "batch_key",
                options=batch_key_options,
                placeholder="None",
                value='None'
            ),
            html.Label('protein_expression_obsm_key'),
            dcc.Dropdown(
                id = "obsm_key",
                options=obsm_key_options,
                placeholder="None",
                value='None'
            ),
            html.Label('protein_names_uns_key'),
            dcc.Dropdown(
                id = "uns_key",
                options=uns_key_options,
                placeholder="None",
                value='None'
            ),
            
            dbc.Button("Submit", id="setup-anndata-submit", n_clicks=0,className="mt-3"),
            

        ], style={})

    )

def train_model_page():
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
                        options=model_options,
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
            progress,
            html.Div(id="train-model-status")
        ])
    )

def visualize_page():

    subprocess.Popen(["cellxgene", "launch", "./data/post_training_data.h5ad"])
    url = read_config("url")
    # subprocess.Popen(["cellxgene", "launch", "https://cellxgene-example-data.czi.technology/pbmc3k.h5ad"])
    return (
        html.Div(
            [
                html.Div([


                ], id="iframe-div"),
                dbc.Button("Visualize", id="visualize-button", className="mr-2"),
                html.A(dbc.Button("Open in browser"), href=url, target='_blank')
                
            ]
        )
    )

iframe_style = {
    'width': 'auto', 'height': '500',
    '-ms-zoom': '0.75',
    '-moz-transform': 'scale(0.75)',
    '-moz-transform-origin': '0 0',
    '-o-transform': 'scale(0.75)',
    '-o-transform-origin': '0 0',
    '-webkit-transform': 'scale(0.75)',
    '-webkit-transform-origin': '0 0',
}

@app.callback(
    Output("iframe-div", "children"),
    Input("visualize-button","n_clicks")
)
def visualize_callback(n):
    url = read_config("url")
    print ("Url for iframe", url)
    if n:
        return html.Iframe(src=url, style=iframe_style)

@app.callback(
    Output("upload-message", "children"),
    Input('submit-datafile', "n_clicks"),
    State("dataset-dropdown", "value")
)
def choose_dataset(n, dataset):
    if n >= 1 and dataset:
        write_config("dataset", dataset)
        return dbc.Alert("Chose "+dataset+" as dataset.", color="success", dismissable=True),

@app.callback(
    Output("loading-output", "children"), 
    [Input("download", "n_clicks")],
    [State("link","value")],
    [State("name","value")],

)
def load_output(n,url,name):
    if n >= 1:
        output = "./data/data.h5ad"
        try:
            gdown.download(url,output, quiet=False)
        except Exception as e:
            print (e)
            return "Error, try again."
        add_path("./data/data.h5ad", name)
        return "Downloaded!"

@app.callback(
    Output("modal", "is_open"),
    Output("dataset-dropdown", "options"),
    [Input("open", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open, get_datasets()
    return is_open, get_datasets()
@du.callback(
    output=Output('callback-output', 'children'),
    id='dash-uploader',
)
def get_a_list(filenames):
    for filename in filenames:
        add_path(filename)
    return html.Ul([html.Li(filenames)])


@app.callback([
    Output("status_dialog","children"),
    Output('n_genes_by_counts', 'figure'),
    Output('total_counts', 'figure'),
    Output('loading-output-preprocess', 'children'),
    Input('submit-button', 'n_clicks'),
    State("filter_genes_min_count", "value"),
    State("num_features_selected", "value"),
    State("subset", "value"),
    State("flavor", "value"),
])
def preprocess_callback(n_clicks, min_counts, num_features, subset, flavor):
    path = read_config("dataset")
    if n_clicks < 1:
        adata = None
        if path not in builtins:
            adata = scvi.data.read_h5ad(path)
        else:
            adata = builtins[path](run_setup_anndata=False)
        sc.pp.filter_genes(adata, min_counts=min_counts)
        adata.layers["counts"] = adata.X.copy()
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=int(num_features),
            subset=subset,
            layer="counts",
            flavor=flavor
        )
        adata.write_h5ad("./data/preprocessed_data.h5ad")
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        alert = dbc.Alert("Successfully preprocessed data.", color="success", dismissable=True)
        return [], px.violin(adata.obs, y="n_genes_by_counts"), px.violin(adata.obs, y="total_counts"), alert

@app.callback(
    Output("setup-anndata-status",  "children"),
    Input("setup-anndata-submit", "n_clicks"),
    State("batch_key", "value"),
    State("obsm_key", "value"),
    State("uns_key", "value"),
)
def setup_anndata_callback(n_clicks, batch_key, obsm_key, uns_key):
    if n_clicks >= 1:
        adata = scvi.data.read_h5ad("data/preprocessed_data.h5ad")
        anndata_config = {

            'batch_key' : batch_key,
            'protein_expression_obsm_key' : obsm_key,
            'protein_names_uns_key' : uns_key,
        }
        anndata_config = {key : anndata_config[key] for key in anndata_config  if anndata_config[key] and anndata_config[key] != 'None'}
        json.dump(anndata_config, open('anndata_config.json','w'))


        return [dbc.Alert("Success! Anndata setup.", dismissable=True, is_open=True, color = "success")]


@app.callback(
    [Output("progress", "value"), Output("progress", "children")],
    [Input("progress-interval", "n_intervals")],
)
def update_progress(n):
    # check progress of some background process, in this example we'll just
    # use n_intervals constrained to be in 0-100

    progress = int(get_status()["epoch"]/400 * 100)
    # only add text after 5% progress to ensure text isn't squashed too much
    return progress, f"{progress} %" if progress >= 5 else ""

@app.callback(
    Output("train-model-status", "children"),
    Input("train-button", "n_clicks"),
    State("n_hidden", "value"),
    State("model-type", "value")
)
def train_model_callback(n_clicks, n_hidden, model_type):
    model = models[model_type]
    if n_clicks >= 1:
        reset_status()
        adata = scvi.data.read_h5ad("data/preprocessed_data.h5ad")
        anndata_config = json.loads(open("anndata_config.json", "r").read())
        scvi.data.setup_anndata(
            adata, 
            **anndata_config
        )
        vae = model(adata, n_hidden=n_hidden)

        vae.train(callbacks=[ProgressCallback()])
        print ('finished training')
        vae.save("my_" + model_type + "_model/")
        print ('saved model')
        latent = vae.get_latent_representation()
        adata.obsm["X_scVI"] = latent
        sc.pp.neighbors(adata, use_rep="X_scVI")
        sc.tl.umap(adata, min_dist=0.2)
        print ('finished tl.umap')
        sc.pl.umap(
            adata, 
            color="cell_type", 
            frameon=False,
        )
        adata.write_h5ad("./data/post_training_data.h5ad")
        print ('written')
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
def toggle_collapse_callback(n, is_open):
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
        return setup_anndata_page()
    elif pathname == '/train-model':
        return train_model_page()
    elif pathname == '/visualize':
        return visualize_page()
    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    )

def run():
    app.run_server(debug=False)

def debug():
    app.run_server(debug=True)


if __name__ == "__main__":
    print ("URL cellxgene will run on",sys.argv[1])
    write_config("url", sys.argv[1])
    run()
