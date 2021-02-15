import json
from pytorch_lightning.callbacks import Callback

def get_paths():
    
    with open("files.json") as f:
        return json.load(f)

def get_datasets():
    paths = get_paths()
    options = []
    for k in paths:
        options.append(make_option(k, paths[k]))
    return options



class ProgressCallback(Callback):

    def on_epoch_end(self, trainer, pl_module):
        with open("status.json") as f:

            status = json.load(f)
            status["epoch"] += 1
        with open("status.json", 'w') as f:
            json.dump(status, f)

    def on_train_end(self, trainer, pl_module):
        with open("status.json") as f:

            status = json.load(f)
            status["trained"] = "True"
        with open("status.json", 'w') as f:
            json.dump(status, f)

def write_config(key, value):
    assert key and value, "key and value must not be None"
    with open("config.json") as f:
        conf = json.load(f)
        conf[key] = value
    with open("config.json",'w') as f:
        json.dump(conf, f)
    return True

def make_option(label, value):
    return {'label':label, 'value':value}

def add_path(path, name=None):
    
    data = get_paths()
    if not name:

        data[path.split('/')[-1]] = path
    else:
        data[name] = path

    with open("files.json", 'w') as f:
        json.dump(data, f)

    return data