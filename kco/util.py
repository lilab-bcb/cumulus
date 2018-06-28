import json
import os
import warnings
from typing import Tuple

from firecloud import api as fapi

warnings.filterwarnings('ignore', 'Your application has authenticated', UserWarning, 'google')


def get_or_create_workspace(workspace_namespace, workspace_name):
    ws = fapi.get_workspace(workspace_namespace, workspace_name)
    if ws.status_code == 404:
        ws = fapi.create_workspace(workspace_namespace, workspace_name)
        if ws.status_code != 201:
            raise ValueError('Unable to create workspace')
        return ws.json()
    else:
        return ws.json()['workspace']


def get_wdl_inputs(wdl_inputs):
    inputs = wdl_inputs
    if type(wdl_inputs) != dict:
        if os.path.exists(wdl_inputs):
            with open(wdl_inputs, 'r') as f:
                inputs = json.loads(f.read())
                # Filter out any key/values that contain #, and escape strings with quotes as MCs need this to not be treated as expressions
                # inputs = {k: "\"{}\"".format(v) for k, v in inputs_json.items() if '#' not in k}
        elif type(wdl_inputs) == str:
            inputs = json.loads(wdl_inputs)
        else:
            print('Unknown input type: ' + str(type(wdl_inputs)))
    else:
        inputs = wdl_inputs
    return inputs


def fs_split(s: str) -> Tuple[str, str, str]:
    """Split a FireCloud namespace/name/version string.

    Args:
        s: namespace/name/version. Version is optional

    Returns:
        Tuple of namespace, name, version
   """
    version = None
    sep = s.find('/')
    if sep == -1:
        return [None, None, None]
    namespace = s[0:sep]
    name = s[sep + 1:]
    # check for version
    sep = name.find('/')
    if sep != -1:
        version = int(name[sep + 1:])
        name = name[0:sep]
    return namespace, name, version


METHOD_HELP = 'Method namespace/name (e.g. regev/cellranger_mkfastq_count). A version can optionally be specified (e.g. regev/cell_ranger_mkfastq_count/4), otherwise the latest version of the method is used.'
