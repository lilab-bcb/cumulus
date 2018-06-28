import argparse
import json
import os
import subprocess
import tempfile

import pkg_resources
from firecloud import api as fapi

import kco


def do_fc_inputs(method, out, config):
    method_namespace, method_name, method_version = kco.fs_split(method)
    if method_namespace is None:
        raise ValueError('Method should be specified as namespace/name (e.g. regev/drop-seq)')

    if method_version is None:
        version = -1
        methods = fapi.list_repository_methods(method_name).json()
        for method in methods:
            if method['namespace'] == method_namespace:
                version = max(version, method['snapshotId'])
        if version == -1:
            raise ValueError(method_name + ' not found')

        method_version = version
        payload = fapi.get_repository_method(method_namespace, method_name, method_version).json()['payload']
        tmp_wdl_file = tempfile.mkstemp(suffix='.wdl')[1]
        with open(tmp_wdl_file, 'w') as w:
            w.write(payload)

        repo_config = None
        if config is not None:
            config_namespace, config_name, config_version = kco.fs_split(config)
            if config_version is not None:
                repo_config = fapi.get_repository_config(config_namespace, config_name, config_version).json()
            else:
                repo_configs = fapi.list_repository_configs().json()
                for config in repo_configs:
                    if config['namespace'] == config_namespace and config['name'] == config_name and config['method'][
                        'name'] == method_name and config['method']['namespace'] == method_namespace:
                        repo_config = config
                        break
                if repo_config is not None:
                    repo_config = fapi.get_repository_config(repo_config['namespace'], repo_config['name'],
                                                             repo_config['snapshotId']).json()

        jar = pkg_resources.resource_filename('kco', 'womtool-32.jar')
        tmp_json_file = tempfile.mkstemp(suffix='.json')[1]
        with open(tmp_json_file, 'w') as tmp_json_out:
            subprocess.check_call(['java', '-jar', jar, 'inputs', tmp_wdl_file], stdout=tmp_json_out)
        with open(tmp_json_file, 'r') as tmp_json_in:
            inputs = json.loads(tmp_json_in.read())
        if repo_config is not None:
            repo_config = json.loads(repo_config['payload'])['inputs']
            for key in repo_config:
                value = repo_config[key]
                if inputs.get(key) is not None and value != '':
                    inputs[key] = value
        if out is None:
            out = method_name + '_inputs.json'
        if not out.lower().endswith('.json'):
            out = out + '.json'
        with open(out, 'w') as json_out:
            json.dump(inputs, json_out, indent=0)
        os.remove(tmp_json_file)
        os.remove(tmp_wdl_file)


def main(argsv):
    parser = argparse.ArgumentParser(
        description='Generate a FireCloud WDL input.json stub.')
    parser.add_argument('-m', '--method', dest='method', action='store', required=True,
                        help=kco.METHOD_HELP)
    parser.add_argument('-c', '--config', dest='config', action='store', required=False,
                        help='Repository config to use for generating input.json stub (e.g. regev/drop-seq-MMUL_8_0_1')
    parser.add_argument('-o', '--out', dest='out', action='store', required=True, help='JSON output file')
    args = parser.parse_args(argsv)
    do_fc_inputs(args.method, args.out, args.config)
