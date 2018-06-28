import argparse
import uuid
from typing import Union

from firecloud import api as fapi

import kco


def do_fc_run(method: str, workspace: str, wdl_inputs: Union[str, dict], upload: bool) -> str:
    """Run a FireCloud method.

    Args:
        method: method namespace/name/version. Version is optional
        workspace: workspace namespace/name
        wdl_inputs: WDL input JSON.
        upload: Whether to upload inputs and convert local file paths to gs:// URLs.

    Returns:
        URL to check submission status
   """
    inputs = kco.get_wdl_inputs(wdl_inputs)
    method_namespace, method_name, method_version = kco.fs_split(method)
    if method_version is None:
        version = -1
        methods = fapi.list_repository_methods(method_name).json()
        for method in methods:
            if method['namespace'] == method_namespace:
                version = max(version, method['snapshotId'])
        if version == -1:
            raise ValueError(method_name + ' not found')
        method_version = version

    root_entity = None
    launch_entity = None
    workspace_namespace, workspace_name, workspace_version = kco.fs_split(workspace)
    if upload:
        kco.do_fc_upload(inputs, workspace, False)
    config_namespace = workspace_namespace
    config_name = method_name + '_' + str(uuid.uuid1())
    method_body = {
        'name': config_name,
        'namespace': workspace_namespace,
        'methodRepoMethod': {'methodNamespace': method_namespace, 'methodName': method_name,
                             'methodVersion': method_version},
        'rootEntityType': root_entity,
        'prerequisites': {},
        'methodConfigVersion': 1,
        'deleted': False,
        'inputs': inputs,
        'outputs': {}
    }
    config_exists = fapi.get_workspace_config(workspace_namespace, workspace_name, config_namespace, config_name)

    if config_exists.status_code == 200:
        config_submission = fapi.update_workspace_config(workspace_namespace, workspace_name, config_namespace,
                                                         config_name, method_body)
        if config_submission.status_code != 200:
            raise ValueError('Unable to update workspace config')

    else:
        config_submission = fapi.create_workspace_config(workspace_namespace, workspace_name, method_body)
        if config_submission.status_code != 201:
            raise ValueError('Unable to create workspace config')

    launch_submission = fapi.create_submission(workspace_namespace, workspace_name, config_namespace, config_name,
                                               launch_entity, root_entity, "")

    if launch_submission.status_code == 201:
        submission_id = launch_submission.json()['submissionId']
        url = 'https://portal.firecloud.org/#workspaces/{}/{}/monitor/{}'.format(workspace_namespace, workspace_name,
                                                                                 submission_id)

        return url
    else:
        raise ValueError('Unable to launch submission - ' + str(launch_submission.json()))


def main(argsv):
    parser = argparse.ArgumentParser(
        description='Run a FireCloud method. Optionally upload files/directories to the workspace Google Cloud bucket.')
    parser.add_argument('-m', '--method', dest='method', action='store', required=True, help=kco.METHOD_HELP)
    parser.add_argument('-w', '--workspace', dest='workspace', action='store', required=True,
                        help='Workspace name (e.g. foo/bar). The workspace is created if it does not exist')
    parser.add_argument('-i', '--input', dest='wdl_inputs', action='store', required=True,
                        help='WDL input JSON.')
    parser.add_argument('-n', '--no-upload', dest='no_upload', action='store_true',
                        help='Do not upload files/directories to the workspace Google Cloud bucket')
    # parser.add_argument('-c', '--config_name', dest='config_name', action='store', required=False,
    #                     help='Method configuration name')
    args = parser.parse_args(argsv)
    url = do_fc_run(args.method, args.workspace, args.wdl_inputs, not args.no_upload)
    print(url)
