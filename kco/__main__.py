#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from kco.commands import *
import argparse


def main():
    command_list = [fc_inputs, fc_run, sample_sheet, fc_upload]
    parser = argparse.ArgumentParser(description='Run a kco command')
    command_list_strings = list(map(lambda x: x.__name__[len('kco.commands.'):], command_list))
    parser.add_argument('command', help='The kco command', choices=command_list_strings)
    parser.add_argument('command_args', help='The command arguments', nargs=argparse.REMAINDER)
    kfc_args = parser.parse_args()
    command_name = kfc_args.command
    command_args = kfc_args.command_args
    cmd = command_list[command_list_strings.index(command_name)]
    sys.argv[0] = cmd.__file__
    cmd.main(command_args)


if __name__ == '__main__':
    main()
