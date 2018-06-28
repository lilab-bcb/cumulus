#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
import sys

parser = argparse.ArgumentParser(
    description='Checks a sample sheet to make sure all files exist and sample ids are unique')
parser.add_argument('file', help='The sample sheet')
parser.add_argument('--remove', help='Optionally remove invalid lines and print the new sample sheet to standard out',
                    action='store_true')

args = parser.parse_args()
sheet = pd.read_table(args.file, header=None)
ids = {}
nerrors = 0
ncols = sheet.shape[1]
file_path_cols = [1, 2] if ncols is 3 else [0, 1]
has_ids = True if ncols is 3 else False
for i in range(sheet.shape[0]):
    line_has_error = False
    if has_ids:
        if ids.get(sheet.iloc[i, 0]) is not None:
            print(sheet.iloc[i, 0] + ' is a duplicate id', file=sys.stderr)
            line_has_error = True
            nerrors = nerrors + 1
        ids[id] = True
    for j in file_path_cols:
        if not os.path.exists(sheet.iloc[i, j]):
            print(sheet.iloc[i, j] + ' not found', file=sys.stderr)
            line_has_error = True
            nerrors = nerrors + 1
        elif os.path.getsize(sheet.iloc[i, j]) is 0:
            print(sheet.iloc[i, j] + ' is empty', file=sys.stderr)
            line_has_error = True
            nerrors = nerrors + 1
    if args.remove and not line_has_error:
        print('\t'.join(sheet.iloc[i]))

if nerrors > 0:
    print(str(nerrors) + ' error' + ('s' if nerrors > 1 else '') + ' found', file=sys.stderr)
