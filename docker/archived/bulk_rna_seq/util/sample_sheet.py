#!/usr/bin/env python

import argparse
import json

import pandas as pd

parser = argparse.ArgumentParser(description='Parse sample sheet')
parser.add_argument('csv', help='Path to sample sheet')
args = parser.parse_args()
input_csv_file = args.csv
df = pd.read_csv(input_csv_file, header=None, dtype=str, sep=None, engine='python', names=['name', 'r1', 'r2'])
for c in df.columns:
    df[c] = df[c].str.strip()

agg_df = df.groupby('name').agg(lambda x: list(x))
agg_df.to_csv('names.txt', columns=[], header=False)
with open('r1.json', 'wt') as f1, open('r2.json', 'wt') as f2:
    r1 = {}
    r2 = {}
    for name in agg_df.index:
        r1[name] = agg_df.loc[name, 'r1']
        r2[name] = agg_df.loc[name, 'r2']
    json.dump(r1, f1)
    json.dump(r2, f2)
