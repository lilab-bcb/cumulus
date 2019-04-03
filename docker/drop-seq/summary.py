#!/usr/bin/env python

import argparse
import os
from subprocess import check_call

import cufflinks as cf
import numpy as np
import pandas as pd
import plotly.offline

cf.go_offline()


def get_file(file):
    if file.startswith('gs://'):
        check_call(['gsutil', '-q', '-m', 'cp', file, '.'])
        file = os.path.basename(file)
    return file


def parse_picard(path):
    path = get_file(path)
    with open(path, 'r') as reader:
        line_number = 0
        current_dict = None
        metrics_info = {'skiprows': None, 'nrows': None}
        histogram_info = {'skiprows': None, 'nrows': None}
        for line in reader:
            line = line.rstrip()
            if line == '':
                continue
            if line.startswith('##') and current_dict is not None:
                current_dict['nrows'] = line_number - current_dict['skiprows'] - 1
            if line.startswith('## METRICS CLASS'):
                current_dict = metrics_info
                current_dict['skiprows'] = line_number + 1
            elif line.startswith('## HISTOGRAM'):
                current_dict = histogram_info
                current_dict['skiprows'] = line_number + 1

            line_number = line_number + 1

    metrics = pd.read_csv(path, comment='#', sep='\t', skiprows=metrics_info['skiprows'],
                          nrows=metrics_info['nrows']) if metrics_info['skiprows'] is not None else None
    histogram = pd.read_csv(path, comment='#', sep='\t', skiprows=histogram_info['skiprows'],
                            nrows=histogram_info['nrows']) if histogram_info['skiprows'] is not None else None
    return metrics, histogram


def star(files):
    df = None
    for i in range(len(files)):
        df_i = pd.read_csv(get_file(files[i]), sep='\t', index_col=0, header=None, engine='python')
        df_i = df_i[df_i[1].isnull() == False]
        df_i.index = df_i.index.str.rstrip(' |')
        df_i.index = df_i.index.str.strip()  # e.g. % of reads mapped to too many loci
        df_i = df_i.T
        df_i.index = [sample_ids[i]]
        df = pd.concat((df, df_i), copy=False) if df is not None else df_i
    fig = df[
        ['Number of input reads', 'Uniquely mapped reads number']].iplot(
        title='Alignment Summary', kind='bar',
        asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def read_picard_files(files):
    df1 = None
    df2 = None
    for i in range(len(files)):
        df1_i, df2_i = parse_picard(files[i])
        if df1_i is not None:
            df1_i['index'] = np.arange(df1_i.shape[0])
            df1_i.index = [sample_ids[i]] * df1_i.shape[0]
            df1 = pd.concat((df1, df1_i), copy=False) if df1 is not None else df1_i
        if df2_i is not None:
            df2_i['index'] = np.arange(df2_i.shape[0])
            df2_i.index = [sample_ids[i]] * df2_i.shape[0]
            df2 = pd.concat((df2, df2_i), copy=False) if df2 is not None else df2_i
    return df1, df2


def dge_summary(files):
    df, hist = read_picard_files(files)
    df['Sample'] = df.index

    fig = df.iplot(title='# Genes Per Cell Barcode',
                   categories='Sample',
                   mode='lines+markers',
                   size=1, kind='scatter',
                   x='index',
                   y='NUM_GENES', xTitle='Cell Barcode', yTitle='# Genes',
                   asFigure=True)
    for i in range(len(fig['data'])):
        fig['data'][i]['mode'] = 'lines'
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def polya_trimmer(files):
    df, hist = read_picard_files(files)
    hist['Sample'] = hist.index
    fig = hist.iplot(title='Poly-A Trim Summary', size=6, mode='lines+markers', kind='scatter', categories='Sample',
                     x='BIN',
                     y='VALUE', xTitle='1st base of polyA tail trimmed', yTitle='# of reads',
                     asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def adapter_trimmer(files):
    df, hist = read_picard_files(files)
    hist['Sample'] = hist.index
    fig = hist.iplot(title='Adapter Trim Summary', size=6, mode='lines+markers', kind='scatter', categories='Sample',
                     x='BIN',
                     y='VALUE', xTitle='# bases trimmed 5\'', yTitle='# of reads',
                     asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def sc_rnaseq_metrics(files):
    df, hist = read_picard_files(files)
    df['Sample'] = df.index
    df = df.groupby('Sample', as_index=False).sum()

    df['genic'] = df['PCT_CODING_BASES'] + df['PCT_INTRONIC_BASES'] + df['PCT_UTR_BASES']
    df['exonic'] = df['PCT_CODING_BASES'] + df['PCT_UTR_BASES']
    df['intronic'] = df['PCT_INTRONIC_BASES']
    df['intergenic'] = df['PCT_INTERGENIC_BASES']
    fig = df[['genic', 'exonic', 'intronic', 'intergenic', 'Sample']].iplot(title='all reads',
                                                                            categories='Sample',
                                                                            kind='bar', asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def synthesis(files):
    df, hist = read_picard_files(files)
    fig = df[['NUM_BEADS', 'NO_ERROR']].iplot(title='Synthesis Error Summary', kind='bar',
                                              asFigure=True,
                                              yTitle='Count')  # asPlot=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


parser = argparse.ArgumentParser()
parser.add_argument('--dge_summary')
parser.add_argument('--sample_id')
parser.add_argument('--star_log')
parser.add_argument('--adapter_trimming_report')
parser.add_argument('--polyA_trimming_report')
parser.add_argument('--sc_rnaseq_metrics_report')
parser.add_argument('--bead_synthesis_summary')

args = parser.parse_args()
sample_ids = args.sample_id.split(',')

writer = open('drop_seq_report.html', 'w')
writer.write(
    '<html><head><script src="https://cdn.plot.ly/plotly-latest.min.js"></script><title>Drop-Seq Report</title></head>')
writer.write('<body>')

if args.dge_summary is not None:
    dge_summary(args.dge_summary.split(','))
if args.star_log is not None:
    star(args.star_log.split(','))
if args.sc_rnaseq_metrics_report is not None:
    sc_rnaseq_metrics(args.sc_rnaseq_metrics_report.split(','))
if args.bead_synthesis_summary is not None:
    synthesis(args.bead_synthesis_summary.split(','))
if args.polyA_trimming_report is not None:
    polya_trimmer(args.polyA_trimming_report.split(','))
if args.adapter_trimming_report is not None:
    adapter_trimmer(args.adapter_trimming_report.split(','))
writer.write('</body>')
writer.write('</html>')
writer.close()
