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


def star_report(files):
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


def lines_only(fig):
    for i in range(len(fig['data'])):
        fig['data'][i]['mode'] = 'lines'


def dge_summary_multi_species_report(file1, file2):
    sample1 = os.path.basename(file1)
    if sample1.endswith('_dge.summary.txt'):
        sample1 = sample1[0:len(sample1) - len('_dge.summary.txt')]
    sample2 = os.path.basename(file2)
    if sample2.endswith('_dge.summary.txt'):
        sample2 = sample2[0:len(sample2) - len('_dge.summary.txt')]
    df1, hist1 = parse_picard(file1)
    df1 = df1.set_index('CELL_BARCODE')

    df2, hist2 = parse_picard(file2)
    df2 = df2.set_index('CELL_BARCODE')
    df = df1.join(df2, how='outer', rsuffix='_2')
    df = df.fillna(0)

    fig = df.iplot(title='Genes Per Cell Barcode',
        mode='markers',
        size=5,
        kind='scatter',
        x='NUM_GENES',
        y='NUM_GENES_2',
        xTitle='{} # Genes'.format(sample1),
        yTitle='{} # Genes'.format(sample2),
        asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))

    fig = df.iplot(title='Genic Reads Per Cell Barcode',
        mode='markers',
        size=5,
        kind='scatter',
        x='NUM_GENIC_READS',
        y='NUM_GENIC_READS_2',
        xTitle='{} # Genic reads'.format(sample1),
        yTitle='{} # Genic reads'.format(sample2),
        asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))

    fig = df.iplot(title='Transcripts Per Cell Barcode',
        mode='markers',
        size=5,
        kind='scatter',
        x='NUM_TRANSCRIPTS',
        y='NUM_TRANSCRIPTS_2',
        xTitle='{} # Transcripts'.format(sample1),
        yTitle='{} # Transcripts'.format(sample2),
        asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def dge_summary_report(files):
    df, hist = read_picard_files(files)
    df['Sample'] = df.index

    fig = df.iplot(title='Genes Per Cell Barcode',
        categories='Sample',
        mode='lines',
        size=1, kind='scatter',
        x='index',
        y='NUM_GENES',
        xTitle='Cell Barcode',
        yTitle='# Genes',
        asFigure=True)

    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))

    fig = df.iplot(title='Genic Reads Per Cell Barcode',
        categories='Sample',
        mode='lines',
        size=1,
        kind='scatter',
        x='index',
        y='NUM_GENIC_READS',
        xTitle='Cell Barcode',
        yTitle='# Genic Reads',
        asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))

    fig = df.iplot(title='Transcripts Per Cell Barcode',
        categories='Sample',
        mode='lines',
        size=1,
        kind='scatter',
        x='index',
        y='NUM_TRANSCRIPTS',
        xTitle='Cell Barcode',
        yTitle='# Transcripts',
        asFigure=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def polya_trimmer_report(files):
    df, hist = read_picard_files(files)
    hist['Sample'] = hist.index
    fig = hist.iplot(title='Poly-A Trim Summary', size=6, mode='lines+markers', kind='scatter', categories='Sample',
        x='BIN',
        y='VALUE', xTitle='1st base of polyA tail trimmed', yTitle='# of reads',
        asFigure=True)
    lines_only(fig)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def adapter_trimmer_report(files):
    df, hist = read_picard_files(files)
    hist['Sample'] = hist.index
    fig = hist.iplot(title='Adapter Trim Summary', size=6, mode='lines+markers', kind='scatter', categories='Sample',
        x='BIN',
        y='VALUE', xTitle='# bases trimmed 5\'', yTitle='# of reads',
        asFigure=True)
    lines_only(fig)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def sc_rnaseq_metrics_report(files):
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


def synthesis_report(files):
    df, hist = read_picard_files(files)
    fig = df[['NUM_BEADS', 'NO_ERROR']].iplot(title='Synthesis Error Summary', kind='bar',
        asFigure=True,
        yTitle='Count')  # asPlot=True)
    writer.write(plotly.offline.plot(fig, auto_open=False, include_plotlyjs=False, output_type='div'))


def get_arg(s):
    if s is None:
        return []
    return list(filter(lambda x: x != '', s.strip().split(',')))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dge_summary')
    parser.add_argument('--dge_summary_multi_species')
    parser.add_argument('--sample_id')
    parser.add_argument('--star_log')
    parser.add_argument('--adapter_trimming_report')
    parser.add_argument('--polyA_trimming_report')
    parser.add_argument('--sc_rnaseq_metrics_report')
    parser.add_argument('--bead_synthesis_summary')

    args = parser.parse_args()
    sample_ids = get_arg(args.sample_id)

    writer = open('drop_seq_report.html', 'w')
    writer.write(
        '<html><head><script src="https://cdn.plot.ly/plotly-latest.min.js"></script><title>Drop-Seq Report</title></head>')
    writer.write('<body>')

    tokens = get_arg(args.dge_summary)
    if len(tokens) > 0:
        dge_summary_report(tokens)
    if args.dge_summary_multi_species is not None and args.dge_summary_multi_species != '':
        # tsv file where each line is an array of species for one sample
        with open(args.dge_summary_multi_species, 'rt') as f:
            for line in f:
                line = line.strip()
                if line != '':
                    sample = line.split('\t')
                    for i in range(len(sample)):
                        for j in range(i):
                            dge_summary_multi_species_report(sample[i], sample[j])

    tokens = get_arg(args.star_log)
    if len(tokens) > 0:
        star_report(tokens)

    tokens = get_arg(args.sc_rnaseq_metrics_report)
    if len(tokens) > 0:
        sc_rnaseq_metrics_report(tokens)

    tokens = get_arg(args.bead_synthesis_summary)
    if len(tokens) > 0:
        synthesis_report(tokens)

    tokens = get_arg(args.polyA_trimming_report)
    if len(tokens) > 0:
        polya_trimmer_report(tokens)

    tokens = get_arg(args.adapter_trimming_report)
    if len(tokens) > 0:
        adapter_trimmer_report(tokens)
    writer.write('</body>')
    writer.write('</html>')
    writer.close()
