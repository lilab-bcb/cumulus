import argparse

import pegasusio as io
import numpy as np
import pandas as pd

from collections import namedtuple
from typing import List, Dict, Tuple


demux_type_dict = {'SNG': 'singlet', 'DBL': 'doublet', 'AMB': 'unknown'}
SNP = namedtuple('SNP', ['CHROM', 'POS', 'REF', 'ALT'])

def check_colnames(fields: List[str]) -> bool:
    template = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for i, key in enumerate(template):
        if fields[i] != key:
            return False
    return True

def parse_denovo_vcf(input_vcf: str) -> Tuple[List[str], Dict[object, object]]:
    sample_names = None
    snp2geno = {}

    file_ext = input_vcf.split('.')[-1]
    if file_ext == 'gz':
        import gzip
        fin = gzip.open(input_vcf, 'rt')
    else:
        fin = open(input_vcf, 'rt')

    for line in fin:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            fields = line.strip()[1:].split('\t')
            assert check_colnames(fields)
            sample_names = fields[9:]
        else:
            fields = line.strip().split('\t')
            assert fields[8].startswith('GT')
            snp = SNP(fields[0], fields[1], fields[3], fields[4])
            snp2geno[snp] = [x.split(':')[0] for x in fields[9:]]
    fin.close()
    return sample_names, snp2geno


def calc_matching(denovo_geno: List[str], ref_geno: List[str], mmat: np.array) -> None:
    for i, r_geno in enumerate(ref_geno):
        for j, d_geno in enumerate(denovo_geno):
            mmat[i, j] += (r_geno == d_geno)

def parse_reference_vcf(reference_vcf: str, snp2geno: dict, sample_names: List[str]) -> Tuple[List[str], np.array]:
    nsample = len(sample_names)

    cnt = 0
    nbingo = 0
    ref_names = None
    mmat = None # mmat: matching matrix

    file_ext = reference_vcf.split('.')[-1]
    if file_ext == 'gz':
        import gzip
        fin = gzip.open(reference_vcf, 'rt')
    else:
        fin = open(reference_vcf, 'rt')

    for line in fin:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            fields = line.strip()[1:].split('\t')
            assert check_colnames(fields)
            ref_names = fields[9:]
            mmat = np.zeros((len(ref_names), nsample), dtype = int)
        else:
            fields = line.strip().split('\t')
            snp = SNP(fields[0], fields[1], fields[3], fields[4])
            d_geno = snp2geno.get(snp, None)
            if d_geno is not None:
                assert fields[8].startswith('GT')
                r_geno = [x.split(':')[0] for x in fields[9:]]
                calc_matching(d_geno, r_geno, mmat)
                nbingo += 1
            cnt += 1
            if cnt % 100000 == 0:
                print("Parsed {0} variants, matched {1} variants.".format(cnt, nbingo))

    fin.close()

    print("\n{0} variants are parsed and {1} SNPs are matched.".format(cnt, nbingo))

    return ref_names, mmat


def replace_ref_names(ref_str: str, ref_names: List[str], has_ref_genotypes: bool) -> List[str]:
    if ref_str is not None:
        res_arr = ref_str.split(',')
        assert len(res_arr) == len(ref_names)
        ref_names = res_arr
    elif not has_ref_genotypes:  # Rename freemuxlet's cluster names from "CLUST<num>" to "Donor<num>" format.
        ref_names = ['Donor' + str(int(x[5:]) + 1) for x in ref_names]
    ref_names = ['_ref_' + x for x in ref_names]
    return ref_names


def find_max_matching(ref_names: List[str], sample_names: List[str], mmat: np.array) -> dict:
    import itertools
    import networkx as nx
    from networkx.algorithms import bipartite

    nref = len(ref_names)
    nsample = len(sample_names)

    G = nx.Graph()
    G.add_nodes_from(ref_names, bipartite = 0)
    G.add_nodes_from(sample_names, bipartite = 1)
    G.add_weighted_edges_from([(ref_names[x], sample_names[y], -mmat[x, y]) for x, y in itertools.product(range(nref), range(nsample))])
    result = bipartite.matching.minimum_weight_full_matching(G)

    for j, sample_name in enumerate(sample_names):
        if sample_name not in result:
            i = np.where(mmat[:, j] == mmat[:, j].max())[0][0]
            result[sample_name] = ref_names[i]
            if isinstance(result[ref_names[i]], str):
                result[ref_names[i]] = [result[ref_names[i]]]
            result[ref_names[i]].append(sample_name)

    ref_n2i = {}
    for i, ref_name in enumerate(ref_names):
        ref_n2i[ref_name] = i

    for j, sample_name in enumerate(sample_names):
        i = ref_n2i[result[sample_name]]
        if mmat[i, j] != mmat[:, j].max():
            k = np.where(mmat[:, j] == mmat[:, j].max())[0][0]
            print("Warning: souporcell donor {} shares most SNPs with ref donor {}, but matches to ref donor {}!".format(sample_name, ref_names[k][5:], ref_names[i][5:]))

    print()
    for sample_name in sample_names:
        print("Souporcell donor {} matches reference donor {}.".format(sample_name, result[sample_name][5:]))
    print()

    return result

def translate_donor_name(inp_str: str, matching: dict) -> str:
    res_str = []
    for donor_id in inp_str.split(','):
        res_str.append(matching[f"CLUST{donor_id}"][5:])
    return ','.join(res_str)

def gen_summary_of_status(df):
    type_counts = df['DROPLET.TYPE'].value_counts()
    for type in demux_type_dict.keys():
        if type not in type_counts.index:
            type_counts[type] = 0
    return type_counts

def set_matching_no_reference(sample_numbers: List[str], ref_names: List[str]) -> dict:
    assert len(sample_numbers) == len(ref_names)
    matching = dict()
    pairs = list(zip(sample_numbers, ref_names))
    for (sample_number, ref_name) in pairs:
        matching[sample_number] = ref_name
        matching[ref_name] = sample_number

    return matching

def write_output(assignment_file: str, input_mat_file: str, output_zarr_file: str, matching: dict = None) -> None:
    algorithm = 'demuxlet' if matching is None else 'freemuxlet'

    df = pd.read_csv(assignment_file, sep='\t', header=0, index_col='BARCODE')
    df.index = pd.Index([x.split('-')[0] for x in df.index])

    type_counts = gen_summary_of_status(df)
    print("\nSinglets = {}, doublets = {}, unknown = {}.".format(type_counts['SNG'], type_counts['DBL'], type_counts['AMB']))

    df['demux_type'] = df['DROPLET.TYPE'].apply(lambda s: demux_type_dict[s])
    df['assignment'] = ''
    df.loc[df['demux_type']=='singlet', 'assignment'] = df.loc[df['demux_type']=='singlet', 'SNG.BEST.GUESS']
    if algorithm == 'demuxlet':
        df.loc[df['demux_type']=='doublet', 'assignment'] = df.loc[df['demux_type']=='doublet', 'DBL.BEST.GUESS'].apply(lambda s: ','.join(s.split(',')[:-1]))
    else:
        df.loc[df['demux_type']=='doublet', 'assignment'] = df.loc[df['demux_type']=='doublet', 'DBL.BEST.GUESS']
        f = np.vectorize(translate_donor_name)
        df['assignment'] = f(df['assignment'].astype('str').values, matching)

    data = io.read_input(input_mat_file)
    data.obs['demux_type'] = ''
    data.obs['assignment'] = ''

    idx = data.obs_names.isin(df.index)
    barcodes = data.obs_names[idx]
    df_valid = df.loc[barcodes, ['demux_type', 'assignment']]
    data.obs.loc[idx, 'demux_type'] = df_valid['demux_type'].values
    data.obs.loc[idx, 'assignment'] = df_valid['assignment'].values

    io.write_output(data, output_zarr_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge popscle result with gene-count matrix, matching donor names if needed.')
    parser.add_argument('demux_res', metavar='(demux_result.best | demux_result.clust1.samples.gz)', help='Demuxlet or Freemuxlet demultiplexing results.')
    parser.add_argument('raw_mat', metavar='raw_feature_bc_matrix.h5', help='Raw gene count matrix in 10x format.')
    parser.add_argument('out_file', metavar='output_result.zarr.zip', help='Output zarr file.')
    parser.add_argument('--cluster-genotypes', metavar = '(demux_result.clust1.vcf.gz | demux_result.clust1.vcf)', dest='cluster_genotypes', help = 'Genotypes detected by freemuxlet from RNA-seq reads.')
    parser.add_argument('--ref-genotypes', metavar = '(reference_genotypes.vcf.gz | reference_genotypes.vcf)', dest = 'ref_genotypes', help = 'Reference genotypes called from exome or genome sequencing data.')
    parser.add_argument('--donor-names', dest = 'ref_names', help = 'A comma-separated list containing donor names that are used to replace the ones in reference_genotypes.vcf.gz. Must match the order in the .vcf.gz file.')
    args = parser.parse_args()

    if args.cluster_genotypes is not None:  # freemuxlet
        sample_names, snp2geno = parse_denovo_vcf(args.cluster_genotypes)
        if args.ref_genotypes is not None:
            ref_names, mmat = parse_reference_vcf(args.ref_genotypes, snp2geno, sample_names)
            ref_names = replace_ref_names(args.ref_names, ref_names, has_ref_genotypes = True)
            matching = find_max_matching(ref_names, sample_names, mmat)
        else:
            ref_names = replace_ref_names(args.ref_names, sample_names, has_ref_genotypes = False)
            matching = set_matching_no_reference(sample_names, ref_names)
        write_output(args.demux_res, args.raw_mat, args.out_file, matching)
    else:  # demuxlet
        write_output(args.demux_res, args.raw_mat, args.out_file)
