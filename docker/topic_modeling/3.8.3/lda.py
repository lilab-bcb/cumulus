#!/bin/usr/python3

import argparse

import gensim
import numpy as np
import pandas as pd
import pegasusio as pio
import pyLDAvis.gensim
from gensim import corpora, models
from gensim.models.coherencemodel import CoherenceModel


def lda_setup(adata, prefix_exclude=('mt-', 'Rpl', 'Rps'), min_percent=2, max_percent=98):
    percents = (adata.X.getnnz(axis=0) / adata.shape[0]) * 100.0
    feature_filters = []
    if prefix_exclude is not None and len(prefix_exclude) > 0:
        regex = '|'.join(['^' + s for s in prefix_exclude])
        feature_filters.append(~(adata.var_names.str.match(regex, case=False)))
    if min_percent is not None:
        feature_filters.append(percents > min_percent)
    if max_percent is not None:
        feature_filters.append(percents < max_percent)

    # subset data to exclude genes expressed above max_percent or below min_percent of cells
    adata = adata[:, np.logical_and.reduce(feature_filters)]
    if adata.shape[1] == 0:
        raise ValueError('Filtered all features')
    print('Using {} features'.format(adata.shape[1]))
    ## Extract Sparse gbm
    # norm_count = 1e5
    # scale = norm_count / adata.X.sum(axis=1).A1
    # adata.X.data *= np.repeat(scale, np.diff(adata.X.indptr))
    # adata.X.data = np.log1p(adata.X.data) # faster than data.X.log1p()
    mat = adata.X.transpose()

    featureids = adata.var_names

    ## Create Vocab list of features
    id_list = featureids.tolist()
    out = [[]]
    for i in id_list: out.append([i])
    ## Turn into dictionary for use in model
    dictionary = corpora.Dictionary(out)

    ## Convert gbm to a corpus format for model
    corpus = gensim.matutils.Sparse2Corpus(mat)
    dictionary.save('dictionary.dict')
    corpora.MmCorpus.serialize('corpus.mm', corpus)
    adata.obs.to_csv('cell_ids.csv', columns=[])


def compute_lda(corpus, dictionary, topics, cell_ids, random_state=0):
    lda = models.LdaModel(corpus=corpus, id2word=dictionary, num_topics=topics, random_state=random_state, alpha='auto')

    ## Cell Topics
    cell_scores = lda.get_document_topics(corpus)
    cell_scores_mat = gensim.matutils.corpus2dense(cell_scores, num_terms=topics)

    ## topic by cell/documents
    cell_topics = pd.DataFrame(cell_scores_mat.T)
    cell_topics.index = cell_ids

    ## every topic by every feature
    topic_scores = pd.DataFrame(lda.get_topics()).T
    topic_scores.index = list(dictionary.values())

    lda.save('lda.model')
    cell_topics.to_csv('cell_scores.csv')
    topic_scores.to_csv('feature_topics.csv')
    coherence_model = CoherenceModel(model=lda, corpus=corpus, coherence='u_mass')
    with open('stats.txt', 'wt') as f:
        f.write('topics\tcoherence\tlog_perplexity\n')
        f.write(str(topics))
        f.write('\t')
        f.write(str(coherence_model.get_coherence()))
        f.write('\t')
        f.write(str(lda.log_perplexity(corpus)))
        f.write('\n')

    prepared_data = pyLDAvis.gensim.prepare(lda, corpus, dictionary)
    with open('report.html', 'wt') as f:
        f.write(pyLDAvis.prepared_data_to_html(prepared_data))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='sub_parser')

    parser_prepare = subparsers.add_parser('prepare', help='Prepare a dataset for input to topic modelling')
    parser_prepare.add_argument('--input', type=str, help='Anndata file to use', required=True)
    parser_prepare.add_argument('--prefix_exclude', type=str,
        help='Comma separated list of features to exclude that start with prefix. For example: mt-,Rpl,Rps')
    parser_prepare.add_argument('--min_percent', type=float, help='Exclude features expressed below min_percent')
    parser_prepare.add_argument('--max_percent', type=float, help='Exclude features expressed above max_percent')

    parser_run = subparsers.add_parser('run', help='Run topic modelling')
    parser_run.add_argument('--dictionary', type=str, required=True)
    parser_run.add_argument('--corpus', type=str, required=True)
    parser_run.add_argument('--cell_ids', type=str, required=True)
    parser_run.add_argument('--topics', type=int, default=20, help='Number of topics for LDA')
    parser_run.add_argument('--random_seed', type=int, help='Random seed', default=0)

    parser_plot = subparsers.add_parser('plot', help='Plot topic modelling stats')
    parser_plot.add_argument('stats', type=str, nargs='+')
    args = parser.parse_args()

    if args.sub_parser == 'prepare':

        prefix_exclude = None
        if args.prefix_exclude is not None:
            prefix_exclude = args.prefix_exclude.split(',')
        input_path = args.input
        d = pio.read_input(input_path)
        lda_setup(adata=d, prefix_exclude=prefix_exclude, min_percent=args.min_percent, max_percent=args.max_percent)
    elif args.sub_parser == 'run':
        dictionary = gensim.corpora.Dictionary.load(args.dictionary)
        corpus = gensim.corpora.MmCorpus(args.corpus)
        compute_lda(corpus=corpus, cell_ids=pd.read_csv(args.cell_ids, index_col=0).index.values, dictionary=dictionary,
            topics=args.topics, random_state=args.random_seed)
    elif args.sub_parser == 'plot':
        stats = []
        from matplotlib import pyplot as plt

        for f in args.stats:
            stats.append(pd.read_csv(f, sep='\t'))
        df = pd.concat(stats)
        plt.plot(df['topics'], df['coherence'], '-o')
        plt.xlabel("# Topics")
        plt.ylabel("Coherence score")
        plt.tight_layout()
        plt.savefig('coherence.png')
        plt.clf()

        plt.plot(df['topics'], df['log_perplexity'], '-o')
        plt.xlabel("# Topics")
        plt.ylabel("Log perplexity")
        plt.tight_layout()
        plt.savefig('perplexity.png')

# dictionary = gensim.corpora.Dictionary.load('dictionary.dict')
# corpus = gensim.corpora.MmCorpus('corpus.mm')
# lda = gensim.models.ldamodel.LdaModel.load('lda.model')
# lda.show_topics()
