#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os.path


def main(argsv):
    parser = argparse.ArgumentParser(description='Generate a sample sheet')
    parser.add_argument(dest='dir', help='Root directory to look for fastq files. Directory is searched recursively',
                        nargs='+')
    parser.add_argument('-f', '--format', dest='output_format', help='Format of sample sheet',
                        choices=['directory', 'r1_r2', 'r1_r2_i1'], default='directory')
    parser.add_argument('-o', '--output', dest='output', help='Output file name', required=True)
    args = parser.parse_args(argsv)
    output_format = args.output_format
    read_suffixes = ['_R1', '_R2']

    if output_format == 'r1_r2_i1':
        read_suffixes.append('_I1')

    sample_name_to_directory_to_fastqs = {}
    for dir_to_search in args.dir:
        for root, dirs, files in os.walk(os.path.abspath(dir_to_search)):
            root_path = os.path.abspath(root)
            for file_name in files:
                file_name_lc = file_name.lower()
                ext_index = file_name_lc.rfind('.fastq.gz')
                if ext_index != -1:
                    if file_name_lc.startswith('undetermined'):
                        print('Skipped ' + os.path.join(root, file_name))
                        continue
                    sample_name = file_name[0:ext_index]
                    is_fastq = False
                    for suffix_index in range(len(read_suffixes)):
                        index = sample_name.rfind(read_suffixes[suffix_index])
                        if index != -1:
                            sample_name = sample_name[0:index]
                            is_fastq = True
                            break
                    if is_fastq:
                        fastq_path = os.path.abspath(os.path.join(root, file_name))
                        directory_to_fastqs = sample_name_to_directory_to_fastqs.get(sample_name)
                        if directory_to_fastqs is None:
                            directory_to_fastqs = {}
                            sample_name_to_directory_to_fastqs[sample_name] = directory_to_fastqs
                        fastqs = directory_to_fastqs.get(root_path)
                        if fastqs is None:
                            fastqs = []
                            directory_to_fastqs[root_path] = fastqs
                        fastqs.append(fastq_path)
    if output_format == 'directory':
        with open(args.output, 'w') as writer:
            for sample_name in sample_name_to_directory_to_fastqs:
                directory_to_fastqs = sample_name_to_directory_to_fastqs[sample_name]
                for directory in directory_to_fastqs:
                    writer.write(sample_name + '\t' + directory + '\n')
    else:
        with open(args.output, 'w') as writer:
            for sample_name in sample_name_to_directory_to_fastqs:
                directory_to_fastqs = sample_name_to_directory_to_fastqs[sample_name]
                for directory in directory_to_fastqs:
                    fastqs = directory_to_fastqs[directory]
                    if len(fastqs) != len(read_suffixes):
                        raise ValueError('Missing fastq files for ' + sample_name + ' in ' + directory)
                    fastqs.sort()
                    writer.write(sample_name + '\t' + '\t'.join(fastqs) + '\n')
