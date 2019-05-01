#include <ctime>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#include "dirent.h"

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"
#include "datamatrix_utils.hpp"

using namespace std;

const int STRLEN = 1005;
const string TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"; // For Perturb-seq


struct InputFile{
	string input_r1, input_r2;

	InputFile(string r1, string r2) : input_r1(r1), input_r2(r2) {}
};

int max_mismatch_cell, max_mismatch_feature, umi_len;
string feature_type, extra_info;

vector<InputFile> inputs; 

Read read1, read2;
iGZipFile gzip_in_r1, gzip_in_r2;

int n_cell, n_feature; // number of cell and feature barcodes
int cell_blen, feature_blen; // cell barcode length and feature barcode length
vector<string> cell_names, feature_names;
HashType cell_index, feature_index;
HashIterType cell_iter, feature_iter;

DataCollector dataCollector;

int f[2][7]; // for banded dynamic programming, max allowed mismatch = 3



void parse_input_directory(char* input_dirs) {
	DIR *dir;
	struct dirent *ent;
	vector<string> mate1s, mate2s;

	string mate1_pattern = string("R1_001.fastq.gz");
	string mate2_pattern = string("R2_001.fastq.gz");
	string dir_name;

	char *input_dir = strtok(input_dirs, ",");
	
	inputs.clear();
	while (input_dir != NULL) {
		assert((dir = opendir(input_dir)) != NULL);
		
		dir_name = string(input_dir) + "/";

		mate1s.clear();
		mate2s.clear();

		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG) {
				string file_name = string(ent->d_name);
				size_t pos;

				pos = file_name.find(mate1_pattern);
				if (pos != string::npos && pos + mate1_pattern.length() == file_name.length()) {
					mate1s.push_back(file_name);
				}

				pos = file_name.find(mate2_pattern);
				if (pos != string::npos && pos + mate2_pattern.length() == file_name.length()) {
					mate2s.push_back(file_name);
				}
			}
		}

		int s = mate1s.size();

		assert(s == mate2s.size());
		sort(mate1s.begin(), mate1s.end());
		sort(mate2s.begin(), mate2s.end());

		for (int i = 0; i < s; ++i) {
			inputs.emplace_back(dir_name + mate1s[i], dir_name + mate2s[i]);
		}

		input_dir = strtok(NULL, ",");
	}
}

// valdiate the BA...A pattern
// lenA, length of A string
inline bool validate_pattern_antibody(const string& tag, int pos, int lenA, int max_mismatch) {
	int nmis = (tag[pos] != 'C' && tag[pos] != 'G' && tag[pos] != 'T');
	++pos;
	for (int i = 0; i < lenA; ++i, ++pos) {
		nmis += (tag[pos] != 'A');
		if (nmis > max_mismatch) return false;
	}
	return true;
}

// return rightmost position + 1
inline int matching(const string& readseq, const string& pattern, int nmax_mis, int pos, int& best_value) {
	int nmax_size = nmax_mis * 2 + 1;
	// f[x][y] : x, pattern, y, readseq
	// f[x][y] = min(f[x - 1][y - 1] + delta, f[x][y - 1] + 1, f[x - 1][y] + 1)
	int rlen = readseq.length(), plen = pattern.length();
	int prev, curr, rpos;
	int value, best_j;

	// init f[-1], do not allow insertion at the beginning
	for (int j = 0; j < nmax_size; ++j) f[1][j] = nmax_mis + 1;
	f[1][nmax_mis] = 0;

	// Dynamic Programming
	prev = 1; curr = 0;
	best_value = 0;
	int i;
	for (i = 0; i < plen; ++i) {
		best_value = nmax_mis + 1; best_j = -1;
		for (int j = 0; j < nmax_size; ++j) {
			value = nmax_mis + 1;
			rpos = pos + i + (j - nmax_mis);
			if (rpos >= 0 && rpos < rlen) value = min(value, f[prev][j] + (pattern[i] != readseq[rpos])); // match/mismatch
			if (j > 0) value = min(value, f[curr][j - 1] + 1); // insertion
			if (j + 1 < nmax_size) value = min(value, f[prev][j + 1] + 1); // deletion
			f[curr][j] = value;
			if (best_value > value) { best_value = value; best_j = j; }
		}
		if (best_value > nmax_mis) break;
		prev = curr; curr ^= 1;
	}

	return best_value <= nmax_mis ? pos + i + (best_j - nmax_mis) : -1;
}

// [start, end]
inline int locate_scaffold_sequence(const string& sequence, const string& scaffold, int start, int end, int max_mismatch) {
	int i, pos, best_value, value;

	for (i = start; i <= end; ++i) {
		pos = matching(sequence, scaffold, max_mismatch, i, best_value);
		if (pos >= 0) break;
	}

	if (best_value > 0) {
		for (int j = i + 1; j <= i + max_mismatch; ++j) {
			pos = matching(sequence, scaffold, max_mismatch, j, value);
			if (best_value > value) best_value = value, i = j;
		}
	}

	return i <= end ? i : -1;
}

// extra_info is the skeleton sequence for crispr and total-A/B/C for antibody
inline bool extract_feature_barcode(const string& sequence, int feature_length, const string& feature_type, const string& extra_info, string& feature_barcode) {
	bool success;
	int start_pos, end_pos, best_value;

	if (feature_type == "antibody") {
		if (extra_info == "TotalSeq-A") {
			success = validate_pattern_antibody(sequence, feature_length, 7, 1);
			if (success) feature_barcode = sequence.substr(0, feature_length);			
		}
		else {
			success = true;
			feature_barcode = sequence.substr(10, feature_length);
		}
	}
	else {
		start_pos = matching(sequence, TSO, 3, 0, best_value); // match template switch oligo
		success = start_pos >= 0;
		if (success) {
			end_pos = locate_scaffold_sequence(sequence, extra_info, start_pos + feature_length - max_mismatch_feature, sequence.length() - (extra_info.length() - 2), 2);
			success = end_pos >= 0;
			if (success) {
				if (end_pos - start_pos >= feature_length) 
					feature_barcode = sequence.substr(end_pos - feature_length, feature_length);
				else 
					feature_barcode = string(feature_length - (end_pos - start_pos), 'N') + sequence.substr(start_pos, end_pos - start_pos);
			}
		}
	}

	return success;
}

void detect_totalseq_type(string& extra_info) {
	int numer, cnt, denom = 1000;

	cnt = numer = 0;
	for (auto&& input_fastq : inputs) {
		gzip_in_r2.open(input_fastq.input_r2.c_str());
		while (gzip_in_r2.next(read2) == 4 && cnt < denom) {
			numer += validate_pattern_antibody(read2.seq, feature_blen, 7, 1);
			++cnt;
		}
		gzip_in_r2.close();
		if (cnt == denom) {
			double ratio = numer * 1.0 / denom;
			printf("Ratio of having poly(A) tails: %.2f.\n", ratio); 
			extra_info = (ratio > 0.5 ? "TotalSeq-A" : "TotalSeq-B or TotalSeq-C");
			return;
		}		
	}
}

int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: generate_count_matrix_ADTs cell_barcodes.txt[.gz] feature_barcodes.csv fastq_folders output_name [--max-mismatch-cell #] [--feature feature_type] [--scaffold-sequence sequence] [--max-mismatch-feature #] [--umi-length len]\n");
		printf("Arguments:\n\tcell_barcodes.txt[.gz]\t10x genomics barcode white list\n");
		printf("\tfeature_barcodes.csv\tfeature barcode file;barcode,feature_name\n");
		printf("\tfastq_folders\tfolder contain all R1 and R2 FASTQ files ending with 001.fastq.gz\n");
		printf("\toutput_name\toutput file name prefix;output_name.csv and output_name.stat.csv\n");
		printf("Options:\n\t--max-mismatch-cell #\tmaximum number of mismatches allowed for cell barcodes [default: 1]\n");
		printf("\t--feature feature_type\tfeature type can be either antibody or crispr [default: antibody]\n");
		printf("\t--scaffold-sequence sequence\tscaffold sequence used to locate the protospacer for sgRNA\n");
		printf("\t--max-mismatch-feature #\tmaximum number of mismatches allowed for feature barcodes [default: 3]\n");
		printf("\t--umi-length len\tlength of the UMI sequence [default: 10]\n");
		printf("Outputs:\n\toutput_name.csv\tfeature-cell count matrix. First row: [Antibody/CRISPR],barcode_1,...,barcode_n;Other rows: feature_name,feature_count_1,...,feature_count_n\n");
		printf("\toutput_name.stat.csv.gz\tgzipped sufficient statistics file. First row: Barcode,UMI,Feature,Count; Other rows: each row describe the read count for one barcode-umi-feature combination\n");
		exit(-1);
	}

	time_t a, b;

	a = time(NULL);

	max_mismatch_cell = 1;
	feature_type = "antibody";
	max_mismatch_feature = 3;
	umi_len = 10;
	extra_info = "";

	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "--max-mismatch-cell")) {
			max_mismatch_cell = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--feature")) {
			feature_type = argv[i + 1];
		}
		if (!strcmp(argv[i], "--scaffold-sequence")) {
			extra_info = argv[i + 1];
		}
		if (!strcmp(argv[i], "--max-mismatch-feature")) {
			max_mismatch_feature = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--umi-length")) {
			umi_len = atoi(argv[i + 1]);
		}
	}

	parse_sample_sheet(argv[1], n_cell, cell_blen, cell_index, cell_names, max_mismatch_cell);
	printf("Time spent on parsing cell barcodes = %.2fs.\n", difftime(time(NULL), a));
	parse_sample_sheet(argv[2], n_feature, feature_blen, feature_index, feature_names, max_mismatch_feature);

	parse_input_directory(argv[3]);

	if (feature_type == "antibody") {
		detect_totalseq_type(extra_info);
		printf("TotalSeq type is automatically detected as %s.\n", extra_info.c_str());
	} else {
		if (feature_type != "crispr") {
			printf("Do not support unknown feature type %s!\n", feature_type.c_str());
			exit(-1);
		}
		if (extra_info == "") {
			printf("Scaffold sequence is required for feature type crispr!\n");
			exit(-1);
		}
	}

	int cnt = 0;
	string cell_barcode, umi, feature_barcode;
	uint64_t binary_cell, binary_umi, binary_feature;
	
	dataCollector.clear();

	for (auto&& input_fastq : inputs) {
		gzip_in_r1.open(input_fastq.input_r1.c_str());
		gzip_in_r2.open(input_fastq.input_r2.c_str());

		while (gzip_in_r1.next(read1) == 4 && gzip_in_r2.next(read2) == 4) {
			++cnt;
			
			cell_barcode = read1.seq.substr(0, cell_blen);
			binary_cell = barcode_to_binary(cell_barcode);
			cell_iter = cell_index.find(binary_cell);

			if (cell_iter != cell_index.end() && cell_iter->second.item_id >= 0) {
				if (extract_feature_barcode(read2.seq, feature_blen, feature_type, extra_info, feature_barcode)) {
					binary_feature = barcode_to_binary(feature_barcode);
					feature_iter = feature_index.find(binary_feature);
					if (feature_iter != feature_index.end() && feature_iter->second.item_id >= 0) {
						umi = read1.seq.substr(cell_blen, umi_len);
						binary_umi = barcode_to_binary(umi);

						dataCollector.insert(cell_iter->second.item_id, binary_umi, feature_iter->second.item_id);
					}
				}
			}

			if (cnt % 1000000 == 0) printf("Processed %d reads.\n", cnt);
		}

		gzip_in_r1.close();
		gzip_in_r2.close();		
	}

	printf("Parsing input data is finished.\n");

	dataCollector.output(argv[4], feature_type, n_feature, cell_names, umi_len, feature_names);

	b = time(NULL);
	printf("Time spent = %.2fs.\n", difftime(b, a));

	return 0;
}
