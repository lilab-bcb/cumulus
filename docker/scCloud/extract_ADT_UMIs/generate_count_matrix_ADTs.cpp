#include <ctime>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>

#include "dirent.h"

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"
#include "datamatrix_utils.hpp"

using namespace std;

const int STRLEN = 1005;

struct InputFile{
	string input_r1, input_r2;

	InputFile(string r1, string r2) : input_r1(r1), input_r2(r2) {}
};

int max_mismatch_cell, max_mismatch_feature, umi_len;
string feature_type;
int min_reads_per_umi;
double min_ratio_per_umi;

vector<InputFile> inputs; 

Read read1, read2;
iGZipFile gzip_in_r1, gzip_in_r2;

int n_cell, n_feature; // number of cell and feature barcodes
int cell_blen, feature_blen; // cell barcode length and feature barcode length
vector<string> cell_names, feature_names;
HashType cell_index, feature_index;
HashIterType cell_iter, feature_iter;

DataCollector dataCollector;

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

inline int locate_feature_start_crispr(const string& sequence, int range, const string& pattern, int max_mismatch) {
	int best_pos = -1, best_nmis = max_mismatch + 1, npat = pattern.length(), n_mis;
	range -= npat - 1;
	for (int i = 0; i < range; ++i) {
		n_mis = 0;
		for (int j = 0; j < npat; ++j) {
			n_mis += (sequence[i + j] != pattern[j]);
			if (n_mis > max_mismatch) break;
		}
		if (best_nmis > n_mis) {
			best_pos = i;
			best_nmis = n_mis;
		}
	}
	return best_pos >= 0 ? best_pos + npat : -1;
}

inline bool extract_feature_barcode(const string& sequence, int feature_length, const string& feature_type, string& feature_barcode) {
	bool success;

	if (feature_type == "antibody") {
		success = validate_pattern_antibody(sequence, feature_length, 7, 1);
		if (success) feature_barcode = sequence.substr(0, feature_length);
	}
	else {
		assert(feature_type == "crispr");
		int pos = locate_feature_start_crispr(sequence, sequence.length() - feature_length, "GGAAAGGACGAAACACCG", 1);
		success = pos >= 0;
		if (success) feature_barcode = sequence.substr(pos, feature_length);
	}

	return success;
}

int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: generate_count_matrix_ADTs cell_barcodes.txt[.gz] feature_barcodes.csv fastq_folders output_name [--max-mismatch-cell #] [--feature feature_type] [--max-mismatch-feature #] [--umi-length len] [--min-reads-per-umi min_reads] [--min-ratio-per-umi ratio]\n");
		printf("Arguments:\n\tcell_barcodes.txt[.gz]\t10x genomics barcode white list\n");
		printf("\tfeature_barcodes.csv\tfeature barcode file;barcode,feature_name\n");
		printf("\tfastq_folders\tfolder contain all R1 and R2 FASTQ files ending with 001.fastq.gz\n");
		printf("\toutput_name\toutput file name prefix;output_name.csv and output_name.stat.csv\n");
		printf("Options:\n\t--max-mismatch-cell #\tmaximum number of mismatches allowed for cell barcodes [default: 1]\n");
		printf("\t--feature feature_type\tfeature type can be either antibody or crispr [default: antibody]\n");
		printf("\t--max-mismatch-feature #\tmaximum number of mismatches allowed for feature barcodes [default: 3]\n");
		printf("\t--umi-length len\tlength of the UMI sequence [default: 10]\n");
		printf("\t--min-reads-per-umi min_reads\tminimum number of reads required to keep one UMI as true signal [default: 1]\n");
		printf("\t--min-ratio-per-umi ratio\tif one barcode-umi combination has multiple features, only features with percentage of reads in the combination > ratio [default: 0.0]\n");
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
	min_reads_per_umi = 1;
	min_ratio_per_umi = 0.0;

	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "--max-mismatch-cell")) {
			max_mismatch_cell = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--feature")) {
			feature_type = argv[i + 1];
		}
		if (!strcmp(argv[i], "--max-mismatch-feature")) {
			max_mismatch_feature = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--umi-length")) {
			umi_len = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--min-reads-per-umi")) {
			min_reads_per_umi = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--min-ratio-per-umi")) {
			min_ratio_per_umi = atof(argv[i + 1]);
		}
	}

	parse_sample_sheet(argv[1], n_cell, cell_blen, cell_index, cell_names, max_mismatch_cell);
	printf("Time spent on parsing cell barcodes = %.2fs.\n", difftime(time(NULL), a));
	parse_sample_sheet(argv[2], n_feature, feature_blen, feature_index, feature_names, max_mismatch_feature);

	parse_input_directory(argv[3]);

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
				if (extract_feature_barcode(read2.seq, feature_blen, feature_type, feature_barcode)) {
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

	dataCollector.output(argv[4], feature_type, n_feature, cell_names, umi_len, feature_names, min_reads_per_umi, min_ratio_per_umi);

	b = time(NULL);
	printf("Time spent = %.2fs.\n", difftime(b, a));

	return 0;
}
