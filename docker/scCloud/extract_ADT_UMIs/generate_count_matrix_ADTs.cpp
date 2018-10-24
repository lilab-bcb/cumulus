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
#include "citeseq_utils.hpp"

using namespace std;

const int STRLEN = 1005;

struct InputFile{
	string input_r1, input_r2;

	InputFile(string r1, string r2) : input_r1(r1), input_r2(r2) {}
};

int max_mismatch, first_n;

vector<InputFile> inputs; 

Read read1, read2;
iGZipFile gzip_in_r1, gzip_in_r2;

int n_cell, n_antibody; // number of cell and antibody barcodes
int cell_blen, antibody_blen; // cell barcode length and antibody barcode length
vector<string> cell_names, antibody_names;
HashType cell_index, antibody_index;
HashIterType cell_iter, antibody_iter;

Cell2Antibody data_matrix;


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
inline bool validate_pattern(const string& tag, int pos, int lenA, int max_mismatch) {
	int nmis = (tag[pos] != 'C' && tag[pos] != 'G' && tag[pos] != 'T');
	++pos;
	for (int i = 0; i < lenA; ++i, ++pos) {
		nmis += (tag[pos] != 'A');
		if (nmis > max_mismatch) return false;
	}
	return true;
}

void produce_output(const char* out_name) {
	int i, total_umi, set_size;
	char outF[STRLEN];
	ofstream fout;
	vector<SortType> sort_arr;
	vector<vector<int> > ADTs;
	vector<int> dummy(n_antibody, 0);
	vector<int> nreads, numis;

	i = 0;
	sort_arr.clear();
	ADTs.clear();
	nreads.clear(); numis.clear();
	for (auto&& kv : data_matrix) {
		ADTs.push_back(dummy);
		total_umi = 0;
		for (auto&& kv2 : kv.second.ADTs) {
			set_size = kv2.second.size();
			ADTs[i][kv2.first] = set_size;
			total_umi += set_size;
		}
		nreads.push_back(kv.second.nreads);
		numis.push_back(total_umi);
		sort_arr.emplace_back(i++, kv.first);
	}

	sort(sort_arr.begin(), sort_arr.end());

	sprintf(outF, "%s.csv", out_name);
	fout.open(outF);
	fout<< "Antibody";
	for (auto&& val : sort_arr) fout<< ","<< cell_names[val.cell_id];
	fout<< endl;
	for (i = 0; i < n_antibody; ++i) {
		fout<< antibody_names[i];
		for (auto&& val : sort_arr) fout<< ","<< ADTs[val.pos][i];
		fout<< endl;
	}
	fout.close();

	int all_cells = sort_arr.size(), all_reads = 0, all_umis = 0;
	sprintf(outF, "%s.stat.csv", out_name);
	fout.open(outF);
	fout<< "Barcode,Total_reads,Total_umis"<< endl;
	for (auto&& val : sort_arr) {
		fout<< cell_names[val.cell_id]<< ","<< nreads[val.pos]<< ","<< numis[val.pos]<< endl;
		all_reads += nreads[val.pos];
		all_umis += numis[val.pos];
	}
	fout.close();

	printf("all_cells = %d, all_reads = %d, all_umis = %d.\n", all_cells, all_reads, all_umis);
}



int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: generate_count_matrix_ADTs cell_barcodes.txt antibody_barcodes.csv fastq_folders output_name [--max-mismatch #] [--first-n n]\n");
		exit(-1);
	}

	time_t a, b;

	a = time(NULL);

	first_n = -1;
	max_mismatch = 2;

	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "--max-mismatch")) {
			max_mismatch = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--first-n")) {
			first_n = atoi(argv[i + 1]);
		}
	}

	parse_sample_sheet(argv[1], n_cell, cell_blen, cell_index, cell_names, 1);
	parse_sample_sheet(argv[2], n_antibody, antibody_blen, antibody_index, antibody_names, max_mismatch);
	
	parse_input_directory(argv[3]);

	int cnt = 0;
	string cell_barcode, umi, antibody_barcode;
	uint64_t binary_cell, binary_umi, binary_antibody;
	
	data_matrix.clear();

	for (auto&& input_fastq : inputs) {
		gzip_in_r1.open(input_fastq.input_r1.c_str());
		gzip_in_r2.open(input_fastq.input_r2.c_str());

		while (gzip_in_r1.next(read1) == 4 && gzip_in_r2.next(read2) == 4) {
			++cnt;
			
			cell_barcode = read1.seq.substr(0, cell_blen);
			binary_cell = barcode_to_binary(cell_barcode);
			cell_iter = cell_index.find(binary_cell);

			if (cell_iter != cell_index.end() && cell_iter->second.item_id >= 0) {
				if (validate_pattern(read2.seq, antibody_blen, 7, 1)) {
					antibody_barcode = read2.seq.substr(0, antibody_blen);
					binary_antibody = barcode_to_binary(antibody_barcode);
					antibody_iter = antibody_index.find(binary_antibody);

					if (antibody_iter != antibody_index.end() && antibody_iter->second.item_id >= 0) {
						umi = read1.seq.substr(cell_blen);
						binary_umi = barcode_to_binary(umi);

						auto& one_cell = data_matrix[cell_iter->second.item_id];
						++one_cell.nreads;
						one_cell.ADTs[antibody_iter->second.item_id].insert(binary_umi);
					}
				}
			}


			if (cnt % 1000000 == 0) printf("Processed %d reads.\n", cnt);

			if (first_n > 0 && cnt >= first_n) break;
		}

		gzip_in_r1.close();
		gzip_in_r2.close();		
	}

	produce_output(argv[4]);

	b = time(NULL);
	printf("Time spent = %.2f.\n", difftime(b, a));

	return 0;
}
