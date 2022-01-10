#include <ctime>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

#include "dirent.h"

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"

using namespace std;

int n_barcode, barcode_len, max_mismatch;
HashType barcode_index;
HashIterType barcode_iter;
vector<string> barcode_names, barcode_seqs;

vector<string> flanking_seqs;

string sample_name, sample_type, output_directory;

struct InputFile{
	static string r1_pattern, r2_pattern, i1_pattern;

	string directory, input_r1, input_r2, input_i1, out_prefix; // out_prefix, adding directory number to distinguish FASTQ files from different directories
	InputFile(string dir, string r1, string r2, string i1, int dir_num) : directory(dir), input_r1(r1), input_r2(r2), input_i1(i1), out_prefix(sample_name + "_d" + to_string(dir_num)) {}
	
	string get_input_name(string choice) {
		if (choice == "r1") return directory + input_r1;
		if (choice == "r2") return directory + input_r2;
		return directory + input_i1;
	}

	string get_output_name(string choice, const string& output_directory) {
		if (choice == "r1") return output_directory + out_prefix + input_r1.substr(sample_name.length(), input_r1.length() - sample_name.length() - r1_pattern.length()) + "R1.fastq.gz";
		if (choice == "r2") return output_directory + out_prefix + input_r2.substr(sample_name.length(), input_r2.length() - sample_name.length() - r2_pattern.length()) + "R2.fastq.gz";
		return output_directory + out_prefix + input_i1.substr(sample_name.length(), input_i1.length() - sample_name.length() - i1_pattern.length()) + "I1.fastq.gz";
	}		
};

string InputFile::r1_pattern = string("R1_001.fastq.gz");
string InputFile::r2_pattern = string("R3_001.fastq.gz");
string InputFile::i1_pattern = string("R2_001.fastq.gz");

vector<InputFile> inputs; 

Read read1, read2, index1;
iGZipFile gzip_in_r1, gzip_in_r2, gzip_in_i1;
oGZipFile gzip_out_r1, gzip_out_r2, gzip_out_i1;

time_t start_time, finish_time;

double percent_mismatch;
int umi_len, check_polyT_len, check_polyT_max_mismatch;


void parse_flanking_csv(const char* flanking_csv) {
	iGZipFile fin(flanking_csv);
	std::string line, index_name, index_seq;
	std::size_t pos;

	flanking_seqs.clear();
	while (fin.getline(line)) {
		pos = line.find_first_of(',');
		assert(pos != string::npos);
		flanking_seqs.push_back(line.substr(pos + 1));
	}
	fin.close();
	printf("%s is parsed, n_flanking_sequences = %d.\n", flanking_csv, (int)flanking_seqs.size());

	assert(flanking_seqs.size() == 3);
}

void parse_input_directory(char* input_dirs, const string& sample_name) {
	DIR *dir;
	struct dirent *ent;
	vector<string> mate1s, mate2s, index1s;

	string dir_name;
	int dir_num = 0;

	char *input_dir = strtok(input_dirs, ",");
	
	inputs.clear();
	while (input_dir != NULL) {
		assert((dir = opendir(input_dir)) != NULL);
		
		dir_name = string(input_dir) + "/";
		++dir_num;

		mate1s.clear();
		mate2s.clear();
		index1s.clear();

		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG) {
				string file_name = string(ent->d_name);
				size_t pos;

				if (file_name.find(sample_name) == 0) {
					pos = file_name.find(InputFile::r1_pattern);
					if (pos != string::npos && pos + InputFile::r1_pattern.length() == file_name.length()) {
						mate1s.push_back(file_name);
					}
					else {
						pos = file_name.find(InputFile::r2_pattern);
						if (pos != string::npos && pos + InputFile::r2_pattern.length() == file_name.length()) {
							mate2s.push_back(file_name);
						}
						else {
							pos = file_name.find(InputFile::i1_pattern);
							if (pos != string::npos && pos + InputFile::i1_pattern.length() == file_name.length()) {
								index1s.push_back(file_name);
							}							
						}						
					}
				}
			}
		}

		int s = mate1s.size();
		assert(s == mate2s.size() && s == index1s.size());

		sort(mate1s.begin(), mate1s.end());
		sort(mate2s.begin(), mate2s.end());
		sort(index1s.begin(), index1s.end());

		for (int i = 0; i < s; ++i) {
			inputs.emplace_back(dir_name, mate1s[i], mate2s[i], index1s[i], dir_num);
		}

		input_dir = strtok(NULL, ",");
	}
}

inline bool check_flanking(const string& sequence, const string& flanking, int start, double percent_mismatch) {
	int end, max_mismatch, n_mismatch;

	end = start + flanking.length();
	if (end > sequence.length()) return false; // if flanking sequence overhang, fail

	max_mismatch = int(flanking.length() * percent_mismatch + 0.5);
	n_mismatch = 0;
	for (int i = start; i < end; ++i) {
		n_mismatch += (sequence[i] != flanking[i - start]);
	}
	return n_mismatch <= max_mismatch;
}

inline bool validate_index_read(const Read& index, string& barcode, string& qual) {
	int pos;
	string barcode_frag;
	uint64_t binary_barcode;

	pos = 0;
	barcode = "";
	qual = "";
	for (int i = 0; i < 3; ++i) {
		if (!check_flanking(index.seq, flanking_seqs[i], pos, percent_mismatch)) return false;
		pos += flanking_seqs[i].length();
		if (pos + barcode_len > index.seq.length()) return false;
		barcode_frag = index.seq.substr(pos, barcode_len);
		binary_barcode = barcode_to_binary(barcode_frag);
		barcode_iter = barcode_index.find(binary_barcode);
		if (barcode_iter == barcode_index.end() || barcode_iter->second.item_id < 0) return false;
		barcode += barcode_seqs[barcode_iter->second.item_id];
		qual += index.qual.substr(pos, barcode_len);
		pos += barcode_len;
	}

	return true;
}

inline bool check_polyT(const string& sequence, int polyT_start, int polyT_len, int polyT_max_mismatch) {
	int polyT_end = polyT_start + polyT_len;
	if (polyT_end > sequence.length()) return false;
	int n_mismatch = 0;
	for (int i = polyT_start; i < polyT_end; ++i) n_mismatch += (sequence[i] != 'T');
	return n_mismatch <= polyT_max_mismatch;
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		printf("Usage: shareseq_reorg_barcodes barcode_index.csv flanking_sequence.csv sample_name sample_type fastq_folders output_directory [--r1-pattern pattern] [--r2-pattern pattern] [--i1-pattern pattern]\n");
		printf("Arguments:\n\tbarcode_index.csv\tSHARE-Seq barcode white list, used by round1 to round3.\n");
		printf("\tflanking_sequence.csv\tFlanking sequences in front of round1 to round3 barcodes\n");
		printf("\tsample_name\tSample name. Only FASTQ files with sample_name as prefix are considered.\n");
		printf("\tsample_type\tSample type, choosing from 'gex' and 'atac'.\n");
		printf("\tfastq_folders\tfolder contain all FASTQ files ending with 001.fastq.gz\n");
		printf("\toutput_directory\tOutput all reorganized FASTQs to this folder. Please do not include slash at the end of directory name.\n");
		printf("\t[--r1-pattern pattern]\tOptional, specify the end of file name for read 1. Default to R1_001.fastq.gz.\n");
		printf("\t[--r2-pattern pattern]\tOptional, specify the end of file name for read 2. Default to R3_001.fastq.gz.\n");
		printf("\t[--i1-pattern pattern]\tOptional, specify the end of file name for index 1. Default to R2_001.fastq.gz.\n");
		exit(-1);
	}

	start_time = time(NULL);

	sample_name = string(argv[3]);
	sample_type = string(argv[4]);
	output_directory = string(argv[6]) + "/";
	max_mismatch = 2;
	percent_mismatch = 0.4; // for check_flanking, allow 40% mismatches
	umi_len = 10; // UMI length is 10
	check_polyT_len = 6; // Follow shareseq paper
	check_polyT_max_mismatch = 1; // Follow shareseq paper

	for (int i = 7; i < argc; ++i) {
		if (!strcmp(argv[i], "--r1-pattern")) { InputFile::r1_pattern = string(argv[i + 1]); ++i; }
		if (!strcmp(argv[i], "--r2-pattern")) { InputFile::r2_pattern = string(argv[i + 1]); ++i; }
		if (!strcmp(argv[i], "--i1-pattern")) { InputFile::i1_pattern = string(argv[i + 1]); ++i; }
	}

	assert(sample_type == "gex" || sample_type == "atac");

	parse_sample_sheet(argv[1], n_barcode, barcode_len, barcode_index, barcode_names, barcode_seqs, max_mismatch);

	parse_flanking_csv(argv[2]);

	parse_input_directory(argv[5], sample_name);

	int cnt = 0;
	int n_valid = 0;
	string barcode, qual;

	for (auto&& input_fastq : inputs) {
		gzip_in_r1.open(input_fastq.get_input_name("r1"));
		gzip_in_r2.open(input_fastq.get_input_name("r2"));
		gzip_in_i1.open(input_fastq.get_input_name("i1"));

		gzip_out_r1.open(input_fastq.get_output_name("r1", output_directory));
		gzip_out_r2.open(input_fastq.get_output_name("r2", output_directory));
		if (sample_type == "atac") gzip_out_i1.open(input_fastq.get_output_name("i1", output_directory));

		while (gzip_in_r1.next(read1) == 4 && gzip_in_r2.next(read2) == 4 && gzip_in_i1.next(index1) == 4) {
			++cnt;

			if (sample_type == "gex" && !check_polyT(read2.seq, umi_len, check_polyT_len, check_polyT_max_mismatch)) continue;

			if (validate_index_read(index1, barcode, qual)) {
				index1.seq = barcode;
				index1.qual = qual;
				++n_valid;

				gzip_out_r1.write(read1);
				if (sample_type == "atac") {
					gzip_out_r2.write(read2);
					gzip_out_i1.write(index1);
				} else {
					read2.seq = index1.seq + read2.seq.substr(0, umi_len);
					read2.qual = index1.qual + read2.qual.substr(0, umi_len);
					gzip_out_r2.write(read2);
				}
			}
			if (cnt % 1000000 == 0) printf("Processed %d reads, %d reads passed filtration.\n", cnt, n_valid);
		}

		gzip_in_r1.close();
		gzip_in_r2.close();
		gzip_in_i1.close();

		gzip_out_r1.close();
		gzip_out_r2.close();
		if (sample_type == "atac") gzip_out_i1.close();
	}

	finish_time = time(NULL);
	printf("Reorg for %s (%s) is finished. Kept %d reads out of %d reads. Time spent = %.2fs\n", sample_name.c_str(), sample_type.c_str(), n_valid, cnt, difftime(finish_time, start_time));

	return 0;
}
