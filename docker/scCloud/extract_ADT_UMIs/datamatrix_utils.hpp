#ifndef DATAMATRIX_UTILS
#define DATAMATRIX_UTILS

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <unordered_map>

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"

typedef std::unordered_map<int, int> Feature2Count;
typedef std::unordered_map<uint64_t, Feature2Count> UMI2Feature;
typedef std::unordered_map<int, UMI2Feature> Cell2UMI;

class DataCollector {
public:
	DataCollector() { clear(); }

	void clear() { data_container.clear(); }

	void insert(int cell_id, uint64_t umi, int feature_id) {
		++data_container[cell_id][umi][feature_id];
	}

	void output(const std::string& output_name, const std::string& feature_type, int n_feature, const std::vector<std::string>& cell_names, int umi_len, const std::vector<std::string>& feature_names, int min_reads_per_umi, double min_ratio_per_umi) {
		std::vector<int> cell_ids;
		std::ofstream fout;
		oGZipFile fstat;

		cell_ids.clear();
		for (auto&& kv : data_container)
			cell_ids.push_back(kv.first);
		if (cell_ids.size() > 1) std::sort(cell_ids.begin(), cell_ids.end());

		std::vector<int> dummy(cell_ids.size(), 0), tot_umis(cell_ids.size(), 0);
		std::vector<std::vector<int> > ADTs(n_feature, dummy);
		double denom;
		int total_reads = 0, total_umis = 0;

		fstat.open(output_name + ".stat.csv.gz");
		fstat() << "Barcode,UMI,Feature,Count"<< std::endl;
		for (int i = 0; i < (int)cell_ids.size(); ++i) {
			auto& one_cell = data_container[cell_ids[i]];
			for (auto&& kv1 : one_cell) {
				denom = 0.0;
				for (auto&& kv2 : kv1.second) denom += kv2.second;
				for (auto&& kv2 : kv1.second) {
					fstat()<< cell_names[cell_ids[i]]<< ","<<  binary_to_barcode(kv1.first, umi_len)<< ","<< feature_names[kv2.first]<< ","<< kv2.second<< std::endl;
					total_reads += kv2.second;
					++total_umis;
					if (kv2.second >= min_reads_per_umi && kv2.second / denom > min_ratio_per_umi) {
						++ADTs[kv2.first][i];
						++tot_umis[i];
					}
				}
			}
		}
		fstat.close();
		printf("%s.stat.csv.gz is written.\n", output_name.c_str());

		fout.open(output_name + ".csv");
		fout<< (feature_type == "antibody" ? "Antibody" : "CRISPR");
		for (int i = 0; i < (int)cell_ids.size(); ++i)
			if (tot_umis[i] > 0) fout<< ","<< cell_names[cell_ids[i]];
		fout<< std::endl;
		for (int i = 0; i < n_feature; ++i) {
			fout<< feature_names[i];
			for (int j = 0; j < (int)cell_ids.size(); ++j) 
				if (tot_umis[j] > 0) fout<< ","<< ADTs[i][j];
			fout<< std::endl;
		}
		fout.close();
		printf("%s.csv is written.\n", output_name.c_str());

		printf("Sequencing saturation = %.2f%%, total_umis = %d, total_reads = %d.\n", (100.0 - total_umis * 100.0 / total_reads), total_umis, total_reads);
	}

private:
	Cell2UMI data_container;
};

#endif
