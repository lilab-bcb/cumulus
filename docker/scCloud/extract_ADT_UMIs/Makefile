.PHONY : all clean

all : generate_count_matrix_ADTs

generate_count_matrix_ADTs : generate_count_matrix_ADTs.cpp barcode_utils.hpp citeseq_utils.hpp gzip_utils.hpp
	g++ --std=c++11 -O3 -lboost_iostreams $< -o $@

clean :
	rm -f generate_count_matrix_ADTs
