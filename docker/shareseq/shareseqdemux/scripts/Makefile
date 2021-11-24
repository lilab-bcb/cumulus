.PHONY : all clean

all : shareseq_reorg_barcodes

shareseq_reorg_barcodes : shareseq_reorg_barcodes.cpp gzip_utils.hpp barcode_utils.hpp
	g++ --std=c++11 -O3 -lboost_iostreams $< -o $@

clean :
	rm -f shareseq_reorg_barcodes
