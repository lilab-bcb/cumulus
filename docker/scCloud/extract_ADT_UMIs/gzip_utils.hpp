#ifndef GZIP_UTILS
#define GZIP_UTILS

#include <string>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

struct Read {
	std::string name, seq, qual;

	std::string toString() {
		return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
	}
};

struct iGZipFile {
	std::ifstream fin;
	boost::iostreams::filtering_istream gin;
	std::string line;

	iGZipFile(const char* input_file = NULL) {
		if (input_file != NULL) open(input_file);
	}
	~iGZipFile() { close(); }

	bool open(const char* input_file) {
		fin.open(input_file, std::ios_base::in | std::ios_base::binary);
		gin.push(boost::iostreams::gzip_decompressor());
		gin.push(fin);

		return true;
	}

	bool close() {
		if (!gin.empty()) gin.reset();
		if (fin.is_open()) fin.close();
		return true;
	}

	int next(Read& aread) {
		bool success;

		success = (bool)getline(gin, aread.name);
		if (!success) return 0;
		aread.name = aread.name.substr(1);
		success = (bool)getline(gin, aread.seq);
		if (!success) return 1;
		success = (bool)getline(gin, line);
		if (!success) return 2;
		success = (bool)getline(gin, aread.qual);
		return (success ? 4 : 3);
	}
};

struct oGZipFile {
	std::ofstream fout;
	boost::iostreams::filtering_ostream gout;

	oGZipFile(const char* output_file = NULL) {
		if (output_file != NULL) open(output_file);
	}

	oGZipFile(const oGZipFile& o) {
	}

	~oGZipFile() { close(); }

	bool open(const char* output_file) {
		fout.open(output_file, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
		gout.push(boost::iostreams::gzip_compressor());
		gout.push(fout);

		return true;
	}

	bool close() {
		if (!gout.empty()) gout.reset();
		if (fout.is_open()) fout.close();
		return true;
	}

	void write(Read& aread) {
		const std::string &outstr = aread.toString();
		gout.write(outstr.c_str(), outstr.length());
	}
};

#endif
