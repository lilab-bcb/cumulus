#ifndef CITESEQ_UTILS
#define CITESEQ_UTILS

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>



typedef std::unordered_map<int, std::unordered_set<uint64_t> > Antibody2UMI;
// typedef std::unordered_map<uint64_t, std::unordered_set<std::string> > UMI2Rest;
// typedef std::unordered_map<int, UMI2Rest> Antibody2UMI;

struct CellType {
	int nreads, numis;
	Antibody2UMI ADTs;

	CellType() { nreads = numis = 0; ADTs.clear(); }
};

typedef std::unordered_map<int, CellType> Cell2Antibody;



struct SortType {
	int pos, cell_id;

	SortType(int pos, int cell_id) : pos(pos), cell_id(cell_id) {}

	bool operator< (const SortType& o) const {
		return cell_id < o.cell_id;
	}
};

#endif
