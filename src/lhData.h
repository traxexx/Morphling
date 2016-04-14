#ifndef LHDATA_H
#define LHDATA_H

#include <vector>

// ref stat single cell
typedef struct
{
	std::vector< std::vector<float> > stat;
	std::vector<float> frequency; 
} RefStat;

// single sample gl record
typedef struct
{
	std::vector<float> gl; // 0 1 2
	std::vector<float> ad; // 0 1 ambiguous
} lhRec;


// single site re-genotype record
typedef struct {
  int position;
  int mtype;
  int subtype;
  std::vector<lhRec> info;
} GtRec;

#endif