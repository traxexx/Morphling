#ifndef REGENOTYPE_H
#define REGENOTYPE_H

#include "VcfRecord.h"
#include "FastaUtil.h"
#include "ProperDeck.h"
#include "DiscPair.h"
#include "RefStats.h"
#include "SamFile.h"
#include "SamFileHeader.h"

using std::vector;
using std::string;

// re-genotype and update a vcf field from reads in bam file
void ResetVcfRecordFromBam( VcfRecord & vcf_rec, RefStats & rstats, vector<RefSeq*> & REF_SEQ, string & chr, int center, SamFile & samIn, SamFileHeader & samHeader );

// set read-type count in a specific window from bam file
void setReadCountInSection( vector<int> & raw_counts, string & chr, int center, SamFile & samIn, SamFileHeader & samHeader, vector<RefSeq*> & REF_SEQ );

// return index to raw counts in proper reads
int getRetrievedIndexToRawCounts( RetrievedIndex & rv );

// return index to raw counts in disc reads 
int getLociToRawCounts( Loci & loci);

#endif