#ifndef BAMUTILITY_H
#define BAMUTILITY_H

#include <string>
#include "SamFile.h"
#include "SamFileHeader.h"

using std::string;

// check if chr exist in bam by checking header
bool ExistChrInBam( SamFileHeader & samHeader, string & focus_chr );

// sort bam by name, return name: (base)disc_name-nsort.bam
string SortBamByName( string & disc_name );

int GetAvrReadLenFromBam( string & bam);

int GetAvrInsSizeFromBam( string & bam );

float EstimateBamDepth( string & bam, int avr_read_len, string & ref_chr );

void SanityCheckBams( SamFile & samIn, SamFile & samOut, bool & bai_status );

void OpenBamAndBai( SamFile & samIn, SamFileHeader & samHeader, string & bam_name );

#endif

