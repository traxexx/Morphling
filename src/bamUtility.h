#ifndef BAMUTILITY_H
#define BAMUTILITY_H

#include "SamFile.h"
#include <string>

void OpenBamOnly( SamFile & bam, SamFileHeader & bam_header, std::string & bam_name );

void OpenBamAndBai( SamFile & bam, SamFileHeader & bam_header, std::string & bam_name );

bool SetDiscMateRec( SamRecord & mate_rec, SamRecord & rec, SamFile & alt_sam, SamFileHeader & alt_sam_header );

//int GetAvrInsSize( std::string & bam_name, std::string & chr );

//int GetAvrReadLength( std::string & bam_name, std::string & chr );

bool IsSupplementary( int flag );

bool qcCheck( int flag );

bool OneEndUnmap( int flag );

std::string GetCorrectMateSeq( SamRecord & rec, SamRecord & mate_rec );

int GetMaxIntervalReadCount( float current_depth, int avr_read_length, int interval );

#endif
