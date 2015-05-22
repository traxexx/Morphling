#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "SamFile.h"
#include "SamFileHeader.h"
#include "MeiCoord.h"
#include "QC.h"

// take ref-chr out + print disc sam
void PreProcessBam( string & inBam, string & outSam, string & disc_name, string & mei_coord_list);

// generate level list from disc.sam
void GenerateLevelListFromDiscBam( string & disc_name, string & list_prefix );

	 
/************ inner functions **************/

// non-ref sections
void processNonRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader, QC & BamQC, MeiCoord & mei_coord );


// ref section (print all to samOut)
void processRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader, QC & BamQC, MeiCoord & mei_coord );

#endif
