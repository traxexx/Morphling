#include "StringArray.h"
#include "StringHash.h"
#include "Parameters.h"
#include "Error.h"

#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <fstream>

#include "Preprocess.h"
#include "bamUtility.h"
#include "Utilities.h"
#include "Globals.h"
#include "ReadMapGlobals.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;

/**
	open raw bam & out-disc-bam & out-ref-bam
	initialize coord MEI bed
	process bam by chr & print
**/
void PreProcessBam( string & inBam, string & outSam, string & disc_name, string & mei_coord_list)
{
	SamFile samIn, samOut;
	SamFileHeader samInHeader;
	OpenBamAndBai( samIn, samInHeader, inBam );
	samOut.OpenForWrite( outSam.c_str() );	
	if ( !samOut.IsOpen() ) {
		cerr << "Unable to open " << outSam << endl;
	}
	samOut.WriteHeader(samInHeader);
	
	SamFile discSam;
	discSam.OpenForWrite( disc_name.c_str() );
	if ( !discSam.IsOpen() ) {
		cerr << "Unable to open " << disc_name << endl;
		exit(1);
	}
	discSam.WriteHeader( samInHeader);

// initialize mei-coord-ref
	MeiCoord mei_coord( mei_coord_list );


	QC BamQC;
// process bam by section
	cout << endl;
	cout << "Pre-process raw bam..." << endl;
	for(int i=0; i< samInHeader.getNumSQs(); i++) {
		String val = samInHeader.getReferenceLabel(i);
		string current_chr = string(val.c_str());
		if ( current_chr.length() > 5 && (!PSEUDO_CHR)) // skip pseudo
			continue;
		bool section_status = samIn.SetReadSection(current_chr.c_str());
		if (!section_status) {
			cerr << "ERROR: Unable to set read section to chr " << current_chr << endl;
			exit(1);
		}
		
		mei_coord.ResetVecPtr( current_chr );
		if ( current_chr.compare(REF_CHR) == 0) { // ref section
			cout << "Working on ctrl chr: " << current_chr << endl;
			processRefSection( samIn, samOut, discSam, samInHeader, BamQC, mei_coord );
		}
		else { // non ref
			cout << "Working on chr: " << current_chr << endl;
			processNonRefSection( samIn, samOut, discSam, samInHeader, BamQC, mei_coord );
		}
	}
	cout << "Finished! QC Summary of raw bam: " << endl;	
	BamQC.PrintQCsummary();
}


/* generate level list from disc.sam
   levels list format: 2 col
   	# evidence disc		# non-evidence ( disc + unmap )
  			0,1,2,3,4,5 etc
   ADD BY MATE (EM is for the current read)
   */
void GenerateLevelListFromDiscBam( string & disc_name, string & list_prefix )
{
	if (DEBUG_MODE)
		cout << "Generating level list from disc bam..." << endl;

	SamFile samIn;
	samIn.OpenForRead( disc_name.c_str() );
	if ( !samIn.IsOpen() ) {
		cerr << "ERROR: Can't open " << disc_name << endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool disc_header_status = samIn.ReadHeader(samHeader);
	if (!disc_header_status) {
		cerr << "ERROR: Fail to read header: " << disc_name << endl;
		exit(1);
	}
	
// generate list for all then print out by chr

  // generate chr map
  	map<string, vector< std::pair<int, int> > > levelMap;
  	for(int i=0; i< samHeader.getNumSQs(); i++) {
  		String chr_Str = samHeader.getReferenceLabel(i);
  		string chr_str = string(chr_Str.c_str());
  		int chr_size = std::stoi( samHeader.getSQTagValue("LN", chr_str.c_str()) );
		int win_count = chr_size / STEP + 1;
		std::pair<int, int> null_pair (-1, -1);
		levelMap[chr_str].resize(win_count, null_pair);
  	}
	
  // add level from bam
	SamRecord sam_rec;
	int line_count = 0;
	while(samIn.ReadRecord(samHeader, sam_rec)) {
		line_count++;
		if (sam_rec.getFlag() & 0x4) // skip unmap reads
			continue;
		string current_chr = sam_rec.getMateReferenceName();
		int current_index_min = (sam_rec.get1BasedMatePosition() + AvrReadLen - WIN) / STEP + 1;
		if (current_index_min < 0)
			current_index_min = 0;
		int current_index_max = sam_rec.get1BasedMatePosition() / STEP;
		if ( current_index_max < 0 ) {
			cerr << "Warning: negateive mapped position at: " << sam_rec.getReadName() << ". Skipped!" << endl;
			continue;
		}
		if ( current_index_max >= int(levelMap[current_chr].size()) ) {
			cerr << "Warning: Mate of : " << sam_rec.getReadName() << " at: " << current_chr << "-" << sam_rec.get1BasedMatePosition() << " out of chr boundary. Skipped!" << endl;
			continue;
		}
		int * tag = sam_rec.getIntegerTag("EM");
		bool em_type;
		if (tag)
			em_type = *tag > 0 ? 1 : 0;
		else // no EM (that's not error, including some disc & unmap)
			em_type = 0;
		vector< std::pair<int, int> >::iterator lmap = levelMap[current_chr].begin() + current_index_min;
		if (em_type > 0) { // add to evidence
			for(int i = current_index_min; i <= current_index_max; i++, lmap++)
				(lmap->first)++;
		}
		else { // add to non-evidence
			for(int i = current_index_min; i <= current_index_max; i++, lmap++)
				(lmap->second)++;
		}
	}


  // print level info
	for( map<string, vector< std::pair<int, int> > >::iterator mit = levelMap.begin(); mit != levelMap.end(); mit++ ) {
		ofstream list_file;
		string list_name = string(list_prefix) + "." + mit->first;
		if ( mit->first.length() > 5 && (!PSEUDO_CHR) ) { // skip pseudo chr
			continue;
		}
		list_file.open(list_name.c_str());
		CheckOutFileStatus(list_file, list_name.c_str());
		vector< std::pair<int, int> >::iterator it = mit->second.begin();
		for( unsigned int dist=0; dist < mit->second.size(); dist++, it++ ) {
			if (it->first > -1 || it->second > -1) {
				list_file << dist*STEP << "\t" << it->first << "\t" << it->second << endl;
			}
		}
		list_file.close();
	}
}

/*************** inner function ***************/

// only print: related to REF-CHR to samOut & disc to discSam
void processNonRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader, QC & BamQC, MeiCoord & mei_coord )
{
	int Counter = 0;
	SamRecord sam_rec;

	while(samIn.ReadRecord(samHeader, sam_rec)) {
		Counter++;

		if ( sam_rec.getFlag() & 0x2 )
			continue;
		bool qc_pass = BamQC.PassQC( sam_rec );
		if ( !qc_pass ) 
			continue;
		if ( sam_rec.getFlag() & 0x8 ) {
			if ( sam_rec.getFlag() & 0x4 )
				continue;
			discSam.WriteRecord(samHeader, sam_rec);
			continue;
		}
		string mate_chr = sam_rec.getMateReferenceName();
		if (mate_chr.compare(REF_CHR) == 0)
			samOut.WriteRecord(samHeader, sam_rec);
		if ( sam_rec.getFlag() & 0x4 ) {
			discSam.WriteRecord(samHeader, sam_rec);
			continue;
		}	
	
	// reset qual if irregular chr? should we do that?
		bool disc_pass = DiscSamPass ( sam_rec );
		if ( disc_pass ) {
			mei_coord.SetEMtag( sam_rec );
			discSam.WriteRecord(samHeader, sam_rec);
		}
	}
	cout << "  Total = " << Counter << " reads processed!" << endl;
}


// print all to samOut & disc to discSam
void processRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader, QC & BamQC, MeiCoord & mei_coord )
{
	int Counter = 0;
	SamRecord sam_rec;
	
	while(samIn.ReadRecord(samHeader, sam_rec)) {
		Counter++;
		if (sam_rec.getFlag() & 0x4) {
			if ( sam_rec.getFlag() & 0x8 )
				continue;
			discSam.WriteRecord(samHeader, sam_rec);
			samOut.WriteRecord(samHeader, sam_rec);
			continue;
		}
		samOut.WriteRecord(samHeader, sam_rec);
		
		if ( sam_rec.getFlag() & 0x2 ) 
			continue;
		bool qc_pass = BamQC.PassQC( sam_rec );
		if ( !qc_pass ) 
			continue;

		if ( sam_rec.getFlag() & 0x8 ) {
			discSam.WriteRecord(samHeader, sam_rec);
			continue;
		}	
	
	// reset qual if irregular chr? should we do that?
		bool disc_pass = DiscSamPass( sam_rec );
		if ( disc_pass ) {
			mei_coord.SetEMtag( sam_rec );
			discSam.WriteRecord(samHeader, sam_rec);
		}
	}
	cout << "  Total = " << Counter << " reads processed!" << endl;
}


















