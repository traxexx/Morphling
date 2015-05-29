#include "bamUtility.h"
#include "StringBasics.h"
#include "Utilities.h"
#include <algorithm>    // std::min_element, std::max_element
#include <iostream>

using std::cerr;
using std::endl;

bool ExistChrInBam( SamFileHeader & samHeader, string & focus_chr )
{
	bool focus_exist = 0;
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String chr_name = samHeader.getReferenceLabel(i);
		if ( focus_chr.compare(chr_name.c_str()) == 0 ) {
			focus_exist = 1;
			break;
		}
	}
	return focus_exist;
}


// sort bam by name, return name: (base)disc_name-nsort.bam
string SortBamByName( string & disc_name )
{
	string nsort_name = disc_name.substr( 0, disc_name.length()-4 ) + "-nsort";
	string cmd = string("samtools sort -n ") + disc_name + " " + nsort_name;
	ExecuteCmd(cmd);
	nsort_name += ".bam";
	
	return nsort_name;
}

// set avr read len based on 1st 201~300st reads
int GetAvrReadLenFromBam( string & bam )
{		
// open bam
	SamFile samIn;
	SamFileHeader samHeader;
	OpenBamAndBai( samIn, samHeader, bam );

// let's process....
	int Counter = 0;
	SamRecord sam_rec;
	int avr_read_len = 0;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
		Counter++;
		if (Counter > 300)
			break;
		else if (Counter < 200)
			continue;
		avr_read_len += sam_rec.getReadLength();
	}
	samIn.Close();

// get avr, return	
	avr_read_len /= (Counter - 200);
	samIn.Close();
	return avr_read_len;
}


// avr read pair insert size: based on 1st 201~300 0x2 reads
int GetAvrInsSizeFromBam( string & bam )
{
// open bam
	SamFile samIn;
	SamFileHeader samHeader;
	OpenBamAndBai( samIn, samHeader, bam );
	
// process
	int Counter = 0;
	SamRecord sam_rec;
	int avr_ins_size = 0;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
		if ( !(sam_rec.getFlag() & 0x2) ) // skip disc
			continue;
		Counter++;
		if (Counter > 300)
			break;
		else if (Counter < 200)
			continue;
		avr_ins_size += abs( sam_rec.getInsertSize() );
	}
	samIn.Close();

// get avr, return	
	avr_ins_size /= (Counter - 200);
	samIn.Close();
	return avr_ins_size;
}

// get 100 regions in chr20 then calculate depth
float EstimateBamDepth( string & bam, int avr_read_len, string & ref_chr )
{
// parameters
	int pad = 200000;
	int cwin = 2000;
	int nwin = 60;
	int minlen = cwin * nwin;

// open bam
	SamFile samIn;
	SamFileHeader samHeader;
	OpenBamAndBai( samIn, samHeader, bam );
// get ref chr length & set regions
	int reflen = atoi( samHeader.getSQTagValue("LN", ref_chr.c_str()) );
	int trimlen = reflen - pad * 2;
	if ( trimlen < minlen ) {
		cerr << "ERROR: ref-chr: " << ref_chr << " is too short to estimate read depth. Are you using the right ref chromosome?" << endl;
		exit(1);
	}
	
// count read number
	int times = trimlen / minlen;
	vector<int> vlc;
	vlc.resize( nwin );
	for( int seg = 0; seg < nwin; seg++ ) {
		bool section_status = samIn.SetReadSection(ref_chr.c_str(), pad + seg*times, pad + seg*times + cwin);
		if (!section_status) {
			std::cerr << "ERROR: Unable to set read section: " << ref_chr << ": " << pad + seg*times << "-" << pad + seg*times + cwin << std::endl;
			exit(1);
		}
		int lc = 0;
		SamRecord sam_rec;
		while( samIn.ReadRecord(samHeader, sam_rec) ) {
			lc++;
		}
		vlc[ seg ] = lc;
	}
	samIn.Close();
	int rcount = 0;
	for( vector<int>::iterator it = vlc.begin(); it != vlc.end(); it++ ) {
		rcount += *it;
	}
	// remove max & min
	rcount -= ( *std::min_element( vlc.begin(), vlc.end() ) + *std::max_element( vlc.begin(), vlc.end() ) );
	
// depth
	float depth;
	int region_size = minlen - 2 * cwin - 2 * avr_read_len * nwin;
	depth = float(rcount) * avr_read_len / region_size;
//	depth = float(rcount) / (region_size / WIN) / 2;
	return depth;
}


// check in & out & bai
void SanityCheckBams( SamFile & samIn, SamFile & samOut, bool & bai_status )
{
	if ( !samIn.IsOpen() ) {
		cerr << "ERROR: Unable to open input bam " << endl;
		exit(1);
	}
	
	if ( !samOut.IsOpen() ) {
		cerr << "ERROR: Unable to open output bam " << endl;
		exit(1);
	}
	
	if ( !bai_status ) {
		cerr << "ERROR: Unable to read bam index, is it indexed?" << endl;
		exit(1);
	}
}

void OpenBamAndBai( SamFile & samIn, SamFileHeader & samHeader, string & bam_name )
{
// file check
	samIn.OpenForRead( bam_name.c_str() );
	if ( !samIn.IsOpen() ) {
		cerr << "ERROR: can't open bam " << bam_name << endl;
		exit(1);
	}
	bool header_status = samIn.ReadHeader(samHeader);
	if (!header_status) {
		cerr << "ERROR: Fail to read header: " << bam_name << endl;
		exit(1);
	}
// bai status
	string bai = bam_name + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai.c_str() );
	if ( !bai_status ) {
		cerr << "ERROR: can't read bam index: " << bai << endl;
		exit(1);
	}
}

void SetChrListFromBamHeader( vector<string> & chr_list, string & bam_name )
{
// open & sanity check
	SamFile samIn;
	samIn.OpenForRead( bam_name.c_str() );
	if ( !samIn.IsOpen() ) {
		cerr << "ERROR: can't open bam " << bam_name << endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool header_status = samIn.ReadHeader(samHeader);
	if (!header_status) {
		cerr << "ERROR: Fail to read header: " << bam_name << endl;
		exit(1);
	}
// add to chr_list
	chr_list.clear();
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String chr_name_str = samHeader.getReferenceLabel(i);
		string chr_name = string(chr_name_str.c_str());
		chr_list.push_back( chr_name );
	}	
	samIn.Close();
}


