#include "bamUtility.h"
#include "seqUtility.h"
#include "morphError.h"

using std::string;

void OpenBamOnly( SamFile & bam, SamFileHeader & bam_header, string & bam_name )
{
	// bam status
	bam.OpenForRead( bam_name.c_str() );
	if ( !bam.IsOpen() )
		morphErrorFile(bam_name);

	// header status
	bool header_status = bam.ReadHeader( bam_header );
	if (!header_status) {
		string str = "fail to read bam header: " + bam_name;
		morphError(str, 11);
	}
}

void OpenBamAndBai( SamFile & bam, SamFileHeader & bam_header, string & bam_name )
{
	OpenBamOnly( bam, bam_header, bam_name );

	// bai status
	string bai = bam_name + ".bai";
	bool bai_status = bam.ReadBamIndex( bai.c_str() );
	if (!bai_status) {
		bai = bam_name.substr(0, bam_name.length()-1) + 'i';
		bai_status = bam.ReadBamIndex( bai.c_str() );
	}
	if (!bai_status)
		morphErrorFile( bai );
}


// theoretically this can find mate of any reads
//	but proper read will be recorded in a vector so don't need to use this
//	use this to find proper reads will dramatically slow down
bool SetDiscMateRec( SamRecord & mate_rec, SamRecord & rec, SamFile & alt_sam, SamFileHeader & alt_sam_header )
{
	bool found = 0;

	string chr = rec.getMateReferenceName();
	int st = rec.get1BasedMatePosition();

	bool section_status = alt_sam.SetReadSection( chr.c_str(), st-1, st+1.5*rec.getReadLength(), 0 ); // overlap=false. only browse reads that completely fall within the region
	if (!section_status)
		return 0;
	while( alt_sam.ReadRecord(alt_sam_header, mate_rec) ) {
		if ( mate_rec.get1BasedPosition() > st )
			break;
//		if ( mate_rec.get1BasedPosition() != st )
//			continue;
//		if ( mate_rec.get1BasedMatePosition() != rec.get1BasedPosition() )
//			continue;
//		found = 1;
//		break;
		string rn = mate_rec.getReadName();
		if ( rn.compare(rec.getReadName()) == 0 ) {
			found = 1;
			break;
		}
	}

// return flag
	if (found)
		return 1;
	else
		return 0;
}

/*

// deprecated because all info can be got from sample list

int GetAvrInsSize( string & bam_name, string & chr )
{
	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam, bam_header, bam_name );
	bool section_status = bam.SetReadSection( chr.c_str() );
	if (!section_status) {
		string str = "[bamUtility: GetAvrInsSize] Unable to set read section to chr " + chr;
		morphError(str, 30);
	}

	int ins_size = 0;
	int counter = 0;
	int lc = 0;
	SamRecord rec;
	while( bam.ReadRecord( bam_header, rec ) ) {
		lc++;
		if (lc < 20)
			continue;
		if ( !(rec.getFlag() & 0x2) ) // skip disc
			continue;
		if ( IsSupplementary( rec.getFlag() ) )
			continue;
		if (counter > 100)
			break;
		counter++;
		ins_size += abs( rec.getInsertSize() );
	}
	bam.Close();

	ins_size /= counter;
	return ins_size;
}


int GetAvrReadLength( string & bam_name, string & chr )
{
	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam, bam_header, bam_name );
	bool section_status = bam.SetReadSection( chr.c_str());
	if (!section_status) {
		string str = "[bamUtility: GetAvrInsSize] Unable to set read section to chr " + chr;
		morphError(str, 30);
	}

	int avr_len = 0;
	int counter = 0;
	int lc = 0;
	SamRecord rec;
	while( bam.ReadRecord( bam_header, rec ) ) {
		lc++;
		if (lc<20)
			continue;
		if ( !(rec.getFlag() & 0x2) ) // skip disc
			continue;
		if ( IsSupplementary( rec.getFlag() ) )
			continue;
		if (counter > 100)
			break;
		counter++;
		avr_len += rec.getReadLength();
	}
	bam.Close();

	avr_len /= counter;
	return avr_len;
}
*/

bool IsSupplementary( int flag )
{
	if ( flag & 0x800 )
		return 1;
	return 0;
}

bool qcCheck( int flag )
{
	if (flag & 0x400)
		return 0;
	if (flag & 0x200)
		return 0;
	if (flag & 0x100)
		return 0;
	return 1;
}

bool OneEndUnmap( int flag )
{
	if (flag & 0x8)
		return 1;
	if (flag & 0x4)
		return 1;
	return 0;
}

// be aware of strand
string GetCorrectMateSeq( SamRecord & rec, SamRecord & mate_rec )
{
	string seq = std::string( mate_rec.getSequence() );
	if (rec.getFlag() & 0x10) { // anchor reverse comp
		if (mate_rec.getFlag() & 0x10) // get it back if automatically rev comped by the mapper
			seq = RevCompSeq(seq);
	}
	else { // anchor plus strand
		if ( !(mate_rec.getFlag() & 0x10) ) // if not reversed then reverse
			seq = RevCompSeq(seq);
	}
	return seq;
}

// get max # read count in window
// if > that, do not perform assembly or genotype
int GetMaxIntervalReadCount( float current_depth, int avr_read_length, int interval )
{
	int max_lc = current_depth / (float)avr_read_length * (float)(interval) * 10;
	return max_lc;
}



