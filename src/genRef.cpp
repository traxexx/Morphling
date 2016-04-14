#include "genRef.h"
#include "Aligner.h"
#include "Globals.h"
#include "morphError.h"
#include "bamUtility.h"
#include "seqUtility.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>

using std::endl;
using std::map;

genRef::genRef( string & ref_bed_name, SingleSampleInfo & si, 
	vector<string> & MEseq, vector<string> & subMEname, bool neg_ctrl )
{
	// initialize
	bam_name = si.bam_name;
	avr_read_length = si.avr_read_length;
	max_ins_size = si.avr_ins_size + 3 * si.var_avr_ins_size;
	min_frac = 0.1;
	psMEseq = &MEseq;

	map<int, string> raw_ref_coord; // position->subtype-name
	loadRefBed(raw_ref_coord, ref_bed_name, neg_ctrl);
	map<int, int> ref_coord;
	// match index
	for( map<int, string>::iterator t = raw_ref_coord.begin(); t != raw_ref_coord.end(); t++ ) {
		if (neg_ctrl) { // neg then use 1st
			ref_coord[t->first] = 0;
			continue;
		}
		vector<string>::iterator vp = find( subMEname.begin(), subMEname.end(), t->second );
		if (vp == subMEname.end()) {
			string str = "[genRef::genRef] cannot find subtype: " + t->second;
			morphError(str, 32);
		}
		int idx = vp - subMEname.begin();
		ref_coord[t->first] = idx;
	}

	int rsize = ref_coord.size();
	rawStat.resize( rsize );
	for(int i=0; i<rsize; i++)
		rawStat[i].resize(NFIELD, 0);

	SamFile bam;
	SamFileHeader bam_header;
	SamFile alt_sam;
	SamFileHeader alt_sam_header;
	OpenBamAndBai( bam, bam_header, bam_name );
	OpenBamAndBai( alt_sam, alt_sam_header, bam_name);
	// read raw
	int n = 0; // trace in rawStat
	for( map<int, int>::iterator t = ref_coord.begin(); t!= ref_coord.end(); t++ ) {
		int center = t->first;
		int st = center - WIN/2 - avr_read_length/2;
		int ed = center + WIN/2 + avr_read_length/2;
		vector<int> stat;
		bool both_strand = 0;
		if (!neg_ctrl)
			both_strand = 1;
		bool stat_set = SetWindowStatFromBam( stat, avr_read_length, max_ins_size,
			(*psMEseq)[t->second], REF_CHR, st, ed, bam, bam_header, alt_sam, alt_sam_header, both_strand );
		if (!stat_set)
			continue;
		float sum = 0;
		int nzero = 0;
		for(int i=0; i<stat.size(); i++) {
			if (stat[i]==0)
				nzero++;
			sum += stat[i];
		}
		if (sum==0) {
			string str = "[genRef::MakeRefStat]: at n=" + std::to_string(n) + ", sum is zero";
			morphError(str, 37);
		}
		sum += nzero * min_frac;  // sometimes disaster to low-coverage
		for(int i=0; i<stat.size(); i++) {
			if (stat[i]>0)
				rawStat[n][i] = (float)stat[i] / sum;
			else
				rawStat[n][i] = min_frac / sum; // here convert 0 to a small value 1/ depth/3
//std::cout << "n=" << n << ";i=" << i << "stat[i]=" << stat[i] << ",sum=" << sum << "rawStat[n][i]=" << rawStat[n][i] << std::endl;	
		}
		n++;
//string str = "n=" + std::to_string(n);
//morphMessage(str);		
	}
	bam.Close();
	alt_sam.Close();

	string msg;
	if (neg_ctrl)
		msg = "   #Neg";
	else
		msg = "   #Pos";
	msg += " ctrl window =  " + std::to_string(n);
	morphMessageNoTime(msg);
	rawStat.resize(n);
	transferRawStatToRefStat();
	msg = "     #De-duped ctrl window = " + std::to_string(refStat.frequency.size());
	morphMessageNoTime( msg );
}


void genRef::transferRawStatToRefStat()
{
	// sort & uniq -c
	std::sort(rawStat.begin(), rawStat.end(), sortRawStat);

	// calculate dup
	int d = 0; // index of current dup
	vector<int> dup;
	dup.resize(1,1);
	for(int i=1; i<rawStat.size(); i++) {
		bool efv = equalFloatVector( rawStat[i], rawStat[i-1] );
		if (efv)
			dup[d]++;
		else {
			d++;
			dup.push_back(1);
		}
	}
	// convert to refStat
	refStat.stat.resize(dup.size());
	refStat.frequency.resize(dup.size());
	int rsi = 0; // rawStat index
	for(int i=0; i<dup.size(); i++) {
		transferSingleRawStatToRefStat( i, rawStat[rsi] );
		refStat.frequency[i] = dup[i];
		rsi += dup[i];
	}
}

void genRef::transferSingleRawStatToRefStat(int index, vector<float> & raw)
{
	refStat.stat[index].resize(NFIELD);
	for(int i=0; i<NFIELD; i++)
		refStat.stat[index][i] = log10(raw[i]);
}

// self defined sort funtion for rawStat
// compare from 1st to last
bool genRef::sortRawStat( const vector<float> & x, const vector<float> & y )
{
	if ( x.size() != y.size() ) {
		string str = "[genRef::sortRawStat] x.size and y.size, " + std::to_string(x.size()) + "," + std::to_string(y.size()) + " are unequal!";
		morphError(str, 60);
	}
	for(size_t i=0; i<x.size(); i++) {
		if (x[i] < y[i])
			return 1;
		else if (x[i] > y[i])
			return 0;
	}
	return 0;
}

/*
	format: frequency stat(log10)
*/
void genRef::Print( string & out_name )
{
	if (refStat.stat.size() != refStat.frequency.size()) {
		string str = "[genRef::Print] refStat frequency size=" + std::to_string(refStat.frequency.size());
		str += ", stat size=" + std::to_string(refStat.stat.size());
		morphError(str, 61);
	}

	std::ofstream output;
	output.open(out_name.c_str());
	for(int i=0; i<refStat.stat.size(); i++) {
		output << refStat.frequency[i];
		for(int j=0; j<NFIELD; j++)
			output << "\t" << refStat.stat[i][j];
		output << endl;
	}
	output.close();
}


// if het_stat is empty, then put all rawStat into that
// if het_stat not empty, then take size n*times, put n rawStat into that and average
void genRef::Export( vector< vector<float> > & het_stat )
{
	// hom
	if (het_stat.empty()) {
		// check if empty exists
		int empty = 0;
		for(int i=0; i<rawStat.size(); i++)
			if (rawStat[i].empty())
				empty++;
		if (empty == 0) { // directly assign
			het_stat = rawStat;
		}
		else { // eleminate empty
			het_stat.resize( rawStat.size() - empty );
			int hi = 0;
			for(int ri=0; ri<rawStat.size(); ri++) {
				if (rawStat[ri].empty())
					continue;
				het_stat[hi] = rawStat[ri];
				hi++;
			}
		}
		return;
	}

	// neg: add up and transform
	int times = 8; // #neg match to each hom window
	int n = het_stat.size();
	int expn = n*times;
	if (rawStat.size() < expn) {
		string str = "[genRef::Export] #ctrl-window are not enough to make het ref! n=";
		str += std::to_string(n) + ",rawStat size=" + std::to_string(rawStat.size());
		morphError(str, 63);
	}
	het_stat.resize(expn); // expand & copy
	for(int i=1; i<times; i++) {
		for(int j=0; j<n; j++)
			het_stat[i*n+j] = het_stat[j];
	}
	int v = rawStat.size() / expn;
	int vi = 0;
	for(int i=0; i<expn; i++) {
		if (rawStat[vi].size()==0) { // if zero, empty it
			het_stat[i].clear();
			continue;
		}
		if (het_stat[i].size() == 0) {
			string str = "at" + std::to_string(i) + ", het size = 0";
			morphError(str, 31);
		}
		for(int j=0; j<NFIELD; j++ ) {
			het_stat[i][j] += rawStat[vi][j];
			het_stat[i][j] /= 2; // average
			het_stat[i][j] = log10(het_stat[i][j]); // then transform to log
		}
		vi += v;
	}
}

bool equalFloatVector( vector<float> & f1, vector<float> & f2 )
{
	if (f1.size() != f2.size())
		morphError("[equalFloatVector] f1 & f2 not the same size");
	for(int i=0; i<f1.size(); i++) {
		if (f1[i] != f2[i])
			return 0;
	}
	return 1;
}

// position -> subtype-name
// neg ctrl only has 3 columns. Use 1st subtype
// pos ctrl has 4 columns
void loadRefBed( map<int, string> & ref_coord, string & ref_bed_name, bool neg_ctrl)
{
	std::ifstream bedfile;
	bedfile.open( ref_bed_name.c_str() );
	if ( !bedfile.is_open() )
		morphErrorFile( ref_bed_name );
	string line;
	while( std::getline(bedfile, line) ) {
		std::stringstream ss;
		ss << line;
		string field;
		std::getline( ss, field, '\t' );
		std::getline( ss, field, '\t' );
		int pos = stoi(field);
		std::getline( ss, field, '\t' );
		if (neg_ctrl)
			field = "";
		else
			std::getline( ss, field, '\t' );
		ref_coord[pos] = field;
	}
	bedfile.close();
}

/*
// if no stat set, return false
bool SetWindowStatFromBam( vector<int> & stat, int avr_read_length, float max_ins_size,
	string & mei_seq, string & chr, int st, int ed, SamFile & bam, SamFileHeader & bam_header,
	SamFile & alt_sam, SamFileHeader & alt_sam_header, bool both_strand )
{
	bool section_status = bam.SetReadSection( chr.c_str(), st, ed, 0 ); // set overlap = false
	if (!section_status)
		return 0;

	stat.resize(NFIELD, 0);

	string rv_mei;
	if (both_strand)
		rv_mei = RevCompSeq(mei_seq);
	SamRecord rec;
//	int lc = 0;
	while( bam.ReadRecord( bam_header, rec ) ) {
		// QC
		if (IsSupplementary(rec.getFlag()))
			continue;
		if (rec.getReadLength() < avr_read_length / 3)
			continue;
		if (rec.getFlag() & 0x4)
			continue;
		if (rec.getFlag() & 0x2) { // only check 1 read of the proper pair
			if (rec.get1BasedPosition() >= rec.get1BasedMatePosition())
				continue;
		}
		else { // skip translocate reads
			string mate_ref = rec.get1MateReferenceName();
			if (mate_ref.compare(rec.getReferenceName())==0) {
				if (abs(rec.getInsertSize())<max_ins_size )
					continue;
			}
		}
		if (rec.getMapQuality() < MIN_QUALITY)
			continue;
//		if ( rec.get1BasedPosition() - st < rec.getReadLength() / 2 )
//			continue;
//		if ( ed - rec.get1BasedAlignmentEnd() < rec.getReadLength() / 2 )
//			continue;

		// get index
		SamRecord mate_rec;
		bool mate_set = SetDiscMateRec(mate_rec, rec, alt_sam, alt_sam_header);
		if (!mate_set)
			continue;
		int index = getSingleTypeIndex( rec, mate_rec, mei_seq ); // rec is definitely the anchor
		if (both_strand && index>=9) { // ref coord doesn't have strand preference
			int index2 = getSingleTypeIndex( rec, mate_rec, rv_mei );
			if (index < 9)
				index = index2;
		}
		stat[index]++;
//		lc++;
	}
//string str = "lc=" + std::to_string(lc);
//morphMessage(str);

	int sum = 0;
	for(int i=0; i<stat.size(); i++)
		sum += stat[i];
	if (sum == 0)
		return 0;

	return 1;
}
*/

// if no stat set, return false
// for proper reads, set a vector to store reads
// for disc reads, find its mate by looking at corresponding section
// count polyA/T disc read as supporting read
bool SetWindowStatFromBam( vector<int> & stat, int avr_read_length, float max_ins_size,
	string & mei_seq, string & chr, int st, int ed, SamFile & bam, SamFileHeader & bam_header,
	SamFile & alt_sam, SamFileHeader & alt_sam_header, bool both_strand )
{
	SamRecord rec;
	/* check read count first
	bool section_status = bam.SetReadSection( chr.c_str(), st-avr_read_length, ed+avr_read_length, 0 ); // set overlap = false
	if (!section_status)
		return 0;
	int lc = 0;
	int max_lc = GetMaxIntervalReadCount( sample_depth[sample_depth.size()-1], avr_read_length, ed - st );
	while( bam.ReadRecord( bam_header, rec ) ) {
		lc++;
		if (lc >= max_lc)
			return 0;
	}
	*/

	// do stat now
	bool section_status = bam.SetReadSection( chr.c_str(), st-avr_read_length, ed+avr_read_length, 0 ); // set overlap = false
	if (!section_status)
		return 0;
	stat.resize(NFIELD, 0);

	struct mate_record {
		int flag;
		string seq;
		string cigar;
	};

	// set vector to store proper reads
	vector< map<string, mate_record> > rstore;
	rstore.resize(ed-st+1); // do not need the last one

	string rv_mei;
	if (both_strand)
		rv_mei = RevCompSeq(mei_seq);
//	int lc = 0;
	while( bam.ReadRecord( bam_header, rec ) ) {
		// QC
		if (!qcCheck(rec.getFlag()))
			continue;
		if (IsSupplementary(rec.getFlag()))
			continue;
		if (rec.getReadLength() < avr_read_length / 3)
			continue;
		if (rec.getFlag() & 0x4)
			continue;
		if (rec.getMapQuality() < MIN_QUALITY)
			continue;
		SamRecord mate_rec;
		if (rec.getFlag() & 0x2) { // proper
			if (rec.get1BasedPosition() < st || rec.get1BasedAlignmentEnd() > ed)
				continue;
			if (rec.get1BasedPosition() < rec.get1BasedMatePosition()) { // add to rstore
				int idx = rec.get1BasedPosition() - st;
//std::cout << rec.get1BasedPosition() << " " << st << " " << idx << std::endl;
				if (idx < 0)
					continue;
				rstore[idx][rec.getReadName()].flag = rec.getFlag();
				rstore[idx][rec.getReadName()].seq = rec.getSequence();
				rstore[idx][rec.getReadName()].cigar = rec.getCigar();
				continue;
			}
			else { // check for mate
				int idx = rec.get1BasedMatePosition() - st;
				if (idx < 0)
					continue;
				if (rstore[idx].empty())
					continue;
				map<string, mate_record>::iterator t = rstore[idx].find(rec.getReadName());
				if (t==rstore[idx].end())
					continue;
				mate_rec.setFlag( t->second.flag );
				mate_rec.setCigar( t->second.cigar.c_str() );
				mate_rec.setSequence( t->second.seq.c_str() );
				rstore[idx].erase( rec.getReadName() );
			}

		}
		else { // disc add directly
			// skip translocate
			string mate_ref = rec.getMateReferenceName();
			if (mate_ref.compare(rec.getReferenceName())==0) {
				if (abs(rec.getInsertSize())<max_ins_size )
					continue;
			}
			// get mate
			bool mate_set = SetDiscMateRec(mate_rec, rec, alt_sam, alt_sam_header);
			if (!mate_set)
				continue;
		}

		// get index
		int index = getSingleTypeIndex( rec, mate_rec, mei_seq ); // rec is definitely the anchor
		if (both_strand && index>=9) { // ref coord doesn't have strand preference
			int index2 = getSingleTypeIndex( rec, mate_rec, rv_mei );
			if (index < 9)
				index = index2;
		}
		stat[index]++;
//		lc++;
	}
//string str = "lc=" + std::to_string(lc);
//morphMessage(str);

	int sum = 0;
	for(int i=0; i<stat.size(); i++)
		sum += stat[i];
	if (sum == 0)
		return 0;

	return 1;
}



/*
	0  	proper
	1 	short-clip

	if not mei, then + 8;
	final index is <below> + 2;

	0	00	left read, begin clip
	1	01	left read, end clip
	2	10	right read, begin clip
	3	11	right read, end clip

	4	100 left anchor (disc)
	5	101	right anchor (disc)
	6	110 left anchor (unmap)
	7	111	right anchor (unmap)
*/
int getSingleTypeIndex( SamRecord & rec, SamRecord & mate_rec, string & mei_seq )
{
	int index;
	if (rec.getFlag() & 0x2) { // proper
		int c1 = GetMaxClipLen(rec);
		int c2 = GetMaxClipLen(mate_rec);
		int abs_cmax;
		int cmax;
		SamRecord * p_clip_rec;
		if ( abs(c1) >= abs(c2) ) {
			cmax = c1;
			abs_cmax = abs(c1);
			p_clip_rec = &rec;
		}
		else {
			cmax = c2;
			abs_cmax = abs(c2);
			p_clip_rec = &mate_rec;
		}
		if (abs_cmax < MIN_CLIP/2)
			return 0;
		if (abs_cmax < MIN_CLIP)
			return 1;
		// then index
		index = cmax > 0 ? 0 : 1;
		if (p_clip_rec->getFlag() & 0x10)
			index += 2;
		// see if it's mei
		string clip_seq = GetMaxClipSeq( *p_clip_rec );
		Aligner al( mei_seq, clip_seq );
		if ( !al.IsMapped() )
			index += 8;
	}
	else { // disc. Be careful about strand
		index = (rec.getFlag() & 0x10) ? 5 : 4;
		if (rec.getFlag() & 0x8)
			index += 2;
		string disc_seq = GetCorrectMateSeq( rec, mate_rec );
		Aligner al(mei_seq, disc_seq);
		if ( !al.IsMapped() ) {
/*		    // polyA status of the disc read
		    bool polyA = ContainRepetative( disc_seq, 'A');
		    if (!polyA)
		        polyA = ContainRepetative( disc_seq, 'T');
		    //polyA status of anchor
		    string anchor_seq;
		    if (!polyA) {
		        seq = rec.getSequence();
		        polyA = ContainRepetative( anchor_seq, 'A');
		    }
		    if (!polyA)
		        polyA = ContainRepetative( anchor_seq, 'T');
		    if (!polyA) */
                index += 8;   
		}
	}

	// finally return
	index += 2;
	return index;
}



