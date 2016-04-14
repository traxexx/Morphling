#include "Likelihood.h"
#include "morphError.h"
#include "genRef.h"
#include "qcCheck.h"
#include "bamUtility.h"
#include "Globals.h"
#include <math.h>
#include <algorithm>
#include <fstream>

Likelihood::Likelihood( vector< vector<string> > & MEseq, SingleSampleInfo & si, string & ref_file_prefix )
{
	string bam_name = si.bam_name;
	avr_read_length = si.avr_read_length;
	max_ins_size = si.avr_ins_size + 3 * si.var_avr_ins_size;
	refStat.resize(NMEI);
	nref.resize(NMEI);
	log10_total.resize(NMEI);
	for(int i=0; i<NMEI; i++) {
		refStat[i].resize(3);
		nref[i].resize(3, -1);
		log10_total[i].resize(3);
	}
	vector<string> suffix_name;
	suffix_name.resize(3);
	suffix_name[2] = ".hom-stat.";
	suffix_name[1] = ".het-stat.";
	suffix_name[0] = ".neg-stat.";
	for (int m=0; m<3; m++) {
		for(int g=0; g<3; g++) {
			string ref_file_name = ref_file_prefix + suffix_name[g] + std::to_string(m);
			// line count first
			nref[m][g] = GetFileLineCount(ref_file_name);
			if (nref[m][g] <= 0) {
				string str = "[Likelihood::Likelihood] " + ref_file_name + "is empty";
				morphError(str, 35);
			}
			// read
			setFromRefFile( m, g, ref_file_name );
			// frequency to log
			int sum = 0;
			for(int i=0; i<refStat[m][g].frequency.size(); i++) {
				sum += refStat[m][g].frequency[i];
				refStat[m][g].frequency[i] = log10(refStat[m][g].frequency[i]);
			}
			log10_total[m][g] = log10(sum);
		}
	}

	// open bam
	OpenBamAndBai( bam, bam_header, bam_name );
	OpenBamAndBai( alt_bam, alt_bam_header, bam_name );
	pMEseq = &MEseq;
}


void Likelihood::SetLhRec( lhRec & lr, string chr, int pos, int mtype, int subtype )
{
	vector<int> stat;
	bool status = setWindowStat( stat, chr, pos, (*pMEseq)[mtype][subtype] );
	if (!status) { // NA. clear everything
		lr.gl.clear();
		lr.ad.clear();
		return;
	}
	// fill in gl & ad
	lr.gl.resize(3);
	lr.ad.resize(3,0);
	for(int i=0; i<3; i++) {
		lr.gl[i] = getLog10Likelihood(mtype, i, stat);
	}
	lr.ad[0] = getFlankingCount( chr, pos ); // #ref read
	for(int i=2; i<10; i++) // #alt read
		lr.ad[1] += stat[i];
	for(int i=10; i<NFIELD; i++) // #ambiguous read
		lr.ad[2] += stat[i];
	// normalize gl
//	int max_index = 0;
	float max_val = lr.gl[0];
	for(int i=1; i<3; i++) {
		if (lr.gl[i] > max_val) {
			max_val = lr.gl[i];
//			max_index = i;
		}
	}
	for(int i=0; i<3; i++)
		lr.gl[i] -= max_val;
}


float Likelihood::getLog10Likelihood( int m, int g, vector<int> & wc )
{
	float sum = 0;
	// initialize
	for(int i=0; i<NFIELD; i++)
		sum += refStat[m][g].stat[0][i] * wc[i];
	if (refStat[m][g].stat.size() == 1) {
		sum -= log10_total[m][g];
		return sum;
	}
	// rest
	for(int i=1; i<refStat[m][g].stat.size(); i++) {
		float lh = 0;
		for(int j=0; j<NFIELD; j++)
			lh += refStat[m][g].stat[i][j] * wc[j];
		lh += refStat[m][g].frequency[i];
		if (sum - lh > 3)
			continue;
		else if (sum - lh < -3)
			sum = lh;
		else
			sum = SumTwoLog10s( sum, lh );
	}
	// avreage by ref window
	sum -= log10_total[m][g];
	return sum;
}


float SumTwoLog10s( float l1, float l2 )
{
	float mid = (l1 + l2) / 2;
	float sum = mid + log10( pow(10, l1-mid) + pow(10, l2-mid) );
	return sum;
}

// format: frequency stat
void Likelihood::setFromRefFile( int m, int g, string & ref_file_name )
{
	refStat[m][g].stat.resize(nref[m][g]);
	refStat[m][g].frequency.resize(nref[m][g]);

	std::ifstream ref_file;
	ref_file.open( ref_file_name.c_str() );
	string line;
	int idx = 0;
	while( std::getline(ref_file, line) ) {
		std::stringstream ss;
		ss << line;
		string field;
		std::getline(ss, field, '\t');
		refStat[m][g].frequency[idx] = stoi(field); // frequency
		refStat[m][g].stat[idx].resize(NFIELD);
		for(int i=0; i<NFIELD; i++) {
			std::getline(ss, field, '\t');
			refStat[m][g].stat[idx][i] = stof(field);
		}
		idx++;
	}
	ref_file.close();
	if (idx != nref[m][g])
		morphError("[Likelihood::setFromRefFile] idx != nref[g]");
}

bool Likelihood::setWindowStat( vector<int> & stat, string & chr, int pos, string & mseq )
{
	int st = pos - WIN/2;
	int ed = pos + WIN/2;
	bool status = SetWindowStatFromBam( stat, avr_read_length, max_ins_size,
		mseq, chr, st, ed, bam, bam_header, alt_bam, alt_bam_header, 0 );
	return status;
}

int Likelihood::getFlankingCount( string & chr, int pos )
{
	bool status = bam.SetReadSection( chr.c_str(), pos - avr_read_length * 2, pos + avr_read_length/2 );
	if (!status)
		return 0;

	int n = 0;
	SamRecord rec;
	while(bam.ReadRecord(bam_header, rec)) {
		if ( !(rec.getFlag() & 0x2) )
			continue;
		if (IsSupplementary(rec.getFlag()))
			continue;
		if (rec.getReadLength() < avr_read_length / 3)
			continue;
		if (rec.getMapQuality() < MIN_QUALITY)
			continue;
		if (rec.get1BasedPosition() + MIN_CLIP/2 <= pos && rec.get1BasedAlignmentEnd() - MIN_CLIP/2 >= pos )
			n++;
	}

	return n;
}

// print gt info in format of GT:PL:AD:OD
void PrintLhRecToVcf( std::ofstream & vcf, int gt, lhRec & lr )
{
	// gt
	if (gt == -1)
		vcf << "./.";
	else if (gt==0)
		vcf << "0/0";
	else if (gt==1)
		vcf << "0/1";
	else if (gt==2)
		vcf << "1/1";
	else
		morphError("[PrintLhRecToVcf] genotype unrecognized");
	vcf << ":";
	// gl
	vector<float> pl;
	pl.resize(3);
	for(int i=0; i<3; i++) {
		pl[i] = round(-lr.gl[i] * 10);
		if (pl[i] <0.0001)
			pl[i] = 0;
	}

	vcf << pl[0] << "," << pl[1] << "," << pl[2];
	vcf << ":";
	// ad & od
	vcf << lr.ad[0] << "," << lr.ad[1] << ":" << lr.ad[2];
}





































