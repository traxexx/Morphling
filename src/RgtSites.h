#ifndef RGTSITES_H
#define RGTSITES_H

#include <string>
#include <map>
#include <vector>
#include "SamFile.h"
#include "SamRecord.h"

using std::string;
using std::map;
using std::vector;

struct GtInfo
{
	int ad_ref;
	int ad_alt;
	int ad_undef;
	float log_gl_ref;
	float log_gl_alt;
	float log_gl_het;
};

struct RgSiteInfo {
	int meitype;
	int subtype;
	int mstart;
	int mend;
	int left_end;
	int right_start;
	bool is_plus_strand;
	map<int, GtInfo> gt_info;
};

class RgtSites
{
public:
	RgtSites( vector< vector<string> > & mei_names, vector< vector<string> > & mei_seqs );
	void LoadSitesFromVcf( string & vcf_name );
	void ReGenotypeFromSingleBam( string & sample_name, string & bam_name );
	void PrintVcf( string & vcf_name );

private:
	void constructMeiNameMap( map<int, map<string, int> > & mmap );
	void addProperRecord( SamRecord & rec, int breakp, RgSiteInfo & rsi, string & mei_seq_left, string & mei_seq_right);
	void addDiscRecord( SamFile & sam, SamFileHeader & sam_header, SamRecord & rec, int breakp, RgSiteInfo & rsi, string & mei_seq_left, string & mei_seq_right);
	void setGinfoFromAltLogP( GtInfo & ginfo, float alt_logp );

	vector< vector<string> > * pMEseqs;
	vector< vector<string> > * pMEnames;

	map<string, map<int, RgSiteInfo> > siteList;
	map<string, int> sampleIndex;
	string stat_chr;
	int nsample;
	int avr_ins_size; // dynamic

	static float P_REF;
	static float P_ALT;
};

#endif
