
#include "RgtSites.h"
#include "Aligner.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "Globals.h"
#include "bamUtility.h"
#include "seqUtility.h"
#include "VcfSiteRecord.h"
#include "MeiSeq.h"

using std::cerr;
using std::endl;

float RgtSites::P_REF = 0.6;
float RgtSites::P_ALT = 0.4;

RgtSites::RgtSites( vector< vector<string> > & mei_names, vector< vector<string> > & mei_seqs )
{
	pMEseqs = &mei_seqs;
	pMEnames = &mei_names;
	nsample = 0;
}

void RgtSites::LoadSitesFromVcf( string & vcf_name )
{

	map<int, map<string, int> > mei_name_map;
	constructMeiNameMap( mei_name_map );

	std::ifstream vcf;
	vcf.open( vcf_name.c_str() );
	if ( !vcf.is_open() ) {
		cerr << "ERROR: unable to open: " << vcf_name << endl;
		exit(20);
	}
	string line;
	while( getline( vcf, line ) ) {
		VcfSiteRecord vrec;
		vrec.SetFromLine( line );
		string chr = vrec.GetChromosome();
		RgSiteInfo rsi;
		string mtype_str = vrec.GetMeiType();
		rsi.meitype = GetMeiIndexFromName( mtype_str );
		string info_name = string("SUB");
		string subtype_str = vrec.GetInfoString( info_name );
		rsi.subtype = mei_name_map[rsi.meitype][subtype_str];
		vector<int> sted;
		info_name = "LA";
		vrec.GetInfoPairValue( sted, info_name );
		rsi.mstart = sted[0];
		rsi.left_end = sted[1];
		sted.clear();
		info_name = "RA";
		vrec.GetInfoPairValue( sted, info_name );
		rsi.right_start = sted[0];
		rsi.mend = sted[1];
	}
	vcf.close();
	// set stat chr
	if (siteList.empty()) {
		cerr << "ERROR: [AsbSites::LoadCandidateSitesFromSites] No available site! Check if vcfs are emtpy" << std::endl;
		exit(30);
	}
	map<string, map<int, RgSiteInfo> >::iterator t = siteList.begin();
	stat_chr = t->first;
}

void RgtSites::constructMeiNameMap( map<int, map<string, int> > & mmap )
{
	for(int m=0; m<NMEI; m++) {
		mmap[m];
		for( int t=0; t<(*pMEnames)[m].size(); t++ ) {
			mmap[m][(*pMEnames)[m][t]] = t;
		}
	}
}


void RgtSites::ReGenotypeFromSingleBam( string & sample_name, string & bam_name )
{
	// add sample index
	int sp;
	if ( sampleIndex.find(sample_name) != sampleIndex.end() ) {
		sp = sampleIndex[sample_name];
		cerr << "Warning: " << sample_name << " has been re-genotyped! Re-do this and replace previous info!" << endl;
	}
	else {
		nsample++;
		sampleIndex[sample_name] = nsample - 1;
		avr_ins_size = GetAvrInsSize( bam_name, stat_chr );
	}

	// add by site
	for( map<string, map<int, RgSiteInfo> >::iterator c = siteList.begin(); c!= siteList.end(); c++ ) {
		for( map<int, RgSiteInfo>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
			SamFile sam;
			SamFileHeader sam_header;
			int breakp = p->first;
			string entire_mei_seq = (*pMEseqs)[p->second.meitype][p->second.subtype];
			string mei_seq_left = entire_mei_seq.substr(p->second.mstart, p->second.left_end - p->second.mstart+1);
			string mei_seq_right = entire_mei_seq.substr(p->second.right_start, p->second.mend - p->second.right_start + 1);
			int st = p->first - WIN;
			int ed = p->first + WIN;
			OpenBamAndBai( sam, sam_header, bam_name );
			bool section_status = sam.SetReadSection( c->first.c_str(), st, ed );
			if (!section_status)
				continue;
			SamRecord rec;
			while( sam.ReadRecord(sam_header, rec) ) {
				if ( rec.getReadLength() < MIN_CLIP )
					continue;
				if (rec.getFlag() & 0x2) {
					addProperRecord( rec, breakp, p->second, mei_seq_left, mei_seq_right);
				}
				else {
					addDiscRecord(sam, sam_header, rec, breakp, p->second, mei_seq_left, mei_seq_right);
				}
			}
			sam.Close();
		}
	}
}


void RgtSites::addProperRecord( SamRecord & rec, int breakp, RgSiteInfo & rsi, string & mei_seq_left, string & mei_seq_right)
{
	GtInfo* ginfo = &rsi.gt_info[nsample-1];

	int st = rec.get1BasedPosition();
	int ed = rec.get1BasedAlignmentEnd();
	int rlen = rec.getReadLength();
	if ( ed < breakp - rlen/2 || st > breakp + rlen/2 )
		return;

	int clip_len = GetMaxClipLen( rec );
	// if no clip, add based on record mapQ
	if (abs(clip_len)==0) {
		if ( st < breakp && ed > breakp ) {
			ginfo->ad_ref++;
			float alt_logp = -(float)rec.getMapQuality() / 10 / log(10);
			float ref_logp = log(1 - exp(alt_logp) );
			float het_logp = log( exp(alt_logp) * P_ALT + exp(het_logp) * P_REF );
			ginfo->log_gl_ref += ref_logp;
			ginfo->log_gl_alt += alt_logp;
			ginfo->log_gl_het += het_logp;
		}
		return;
	}

	// check if too short clip
	if ( abs(clip_len) < MIN_CLIP ) {
		ginfo->ad_undef++;
		return;
	}

	// map clip to subtype
	string clip_seq = GetMaxClipSeq( rec );
	if ( !rsi.is_plus_strand )
		clip_seq = RevCompSeq(clip_seq);
	string *pseq;
	if ( clip_len > 0 )
		pseq = &mei_seq_left;
	else
		pseq = &mei_seq_right;

	Aligner al( clip_seq, *pseq );
	float alt_logp = al.GetLogAlignProb();
	if (alt_logp == 1) { // unable to map
		ginfo->ad_undef++;
		return;
	}
	// set alt based
	setGinfoFromAltLogP( *ginfo, alt_logp);
}


void RgtSites::addDiscRecord( SamFile & sam, SamFileHeader & sam_header, SamRecord & rec, int breakp, RgSiteInfo & rsi, string & mei_seq_left, string & mei_seq_right)
{
	GtInfo* ginfo = &rsi.gt_info[nsample-1];

	int st = rec.get1BasedPosition();
	int ed = rec.get1BasedAlignmentEnd();

	// skip unmap
	if ( rec.getFlag() & 0x4 )
		return;
	if ( rec.getFlag() & 0x8 )
		return;
	if ( rec.getMapQuality() == 0 )
		return;

	if ( st < breakp && ed > breakp ) {
		ginfo->ad_undef++;
		return;
	}
	string mate_name = rec.getMateReferenceName();
	if ( mate_name.compare(rec.getReferenceName()) == 0 && abs(rec.getInsertSize()) < 3*avr_ins_size ) {
		ginfo->ad_undef++;
		return;
	}

	SamRecord mate_rec;
	SetDiscMateRec( mate_rec, rec, sam, sam_header );
	string mate_seq = mate_rec.getSequence();
	if (!rsi.is_plus_strand)
		mate_seq = RevCompSeq(mate_seq);
	string * pseq;
	if ( !rec.getFlag() & 0x10 ) // left anchor
		pseq = &mei_seq_left;
	else
		pseq = &mei_seq_right;

	float anchor_qual = rec.getMapQuality();
	Aligner al( mate_seq, *pseq );
	float alt_logp = al.GetLogAlignProb() - anchor_qual/10*log(10);
	if (alt_logp == 1) { // unable to map
		ginfo->ad_undef++;
		return;
	}
	setGinfoFromAltLogP( *ginfo, alt_logp );
}

void RgtSites::setGinfoFromAltLogP( GtInfo & ginfo, float alt_logp )
{
	float ref_logp = log( 1 - exp(alt_logp) );
	float het_logp = log( exp(alt_logp)*P_ALT + exp(het_logp)*P_REF );
	ginfo.log_gl_ref += ref_logp;
	ginfo.log_gl_alt += alt_logp;
	ginfo.log_gl_het += het_logp;
	ginfo.ad_alt++;
}

/* AF/AC were also calculated here
*/
void RgtSites::PrintVcf( string & vcf_name )
{

}



