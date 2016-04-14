#include "AsbSites.h"
#include "MeiSeq.h"
#include "Globals.h"
#include "bamUtility.h"
#include "morphError.h"
#include "seqUtility.h"
#include <math.h>
#include <algorithm>
#include <iomanip>      // std::setprecision

using std::endl;

AsbSites::AsbSites( vector< vector<string> > & mei_names, vector< vector<string> > & mei_seqs )
{
	MIN_SPAN = 10;
	MIN_QUALITY = 20;
	pCandidateSites = new AsbSiteList [NMEI];
	sample_count = 0;
	sample_depth.clear();
	avr_ins_size_across_all = 0;
	pMEseqs = &mei_seqs;
	pMEnames = &mei_names;
	first_add = 1;
}


void AsbSites::Assembly( map<int, map<string, map<int, int> > > & candidate_sites, vector<SingleSampleInfo> & msinfo )
{
	LoadCandidateSitesFromSites(candidate_sites);

	// add from single bam
	for( int s=0; s<msinfo.size(); s++ )
		AddEvidenceFromSingleBam( msinfo[s] );
	first_add = 0;

	// first assembly
	assemblySites();
	int sum = _getCompSiteCount();
	string str = "  Identified " + std::to_string(sum) + " sites in first round";
	morphMessage( str );	

	// rescue
	rescueCsi( msinfo );
	sum = _getCompSiteCount();
	str = "  Identified " + std::to_string(sum) + " sites after rescue";
	morphMessage( str );

	// remove nearby
	removeNearbySites();

	// refine breakpoint
	for( int s=0; s<msinfo.size(); s++ )
		ImplementBreakpFromSingleBam( msinfo[s] );
}


// mtype -> chr-> breakp -> extend
void AsbSites::LoadCandidateSitesFromSites( map<int, map<string, map<int, int> > > & candidate_sites )
{
	for( int m=0; m<NMEI; m++ ) { // mtype
		for( map<string, map<int, int> >::iterator c=candidate_sites[m].begin(); c!= candidate_sites[m].end(); c++) { // chr
			pCandidateSites[m][c->first];
			for( map<int, int>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
				pCandidateSites[m][c->first][p->first].nclip = 0;
				pCandidateSites[m][c->first][p->first].ndisc = 0;				
				pCandidateSites[m][c->first][p->first].raw_breakp = p->first;
				pCandidateSites[m][c->first][p->first].true_breakp = -1;
				pCandidateSites[m][c->first][p->first].disc_breakp = -1;
				pCandidateSites[m][c->first][p->first].extend = p->second;
//				pCandidateSites[m][c->first][p->first].doing_rescue = 0;
				pCandidateSites[m][c->first][p->first].polyA_left = 0;
				pCandidateSites[m][c->first][p->first].polyA_right = 0;
				pCandidateSites[m][c->first][p->first].polyT_left = 0;
				pCandidateSites[m][c->first][p->first].polyT_right = 0;
			}
		}
	}
	/* set stat chr
	for(int m=0; m<NMEI; m++) {
		if (candidate_sites[m].size() == 0)
			continue;
		map<string, map<int, int> >::iterator c=candidate_sites[m].begin();
		stat_chr = c->first;
		break;
	}
	if (stat_chr.empty())
		morphError("[AsbSites::LoadCandidateSitesFromSites] No available site! Check if vcfs are emtpy");
	*/
}


void AsbSites::AddEvidenceFromSingleBam( SingleSampleInfo & si)
{
	// set in class
	current_bam_name = si.bam_name;
	AsbSiteList* p = pCandidateSites;
	avr_ins_size = si.avr_ins_size;
	avr_read_length = si.avr_read_length;
	min_read_length = std::min( avr_read_length * 0.33, MIN_CLIP*1.5);
//	string str = "   Adding " + si.sample_name;
//	morphMessage(str);	
	// update inner stats if needed
	if (first_add) {
		sample_count++;
		sample_depth.push_back( si.depth );
		float a1 = (float)si.depth / (float)(si.depth + total_depth);
		avr_ins_size_across_all = avr_ins_size*(1-a1) + si.avr_ins_size*a1;
		avr_read_length_across_all = avr_read_length*(1-a1) + si.avr_read_length*a1;
		total_depth += si.depth;
	}

	// adding from bam
	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam, bam_header, current_bam_name );
	// alt ptr to set disc mate
	SamFile alt_bam;
	SamFileHeader alt_bam_header;
	OpenBamAndBai( alt_bam, alt_bam_header, current_bam_name );
	for( int m=0; m<NMEI; m++) {
		for( AsbSiteList::iterator c=p->begin(); c!= p->end(); c++ ) {
			for( map<int, AsbSiteListElement>::iterator it=c->second.begin(); it!=c->second.end(); it++ ) {
				int st = it->first - it->second.extend;
				int ed = it->first + it->second.extend;
				if (st<0 || ed <0)
					morphError("Site cluster < 0");
//std::cout << "m=" << m << ", " << c->first << ":" << st << "-" << ed << "; " << it->first << std::endl;
				bool section_status = bam.SetReadSection( c->first.c_str(), st, ed ); // set as overlap here!
				if (!section_status)
					continue;
				int max_lc = GetMaxIntervalReadCount( sample_depth[sample_depth.size()-1], avr_read_length, ed - st );
				SamRecord sam_rec;
				int lc = 0;
				while( bam.ReadRecord(bam_header, sam_rec) ) {
					lc++;
					if (lc >= max_lc)
						break;
				}
				if (lc >= max_lc)
					continue;
				section_status = bam.SetReadSection( c->first.c_str(), st, ed ); // set as overlap here!
				if (!section_status)
					continue;
//str = "adding " + std::to_string(st) + "-" + std::to_string(ed) + "; lc=" + std::to_string(lc) + ";mlc=" + std::to_string(max_lc);				
//morphMessage(str);
				addSectionEvidence( m, it, bam, bam_header, alt_bam, alt_bam_header );
			}
		}
		p++;
	}
	bam.Close();
	alt_bam.Close();
}


// add single section
void AsbSites::addSectionEvidence( int m, map<int, AsbSiteListElement>::iterator & it, SamFile & bam, SamFileHeader & bam_header, SamFile & alt_bam, SamFileHeader & alt_bam_header)
{
	SamRecord sam_rec;	
	while( bam.ReadRecord(bam_header, sam_rec) ) {
		if (!qcCheck(sam_rec.getFlag()))
			continue;
		if (IsSupplementary(sam_rec.getFlag()))
			continue;
		if (sam_rec.getMapQuality() < MIN_QUALITY)
			continue;
		if ( sam_rec.getReadLength() <= min_read_length )
			continue;
//std::cout << sam_rec.getReadName() << std::endl;		
		if ( sam_rec.getFlag() & 0x2 )
			addProperRead( m, it, sam_rec );
		else
			addDiscRead( m, it, sam_rec, alt_bam, alt_bam_header );
	}
}

void AsbSites::addProperRead( int m, map<int, AsbSiteListElement>::iterator & it, SamRecord & sam_rec )
{
	int breakpoint = it->second.raw_breakp;

	Cigar * myCigar = sam_rec.getCigarInfo();
	int begin_clip = myCigar->getNumBeginClips();
	int end_clip = myCigar->getNumEndClips();

/*	int max_clip = GetMaxClipLen( sam_rec );
	if (abs(max_clip)<MIN_CLIP) {
		return;
	}
*/
//	addToSpanning(sam_rec.get1BasedPosition(), sam_rec.get1BasedAlignmentEnd(), it);
//	addToDepth(sam_rec.get1BasedPosition(), sam_rec.get1BasedAlignmentEnd(), it);

    // check overlap with break point & re-map
	string seq = sam_rec.getSequence();
	 // right cluster (left clip)
	if (begin_clip>MIN_CLIP) {
		if( sam_rec.get1BasedPosition() < breakpoint - sam_rec.getReadLength()*2 || sam_rec.get1BasedPosition() > breakpoint + sam_rec.getReadLength()*2 )
			return;
		string clip_str = seq.substr(0, begin_clip);
		bool mapped = addToSubtypeInfo( m, clip_str, it, 0 );
		bool pstatus = ContainRepetative( seq, 'A');
//std::cout << "proper right: " << sam_rec.getReadName() << " " << begin_clip << " status=" << pstatus << std::endl;		
		if (mapped)
			updateTrueBreakp( it->second, sam_rec.get1BasedPosition() );
		else {
			if (pstatus)
				updateTrueBreakp( it->second, sam_rec.get1BasedPosition() );
		}
		if (pstatus)
			it->second.polyA_right++;
	}
	
	// left cluster
	if (end_clip>MIN_CLIP) {
		if ( sam_rec.get1BasedAlignmentEnd() < breakpoint - sam_rec.getReadLength() || sam_rec.get1BasedAlignmentEnd() > breakpoint + sam_rec.getReadLength() )
			return;
		string clip_str = seq.substr(seq.length() -end_clip, end_clip);
		bool mapped = addToSubtypeInfo( m, clip_str, it, 1);
		bool pstatus = ContainRepetative( seq, 'T');
//std::cout << "proper left: " << sam_rec.getReadName() << " " << end_clip << " status=" << pstatus << std::endl;
		if (mapped)
			updateTrueBreakp( it->second, sam_rec.get1BasedAlignmentEnd() );
		else {
			if (pstatus)
				updateTrueBreakp( it->second, sam_rec.get1BasedAlignmentEnd() );
		}
		if (pstatus)
			it->second.polyT_left++;
	}
}

void AsbSites::addDiscRead( int m, map<int, AsbSiteListElement>::iterator & it, SamRecord & sam_rec, SamFile & alt_bam, SamFileHeader & alt_bam_header )
{
	if (sam_rec.getFlag() & 0x4)
		return;

	int breakpoint = it->second.raw_breakp;
	if( sam_rec.get1BasedPosition() > breakpoint + avr_ins_size*2 || sam_rec.get1BasedAlignmentEnd() < breakpoint - avr_ins_size*2 )
		return;	

	// see if anchor end needs to be added to proper
	addProperRead( m, it, sam_rec );

	if ( !(sam_rec.getFlag() & 0x8) ) { // disc qc check
		string ref_name = sam_rec.getReferenceName();
		string mate_ref = sam_rec.getMateReferenceName();
		if (ref_name.compare(mate_ref)==0 && abs(sam_rec.getInsertSize()) < 3*avr_ins_size)
			return;
	}

	SamRecord mate_rec;
	bool mate_set = SetDiscMateRec( mate_rec, sam_rec, alt_bam, alt_bam_header );
	if (!mate_set)
		return;
	string seq = GetCorrectMateSeq(sam_rec, mate_rec);
	bool mapped;
	if (sam_rec.getFlag() & 0x10) { // mate reverse: right anchor
		mapped = addToSubtypeInfo( m, seq, it, 0);
		bool pstatus = ContainRepetative( seq, 'A');
	    if (pstatus)
			it->second.polyA_right++;
		pstatus = ContainRepetative( seq, 'T');
		if (pstatus)
			it->second.polyT_right++;
	}
	else { // left anchor
		mapped = addToSubtypeInfo( m, seq, it, 1);
		bool pstatus = ContainRepetative( seq, 'A');
	    if (pstatus)
			it->second.polyA_left++;
		pstatus = ContainRepetative( seq, 'T');
		if (pstatus)
			it->second.polyT_left++;
	}

	// set breakp back up
	if (mapped)
		updateDiscBreakp( it->second, sam_rec );

}


bool AsbSites::addToSubtypeInfo( int m, string & seq, map<int, AsbSiteListElement>::iterator & it, bool is_left_cluster )
{
	bool mapped = 0;
	if ( it->second.p_cluster.empty() ) {
		it->second.p_cluster.resize(4);
		for( int i=0; i<4; i++ ) {
			it->second.p_cluster[i].resize((*pMEseqs)[m].size());
			for( int j=0; j<(*pMEseqs)[m].size(); j++ )
				it->second.p_cluster[i][j].coverage.resize( (*pMEseqs)[m][j].length(), 0 );
		}
	}

	vector< vector<SubtypeInfo> >::iterator p = it->second.p_cluster.begin();
	if (!is_left_cluster)
		p++;
	// plus strand
	bool new_map = _realignAndAddSingleSide( m, *p, it, seq);
	if (!mapped)
		if (new_map)
			mapped = 1;
//if (new_map)
//std::cout << "plus: " << new_map << ", left : " << is_left_cluster << " " << seq << std::endl;

	// minus strand
	p++; p++;
	string rv_seq = RevCompSeq(seq);
	new_map = _realignAndAddSingleSide( m, *p, it, rv_seq);
	if (!mapped)
		if (new_map)
			mapped = 1;
//if (new_map)
//std::cout << "minus: " << new_map << ", left : " << is_left_cluster << " " << seq << std::endl;	
	return mapped;
}

bool AsbSites::_realignAndAddSingleSide( int m, vector<SubtypeInfo> & p, map<int, AsbSiteListElement>::iterator & it, string & seq )
{
	bool mapped = 0;
	for(int s=0; s<(*pMEseqs)[m].size(); s++) {
//std::cout << "size=" << (*pMEseqs)[m].size() << ", m=" << m << ", s=" << s << endl;
		Aligner al( (*pMEseqs)[m][s], seq );
		if ( !al.IsMapped() )
			continue;
		if (!mapped)
			mapped = 1;
//std::cout << "add s= " << s << std::endl;
		_updateBaseCoverageWithRead(p[s].coverage, al );
		p[s].evidence++;
	}
	return mapped;
}


void AsbSites::_updateBaseCoverageWithRead( vector<float> & coverage, Aligner & al )
{
	int st = al.GetLeftRefPosition();
	int ed = al.GetRightRefPosition();
//std::cout << "st=" << st << ", ed=" << ed << std::endl;
	float base_score = (float)al.GetScore() / (ed - st + 1) / MATCH;
	for(int i=st; i<=ed; i++)
		coverage[i] += base_score;
}

/******************** assembly ********************
	assembly using recorded info
	convert to AsbSite to  CompleteSites (merge different families)
	try to rescue failed sites
	merge nearby
	*/

void AsbSites::assemblySites()
{
	// convert & delete redundant keys
	for(int m=0; m<NMEI; m++) {
		for( AsbSiteList::iterator pchr = pCandidateSites[m].begin(); pchr != pCandidateSites[m].end(); pchr++ ) {
			compSites[pchr->first];
			for( map<int, AsbSiteListElement>::iterator p = pchr->second.begin(); p != pchr->second.end(); p++ ) {
				CompleteSiteInfo new_csi;
				setCsiBreakPoint( new_csi, p->second );
				new_csi.mei_type = m;
				for(int j=0; j<NMEI; j++)
					new_csi.merged[j] = 0;
				new_csi.merged[m] = 1;
//std::cout << "doing " << pchr->first << ":" << new_csi.breakpoint << ", m = " << m << ", size= " << compSites[pchr->first].size() << std::endl;
				bool assembled = setCompleteSiteInfoFromAsbSiteListElement( m, new_csi, p->second );
//				if ( !assembled ) // need to do recue. Cannot skip unassembled!
//					continue;
				int new_key = round(new_csi.breakpoint / WIN);
				map<int, CompleteSiteInfo>::iterator t = compSites[pchr->first].find(new_key);
				if ( t != compSites[pchr->first].end() ) {
					bool new_csi_better = isNewCsiBetter( new_csi, t );
					if (new_csi_better) { // assign and copy merge info
						for(int i=0; i<NMEI; i++)
							new_csi.merged[i] = (new_csi.merged[i] || t->second.merged[i]);
						t->second = new_csi;
					}
					else 
						t->second.merged[m] = 1;
				}
				else {
					compSites[pchr->first][new_key] = new_csi;
//					std::cout << compSites[pchr->first].size() << endl;
				}
			}
		}
	}
}


void AsbSites::rescueCsi( vector<SingleSampleInfo> & msinfo )
{
	// transfer rescue list to pCandidate->AsbSiteList
	constructRescueAsbList();

	// add single sample info
	for( int s=0; s<msinfo.size(); s++ )
		AddEvidenceFromSingleBam( msinfo[s] );

	// assembly again
	assemblySites();

	// remove unassembled
	removeUnAssembledSites();
}


void AsbSites::removeUnAssembledSites()
{
	for( map<string, map<int, CompleteSiteInfo> >::iterator c=compSites.begin(); c!=compSites.end(); c++ ) {
		vector<int> keys;
		for( map<int, CompleteSiteInfo>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
			if (!p->second.assembled)
				keys.push_back(p->first);
		}
		for(int i=0; i<keys.size(); i++)
			c->second.erase(keys[i]);
	}

}



void AsbSites::removeNearbySites()
{
	// check & remove nearby key in csi list
	for( CompleteSiteList::iterator c = compSites.begin(); c!= compSites.end(); c++ ) {
		// get key list
		vector<int> exist_keys;
		exist_keys.resize( c->second.size() );
		int k = 0;
		for( map<int, CompleteSiteInfo>::iterator p = c->second.begin(); p != c->second.end(); p++ ) {
			exist_keys[k] = p->first;
			k++;
		}
		// then merge nearby +/- 1
		for( k=0; k<exist_keys.size(); k++ ) {
			// +1
			map<int, CompleteSiteInfo>::iterator t2 = c->second.find(k+1);
			if ( t2 != c->second.end() ) {
				map<int, CompleteSiteInfo>::iterator t = c->second.find(k);
				if ( t != c->second.end() ) {
					bool new_csi_better = isNewCsiBetter( t2->second, t );
					if (new_csi_better)
						c->second.erase( t->first );
					else
						c->second.erase( t2->first );
				}
			}
			// -1
			t2 = c->second.find(k-1);
			if ( t2 != c->second.end() ) {
				map<int, CompleteSiteInfo>::iterator t = c->second.find(k);
				if ( t != c->second.end() ) {
					bool new_csi_better = isNewCsiBetter( t2->second, t );
					if (new_csi_better)
						c->second.erase( t->first );
					else
						c->second.erase( t2->first );
				}
			}
		}
	}
}

// clear asb list first
// add from compSite
void AsbSites::constructRescueAsbList()
{
	// clear
	for(int i=0; i<NMEI; i++)
		pCandidateSites[i].clear();

	// add from assembled list
	for( CompleteSiteList::iterator c=compSites.begin(); c!=compSites.end(); c++ ) {
		for( map<int, CompleteSiteInfo>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
			bool to_rescue = 0;
			if (p->second.assembled) {
				if (p->second.stringent_filter)
					to_rescue = 1;
			}
/*			if ( !p->second.assembled)
				to_rescue = 1;
			else {
				if (p->second.stringent_filter)
					to_rescue = 1;
			}
*/			
			// then prepare for rescue
			if ( to_rescue ) {
				for(int m=0; m<NMEI; m++) {
					if ( !(p->second.merged[m]) ) { // not merged
						pCandidateSites[m][c->first][p->second.breakpoint].raw_breakp = p->second.breakpoint;
						pCandidateSites[m][c->first][p->second.breakpoint].true_breakp = -1;
						pCandidateSites[m][c->first][p->second.breakpoint].nclip = 0;
						pCandidateSites[m][c->first][p->second.breakpoint].extend = avr_ins_size_across_all*1.5;
//						pCandidateSites[m][c->first][p->second.breakpoint].doing_rescue = 1;
						pCandidateSites[m][c->first][p->second.breakpoint].polyA_left = 0;
						pCandidateSites[m][c->first][p->second.breakpoint].polyA_right = 0;
						pCandidateSites[m][c->first][p->second.breakpoint].polyT_left = 0;
						pCandidateSites[m][c->first][p->second.breakpoint].polyT_right = 0;
					}
				}
			}
		}
	}
}


void AsbSites::ImplementBreakpFromSingleBam( SingleSampleInfo & si )
{
	min_read_length = std::min( avr_read_length * 0.33, MIN_CLIP*1.5);
	// adding
	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam, bam_header, current_bam_name );
	
	for( CompleteSiteList::iterator c = compSites.begin(); c!= compSites.end(); c++ ) {
		string chr = c->first;
		for(map<int, CompleteSiteInfo>::iterator t=c->second.begin(); t!= c->second.end(); t++) {
			int st = t->second.breakpoint - WIN/2;
			int ed = t->second.breakpoint + WIN/2;
			bool section_status = bam.SetReadSection( c->first.c_str(), st, ed ); // set as overlap here!
			if (!section_status) {
				string str = "At section " + chr + ":" + std::to_string(st) + "-" + std::to_string(ed);
				str += ", no read in assembled section!";
				morphWarning(str);
				continue;
			}
			implementSection( bam, bam_header,t->second );

		}
	}
	bam.Close();
}




void AsbSites::implementSection( SamFile & bam, SamFileHeader & bam_header, CompleteSiteInfo & csi )
{
	SamRecord sam_rec;
	while( bam.ReadRecord(bam_header, sam_rec) ) {
		if (!qcCheck(sam_rec.getFlag()))
			continue;
		if (IsSupplementary(sam_rec.getFlag()))
			continue;
		if (sam_rec.getMapQuality() < MIN_QUALITY)
			continue;
		if ( sam_rec.getReadLength() <= min_read_length )
			continue;
		addToSpanning( sam_rec, csi );
		addToDepth(sam_rec, csi);
	}
}


void AsbSites::addToSpanning( SamRecord & rec, CompleteSiteInfo & csi )
{
	int st = rec.get1BasedPosition();
	int ed = rec.get1BasedAlignmentEnd();
	if (ed - st < MIN_SPAN)
		return;

	int breakpoint = csi.breakpoint;
	if (st<breakpoint && ed>breakpoint) {
		if (st + MIN_SPAN > breakpoint)
			return;
		if ( ed - MIN_SPAN < breakpoint)
			return;
		csi.n_spanning++;
	}
}

void AsbSites::addToDepth( SamRecord & rec, CompleteSiteInfo & csi )
{
	int st = rec.get1BasedPosition();
	int ed = rec.get1BasedAlignmentEnd();
	int breakpoint = csi.breakpoint;
	int new_st, new_ed;
	new_st = (st<breakpoint-WIN/2) ? breakpoint-WIN/2 : st;
	new_ed = (ed>breakpoint+WIN/2) ? breakpoint+WIN/2 : ed;
	if (new_ed > new_st)
		csi.avr_flank_depth += (new_ed - new_st);
}




/*
	if unable to assemble, return FALSE
*/
bool AsbSites::setCompleteSiteInfoFromAsbSiteListElement( int m, CompleteSiteInfo & csi, AsbSiteListElement & ase )
{
	csi.merged[m] = 1;
	if ( ase.p_cluster.size() == 0 ) {
		csi.assembled = 0;
		return 0;
	}

	// get most likely strand & subtype first
	int n = ase.p_cluster[0].size(); // # subtypes
	bool is_plus_strand;
	int optimized_subtype;

	// refine first
	vector< vector<SubtypeInfo> >::iterator p = ase.p_cluster.begin();
	for( int s=0; s<4; s++ ) {
		for(int t=0; t<p->size(); t++)
			refineAnchorCluster( (*p)[t] );
		p++;
	}

	// then add up score and compare
	vector< vector<float> > score_sum; // 0 is plus, 1 is minus strand -> then every subtype
	score_sum.resize(2);
	score_sum[0].resize( n );
	score_sum[1].resize( n );
	for(int i=0; i<ase.p_cluster[0].size(); i++) {
		score_sum[0][i] = _sumFloatVector( ase.p_cluster[0][i].coverage );
		score_sum[0][i] += _sumFloatVector( ase.p_cluster[1][i].coverage );
	}
	for(int i=0; i<ase.p_cluster[1].size(); i++) {
		score_sum[1][i] = _sumFloatVector( ase.p_cluster[2][i].coverage );
		score_sum[1][i] += _sumFloatVector( ase.p_cluster[3][i].coverage );
	}

	vector<float>::iterator plus_pmax = std::max_element(score_sum[0].begin(), score_sum[0].end());
	vector<float>::iterator minus_pmax = std::max_element(score_sum[1].begin(), score_sum[1].end());

	if ( *plus_pmax + (*minus_pmax) == 0 ) { // unable to assemble
		csi.assembled = 0;
		return 0;
	}

	if (*plus_pmax > *minus_pmax) {
		is_plus_strand = 1;
		optimized_subtype = plus_pmax - score_sum[0].begin();
	}
	else if (*plus_pmax < *minus_pmax) {
		is_plus_strand = 0;
		optimized_subtype = minus_pmax - score_sum[1].begin();
	}
	else {
	// compare coverage
		int plus_max_index = plus_pmax - score_sum[0].begin();
		int minus_max_index = minus_pmax - score_sum[1].begin();
		float plus_cov_sum = 0;
		for(int i=0; i<(*pMEseqs)[m][plus_max_index].size(); i++) {
			plus_cov_sum += ase.p_cluster[0][plus_max_index].coverage[i];
			plus_cov_sum += ase.p_cluster[2][plus_max_index].coverage[i];
		}
		int minus_cov_sum = 0;
		for(int i=0; i<(*pMEseqs)[m][minus_max_index].size(); i++) {
			minus_cov_sum += ase.p_cluster[1][minus_max_index].coverage[i];
			minus_cov_sum += ase.p_cluster[3][minus_max_index].coverage[i];
		}
		if (plus_cov_sum>minus_cov_sum) {
			is_plus_strand = 1;
			optimized_subtype = plus_max_index;
		}
		else if (plus_cov_sum<minus_cov_sum) {
			is_plus_strand = 0;
			optimized_subtype = minus_max_index;
		}
		else {
			is_plus_strand = 1;
			optimized_subtype = plus_max_index;
			morphWarning("[AsbSites::setCompleteSiteInfoFromAsbSiteListElement] equal rank for +/- strand!");
		}
	}

	// then set to csi from ase
	bool convert_success = _convertMostLikelySubtypeToCsiInfo( csi, m, is_plus_strand, optimized_subtype, ase );	
	if (convert_success) {
		csi.assembled = 1;
		return 1;
	}
	else {
		csi.assembled = 0;
		return 0;
	}
}


/* 
previous;
	refine left and right anchor
current: get the width window with maximized value

finally:
	set middle part as coverage = 0
*/
void AsbSites::refineAnchorCluster( SubtypeInfo & sbi )
{
	if ( sbi.coverage.size() < avr_ins_size_across_all * 3 )// no refine
		return;
	int width = avr_ins_size_across_all * 1.5;
	int last_start = -1; // mark start of last window

	// first check partially
	vector<float> cum1; // recording vec
	cum1.resize( (int)sbi.coverage.size() / width + 1, 0);
	int i1 = 0; // starting base
	int j1 = 0; // starting base's index in cum1
	while( i1 < sbi.coverage.size() ) {
		for( int i=i1; i<std::min(i1+width, (int)sbi.coverage.size()-1); i++ )
			cum1[j1] += sbi.coverage[i];
		i1 += width;
		j1++;
	}
	// correct the last window by adding the positions before that
	if (cum1.size() * width > sbi.coverage.size()) {
		j1 = cum1.size() - 1;
		for(int i=sbi.coverage.size() - width; i<width*(cum1.size()-1); i++)
			cum1[j1] += sbi.coverage[i];
		last_start = sbi.coverage.size() - width;
	}

	// then check top 2 nearby
	vector<float>::iterator cum1_max1 = std::max_element( cum1.begin(), cum1.end() );
	if (*cum1_max1 == 0)
		return;
	int max_val1 = *cum1_max1;
	*cum1_max1 = -1;
	vector<float>::iterator cum1_max2 = std::max_element( cum1.begin(), cum1.end() );

	// 1st: try to find maximized cluster
	int p1 = cum1_max1 - cum1.begin();
	int val1;
	int st1;
	if (p1 == cum1.size() - 1 && last_start != -1) { // last one & corrected
		val1 = max_val1;
		st1 = last_start;
	}
	else {
		p1 *= width;
		st1 = std::max( 0,p1 - width/2 + 1);
		int ed1 = std::min( (int)sbi.coverage.size() - width - 1, p1 + width/2 - 1);
		if (st1 > ed1)
			morphError("[AsbSites::refineAnchorCluster] st1 > ed1");
		vector<float> cum2;
		cum2.resize( ed1-st1+1, 0 );
		for(int i=st1; i<=ed1; i++) {
			for(int j=i; j<i+width; j++)
				cum2[i-st1] += sbi.coverage[j];
		}
		vector<float>::iterator t = std::max_element( cum2.begin(), cum2.end() );
		int idx1 = t - cum2.begin();
		st1 += idx1;
		val1 = *t;
	}

	int left_index = 0;
	if ( (*cum1_max2) == 0 || (max_val1) > 2 * (*cum1_max2) ) { // skip 2nd
		left_index = st1;
	}
	else { // test 2nd
		// also check if 2nd is last
		int p2 = cum1_max2 - cum1.begin();
		int max_val2;
		int st2;
		if (p2 == cum1.size() -1 && last_start != -1) {
			max_val2 = *cum1_max2;
			st2 = last_start;
		}
		else {
			p2 *= width;
			st2 = std::max( 0,p2 - width/2 + 1);
			int ed2 = std::min( (int)sbi.coverage.size() - width - 1, p2 + width/2 - 1);
			if (st2 > ed2)
				morphError("[AsbSites::refineAnchorCluster] st2 > ed2");
			vector<float> cum2;
			cum2.resize( ed2 - st2 + 1, 0 );
			for(int i=st2; i<=ed2; i++) {
				for(int j=i; j<i+width; j++)
					cum2[i-st2] += sbi.coverage[j];
			}
			vector<float>::iterator t = std::max_element( cum2.begin(), cum2.end() );
			max_val2 = *t;
			int idx2 = t - cum2.begin();
			st2 += idx2;		
		}
		if (max_val2 > val1) // prefer 5' truncate
			left_index = st2;
		else
			left_index = st1;
	}
	
	if ( left_index > 0 )
		std::fill( sbi.coverage.begin(), sbi.coverage.begin() + left_index, 0 );
	std::fill( sbi.coverage.begin() + left_index + width, sbi.coverage.end(), 0 );
}


bool AsbSites::isNewCsiBetter( CompleteSiteInfo & new_csi, map<int, CompleteSiteInfo>::iterator & t )
{
	if ( new_csi.filter + t->second.filter != 0 && new_csi.filter*t->second.filter == 0 ) {
		if (new_csi.filter == 0)
			return 1;
		else
			return 0;
	}
	if (new_csi.score_sum > t->second.score_sum)
		return 1;
	else if (new_csi.score_sum < t->second.score_sum)
		return 0;
	if (new_csi.mei_type<t->second.mei_type)
		return 1;
	else
		return 0;
}


void AsbSites::PrintToVcf(string & vcf_name)
{
	std::ofstream vcf;
	vcf.open(vcf_name.c_str());
	if (!vcf.is_open())
		morphErrorFile(vcf_name);

	// print header
	printHardFilterSiteVcfHeader( vcf );

	// set the statisc
	_setRecordStatics();

	// print record
	for( CompleteSiteList::iterator c = compSites.begin(); c!= compSites.end(); c++ ) {
		string chr = c->first;
		for(map<int, CompleteSiteInfo>::iterator t=compSites[chr].begin(); t!= compSites[chr].end(); t++) {
			printSingleCsiToVcf(chr, vcf, t->second);
		}
	}
	vcf.close();
}

void AsbSites::printSingleCsiToVcf( string & chr, std::ofstream & vcf, CompleteSiteInfo & csi )
{
	int avr_score = round( csi.score_sum/(csi.n_left_anchor + csi.n_right_anchor) );
	string mname = GetMeiNameFromIndex(csi.mei_type);

	vcf << chr << "\t" << csi.breakpoint << "\t.\t.\t<INS:ME:" << mname;
	vcf << ">\t" << avr_score/2 << "\t" << filter_name[csi.filter] << "\tSUB=";
	vcf << (*pMEnames)[csi.mei_type][csi.subtype] << ";STR=" << strand_name[csi.strand] << ";SVLEN=" << csi.svlen;
	vcf << ";DP=" << std::setprecision(3) << csi.avr_mei_depth;
	vcf << ";SCORE=" << avr_score << ";MPOS=" << csi.left_most << "," << csi.right_most;
	vcf << ";ANCLEN=" << csi.left_anchor_length << "," << csi.right_anchor_length;
	vcf << ";MISS=" << csi.n_missing << ";SPAN=" << csi.n_spanning;
	vcf << ";FLANK=" << std::setprecision(3) << csi.avr_flank_depth/WIN;
	vcf << ";NANC=" << csi.n_left_anchor << "," << csi.n_right_anchor;
	vcf << ";POLYA=" << csi.n_left_poly << "," << csi.n_right_poly;
	if (csi.imprecise)
		vcf << ";IMPRECISE=1";
	else
		vcf << ";IMPRECISE=0";
	vcf << endl;
}


/**** lower level functions for assembly *****/

void AsbSites::_setRecordStatics()
{
 	filter_name[0] = "PASS";
 	filter_name[1] = "SVLEN";
 	filter_name[2] = "SR";
 	filter_name[3] = "MISS";
 	filter_name[4] = "S_ANCHOR_LEN";
 	filter_name[5] = "ANCHOR_LEN";
	strand_name[1] = '+';
	strand_name[0]= '-';
}


// if not assembled, return FALSE
bool AsbSites::_convertMostLikelySubtypeToCsiInfo( CompleteSiteInfo & csi, int mei_type, bool is_plus_strand, int optimized_subtype, AsbSiteListElement & ase )
{
	int subtype_length = (*pMEseqs)[mei_type][optimized_subtype].size();
	int pbase;
	if (is_plus_strand)
		pbase = 0;
	else
		pbase = 2;
	vector<float> * pleft = &ase.p_cluster[pbase][optimized_subtype].coverage;
	vector<float> * pright = &ase.p_cluster[pbase+1][optimized_subtype].coverage;
	vector<float> new_cov;
	new_cov.resize( subtype_length );
	for(int i=0; i<subtype_length; i++) {
		new_cov[i] = (*pleft)[i] + (*pright)[i];
	}

	csi.score_sum = _sumFloatVector( *pleft ) + _sumFloatVector( *pright );
	csi.strand = is_plus_strand;
	csi.mei_type = mei_type;
	csi.subtype = optimized_subtype;
	csi.right_most = _getRightMost( new_cov );
	if ( csi.right_most == 0 ) // unable to assemble
		return 0;
	csi.left_most = _getLeftMost( new_cov );
	csi.left_anchor_length = _getAnchorLength(*pleft);
	csi.right_anchor_length = _getAnchorLength(*pright);
	csi.svlen = csi.right_most - csi.left_most + 1;
	// correct polyA/T count if sv is long
	if (csi.svlen > 1000) {
		if (csi.strand)
			ase.polyA_left = 0;
		else
			ase.polyT_right = 0;
	}
	if (csi.strand) { // plus only has polyA
		csi.n_left_poly = ase.polyA_left;
		csi.n_right_poly = ase.polyA_right;
	}
	else {
		csi.n_left_poly = ase.polyT_left;
		csi.n_right_poly = ase.polyT_right;
	}
	csi.n_left_anchor = ase.p_cluster[pbase][optimized_subtype].evidence;
	csi.n_right_anchor = ase.p_cluster[pbase+1][optimized_subtype].evidence;
	csi.n_missing = _getMissingBaseCount(new_cov);
	csi.avr_mei_depth = _getSumDepth( new_cov ) / (float)(csi.svlen - csi.n_missing);
	csi.n_spanning = 0;
	csi.avr_flank_depth = 0;
	_setCsiFilter(csi);
	_setCsiStringentFilter(csi);
	return 1;
}

int AsbSites::_getLeftMost( vector<float> & coverage )
{
	int st = -1;
	for(int i=0; i<coverage.size(); i++) {
		if (coverage[i]>0) {
			st = i;
			break;
		}
	}
	if (st<0)
		return 0;
	else
		return st;
}

int AsbSites::_getRightMost( vector<float> & coverage )
{
	int ed = -1;
	for(int i=coverage.size()-1; i>=0; i--) {
		if (coverage[i]>0) {
			ed = i;
			break;
		}
	}
	if (ed<0)
		return 0;
	else
		return ed;
}

// when no coverage, return 0
int AsbSites::_getAnchorLength( vector<float> & coverage )
{
	int st = _getLeftMost( coverage );
	int ed = _getRightMost( coverage );
	if (st==ed)
		return 0;
	else
		return ed - st + 1;
}


int AsbSites::_getMissingBaseCount( vector<float> & coverage )
{
	int st = -1;
	int ed = -1;
	for(int i=0; i<coverage.size(); i++) {
		if(coverage[i]>0) {
			st=i;
			break;
		}
	}
	for(int i=coverage.size()-1; i>=0; i--) {
		if (coverage[i]>0) {
			ed = i;
			break;
		}
	}
	if ( st == -1 || ed == -1 ) {
		string str = "[AsbSites::_getMissingBaseCount] st=" + std::to_string(st) + ", ed=" + std::to_string(ed);
		morphError(str, 60);
	}
	int nmiss = 0;
	for(int i=st; i<=ed; i++) {
		if (coverage[i]==0)
			nmiss++;
	}
	return nmiss;
}

int AsbSites::_getSumDepth( vector<float> & coverage )
{
	float sum = 0;
	for(int i=0; i<coverage.size();i++)
		sum += coverage[i];
	return sum;
}


// this is to decide if need to do rescue or not
// if true, then do rescue
void AsbSites::_setCsiStringentFilter( CompleteSiteInfo & csi )
{
	if ( csi.svlen < avr_ins_size_across_all * 2) { // for shorter MEI as ALU
		if (float(csi.n_missing/csi.svlen)>0.2) {
			csi.stringent_filter = 1;
			return;
		}

	}
	if (csi.n_left_poly + csi.n_right_poly < 2) {
		csi.stringent_filter = 1;
		return;
	}
	if (csi.filter > 0) {
		csi.stringent_filter = 1;
		return;
	}
	// more stringent length
	if (csi.svlen<200) {
		csi.stringent_filter = 1;
		return;
	}
	csi.stringent_filter = 0;
}


/* csi filter
1. svlen < 90
2. SR > 3
(for svlen < avr_ins * 2)
3. %missing > 0.4
4. either anchor < svlen/3
5. anchor sum < svlen/2
(for svlen > avr_ins * 2)
*/

void AsbSites::_setCsiFilter( CompleteSiteInfo & csi )
{
	if (csi.svlen<90) {
		csi.filter = 1;
		return;
	}
	// sr
	if ( csi.n_right_anchor == 0 || csi.n_left_anchor == 0 ) {
		csi.filter = 2;
		return;
	}
	float sr = (float)csi.n_left_anchor / (float)csi.n_right_anchor;
	float sr_alt = (float)(csi.n_left_anchor + csi.n_left_poly) / (float)(csi.n_right_anchor + csi.n_right_poly); // include polyA/T
	float sr_thred;
	if (total_depth>=10)
		sr_thred = 0.3;
	else if (total_depth<2)
		sr_thred = 0.2;
	else
		sr_thred = 0.2 + (float)(total_depth - 2)*0.01;
	float sr_thred1 = float(1.0) / sr_thred;
	if (sr > sr_thred1 || sr < sr_thred) {
		if ( sr_alt > sr_thred1 || sr_alt < sr_thred ) {
			csi.filter = 2;
			return;
		}
	}
	if ( csi.svlen < avr_ins_size_across_all * 2) { // for shorter MEI as ALU
		if (float(csi.n_missing/csi.svlen)>0.4) {
			csi.filter = 3;
			return;
		}
		if (std::min(csi.left_anchor_length, csi.right_anchor_length) < avr_read_length_across_all * 0.8) {
			csi.filter = 4;
			return;
		}
		int min_sum = round(avr_ins_size_across_all * 0.8);
		if (csi.left_anchor_length+csi.right_anchor_length < std::min(csi.svlen/2, min_sum)) {
			csi.filter = 5;
			return;
		}
	}
	else { // for longer MEI as L1
		if (std::min(csi.left_anchor_length, csi.right_anchor_length) < avr_read_length_across_all) {
			csi.filter = 4;
			return;
		}
		if (csi.left_anchor_length+csi.right_anchor_length < avr_ins_size_across_all) {
			csi.filter = 5;
			return;
		}
	}
	csi.filter = 0;
}



void AsbSites::printHardFilterSiteVcfHeader( std::ofstream & vcf )
{
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << endl;
}

float AsbSites::_sumFloatVector( vector<float> & v )
{
	float sum = 0;
	for( int i=0; i<v.size(); i++ )
		sum += v[i];
	return sum;
}


void AsbSites::updateTrueBreakp( AsbSiteListElement & ase, int breakp )
{
	float a1 = 1 / (float)(ase.nclip + 1);
	if (ase.true_breakp != -1)
		ase.true_breakp = round( (float)ase.true_breakp * (1-a1) + (float)breakp * a1 );
	else
		ase.true_breakp = breakp;
	ase.nclip++;
}


// if clip , then set clip position
// then as avr_ins_size/3
void AsbSites::updateDiscBreakp( AsbSiteListElement & ase, SamRecord & rec )
{
	int breakp;
	Cigar * myCigar = rec.getCigarInfo();
	if (rec.getFlag() & 0x10) { // right anchor
		int clen = myCigar->getNumBeginClips();
		if (clen > MIN_CLIP/2)
			breakp = rec.get1BasedPosition();
		else
			breakp = rec.get1BasedAlignmentEnd() - avr_ins_size*0.35;
	}
	else { // left anchor
		int clen = myCigar->getNumEndClips();
		if (clen < -MIN_CLIP/2)
			breakp = rec.get1BasedAlignmentEnd();
		else
			breakp = rec.get1BasedPosition() + avr_ins_size*0.35;
	}

	float a1 = 1.0 / (float)(ase.ndisc + 1);
	if (ase.disc_breakp != -1)
		ase.disc_breakp = round( (float)ase.disc_breakp * (1-a1) + (float)breakp * a1 );
	else
		ase.disc_breakp = breakp;
	ase.ndisc++;
}


int AsbSites::_getCompSiteCount()
{
	int sum = 0;
	for( map<string, map<int, CompleteSiteInfo> >::iterator c=compSites.begin(); c!=compSites.end(); c++ ) {
		for( map<int, CompleteSiteInfo>::iterator p=c->second.begin(); p!=c->second.end(); p++ )
			sum++;
	}
	return sum;
}

void AsbSites::setCsiBreakPoint( CompleteSiteInfo & csi, AsbSiteListElement & ase )
{
	// if no clip, then imprecise
	if (ase.true_breakp == -1) {
		if (ase.disc_breakp != -1)
			csi.breakpoint = ase.disc_breakp;
		else
			csi.breakpoint = ase.raw_breakp;
		csi.imprecise = 1;
		return;
	}

	// if 2 estimates are near, then return clip estimate
	if ( abs( ase.true_breakp - ase.disc_breakp ) < avr_ins_size_across_all * 0.8 ) {
		csi.breakpoint = ase.true_breakp;
		csi.imprecise = 0;
		return;
	}

	// evaluate
	if ( ase.ndisc > ase.nclip * 10 ) {
		csi.breakpoint = ase.disc_breakp;
		csi.imprecise = 1;
	}
	else {
		csi.breakpoint = ase.true_breakp;
		csi.imprecise = 0;
	}
}


/* 
	move to discover phase

check if too many reads in a section. If so, return false
bool AsbSites::isGoodSectionDepth( const char* chr, int st, int ed, SamFile & bam, SamFileHeader & bam_header )
{
	bool section_status = bam.SetReadSection( chr, st, ed, 0 );
	if (!section_status)
		return 0;

	int lc = 0;
	int max_lc = (float)sample_depth[sample_depth.size()-1] / (float)avr_read_length * (float)(ed - st) * 10;
	SamRecord rec;
	while(bam.ReadRecord(bam_header, rec)) {
		lc++;
		if (lc > max_lc)
			return 0;
	}

	return 1;
}
*/

