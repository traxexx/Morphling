#include "AsbSite.h"
#include "Globals.h"
#include "MapParameters.h"
#include "ssw_cpp.h"
#include "Utilities.h"
#include <iostream>
#include <utility>
#include <algorithm> // max_element

using std::cout;
using std::cerr;
using std::endl;

// constructor: prepare all necessary members
AsbSite::AsbSite( vector< vector<string> > & Seqs, subCluster & sc1, subCluster & sc2, bool strand, string & ref_seq )
{
	rSeq = ref_seq;
// set default
	corrected = 0;
	consecutive_miss = -1;
	left_length = -1;
	right_length = -1;
	if ( !sc1.empty() )
		findMostLikelyCenter( sc1, Seqs );
	if ( !sc2.empty() )	
		findMostLikelyCenter( sc2, Seqs );
	setAssemblyInfo( sc1, sc2 );
}

void AsbSite::findMostLikelyCenter( subCluster & sc, vector< vector<string> > & Seqs )
{
	int center_index = 0;
	vector<subCluster> AltVec; // evi.Boundary here stores evi-vector index
	vector<subCluster> NeutVec;
	vector<int> centerScore; // score of each center
	AltVec.resize( sc.size() );
	NeutVec.resize( sc.size() );
	centerScore.resize( sc.size(), -1 );
	for( subCluster::iterator pcenter = sc.begin(); pcenter != sc.end(); pcenter++, center_index++ ) {
	// score = 0 key cannot be used as center
		bool score0 = 1;
		for( int i=0; i<(int)pcenter->second.size(); i++ ) {
			if ( pcenter->second[i].Score > 0 ) {
				score0 = 0;
				break;
			}
		}
		if ( score0 ) {
			centerScore[center_index] = 0;
			continue;
		}
	// see if need boundary
		bool lb = 0;
		bool rb = 0;
		for( int i=0; i<(int)pcenter->second.size(); i++ ) {
			if ( pcenter->second[i].Boundary == 2 )
				lb = 1;
			else if ( pcenter->second[i].Boundary == 1 )
				rb = 1;
		}
		if ( lb & rb ) {
			cerr << "ERROR: [refine-score & position] left/right boundary cannot exist simultaneously in one cluster!" << endl;
			exit(1);
		}
	// set range
		int range_st;
		int range_ed;
		if ( !lb && !rb ) { // not truncated
			range_st = pcenter->first - WIN / 2 - 20;
			range_ed = pcenter->first + WIN / 2 + 20;
		}
		else {
			if ( lb ) { // left boundary
				range_st = pcenter->first - 20;
				range_ed = pcenter->first + WIN;
			}
			else { // right boundary
				range_st = pcenter->first - WIN;
				range_ed = pcenter->first + 20;
			}
		}	
	// adjust	
		if (range_st < 0)
			range_st = 0;
		if ( range_ed > (int)rSeq.length() )
			range_ed = rSeq.length();
		if ( range_ed - range_st < 100 ) { // skip short refs (do not use this center
			centerScore[center_index] = 0;
			continue;
		}
	// sum score. remap if necessary
		int score_sum = 0;
		int read_index = 0;
		string partial_ref;
		for( subCluster::iterator pread = sc.begin(); pread != sc.end(); pread++,read_index++ ) { // each read
			if ( read_index == center_index ) {
				for( int i=0; i<(int)pread->second.size(); i++ )
					score_sum += pread->second[i].Score;
				continue;
			}
			for( int i=0; i<(int)pread->second.size(); i++ ) {
				if ( pread->second[i].LAlign < range_st || pread->second[i].RAlign > range_ed ) { // need alt mapping
					if ( partial_ref.empty() ) { // if ref not set, do it
						partial_ref = rSeq.substr( range_st, range_ed - range_st + 1 );
					}
					StripedSmithWaterman::Aligner aligner;
					StripedSmithWaterman::Filter filter;
					StripedSmithWaterman::Alignment alignment;
					int sp = pread->second[i].SampleKey;
					if ( sp >= (int)Seqs.size() ) {
						cerr << "ERROR: [refine-score & position] at " << i << "th read, sp = " << sp << " > Seqs.size()=" << Seqs.size() << endl;
						exit(1);
					}
					int qk = pread->second[i].SeqKey;
					if ( qk >= (int)Seqs[sp].size() ) {
						cerr << "ERROR: [refine-score & position] at " << i << "th read, qk=" << qk << " > Seqs[sp].size(). Something goes wrong with preAssemble!" << endl;
						exit(1);
					}
					aligner.Align( Seqs[sp][qk].c_str(), partial_ref.c_str(), partial_ref.length(), filter, & alignment);
					int SR = Seqs[sp][qk].length() * MATCH * MAP_RATIO;
					EviInfo evi;
					evi.Boundary = i;
					if ( alignment.sw_score >= SR ) { // add to alternative map
						evi.LAlign = alignment.ref_begin;
						evi.RAlign = alignment.ref_end;
						evi.Score = alignment.sw_score;
						AltVec[center_index][pread->first].push_back(evi);
						score_sum += alignment.sw_score;
					}
					else // add to neutralize map. score = 0
						NeutVec[center_index][pread->first].push_back( evi );
				}
				else
					score_sum += pread->second[i].Score;
			}
		}
		centerScore[center_index] = score_sum;
	}

// get max center
	vector<int>::iterator pmax = std::max_element( centerScore.begin(), centerScore.end() );
	int max_index = pmax - centerScore.begin();
	if ( max_index >= (int)AltVec.size() || max_index >= (int)NeutVec.size() ) {
		cerr << "ERROR: [findMostLikelyCenter] max_index = " << max_index << ", but AltMap = " << AltVec.size() << " & Neut = " << NeutVec.size() << endl;
		exit(1);
	}
	
// adjust sc
  // alt map
	for( subCluster::iterator ap = AltVec[max_index].begin(); ap != AltVec[max_index].end(); ap++ ) {
		subCluster::iterator ep = sc.find(ap->first);
		if ( ep == sc.end() ) {
			cerr << "ERROR: [findMostLikelyCenter] Can't find Alt key: " << ap->first << "in psub!" << endl;
			exit(1);
		}
		for( int i=0; i<(int)ap->second.size(); i++  ) {
			int vindex = ap->second[i].Boundary;
			ep->second[vindex].LAlign = ap->second[i].LAlign;
			ep->second[vindex].RAlign = ap->second[i].RAlign;
			ep->second[vindex].Score = ap->second[i].Score;
		}
	}
 // neutral map
	for( subCluster::iterator np = NeutVec[max_index].begin(); np != NeutVec[max_index].end(); np++ ) {
		subCluster::iterator ep = sc.find(np->first);
		if ( ep == sc.end() ) {
			cerr << "ERROR: [findMostLikelyCenter] Can't find Neut key: " << np->first << "in psub!" << endl;
			exit(1);
		}
		for( int i=0; i<(int)np->second.size(); i++  ) {
			int vindex = np->second[i].Boundary;
			ep->second[vindex].Score = 0;
		}
	}	
}

void AsbSite::setAssemblyInfo( subCluster & cluster1, subCluster & cluster2 )
{
	if ( cluster1.empty() && cluster2.empty() ) {
		this->assembled = 0;
		sample_count = 0;
		return;
	}
	this->assembled = 1;
	setSampleCount( cluster1, cluster2 );
	vector<int> basecov( rSeq.size(), 0);

// cluster 1
	for( map<int, vector<EviInfo> >::iterator mp = cluster1.begin(); mp != cluster1.end(); mp++) {
		for( vector<EviInfo>::iterator ev = mp->second.begin(); ev != mp->second.end(); ev++ ) {
			if ( ev->Score <= 0 )
				continue;
			for( int i=ev->LAlign; i<=ev->RAlign; i++ )
				basecov[i]++;
		}
	}
// cluster 2
	for( map<int, vector<EviInfo> >::iterator mp = cluster2.begin(); mp != cluster2.end(); mp++) {
		for( vector<EviInfo>::iterator ev = mp->second.begin(); ev != mp->second.end(); ev++ ) {
			if ( ev->Score <= 0 )
				continue;
			for( int i=ev->LAlign; i<=ev->RAlign; i++ )
				basecov[i]++;
		}
	}
// set info
	left_most = 0;
	right_most = (int)basecov.size() - 1;
  // left most
	for( int i=0; i<(int)basecov.size();i++ ) {
		if ( basecov[i] > 0 ) {
			left_most = i;
			break;
		}
	}
  // right most
	for( int i=(int)basecov.size()-1;i>=0 ;i-- ) {
		if ( basecov[i] > 0 ) {
			right_most = i;
			break;
		}
	}
  // missing & basecount
	basecount = calculateBaseCount( basecov );
  	missing_base = calculateMissingBase( basecov );
	
  // if missing base > 500, see #missing bases in the middle
  	if ( missing_base > 500 )
  		setMiddleBase( basecov );
}


void AsbSite::setMiddleBase( vector<int> & basecov )
{
	int max_miss, max_left, max_right;
	bool side_fixed = 0; // if both side contribute >20% dp to all
	while( !side_fixed ) { // do this iteratively until either side consists of >20% depth
		max_miss = 0;
		max_left = left_most+1; // first 0
		max_right = right_most-1; // last 0
	// find max 0 segment length
		for( int i=left_most+1; i<right_most-1; i++ ) {
			if ( basecov[i] == 0 ) {
				int left = i;
				int right = i;
				for( int j=i+1; j<(int)basecov.size(); j++ ) {
					if ( basecov[j] != 0 ) {  // summary this 0 segment
						right = j-1;
						i = j+1<right_most-1 ? j+1 : right_most-2;
						break;
					}
				}
				int current_miss = right - left + 1;
				if ( current_miss > max_miss ) {
					max_miss = current_miss;
					max_left = left;
					max_right = right;
				}
			}
		}
	// correct length: 	 if >80% coverage fall on one side anchor, discard the other side and use this side as MEI length.
		int ldp = 0;
		for( int i=left_most; i<max_left; i++ )
			ldp += basecov[i];
		int rdp = 0;
		for( int i=max_right+1; i<=right_most; i++ )
			rdp += basecov[i];
		int dpsum = ldp + rdp;
		if ( (float)ldp/dpsum < 0.2 || (float)rdp/dpsum < 0.2 ) {
			if ( (float)ldp/dpsum < 0.2 )
				left_most = max_right + 1;
			else
				right_most = max_left - 1;
			corrected++;
			int current_miss_base = calculateMissingBase( basecov );
			if ( current_miss_base <= 500 )
				side_fixed = 1;
		}
		else
			side_fixed = 1;
	}

	int current_miss_base = calculateMissingBase( basecov );
	this->missing_base = current_miss_base;
	this->basecount = calculateBaseCount( basecov );
	if ( current_miss_base > 500 ) {
		consecutive_miss = max_miss;
		left_length = max_left - left_most;
		right_length = right_most - max_right;
	}
	else {
		consecutive_miss = -1;
		left_length = -1;
		right_length = -1;	
	}
}

int AsbSite::calculateMissingBase( vector<int> & basecov )
{
	map<int, bool> miss_map;
	for( int i=left_most+1; i<right_most; i++ ) {
  		if( basecov[i] == 0 )
  			miss_map[i];
  	}
  	int mb = miss_map.size();
  	return mb;
}

int AsbSite::calculateBaseCount( vector<int> & basecov )
{
	int sum = 0;;
  	for( int i=left_most; i<=right_most; i++ )
  		sum += basecov[i];
  	return sum;
}

// set asb sample count
// if sample key exist, add to sample count if not added
// if exists score > 0, add to valid sample count if not added
void AsbSite::setSampleCount( subCluster & cluster1, subCluster & cluster2 )
{
	bool scount[NSAMPLE]; // sample with candidate reads
	bool vs[NSAMPLE]; // sample with score > 0
	
	for( int i=0; i<NSAMPLE; i++ ) {
		scount[i] = 0;
		vs[i] = 0;
	}
	
	for( subCluster::iterator p = cluster1.begin(); p != cluster1.end(); p++ ) {
		for( vector<EviInfo>::iterator ev = p->second.begin(); ev != p->second.end(); ev++ ) {
			if ( !scount[ev->SampleKey] )
				scount[ev->SampleKey] = 1;
			if ( ev->Score > 0 ) {
				if ( !vs[ev->SampleKey] )
					vs[ev->SampleKey] = 1;
			}
		}
	}
	for( subCluster::iterator p = cluster2.begin(); p != cluster2.end(); p++ ) {
		for( vector<EviInfo>::iterator ev = p->second.begin(); ev != p->second.end(); ev++ ) {
			if ( !scount[ev->SampleKey] )
				scount[ev->SampleKey] = 1;
			if ( ev->Score > 0 ) {
				if ( !vs[ev->SampleKey] )
					vs[ev->SampleKey] = 1;
			}
		}	
	}
	
	int sum = 0;
	for( int i=0; i<NSAMPLE; i++ ) {
		if ( scount[i] )
			sum++;
	}
	sample_count = sum;
	
	sum = 0;
	for( int i=0; i<NSAMPLE; i++ ) {
		if ( vs[i] )
			sum++;
	}
	valid_sample_count = sum;
}


// get methods
bool AsbSite::IsAssembled()
{
	return this->assembled;
}

int AsbSite::GetSVlength()
{
	return right_most - left_most;
}

float AsbSite::GetSVdepth()
{
	if ( basecount <= 0 )
		return 0;
	else
		return ((float)basecount / ( right_most - left_most - missing_base ) / sample_count);
}

int AsbSite::GetMissingBaseCount()
{
	return this->missing_base;
}

int AsbSite::GetCorrection()
{
	return this->corrected;
}

int AsbSite::GetLeftMost()
{
	return this->left_most;
}

int AsbSite::GetRightMost()
{
	return this->right_most;
}

int AsbSite::GetConsecutiveMiss()
{
	return this->consecutive_miss;
}

int AsbSite::GetLeftLength()
{
	return this->left_length;
}

int AsbSite::GetRightLength()
{
	return this->right_length;
}

int AsbSite::GetSampleCount()
{
	return sample_count;
}

int AsbSite::GetValidSampleCount()
{
	return valid_sample_count;
}


































