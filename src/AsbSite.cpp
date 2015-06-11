#include "AsbSite.h"
#include "Globals.h"
#include "ReadMapGlobals.h"
#include "MapParameters.h"
#include "ProperDeckUtility.h"
#include "QC.h"
#include "ssw_cpp.h"
#include "Utilities.h"
#include <iostream>
#include <utility>
#include <algorithm> // max_element

using std::cout;
using std::cerr;
using std::endl;

// constructor: prepare all necessary members
AsbSite::AsbSite( string & chr, int position, bool* gtList, vector<SamFile> & BamFiles, vector<SamFileHeader> & BamFileHeaders, vector<string> & MEseqs )
{
	Chr = chr;
	Position = position;
	GtList = gtList;
	pBamFiles = &BamFiles;
	pBamFileHeaders = &BamFileHeaders;
	pMEseqs = &MEseqs;
	assembled = 0;
	Clusters.resize( 4 );
	for( int i=0; i<4; i++ )
		Clusters[i].resize( pMEseqs->size() );
	sample_count = 0;
	subtype = -1;
	plusstr = 0;
	left_most = -1;
	right_most = -1;
	missing_base = -1;
	basecount = -1;
}


// do assembly
void AsbSite::Assembly()
{
// get cluster info from bams
	setClusterInfoFromBams();
	if ( sample_count == 0 ) {
		cerr << "Warning: no available samples at site: " << Position << ". Skip assembly" << endl;
		return;
	}

// adjust score with proper remap
	findMostLikelyCluster();
}

// get methods
bool AsbSite::IsAssembled()
{
	return this->assembled;
}

int AsbSite::GetPosition()
{
	return this->Position;
}

int AsbSite::GetSubtype()
{
	return subtype;
}

int AsbSite::GetSVlength()
{
	return right_most - left_most;
}

float AsbSite::GetSVdepth()
{
	return ((float)basecount / ( right_most - left_most ) / sample_count);
}

int AsbSite::GetMissingBaseCount()
{
	return this->missing_base;
}

int AsbSite::GetLeftMost()
{
	return this->left_most;
}

int AsbSite::GetRightMost()
{
	return this->right_most;
}

bool AsbSite::GetStrand()
{
	return this->plusstr;
}

int AsbSite::GetSampleCount()
{
	return sample_count;
}


/*** add read info from cluster ****/

void AsbSite::setClusterInfoFromBams()
{
// add from each bam
	int cluster_start = Position - WIN / 2;
	int cluster_end = Position + WIN / 2;
	int ns = 0;
	for( int sp=0; sp < NSAMPLE; sp++ ) {
		if ( GtList[sp] ) { // this bam contain MEI info
			add2ClusterFromSingleBam( sp, cluster_start, cluster_end );
			ns++;
		}
	}
	sample_count = ns;
}

void AsbSite::add2ClusterFromSingleBam( int sp, int cluster_start, int cluster_end )
{
 	vector<SamFile>::iterator pBam = pBamFiles->begin(); // point to opened bam files
 	pBam += sp;
 	vector<SamFileHeader>::iterator pHeader = pBamFileHeaders->begin();
 	pHeader += sp;
// sanity check
	if ( !pBam->IsOpen() ) {
		cerr << "ERROR: [add2ClusterFromSingleBam] index = " << sp << ", bamFile is not open!" << endl;
		exit(1);
	}
	bool section_status = pBam->SetReadSection( Chr.c_str(), cluster_start, cluster_end );
	if ( !section_status ) {
		cerr << "ERROR: [add2ClusterFromSingleBam] cannot set section: chr" << Chr << " " << cluster_start << "-" << cluster_end << endl;
		exit(1);
	}

// find disc or clip reads
	vector<DiscInfo> discVec;
	SamRecord sam_rec;
	while( pBam->ReadRecord( *pHeader, sam_rec ) ) {
		if ( !PassQC(sam_rec) )
			continue;
		if ( sam_rec.getFlag() & 0x2 ) // proper
			addProper2Cluster( sam_rec );
		else { // save disc info
			DiscInfo di;
			di.Chr = sam_rec.getMateReferenceName();
			di.Position = sam_rec.get1BasedMatePosition();
			di.MatePosition = sam_rec.get1BasedPosition();
			discVec.push_back( di );
		}
	}
	
// add from disc info
	for( vector<DiscInfo>::iterator pd = discVec.begin(); pd != discVec.end(); pd++ ) {
		bool section_status = pBam->SetReadSection( pd->Chr.c_str(), pd->Position, pd->Position + 1 );
		if ( !section_status )
			cerr << "Warning: [add2ClusterFromSingleBam] Unable to set read section. Skip this DiscInfo record!" << endl;
		bool got_mate = 0;
		while(pBam->ReadRecord(*pHeader, sam_rec)) {
			if ( !PassQC(sam_rec) )
				continue;
			if ( pd->Position == sam_rec.get1BasedPosition() &&
				 pd->MatePosition == sam_rec.get1BasedMatePosition() &&
				 Chr.compare( sam_rec.getMateReferenceName() ) == 0 &&
				 pd->Chr.compare( sam_rec.getReferenceName() ) == 0 ) {
				addDisc2Cluster( sam_rec );
				got_mate = 1;
				break;
			}
		}
//		if ( !got_mate )
//			cerr << "Warning: [add2ClusterFromSingleBam] Unable to find its mate. Skip this DiscInfo record!" << endl;
	}
}

void AsbSite::addProper2Cluster( SamRecord & sam_rec )
{
// sanity check: skip reads with no clip
	int max_clip = getMaxClipLen( sam_rec );
	if ( max_clip == 0 )
		return;
	if ( max_clip > 0 )
		if ( max_clip < MinClip )
			return;
	if ( max_clip < 0 )
		if ( -max_clip < MinClip )
			return;

// re-map
	string seq = sam_rec.getSequence();
	string clip_str;
	if ( max_clip > 0 ) { // left clip: make right cluster
		clip_str = seq.substr( 0, max_clip );
		setEviInfoByRemap( clip_str, Clusters[2], 1, 0 ); // bool: boundary?, left-bound?
		string rv_str = ReverseCompSeq( clip_str );
		setEviInfoByRemap( rv_str, Clusters[3], 1, 0 );
	}
	else { // right clip: make left cluster
		clip_str = seq.substr( seq.length() + max_clip, -max_clip );
		setEviInfoByRemap( clip_str, Clusters[0], 1, 1 );
		string rv_str = ReverseCompSeq( clip_str );
		setEviInfoByRemap( rv_str, Clusters[1], 1, 1 );
	}
}

void AsbSite::addDisc2Cluster( SamRecord & sam_rec )
{
	string seq = sam_rec.getSequence();
	if ( sam_rec.getFlag() & 0x20 ) { // mate(anchor) reverse: right cluster
		setEviInfoByRemap( seq, Clusters[2], 0, 0 );
		string rv_seq = ReverseCompSeq( seq );
		setEviInfoByRemap( rv_seq, Clusters[3], 0, 0 );
	}
	else { // left cluster
		setEviInfoByRemap( seq, Clusters[0], 0, 0 );
		string rv_seq = ReverseCompSeq( seq );
		setEviInfoByRemap( rv_seq, Clusters[1], 0, 0 );	
	}
}

void AsbSite::setEviInfoByRemap( string & seq, vector< map<int, vector<EviInfo> > > & evec, bool boundary, bool lbound )
{
	int SR = seq.length() * MATCH * MAP_RATIO;

	vector< map<int, vector<EviInfo> > >::iterator pe = evec.begin();
	for( vector<string>::iterator pseq = pMEseqs->begin(); pseq != pMEseqs->end(); pseq++, pe++ ) {
		StripedSmithWaterman::Aligner aligner;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment alignment;
		aligner.Align(seq.c_str(), pseq->c_str(), pseq->length(), filter, & alignment);
		int center = alignment.ref_begin + seq.length() / 2;
		map<int, vector<EviInfo> >::iterator mp = pe->find( center );
		if ( mp == pe->end() ) { // add key if not exist
			(*pe)[center].clear();
			mp = pe->find( center );
		}
		EviInfo ei;
		if ( boundary )
			ei.Boundary = lbound? -alignment.ref_begin : alignment.ref_end;
		else
			ei.Boundary = 0;
		if ( alignment.sw_score < SR ) {
			ei.LAlign = 0;
			ei.RAlign = 0;
			ei.Score = 0;
		}
		else {
			ei.LAlign = alignment.ref_begin;
			ei.RAlign = alignment.ref_end;
			ei.Score = alignment.sw_score;		
		}
		ei.Seq = seq;	
		mp->second.push_back( ei ); // add to sub vector
	}
}


/*** calculate scores ***/

void AsbSite::findMostLikelyCluster() // 0->2, 1->3
{
// calculate scores
	vector< vector<int> > SumScoreVec;
	SumScoreVec.resize(4);
	for( int cl=0; cl<4; cl++ ) {
		SumScoreVec[cl].resize( pMEseqs->size(), 0 );
		int sub_index = 0;
		for( vector< map<int, vector<EviInfo> > >::iterator psub = Clusters[cl].begin(); psub != Clusters[cl].end(); psub++, sub_index++ ) { // subtype
			if ( psub->size() == 0 ) {
				SumScoreVec[cl][sub_index] = 0;
				continue;
			}
			int center_index = 0;
			vector< map<int, vector<EviInfo> > > AltMapVec; // evi.Boundary here stores evi-vector index
			vector< map<int, vector<EviInfo> > > NeutMapVec;
			AltMapVec.resize( psub->size() );
			NeutMapVec.resize( psub->size() );
			vector< map<int, vector<EviInfo> > >::iterator altp = AltMapVec.begin();
			vector< map<int, vector<EviInfo> > >::iterator neutp = NeutMapVec.begin();
			vector<int> centerScore; // score of each center
			centerScore.resize( psub->size() );
			for( map<int, vector<EviInfo> >::iterator pcenter = psub->begin(); pcenter != psub->end(); pcenter++, center_index++, altp++, neutp++ ) { // center
				if ( altp == AltMapVec.end() || neutp == NeutMapVec.end() ) {
					cerr << "ERROR: At position " << this->Position << ", AltMap / NeutMap reached end!" << endl;
					exit(1);
				}
			// see if we need boundary
				bool lb = 0;
				bool rb = 0;
				for( int i=0; i<(int)pcenter->second.size(); i++ ) {
					if ( pcenter->second[i].Boundary > 0 )
						lb = 1;
					else if ( pcenter->second[i].Boundary < 0 )
						rb = 1;
				}
				if ( lb & rb ) {
					cerr << "ERROR: [findMostLikelyCluster] At " << this->Position << ", left/right boundary cannot exist simultaneously in one cluster!" << endl;
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
				if ( range_ed > (int)(*pMEseqs)[sub_index].length() )
					range_ed = (int)(*pMEseqs)[sub_index].length();
				if ( range_ed - range_st < 100 ) { // skip short refs (do not use this center
					centerScore[center_index] = 0;
					continue;
				}
			// sum score. remap if necessary
				int score_sum = 0;
				int read_index = 0;
				string partial_ref;
				for( map<int, vector<EviInfo> >::iterator pread = psub->begin(); pread != psub->end(); pread++,read_index++ ) { // each read
					if ( read_index == center_index )
						continue;
					for( int i=0; i<(int)pread->second.size(); i++ ) {
						if ( pread->second[i].LAlign < range_st || pread->second[i].RAlign > range_ed ) { // need alt mapping
							if ( partial_ref.empty() ) { // if ref not set, do it
								partial_ref = (*pMEseqs)[sub_index].substr( range_st, range_ed - range_st + 1 );
							}
							StripedSmithWaterman::Aligner aligner;
							StripedSmithWaterman::Filter filter;
							StripedSmithWaterman::Alignment alignment;
							aligner.Align(pread->second[i].Seq.c_str(), (*pMEseqs)[sub_index].c_str(), partial_ref.length(), filter, & alignment);
							int SR = pread->second[i].Seq.length() * MATCH * MAP_RATIO;
							EviInfo evi;
							evi.Boundary = i;
							if ( alignment.sw_score < SR ) { // add to alternative map
								evi.LAlign = alignment.ref_begin;
								evi.RAlign = alignment.ref_end;
								evi.Score = alignment.sw_score;
								(*altp)[pread->first].push_back(evi);
								score_sum += alignment.sw_score;
							}
							else // add to neutralize map. score = 0
								(*neutp)[pread->first].push_back( evi );
						}
						else
							score_sum += pread->second[i].Score;
					}
				}
				centerScore[center_index] = score_sum;
			}
		// find score sum for each center
			vector<int>::iterator pmax = std::max_element( centerScore.begin(), centerScore.end() );
			int max_index = pmax - centerScore.begin();
			if ( max_index >= (int)AltMapVec.size() || max_index >= (int)NeutMapVec.size() ) {
				cerr << "ERROR: at " << this->Position << ", max_index = " << max_index << ", but AltMap = " << AltMapVec.size() << " & Neut = " << NeutMapVec.size() << endl;
				exit(1);
			}
		// adjust evi-info
		  // alt map
			for( map<int, vector<EviInfo> >::iterator ap = AltMapVec[max_index].begin(); ap != AltMapVec[max_index].end(); ap++ ) {
				map<int, vector<EviInfo> >::iterator ep = psub->find(ap->first);
				if ( ep == psub->end() ) {
					cerr << "ERROR: [findMostLikelyCluster] at " << this->Position << ", can't find key: " << ap->first << "in psub!" << endl;
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
			for( map<int, vector<EviInfo> >::iterator np = NeutMapVec[max_index].begin(); np != NeutMapVec[max_index].end(); np++ ) {
				map<int, vector<EviInfo> >::iterator ep = psub->find(np->first);
				if ( ep == psub->end() ) {
					cerr << "ERROR: [findMostLikelyCluster] at " << this->Position << "can't find key: " << np->first << "in psub!" << endl;
					exit(1);
				}
				for( int i=0; i<(int)np->second.size(); i++  ) {
					int vindex = np->second[i].Boundary;
					ep->second[vindex].LAlign = 0;
					ep->second[vindex].RAlign = 0;
					ep->second[vindex].Score = 0;
				}
			}
			SumScoreVec[cl][sub_index] = *pmax;
		}
	}
	
// find most likely cluster by 2-end score
	vector<int> plus_strand;
	vector<int> minus_strand;
	plus_strand.resize( (int)SumScoreVec[0].size() );
	minus_strand.resize( (int)SumScoreVec[0].size() );
	for( int i=0; i<(int)SumScoreVec[0].size(); i++ ) {
		plus_strand[i] = SumScoreVec[0][i] + SumScoreVec[2][i];
		minus_strand[i] = SumScoreVec[1][i] + SumScoreVec[3][i];
	}
	vector<int>::iterator pmax_plus = std::max_element( plus_strand.begin(), plus_strand.end() );
	vector<int>::iterator pmax_minus = std::max_element( minus_strand.begin(), minus_strand.end() );
	if ( *pmax_plus == *pmax_minus ) {
		if ( *pmax_plus == 0 )
			return;
		else
		cerr << "Warning: [findMostLikelyCluster] +/- cluster has same SW score = " << *pmax_plus << ", use (+) strand..." << endl;
	}
	if ( *pmax_plus >= *pmax_minus ) {
		this->subtype = pmax_plus - plus_strand.begin();
		this->plusstr = 1;
		setAssemblyInfo( Clusters[0][subtype], Clusters[2][subtype] );
	}
	else {
		this->subtype = pmax_minus - minus_strand.begin();
		this->plusstr = 0;
		setAssemblyInfo( Clusters[1][subtype], Clusters[3][subtype] );
	}
}


void AsbSite::setAssemblyInfo( map<int, vector<EviInfo> > & cluster1, map<int, vector<EviInfo> > & cluster2 )
{
	this->assembled = 1;
	vector<int> basecov( (*pMEseqs)[subtype].size(), 0);

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
	missing_base = 0;
	basecount = getSumOfVector( basecov );
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
  // missing basecount
  	map<int, bool> miss_map;
  	for( int i=left_most; i<=right_most; i++ ) {
  		if( basecov[i] == 0 )
  			miss_map[i];
  	}
  	missing_base = miss_map.size();
}

string ReverseCompSeq( string & seq ) {
	string rv;
	rv.resize( seq.size() );
	for( int i=0; i<(int)seq.size(); i++ )
		rv[i] = GetCompNt( seq[i] );
	return rv;
}

char GetCompNt( char nt)
{
	char c;
	switch( nt ) {
		case('A'):
			c = 'T'; break;
        case ('T'):
			c = 'A'; break;
		case ('C'):
			c = 'G'; break;
		case ('G'):
			c = 'C'; break;
		case ('a'):
			c = 'T'; break;
		case ('t'):
			c = 'A'; break;
		case ('c'):
			 c = 'G'; break;
		case ('g'):
			c = 'C'; break;
		default:
			c = 'N'; break; 
	}
	return c;
}






































