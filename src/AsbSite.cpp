#include "AsbSite.h"


// constructor
AbsSite()
{

}

// do assembly
void Assembly()
{
// get cluster info from bams
	setClusterInfoFromBams();

// adjust score with proper remap
	findMostLikelyCluster();
}


/*** add read info from cluster ****/

void setClusterInfoFromBams()
{
// add from each bam
	int cluster_start = Position - WIN / 2;
	int cluster_end = Position + WIN / 2;
	for( int sp=0; sp < NSAMPLE; sp++ ) {
		if ( GtList[sp] )
			add2ClusterFromSingleBam( bFiles[sp], cluster_start, cluster_end );
	}
}

void add2ClusterFromSingleBam( SamFile & samFile, SamFileHeader & samHeader, string chr, int cluster_start, int cluster_end )
{
// sanity check
	if ( samFile.IsOpen() ) {
		cerr << "ERROR: [add2ClusterFromSingleBam] samFile is not open!" << endl;
		exit(1):
	}
	bool section status = samFile.SetReadSection( chr, cluster_start, cluster_end );
	if ( !section_status ) {
		cerr << "ERROR: [add2ClusterFromSingleBam] cannot set section: chr" << chr << " " << cluster_start << "-" << cluster_end << endl;
		exit(1);
	}

// find disc or clip reads
	vector<DiscInfo> discVec;
	SamRecord sam_rec;
	while( samFile.ReadRecord( samHeader, sam_rec ) ) {
		if ( !PassQC(sam_rec) )
			continue;
		if ( sam_rec.getFlag() & 0x2 ) // proper
			addProperToCluster( sam_rec );
		else { // save disc info
			DiscInfo di;
			di.Chr = sam_rec.GetMateReferenceName();
			di.MateChr = sam_rec.GetReferenceName();
			di.Position = sam_rec.Get1BasedMatePosition();
			di.MatePosition = sam_rec.Get1BasedPosition();
			discVec.push_back( di );
		}
	}
	
// add from disc info
	for( vector<DiscInfo>::iterator pd = discVec.begin(); pd != discVec.end(); pd++ ) {
		bool section_status = samFile.SetReadSection( pd->Chr, pd->Position, pd->Position + 1 );
		if ( !section_status )
			cerr << "Warning: [add2ClusterFromSingleBam] Unable to set read section. Skip this DiscInfo record!" << endl;
		bool got_mate = 0;
		while(samIn.ReadRecord(samHeader, sam_rec)) {
			if ( !PassQC(sam_rec) )
				continue;
			if ( di->Position == sam_rec.get1BasedPosition() &&
				 di->MatePosition == sam_rec.get1BasedMatePosition() &&
				 di->Mate.compare( sam_rec.getMateReferenceName() )) {
				addDiscToCluster( sam_rec );
				got_mate = 1;
				break;
			}
		}
		if ( !got_mate )
			cerr << "Warning: [add2ClusterFromSingleBam] Unable to find its mate. Skip this DiscInfo record!" << endl;
	}
}

void addProperToCluster( SamRecord & sam_rec )
{
// sanity check
	int max_clip = getMaxClipLen( sam_rec ):
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
	if ( max_clip > 0 ) { // left clip
		clip_str = seq.substr( 0, max_clip );
		setEviInfoByRemap( clip_str, Cluster[0], 1, 0 );
		string rv_str = ReverseCompSeq( clip_str );
		setEviInfoByRemap( rv_str, Cluster[1], 1, 0 );
	}
	else {
		clip_str = seq.substr( seq.length() + max_clip, -max_clip );
		setEviInfoByRemap( clip_str, Cluster[2], 1, 1 );
		string rv_str = ReverseCompSeq( clip_str );
		setEviInfoByRemap( rv_str, Cluster[3], 1, 1 );
	}
}

void addDiscToCluster( SamRecord & sam_rec )
{
	string seq = sam_rec.getSequence();
	if ( sam_rec.getFlag() & 0x20 ) { // mate(anchor) reverse: right cluster
		setEviInfoByRemap( seq, Cluster[2], 0, 0 );
		string rv_seq = ReverseCompSeq( seq );
		setEviInfoByRemap( rv_seq, Cluster[3], 0, 0 );
	}
	else { // left cluster
		setEviInfoByRemap( seq, Cluster[0], 0, 0 );
		string rv_seq = ReverseCompSeq( seq );
		setEviInfoByRemap( rv_seq, Cluster[1], 0, 0 );	
	}
}

void setEviInfoByRemap( string & seq, vector< map<int, vector<EviInfo> > > & evec, bool boundary, bool lbound )
{
	int SR = seq.length() * MATCH * MAP_RATIO;

	vector< map<int, vector<EviInfo> > >::iterator pe = evec.begin();
	for( map<string, string>::iterator pseq = SeqH.begin(); pseq != SeqH.end(); pseq++, pe++ ) {
		StripedSmithWaterman::Aligner aligner;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment alignment;
		aligner.Align(seq.c_str(), pseq->second.c_str(), pseq->length(), filter, & alignment);
		int center = alignment.ref_begin + seq.length() / 2;
		map<int, EviInfo>::iterator mp = (*pe).find( center );
		if ( mp == pe->end() ) { // add key if not exist
			(*pe)[center].clear();
			mp = (*pe).find( center );
		}
		EviInfo ei;
		if ( boundary )
			ei.Boundary = lbound? -alignment.ref_begin : alignment.ref_end;
		else
			ei.Boundary = 0;
		if ( alignment.sw_score < SR ) {
			ei.LAlign = 0;
			ei.Length = 0;
			ei.Score = 0;
		}
		else {
			ei.LAlign = alinment.ref_begin;
			ei.Length = seq.length();
			ei.Score = alignment.score;		
		}
		ei.Seq = seq;	
		mp->push_back( ei );
	}
}


/*** calculate scores ***/

void findMostLikelyCluster() // 0->2, 1->3
{
// calculate scores
	vector< vector<int> > SumScoreVec;
	for( cl=0; cl<3; cl++ ) {
		SumScoreVec[cl].resize( Clusters[cl].size() );
		vector< map<int, int> >::iterator altp = AltMapVec[cl].begin();
		vector< map<int, int> >::iterator neutp = NeutMapVec[cl].begin();
		altp += cl;
		neutp += cl;
		int sub_index = 0;
		for( vector< map<int, EviInfo> >::iterator psub = Clusters[cl].begin(); psub != Clusters[cl].end(); psub++, sub_index++ ) { // subtype
			int center_index = 0;
			vector< map<int, int> > AltMapVec;
			vector< map<int, int> > NeutMapVec;
			AltMapVec.resize( Clusters[cl].size() );
			NeutMapVec.resize( Clusters[cl].size() );
			vector< map<int, EviInfo> >::iterator altp = AltMapVec.begin();
			vector< map<int, EviInfo> >::iterator neutp = NeutMapVec.begin();
			for( map<int, vector<EviInfo> >::iterator pcenter = psub->begin(); pcenter != psub->end(); pcenter++, center_index++, altp++, neutp++ ) { // center
				vector<int> centerScore; // score of each center
				centerScore.resize( psub->size() );
			// see if we need boundary
				bool lb = 0;
				bool rb = 0;
				for( int i=0; i<(int)pcenter->second.size(); i++ ) {
					if ( pcenter->second[i].Boundary > 0 )
						lb = 1;
					else ( pcenter->second[i].Boundary < 0 )
						rb = 1;
				}
				if ( lb & rb ) {
					cerr << "ERROR: [findMostLikelyCluster] left/right boundary cannot exist simultaneously in one cluster!" << endl;
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
						range_st = pcenter->first + 20;
						range_ed = pcenter->first - WIN;
					}
				}
				string partial_ref;
			// sum score. remap if necessary
				int score_sum 0;
				for( int i=0; i<(int)pread->second.size(); i++ )
					score_sum += pread->second[i].Score;
				int read_index = 0;
				for( map<int, vector<EviInfo> >::iterator pread = psub->begin(); pread != psub->end(); pread++,read_index++ ) { // each read
					if ( read_index = center_index )
						continue;
					for( int i=0; i<(int)pread->second.size(); i++ ) {
						if ( pread->second[i].LAlign < range_st || pread->second[i].RAlign > range_ed ) { // need alt mapping
							if ( partial_ref.empty() ) { // if ref not set, do it
								int adjusted_st = range_st > 0 ? range_st : 0;
								int adjusted_ed < SeqH[MeiType][sub_index].length() ? range_ed - range_st : SeqH[MeiType][sub_index].length() - range_st;
								partial_ref = SeqH[MeiType][sub_index].substr( adjusted_st, adjusted_range_ed );
							}
							StripedSmithWaterman::Aligner aligner;
							StripedSmithWaterman::Filter filter;
							StripedSmithWaterman::Alignment alignment;
							aligner.Align(pread->second[i].Seq.c_str(), SeqH[MeiType][sub_index].c_str(), partial_ref.length(), filter, & alignment);
							int SR = pread->second.Seq.length() * MATCH * MAP_RATIO;
							EviInfo evi;
							evi.Boundary = i;
							if ( alignment.sw_score < SR ) { // add to alternative map
								evi.LAlign = alignment.ref_begin;
								evi.RAlign = alignment.ref_end;
								evi.Score = alignment.sw_score;
								altp[read_index].push_back(evi);
							}
							else // add to neutralize map
								neutp[read_index].push_back( evi );
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
		// adjust evi-info
		  // alt map
			for( map<int, EviInfo >::iterator ap = AltMapVec[max_index].begin(); ap != AltMapVec[max_index].end(); ap++ ) {
				map<int, EviInfo>::iterator ep = psub->find(ap->first);
				if ( ep == psub->end() ) {
					cerr << "ERROR: [findMostLikelyCluster] can't find key: " << ap->fist << "in psub!" << endl;
					exit(1);
				}
				int vindex = ap->second.Boundary;
				ep->second[vindex].LAlign = ap->second.LAlign;
				ep->second[vindex].RAlign = ap->second.RAlign;
				ep->second[vindex].Score = ap->second.Score;
			}
		  // neutral map
				for( map<int, EviInfo >::iterator ap = NeutMapVec[max_index].begin(); ap != NeutMapVec[max_index].end(); ap++ ) {
				map<int, EviInfo>::iterator ep = psub->find(ap->first);
				if ( ep == psub->end() ) {
					cerr << "ERROR: [findMostLikelyCluster] can't find key: " << ap->fist << "in psub!" << endl;
					exit(1);
				}
				ep->second[vindex].LAlign = 0;
				ep->second[vindex].RAlign = 0;
				ep->second[vindex].Score = 0;
			}
			SumScoreVec[cl][sub_index] = *pmax;
		}
	}
	
// find most likely cluster by 2-end score
	plus_strand.resize( (int)SumScoreVec[0].size() );
	minus_strand.resize( (int)SumScoreVec[0].size() );
	vector<int> plus_strand;
	for( int i=0; i<(int)SumScoreVec[0].size(); i++ ) {
		plus_strand[i] = SumScoreVec[0][i] + SumScoreVec[2][i];
		minus_strand[i] = SumScoreVec[1][i] + SumScoreVec[3][i];
	}
	vector<int>::iteratpr pmax_plus = std::max_element( plus_strand.begin(), plus_strand.end() );
	vector<int>::iteratpr pmax_minus = std::max_element( minus_strand.begin(), minus_strand.end() );
	if ( *pmax_plus == *pmax_minus )
		cerr << "Warning: [findMostLikelyCluster] +/- cluster has same SW score = " << *pmax_plus << ", use (+) strand..." << endl;
	int max_subtype_index;
	if ( *pmax_plus >= *pmax_minus ) {
		max_subtype_index = pmax_plus - plus_strand.begin();
		setAssemblyInfo( max_subtype_index, Clusters[0], Cluster[2] );
	}
	else {
		max_subtype_index = pmax_minus - minus_strand.begin();
		setAssemblyInfo( max_subtype_index, Clusters[1], Cluster[3] );
	}
}


void setAssemblyInfo( int max_subtype_index, vector< map<int, vector<EviInfo> > & cluster1, vector< map<int, vector<EviInfo> > & cluster2 )
{
	vector<int> basecov( (int)SeqH[max_subtype_index].size(), 0);
// cluster 1
	for( vector< map<int, vector<EviInfo> >::iterator pv = cluster1.begin(); pv != cluster1.end(); pv++ ) {
		for( vector<EviInfo>::iterator ev = pv->second.begin(); ev != pv->second.end(); ev++ ) {
			for( int i=ev->LAlign; i<=ev->RAlign; i++ )
				basecov[i]++;
		}
	}
// cluster 2
	for( vector< map<int, vector<EviInfo> >::iterator pv = cluster1.begin(); pv != cluster1.end(); pv++ ) {
		for( vector<EviInfo>::iterator ev = pv->second.begin(); ev != pv->second.end(); ev++ ) {
			for( int i=ev->LAlign; i<=ev->RAlign; i++ )
				basecov[i]++;
		}
	}
// set info
	left_most = 0;
	right_most = (int)basecov.size();
	missing_basecount = 0;
	basecount = GetVecSum( basecov );
  // left most
	for( int i=0; i<(int)basecov.size();i++ ) {
		if ( basecov[i] > 0 ) {
			left_most = i;
			break;
		}
	}
  // right most
	for( int i=(int)basecov.size();i>=0 ;i-- ) {
		if ( basecov[i] > 0 ) {
			right_most = i;
			break;
		}
	}
  // missing basecount
  	map<int, bool> miss_map;
  	for( int i=left_most; i<=right_most; i++ ) {
  		if( baseconv[i] == 0 )
  			miss_map[i];
  	}
  	missing_base = miss_map.size();
}


void Print( ofstream & out_vcf )
{
	int svlen = right_most - left_most;
	float svcov = (float)basecount / svlen;
	out_vcf << str_before_format;
	out_vcf << "SUB=" << SeqName[subtype] << ";SVLEN=" << svlen << ";SVCOV=" << setprecision(2) << svcov << ";MPOS=" << left_most << "," << right_most << ";MISSING=" << missing_base;
	out_vcf << "\t" << str_after_format << endl;
}










































