#include "Cluster.h"
#include "QC.h"
#include "ProperDeckUtility.h"
#include "ssw_cpp.h"
#include "MapParameters.h"

using std::endl;

Cluster::Cluster( int mei_index, vector<string> & MEseqs )
{
	this->mei_index = mei_index;
	pMEseqs = &MEseqs;
	rClusters.resize(4);
	for( int i=0; i<4; i++ )
		rClusters[i].resize( pMEseqs->size() );
}

void Cluster::AddProper( SamRecord & sam_rec )
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
	Seqs.push_back( seq );
	int key = Seqs.size() - 1;
	string clip_str;
	if ( max_clip > 0 ) { // left clip: make right cluster
		clip_str = seq.substr( 0, max_clip );
		setEviInfoByRemap( key, clip_str, rClusters[2], 1, 0 ); // bool: boundary?, left-bound?
		string rv_str = ReverseCompSeq( clip_str );
		setEviInfoByRemap( key, rv_str, rClusters[3], 1, 0 );
	}
	else { // right clip: make left cluster
		clip_str = seq.substr( seq.length() + max_clip, -max_clip );
		setEviInfoByRemap( key, clip_str, rClusters[0], 1, 1 );
		string rv_str = ReverseCompSeq( clip_str );
		setEviInfoByRemap( key, rv_str, rClusters[1], 1, 1 );
	}
}

void Cluster::AddDisc( SamRecord & sam_rec )
{
	string seq = sam_rec.getSequence();
	Seqs.push_back( seq );
	int key = Seqs.size() - 1;
	if ( sam_rec.getFlag() & 0x20 ) { // mate(anchor) reverse: right cluster
		setEviInfoByRemap( key, seq, rClusters[2], 0, 0 );
		string rv_seq = ReverseCompSeq( seq );
		setEviInfoByRemap( key, rv_seq, rClusters[3], 0, 0 );
	}
	else { // left cluster
		setEviInfoByRemap( key, seq, rClusters[0], 0, 0 );
		string rv_seq = ReverseCompSeq( seq );
		setEviInfoByRemap( key, rv_seq, rClusters[1], 0, 0 );	
	}
}


void Cluster::Print( ofstream & out_info, ofstream & seq_info )
{
	out_info << mei_index;
// print cluster
	for( int cl=0; cl<4; cl++ ) {
		out_info << ":";
		for( vector< subCluster >::iterator vp = rClusters[cl].begin(); vp != rClusters[cl].end(); vp++ ) {
			if ( vp !=  rClusters[cl].begin() )
				out_info << ";";
			bool first_pass = 0;
			for( subCluster::iterator sub = vp->begin(); sub != vp->end(); sub++ ) {
				for( vector<EviInfo>::iterator evi = sub->second.begin(); evi != sub->second.end(); evi++ ) {
					if ( !first_pass )
						first_pass = 1;
					else
						out_info << ",";
					out_info << sub->first << "," << evi->Boundary << "," << evi->LAlign << "," << evi->RAlign << ",";
					out_info << evi->Score << "," << evi->SeqKey;
				}
			}
		}
	}
// print seq after cluster
	for( vector<string>::iterator str = Seqs.begin(); str != Seqs.end(); str++ ) {
		if ( str != Seqs.begin() )
			seq_info << ",";
		seq_info << *str;
	}
// end this line
	seq_info << endl;
	out_info << endl;
}	
	
void Cluster::setEviInfoByRemap( int key, string & seq, vector< subCluster > & evec, bool boundary, bool lbound )
{
	int SR = seq.length() * MATCH * MAP_RATIO;

	vector< subCluster >::iterator pe = evec.begin();
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
			ei.Boundary = lbound? 2 : 1;
		else
			ei.Boundary = 0;
		if ( alignment.sw_score < SR ) { // unmap
			ei.LAlign = 0;
			ei.RAlign = 0;
			ei.Score = 0;
			if ( seq.length() > 30 ) { // let's try partial map for long reads
				vector< string > vec_str;
				int part_count = seq.length() / 30 + 1;
				vec_str.resize( part_count );
				for( int i = 0; i < part_count - 1; i++ )
					vec_str[i] = seq.substr( i*30, 30 );
				vec_str[ part_count - 1 ] = seq.substr( seq.length() - 30 );
				int tail_length = part_count * 30 - seq.length();
				int map_count;
				if ( tail_length < 15 && part_count > 2 )
					map_count = part_count - 1;
				else
					map_count = part_count;
				for( int i = 0; i < map_count; i++ ) {
					StripedSmithWaterman::Aligner aligner2;
					StripedSmithWaterman::Filter filter2;
					StripedSmithWaterman::Alignment alignment2;
					aligner2.Align(seq.c_str(), pseq->c_str(), pseq->length(), filter2, & alignment2);
					if ( alignment2.sw_score >= 48 ) {
						ei.LAlign = alignment2.ref_begin;
						ei.RAlign = alignment2.ref_end;
						ei.Score = alignment2.sw_score;
						break;
					}
				}
			}
		}
		else { // mapped
			ei.LAlign = alignment.ref_begin;
			ei.RAlign = alignment.ref_end;
			ei.Score = alignment.sw_score;		
		}
		ei.SeqKey = key;	
		mp->second.push_back( ei ); // add to sub vector
	}
}
	
	