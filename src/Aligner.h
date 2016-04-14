#ifndef ALIGNER_H
#define ALIGNER_H

#include "ssw_cpp.h"
#include <string>

using std::string;

class Aligner {
public:
	Aligner( string & subj, string & query );
	int GetMapLength();
	bool IsMapped();
	int GetScore();
	int GetLeftRefPosition();
	int GetRightRefPosition();
//	float GetLogAlignProb();

private:
	string * subject_seq;
	string * query_seq;
	int subject_len;
	int query_len;

//	float logAlignProb;
	bool is_log_align_prob_set;
	int map_length;

	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;

//	static float SCORE_DENOM;
	static float MAP_RATIO;
};

#endif