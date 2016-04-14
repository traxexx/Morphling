#include "Aligner.h"
#include "Cigar.h"
#include "seqUtility.h"
#include <math.h>
#include "Globals.h"

//float Aligner::SCORE_DENOM = 3.4; // 1/lambda in e value
float Aligner::MAP_RATIO = 0.6;

Aligner::Aligner( string & subj, string & query )
{
	if (subj.length() < query.length())
		std::cerr << "Warning: [Aligner] subject sequence length < query seq length! Are you using the right order?\n" << std::endl;

	subject_seq = &subj;
	query_seq = &query;
	subject_len = subject_seq->length();
	query_len = query_seq->length();

	is_log_align_prob_set = 0;

	aligner.Align( query_seq->c_str(), subject_seq->c_str(), subject_len, filter, &alignment );
	map_length = query_seq->length() - GetTotalClipLength( alignment.cigar_string );
}

int Aligner::GetMapLength()
{
	return map_length;
}

bool Aligner::IsMapped()
{
	if (map_length<MIN_CLIP*0.8)
		return 0;
	int score_thred = map_length * MATCH * MAP_RATIO;
	if ( alignment.sw_score < score_thred )
		return 0;
	else
		return 1;
}

int Aligner::GetScore()
{
	return alignment.sw_score;
}

int Aligner::GetLeftRefPosition()
{
	return alignment.ref_begin;
}

int Aligner::GetRightRefPosition()
{
	return alignment.ref_end;
}

/*
float Aligner::GetLogAlignProb() {
	if (is_log_align_prob_set)
		return logAlignProb;

	float eval = log(subject_len) + log(query_len) + (float)alignment.sw_score / SCORE_DENOM;
	if (eval + log(2) >= 0)
		logAlignProb = 1;
	else
		logAlignProb = eval;
	is_log_align_prob_set = 1;
	return logAlignProb;
}
*/

