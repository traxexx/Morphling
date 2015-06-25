#ifndef CONSENSUSVCFRECORD_H
#define CONSENSUSVCFRECORD_H

#include <string>
#include <vector>
#include <fstream>
#include <GenomeSequence.h> // get ref allele

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;

class ConsensusVcfRecord
{
	public:
		ConsensusVcfRecord( int position );
		~ConsensusVcfRecord();
		
		bool CheckPosition( int & location );
		
		int GetPosition();
		int GetRawVariantQuality();
		int GetEvidence();
		int GetWinCount();
		int GetSumDosage();
		
		void SetSampleSize( int ns );
		void SetAltAllele( string & alt );
		void SetInfoFields( int sp_index, string & info_str );
		void SetGLFields( int sp_index, string & gl_str );
		
		void SetRefAllele( GenomeSequence* gs, string & chr_name );
		void EstimateAFviaEM();
		void SetFilter( vector<float> & dps, vector<int> & ref_dp_cuts, vector<int> & alt_dp_cuts );
		
		void Print( string & chr, ofstream & out_vcf );
		
	private:
		void setProbViaGLs( vector< vector<float> > & prob_gls );
		void printInfoStr( ofstream & ovcf );
		void printGLStr( int sp, ofstream & ovcf );
		string getGenotypeFromGLandAF( int sp );
		int getDosageFromGenotype( string & gt );
		int getVarintQualityWithAF();
		void setPLsFromGL( int sp, int* vp); // convert GL to PL
		void estimateRawVariantQuality();
		void updateCIPOS( int pos );
		void updateVariantEnd( int current_end );
		void updateCIEND( int current_end );
		int countSampleWithVariant();
	
		int position;
		char ref_allele;
		string alt_allele;
		int variant_quality;
		string filter;
		float gt_freq[2];
		vector< int > Dosages;
		vector< int > DPs;
		vector< int > GQs;
		vector< vector<float> > GLs;
	// info fields
		int lanchor;
		int ranchor;
		int variant_end_baseline;
		int variant_end_diff_sum;
		int left_most;
		int right_most;
		int left_most_end;
		int right_most_end;
		int evidence;
		int wcount;
		float ref_af;
		float hwe;
		float fic;	
};


#endif
