#ifndef CONSENSUSVCF_H
#define CONSENSUSVCF_H

#include <string>
#include <fstream>
#include <map>
#include "ConsensusVcfRecord.h"

using std::ifstream;
using std::ofstream;
using std::map;

// contain info for a whole vcf. Become memory-intensive for variants with huge site list
class ConsensusVcf
{
	public:
		ConsensusVcf( int n_sample, int dist, string & chr_name );
		~ConsensusVcf();
		
		void SetSampleList( vector< vector<string> > sampleList );
		void SetAltAlleleByMeiType( int mei_index );
		void SetFasta( string & fasta_name );
		void InitializeSdataFromSiteList( string & site_list_name );
		
		void AddFromSingleVcf( int sample_index, string & vcf_name );
		void MergeData();
		void Polish(); // add Morphling-specific fields
		void Print( string & out_name );
		
	private:
		void appendToCvcf( ConsensusVcfRecord * cv, int & sp_index, string & line );
		bool NewCvcfHigherRank( ConsensusVcfRecord* old_cv, ConsensusVcfRecord* new_cv );
		void printFinalHeader( ofstream & out_vcf );
	
		vector< ConsensusVcfRecord* > sdata; // sinlge MEI vcf info
		map<int, ConsensusVcfRecord* > Data; // store all vcf info
		vector< string > SampleNames; // start from 0
		
		string chr;
		string alt_allele;
		const int _nsample;
		const int _max_dist;
		
	// something useful in final polish
		string _ref_fasta_name;
};

#endif