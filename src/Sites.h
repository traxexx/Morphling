#ifndef SITES_H
#define SITES_H

#include "MeiSeqs.h"
#include "AsbSite.h"

#include <string>
#include <vector>
#include <fstream>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;


// contain all site info from vcf
class Sites
{
  public:
	Sites( string & vcf_name, vector<string> & PreAsb, string & out_vcf_name, string & me_list_name );
	void AssemblySubtypes();

  private:
  	void initializeSiteInfo();
  	void loadPreAsb( vector<string> & PreAsb );
  	void loadSinglePreAsb( string & pre_name, int sp );
  	void keepMaxSubtype();
  	void setSinlgeSiteSubtypeAndStrand( vector< vector< subCluster > > & sclust, int sub );
  	void clearOtherSubcluster( vector< subCluster > & sc, int sp );
  	void loadPreSeq( vector<string> & PreAsb );
  // set site
  	void setSiteAndSamples();
  	void printAddedVcfHeaders();
  	void printUnAssembledRecord( string & vline, int site );
  	void printSingleRecord( string & vline, int site, AsbSite & cs );
  	void setAssemblyFilter( string & filter, int & site, AsbSite & cs );
 
 	string InVcfName;
  	ofstream OutVcf;
  	vector<int> MeiType;
  	vector<int> NpolyA;
  	vector<int> NpolyT;
	int Nsite;
	vector<string> PreNames;
	
	vector< vector<string> > MEseqs;
	vector< vector<string> > MEnames;	
	vector< vector< vector<string> > > Seqs; // site->sample->seqs
	vector<int> Subtypes;
	vector<bool> Strands;

// data record
	vector< vector< vector< subCluster > > > Clusters; // sites-> 0~4 -> subtype -> read map
};


#endif


