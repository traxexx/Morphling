#ifndef SITES_H
#define SITES_H

#include "SamFile.h"
#include "SamFileHeader.h"
#include "SamRecord.h"
#include "AsbSite.h"

#include <string>
#include <vector>
#include <fstream>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;

// store site info. if gt == 0, do not read from this bam, in GtList mark as FALSE
struct PotSite {
	int Position;
	int MEtype;
	bool* GtList;
};

// contain all site info from vcf
class Sites
{
  public:
	Sites( string & vcf_name, string & sample_list_name, string & out_vcf_name, string & me_list_name );
	void AssemblySubtypes();

  private:
  	void initializeSiteInfo();
  // set site
  	void setSiteAndSamples();
  	void setBamFilesFromSampleList( string & sample_list_name );
  	void printAddedVcfHeaders();
  	void loadMEsequence( string & me_list_name );
  	void printSingleRecord( int pend, int mei_index, string & vline, AsbSite & cs );
 
 	string InVcfName;
  	ofstream OutVcf;
	vector<string> SampleNames;
	vector< PotSite > SiteInfo;
	vector<string> BamFileNames;
	vector<SamFile> BamFiles;
	vector<SamFileHeader> BamFileHeaders;
	vector< vector<string> > MEseqs;
	vector< vector<string> > MEnames;
};

int GetMEtypeFromAlt( string & field);
int GetTabLocation( int search_start, int noccur, string & line );

#endif