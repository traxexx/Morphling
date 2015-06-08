#ifndef SITES_H
#define SITES_H

#include <string>

struct PosSite {
	int Position;
	bool* GtList;
}

class Sites
{
  public:
	Sites( string & vcf_name, string & sample_list_name );
	AssemblySubtypes();

  private:
  	void setSampleNamesFromVcf(string & vcf_name);
  	void reorganizeBamsFromSampleList( string & sample_list_name); // to match the one from vcf
  	void openAllBamFiles();
  	
	vector<string> SampleNames;
	vector<string> BamFileNames;
	vector<SamFile> bamFiles;
	vector< PotSite > SiteInfo; // if 1/1 or 1/0, mark as TRUE

};

#endif