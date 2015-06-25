#ifndef SITELISTMAP_H
#define SITELISTMAP_H

#include <map>
#include <string>
#include <vector>
#include "VcfRecord.h"

using std::map;
using std::string;
using std::vector;

struct smComp // SiteListMap per-chr order
{
	bool operator()( const int & a1, const int & a2 ) const
	{
		return a1 < a2;
	}
};

typedef map<int, vector<int>, smComp > SingleSiteListMap; // position --> <all evidence-depth>
typedef map<int, int, smComp > SiteListMap; // chr ->[ center -> evidence depth ]

class SiteList
{
  public:
	SiteList( vector< vector<string> > & SampleList, string & sample_suffix );
	void Print( string & out_list_name );
	bool IsCandidateSite( VcfRecord & vcf_rec );
	
  private:
  	void addFromVcfRecord( SingleSiteListMap & singleList, VcfRecord & vcf_rec );
  	void addClusterToSiteList( SingleSiteListMap::iterator & start_anchor, SingleSiteListMap::iterator & next_of_end);
  	void addSingleListToSiteList( SingleSiteListMap & singleList );
  	void updateVariantQuality( string & dname );

// data structure  	
  	SiteListMap siteList;
// dynamics
	int MinVariantQuality;  
};

// other site list constructing related functions

// if can be added as candidate site, return TURE
bool IsCandidateSite( VcfRecord & vcf_rec );

// read all site list from files
void LoadSiteList( vector<int> & siteVec, string & site_list_name );

#endif
