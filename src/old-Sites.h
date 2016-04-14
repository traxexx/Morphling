#ifndef SITES_H
#define SITES_H

#include "VcfSiteRecord.h"

using std::string;
using std::map;
using std::vector;

typedef struct {
	int position;
	int evidence;
	int nsample;
} SiteListElement;

typedef map<string, map<int, SiteListElement> > SiteList; // chr-> round(position/WIN)->element

class Sites
{
public:
	Sites( int mtype, vector<string> & chr_list);
	void SetSiteList( vector<string> & dirs );
	void ExportPositions( map<string, map<int, int> > & esites );
//	void Print();

private:
	bool isCandidateSite( VcfSiteRecord & vrec );
	void addToSiteListFromVcfSiteRecord( VcfSiteRecord & vrec );
	void addKeyToSiteList( int new_key, VcfSiteRecord & vrec );
	void addKeyToSiteListExistOnly( int new_key, VcfSiteRecord & vrec );

	SiteList siteList;

	int mei_type;
	vector<string> * pchr_list;
};

#endif
