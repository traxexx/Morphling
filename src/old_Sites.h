#ifndef SITES_H
#define SITES_H

#include <vector>
#include <string>
#include <map>
#include "SamRecord.h"
#include "singleSample.h"


using std::string;
using std::map;
using std::vector;

typedef struct {
	int start; // cluster start. Also the key in map
	int end; // cluster end
	int breakp; // estimated raw breakpoint
//	int break_point;
	int left;
	int right;
	int left_all;
	int right_all;
	bool left_clip_only; // indicate if left anchor only covers left clip. adjusted with other mtype
	bool right_clip_only; // same
	int depth; // sum depth of samples that contain this sites
//	bool dp_update; // indicate if this has been found in single sample. Clear before every sample adding.
} SingleSite;

class Sites
{
public:
	Sites( vector<string> & mei_coord_vec, vector<string> & chr_list);
	void MakePreliminarySiteList( vector<SingleSampleInfo> & msinfo, string & range, map<int, map<string, map<int, int> > > & candidate_sites );

private:
	void AddDiscoverSiteFromSingleBam( int m, SingleSampleInfo & si);
	void AdjustSiteList();
	void ExportPositions( map<int, map<string, map<int, int> > > & esites );

	void setSingleMeiCoord( int mtype, string & filename );
	void exportSinglePositions( int m, map<string, map<int, int> > & esites );

	void setMeiCoord( vector<string> & mei_coord_vec );
	int getEstimatedBreakPoint( SamRecord & rec );
	bool isWithinCoord( int position, map<int, bool> & coord );

	void setNewCluster( SingleSite & new_site, SamRecord & rec );
	void addToCurrentCluster( SingleSite & new_site, SamRecord & rec );
	void addClusterToMap( SingleSite & new_site, map<int, SingleSite> & smap );
	int getMinEvidenceFromDepth( int depth );

	int range_start;
	int range_end;

	vector<string> * pchr_list;
	map<int, map<string, map<int, SingleSite> > > siteList; // site list: mtype -> chr -> round(cluster_mid/WIN) -> info
	map<int, map<string, map<int, bool> > > meiCoord;
	int avr_read_length;
	int avr_ins_size;
	int min_read_length;
	int current_depth;

	// dynamics in adding sites
	map<int, SingleSite>::iterator pnearest;
	bool p_reach_last;
};

#endif
