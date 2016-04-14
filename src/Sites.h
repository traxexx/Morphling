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
	int rcount;
//	int break_point;
	int left[3];
	int right[3];
	int evidence;
	bool left_clip_only; // indicate if left anchor only covers left clip. adjusted with other mtype
	bool right_clip_only; // same
	int depth; // sum depth of samples that contain this sites
	bool depth_add; // indicate if current depth is added to depth on this sample
//	bool dp_update; // indicate if this has been found in single sample. Clear before every sample adding.
} SingleSite;

class Sites
{
public:
	Sites( vector<string> & mei_coord_vec, vector<string> & chr_list);
	void MakePreliminarySiteList( vector<SingleSampleInfo> & msinfo, string & range, map<int, map<string, map<int, int> > > & candidate_sites );

private:
	void AddDiscoverSiteFromSingleBam( SingleSampleInfo & si);
//	void AdjustSiteList();
	void ExportPositions( map<int, map<string, map<int, int> > > & esites );

	void setSingleMeiCoord( int mtype, string & filename );

	void setMeiCoord( vector<string> & mei_coord_vec );
	int getEstimatedBreakPoint( SamRecord & rec );
	bool isWithinCoord( int position, map<int, bool> & coord );

	void setNewCluster( vector<bool> & is_in_coord, SingleSite & new_site, SamRecord & rec );
	void addToCurrentCluster( vector<bool> & is_in_coord, SingleSite & new_site, SamRecord & rec );
	void addClusterToMap( SingleSite & new_site, map<int, SingleSite> & smap );
	int getMinEvidenceFromDepth( int depth );

	void mergeTwoKeys( map<int, SingleSite>::iterator & pkeep, SingleSite & pmerge );
	void resetDepthAddFlag();

	int range_start;
	int range_end;

	vector<string> * pchr_list;
	map< string, map<int, SingleSite> > siteList; // site list: chr -> start -> info
	map<int, map<string, map<int, bool> > > meiCoord;
	int avr_read_length;
	int avr_ins_size;
	int min_read_length;
	int current_depth;
//	int total_depth;

	// dynamics in adding sites
	map<int, SingleSite>::iterator pnearest;
//	bool p_reach_last;
};

#endif
