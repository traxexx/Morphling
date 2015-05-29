#include "SiteListMap.h"
#include "Globals.h"
#include "Utilities.h"
#include "GLs.h"
#include <math.h>
#include <fstream>
#include <sstream>

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;


SiteList::SiteList( vector< vector<string> > & SampleList, string & sample_suffix )
{
	siteList.clear();
	SingleSiteListMap singleList;
	for( vector< vector<string> >::iterator info = SampleList.begin(); info != SampleList.end(); info++ ) {
		ifstream current_vcf;
		string vcf_name = (*info)[2] + sample_suffix;
		current_vcf.open( vcf_name.c_str() );
		if ( !current_vcf.is_open() ) {
			cerr << "Warning: " << vcf_name << " does not exist. Skipped!" << endl;
			continue;
		}
		string line;
		while( getline( current_vcf, line ) ) {
			bool is_info_line = IsInfoLine( line );
			if ( is_info_line )
				continue;
			VcfRecord vcf_rec;
			vcf_rec.SetFromLine( line );
			bool is_candidate_site = IsCandidateSite( vcf_rec );
			if ( !is_candidate_site )
				continue;
			addFromVcfRecord( singleList, vcf_rec );			
		}
		current_vcf.close();
	}
	addSingleListToSiteList( singleList );
	singleList.clear();
}

// print by chr: out_vcf_name.chrname
// format (2 cols): round(voted-position) total-depth
// site map key is already numerically ordered when adding

void SiteList::Print( string & out_list_name )
{
	ofstream out_list;
	out_list.open( out_list_name.c_str() );
	CheckOutFileStatus( out_list, out_list_name.c_str() );
	for( SiteListMap::iterator loc_it = siteList.begin(); loc_it != siteList.end(); loc_it++ ) {
		int position = round(loc_it->first);
		out_list << position << "\t" << loc_it->second << endl;
	}
	out_list.close();
}


// add to single list from vcf rec
// position --> #evidence
void SiteList::addFromVcfRecord( SingleSiteListMap & singleList, VcfRecord & vcf_rec )
{
// get necessary info from rec
	int position = vcf_rec.GetPosition();
	int evidence = vcf_rec.GetEvidenceDepth();
	
	SingleSiteListMap::iterator loc_it = singleList.find( position );
	if ( loc_it == singleList.end() ) { // position not exist, add
		singleList[ position ].clear();
		singleList[ position ].push_back( evidence );
	}
	else { // position exists
		singleList[ position ].push_back( evidence );
	}
}

// use 1 as anchor:
// 		if ( cluster start <----> current anchor  > max dist)
//				add cluster to siteList (exclude current anchor)
//				cluster start ---> to find an anchor with anchor <-----> current <= max_disc, set anchor start
//						if not found, anchor start is current anchor
//				update cluster with anchor start
//		(now it must be disc <= max dist)
//		nothing do to do since new anchor is set
//		
void SiteList::addSingleListToSiteList( SingleSiteListMap & singleList )
{
	siteList.clear();
// add from single
	SingleSiteListMap::iterator start_anchor = singleList.begin();
	SingleSiteListMap::iterator current_anchor = singleList.begin();
	current_anchor++;
	for( ; current_anchor != singleList.end(); current_anchor++ ) {
		int dist = current_anchor->first - start_anchor->first;
		// update
		if ( dist > WIN ) {
		// add to siteList
			addClusterToSiteList( start_anchor, current_anchor );
		// set new start anchor
			start_anchor++;
			bool is_self_anchor = 1;
			for( SingleSiteListMap::iterator ptr = start_anchor; ptr != current_anchor; ptr++ ) {
				int dist = current_anchor->first - ptr->first;
				if ( dist <= WIN ) {
					start_anchor = ptr;
					is_self_anchor = 0;
					break;
				}
			}
			if ( is_self_anchor )
				start_anchor = current_anchor;
			dist = current_anchor->first - start_anchor->first;
		}
	}
// check last anchor
	addClusterToSiteList( start_anchor, current_anchor);
}

void SiteList::addClusterToSiteList( SingleSiteListMap::iterator & start_anchor, SingleSiteListMap::iterator & next_of_end )
{
	int total_depth = 0;
	int center = 0;
	for( SingleSiteListMap::iterator ptr = start_anchor; ptr != next_of_end; ptr++ ) {
		int sum = GetVecSum( ptr->second );
		total_depth += sum;
		center += ptr->first * sum;
	}
	center = round( center/ float(total_depth) );
	siteList[ center ] = total_depth;
}

// check multiple fields in vcf record. If pass all, return TRUE
bool IsCandidateSite( VcfRecord & vcf_rec )
{
// check filter	
	bool is_pass = vcf_rec.GetFilterPassOrNot();
	if ( !is_pass )
		return 0;

// check variant quality
	int qual = vcf_rec.GetVariantQuality();
	if ( qual < MIN_VARIANT_QUALITY )
		return 0;

// all checkpoint passed, return true
	return 1;			
}


// read all site list from files
void LoadSiteList( vector<int> & siteVec, string & site_list_name )
{
	ifstream site_list;
	site_list.open( site_list_name.c_str() );
	CheckInputFileStatus( site_list, site_list_name.c_str() );
	string line;
	while( getline( site_list, line ) ) {
		std::stringstream ss;
		ss << line;
		string pos_str;
		getline( ss, pos_str, '\t');
		int position = std::stoi( pos_str );
		siteVec.push_back( position );
	}
	site_list.close();
}


