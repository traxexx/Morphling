#include "Sites.h"
#include "SamFile.h"
#include "Globals.h"
#include "bamUtility.h"
#include "morphError.h"
#include "seqUtility.h"
#include <fstream>
#include <map>
#include <sstream>
#include <math.h>

using std::stringstream;

// constructor: prepare for asembly
Sites::Sites( vector<string> & mei_coord_vec, vector<string> & chr_list ) {
	range_start = -1;
	range_end = -1;
//	total_depth = 0;

	pchr_list = &chr_list;
	setMeiCoord( mei_coord_vec );
}


void Sites::MakePreliminarySiteList( vector<SingleSampleInfo> & msinfo, string & range, map<int, map<string, map<int, int> > > & candidate_sites )
{
	//check range
	if (range.compare(".") != 0) {
		std::stringstream ss;
		ss << range;
		string field;
		std::getline(ss, field, '-');
		range_start = stoi(field);
		std::getline(ss, field, '-');
		range_end = stoi(field);
	}
	// add by file
	for( int s=0; s<msinfo.size(); s++ )
		AddDiscoverSiteFromSingleBam( msinfo[s] );
	// clean & merge
//	AdjustSiteList();
	// export
	ExportPositions( candidate_sites );
}


/*
	if a discordant read is mapped to MEI (overlap with MEI coord)
		add centor of ( anchor_end + 3*avr_ins_var )
	skip unmap & clip
	check 3 types at the same time
*/
void Sites::AddDiscoverSiteFromSingleBam( SingleSampleInfo & si )
{
/*for(int i=0;i<NMEI; i++) {
std::cout << "m=" << i << ": ";
for(map<string, map<int, bool> >::iterator t=meiCoord[m].begin(); t!=meiCoord[m].end(); t++)
std::cout << t->first << "->" << t->second.size() << " ";
std::cout << std::endl;
}*/

	avr_read_length = si.avr_read_length;
	avr_ins_size = si.avr_ins_size;
	min_read_length = avr_read_length / 3;
	current_depth = si.depth;
//	total_depth += current_depth;

	resetDepthAddFlag();
	
	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam,bam_header, si.bam_name );
	for( int i=0; i<pchr_list->size(); i++ ) {
		string chr = (*pchr_list)[i];
//		if ( !single_chr.empty() && chr.compare(single_chr)!=0 )
//			continue;
		if ( siteList.find(chr) == siteList.end() )
			siteList[chr].clear();
//		map<string, map<int, SingleSite> >::iterator pchr = siteList[m].find(chr);
//		map<string, map<int, bool> >::iterator coord_chr_ptr = meiCoord[m].find(chr);
//		if (coord_chr_ptr == meiCoord[m].end())
//			continue;
		bool section_status;
		if (range_start<0) { // no range
			section_status = bam.SetReadSection( chr.c_str() );
			if (!section_status) {
				string str = "Cannot set read section at chr " + chr;
				morphWarning( str );
			}
		}
		else { // set range
			section_status = bam.SetReadSection( chr.c_str(), range_start, range_end );
			if (!section_status) {
				string str = "Cannot set read section at chr " + chr + " " + std::to_string(range_start) + "-" + std::to_string(range_end); 
				morphWarning( str );
			}			
		}
		
		// DO ADDING
//		if (siteList[chr].empty())
//			p_reach_last = 1;
//		else {
//			p_reach_last = 0;
		pnearest = siteList[chr].begin();
//		}
		SingleSite new_site; // temporary cluster. will be added to map later.
		new_site.depth = current_depth;
		bool start_new = 1; // check if need to add new_site to map and start new new_site
		SamRecord rec;
		int between = 0; // count #reads after new_site.end. If end changed, add it to rcount and reset to zero
		while( bam.ReadRecord( bam_header, rec ) ) {
			if (!start_new) {
				if (rec.get1BasedPosition() >= new_site.end)
					between++;
				else
					new_site.rcount++;
			}
			if (rec.getFlag() & 0x2)
				continue;
			if ( OneEndUnmap( rec.getFlag() ) )
				continue;
			if ( IsSupplementary(rec.getFlag()) )
				continue;
			if ( rec.getReadLength() < min_read_length )
				continue;
			if ( rec.getMapQuality() < MIN_QUALITY )
				continue;
			if (chr.compare(rec.getMateReferenceName())==0 && rec.getInsertSize() < abs(avr_ins_size*2))
				continue;
			bool is_mei = 0;
			vector<bool> is_in_coord;
			is_in_coord.resize(3, 0);
			for(int m=0; m<NMEI; m++) {
				map<string, map<int, bool> >::iterator coord_chr_ptr = meiCoord[m].find(rec.getMateReferenceName());
				if (coord_chr_ptr == meiCoord[m].end())
					is_in_coord[m] = 0;
				else
					is_in_coord[m] = isWithinCoord( rec.get1BasedMatePosition(), coord_chr_ptr->second ); // within MEI coord
				if (is_in_coord[m])
					is_mei = 1;
			}
			if (!is_mei)
				continue;
			if (start_new) {
				setNewCluster( is_in_coord, new_site,rec);
				start_new = 0;
				between = 0;
			}
			else { // add to existing cluster
				if ( rec.get1BasedPosition() > new_site.end + avr_ins_size ) { // start new coord
					addClusterToMap(new_site, siteList[chr]);
					setNewCluster( is_in_coord, new_site, rec);
					start_new = 0;
					between = 0;
				}
				else {
					addToCurrentCluster( is_in_coord, new_site, rec);
					new_site.rcount += between;
					between = 0;
				}
			}
		}
		// add last one
		if (!start_new)
			addClusterToMap(new_site, siteList[chr]);
	}
	bam.Close();
}

void Sites::setMeiCoord( vector<string> & mei_coord_vec)
{
	for(int m=0; m<NMEI; m++) {
		meiCoord[m].clear();
		setSingleMeiCoord( m, mei_coord_vec[m] );
	}
}

void Sites::setSingleMeiCoord( int mtype, string & filename )
{
	std::ifstream file;
	file.open(filename.c_str());
	string line;
	string last_chr = "";
	map<string, map<int, bool> >::iterator pchr;
	while( std::getline( file, line ) ) {
		std::stringstream ss;
		ss << line;
		string chr;
		std::getline(ss, chr, '\t');
		if ( chr.compare(last_chr) != 0 ) {
			last_chr = chr;
			pchr = meiCoord[mtype].find(last_chr);
			if (pchr == meiCoord[mtype].end()) {
				meiCoord[mtype][last_chr];
				pchr = meiCoord[mtype].find(last_chr);
			}
		}
		string field;
		std::getline(ss, field, '\t');
		int st = stoi(field);
		std::getline(ss, field, '\t');
		int ed = stoi(field);
		if (ed - st > WIN) {
			int n = (ed - st) / WIN + 1;
			for(int i=0; i<n; i++) {
				if ( (i+1)*WIN > ed ) {
					if ( ed - i*WIN > WIN/2 ) {
						int key = round( float(st) / WIN + 0.49 + i);
						pchr->second[key];
					}
					break;
				}
				int key = round( float(st) / WIN + 0.49 + i);
				pchr->second[key] = 1;
			}
		}
		else {
			int key = round( float(st + ed / 2) / WIN + 0.49 );
			pchr->second[key] = 1;
		}
	}
	file.close();
}

// given mapping start position
// if anchor on left
//		expected_breakp = position + avr_ins_size / 2;
// if anchor on right
//		expected_breakp = position + read_length - avr_ins_size / 2;
// key = round( expected_breakp / WIN )
int Sites::getEstimatedBreakPoint( SamRecord & rec )
{
	int ep;
	int clen = GetMaxClipLen(rec);

	if ( !rec.getFlag() & 0x10 ) { // left anchor
		if (clen < -MIN_CLIP/2) // end clip of anchor decides position
			ep = rec.get1BasedAlignmentEnd();
		else
			ep = rec.get1BasedPosition() + avr_ins_size / 3;
	}
	else { // right anchor
		if (clen > MIN_CLIP/2)
			ep = rec.get1BasedPosition();
		else
			ep = rec.get1BasedPosition() + avr_read_length - avr_ins_size / 3;
	}
	return ep;
}


bool Sites::isWithinCoord( int position, map<int, bool> & coord )
{
	int key = round( (float)position / WIN );
	if (coord.find(key) != coord.end())
		return 1;
	else
		return 0;
}


void Sites::ExportPositions( map< int, map<string, map<int, int> > > & esites )
{
	int count = 0;
	int exclude = 0;
	for( map<string, map<int, SingleSite> >::iterator c=siteList.begin(); c!= siteList.end(); c++ ) {
		for(int m=0; m<NMEI; m++)
			esites[m][c->first].clear();
		for( map<int, SingleSite>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
			int cutoff = getMinEvidenceFromDepth( p->second.depth );
			bool left_exist = 0;
			bool right_exist = 0;
			for(int m=0; m<NMEI; m++) {
				if (!left_exist) {
					if (p->second.left[m]>0)
						left_exist = 1;
				}
				if (!right_exist) {
					if (p->second.right[m]>0)
						right_exist = 1;
				}
			}
			if ( cutoff>=3 && ( !left_exist || !right_exist )) {
				count++;
				exclude++;
				continue;
			}
			if (p->second.evidence < cutoff) {
				count++;
				exclude++;
				continue;
			}
			int st = p->second.start;
			int ed = p->second.end;
			if ( p->second.left_clip_only )
				st -= WIN/2;
			if (st<0) {
				morphWarning("negative st in Sites::ExportPositions");
				st = 0;
			}
			if ( p->second.right_clip_only ) {
				if (ed<0)
					morphError("negative ed in Sites::ExportPositions");
				ed += WIN/2;
			}
			int extend = std::max(ed - p->second.breakp, p->second.breakp - st);
			extend = std::max(extend, WIN);
			int n;
			if (extend > WIN*2) // extreme long extend, break into chuncks in WIN*2 size
				n = extend / WIN*2;
			else
				n = 1;
			for(int i=0; i<n; i++) {
				bool exported = 0;
				int breakp;
				if (n==1) {
					breakp = p->second.breakp;
				}
				else {
					breakp = st + (float)(ed - st) / (float)n * i + WIN*2;
					extend = WIN*1.6;
				}
				for(int m=0; m<NMEI; m++) {
					if (p->second.left[m] + p->second.right[m] >= cutoff) {
						esites[m][c->first][breakp] = extend;
						exported = 1;
					}
				}
				if (!exported) // if not then assume it is Alu. If not, will be rescued in assembly
					esites[0][c->first][breakp] = extend;
			}
//std::cout << p->second.breakp << " " << extend << " " << st << " " << ed << std::endl;
			count+=n;
		}
	}
	string str = "Finished Preliminary Discovery. \n    Discovered " + std::to_string(count) + " preliminary sites. Excluded " + std::to_string(exclude);
	morphMessage(str);
}


/* add other type info
void Sites::AdjustSiteList()
{
	for( int m=0; m<NMEI; m++ ) {
		for(map<string, map<int, SingleSite> >::iterator c=siteList[m].begin(); c!= siteList[m].end(); c++) {
			// find in other 2 maps
			for( map<int, SingleSite>::iterator p = c->second.begin(); p!= c->second.end(); p++ ) {
				int key = p->first;
				p->second.left_all = p->second.left;
				p->second.right_all = p->second.right;
				int msub [2];
				if (m==0) {
					msub[0] = 1;
					msub[1] = 2;
				}
				else if (m==1) {
					msub[0] = 0;
					msub[1] = 2;
				}
				else {
					msub[0] = 0;
					msub[1] = 1;
				}
				for( int i=0; i<=1; i++ ) {
					int sm = msub[i];
					map<int, SingleSite>::iterator t = siteList[sm][c->first].find(key);
					if ( t != siteList[sm][c->first].end() ) {
						p->second.left_all += t->second.left;
						p->second.right_all += t->second.right;
						if (!p->second.left_clip_only || !t->second.left_clip_only)
							p->second.left_clip_only = 0;
						if (!p->second.right_clip_only || !t->second.right_clip_only)
							p->second.right_clip_only = 0;
						float a1 = (float)(p->second.left_all + p->second.right_all) / (float)(p->second.left_all + p->second.right_all + t->second.left_all + t->second.right_all);
						p->second.breakp = p->second.breakp * a1 + t->second.breakp * (1-a1);
					}
				}
			}	
		}
	}
}
*/


int Sites::getMinEvidenceFromDepth( int depth )
{
	int evi;
	if (depth<20)
		evi = 1;
	else
		evi = depth / 20;
/*	
	if (depth < 60)
		evi = 1;
	else
		evi = depth / 60 + 1;
	*/
/*	
	if (depth <= 10)
		evi = 1;
	else {
		int times = (depth - 10) / 23;
		evi = times + 1;
	}
*/	
	return evi;
}

/* functions for updating site key */

void Sites::setNewCluster( vector<bool> & is_in_coord, SingleSite & new_site, SamRecord & rec )
{
	if (is_in_coord.size() != NMEI)
		morphError("[Sites::setNewCluster] is_in_coord size error");

	// set info
	new_site.breakp = getEstimatedBreakPoint(rec);
	new_site.rcount = 1;
	new_site.evidence = 1;
	for(int m=0; m<NMEI; m++) {
		new_site.left[m] = 0;
		new_site.right[m] = 0;
	}
	new_site.left_clip_only = 1;
	new_site.right_clip_only = 1;
	new_site.depth = current_depth;
	new_site.depth_add = 1;

	// set position & mtype
	if ( rec.getFlag() & 0x10 )  { // right anchor
		new_site.start = rec.get1BasedPosition();
		new_site.end = rec.get1BasedAlignmentEnd();
		Cigar * myCigar = rec.getCigarInfo();
		int begin_clip = myCigar->getNumBeginClips();
		if ( begin_clip < MIN_CLIP/2)
			new_site.right_clip_only = 0;
		for(int m=0; m<NMEI; m++) {
			if (is_in_coord[m])
				new_site.right[m] = 1;
		}
	}
	else {
		new_site.start = rec.get1BasedPosition();
		new_site.end = rec.get1BasedAlignmentEnd();
		Cigar * myCigar = rec.getCigarInfo();
		int end_clip = myCigar->getNumEndClips();
		if (end_clip < MIN_CLIP/2)
			new_site.left_clip_only = 0;
		for(int m=0; m<NMEI; m++) {
			if (is_in_coord[m])
				new_site.left[m] = 1;
		}
	}	
}

void Sites::addToCurrentCluster( vector<bool> & is_in_coord, SingleSite & new_site, SamRecord & rec )
{
	if (is_in_coord.size() != NMEI)
		morphError("[Sites::setNewCluster] is_in_coord size error");

	// update breakpoint
	int old_evi = new_site.evidence;
	float a1 = (float)1 / float(old_evi+1);
	int ep = getEstimatedBreakPoint(rec);
	new_site.breakp = round( a1 * (float)ep + (float)new_site.breakp * (1-a1));
	new_site.evidence++; 

	// update position
	if (rec.get1BasedPosition() < new_site.start)
		new_site.start = rec.get1BasedPosition();
	else if (rec.get1BasedAlignmentEnd() > new_site.end)
		new_site.end = rec.get1BasedAlignmentEnd();

	// update info
	if (rec.getFlag() & 0x10) {
		if (new_site.right_clip_only) {
			Cigar * myCigar = rec.getCigarInfo();
			int begin_clip = myCigar->getNumBeginClips();
			if ( begin_clip < MIN_CLIP/2)
				new_site.right_clip_only = 0;	
		}
		for(int m=0; m<NMEI; m++) {
			if (is_in_coord[m])
				new_site.right[m]++;
		}	
	}
	else {
		if (new_site.left_clip_only) {
			Cigar * myCigar = rec.getCigarInfo();
			int end_clip = myCigar->getNumEndClips();
			if (end_clip < MIN_CLIP/2)
				new_site.left_clip_only = 0;			
		}
		for( int m=0; m<NMEI; m++) {
			if (is_in_coord[m])
				new_site.left[m]++;
		}
	}
}


// add new_site to siteList map
void Sites::addClusterToMap( SingleSite & new_site, map<int, SingleSite> & smap )
{
	int interval = new_site.end - new_site.start;
	if (interval < WIN*2)
		interval = WIN*2;
	int max_lc = GetMaxIntervalReadCount( current_depth, avr_read_length, interval );
	if (new_site.rcount > max_lc)
		return;
/*std::cout << "\nns=" << new_site.start << std::endl;
if(!smap.empty())
std::cout << "pnearest=" << pnearest->first << std::endl;
for(map<int, SingleSite>::iterator t=smap.begin(); t!=smap.end(); t++)
std::cout << t->first << std::endl;*/

	bool is_new_key = 0;
	if (smap.empty())
		is_new_key = 1;
	else {
		if (pnearest->second.start > new_site.end)
			is_new_key = 1;
		else {
			while(new_site.start > pnearest->second.end) {
				pnearest++;
				if (pnearest == smap.end()) {
					pnearest--;
					is_new_key = 1;
					break;
				}
			}
			if (!is_new_key) {
				if (pnearest->second.start > new_site.end)
					is_new_key = 1;
			}
		}
	}

	// add new key
	if (is_new_key) {
		// add to map
		smap[new_site.start] = new_site;
		pnearest = smap.find(new_site.start);
	}
	
	// merge with existing key
	if (!is_new_key)
		mergeTwoKeys(pnearest, new_site);

	// check if need to merge nearby existing keys
	map<int, SingleSite>::iterator t = pnearest;
	// check before
	vector<int> key_to_del;
	int keep_key = -1;
	if (pnearest != smap.begin()) {
		t--;
		while (t->second.end > pnearest->second.start) {
			int small_key = t->second.start;
			int large_key = pnearest->second.start;
			mergeTwoKeys( t, pnearest->second );
			keep_key = small_key;
			key_to_del.push_back(large_key);
			if (t==smap.begin())
				break;
			else
				t--;
		}
	}
	// clear
	if (keep_key != -1) {
		for(int i=0; i<key_to_del.size(); i++)
			smap.erase(key_to_del[i]);
		t = smap.find(keep_key);
		if (t==smap.end()) {
			string str = "[Sites::addClusterToMap] (1) Cannot find keep key " + std::to_string(keep_key);
			morphError(str, 33);
		}
		pnearest = t;
		keep_key = -1;
		key_to_del.clear();
	}
	else
		t = pnearest;
	// check after
	t++;
	if (t != smap.end()) {
		while (t->second.end > pnearest->second.start) {
			int small_key = pnearest->second.start;
			int large_key = t->second.start;
			mergeTwoKeys( pnearest, t->second );
			keep_key = small_key;
			key_to_del.push_back(large_key);
			t++;
			if (t==smap.end())
				break;		
		}
	}
	if (keep_key != -1) {
		for(int i=0; i<key_to_del.size(); i++)
			smap.erase(key_to_del[i]);
		pnearest = smap.find(keep_key);
		if (pnearest==smap.end()) {
			string str = "[Sites::addClusterToMap] (2) Cannot find keep key " + std::to_string(keep_key);
			morphError(str, 33);
		}
	}	
}


void Sites::mergeTwoKeys( map<int, SingleSite>::iterator & pkeep, SingleSite & pmerge )
{
	int old_evi = pkeep->second.evidence;
	int new_evi = pmerge.evidence;
	int esum = old_evi + new_evi;
	float a1 = (float)old_evi / (float)esum;
	int new_ep = round( a1 * (float)pkeep->second.breakp + (1-a1) * (float)pmerge.breakp);
	pkeep->second.breakp = new_ep;
	pkeep->second.end = std::max(pkeep->second.end, pmerge.end);
	for(int m=0; m<NMEI; m++) {
		pkeep->second.left[m] += pmerge.left[m];
		pkeep->second.right[m] += pmerge.right[m];
	}
	pkeep->second.evidence += pmerge.evidence;
	if (pkeep->second.left_clip_only) {
		if (!pmerge.left_clip_only)
			pkeep->second.left_clip_only = 0;
	}
	if (pkeep->second.right_clip_only) {
		if (!pmerge.right_clip_only)
			pkeep->second.right_clip_only = 0;
	}
	if (!pkeep->second.depth_add) {
		pkeep->second.depth += current_depth;
		pkeep->second.depth_add = 1;
	}
}

// reset all depth_add flags in existing map
void Sites::resetDepthAddFlag()
{
	// site list: chr -> start -> info )
	for( map< string, map<int, SingleSite> >::iterator c=siteList.begin(); c!= siteList.end(); c++) {
		if (c->second.empty())
			continue;
		for( map<int, SingleSite>::iterator p=c->second.begin(); p!=c->second.end(); p++ )
			p->second.depth_add = 0;
	}
}


