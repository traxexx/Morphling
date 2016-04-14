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
	for( int m=0; m<NMEI; m++ ) {
		for( int s=0; s<msinfo.size(); s++ )
			AddDiscoverSiteFromSingleBam( m, msinfo[s] );
	}
	// clean & merge
	AdjustSiteList();
	// export
	ExportPositions( candidate_sites );
}


/*
	if a discordant read is mapped to MEI (overlap with MEI coord)
		add centor of ( anchor_end + 3*avr_ins_var )
	skip unmap & clip
*/
void Sites::AddDiscoverSiteFromSingleBam( int m, SingleSampleInfo & si )
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

	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam,bam_header, si.bam_name );
	for( int i=0; i<pchr_list->size(); i++ ) {
		string chr = (*pchr_list)[i];
//		if ( !single_chr.empty() && chr.compare(single_chr)!=0 )
//			continue;
		if ( siteList[m].find(chr) == siteList[m].end() )
			siteList[m][chr];
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
		if (siteList[m][chr].empty())
			p_reach_last = 1;
		else {
			p_reach_last = 0;
			pnearest = siteList[m][chr].begin();
		}
		SingleSite new_site; // temporary cluster. will be added to map later.
		new_site.depth = current_depth;
		bool start_new = 1; // check if need to add new_site to map and start new new_site
		SamRecord rec;
		while( bam.ReadRecord( bam_header, rec ) ) {
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
			map<string, map<int, bool> >::iterator coord_chr_ptr = meiCoord[m].find(rec.getMateReferenceName());
			bool is_in_coord = isWithinCoord( rec.get1BasedMatePosition(), coord_chr_ptr->second ); // within MEI coord
			if ( !is_in_coord )
				continue;
			if (start_new) {
				setNewCluster(new_site,rec);
				start_new = 0;
			}
			else {
				if ( rec.get1BasedPosition() > new_site.end + avr_ins_size ) { // start new coord
					addClusterToMap(new_site, siteList[m][chr]);
					setNewCluster(new_site, rec);
					start_new = 0;
				}
				else
					addToCurrentCluster(new_site, rec);
			}
		}
		// add last one
		if (!start_new)
			addClusterToMap(new_site, siteList[m][chr]);
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
	for(int m=0; m<NMEI; m++)
		exportSinglePositions(m, esites[m]);
}


// convert from
// chr-> start -> info
//to
// chr -> raw_breakp -> extend
void Sites::exportSinglePositions( int m, map< string, map<int, int> > & esites )
{
	int count = 0;
	int exclude = 0;
	for( map<string, map<int, SingleSite> >::iterator c=siteList[m].begin(); c!= siteList[m].end(); c++ ) {
		esites[c->first];
		for( map<int, SingleSite>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
			int cutoff = getMinEvidenceFromDepth( p->second.depth );
			if ( cutoff>1 && p->second.left_all * p->second.right_all == 0) {
				count++;
				exclude++;
				continue;
			}
			if (p->second.left_all + p->second.right_all < cutoff) {
				count++;
				exclude++;
				continue;
			}
			int st = p->second.start;
			int ed = p->second.end;
			if ( p->second.left_clip_only || p->second.left == 0)
				st -= WIN/2;
			if (st<0) {
				morphWarning("negative st in Sites::ExportPositions");
				st = 0;
			}
			if ( p->second.right_clip_only || p->second.right == 0) {
				if (ed<0)
					morphError("negative ed in Sites::ExportPositions");
				ed += WIN/2;
			}
			int extend = std::max(ed - p->second.breakp, p->second.breakp - st);
			extend = std::max(extend, WIN);
			esites[c->first][p->second.breakp] = extend;
//std::cout << p->second.breakp << " " << extend << " " << st << " " << ed << std::endl;
			count++;
		}
	}
	string str = "Finished Preliminary Discovery. \n    Discovered " + std::to_string(count) + " preliminary sites. Excluded " + std::to_string(exclude);
	morphMessage(str);
}


// add other type info
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


int Sites::getMinEvidenceFromDepth( int depth )
{
	int evi;
	if (depth < 60)
		evi = 1;
	else
		evi = depth / 60 + 1;
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

void Sites::setNewCluster( SingleSite & new_site, SamRecord & rec )
{
	new_site.left_clip_only = 1;
	new_site.right_clip_only = 1;
	new_site.breakp = getEstimatedBreakPoint(rec);
	if ( rec.getFlag() & 0x10 )  { // right anchor
		new_site.right = 1;
		new_site.left = 0;
		new_site.start = rec.get1BasedPosition();
		new_site.end = rec.get1BasedAlignmentEnd();
		Cigar * myCigar = rec.getCigarInfo();
		int begin_clip = myCigar->getNumBeginClips();
		if ( begin_clip < MIN_CLIP/2)
			new_site.right_clip_only = 0;
	}
	else {
		new_site.left = 1;
		new_site.right = 0;
		new_site.start = rec.get1BasedPosition();
		new_site.end = rec.get1BasedAlignmentEnd();
		Cigar * myCigar = rec.getCigarInfo();
		int end_clip = myCigar->getNumEndClips();
		if (end_clip < MIN_CLIP/2)
			new_site.left_clip_only = 0;
	}
	new_site.depth = current_depth;
}

void Sites::addToCurrentCluster( SingleSite & new_site, SamRecord & rec )
{
	// update breakpoint
	int old_evi = new_site.left + new_site.right;
	float a1 = (float)1 / float(old_evi+1);
	int ep = getEstimatedBreakPoint(rec);
	new_site.breakp = round( a1 * (float)ep + (float)new_site.breakp * (1-a1)); 

	// updae rest
	if (rec.get1BasedPosition() < new_site.start)
		new_site.start = rec.get1BasedPosition();
	else if (rec.get1BasedAlignmentEnd() > new_site.end)
		new_site.end = rec.get1BasedAlignmentEnd();
	if (rec.getFlag() & 0x10) {
		new_site.right++;
		if (new_site.right_clip_only) {
			Cigar * myCigar = rec.getCigarInfo();
			int begin_clip = myCigar->getNumBeginClips();
			if ( begin_clip < MIN_CLIP/2)
				new_site.right_clip_only = 0;			
		}
	}
	else {
		new_site.left++;
		if (new_site.left_clip_only) {
			Cigar * myCigar = rec.getCigarInfo();
			int end_clip = myCigar->getNumEndClips();
			if (end_clip < MIN_CLIP/2)
				new_site.left_clip_only = 0;			
		}
	}
}

void Sites::addClusterToMap( SingleSite & new_site, map<int, SingleSite> & smap )
{
	bool is_new_key = 0;
	if (p_reach_last)
		is_new_key = 1;
	else {
		if (pnearest->second.start > new_site.end) // far away from nearest
			is_new_key = 1;
		else {
			while(new_site.start > pnearest->second.end) {
				pnearest++;
				if (pnearest == smap.end()) {
					p_reach_last = 1;
					is_new_key = 1;
				}
			}
			if (pnearest->second.start > new_site.end)
				is_new_key = 1;		
		}
	}

	// add new key
	if (is_new_key) {
		// add to map
//		int key = (new_site.end - new_site.start) / 2 + new_site.start; // key is the average
		smap[new_site.start] = new_site;
		// set dynamics
		pnearest = smap.find(new_site.start);
		return;
	}
	
	// merge with existing key
	int old_evi = pnearest->second.left + pnearest->second.right;
	int new_evi = new_site.left + new_site.right;
	int esum = old_evi + new_evi;
	float a1 = (float)old_evi / (float)esum;
	int new_ep = round( a1 * (float)pnearest->second.breakp + (1-a1) * (float)new_site.breakp);
	pnearest->second.breakp = new_ep;
	pnearest->second.start = std::min(pnearest->second.start, new_site.start);
	pnearest->second.end = std::max(pnearest->second.end, new_site.end);
	pnearest->second.left += new_site.left;
	pnearest->second.right += new_site.right;	
	pnearest->second.depth += current_depth;
	int prev_start = pnearest->second.start;
	int prev_end = pnearest->second.end;
	// set dynamics
	pnearest++;
	if (pnearest == smap.end()) {
		p_reach_last = 1;
		return;
	}
	if (pnearest->second.start <= prev_end) {
		string str = "Added " + std::to_string(prev_start) + "-" + std::to_string(prev_end);
		str += "; but already exists " + std::to_string(pnearest->second.start) + "-" + std::to_string(pnearest->second.end);
		morphError(str, 20);
	}
}



