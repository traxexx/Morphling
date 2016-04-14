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
Sites::Sites( string mei_coord_file, vector<string> & chr_list ) {
	pchr_list = &chr_list;
	setMeiCoord( mei_coord_file );
	min_supporting = 1;
}

/*
	if a discordant read is mapped to MEI (overlap with MEI coord)
		add centor of ( anchor_end + 3*avr_ins_var )
	skip unmap & clip
*/
void Sites::AddDiscoverSiteFromSingleBam( SingleSampleInfo & si, string & single_chr, string & range )
{
	avr_read_length = si.avr_read_length;
	avr_ins_size = si.avr_ins_size;
	min_read_length = avr_read_length / 3;
	//check range
	int range_st = -1;
	int range_ed = -1;
	if (range.compare(".") != 0) {
		std::stringstream ss;
		ss << range;
		string field;
		std::getline(ss, field, '-');
		range_st = stoi(field);
		std::getline(ss, field, '-');
		range_ed = stoi(field);
	}

	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam,bam_header, si.bam_name );
	for( int i=0; i<pchr_list->size(); i++ ) {
		string chr = (*pchr_list)[i];
		if ( !single_chr.empty() && chr.compare(single_chr)!=0 )
			continue;
		if ( siteList.find(chr) == siteList.end() )
			siteList[chr];
		map<string, map<int, SingleSite> >::iterator pchr = siteList.find(chr);
		map<string, map<int, bool> >::iterator coord_chr_ptr = meiCoord.find(chr);
		if (coord_chr_ptr == meiCoord.end())
			continue;
		bool section_status;
		if (range_st<0) { // no range
			section_status = bam.SetReadSection( chr.c_str() );
			if (!section_status) {
				string str = "Cannot set read section at chr " + chr;
				morphWarning( str );
			}
		}
		else { // set range
			section_status = bam.SetReadSection( chr.c_str(), range_st, range_ed );
			if (!section_status) {
				string str = "Cannot set read section at chr " + chr + " " + range;
				morphWarning( str );
			}			
		}
		
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
			bool is_in_coord = isWithinCoord( rec.get1BasedMatePosition(), coord_chr_ptr->second ); // within MEI coord
			if ( is_in_coord ) {
				int ep = getEstimatedBreakPoint( rec );
				int key = round( (float)ep / WIN );
				map<int, SingleSite> ::iterator sptr = pchr->second.find(key);
				if (sptr != pchr->second.end()) { // key exist
					sptr->second.break_point = ((float)sptr->second.break_point * sptr->second.evidence + ep) / (sptr->second.evidence +1);
					sptr->second.evidence++;
					if (!sptr->dp_update)
						sptr->depth += si.depth;
					// check if need to update key
					int new_key = round( (float)sptr->second.break_point / WIN );
					if ( new_key != key )
						updateSiteKey( chr, sptr, new_key );
				}
				else { // key does not exist
					(pchr->second)[key].evidence = 1;
					(pchr->second)[key].break_point = ep;
					(pchr->second)[key].dp_update = 1;
					(pchr->second)[key].depth = si.depth;
				}
			}
		}
	}
	bam.Close();
}

void Sites::setMeiCoord( string & filename )
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
			pchr = meiCoord.find(last_chr);
			if (pchr == meiCoord.end()) {
				meiCoord[last_chr];
				pchr = meiCoord.find(last_chr);
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
		if (clen < -MIN_CLIP) // end clip of anchor decides position
			ep = rec.get1BaseAlignmentEnd();
		else
			ep = rec.get1BasedPosition() + avr_ins_size / 2;
	}
	else { // right anchor
		if (clen > MIN_CLIP)
			ep = rec.get1BasedPosition();
		else
			ep = rec.get1BasedPosition() + avr_read_length - avr_ins_size / 2;
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


// if key changes
//		move to new key
//		delete old key
void Sites::updateSiteKey( string & chr, map<int, SingleSite>::iterator & sptr, int new_key )
{
	map<int, SingleSite>::iterator prev_t = sptr;
	int original_key = sptr->first;
	int new_position = sptr->second.break_point;
	int new_evidence = sptr->second.evidence;
	int new_depth = sptr->second.depth;

	while( new_key != original_key ) {
		map<int, SingleSite>::iterator t = siteList[chr].find(new_key);
		if ( t != siteList[chr].end() ) { // be careful about int limit
			int esum = new_evidence + t->second.evidence;
			float a1 = (float)new_evidence / esum;
			float a2 = (float)t->second.evidence / esum;
			new_position = round( (float)new_position * a1 + (float)t->second.break_point * a2);
			new_evidence += t->second.evidence;
			int next_key = round( (float)new_position / WIN );
			if (next_key == new_key) {
				t->second.break_point = new_position;
				t->second.evidence = new_evidence;
				t->second.dp_update = 1;
				t->second.depth = new_depth;
				siteList[chr].erase(original_key);
				break;
			}
			else { // continue the loop
				prev_t = t;
				siteList[chr].erase(original_key);
				original_key = new_key;
				new_key = next_key;
			}
		}
		else { // move to new key, get out of loop
			siteList[chr][new_key].break_point = new_position;
			siteList[chr][new_key].evidence = new_evidence;
			siteList[chr][new_key].dp_update = 1;
			siteList[chr][new_key].depth = new_depth;
			break;
		}
	}
}

// convert to chr -> key -> brekp version
void Sites::ExportPositions( map<string, map<int, int> > & esites )
{
	for( map<string, map<int, SingleSite> >::iterator c=siteList.begin(); c!= siteList.end(); c++ ) {
		esites[c->first];
//int i=0;
		for( map<int, SingleSite>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
//if (p->second.break_point<0) {
//std::cout << i << " " << p->second.break_point << std::endl;
//}
//i++;
			int cutoff = getMinEvidenceFromDepth( p->second.depth );
			if (p->second.evidece < cutoff)
				continue;
			esites[c->first][p->first] = p->second.break_point;
		}
	}
}

int Sites::getMinEvidenceFromDepth( int depth )
{
	int evi;
	if (depth < 10)
		evi = 1;
	else {
		int times = (depth - 10) / 23;
		evi = times + 1;
	}
	return evi;
}




