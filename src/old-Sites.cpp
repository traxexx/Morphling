#include "Sites.h"
#include "Globals.h"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <math.h>

using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;

// constructor: prepare for asembly
Sites::Sites( int mtype, vector<string> & chr_list ) {
	mei_type = mtype;
	pchr_list = &chr_list;
}


void Sites::SetSiteList( vector<string> & dirs )
{
	for( int c=0; c<pchr_list->size(); c++) {
		for( vector<string>::iterator d = dirs.begin(); d!= dirs.end(); d++ ) {
			string vcf_name = (*d) + "split/Hits-" + (*pchr_list)[c] + "." + std::to_string(mei_type) + ".vcf";
			std::ifstream current_vcf;
			current_vcf.open( vcf_name.c_str() );
			if( !current_vcf.is_open() ) {
				cerr << "Warning: " << vcf_name << "does not exist. Skipped!" << endl;
				continue;
			}
			string line;
			while( getline(current_vcf, line) ) {
				if ( IsHeaderLine(line) )
					continue;
				VcfSiteRecord vrec;
				vrec.SetFromLine( line );
				if( isCandidateSite(vrec) )
					addToSiteListFromVcfSiteRecord( vrec );
			}
			current_vcf.close();
		}
	}
}

bool Sites::isCandidateSite( VcfSiteRecord & vrec )
{
	return 1;
}

void Sites::addToSiteListFromVcfSiteRecord( VcfSiteRecord & vrec )
{
	int new_key = round( vrec.GetPosition() / WIN);
	int new_key1 = new_key - 1;
	int new_key2 = new_key - 2;

	addKeyToSiteList( new_key, vrec );
	addKeyToSiteListExistOnly( new_key1, vrec );
	addKeyToSiteListExistOnly( new_key2, vrec );
}

void Sites::addKeyToSiteList( int new_key, VcfSiteRecord & vrec )
{
	string chr = vrec.GetChromosome();
	int position = vrec.GetPosition();
	int evidence = vrec.GetEvidenceDepth();
	
	map<int, SiteListElement>::iterator it = siteList[chr].find(new_key);
	if ( it != siteList[chr].end() ) { // merge and adjust
//		int dist = abs(it->second.position - position);
		int new_position = (position * evidence + it->second.position * it->second.evidence) / (evidence + it->second.evidence);
		int new_evidence = evidence + it->second.evidence;
		it->second.position = new_position;
		it->second.evidence = new_evidence;
		it->second.nsample++;
	}
	else { // new key in hash
		SiteListElement new_element;
		new_element.position = position;
		new_element.evidence = evidence;
		new_element.nsample = 1;
		siteList[chr][new_key] = new_element;
	}
}

void Sites::addKeyToSiteListExistOnly( int new_key, VcfSiteRecord & vrec )
{
	string chr = vrec.GetChromosome();
	int position = vrec.GetPosition();
	int evidence = vrec.GetEvidenceDepth();

	map<int, SiteListElement>::iterator it = siteList[chr].find(new_key);
	if ( it == siteList[chr].end() )
		return;

	int dist = abs(it->second.position - position);
	if ( dist > WIN )
		return;

// begin add	
	int new_position = (position * evidence + it->second.position * it->second.evidence) / (evidence + it->second.evidence);
	int new_evidence = evidence + it->second.evidence;
	int new_nsample = it->second.nsample;
	int adjust_key = round(new_position/WIN);
	if (adjust_key != new_key) {
		int original_key = new_key;
		map<int, SiteListElement>::iterator prev_t = it;
		while(adjust_key != original_key) {
			map<int, SiteListElement>::iterator t = siteList[chr].find(adjust_key);
			if ( t != siteList[chr].end() ) {
				new_position = (new_position*new_evidence + t->second.position*t->second.evidence) / (new_evidence+t->second.evidence);
				new_evidence = new_evidence + t->second.evidence;
				if ( t->second.nsample > new_nsample )
					new_nsample = t->second.nsample;
				original_key = adjust_key;
				adjust_key = round(new_position/WIN);
				if ( original_key == adjust_key ) {
					t->second.position = new_position;
					t->second.evidence = new_evidence;
					t->second.nsample = new_nsample + 1;
					siteList[chr].erase(original_key);
					break;
				}
				else { // continue the loop
					prev_t = t;
					siteList[chr].erase(original_key);
				}
			}
			else { // move to new key. get out of loop
				SiteListElement new_element;
				new_element.position = new_position;
				new_element.evidence = new_evidence;
				new_element.nsample = new_nsample + 1;
				siteList[chr].erase(original_key);
				siteList[chr][adjust_key] = new_element;
				break;
			}	
		}
	}
	else { // same key
		it->second.position = new_position;
		it->second.evidence = new_evidence;
		it->second.nsample++;
	}
}


void Sites::ExportPositions( map<string, map<int, int> > & esites )
{
	for( SiteList::iterator c=siteList.begin(); c!= siteList.end(); c++ ) {
		esites[c->first];
		for( map<int, SiteListElement>::iterator p=c->second.begin(); p!=c->second.end(); p++ ) {
			esites[c->first][p->first] = p->second.position;
		}
	}
}






