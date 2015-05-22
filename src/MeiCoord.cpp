#include "MeiCoord.h"
#include "Utilities.h"

#include <fstream>
#include <sstream>

// constructor
MeiCoord::MeiCoord( string & mei_coord_list )
{
	std::ifstream list_file;
	list_file.open( mei_coord_list.c_str() );
	CheckInputFileStatus(list_file, mei_coord_list.c_str());
	
	string line_str;
	int mei_type = 0;
	string line1;
	while( std::getline( list_file, line1) ) {
		std::stringstream ss1; ss1 << line1;
		string field1;
		std::getline(ss1, field1, '\t');
		std::getline(ss1, field1, '\t');
		std::ifstream coord_bed;
		coord_bed.open( field1.c_str() );
		CheckInputFileStatus( coord_bed, field1.c_str() );
		string line;
		Coord new_coord;
		string last_chr;
		MultiCoordMap::iterator chr_it;
		while(std::getline( coord_bed, line )) {
			std::stringstream ss; ss << line;
			string chr;
			std::getline(ss, chr, '\t');
			if (chr.compare(last_chr) != 0) {
				last_chr = chr;
				RefCoord[last_chr].resize(3);
				chr_it = RefCoord.find(last_chr);
			}
			string field;
			std::getline(ss, field, '\t');
			new_coord.start = stoi(field);
			std::getline(ss, field, '\t');
			new_coord.end = stoi(field);
//std::cout << last_chr << " " << new_coord.start << " " << new_coord.end << std::endl;
			MultiCoord::iterator vvit = chr_it->second.begin();
			vvit += mei_type;
			vvit->push_back(new_coord);
		}
		coord_bed.close();
		mei_type++;		
	}
	
	list_file.close();
	No_rec_on_this_chr = 0; // set zero first
}

// destructor
MeiCoord::~MeiCoord() {}

// call it when change chr in process disc
void MeiCoord::ResetVecPtr( string & current_chr )
{
	CordPtrVec.clear();
	EndPtrVec.clear();

	if ( RefCoord.find(current_chr) == RefCoord.end() ) {
		No_rec_on_this_chr = 1;
		return;
	}
// then do regular stuff
	No_rec_on_this_chr = 0;
	CordPtrVec.reserve(3);
	EndPtrVec.reserve(3);
	MultiCoord::iterator it = RefCoord[current_chr].begin();
	for(int i=0; i<3; i++, it++) {
		CordPtrVec[i] = it->begin();
		EndPtrVec[i] = it->end();
	}
}


// set EM tag in a single sam record
// 1  2  3  A  L  S
void MeiCoord::SetEMtag( SamRecord & sam_rec )
{
// all return em = 0 if no rec
	if (No_rec_on_this_chr) {
		sam_rec.addIntTag("EM", 0);
		return;
	}
	
// regular
	if (sam_rec.getFlag() & 0x4 || sam_rec.getFlag() & 0x8) {
		sam_rec.addIntTag("EM", 0);
		return; // nothing else to do
	}
// check if inside RefCoord
	int tag = 0;
	int align_tag = 1;
	for( int i=0; i<3; i++ ) {
		if (i != 0)
			align_tag *= 2;
		CordPtr cp = CordPtrVec[i];
		if (cp == EndPtrVec[i]) 
			continue;
	// if read has not reach first mei coord....
		if ( sam_rec.get1BasedPosition() < cp->start )
			continue;
	// now see if need to adjust mei coord if coord too small
//std::cout << sam_rec.getReferenceName() << ": " << sam_rec.get1BasedPosition() << " " << cp->start << " " << cp->end << std::endl;
		while ( cp->end < sam_rec.get1BasedPosition() ) {
			cp++;
			if (cp == EndPtrVec[i]) {
				break;
			}
		}
		CordPtrVec[i] = cp;
		if (cp == EndPtrVec[i]) // for the break result
			continue;		
		if ( sam_rec.get1BasedUnclippedEnd() <= cp->end )
			tag |= align_tag;
	}
	sam_rec.addIntTag("EM", tag);
}



