#include "ReadMap.h"
#include "Utilities.h" // CheckFileStatus
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <algorithm> // find

ReadMap::ReadMap( std::string & sample_name, int & winlen, int & steplen, std::string & ref_chr, std::string & ref_fasta, std::string & mei_list)
{
	Supplementary_Info = 0;
	PCR_Duplicates = 0;
	QC_Fail = 0;
	Secondary_Alignment = 0;
	
	MinQuality = 10;
	MinClip = 20;
	ShortClip = 10;
	
	ReadMapCellSize_ = 50; //2 + 8 * 3 + 8*3
	SampleName = sample_name;
	win_len = winlen;
	step_len = steplen;
	
	REF_CHR = ref_chr;
	
	initializeMeiSeqRef( mei_list );
	initializeChrLenTable( ref_fasta );
}

ReadMap::~ReadMap()
{
}

void ReadMap::initializeMeiSeqRef(std::string & mei_list)
{
// initialize mei ref seq from mei_list
	std::ifstream mei_list_file; mei_list_file.open(mei_list);
	CheckFileStatus(mei_list_file);
	std::string alu_fasta;
	std::string line_fasta;
	std::string sva_fasta;
	std::string line;
	while(std::getline(mei_list_file, line)) {
		std::stringstream ss; ss << line;
		std::string field;
		std::getline(ss, field, '\t');
		std::string mei_field = field;
		if (mei_field.compare("ALU") == 0) {
			std::getline(ss, field, '\t');
			alu_fasta = field;
		}
		else if (mei_field.compare("L1") == 0) {
			std::getline(ss, field, '\t');
			line_fasta = field;		
		}
		else if (mei_field.compare("SVA") == 0) {
			std::getline(ss, field, '\t');
			sva_fasta = field;		
		}
		else {
			std::cerr << "ERROR: Unidentified MEI name: " << mei_field << std::endl; exit(2);
		}
	}
	mei_list_file.close();
	if (alu_fasta.length() == 0 || line_fasta.length() == 0 || sva_fasta.length() == 0) {
		std::cerr << "ERROR: Invalid MEI fasta. Check " << mei_list << std::endl; exit(2);
	}
	float map_ratio = 0.6; float mismatch_ratio = 0.2;
	REF_SEQ.resize(3);
	REF_SEQ[0] = new RefSeq;
	REF_SEQ[0]->ReadSeq(alu_fasta.c_str(), map_ratio, mismatch_ratio);
	REF_SEQ[0]->addPolyAtail();	
	REF_SEQ[1] = new RefSeq;
	REF_SEQ[1]->ReadSeq(line_fasta.c_str(), map_ratio, mismatch_ratio);
	REF_SEQ[1]->addPolyAtail();
	REF_SEQ[2] = new RefSeq;
	REF_SEQ[2]->ReadSeq(sva_fasta.c_str(), map_ratio, mismatch_ratio);
	REF_SEQ[2]->addPolyAtail();
}

void ReadMap::initializeChrLenTable(std::string & ref_fasta)
{
	std::string FaiName = std::string(ref_fasta) + ".fai";
	std::ifstream FaiFile; FaiFile.open(FaiName.c_str());
	std::string line, field;
	CheckFileStatus(FaiFile);
	while(std::getline(FaiFile, line)) {
		std::stringstream ss; ss << line;
		std::getline(ss, field, '\t');
		std::string chr = field;
		if (field.length() > 3) {
			std::string subfield = field.substr(0,3);
			if (subfield.compare("chr") == 0)
				chr = field.substr(3);
			else continue;
		}	
		if (!std::all_of(chr.begin(), chr.end(), ::isdigit)) {
			if (chr[0] != 'X' && chr[0] != 'Y')
				continue;
		}
		std::getline(ss, field, '\t');
		ChrLenTable[chr] = (stoi(field) - win_len) / step_len + 2;
	}
	FaiFile.close();
}


// do this in set map from bam since i don't want to do TE lift-over for ctrl bam :(
void ReadMap::initializeMeiCoordRef( std::string & mei_coord_list)
{
// struct:	std::map<std::string, std::vector< std::vector<Coord> > > MeiCoordRef;
  // get chr name
	for(std::map<std::string, int>::iterator ref_it = ChrLenTable.begin(); ref_it != ChrLenTable.end(); ref_it++)
		MeiCoordRef[ ref_it->first ].resize(3);
		
  // read mei coord list
	std::ifstream list_file; list_file.open( mei_coord_list.c_str() );
	CheckFileStatus( list_file );
	std::string line_str;
	int mei_type = 0;
	std::string line1;
	while( std::getline( list_file, line1) ) {
		std::stringstream ss1; ss1 << line1;
		std::string field1;
		std::getline(ss1, field1, '\t');
		std::getline(ss1, field1, '\t');
		std::ifstream coord_bed; coord_bed.open( field1.c_str() );
		CheckFileStatus( coord_bed );
		std::string line;
		Coord new_coord;
		std::string last_chr;
		std::map<std::string, std::vector< std::vector<Coord> > >::iterator chr_it;
		bool skip_this_chr;
		while(std::getline( coord_bed, line )) {
			std::stringstream ss; ss << line;
			std::string chr;
			std::getline(ss, chr, '\t');
			if (chr.compare(last_chr) != 0) {
				last_chr = chr;
				chr_it = MeiCoordRef.find(last_chr);
				skip_this_chr = chr_it == MeiCoordRef.end() ? 1 : 0;
			}
			if (skip_this_chr)  continue;
			std::string field;
			std::getline(ss, field, '\t');
			new_coord.start = stoi(field);
			std::getline(ss, field, '\t');
			new_coord.end = stoi(field);
//std::cout << last_chr << " " << new_coord.start << " " << new_coord.end << std::endl;
			std::vector< std::vector<Coord> >::iterator vvit = chr_it->second.begin();
			vvit += mei_type;
			vvit->push_back(new_coord);
		}
		coord_bed.close();
		mei_type++;
	}
	list_file.close();
}

void ReadMap::SetControls( std::string & ctrl_bam, std::string & work_dir )
{
	std::string disc_name = std::string(work_dir) + "/ctrl.disc.sam";
	if ( !ExistDoneFile( work_dir, "proper" ) ) {
		processProperReads(ctrl_bam, work_dir, disc_name, REF_CHR);
		GenerateDoneFile( work_dir, "proper" );
	}
	
	if ( !ExistDoneFile( work_dir, "disc" ) ) {
	// sort
		std::string sort_cmd = std::string("samtools sort -n ") + std::string(disc_name) + " " +  std::string(work_dir) + "/disc.nsort";
		int sort_disc_status = system(sort_cmd.c_str());
		CheckCmdStatus(sort_cmd, sort_disc_status);
		// disc
		std::string disc_bam = work_dir + "/ctrl.disc.nsort.bam";
		processDiscReads(disc_bam, work_dir, REF_CHR);
		GenerateDoneFile( work_dir, "disc" );
	}
}

void ReadMap::SetMapFromBam 
	(std::string & bam, std::string & ref_fasta, std::string & work_dir, std::string & mei_list, std::string & mei_coord_list)
{

	initializeMeiCoordRef( mei_coord_list );

	std::string disc_name = std::string(work_dir) + "/disc.sam";
	std::string empty_focus_chr = std::string("");
	if ( !ExistDoneFile( work_dir, "proper" ) ) {
		processProperReads(bam, work_dir, disc_name, empty_focus_chr);
		GenerateDoneFile( work_dir, "proper" );
	}
	
	if ( !ExistDoneFile( work_dir, "disc" ) ) {
	// sort disc (system)
		std::string sort_cmd = std::string("samtools sort ") + disc_name + " " + work_dir + "/disc.nsort";
		int sort_disc_status = system(sort_cmd.c_str());
		CheckCmdStatus(sort_cmd, sort_disc_status);
	// do disc
		processDiscReads(disc_name, work_dir, empty_focus_chr);
		GenerateDoneFile( work_dir, "disc" );
	}
}

void ReadMap::processProperReads(std::string & bam, std::string & work_dir, std::string & disc_name, std::string & focus_chr)
{
// if focus_chr length > 0, skip chrs != focus_chr
	bool single_chr = focus_chr.length() < 1 ? 0 : 1;

	SamFile sam_file, DiscSam;
	sam_file.OpenForRead(bam.c_str());
	if ( !sam_file.IsOpen() ) {
		std::cerr << "ERROR: Unable to open input bam file: " << bam << std::endl; exit(1);
	}
	
	SamFileHeader sam_header;
	sam_file.ReadHeader(sam_header);
	SamRecord sam_rec;
	
	DiscSam.OpenForWrite(disc_name.c_str());
	
	if ( !DiscSam.IsOpen() ) {
		std::cerr << "ERROR: Unable to open output disc bam: " << disc_name << std::endl; exit(1);
	}
	
	DiscSam.WriteHeader(sam_header);
	
	properDeck.resize(win_len);
	
	int line_count = 0; int line_times = 0;
	int MinDiscIns = win_len * 3;
	bool new_deck = 1; // empty deck when go into a new chr
	bool skipThisChr = 1; // skip abnormal chr, like hs37d5
	std::deque < std::vector<ProperDeckCell> >::iterator it_add_deck = properDeck.begin();
	std::vector<ProperMapCell> properRawMap;
	LastChr = std::string("");
	std::vector<CordPtr> CordPtrVec; // for set MC tag of disc
	std::vector<CordPtr> EndPtrVec; // for end purpose
	while (sam_file.ReadRecord(sam_header, sam_rec)) {
		line_count++;
/*if (line_count > 3722000)
	std::cout << line_count << std::endl;
else continue;
*/
		if (line_count > 10000) {
			line_times++; line_count = 0;
			time_t raw_time;
			time(&raw_time);
			std::cout << "Processed " << line_times << "0 k lines at: " << ctime(&raw_time) << std::endl;
		}

	// sanity check	
		int flag = sam_rec.getFlag();		
		if (flag & 0x800) {
			Supplementary_Info++; continue; // skip supplementary info
		}
		if (flag & 0x400) {
			PCR_Duplicates++; continue; // skip PCR duplicates;
		}
		if (flag & 0x200) {
			QC_Fail++; continue; // skip reads with QC fail;
		}
		if (flag & 0x100) {
			Secondary_Alignment++; continue; // skip secondary alignment
		}
		
// when change chr, rest properDeck		
		std::string current_chr = sam_rec.getReferenceName();
		if (current_chr.compare(LastChr) != 0) {
			if (!skipThisChr) { // check if we need to print the previous chr
				printProperRawMap(properRawMap, work_dir, current_chr);
			}
			properRawMap.clear();
			if (single_chr) {
				skipThisChr = current_chr.compare(focus_chr) == 0 ? 1 : 0;
			}
			else {
				skipThisChr = ChrLenTable.find(current_chr) != ChrLenTable.end() ? 0 : 1;
			}
			new_deck = 1;
			properDeck.clear(); properDeck.resize(win_len);
			it_add_deck = properDeck.begin();
			LastChr = current_chr;
			if (!skipThisChr) {
				properRawMap.resize(ChrLenTable[current_chr]);
			}
// set CordPtrVec for adding disc
			CordPtrVec.clear();
			EndPtrVec.clear();
			if ( MeiCoordRef.find(LastChr) != MeiCoordRef.end() ) {
				CordPtrVec.reserve(3);
				EndPtrVec.reserve(3);
				for(int i=0; i<3; i++) {
					std::vector< std::vector<Coord> >::iterator it = MeiCoordRef[LastChr].begin();
					CordPtrVec[i] = it->begin();
					EndPtrVec[i] = it->end();
					it++;
				}
			}
		}
		
// output disc reads
		if (skipThisChr) { // read on irregular-chr
			if (flag & 0x2) continue; // skip proper
			if (flag & 0x8) continue; // mate not mapped
			std::string mate_chr = sam_rec.getMateReferenceName();
			if (mate_chr.compare("=") == 0 || mate_chr.compare(current_chr) == 0)  continue; // on same chr, skip
			if ( ChrLenTable.find( mate_chr ) != ChrLenTable.end()) { // mate is regular, change this quality to zero
				bool success_change = sam_rec.setMapQuality(0);
				if (!success_change) {
					std::cerr << "ERROR: Unable to change map quality: " << sam_rec.getReadName() << std::endl; exit(1);
				}
				simplifySamRec(sam_rec);
				setEMtag(sam_rec, CordPtrVec, EndPtrVec);
				DiscSam.WriteRecord(sam_header, sam_rec);
			}
			continue;
		}
		else { // this chr is available
			if (!(flag & 0x2)) {
				if (flag & 0x4) {
					if (flag & 0x8) continue; // skip orphan
				
				}
				else if (!(flag & 0x8)) {
					std::string mate_chr = sam_rec.getMateReferenceName();
					if (mate_chr.compare("=") == 0 || mate_chr.compare(current_chr) == 0) {
						int ins_size = sam_rec.getInsertSize();
						if (ins_size < MinDiscIns)
							continue;
					}
				}
				simplifySamRec(sam_rec);
				setEMtag(sam_rec, CordPtrVec, EndPtrVec);
				DiscSam.WriteRecord(sam_header, sam_rec);
				continue;
			}
		}

	
	// now do proper
		int qual = sam_rec.getMapQuality();
		if (qual < MinQuality) continue;
	
		int ins_size = sam_rec.getInsertSize();
		if (new_deck) {
			new_deck = 0;
			deck_start = ins_size > 0 ? sam_rec.get1BasedPosition() : 0;
			properDeck.clear();
			properDeck.resize(win_len);
		}
		
		if (ins_size > 0) { // left in pair
			addToDeck(sam_rec);
		}
		else { // right in pair
			addPairToProperRawMap(properRawMap, sam_rec);
		}	
	}
// the last chr
	if (!skipThisChr) { // check if we need to print the previous chr
		printProperRawMap(properRawMap, work_dir, LastChr);
	}
	properRawMap.clear();
}

/*
	<--front		back-->
	0  1  2  ...	3  4 5
*/ 
void ReadMap::addToDeck (SamRecord & sam_rec)
{
	std::deque < std::vector<ProperDeckCell> >::iterator it_deck; // use for adding
	int current_loc = sam_rec.get1BasedPosition();
	int dist = current_loc - deck_start;
  // locate it_deck
  	if (dist >= win_len * 2 - 1) { // start new deck
  		deck_start = current_loc;
  		properDeck.clear();
  		properDeck.resize(win_len);
  		it_deck = properDeck.begin();
  	}
	else if ( dist >= win_len ) { // do pop & push
		int move_len = dist - win_len + 1;
		deck_start += move_len;
		for(int i=0; i<move_len; i++) {
			properDeck.pop_front();
		}
		properDeck.resize(win_len); // add empty element to the back
		it_deck = properDeck.end(); it_deck--;
	}
	else { // deck not full
		it_deck = properDeck.begin(); it_deck += dist;
	}
  // add to it_deck
  	int old_size = it_deck->size();
	it_deck->resize(old_size + 1);
	std::vector<ProperDeckCell>::iterator it_vec = it_deck->end(); it_vec--;
	it_vec->MatePosition = sam_rec.get1BasedMatePosition();
	it_vec->ReadID = sam_rec.getReadName();
	it_vec->ClipType = getClipType(sam_rec, REF_SEQ);
	it_vec->ClipLength = it_vec->ClipType & 2 ? abs(getMaxClipLength(sam_rec)) : 0;
}


void ReadMap::addPairToProperRawMap (std::vector<ProperMapCell> & properMap, SamRecord & sam_rec)
{
// find its mate
	std::deque < std::vector<ProperDeckCell> >::iterator it_deck = properDeck.begin();
	std::vector<ProperDeckCell>::iterator it_mate;
	int pos1 = sam_rec.get1BasedMatePosition();
	int deck_dist = pos1 - deck_start;
// some additional sanity check
	if (deck_dist < 0 || deck_dist >= win_len) {
		std::cerr << "Mate not found in deck (possibly low-qual): " << sam_rec.getReadName() << std::endl;
		return;
	}

// then search for mate
	it_deck += deck_dist;
	if (it_deck->size() == 1) { // if only one then that's possibly what we want (but still, a dropped zero-qual + non-dropped-non-mate-full-qual....)
		it_mate = it_deck->begin();
		if ( it_mate->MatePosition != sam_rec.get1BasedPosition() ) // not the same
			it_mate = it_deck->end();
	}
	else if (it_deck->size() > 1){ // if not search for mate pos match; if same, search ID;
		it_mate = it_deck->end();
		for( std::vector<ProperDeckCell>::iterator it_vec = it_deck->begin(); it_vec != it_deck->end(); it_vec++) {
			if (it_vec->MatePosition == sam_rec.get1BasedPosition()) {
				if (it_vec->ReadID.compare(sam_rec.getReadName()) == 0)
					it_mate = it_vec;
			}
		}
	}	
	
// decide type
	char pairType = 0;
	bool meiReadDirection = 0; // 8 Alu 16 LINE 32 SVA
	if (it_mate == it_deck->end()) {
		std::cerr << "Mate in Deck not found: " << sam_rec.getReadName() << std::endl; return;
	}

	if (it_mate->ClipType & 1) { // clip in read-1
		pairType |= 1;
		if (!(it_mate->ClipType & 2)) {// only short clip in read-1: get type of 2
			char type2 = getClipType(sam_rec, REF_SEQ);
			if (type2 & 1) { // also clip in read-2
				if (type2 & 2) { // 2 is not short clip: USE 2
					pairType = type2;
					meiReadDirection = sam_rec.getFlag() & 0x6 ? 0 : 1;
				}
			}
			else { // no clip in 2: USE 1
				pairType = it_mate->ClipType;
				meiReadDirection = sam_rec.getFlag() & 0x6 ? 1 : 0;
			}
		}
		else { // real clips in read-1
			int maxClipLen = abs(getMaxClipLength(sam_rec));
			if (maxClipLen > it_mate->ClipLength) { // USE 2
				pairType = getClipType(sam_rec, REF_SEQ);
				meiReadDirection = sam_rec.getFlag() & 0x6 ? 0 : 1;
			}
			else { // USE 1
				pairType = it_mate->ClipType;
				meiReadDirection = sam_rec.getFlag() & 0x6 ? 1 : 0;
			}
		}
	}
	else { // no clip in read-1
		char type2 = getClipType(sam_rec, REF_SEQ);
		if (type2 > 0) { // something in 2, USE 2
			pairType = type2;
			meiReadDirection = sam_rec.getFlag() & 0x6 ? 0 : 1;
		}
//		else 0, whatever
	}
	
/* add to Proper Raw Map: pairType meiReadDirection */
	int right_most = sam_rec.get1BasedAlignmentEnd();
	int first_move = ((right_most - win_len) / step_len + 1);
	std::vector<ProperMapCell>::iterator it_add = properMap.begin(); it_add += first_move;
	int off_set = pos1 / step_len - first_move;
	for(int i=0; i<= off_set; i++, it_add++) {
		if (pairType & 1) { // clip is here
			if (pairType & 2) {			
			// alu
				if (it_add->Alu.size() == 0) {
					it_add->Alu.resize(8,0);
				}
				int index = (pairType >> 2) & 3;
				if (meiReadDirection)  index += 4;
				it_add->Alu[index]++;
			 // L1	
				if (it_add->L1.size() == 0) {
					it_add->L1.resize(8,0);
				}
				bool mei = pairType & 16;
				bool isBegin = pairType & 4;
				index = ((mei << 1) | isBegin) & 3;
				if (meiReadDirection)  index += 4;
				it_add->L1[index]++;
			// SVA	
				if (it_add->SVA.size() == 0) {
					it_add->SVA.resize(8,0);
				}
				mei = pairType & 32;
				isBegin = pairType & 4;
				index = ((mei << 1) | isBegin) & 3;
				if (meiReadDirection)  index += 4;
				it_add->SVA[index]++;				
			}
			else { // short clip
				if (it_add->regular.size() == 0) {
					it_add->regular.resize(2,0);
				}
				it_add->regular[1]++;			
			}
		}
		else { // no clip
			if (it_add->regular.size() == 0) {
				it_add->regular.resize(2,0);
			}
			it_add->regular[0]++;
		}
	}
}

//
void ReadMap::printProperRawMap(std::vector<ProperMapCell> & properMap, std::string & work_dir, std::string & current_chr)
{
	std::string new_work_dir = work_dir + "/";
	std::string out_name = new_work_dir + "chr" + current_chr + ".proper.regular";
	std::ofstream regular_file; regular_file.open(out_name.c_str());
	CheckOutFileStatus(regular_file);
	
	out_name = new_work_dir + "chr" + current_chr + ".proper.Alu";
	std::ofstream Alu_file; Alu_file.open(out_name.c_str());
	CheckOutFileStatus(Alu_file);
	
	out_name = new_work_dir + "chr" + current_chr + ".proper.L1";
	std::ofstream L1_file; L1_file.open(out_name.c_str());
	CheckOutFileStatus(L1_file);
	
	out_name = new_work_dir + "chr" + current_chr + ".proper.SVA";
	std::ofstream SVA_file; SVA_file.open(out_name.c_str());
	CheckOutFileStatus(SVA_file);
	for(std::vector<ProperMapCell>::iterator it_print = properMap.begin(); it_print != properMap.end(); it_print++) {
		if (it_print->regular.size() > 0) {
			regular_file << it_print->regular[0] << " " << it_print->regular[1];
		}
		regular_file << std::endl;
		
		if (it_print->Alu.size() > 0) {
			int sum = it_print->Alu[2] + it_print->Alu[3] + it_print->Alu[6] + it_print->Alu[7];
			if (sum > 0) {
				Alu_file << "1";
			}
			else {
				Alu_file << "0";
			}
			for(std::vector<int>::iterator mit = it_print->Alu.begin(); mit != it_print->Alu.end(); mit++) {
				Alu_file << " " << *mit;
			}
		}
		Alu_file << std::endl;
		
		if (it_print->L1.size() > 0) {
			int sum = it_print->L1[2] + it_print->L1[3] + it_print->L1[6] + it_print->L1[7];
			if (sum > 0) {
				L1_file << "1";
			}
			else {
				L1_file << "0";
			}
			for(std::vector<int>::iterator mit = it_print->L1.begin(); mit != it_print->L1.end(); mit++) {
				L1_file << " " << *mit;
			}
		}
		L1_file << std::endl;
		
		if (it_print->SVA.size() > 0) {
			int sum = it_print->SVA[2] + it_print->SVA[3] + it_print->SVA[6] + it_print->SVA[7];
			if (sum > 0) {
				SVA_file << "1";
			}
			else {
				SVA_file << "0";
			}
			for(std::vector<int>::iterator mit = it_print->SVA.begin(); mit != it_print->SVA.end(); mit++) {
				SVA_file << " " << *mit;
			}
		}
		SVA_file << std::endl;
		
	}
	regular_file.close(); Alu_file.close(); L1_file.close(); SVA_file.close();
}

/* Let's do this to avoid binary search of MeiRefCoord:
	When doing proper, set coord tag "EM" to record MeiRefCoord. 1 Alu 2 L1 4 SVA
	Also do some modification to ProperMap so do not output QO tags (to reduce I/O)
*/

void ReadMap::processDiscReads(std::string & disc_bam, std::string & work_dir, std::string & focus_chr )
{
	SamFile sam_file;
	sam_file.OpenForRead(disc_bam.c_str());
	
	if ( !sam_file.IsOpen() ) {
		std::cerr << "ERROR: Unable to open input disc bam: " << disc_bam << std::endl; exit(1);
	}
	
	SamFileHeader sam_header;
	sam_file.ReadHeader(sam_header);
	SamRecord sam_rec;
	
	int line_count = 0; int line_times = 0;
// initialize disc map
	if ( focus_chr.length() >= 1 ) {
		std::map<std::string, int>::iterator mit = ChrLenTable.find(focus_chr);
		if ( mit == ChrLenTable.end() ) {
			std::cerr << "ERROR: processDiscReads: Can't find focus_chr " << focus_chr << " in ChrLenTable!" << std::endl;
			exit(1);
		}
		DiscRawMap[mit->first].resize(mit->second);
	}
	else { // all chrs
		for(std::map<std::string, int>::iterator mit = ChrLenTable.begin(); mit != ChrLenTable.end(); mit++)
			DiscRawMap[mit->first].resize(mit->second);
	}
	
	std::cout << std::endl;
	std::cout << "Processing Disc Bam..." << std::endl;
	DiscReadPair * rp; rp = new DiscReadPair;
	rp->is_empty = 1;
	int SanityCheckFailPairs = 0;
	while (sam_file.ReadRecord(sam_header, sam_rec)) {
		line_count++;
		if (line_count > 1000000) {
			line_times++; line_count = 0;
			time_t raw_time;
			time(&raw_time);
			std::cout << "    Processed " << line_times << " million lines at: " << ctime(&raw_time) << std::endl;
		}
		
		if (rp->is_empty) { // record info of 1st read
			setDiscReadPair(rp, sam_rec);
		}
		else { // check if second -> add disc pair info
			if (rp->id.compare(sam_rec.getReadName()) == 0) { // is second: sanity check -> re-map -> add to map
				if ( sanityCheckDiscReadPair(rp, sam_rec) ) {
					addDiscPairToRawMap(rp, sam_rec);
				}
				else {
					SanityCheckFailPairs++;
				}
				clearDiscReadPair(rp);
			}
			else { // not match...get a new one
				setDiscReadPair(rp, sam_rec);
			}
		}
		
	}
	sam_file.Close();
	std::cout << "Discordant reads processed. #Sanity check failed pairs = " << SanityCheckFailPairs << "." << std::endl;
// print out
	if (focus_chr.length() < 1) {
		for( std::map<std::string, std::vector<DiscMapCell> >::iterator disc_it = DiscRawMap.begin(); disc_it != DiscRawMap.end(); disc_it++ ) {
			if ( isEmptyDiscRawMapUnit( disc_it->second ) )  continue; // skip those with no disc stats. This is more convenient in test mode
			printDiscRawMap( work_dir, disc_it->first);
		}
	}
	else {
		printDiscRawMap(work_dir, focus_chr);
	}
}


// add from read-pair
void ReadMap::addDiscPairToRawMap (DiscReadPair * rp, SamRecord & sam_rec)
{	
// remap
	if (rp->first_map && rp->second_map) {// both disc
		if (rp->first_high_qual) {
			addDiscSingleToRawMap(rp, sam_rec, 1, 1);
		}
		if (rp->second_high_qual) {
			addDiscSingleToRawMap(rp, sam_rec, 1, 0);
		}	
	}
	else { // mate unmap
		if (!rp->first_map) { // 1 as anchor
			addDiscSingleToRawMap(rp, sam_rec, 0, 1);
		}	
		else if (!rp->second_map) {
			addDiscSingleToRawMap(rp, sam_rec, 0, 0);
		}	
	}
}


void ReadMap::addDiscSingleToRawMap(DiscReadPair * rp, SamRecord & sam_rec, bool mate_mapped, bool first_anchor)
{
	int base = mate_mapped ? 0 : 4;
	// 1 Alu 2 L1 4 SVA
	int disc_type;
	if (first_anchor) {
		rescueSingleMCinReadPair( rp->mc1, rp->seq1, REF_SEQ );
		disc_type = rp->mc1;
	}
	else {
		rescueSingleMCinReadPair( rp->mc2, rp->seq2, REF_SEQ );
		disc_type = rp->mc2;
	}
	
	bool direction; // L 1 R 0
	int flag = sam_rec.getFlag();
	if (first_anchor)
		direction = flag & 0x20 ? 0 : 1;
	else
		direction = flag & 0x10 ? 0: 1;
		
	int left_most;
	int right_most;
	std::string last_chr;
	if (first_anchor) {
		left_most = sam_rec.get1BasedMatePosition();
		right_most = rp->first_alignment_end;
		last_chr = sam_rec.getMateReferenceName();
	}
	else {
		left_most = sam_rec.get1BasedPosition();
		right_most = sam_rec.get1BasedAlignmentEnd();
		last_chr = sam_rec.getReferenceName();
	}	
	std::vector<DiscMapCell>::iterator it_add = DiscRawMap[last_chr].begin();
	int first_move = ((right_most - win_len) / step_len + 1);
	it_add += first_move;
	int off_set = left_most / step_len - first_move;
	for(int i=0; i<=off_set; i++, it_add++) {
		if (it_add->Alu.size() == 0) {
			it_add->Alu.resize(8,0);
		}
		int index = disc_type & 1 ? 1 : 0;
		if (direction)  index += 2;
		index += base;
		it_add->Alu[index]++;
		
		if (it_add->L1.size() == 0) {
			it_add->L1.resize(8,0);
		}
		index = disc_type & 2 ? 1 : 0;
		if (direction)  index += 2;
		index += base;
		it_add->L1[index]++;
		
		if (it_add->SVA.size() == 0) {
			it_add->SVA.resize(8,0);
		}
		index = disc_type & 4 ? 1 : 0;
		if (direction)  index += 2;
		index += base;
		it_add->SVA[index]++;
	}	
}


void ReadMap::printDiscRawMap( std::string & work_dir, std::string current_chr )
{
	std::string new_work_dir = work_dir + "/";
	std::string out_name = new_work_dir + "chr" + current_chr + ".disc.Alu";
	std::ofstream Alu_file; Alu_file.open(out_name.c_str());
	CheckOutFileStatus(Alu_file);
	out_name = new_work_dir + "chr" + current_chr + ".disc.L1";
	std::ofstream L1_file; L1_file.open(out_name.c_str());
	CheckOutFileStatus(L1_file);
	out_name = new_work_dir + "chr" + current_chr + ".disc.SVA";
	std::ofstream SVA_file; SVA_file.open(out_name.c_str());
	CheckOutFileStatus(SVA_file);
	for(std::vector<DiscMapCell>::iterator it_print = DiscRawMap[current_chr].begin(); it_print != DiscRawMap[current_chr].end(); it_print++) {
		if (it_print->Alu.size() > 0) {
			int sum = it_print->Alu[1] + it_print->Alu[3] + it_print->Alu[5] + it_print->Alu[7];
			if (sum > 0) {
				Alu_file << "1";
			}
			else {
				Alu_file << "0";
			}
			for(std::vector<int>::iterator mit = it_print->Alu.begin(); mit != it_print->Alu.end(); mit++) {
				Alu_file << " " << *mit;
			}
		}
		Alu_file << std::endl;
		
		if (it_print->L1.size() > 0) {
			int sum = it_print->L1[1] + it_print->L1[3] + it_print->L1[5] + it_print->L1[7];
			if (sum > 0) {
				L1_file << "1";
			}
			else {
				L1_file << "0";
			}
			for(std::vector<int>::iterator mit = it_print->L1.begin(); mit != it_print->L1.end(); mit++) {
				L1_file << " " << *mit;
			}
		}
		L1_file << std::endl;
		
		if (it_print->SVA.size() > 0) {
			int sum = it_print->SVA[1] + it_print->SVA[3] + it_print->SVA[5] + it_print->SVA[7];
			if (sum > 0) {
				SVA_file << "1";
			}
			else {
				SVA_file << "0";
			}
			for(std::vector<int>::iterator mit = it_print->SVA.begin(); mit != it_print->SVA.end(); mit++) {
				SVA_file << " " << *mit;
			}
		}
		SVA_file << std::endl;
	}
	Alu_file.close(); L1_file.close(); SVA_file.close();
}


/****************************************************************************************************************/
/**************************** Do Likelihood ***************************/

void ReadMap::setGenomicIndexToControl() // set merge & set link
{
	if ( ChrLenTable.find( REF_CHR ) == ChrLenTable.end() ) {
		std::cerr << "ERROR: Can't find REF_CHR: " << REF_CHR << " in ChrLenTable!" << std::endl;
		exit(1);
	}

	std::vector<int> duplicates;
	setDupVecForReadMapCellVec(duplicates, RawCtrlMap);

  // set link & merge cell
	MergeCtrlMap.resize(duplicates.size() - 1);
	CtrlLocationLink.resize(ChrLenTable[REF_CHR], MergeCtrlMap.end()); // so all invalid iterators in clink points to merge map end
	std::vector<StatCell>::iterator merge_it = MergeCtrlMap.begin();
	std::vector<ReadMapCell>::iterator raw_it = RawCtrlMap.begin();
	raw_it += duplicates[0];
	for(std::vector<int>::iterator dp_it = ++duplicates.begin(); dp_it != duplicates.end(); dp_it++) {
		merge_it->count = *dp_it;
		merge_it->log_table = raw_it->stats;
		for(int ct = 0; ct < *dp_it; ct++) {
			CtrlLocationLink[raw_it->win_index] = merge_it;
		}
		raw_it += *dp_it;
		merge_it++;
	}

// clear raw ctrl
	RawCtrlMap.clear();	
}

/* het ctrl index:
	mei-type	position	ref-type	lift-over-length
0 Alu 1 L1 2 SVA

  set lift over before set ctrl stat:
	lift-over all CtrlLocationLink
	When SetCtrlStat after that, set stats based on hs37d5 (using lift-length) instead of setting with lift-over hs37d5
	method: read lift-over mei index into sorted-vector --> (from back) lift-over every CtrlLocationLink */
void ReadMap::setCtrlLiftOver(std::string & het_index_name)
{
	std::vector< std::pair<int, int> > LiftOver; // for doing lift over: index lift-LENGTH(need to adjust)
	std::ifstream het_index_file; het_index_file.open(het_index_name.c_str());
	CheckFileStatus(het_index_file);
	std::string line;
	while(std::getline(het_index_file, line)) {
		std::stringstream ss; ss << line;
		std::string field;
    	std::getline(ss, field, '\t');
    	std::string str_st;
    	std::getline(ss, str_st, '\t');
    	std::getline(ss, field, '\t');
    	if (field.compare("0") != 0)  continue;
    	int st = stoi(str_st);
    	std::getline(ss, field, '\t');
    	int lift_length = stoi(field);
    	int index = (st - win_len ) / step_len;
    	LiftOver.push_back( std::make_pair(index, lift_length) );
	}
	het_index_file.close();
	
  // sort pair & adjust length to length-index
	std::sort(LiftOver.begin(), LiftOver.end(), sortLiftOver);
	for( std::vector< std::pair<int, int> >::iterator lit = LiftOver.begin(); lit != LiftOver.end(); lit++) {
		lit->second = round(float(lit->second) / step_len);
	}
	
// do lift-over on CtrlLocationLink
	std::vector< std::pair<int, int> >::iterator lit = LiftOver.begin();
	int i = CtrlLocationLink.size();
	for( std::vector< std::vector<StatCell>::iterator >::iterator cit = (--CtrlLocationLink.end()); cit != CtrlLocationLink.begin(); cit--) {
		if (i <= lit->first)	lit++;
		std::vector< std::vector<StatCell>::iterator >::iterator new_it = cit + lit->second;
		new_it = cit;
		if (lit == LiftOver.end())	break;
		i--;
	}
}


// when set ctrl-stat, read position --> transition to original index --> then add
//  mei-type	position	ref-type	lift-over-length
void ReadMap::setCtrlStat(int & mei_type, std::string & het_index) // set hom + het + neg
{
	setGenomicIndexToControl();
	setCtrlLiftOver( het_index );

	std::string mei_str = getMeiTypeString( mei_type );
// set according to ref
  // hom
  	std::vector< std::vector<StatCell>::iterator > hom_nearby; // nearby hom & hom itself of ALL ctrl-mei
	std::ifstream het_index_file; het_index_file.open(het_index);
	CheckFileStatus(het_index_file);
    std::string line;
    int total_hom_skip = 0;
    while(std::getline(het_index_file, line)) {
		std::stringstream ss; ss << line;
		std::string current_mei_type;
    	std::getline(ss, current_mei_type, '\t');
    	std::string field;
    	std::getline(ss, field, '\t');
    	int st = stoi(field);
    	std::getline(ss, field, '\t');
    	bool type = stoi(field) > 0 ? 1 : 0;
    	if (type)  continue; // skip het
    	std::getline(ss, field, '\t');
    	int index = round( (stof(field) + st - win_len) / step_len );
    // do nearby
    	for(int i= index - win_len/step_len; i<= index + win_len/step_len; i++) {
    		if (CtrlLocationLink[i] != MergeCtrlMap.end())
    			hom_nearby.push_back(CtrlLocationLink[i]);
    	}
    // add to hom if it's the right mei_type
    	if (current_mei_type.compare(mei_str) == 0) {
    		if (CtrlLocationLink[index] == MergeCtrlMap.end()) {
    			total_hom_skip++;
    		}
    		else {
    			HomCtrl.push_back(CtrlLocationLink[index]);
    		}
    	}
    }
    het_index_file.close();
    
  // neg
	NegCtrl.resize(MergeCtrlMap.size());
  //copy first
	std::vector< std::vector<StatCell>::iterator >::iterator st_it = NegCtrl.begin();
	for( std::vector<StatCell>::iterator mit = MergeCtrlMap.begin(); mit != MergeCtrlMap.end(); mit++ ) {
		if ( std::find( hom_nearby.begin(), hom_nearby.end(), mit ) == hom_nearby.end()) { // not find in hom_nearby, add
			*st_it = mit; st_it++;
		}
	}
	NegCtrl.resize(st_it - NegCtrl.begin() + 1);

  // het: need to convert lift-over index to original index
	het_index_file.open(het_index);
	CheckFileStatus(het_index_file);
	int line_count = std::count(std::istreambuf_iterator<char>(het_index_file), std::istreambuf_iterator<char>(), '\n');
	het_index_file.close();
	het_index_file.open(het_index); CheckFileStatus(het_index_file);	
	HetCtrl.resize(line_count - HomCtrl.size());
	std::vector<HetCtrlCell>::iterator het_it = HetCtrl.begin();
	bool skip_this_rec = 0;
	std::vector<StatCell>::iterator current_hom_it;
	int skip_het_due_to_bad_hom = 0;
	int skip_het_single = 0;
	while(std::getline(het_index_file, line)) {
		std::stringstream ss; ss << line;
		std::string current_mei_type;
    	std::getline(ss, current_mei_type, '\t');
    	std::string field;
    	std::getline(ss, field, '\t');
    	int st = stoi(field);
    	std::getline(ss, field, '\t');
    	bool type = stoi(field) > 0 ? 1 : 0;
    	std::getline(ss, field, '\t');
    	if (type) { // het
    		if (skip_this_rec)	continue;
    		int index = round( (stof(field) + st - win_len) / step_len );
    		if (CtrlLocationLink[index] != MergeCtrlMap.end()) { // valid
    			if ( std::find( hom_nearby.begin(), hom_nearby.end(), CtrlLocationLink[index] ) != hom_nearby.end()) // hom nearby skip
    				skip_het_single++;
    			else { // add
    				int index = round( (stof(field) + st - win_len) / step_len );
    				std::vector< StatCell >::iterator current_neg_it = CtrlLocationLink[index];
    				het_it->log_table.resize(FIELD_COUNT);
    				std::vector<int>::iterator ot = current_hom_it->log_table.begin();
    				std::vector<int>::iterator nt = current_neg_it->log_table.begin();				
    				for( std::vector<int>::iterator it = het_it->log_table.begin(); it != het_it->log_table.end(); it++ ) {
    					*it = *ot + *nt; // 2 times but ok...need to convert later
    					ot++; nt++;
    				}  				    				
    				convertIntVecFromCountToLog( het_it->log_table );
    				het_it->hom_trace = current_hom_it;
    				het_it->neg_trace = current_neg_it;
    				het_it++;
    			}	
    		}
    		else // empty
				skip_het_single++;
    	}
    	else { // hom
    		if (current_mei_type.compare(mei_str) == 0) {
    			int index = round( (stof(field) + st - win_len) / step_len );
    			if (CtrlLocationLink[index] != MergeCtrlMap.end()) {
					current_hom_it = CtrlLocationLink[index];
    			}
    			else {
    				skip_this_rec = 1;
    				skip_het_due_to_bad_hom++;
    			}
    		}
    		else // wrong mei type
    			skip_this_rec = 1;
    	}
    	
    }
    HetCtrl.resize(het_it - HetCtrl.begin() + 1);
    het_index_file.close();
  
// finally clear hom_nearby
	hom_nearby.clear();
	
// and transform hom & neg stats from count to -log()
	convertStatCellVecFromCountToLog(HomCtrl);
	convertStatCellVecFromCountToLog(NegCtrl);
}

void ReadMap::setPerMergeCellPL (std::vector<MergeCell>::iterator & current_merge_it)
{
	std::map< std::vector<StatCell>::iterator, int > empty_exclude;
	current_merge_it->PL.resize(3);
	current_merge_it->PL[0] = getPerRefPerMergeCellPL(current_merge_it, HomCtrl, empty_exclude, 0);
	current_merge_it->PL[2] = getPerRefPerMergeCellPL(current_merge_it, NegCtrl, empty_exclude, 1);
	current_merge_it->PL[1] = getPerRefPerMergeCellHetPL(current_merge_it, HetCtrl, empty_exclude);
}


// std::map<std::string, std::vector<ReadMapCell> > RawReadMap --> std::vector<MergeCell> MergeMap & genomeLocationLink GenomeLocationLink;
void ReadMap::setDataFromRawCounts()
{
// initialize genomeLink
	for( std::map< std::string, std::vector<ReadMapCell> >::iterator it = RawReadMap.begin(); it != RawReadMap.end(); it++ ) {
		GenomeLocationLink[it->first].clear();
	}	
	
// merge RawReadMap to a whole ReadMapCell vector
	std::vector<ReadMapCell> sortedRawVec;
	int sum = 0;
	for( std::map< std::string, std::vector<ReadMapCell> >::iterator it = RawReadMap.begin(); it != RawReadMap.end(); it++ ) {
		sum += it->second.size();
	}
	sortedRawVec.resize(sum);
	std::vector<ReadMapCell>::iterator sort_it = sortedRawVec.begin();
	for( std::map< std::string, std::vector<ReadMapCell> >::iterator it = RawReadMap.begin(); it != RawReadMap.end(); it++ ) {
		for( std::vector<ReadMapCell>::iterator vit = it->second.begin(); vit != it->second.end(); vit++ ) {
			*sort_it = *vit;
			vit->stats.clear();
			sort_it++;
		}
	}
	RawReadMap.clear();	

// set duplicates
	std::vector<int> dup;
	setDupVecForReadMapCellVec( dup, sortedRawVec );
	
// set GenomicLocationLink
	MergeMap.resize( dup.size() - 1 );	
	std::vector<ReadMapCell>::iterator raw_it = sortedRawVec.begin();
	raw_it += dup[0]; // skip empty
	std::vector<MergeCell>::iterator merge_it = MergeMap.begin();
	for(std::vector<int>::iterator dp_it = ++dup.begin(); dp_it != dup.end(); dp_it++) {
		merge_it->count_table = raw_it->stats;
		GenomeLocationCell new_cell;
		for( int ct=0; ct < *dp_it; ct++ ) {
			new_cell.win_index = raw_it->win_index;
			new_cell.data = merge_it;
			GenomeLocationLink[raw_it->chr].push_back(new_cell);
		}
		merge_it++;
	}
	
// sort GenomeLocationLink & clean
	for( genomeLocationLinkIterator it = GenomeLocationLink.begin(); it != GenomeLocationLink.end(); it++ ) {
		std::sort( it->second.begin(), it->second.end(), sortGenomeLocationLink );
	}
	sortedRawVec.clear();
}



// let's do likelihood.......
// only calculate posterior when print likelihood out
// umbrella function
void ReadMap::SetDataPL (int & mei_type, std::string & work_dir, std::string & ctrl_dir, std::string & het_index)
{
// read ctrl-data first
	std::string meiTypeStr = getMeiTypeString(mei_type);
	std::string regular = ctrl_dir + "chr" + REF_CHR + ".proper.regular";
	std::string clip = ctrl_dir + "chr" + REF_CHR + ".proper." + meiTypeStr;
	std::string disc = ctrl_dir + "chr" + REF_CHR + ".disc." + meiTypeStr;
	bool status = setPerChrRawStats(RawCtrlMap, regular, clip, disc, 0);
	if (!status) {
		std::cerr << "ERROR: No regular or clip counts in REF_CHR " << REF_CHR << "!" << std::endl;
		exit(1);
	}

// set this mei_type ctrl first
	setCtrlStat( mei_type, het_index );
// needs to do the same thing for neg for 3 times but it's fast so should be fine...

// read raw-data
	for( std::map<std::string, int>::iterator mit = ChrLenTable.begin(); mit != ChrLenTable.end(); mit++ ) {
		std::string current_chr = mit->first;
		RawReadMap[current_chr].clear();
		regular = work_dir + "chr" + current_chr + ".proper.regular";
		clip = work_dir + "chr" + current_chr + ".proper." + meiTypeStr;
		disc = work_dir + "chr" + current_chr + ".disc." + meiTypeStr;
		bool status = setPerChrRawStats(RawReadMap[current_chr], regular, clip, disc, 1);
		if (!status) {
			std::cerr << "Skipped chr " << current_chr << " due to no regular & clip counts! Also removed it from ChrLenTable..." << std::endl;
			ChrLenTable.erase(current_chr);
		}
	}

// merge & sort	
	setDataFromRawCounts();
	
// set PL in merge without excluding any refs
	for( std::vector<MergeCell>::iterator mit = MergeMap.begin(); mit != MergeMap.end(); mit++ ) {
		setPerMergeCellPL (mit);
	}
/* do excluding for some cells:
	- scan REF_CHR
	- calculate excludes
*/
	std::vector< SpecialPLcell >::iterator sp_it = SpecialPLs.begin();
	int times = win_len / step_len;
	for( std::vector< GenomeLocationCell >::iterator gn_it = GenomeLocationLink[REF_CHR].begin(); gn_it != GenomeLocationLink[REF_CHR].end(); gn_it++ ) {
		int win_index = gn_it->win_index;
		std::vector<StatCell>::iterator ctrl_it = CtrlLocationLink[win_index];
		ctrl_it -= times;
		std::map< std::vector<StatCell>::iterator, int > exclude_stats;
		for(int i=0; i<=times; i++, ctrl_it++) {
			if ( ctrl_it->count < 10 ) {
				if( exclude_stats.find(ctrl_it) != exclude_stats.end())
					exclude_stats[ctrl_it]++;
				else
					exclude_stats[ctrl_it] = 1;
			}
		}
		if (exclude_stats.empty()) continue; // if empty, nothing should be done

		bool redo_hom = ExistOverlap(HomCtrl, exclude_stats);
		bool redo_het = ExistHetOverlap(HetCtrl, exclude_stats);
		bool redo_neg = ExistOverlap(NegCtrl, exclude_stats);
		if (!redo_hom && !redo_het && !redo_neg) continue; // no need to redo
		
		SpecialPLs.resize( SpecialPLs.size() + 1 );
		sp_it->PL.resize(3, -1);
				
		sp_it->PL[0] = redo_hom ? getPerRefPerMergeCellPL(gn_it->data, HomCtrl, exclude_stats, 0) : gn_it->data->PL[0];
		sp_it->PL[1] = redo_het ? getPerRefPerMergeCellHetPL(gn_it->data, HetCtrl, exclude_stats) : gn_it->data->PL[1];
		sp_it->PL[2] = redo_neg ? getPerRefPerMergeCellPL(gn_it->data, NegCtrl, exclude_stats, 1) : gn_it->data->PL[2];
	}
	
// convert special PL to merge cell. Redirect related GenomeLink to new merge cell
	specialMergeMap.resize(SpecialPLs.size());
	std::vector< MergeCell >::iterator ms_it = specialMergeMap.begin();
	for( sp_it = SpecialPLs.begin(); sp_it != SpecialPLs.end(); sp_it++ ) {
		ms_it->PL = sp_it->PL;
		std::vector< GenomeLocationCell >::iterator g_it = GenomeLocationLink[REF_CHR].begin();
		g_it += sp_it->index;
		g_it->data = ms_it;
		ms_it++;
	}
	SpecialPLs.clear();
}


using std::endl;
using std::vector;
void ReadMap::PrintToVcf(std::string & outVcfName, int & mei_type)
{
	std::string MeiType = getMeiTypeString(mei_type);
	
	std::ofstream outVcf; outVcf.open(outVcfName.c_str());
	CheckOutFileStatus(outVcf);

	std::string maxDepth = "NA";
	std::string minPosterior = "NA";
	
	time_t raw_time;
	time(&raw_time);
	
	outVcf << "##fileformat=VCFv4.1" << endl;
	outVcf << "##filedate=" << ctime(&raw_time) << endl;
	outVcf << "##source=LHMEI" << endl;
	outVcf << "##minDepth=1" << endl;
	outVcf << "##maxDepth=" << maxDepth << endl;
	outVcf << "##minMapQuality=0" << endl;
	outVcf << "##minPosterior=" << minPosterior << endl;
	
	outVcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Most Likely Genotype\">" << endl;
	outVcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
	outVcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Call Quality\">" << endl;
	outVcf << "##FORMAT=<ID=PL,Number=.,Type=Integer,Description=\"Genotype Likelihoods for Genotypes in Phred Scale, for 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, ...\">" << endl;
	
	outVcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << SampleName << endl;
	
	// chr pos id(.) type qual info gt:dp:gq:pl (samples)
	// need to implement info
	int half_win = win_len / 2;
	for( genomeLocationLinkIterator map_it = GenomeLocationLink.begin(); map_it != GenomeLocationLink.end(); map_it++) {
		for( vector<GenomeLocationCell>::iterator g_it = map_it->second.begin(); g_it != map_it->second.end(); g_it++) {
			std::string out_string;
			int position = (g_it->win_index * step_len + half_win);
			out_string = std::string("\t.\t.\t<") + MeiType + ">\t99\t" ;
			outVcf << map_it->first << "\t" << position << out_string << endl;
			outVcf << "."; // info
			outVcf << "\tGT:DP:GQ:PL\t";
			std::string gt = getGenotypeFromPL(g_it->data->PL);
			int dp = getSumOfVector(g_it->data->count_table); // #total reads
			int gq = getSumOfMeiReads(g_it->data->count_table); // #total mei reads
			vector<int> pl; pl.resize(3); setNormalizedPLfromPL(pl, g_it->data->PL);
			outVcf << gt << ":" << dp << ":" << gq << ":" << pl[0] << "," << pl[1] << "," << pl[2];
			
			outVcf << endl;
		}
	}
	
	outVcf.close();
}







































