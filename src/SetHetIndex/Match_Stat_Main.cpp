#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include "Match_Stat.h"
     
/*
get het based on b-dist to hom from ctrl
required:
hom-stats
print out:
.index: index, corresponding hom
.het.stat: index, het stats
*/

/* g++ src/Match_Stat_Main.cpp src/Match_Stat.h src/Match_Stat.cpp -o bin/SetReference -std=c++11 -g -O0
-slice_bed /net/assembly/saichen/MEI/simus/CEUTrio/NA12878/chr20/20.alu.bed -mappability /net/assembly/saichen/MEI/ref/smooth.wgEncode.100mer.bedGraph -reference /net/wonderland/home/saichen/reference/archive/hs37d5.fa -win 600 -step 100 -out test/match.index
*/

void MatchStatHelp();

int main(int argc, char * argv[])
{
	std::map<std::string, std::string> arg_map;
/* Required */
	arg_map["-slice_bed"] = std::string(""); // lift-overed
	arg_map["-mappability"] = std::string("");
	arg_map["-reference"] = std::string("");
	arg_map["-win"] = std::string("");
	arg_map["-step"] = std::string("");
	arg_map["-out"] = std::string("");
	arg_map["-chr"] = std::string("");
	
	if (argc <= 1) {
		MatchStatHelp(); return 0;
	}

	int i=1;
	while( i < argc) {
		std::string argv_str = std::string(argv[i]);
		if (arg_map.find(argv_str) != arg_map.end()) {
			i++;
			arg_map[argv_str] = argv[i];
			i++;
		}
		else {
			std::cerr << "Invalid option: " << argv[i] << "!" << std::endl;
			exit(1);
		}
	}

	const char * mappability = arg_map["-mappability"].c_str();
	const char * reference = arg_map["-reference"].c_str();
	const char * slice_bed = arg_map["-slice_bed"].c_str();
	int chr;
	if (arg_map["-chr"].length() == 0)
		chr = 0;
	else
		chr = stoi(arg_map["-chr"]);
	int win = stoi(arg_map["-win"]);
	int step = stoi(arg_map["-step"]);
	const char * out = arg_map["-out"].c_str();
	
	CtrlTable * ctrl_table = new CtrlTable(reference, mappability, slice_bed, win, step, chr);
	ctrl_table->PrintHetIndex(out);
}

void MatchStatHelp() {
	std::cerr << "Procrastination. Will do it soon..." << std::cerr;
}
