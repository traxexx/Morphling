#include "Match_Stat.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>	// distance
#include <algorithm>	// max_element, sort
#include <math.h>	// round, ceil


CtrlTable::CtrlTable(const char * reference, const char * mappability, const char * slice_bed, int win, int step, int chr)
{
	win_size = win;
	step_size = step;
	het_element_size_ = 10;
	
	ConstructCtrlTable(reference);
	ConstructMappability(mappability);
	ConstructGcFromFasta(reference);
	LabelCtrlTable(slice_bed);
	setLiftLengthInCtrlTable( chr );

/*	const char * out = "test/ptable";
	PrintMappabilityAndGc(out);
	exit(0);	
*/	
	SetHetTable(chr);
}

void CtrlTable::ConstructCtrlTable(const char * reference)
{
std::cout << "Constructing Control Table..." << std::endl;
/* construct raw map */
	std::string fai_name = std::string(reference) + ".fai";
	std::ifstream fai; fai.open(fai_name.c_str());
	
	if (!fai.is_open()) {
		std::cerr << "ERROR: Can't open " << reference << ".fai !" << std::endl; exit(1);
	}
	
	std::string line;
	
	while(std::getline(fai, line)) {
		std::stringstream ss; ss << line;
		std::string field;
		std::getline(ss, field, '\t');
		int chr = -1;
		if (std::all_of(field.begin(), field.end(), ::isdigit))
			chr = std::stoi(field);
		else {
			if (field.length() > 3) {
				std::string sub_field = field.substr(0,3);
				if (sub_field.compare("chr") == 0) {
					sub_field = field.substr(3);
					if (std::all_of(sub_field.begin(), sub_field.end(), ::isdigit))
						chr = std::stoi(sub_field);
				}
			}
		}	
		if (chr > 0) {
			std::getline(ss, field, '\t');
			int chr_len = std::stoi(field);
			ctrl_table_[chr].resize(int(ceil((chr_len - win_size)/ step_size)));
		}
	}
	
	fai.close();
	
std::cout << "\t" << "Map Size = " << ctrl_table_.size() << std::endl;
}

void CtrlTable::ConstructMappability(const char * mappability)
{
std::cout << "Constructing Mappability..." << std::endl;
/* Construct mappability by Step & Merge Step to Window*/
	std::ifstream mappability_bed; mappability_bed.open(mappability);
	
	if (!mappability_bed.is_open()) {
		std::cerr << "ERROR: Can't open " << mappability << std::endl; exit(1);
	}
	
	
	std::string line;
	int last_chr = -1;
	int diff;
	float last_align;
	std::vector<CtrlTableElement>::iterator it;
	while(std::getline(mappability_bed, line)) {
		std::stringstream ss; ss << line;
		std::string field;
		std::getline(ss, field, '\t');
		int chr = -1;
		if (std::all_of(field.begin(), field.end(), ::isdigit))
			chr = std::stoi(field);
		else {
			if (field.length() > 3) {
				std::string sub_field = field.substr(0,3);
				if (sub_field.compare("chr") == 0) {
					sub_field = field.substr(3);
					if (std::all_of(sub_field.begin(), sub_field.end(), ::isdigit))
						chr = std::stoi(sub_field);
				}
			}
		}
		
		if (chr < 0)	continue;
		std::getline(ss, field, '\t');
		int cord_start = std::stoi(field);
		std::getline(ss, field, '\t');
		int cord_end = std::stoi(field);
		std::getline(ss, field, '\t');
		float align = std::stof(field);
		if (chr != last_chr) {
			last_chr = chr;
			it = ctrl_table_[chr].begin();
			diff = cord_end;
			last_align = 0;
			while(diff >= step_size) {
				it->align = align;
				it++; diff -= step_size;
			}
		}
		else {
			int current_cord_len = cord_end - cord_start;
			int new_diff = diff + current_cord_len;
			if (new_diff >= step_size) {
				it->align = (diff * last_align + (step_size - diff) * align) / step_size;
				it++;
				new_diff -= step_size;
				while(new_diff >= step_size) {
					it->align = align;
					it++; new_diff -= step_size;
				}
			}
			else {
				last_align = (diff * last_align + current_cord_len * align) / new_diff;
			}
			diff = new_diff;
		}
	}
	mappability_bed.close();

 /* Merge */
 	int times = win_size / step_size;
	for(std::map<int, std::vector<CtrlTableElement> >::iterator chr_it = ctrl_table_.begin(); chr_it != ctrl_table_.end(); chr_it++) {
		std::vector<CtrlTableElement>::iterator vit_end = chr_it->second.end();
		vit_end -= times;
		std::vector<CtrlTableElement>::iterator vit = chr_it->second.begin();
		for(; vit != vit_end; vit++) {
			std::vector<CtrlTableElement>::iterator current_end = vit; current_end += times;
			float align = 0;
			for (std::vector<CtrlTableElement>::iterator current = vit; current != current_end; current++)
				align += current->align;
			align = align / times;
			vit->align = align;
		}
	}
}



void CtrlTable::ConstructGcFromFasta(const char * reference)
{
std::cout << "Constructing GC..." << std::endl;
/* Construct GC from fasta */
	std::ifstream fasta; fasta.open(reference);
	
	if (!fasta.is_open()) {
		std::cerr << "ERROR: Can't open " << reference << std::endl; exit(1);
	}
	
	std::string line;
	bool skip_this_chr = 0;
	
	std::string last_seq;
	std::vector<CtrlTableElement>::iterator it;
	
	while(std::getline(fasta, line)) {
		if (line.length() == 0) continue; // skip empty line
			
		if (line[0] == '>') {
			std::stringstream ss; ss << line;
			std::string pre_field;
			std::getline(ss, pre_field, ' ');
			int chr = -1;

			std::string field = pre_field.substr(1);
			if (std::all_of(field.begin(), field.end(), ::isdigit))
				chr = std::stoi(field);
			else {
				if (field.length() > 3) {
					std::string sub_field = field.substr(0,3);
					if (sub_field.compare("chr") == 0) {
						sub_field = field.substr(3);
						if (std::all_of(sub_field.begin(), sub_field.end(), ::isdigit))
							chr = std::stoi(sub_field);
					}
				}
			}
			if (chr > 0) {
				it =  ctrl_table_[chr].begin();
				skip_this_chr = 0;
				last_seq = std::string("");
			}
			else
				skip_this_chr = 1;
		}
		else {
			if (skip_this_chr) continue;
			last_seq = last_seq + line;
			while(last_seq.length() >= step_size) {
				int count = 0;
				for(int i=0; i<step_size; i++) {
					if (last_seq[i] == 'g' || last_seq[i] == 'G' || last_seq[i] == 'c' || last_seq[i] == 'C')
						count++;
				}
				it->gc = float(count) / step_size;
				last_seq = last_seq.substr(step_size); it++;
			}
		}
	}
	
	fasta.close();
}


void CtrlTable::LabelCtrlTable(const char * slice_bed)
{
std::cout << "Label & Init Het Table..." << std::endl;
	std::ifstream slice_bed_file; slice_bed_file.open(slice_bed);
	int line_count = std::count(std::istreambuf_iterator<char>(slice_bed_file), std::istreambuf_iterator<char>(), '\n');
	slice_bed_file.close(); slice_bed_file.open(slice_bed);
	
	Het_Table.resize(line_count);
	std::vector<HetTableElement>::iterator het_it = Het_Table.begin();
	std::string line;
	int half_win = win_size / 2;
	int invalid_hom = 0;
	int times = win_size / step_size;
	while(std::getline(slice_bed_file, line)) {
		std::stringstream ss; ss << line;
		std::string field;
		std::getline(ss, field, '\t');
		int chr = std::stoi(field);
		std::getline(ss, field, '\t');
		int cord = std::stoi(field);
		int win_index = int(round((cord - half_win) / step_size));
		std::getline(ss, field, '\t');
		int lift = stoi(field);
		std::string mei_type_str;
		std::getline(ss, mei_type_str, '\t');
		
		std::vector<CtrlTableElement>::iterator it = ctrl_table_[chr].begin() + win_index;
		if (win_index >= ctrl_table_[chr].size()) {
			std::cerr << "Warning: chr" << chr << ": " << cord << " >= chr" << chr << " length, skipped!" << std::endl;
			invalid_hom++; continue;
		}
		if (it->type < 0) {
			invalid_hom++; continue;
		}
		it->type = 1;
		it -= half_win;
		for (int i=0; i<=times; i++) {
			if(it->type == 0)
				it->type = -1;
			it++;
		}
		het_it->chr = chr;
		het_it->hom_index = win_index;
		het_it->hom_cord = cord;
		het_it->mei_type = mei_type_str;
		het_it->lift = lift;
		het_it++;
	}
	Het_Table.resize(line_count - invalid_hom);
	slice_bed_file.close();
}

// also set used = 0;
// notice here ctrl1 <= win < ctrl2, use lift-length of ctrl1
void CtrlTable::setLiftLengthInCtrlTable( int & chr )
{
	if (chr == 0) {
		std::cerr << "ERROR: Do not support multiple chr in this version!" << std::endl; exit(1);
	}
	
	std::vector<HetTableElement>::iterator het_it = Het_Table.begin();
	int lift = 0; // prev lift of het_it->hom_index
	bool reach_end = 0;
	for(std::vector<CtrlTableElement>::iterator it = ctrl_table_[chr].begin(); it != ctrl_table_[chr].end(); it++) {
		it->used = 0;
		int win_index = it - ctrl_table_[chr].begin();
		if ( win_index < het_it->hom_index || reach_end)
			it->lift = lift;
		else {
			while( win_index > het_it->hom_index && het_it != Het_Table.end() ) {
				het_it++;
			} // after this win_index <= het_it
			if (het_it == Het_Table.end())  reach_end = 1;
			std::vector<HetTableElement>::iterator prev_het_it = het_it;
			prev_het_it--;
			lift = prev_het_it->lift;
			it->lift = lift;
		}
	}
	
}




// find #HetSize similar ctrl window to the hom window
// implement:
// let's do not use a reference window more than once by labeling it after use. Anyway we have enough windows.
// 			These exclusion will include nearby windows

void CtrlTable::SetHetTable(int & chr)
{
std::cout << "Set Het Table, Length = " << Het_Table.size() << std::endl;
	bool single_chr = 0;
	if (chr > 0)
		single_chr = 1;
	else {
		std::cerr << "ERROR: Multi-chr not implemented!" << std::endl; exit(1);
	}
	
	int flank_steps = (win_size / step_size) * 16; // label around used (600 + 1kb)
	int jump_steps = (win_size / step_size) * 106; // make sure the distance >= 10kb between het windows
	for(std::vector<HetTableElement>::iterator het_it = Het_Table.begin(); het_it != Het_Table.end(); het_it++) {
		int fill = 0;
		het_it->match_index = new CordElement [het_element_size_];
		std::vector<CtrlTableElement>::iterator hom_it = ctrl_table_[het_it->chr].begin() + het_it->hom_index;
		float max_dist;
		int max_index;
		std::vector<float> dist_vec; dist_vec.resize(het_element_size_);
		std::vector<CordElement> index_vec; index_vec.resize(het_element_size_);
		
		std::vector<std::vector<CtrlTableElement>::iterator> it_vec; index_vec.resize(het_element_size_);
		
		std::map<int, std::vector<CtrlTableElement> >::iterator map_it = ctrl_table_.begin();
		for(; map_it != ctrl_table_.end(); map_it++) { // do not consider overlapping windows
			if (single_chr) {
				if (map_it->first != chr) continue;
			}
			int Times = 0;
			int LowBound = 2000000 / step_size;
			int HighBound = (ctrl_table_[map_it->first].size() * step_size - 2000000) / step_size;
			for(std::vector<CtrlTableElement>::iterator ctrl_it = map_it->second.begin(); ctrl_it != map_it->second.end(); ctrl_it++) {
				Times++;
				if (Times <= LowBound) continue;
				if (Times >= HighBound) break;
				if (ctrl_it->type == 0 && !ctrl_it->used) {
					float dist = fabs(ctrl_it->align - hom_it->align) + fabs(ctrl_it->gc - hom_it->gc);
					if (fill < het_element_size_) {
						dist_vec[fill] = dist;
						index_vec[fill].chr = map_it->first;
						index_vec[fill].win_index = (ctrl_it - map_it->second.begin());
						if (fill == 0) {
							max_dist = dist; max_index = 0;
						}
						else {
							if (dist > max_dist) {
								max_dist = dist; max_index = fill;
							}
						}
						fill++;
					}
					else {
						if (dist < max_dist) {
							dist_vec[max_index] = dist;
							index_vec[max_index].chr = map_it->first;
							index_vec[max_index].win_index = (ctrl_it - map_it->second.begin());
							std::vector<float>::iterator dist_vec_it = std::max_element(dist_vec.begin(), dist_vec.end());
							max_dist = *dist_vec_it;
							max_index = std::distance(dist_vec.begin(), dist_vec_it);
							if (ctrl_it - map_it->second.begin() < HighBound) {
								ctrl_it += jump_steps;
								Times += jump_steps;
							}
							else break;
						}
					}						
				}
			}
		}
	/* write to het table */
		for (int i=0; i < het_element_size_; i++) {
			het_it->match_index[i].chr = index_vec[i].chr;
			het_it->match_index[i].win_index = index_vec[i].win_index;
			std::vector<CtrlTableElement>::iterator ctrl_it = ctrl_table_[het_it->chr].begin();
			ctrl_it += index_vec[i].win_index;
			ctrl_it -= (flank_steps/2 + 1);
			for(int i=0; i<= flank_steps; i++) {
				ctrl_it->used = 1;
				ctrl_it++;
			}
		}
	}

}


void CtrlTable::PrintHetIndex(const char * out)
{
std::cout << "Printing Het Index..." << std::endl;
	std::ofstream het_index_file; het_index_file.open(out);
	for(std::vector<HetTableElement>::iterator het_it = Het_Table.begin(); het_it != Het_Table.end(); het_it++) {
		het_index_file << het_it->mei_type << "\t" << het_it->hom_cord << "\t0\t" << het_it->lift << std::endl;
		for(int i=0; i<het_element_size_; i++) {
			het_index_file << het_it->mei_type << "\t" << (het_it->match_index[i].win_index * step_size) << "\t1\t";
			std::vector<CtrlTableElement>::iterator vit = ctrl_table_[het_it->chr].begin();
			vit += het_it->match_index[i].win_index;
			het_index_file << vit->lift << std::endl;
		}
	}
	het_index_file.close();
}


/*********************** debug functions **************************/
void CtrlTable::PrintMappabilityAndGc(const char * out)
{
std::cout << "Printing Mappability & GC for debug..." << std::endl;
	std::string out_debug = std::string(out) + ".table";
	std::ofstream out_debug_file; out_debug_file.open(out_debug);

	for (std::map<int, std::vector<CtrlTableElement> >::iterator map_it = ctrl_table_.begin(); map_it != ctrl_table_.end(); map_it++) {
		for (std::vector<CtrlTableElement>::iterator vit = map_it->second.begin(); vit != map_it->second.end(); vit++) {
			out_debug_file << map_it->first << "\t" << vit->align << "\t" << vit->gc << "\t" << vit->type << std::endl;
		}
	}

	out_debug_file.close();
}







