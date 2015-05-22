#include "Debugs.h"
#include <iostream>

void PrintValidRawMapCounts(ReadMap * Rmap)
{
	for(std::map<std::string, std::vector<ReadMapCell> >::iterator im = Rmap->RawMap.begin(); im != Rmap->RawMap.end(); im++) {
		for(std::vector<ReadMapCell>::iterator it = im->second.begin(); it != im->second.end(); it++) {
			if (it->second.size() > 0) {
				std::cout << "chr" << im->first << ", ";
				std::cout << it->first << " -->";
				for(std::vector<int>::iterator ia = it->second.begin(); ia != it->second.end(); ia++) {
					std::cout << " " << *ia;
				}
				std::cout << std::endl;
			}
		}
	}
}