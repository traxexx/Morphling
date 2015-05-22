#include "Options.h"

#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <iterator>

Options::Options(int argc, char * argv[], std::string & arg_string, std::string & dummies)
{
  // set options first from ArgString & Dummies
	OptMap[std::string("verbose")] = 0; // verbose is always here
  
	std::stringstream ss; ss  << arg_string;
	std::string token;
	while(std::getline(ss, token, ';')) {	
		std::stringstream tss; tss << token;
		std::string raw_option;
		std::string val_default;
		std::getline(tss, raw_option, '=');
		std::getline(tss, val_default, '=');
		if (raw_option[0] == '-') {
			if (raw_option[1] == '-') {
				std::cerr << "No --option allowed in this setting. Exit!" << std::endl; exit(1);
			}
			std::string option = raw_option.substr(1);
			if (val_default.length() == 1 && val_default[0] == ' ')
				ArgMap[option].clear();
			else
				ArgMap[option] = val_default;
		}
		else {
			std::cerr << "Invalid inline argument: " << raw_option << std::endl; exit(1);
		}
	}
	ss.clear(); token.clear();
	
	ss  << dummies;
	while(std::getline(ss, token, ';')) {	
		std::stringstream tss; tss << token;
		std::string raw_option;
		std::getline(tss, raw_option, ':');
		if (raw_option[0] == '-' && raw_option[1] == '-') {
			std::string option = raw_option.substr(2);
			OptMap[option] = 0;
		}
		else {
			std::cerr << "Invalid inline option: " << raw_option << std::endl; exit(1);
		}
	}
	ss.clear(); token.clear();

  // read OptMap first
	int current_argc = 1;
	while(current_argc < argc) {
		std::string argv_str = std::string(argv[current_argc]);
		if (argv_str[0] == '-' && argv_str[1] == '-') {
			std::string option = argv_str.substr(2);
			if (OptMap.find(option) != OptMap.end())
				OptMap[option] = 1;
			else {
				std::cerr << "Invalid OptMap option: " << argv_str << std::endl; exit(1);
			}
		}
		current_argc++;
	}

	// set dummies
	ss << dummies;
	while(std::getline(ss, token, ';')) {
		std::stringstream tss; tss << token;
		std::string raw_mother_option;
		std::getline(tss, raw_mother_option, ':');
		if ( !(raw_mother_option[0] == '-' && raw_mother_option[1] == '-') ) {
			std::cerr << "When set dummies, mother option should be --option. Exit!" << std::endl;
			exit(1);
		}
		std::string mother_option = raw_mother_option.substr(2);
		if (OptMap.find(mother_option) != OptMap.end()) {
			if (!OptMap[mother_option]) continue; // if --option does not exist, skip
		}
		else {
			std::cerr << "Mother option doesn't exist in OptMap: --" << mother_option << std::endl; exit(1); 
		}
		std::string field;
		while(std::getline(tss, field, ',')) {
			if (field[0] == '-') {
				if (field[1] != '-') {
					std::stringstream field_ss; field_ss << field;
					std::string raw_child_option;
					std::getline(field_ss, raw_child_option, '=');
					std::string child_option = raw_child_option.substr(1);
					std::string val_default;
					std::getline(field_ss, val_default, '=');
					if (val_default.length() == 1 && val_default[0] == ' ')
						ArgMap[child_option].clear();
					else
						ArgMap[child_option] = val_default;
				}
				else {
					std::string child_option = field.substr(2);
					if (OptMap.find(child_option) != OptMap.end()) {
						std::stringstream field_ss; field_ss << field;
						std::string raw_child_option;
						std::getline(field_ss, raw_child_option, '=');
						std::string child_option = raw_child_option.substr(1);
						std::string val_default;
						std::getline(field_ss, val_default, '=');
						if (val_default.length() == 1 && (val_default[0] == '1' || val_default[0] == '0'))
							OptMap[child_option] = val_default[0] == '1' ? 1 : 0;
						else {
							std::cerr << "Invalid OptMap child default value: " << field << std::endl; exit(1);
						}
					}
					else {
						std::cerr << "Invalid OptMap child option: " << field << std::endl; exit(1);
					}
				}
			}
			else {
				std::cerr << "Invalid child option: " << field << std::endl; exit(1);
			}
		}
	}
	ss.clear(); token.clear();
	
// read ArgMap
	current_argc = 1;
	while(current_argc < argc) {
		std::string argv_str = std::string(argv[current_argc]);
		if (argv_str[0] == '-') {
			if (argv_str[1] != '-') {
				std::string option = argv_str.substr(1);
				current_argc++;
				if (ArgMap.find(option) != ArgMap.end())
					ArgMap[option] = argv[current_argc];
				else {
					std::cerr << "Argument do not exist: " << argv_str << std::endl; exit(1);
				}
			}
		}
		else {
			std::cerr << "Invalid argument: " << argv_str << std::endl; exit(1);
		}
		current_argc++;
	}
	ss.clear(); token.clear();
	
	verifyOptions();
	if (OptMap[std::string("verbose")])
		displayOptions();
}


void Options::verifyOptions()
{
	for( std::map<std::string, std::string>::iterator it_ArgMap = ArgMap.begin(); it_ArgMap != ArgMap.end(); it_ArgMap++) {
		if (it_ArgMap->second.length() == 0) {
			std::cerr << "Missing required argument: -" << it_ArgMap->first << ", exit! " << std::endl; exit(1);
		}
	}
}


/* check validity of options */
void Options::displayOptions()
{
	std::cout << std::endl;
	std::cout << "Display Options: " << std::endl;
	std::cout << " Options Map: " << std::endl;
	for( std::map<std::string, bool>::iterator it_OptMap = OptMap.begin(); it_OptMap != OptMap.end(); it_OptMap++) {
		std::cout << "\t" << it_OptMap->first << ": ";
		if (it_OptMap->second)
			std::cout << "Yes" << std::endl;
		else
			std::cout << "No" << std::endl;
	}
	std::cout << " Arguments Map: " << std::endl;
	for( std::map<std::string, std::string>::iterator it_ArgMap = ArgMap.begin(); it_ArgMap != ArgMap.end(); it_ArgMap++) {
		std::cout << "\t" << it_ArgMap->first << ": " << it_ArgMap->second << std::endl;
	}
	std::cout << std::endl;
}
