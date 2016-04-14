#include "VcfRecord.h"
#include "morphError.h"
#include "MeiSeq.h"
#include <sstream>
#include <fstream>
#include <iostream>

using std::stringstream;

VcfRecord::VcfRecord( string & line, bool contain_gt )
{
	if (contain_gt)
		common_fields.resize(10);
	else
		common_fields.resize(8);
	stringstream ss;
	ss << line;
	for(int i=0; i<common_fields.size(); i++)
		std::getline( ss, common_fields[i], '\t' );

	// assign
	position = stoi(common_fields[1]);
	string alt_name = _getAltName( common_fields[4]  );
	mei_type = GetMeiIndexFromName( alt_name );
	// parse info to info_map
	ss.clear();
	ss << common_fields[7];
	string field;
	while( getline(ss, field, ';') ) {
		int equal_sign = -1;
		for(int i=0; i<field.size(); i++) {
			if ( field[i] == '=' ) {
				equal_sign = i;
				break;
			}
		}
		if (equal_sign == -1) {
			string str = "info field: " + field + " doesn't have '=' sign!";
			morphError(str, 30);
		}
		string f1 = field.substr(0, equal_sign);
		string f2 = field.substr(equal_sign + 1);
		info_map[f1] = f2;
	}
}

string VcfRecord::GetChr()
{
	return common_fields[0];
}

int VcfRecord::GetPosition()
{
	return position;
}

int VcfRecord::GetMeiType()
{
	return mei_type;
}

int VcfRecord::GetMeiSubtypeIndex( vector< vector<string> > & MEnames )
{
	string str = "SUB";
	string mname = GetInfoString( str );
	int index = -1;
	for(int i=0; i<MEnames[mei_type].size(); i++) {
		if ( mname.compare( MEnames[mei_type][i] ) == 0 ) {
			index = i;
			break;
		}
	}
	if (index == -1) {
		string str = "subtype" + mname + " not valid!";
		morphError(str, 21);
	}
	return index;
}

string VcfRecord::GetInfoString( string & info_name )
{
	_checkInfoFieldExists( info_name );
	string istring = info_map[info_name];
	return istring;
}

// set af, ac
void VcfRecord::SetAFinfo( vector<float> & prior, vector<int> & gt )
{
	if (prior.size() != 3 || gt.empty())
		morphError("[VcfRecord::SetAFinfo] Irregular prior or gt info");
	// af
	float af = prior[1]/2 + prior[2];
	setFloatInfoField( "AF", af );
	// ac
	int ac = 0;
	for(int i=0; i<gt.size(); i++) {
		if (gt[i] == -1)
			continue;
		ac += gt[i];
	}
	setIntInfoField("AC", ac);
}


void VcfRecord::PrintNoEnding( std::ofstream & out )
{
	for( int i=0; i<8; i++ ) {
		if (i>0)
			out << '\t';
		out << common_fields[i];
	}
}


void VcfRecord::setFloatInfoField( const char* str, float val )
{
	info_map[string(str)] = std::to_string(val);
	common_fields[7] += ";" + string(str) + "=" + std::to_string(val);
}

void VcfRecord::setIntInfoField( const char* str, int val )
{
	info_map[string(str)] = std::to_string(val);
	common_fields[7] += ";" + string(str) + "=" + std::to_string(val);
}


void VcfRecord::_checkInfoFieldExists( string & name )
{
	if (info_map.find(name) != info_map.end())
		return;
	string str = "info field: " + name + " does not exits";
	morphError(str, 41);
}

string VcfRecord::_getAltName( string & field )
{
	std::stringstream ss;
	ss << field;
	string sub;
	std::getline( ss, sub, ':' );
	std::getline( ss, sub, ':' );
	std::getline( ss, sub, ':' );
	sub = sub.substr( 0,sub.size()-1 );
	return sub;
}

bool IsHeaderLine(string & line)
{
	if (line.empty())
		morphError("vcf contains empty line");
	if (line[0] == '#')
		return 1;
	return 0;
}


