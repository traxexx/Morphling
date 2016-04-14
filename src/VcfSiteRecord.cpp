#include <sstream>
#include <algorithm> // std::all_of
#include "VcfSiteRecord.h"
#include "Globals.h"
#include "SamFile.h"
#include "bamUtility.h"

using std::stringstream;
using std::getline;
using std::cerr;
using std::endl;
using std::to_string;

// constructor
VcfSiteRecord::VcfSiteRecord():
	chr( string("") ),
	position( -1 ),
	ref_allele( string("") ),
	alt_allele( string("") ),
	variant_quality( -1 ),
	filter( string("") ),
	info_field( string("") ),
	format_field( string("") ),
	gl_field( string("") )
	{}


/*****set methods******/

// from file
void VcfSiteRecord::SetFromLine( string & line )
{
	if ( line.size() < 1 ) {// create an empty class
		cerr << "Warning: use an empty line to set vcf record!" << endl;
		return;
	}

	stringstream ss;
	ss << line;
	getline(ss, chr, '\t');
	string position_str;
	getline(ss, position_str, '\t');
	position = stoi( position_str );
	getline(ss, id, '\t'); // ref
	getline(ss, ref_allele, '\t'); // ref
	getline(ss, alt_allele, '\t'); // alt
	string qual_str;
	getline(ss, qual_str, '\t'); // qual
	variant_quality = stoi( qual_str );
	getline(ss, filter, '\t'); // filter
	getline(ss, info_field, '\t'); // info
	getline(ss, format_field, '\t'); // format
	getline(ss, gl_field, '\t'); // gl field
}

/*** set methods members ***/
void VcfSiteRecord::SetChrName( string & chr_name) {
	if ( !chr.empty() )
		cerr << "Warning: replace existing chr: " << chr << " to " << chr_name << "?" << endl;
	chr = chr_name;
}

void VcfSiteRecord::SetPosition( int center )
{
	position = center;
}

void VcfSiteRecord::SetAltAllele( string & alt_str )
{
	alt_allele = alt_str;
}

/** get methods for members **/
string VcfSiteRecord::GetInfoString( string & info_name )
{
	_checkInfoFieldExists( info_name );
	string istring = info_map[info_name];
	return istring;
}

int VcfSiteRecord::GetInfoValue( string & info_name )
{
	_checkInfoFieldExists( info_name );
	string val_str = info_map[ info_name ];
	_checkInfoFieldNumeric( val_str );
	int val = stoi(val_str);
	return val;
}

void VcfSiteRecord::GetInfoPairValue( vector<int> & val, string & info_name )
{
	_checkInfoFieldExists( info_name );
	string val_str = info_map[ info_name ];
	int comma = getFirstCharPosition( val_str, ',' );
	string f1 = val_str.substr(0, comma);
	string f2 = val_str.substr(comma + 1);
	_checkInfoFieldNumeric( f1 );
	_checkInfoFieldNumeric( f2 );
	val.resize(2);
	val[0] = stoi(f1);
	val[1] = stoi(f2);
}


string VcfSiteRecord::GetChromosome()
{
	if ( chr.empty() ) {
		cerr << "ERROR: chr name not defined at: " << position << endl;
		exit(1);
	}
	return chr;
}	
	
int VcfSiteRecord::GetPosition()
{
	if ( position < 0 ) {
		cerr << "ERROR: [VcfSiteRecord::GetPosition()] position out of range at: " << position << endl;
		exit(1);
	}
	return position;
}

string VcfSiteRecord::GetMeiType()
{
	int bk = -1;
	for( int i=alt_allele.size()-1; i>0; i-- ) {
		if ( alt_allele[i] == ':' ) {
			bk = i;
			break;
		}
	}
	if ( bk == -1 ) {
		cerr << "ERROR: [VcfSiteRecord::GetMeiType] vcf alt allele is not <INS:ME:XX>!\n" << endl;
		exit(40);
	}
	string str = alt_allele.substr( bk + 1, alt_allele.size() - bk - 2 );
	return str;
}

int VcfSiteRecord::GetVariantQuality()
{
	return variant_quality;
}

bool VcfSiteRecord::GetFilterPassOrNot()
{
	if ( filter.compare("PASS") == 0 )
		return 1;
	else
		return 0;
}

int VcfSiteRecord::GetEvidenceDepth()
{
	string str = "CLIP";
	int evi = GetInfoValue(str);
	str = "DISC";
	evi += GetInfoValue( str );
	str = "UNMAP";
	evi += GetInfoValue(str);
}


/*** inner parsers ***/
void VcfSiteRecord::parseGLfield()
{
	if ( format_field.empty() || gl_field.empty() ) {
		cerr << "ERROR: no format or gl field recorded: " << chr << ": " << position << endl;
		exit(1);
	}
// parse format field
	stringstream ss;
	ss << format_field;
	int gt_index = -1;
	int dp_index = -1;
	int fc = 0;
	string field;
	while( getline(ss, field, ':') ) {
		if ( field.compare("GT") == 0 )
			gt_index = fc;
		else if ( field.compare("DP") == 0 )
			dp_index = fc;
		fc++;
	}
	if ( fc == 0 ) {
		cerr << "ERROR: can't parse gl field: " << gl_field << ", maybe not separated by ':' ?" << endl;
		exit(1);
	}

	stringstream glss;
	glss << gl_field;
	fc = 0;
	while( getline( glss, field, ':' ) ) {
		if ( fc == gt_index ) {
			dosage = dosageFromGenotype( field );
		}
		else if ( fc == dp_index ) {
			depth = stoi( field );
		}
		fc++;
	}
}


void VcfSiteRecord::_parseInfoField()
{
	if ( info_field.empty() ) {
		cerr << "ERROR: no info field to parse at: " << chr << ": " << position << endl;
		exit(1);
	}
	// parse
	stringstream ss;
	ss << info_field;
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
			cerr << "ERROR: " << info_field << ", " << field << " doesn't have '=' sign!\n" << endl;
			exit(30);
		}
		string f1 = field.substr(0, equal_sign);
		string f2 = field.substr(equal_sign + 1);
		info_map[f1] = f2;
	}
}


void VcfSiteRecord::PrintRecord( std::ofstream & out_vcf )
{
// update all empty fields to dot	
	out_vcf << chr << "\t" << position << "\t.\t";
  // ref allele
	if ( ref_allele.empty() )
		out_vcf << ".\t";
	else
		out_vcf << ref_allele << "\t";
  // alt allele
	if ( alt_allele.empty() )
		out_vcf << ".\t";
	else
		out_vcf << alt_allele << "\t";
  // int
	out_vcf << variant_quality << "\t";
  // filter
	if ( filter.empty() )
		out_vcf << ".\t";
	else
		out_vcf << filter << "\t";
  // info
	if ( info_field.empty() )
		out_vcf << ".\t";
	else
		out_vcf << info_field << "\t";
  // these fields shouldn't be empty
	out_vcf << format_field << "\t" << gl_field << "\t" << endl;
}

int VcfSiteRecord::getFirstCharPosition( string & str, char c )
{
	int n = -1;
	for( int i=0; i<str.length(); i++ ) {
		if ( str[i] == c ) {
			n = i;
			break;
		}
	}
	if ( n==-1 ) {
		cerr << "ERROR: [VcfSiteRecord::getFirstCharPosition] cannot find " << c << " in " << str << endl;
		exit(60);
	}
	return n;
}

int VcfSiteRecord::dosageFromGenotype( string & field )
{
	bool fail = 0;
	if ( field[1] != '/' || field[1] != '|' )
		fail = 1;
	if (field.length() != 3)
		fail = 1;
	int dos = 0;
	if ( !fail ) {
		if ( field[0] == '1' )
			dos = 1;
		else if ( field[0] != '0' )
			fail = 1;
		if ( field[2] == '1' )
			dos++;
		else if ( field[2] != '0' )
			fail = 1;
	}
	if ( fail ) {
		cerr << "ERROR: GL field " << field << " is abnormal!\n" << endl;
		exit(42);
	}
	return dos;
}

void VcfSiteRecord::_checkInfoFieldExists( string & info_name )
{
	if ( info_map.empty() )
		_parseInfoField();
	if ( info_map.find( info_name ) == info_map.end() ) {
		cerr << "ERROR: [VcfSiteRecord::_checkInfoFieldExists] " << info_name << " doesn't exist in info field!\n" << endl;
		exit(50);
	}
}

void VcfSiteRecord::_checkInfoFieldNumeric( string & f )
{
	for( int i=0; i<f.size(); i++ ) {
		if (!isdigit(f[i])) {
			cerr << "ERROR: [VcfSiteRecord::_checkInfoFieldNumeric] numeric info field " << f << " is not numeric!\n" << endl;
			exit(51);
		}
	}
}


// if start with '#', then is header line
bool IsHeaderLine( string & line )
{
	if ( line[0] == '#' )
		return 1;
	else
		return 0;
}
