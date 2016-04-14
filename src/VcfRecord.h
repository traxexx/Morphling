#ifndef VCFRECORD_H
#define VCFRECORD_H

#include <string>
#include <map>
#include <vector>

using std::vector;
using std::string;

class VcfRecord
{
  public:
  	VcfRecord( string & line, bool contain_gt );
  	string GetChr();
  	int GetPosition();
  	int GetMeiType();
  	int GetMeiSubtypeIndex( vector< vector<string> > & MEnames );
  	string GetInfoString( string & info_name );
    void SetAFinfo( vector<float> & prior, vector<int> & gt );
    void PrintNoEnding( std::ofstream & out );

  private:
    void setFloatInfoField( const char* str, float val );
    void setIntInfoField( const char* str, int val );
    void _checkInfoFieldExists( string & name );
    string _getAltName( string & field );

  	vector<string> common_fields;
  	int position;
  	int mei_type;
  	std::map<string, string> info_map;
};

bool IsHeaderLine(string & line);

#endif
