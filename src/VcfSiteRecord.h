#ifndef VCFSITERECORD_H
#define VCFSITERECORD_H

#include <string>
#include <fstream>
#include <map>
#include "SamFile.h"

using std::string;

class VcfSiteRecord
{
public:
  VcfSiteRecord();
  	
  /** set methods ***/	
  void SetFromLine( string & line );

  // for memebrs
  void SetChrName( string & chr_name );
  void SetPosition( int center );
  void SetAltAllele( string & alt_str );
  
  /** get methods ***/
	string GetChromosome();
	int GetPosition();
  string GetMeiType();
	int GetVariantQuality();
	bool GetFilterPassOrNot();
  int GetEvidenceDepth();

  string GetInfoString( string & info_name );
  int GetInfoValue( string & info_name );
  void GetInfoPairValue( vector<int> & val, string & info_name );
  	  	
  /** print ***/	
	void PrintRecord( std::ofstream & out_vcf );

  private:
  int getFirstCharPosition( string & str, char c );
  int dosageFromGenotype( string & field );
  void parseGLfield();
  void _parseInfoField();

  void _checkInfoFieldExists( string & info_name );
  void _checkInfoFieldNumeric( string & f );

// regular fields
  string chr;
  int position;
  string id;
  string ref_allele;
  string alt_allele;
  int variant_quality;
  string filter;
  string info_field;
  std::map<string, string> info_map;
  string format_field;
  string gl_field; // store unparsed gl field
  int dosage;
  int depth;
};

/** other vcf record related functions **/
// if it's a header line in vcf, return TRUE
bool IsHeaderLine( string & line );

#endif
