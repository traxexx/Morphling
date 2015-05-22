#ifndef VCFRECORD_H
#define VCFRECORD_H

#include <string>
#include <fstream>
#include "SamFile.h"
#include "OriginalStats.h"

using std::string;

class VcfRecord
{
  public:
  	VcfRecord();
  	~VcfRecord();
  	
  /** set methods ***/	
  	void SetFromLine( string & line );
//  void SetFromInputParameters( string & chr_name, int breakp, string & ref_al, string & alt_al, int quality, string & filter_str, string & info_str, string & format_str, string & gl_str );
	void SetFromGenomeLocationCell( GenomeLocationCell & glc );
  // for memebrs
  	void SetChrName( string & chr_name );
  	void SetPosition( int center );
  
  /** get methods ***/
	string GetChromosome();
	int GetPosition();
	int GetVariantQuality();
	bool GetFilterPassOrNot();
	int GetVariantEnd();
	int GetVariantDepth();
	int GetDosage();
	int GetEvidenceDepth();
	
  /** update from merge cell (data ) ***/
  	void UpdateFromMergeCellPtr( MergeCellPtr & ptr );	
	
  /** anchor update methods ***/
  	bool UpdateByRankAnchor( GlcPtr & anchor, GlcPtr & new_anchor );	
  	
  /** break point ***/
  	void SetBreakPointAndCIFromBam( SamFile & currentSam, SamFileHeader & currentSamHeader );	
  	
  /** print ***/	
	void PrintRecord( ofstream & out_vcf );

  private:
  	void parseGLfield();
  	void parseInfoField();
  	
  /** used in update-by-rank-anchor ***/	
  	void updateRecWithNewAnchorIndex( int win_index );
  	void updateRecWithEqualAnchorIndex( int win_index );
  	void updateBothEnd();
  	void updateFilter();
  	void updateInfoField();
  	void updateEvidenceFromReadCount();
  	void updateEvidence();
  	void updateGLfield();
  	
  /** filters **/
  	bool depthQC();	
 
// regular fields
  	string chr;
  	int position;
  	string ref_allele;
  	string alt_allele;
  	int variant_quality;
  	string filter;
  	string info_field;
  	string format_field;
  	string gl_field; // store unparsed gl field

// indicate if field parsed 	
  	bool info_parsed;
  	bool gl_parsed; // if gl parsed, format also parsed
// field after parse	
  	int depth;
  	int dosage;
  	bool breakp_refined;
  	
// LHMEI fields
	bool both_end; // a filter
	int variant_end; // info
	int ci_low;
	int ci_high;
	int clip_count;
	int disc_count;
	int unmap_count;
	int evidence;
	int win_count; // info
	vector<int> read_counts;
	vector<float > GL;
};

/** other vcf record related functions **/
// if it's a header line in vcf, return TRUE
bool IsInfoLine( string & line );

#endif
