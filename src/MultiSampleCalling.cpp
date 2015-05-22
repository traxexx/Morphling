#include <iostream>
#include <fstream>
#include <math.h> // round

#include "MultiSampleCalling.h"
#include "SiteListMap.h"
#include "ReadMap.h"
#include "Utilities.h"
#include "Globals.h"
#include "VcfRecord.h" // VcfRecord
#include "ReGenotype.h"
#include "bamUtility.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::getline;

/* get consensus call set from multiple samples
 steps:
 	1. read through all vcfs, get candidate site list (q>=10 ? )
 		print list out
 	2. read through each vcf, get GL from bam on these sites.
 	3. bcftools -e for AF
 		make final calling based on posterior-likelihood
 	4. (optional) ad-hoc filtering?
*/

// main function for doing multi-calling
void MultiSampleCalling( Options * ptrMainOptions )
{	
// mei type	
	vector<string> mei_type;
	if ( ptrMainOptions->ArgMap["MeiType"].compare("-1") == 0 ) { // all mei type
		for( int i = 0; i <= 2; i++ )
			mei_type.push_back( std::to_string(i) );
	}
	else
		mei_type.push_back( ptrMainOptions->ArgMap["MeiType"] );

// chr
	vector<string> chrs;	
	if ( ptrMainOptions->ArgMap["Chr"].compare("-1") == 0 ) { // whole genome
		for( int i = 1; i <= 22; i++ )
			chrs.push_back( std::to_string(i) );
	}
	else
		chrs.push_back( ptrMainOptions->ArgMap["Chr"] );

// load sample list
	vector< vector<string> > SampleList;
	LoadSampleList( ptrMainOptions->ArgMap["SampleList"], SampleList );

// generate site list
	string work_dir = ptrMainOptions->ArgMap["WorkDir"];
	if (work_dir[ work_dir.size() - 1 ] != '/')
		work_dir += '/';
	string site_list_prefix = work_dir + "/Site-list";
	
	for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) {
		for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
			string vcf_suffix = string("-") + *current_chr + "." + *mt + ".vcf";
			SiteList sList( SampleList, vcf_suffix );
			string site_list_out_name = site_list_prefix + vcf_suffix;
			sList.Print( site_list_out_name );
		}
	}

/* generate parallel command if needed
	if ( ptrMainOptions->OptMap["Parallel"] ) { // print commands to work_dir/sub_commands/sample.gt.cmd
		string ExeFullName = GetExeFullName(); // secured last is '/'
		if ( ExeFullName.length() <= 4 ) {
			std::cerr << "ERROR: can't get full name of LHMEI!" << std::endl;
			exit(1);
		}
		for( vector< vector<string> >::iterator item_ptr = vcf_list.begin(); item_ptr != vcf_list.end(); item_ptr++ ) {
			string out_cmd_name = work_dir + "/" + (*item_ptr)[0] + ".gt.cmd";
			ofstream out_cmd;
			out_cmd.open( out_cmd_name.c_str() );
			CheckOutFileStatus( out_cmd, out_cmd_name );
			out_cmd << ExeFullName << " ReGenotype -Sample " << (*item_ptr)[0] << " -Bam " << (*item_ptr)[1] << " -Vcf " << (*item_ptr)[2];
			out_cmd << " -CtrlDir " << (*item_ptr)[3] << " -Win " << WIN << " -SiteList " << site_list_name << endl;
			out_cmd.close();
		}
		return;
	}
*/
	
// do single-thread: load site list first
	vector<RefSeq*> REF_SEQ;
	InitializeMeiSeqRef( REF_SEQ, ptrMainOptions->ArgMap["MElist"] );
	vector< map<string, vector<int> > > SimplifiedSiteList; // mt -> chr -> siteListVec
	for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) {
		int mt_index = mt - mei_type.begin();
//		vector< map<string, vector<int> > >::iterator ms_it = SimplifiedSiteList.begin();
//		ms_it += mt_index;
		for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
			string site_list_name = site_list_prefix + "-" + *current_chr + "." + *mt;
			LoadSiteList( SimplifiedSiteList[ mt_index ][ *current_chr ], site_list_name );
		}
	}

// load het index
	AllHetIndex allHetIndex;
	SetAllHetIndex( ptrMainOptions->ArgMap["HetIndex"], allHetIndex );

// re-genotype each vcf
	for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) {
		vector<RefSeq*> single_REF_SEQ;
		single_REF_SEQ.resize(3, NULL);
		int mei_index = std::stoi( *mt );
		single_REF_SEQ[ mei_index ] = REF_SEQ[ mei_index ];
		string mt_suffix = string(".") + *mt + ".vcf";
		for( vector< vector<string> >::iterator item_ptr = SampleList.begin(); item_ptr != SampleList.end(); item_ptr++ ) {
			string sample_name = (*item_ptr)[0];
			string bam_name = (*item_ptr)[1];
			string vcf_prefix = (*item_ptr)[2];
			string ctrl_dir = (*item_ptr)[3];		
			
			string proper_prefix = ctrl_dir + "/proper-" + REF_CHR;
			string disc_prefix = ctrl_dir + "/disc-" + REF_CHR;
			RefStats rstats( proper_prefix, disc_prefix, mei_index, allHetIndex );
			
			// open bam & bai
			SamFile samIn;
			SamFileHeader samHeader;
			OpenBamAndBai( samIn, samHeader, bam_name );
			for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
				bool section_status = samIn.SetReadSection( current_chr->c_str() );
				if ( !section_status ) {
					cerr << "ERROR: Unable to set read section at: " << *current_chr << ", in file: " << bam_name << endl;
					exit(1);
				}
				string vcf_suffix = string("-") + *current_chr + mt_suffix;
				string in_vcf_name = vcf_prefix + vcf_suffix;
				string out_vcf_name = work_dir + "/refined-" + sample_name + vcf_suffix;
				ReGenotypeSingleVcf(  SimplifiedSiteList[ mei_index ][*current_chr], single_REF_SEQ, rstats, samIn, samHeader, in_vcf_name, ctrl_dir, out_vcf_name);
			}
			samIn.Close();
		}
	}
}

// get GL for single sample based on site list, then print out to vcf
// do it by chr
void ReGenotypeSingleVcf( vector<int> & siteVec, vector<RefSeq*> & REF_SEQ, RefStats & rstats, SamFile & samIn, SamFileHeader & samHeader, string & in_vcf_name, string & ctrl_dir, string & out_vcf_name )
{
	int MaxKeepDist = WIN / 3; // max dist to use the GL instead of getting a new one
	ifstream in_vcf;
	in_vcf.open( in_vcf_name.c_str() );
	CheckInputFileStatus( in_vcf, in_vcf_name.c_str() );
	ofstream out_vcf;
	out_vcf.open( out_vcf_name );
	CheckOutFileStatus( out_vcf, out_vcf_name.c_str() );
	
	string line;
	vector<int>::iterator site_ptr = siteVec.begin();
	while( getline( in_vcf, line ) ) {
		VcfRecord vcf_rec;
		vcf_rec.SetFromLine( line );
		int position = vcf_rec.GetPosition();
		int dist = *site_ptr - position; // ref - current_read_position
		while ( dist <= MaxKeepDist ) {
			while ( dist < -MaxKeepDist ) { // We've passed some ref positions without genotyping them
				VcfRecord implemented_rec;
				string chr = vcf_rec.GetChromosome();
				ResetVcfRecordFromBam( implemented_rec, rstats, REF_SEQ, chr, *site_ptr, samIn, samHeader );
				implemented_rec.PrintRecord( out_vcf );
				site_ptr++;
				if ( site_ptr == siteVec.end() ) { // no more site to be re-genotyped
					out_vcf.close();
					in_vcf.close();
					return;
				}
				dist = *site_ptr - position;
			}
			vcf_rec.PrintRecord( out_vcf );
			site_ptr++;
			if ( site_ptr == siteVec.end() ) { // no more site to be re-genotyped
				out_vcf.close();
				in_vcf.close();
				return;
			}
			dist = *site_ptr - position;
		}
		if ( dist > MaxKeepDist ) // haven't reached ref site, continue reading
			continue;
	}
}

// read sample list & parse
void LoadSampleList( string & sample_list_name, vector<vector<string> > & SampleList )
{
	string line;
	ifstream sample_list;
	sample_list.open( sample_list_name.c_str() );
	CheckInputFileStatus( sample_list, sample_list_name.c_str() );
	while( getline( sample_list, line ) ) {
		std::stringstream ss;
		ss << line;
		vector<string> current_info;
		string field;
		while( getline( ss, field, '\t') )
			current_info.push_back(field);
		if ( current_info.size() != 4 ) {
			cerr << "ERROR: current line:\n " << line << "is not a regular sample list line!\n sample list should be 4 fields: sample-name bam vcf-prefix ctrl-dir." << endl;
			exit(1);
		}
		SampleList.push_back( current_info );
	}
	sample_list.close();
}


/* main function for doing re-genotype
ReGenotypeSingleSample( Options * ptrMainOptions )
{



}

*/
// invoke bcftools for AF

// ad hoc filter
