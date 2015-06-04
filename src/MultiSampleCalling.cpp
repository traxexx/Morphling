#include <iostream>
#include <fstream>
#include <math.h> // round
#include <sstream>

#include "MultiSampleCalling.h"
#include "SiteListMap.h"
#include "ReadMap.h"
#include "Utilities.h"
#include "SetReadMapGlobals.h"
#include "Globals.h"
#include "VcfRecord.h" // VcfRecord
#include "ReGenotype.h"
#include "bamUtility.h"
#include "ConsensusVcf.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::getline;

/* multi-sample calling steps:
 	1. read through all vcfs, get candidate site list (q>=10 ? )
 		print list out
 	2. read through each vcf, get GL from bam on those un-genotyped sites. ---> ReGenotype()
 	3. make final calling based on posterior-likelihood
 	4. (optional) ad-hoc filtering?
*/

// main function for doing multi-calling
void MultiSampleCalling( Options * ptrMainOptions )
{	
// set globals
	SetGenotypeGlobalOptions( ptrMainOptions );
	SetGenotypeGlobalParameters( ptrMainOptions );

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
	
// set global n_sample
	int NSAMPLE = (int)SampleList.size();
	cout << "Toal sample size = " << NSAMPLE << "." << endl;

// generate site list
	string work_dir = ptrMainOptions->ArgMap["WorkDir"];
	string mkp_cmd = string("mkdir -p ") + work_dir;
	ExecuteCmd( mkp_cmd );
	if (work_dir[ work_dir.size() - 1 ] != '/')
		work_dir += '/';
	string rg_dir = work_dir + "re-genotype/";
	string sl_dir = work_dir + "sites/";
	string final_dir = work_dir + "final/";
	mkp_cmd = string("mkdir -p ") + rg_dir + " " + sl_dir + " " + final_dir;
	ExecuteCmd( mkp_cmd );
	
	string chr_str;
	for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
		chr_str += " " + (*current_chr);
	}
	cout << "Generating site list for chr:" << chr_str << " ..." << endl;
	for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) {
		for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
			string site_list_base_name = "SiteList-" + *current_chr + "." + *mt;
			if ( ExistDoneFile( sl_dir, site_list_base_name.c_str() ) ) {
				cerr << "Warning: [MultiSampleCalling()] exist done file: " << site_list_base_name << ", skip this chr & mei-type!" << endl;
				continue;
			}
			string vcf_suffix = string("Hits-") + *current_chr + "." + *mt + ".vcf";
			string site_list_out_name = sl_dir + site_list_base_name;
			SiteList sList( SampleList, vcf_suffix );
			sList.Print( site_list_out_name );
			GenerateDoneFile( sl_dir, site_list_base_name.c_str() );
		}
	}
	cout << "  All site lists completed!" << endl;

// check if parallel option toggled
	ofstream paraFile;
	vector<RefSeq*> REF_SEQ;
	vector< map<string, vector<int> > > SimplifiedSiteList; // mt -> chr -> siteListVec
	
	string paraFileName = work_dir + "Parallel-Morphling.cmd";
	if ( PARALLEL ) { // open command file, ready to print cmd
		cout << "--Parallel option toggled..." << endl;
		paraFile.open( paraFileName.c_str() );
		CheckOutFileStatus( paraFile, paraFileName.c_str() );
	}
	else { // load site list for making final vcf
		cout << "Doing single-thread re-genotyping ..." << endl;
	  // load site list
		InitializeMeiSeqRef( REF_SEQ, ptrMainOptions->ArgMap["MElist"] );
		SimplifiedSiteList.resize( mei_type.size() );
		for( int i = 0; i < int(mei_type.size()); i++ ) {
			for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ )
				SimplifiedSiteList[i][ *current_chr ].clear();
		}
		for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) {
			int mt_index = mt - mei_type.begin();
			for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
				string site_list_name = sl_dir + "SiteList-" + *current_chr + "." + *mt;
				LoadSiteList( SimplifiedSiteList[ mt_index ][ *current_chr ], site_list_name );
			}
		}
		cout << "  Site list loaded!" << endl;
	}

	
// do re-genotype or print to parallel cmd
	for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) { // per mei type
		int mei_index = std::stoi( *mt );
		if ( !PARALLEL)
			cout << "Re-genotyping MEI type: " << *mt << " ..." << endl;
		for( vector< vector<string> >::iterator item_ptr = SampleList.begin(); item_ptr != SampleList.end(); item_ptr++ ) { // per sample
			if ( !PARALLEL )
				cout << "  Working on sample " << (*item_ptr)[0] << " ..." << endl;
			for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) { // per chr
				if ( PARALLEL ) {
					string site_list_name = sl_dir + "SiteList-" + *current_chr + "." + *mt;
					PrintParallelCommand( paraFile, *item_ptr, site_list_name, ptrMainOptions, *mt, *current_chr, rg_dir );
				}
				else {
					string sm = (*item_ptr)[0] + "." + *mt;
					if ( !ExistDoneFile( rg_dir, sm.c_str() ) ) {
						ReGenotypeSingleVcf( REF_SEQ, SimplifiedSiteList[ mei_index ][*current_chr], *item_ptr, rg_dir, *mt, *current_chr );
						GenerateDoneFile( rg_dir, sm.c_str() );
					}
				}
			}
			if ( !PARALLEL )
				cout << "  Sample " << (*item_ptr)[0] << " re-genotype finished." << endl;
		}
	}

// if parallel, do not proceed unitl all commands are executed
	if ( PARALLEL ) {
		paraFile.close();
		cout << "Parallel Morphling commands printed to " << paraFileName << ". Finish all these commands in it before running Morphling Genotypeagain!" << endl;
		return;
	}


// generate consensus vcf per chr: read each sample all mei-type then print the whole chr out
	cout << "Generating final vcf by merging all re-genotyped vcfs..." << endl;
	for( vector<string>::iterator current_chr = chrs.begin(); current_chr != chrs.end(); current_chr++ ) {
		ConsensusVcf FinalVcf( NSAMPLE, WIN, *current_chr );
		for( vector<string>::iterator mt = mei_type.begin(); mt != mei_type.end(); mt++ ) {
			string site_list_name = sl_dir + "SiteList-" + *current_chr + "." + *mt;
			FinalVcf.InitializeSdataFromSiteList( site_list_name );
			FinalVcf.SetSampleList( SampleList );
			int mei_index = stoi( *mt );
			FinalVcf.SetAltAlleleByMeiType( mei_index );
			FinalVcf.SetFasta( ptrMainOptions->ArgMap["GenomeFasta"] );
			int sp_index = 0;
			for( vector< vector<string> >::iterator item_ptr = SampleList.begin(); item_ptr != SampleList.end(); item_ptr++, sp_index++ ) {
				string sample_name = (*item_ptr)[0];
				string vcf_prefix = (*item_ptr)[2];
				string vcf_suffix = string("-") + *current_chr + string(".") + *mt + ".vcf";
				string in_vcf_name = vcf_prefix + vcf_suffix;
				string rg_vcf_name = rg_dir + "refined-" + sample_name + vcf_suffix;
				FinalVcf.AddFromSingleVcf( sp_index, rg_vcf_name );
			}
			FinalVcf.MergeData();
		}
		FinalVcf.Polish();
		string final_vcf_name = work_dir + "final/final." + (*current_chr) + ".vcf";
		FinalVcf.Print( final_vcf_name );
	}

// finish	
	cout << "Morphling Genotype finished with no error reported. Check final output at: " << final_dir << endl;
}

// child-function for Morphling reGenotype
void ReGenotype( Options * ptrMainOptions )
{
// set globals
	SetGenotypeGlobalOptions( ptrMainOptions );
	SetGenotypeGlobalParameters( ptrMainOptions );

// load site list
	vector<int> siteVec;
	LoadSiteList( siteVec, ptrMainOptions->ArgMap["SiteList"] );
	cout << "site list: " << ptrMainOptions->ArgMap["SiteList"] << " loaded!" << endl;

// load ref seq
	vector<RefSeq*> REF_SEQ;
	InitializeMeiSeqRef( REF_SEQ, ptrMainOptions->ArgMap["MElist"] );
	
	vector<string> subinfo;
	subinfo.resize(3);
	subinfo[0] = ptrMainOptions->ArgMap["Sample"];
	subinfo[1] = ptrMainOptions->ArgMap["Bam"];
	subinfo[2] = ptrMainOptions->ArgMap["DiscoverDir"];

// re-genotype single sample
	cout << "Re-genotype sample " << subinfo[0] << " at chr: " << ptrMainOptions->ArgMap["Chr"] << ", mei-type = " << ptrMainOptions->ArgMap["MeiType"] << "..." << endl;
	string sm = subinfo[0] + "." + ptrMainOptions->ArgMap["MeiType"];
	ReGenotypeSingleVcf( REF_SEQ, siteVec, subinfo, ptrMainOptions->ArgMap["rgDir"], ptrMainOptions->ArgMap["MeiType"], ptrMainOptions->ArgMap["Chr"] );
	GenerateDoneFile( ptrMainOptions->ArgMap["rgDir"], sm.c_str() );
	string out_vcf_name = ptrMainOptions->ArgMap["rgDir"] + "refined-" + subinfo[0] + string("-") + ptrMainOptions->ArgMap["Chr"] + "." + ptrMainOptions->ArgMap["MeiType"] + ".vcf";
	cout << "Re-genotype finished with no error reported. Check re-genotype vcf at: " << out_vcf_name << endl;
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
		if ( current_info.size() != 3 ) {
			cerr << "ERROR: [LoadSampleList()] current line:\n " << line << "is not a regular sample list line!\n sample list should be 3 fields: sample-name bam discover-out-directory." << endl;
			exit(1);
		}
		SampleList.push_back( current_info );
	}
	sample_list.close();
// add slash to discover directory
	for( vector<vector<string> >::iterator ptr = SampleList.begin(); ptr != SampleList.end(); ptr++ ) {
		int last = (*ptr)[2].length() - 1;
		if ( (*ptr)[2][last] != '/' )
			(*ptr)[2] += '/';
	}
}

// binary options
void SetGenotypeGlobalOptions( Options * ptrMainOptions )
{
// set parallel
	if (ptrMainOptions->OptMap["Parallel"])
		PARALLEL = 1;

// set debug mode
	if (ptrMainOptions->OptMap["debug"])
		DEBUG_MODE = 1;

// set single end
	if ( ptrMainOptions->OptMap["includeSingleAnchor"] )
		SINGLE_SIDE = 1;
		
// pseodu-chr
	if ( ptrMainOptions->OptMap["pseudoChr"] )
		PSEUDO_CHR = 1;
		
// non-variant
	if ( ptrMainOptions->OptMap["printNonVariant"] )
		PRINT_NON_VARIANT = 1;

// dp filter
	if ( ptrMainOptions->OptMap["disableDPfilter"] )
		APPLY_DEPTH_FILTER = 0;
		
// no ref allele base?		
	if ( ptrMainOptions->OptMap["noRefAllele"] )
		REF_ALLELE = 0;

// no break point refine?
	if ( ptrMainOptions->OptMap["noBreakPoint"] )
		REFINE_BREAK_POINT = 0;
}

// set parameters used in both read map & GL calculation
void SetGenotypeGlobalParameters( Options * ptrMainOptions )
{
	WIN = stoi(ptrMainOptions->ArgMap["Win"]);
	STEP = stoi(ptrMainOptions->ArgMap["Step"]);
	REF_CHR = ptrMainOptions->ArgMap["CtrlChr"];
	LEVEL = 1;
//	NON_OFFSET = LEVEL + 1;
}

// set bam-related parameters
void SetGenotypeReadMapGlobals( string & qinfo_name )
{
	ifstream qinfo;
	qinfo.open( qinfo_name.c_str() );
	string line;
	vector<string> vc;
	while( getline( qinfo, line ) ) {
		std::stringstream ss;
		ss << line;
		string item;
		getline( ss, item, '\t');
		getline( ss, item, '\t');
		vc.push_back( item );
	}
	if ( vc.size() != 3 ) {
		cerr << "ERROR: [SetGenotypeReadMapGlobals] " << qinfo_name << " doesn't have 3 lines. Please check!" << endl;
		exit(1);
	}
// avr read length & minimum read length
	int avr_read_len = stoi( vc[0] );
//	cout << "    Average length = " << avr_read_len << " bp." << endl;
	SetRMGReadLength( avr_read_len );

// avr ins size & minimum ins size	
	int avr_ins_size = stoi( vc[1] );
//	cout << "    Average insert size = " << avr_ins_size << " bp." << endl;
	SetRMGinsSize( avr_ins_size );
	
//	cout << "  Estimating bam depth: " << endl;
	float dp = stof ( vc[2] );
//	cout << "    Rough depth = " << dp << "x." << endl;
	SetRMGdepth( dp );
	
// close
	qinfo.close();	
}


// print re-genotype command to file
void PrintParallelCommand( ofstream & paraFile, vector<string> & subsample, string & site_list_name, Options* ptrMainOptions, string & mt, string & current_chr, string & rgdir )
{
	paraFile << MPATH << "bin/Morphling reGenotype -Win " << WIN << " -Step " << STEP << " -CtrlChr " << REF_CHR << " -MElist " << ptrMainOptions->ArgMap["MElist"];
	paraFile << " -DiscoverDir " << subsample[2] << " -Bam " << subsample[1] << " -rgDir " << rgdir << " -MeiType " << mt << " -Chr " << current_chr;
	paraFile << " -SiteList " << site_list_name << " -Sample " << subsample[0] << endl;
}


// re-genotype a single sample, must specify chr & mei-type
void ReGenotypeSingleVcf( vector<RefSeq*> REF_SEQ, vector<int> & siteVec, vector<string> subinfo, string & rg_dir, string & mt, string & current_chr )
{
	string sample = subinfo[0];
	string sbam = subinfo[1];
	string discover_dir = subinfo[2];

// use a dummy to make sure this event won't mapped to other type of MEI
	RefSeq * Dummy_ref_seq = new RefSeq;
	vector<RefSeq*> single_REF_SEQ;
	single_REF_SEQ.resize(3, NULL);
	int mei_index = std::stoi( mt );
	for( int rs = 0; rs <= 2; rs++ ) {
		if ( rs == mei_index )
			single_REF_SEQ[ rs ] = REF_SEQ[ mei_index ];
		else
			single_REF_SEQ[ rs ] = Dummy_ref_seq;
	}
	
// set numeric threshold according to qinfo generated by Morphling-discover
	string qinfo_name = discover_dir + "QC/QC.info";
	SetGenotypeReadMapGlobals( qinfo_name );
	
// open bam & ref stat
	string refs_prefix = discover_dir + "QC/sref." + mt;
	RefStats rstats( refs_prefix, mei_index );
	SamFile samIn;
	SamFileHeader samHeader;
	OpenBamAndBai( samIn, samHeader, sbam );
	bool section_status = samIn.SetReadSection( current_chr.c_str() );
	if ( !section_status ) {
		cerr << "ERROR: [ReGenotypeSingleVcf] Unable to set read section at: " << current_chr << ", in file: " << sbam << endl;
		exit(1);
	}
	
// re-genotype
	string mt_suffix = string(".") + mt + ".vcf";
	string vcf_suffix = string("-") + current_chr + mt_suffix;
	string in_vcf_name = discover_dir + "split/Hits" + vcf_suffix;
	string out_vcf_name = rg_dir + "refined-" + sample + vcf_suffix;
	implementSingleVcf( siteVec, single_REF_SEQ, rstats, samIn, samHeader, in_vcf_name, out_vcf_name);
	
// clear
	delete Dummy_ref_seq;	
}

// get GL for single sample based on site list, then print out to vcf
// do it by chr
void implementSingleVcf( vector<int> & siteVec, vector<RefSeq*> & REF_SEQ, RefStats & rstats, SamFile & samIn, SamFileHeader & samHeader, string & in_vcf_name, string & out_vcf_name )
{
	int MaxKeepDist = WIN / 3; // max dist to use the GL instead of getting a new one
	ifstream in_vcf;
	in_vcf.open( in_vcf_name.c_str() );
	CheckInputFileStatus( in_vcf, in_vcf_name.c_str() );
	ofstream out_vcf;
	out_vcf.open( out_vcf_name );
	CheckOutFileStatus( out_vcf, out_vcf_name.c_str() );
	
	string line;
	string chr;
	vector<int>::iterator site_ptr = siteVec.begin();
	while( getline( in_vcf, line ) ) {
		VcfRecord vcf_rec;
		vcf_rec.SetFromLine( line );
		int position = vcf_rec.GetPosition();
		int dist = *site_ptr - position; // ref - current_read_position
		while ( dist < -MaxKeepDist ) { // We've passed some ref positions without genotyping them
			VcfRecord implemented_rec;
			if ( chr.empty() )
				chr = vcf_rec.GetChromosome();
			ResetVcfRecordFromBam( implemented_rec, rstats, REF_SEQ, chr, *site_ptr, samIn, samHeader );
// add the break point estimated from single sample as an Info field
			if ( implemented_rec.GetDosage() > 0 && (implemented_rec.GetPosition() - *site_ptr != 0) )
				implemented_rec.AddIntegerInfoField( "Breakp", implemented_rec.GetPosition() );
//then set the position from site list as position
			implemented_rec.SetPosition( *site_ptr );
			implemented_rec.PrintRecord( out_vcf );
			site_ptr++;
			if ( site_ptr == siteVec.end() ) { // no more site to be re-genotyped
				out_vcf.close();
				in_vcf.close();
				return;
			}
			dist = *site_ptr - position;
		}
  // now current vcf is either a candidate or not candidate before the next candidate
  		while ( dist <= MaxKeepDist ) {
  		// set original position as breakp
  			if ( vcf_rec.GetDosage() > 0 && dist != 0)
  				vcf_rec.AddIntegerInfoField( "Breakp", vcf_rec.GetPosition() );
  		// then set cross-files vcf position
  			vcf_rec.SetPosition( *site_ptr );
  			vcf_rec.PrintRecord( out_vcf );
  			site_ptr++;
  			if ( site_ptr == siteVec.end() ) { // no more site to be re-genotyped
				out_vcf.close();
				in_vcf.close();
				return;
			}
			dist = *site_ptr - position;
  		}
  	// should do a continue for dist > MaxKeepDist. But no code below this. not necessary.
	}
// implement candidate sites after last vcf record
	if ( chr.empty() ) {
		cerr << "ERROR: [ReGenotypeSingleVcf()] No vcf record at " << in_vcf_name << endl;
		exit(1);
	}
	while( site_ptr != siteVec.end() ) {
		VcfRecord implemented_rec;
		ResetVcfRecordFromBam( implemented_rec, rstats, REF_SEQ, chr, *site_ptr, samIn, samHeader );
		if ( implemented_rec.GetDosage() > 0 && (implemented_rec.GetPosition() - *site_ptr != 0) )
			implemented_rec.AddIntegerInfoField( "Breakp", implemented_rec.GetPosition() );
//then set the position from site list as position
		implemented_rec.SetPosition( *site_ptr );
		implemented_rec.PrintRecord( out_vcf );
		site_ptr++;
	}
	out_vcf.close();
}










































