#include <string> 
#include <dirent.h>
#include <iostream>
#include <utility>

#include "ComputeLHMEI.h"
#include "Utilities.h" // file check , done file generate etc
#include "Preprocess.h" // generate ref-b am & disc_bam
#include "ReadMap.h"
#include "bamUtility.h" // getAvrReadLen;
#include "Globals.h" // DEBUG_MODE, WIN, STEP
#include "SetReadMapGlobals.h"
#include "RefStats.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

void ComputeLHMEI (Options * ptrMainOptions)
{
	SetGlobalOptions( ptrMainOptions );
	SetGlobalParameters( ptrMainOptions );
	SetReadMapGlobals( ptrMainOptions );

	string work_dir = ptrMainOptions->ArgMap["WorkDir"];
	if (work_dir[ work_dir.size() - 1 ] != '/')
		work_dir += '/';
	string bam_dir = work_dir + "bam_tmps/";
	string ctrl_dir = work_dir + "ctrl_tmps/";
	string pre_dir = work_dir + "preprocess/";
	string split_dir = work_dir + "split/";
	string vcf_name = split_dir + "Hits";

// prepare directories	
	string cmd = "mkdir -p " + bam_dir + " " + ctrl_dir + " " + pre_dir + " " + split_dir;
	ExecuteCmd(cmd);	
	
/*** pre-process: do not focus on any chr ***/
	string outSam = pre_dir + REF_CHR + ".sam";
	string discNsortBam = pre_dir + "disc-nsort";
	string discSam = pre_dir + "disc.sam";
	if (!ExistDoneFile( pre_dir, "PreProcess" )) {
		PreProcessBam( ptrMainOptions->ArgMap["Bam"], outSam, discSam, ptrMainOptions->ArgMap["MEcoord"]);
		string disc_nsort_cmd = string("samtools sort -n ") + discSam + " " + discNsortBam;
		ExecuteCmd(disc_nsort_cmd);
		GenerateDoneFile( pre_dir, "PreProcess" );
	}
	discNsortBam += ".bam";

/*** Generate level list if >= 0 ***/
  // generate bam index of discSam first
	string level_file_prefix = pre_dir + "level-list";
	if (LEVEL >= 0) {
		if ( !ExistDoneFile( pre_dir, "LevelList" ) ) {
			GenerateLevelListFromDiscBam( discNsortBam, level_file_prefix );
			GenerateDoneFile( pre_dir, "LevelList" );	
		}
	}

/*** re-map ***/
	string ctrl_bam = ctrl_dir + REF_CHR + "-remap-sort-recal.bam";
	if ( !ExistDoneFile(ctrl_dir, "Remap") ) {
	  // nsort
		string cmd = string("samtools sort -n ") + outSam + " " + ctrl_dir + REF_CHR + "-nsort";
		ExecuteCmd(cmd);
	  // fastq
	  	cmd = string("bam bam2FastQ --in ") + ctrl_dir + REF_CHR + "-nsort.bam --readName --outBase " + ctrl_dir + REF_CHR;
	  	ExecuteCmd(cmd);
	  // re-map
	  	string fastq_prefix = ctrl_dir + REF_CHR;
	  	string remapSam = ctrl_dir + "align-pe.sam";
	  	cmd = GetRemapCmd(ptrMainOptions->ArgMap["Mapper"], fastq_prefix, ptrMainOptions->ArgMap["SliceFA"], remapSam);
	  	ExecuteCmd(cmd);
	  // sort by coord
	  	cmd = string("samtools sort ") + remapSam + " " + ctrl_dir + REF_CHR + "-remap-sort";
	  	ExecuteCmd(cmd);
	  // dedup
	  	cmd = string("bam dedup --recab --in ") + ctrl_dir + REF_CHR + "-remap-sort.bam --out ";
	  	cmd += ctrl_bam + " --force --refFile " + ptrMainOptions->ArgMap["SliceFA"] + " --storeQualTag OQ --maxBaseQual 40";
	  	ExecuteCmd(cmd);
	  // generate .bai
	  	cmd = string("samtools index ") + ctrl_bam;
	  	ExecuteCmd(cmd);
	  // generate done
	  	GenerateDoneFile( ctrl_dir, "Remap" );
	}

	string focus_chr_str = ptrMainOptions->ArgMap["Chr"];
/*** Read Map (target + ctrl )****/
	string bam_done_flag = string("BamMap-");
	bool bam_counted;
	if ( focus_chr_str.compare("-1")  == 0 ) {
		bam_done_flag += string("all");
		bam_counted = ExistDoneFile( bam_dir, bam_done_flag.c_str() );
	}
	else {
		string all_done_flag = bam_done_flag; // if all done, then consider single-chr is done
		all_done_flag += string("all");
		bam_done_flag += focus_chr_str;
		bam_counted = ( ExistDoneFile( bam_dir, bam_done_flag.c_str() ) | ExistDoneFile( bam_dir, all_done_flag.c_str() ) );
	}
	bool ctrl_counted = ExistDoneFile( ctrl_dir, "CtrlBamMap" );
	if ( !bam_counted || !ctrl_counted ) {
		ReadMap * Rmap = new ReadMap( ptrMainOptions->ArgMap["MElist"] );
		string focus_chr = focus_chr_str.compare("-1") == 0 ? string("") : focus_chr_str;
		if ( !bam_counted ) {
			cout << "Setting read types from raw bam..." << endl;
			Rmap->SetMapFromBam( ptrMainOptions->ArgMap["Bam"], discNsortBam, bam_dir, focus_chr, level_file_prefix );
			GenerateDoneFile( bam_dir, bam_done_flag.c_str() );
		}
		if ( !ctrl_counted ) {
			cout << "Setting read types from ctrl bam..." << endl;
			string ctrl_coord = ptrMainOptions->ArgMap["MEcoord"] + ".liftOver";
			Rmap->SetMapFromCtrlBam( ctrl_bam, ctrl_dir, ctrl_coord );
			GenerateDoneFile( ctrl_dir, "CtrlBamMap" );
		}
		delete Rmap;
	}
	
/*** Calculate LH ****/
// here do not set done file
	AllHetIndex allHetIndex;
	SetAllHetIndex( ptrMainOptions->ArgMap["HetIndex"], allHetIndex );
	string sample_name = ptrMainOptions->ArgMap["Sample"];
	if ( sample_name.compare(".") == 0 ) { // sample name not specified. Use base name of -Bam with no .bam prefix
		string base_name = GetFileBaseName( ptrMainOptions->ArgMap["Bam"] );
		string suffix = base_name.substr( base_name.length() - 4 );
		if ( suffix.compare(".bam") == 0 && base_name.length() > 4 )
			base_name = base_name.substr(0, base_name.length() - 4);
		sample_name = base_name;
	}
// do by mei type
	int focus_type = stoi(ptrMainOptions->ArgMap["MeiType"]);
	for( int mei_type = 0; mei_type <= 2; mei_type++ ) {
		if ( focus_type != -1 && mei_type != focus_type ) {
			cout << "Skip mei type " << mei_type << "!" << endl;
			continue;
		}
		cout << "Discovering mei-type: " << mei_type << " ..." << endl;
		string lh_flag = string("LH-") + focus_chr_str + "." + std::to_string(mei_type);
		if ( !ExistDoneFile( bam_dir, lh_flag.c_str() ) ) {	
			string ctrl_proper_prefix = ctrl_dir + "proper-" + REF_CHR;
			string ctrl_disc_prefix = ctrl_dir + "disc-" + REF_CHR;
		// construct ref stats
			RefStats* rStats = new RefStats( ctrl_proper_prefix, ctrl_disc_prefix, mei_type, allHetIndex );
			rStats->SetCtrlGLs();
			string outRecord = work_dir + "refLH." + std::to_string(mei_type) + ".report";
		// ref LH
			if ( !( ptrMainOptions->OptMap["noCtrlVcf"] ) ) {
				string ref_lh_done = string("RefLH-") + std::to_string( mei_type );
				if ( ExistDoneFile( work_dir, ref_lh_done.c_str() ) )
					cerr << "Warning: exist " << ref_lh_done << ". Skip RefLH part!" << endl;
				else {
					rStats->PrintCtrlGLasRecord( outRecord, ctrl_bam, ptrMainOptions->ArgMap["SliceFA"] );
					GenerateDoneFile( work_dir, ref_lh_done.c_str() );	
				}
			}
			if ( ptrMainOptions->OptMap["printRefStats"] ) { // debug: print refStats
				string refPrefix = ptrMainOptions->ArgMap["refPrefix"] + "." + std::to_string(mei_type);
				rStats->PrintRefStats( refPrefix );
			}
			rStats->MarkRefLHasDone();
			rStats->ReAdjustSelfsWithLiftOver();
		// data LH: loop-through bam dir ---> re-organize --> 
			OriginalStats* dataOsPtr = new OriginalStats( mei_type, sample_name );
			
		  //add to memory
		  	vector<string> chr_list; // also used in LH if needed
			if ( focus_chr_str.compare("-1") != 0 ) { // single chr
				cout << "Discovering single chr: " << focus_chr_str << " ..." << endl;
				string data_proper_name = bam_dir + "proper-" + focus_chr_str;
				string data_disc_name = bam_dir + "disc-" + focus_chr_str;
				bool add_success = dataOsPtr->Add( focus_chr_str, data_proper_name, data_disc_name );
				if ( !add_success ) {
					cerr << "Warning: no available proper reads in chromesome " << focus_chr_str << ", skipped this chr on mei_type = " << mei_type << "!" << endl;
					continue;
				}
			}
			else { // get chr list from bam header. Then add to memory. Skip pseudo chr
				SetChrListFromBamHeader( chr_list, ptrMainOptions->ArgMap["Bam"] );
				string out_str;
				for( vector<string>::iterator current_chr = chr_list.begin(); current_chr != chr_list.end(); current_chr++ ) {			
					if ( current_chr->length() > 5 && (!PSEUDO_CHR)) // skip pseudo chr
						continue;
					out_str += (" " + *current_chr);
				}
				cout << "Discovering whole genome, chr: " << out_str << " ..." << endl;;
				for( vector<string>::iterator current_chr = chr_list.begin(); current_chr != chr_list.end(); current_chr++ ) {			
					if ( current_chr->length() > 5 && (!PSEUDO_CHR)) // skip pseudo chr
						continue;
					string proper_prefix = bam_dir + "proper-" + *current_chr;
					string disc_prefix = bam_dir + "disc-" + *current_chr;
					bool single_success = dataOsPtr->Add( *current_chr, proper_prefix, disc_prefix );
					if (!single_success)
						cerr << "Warning: no available proper reads in chr: " << *current_chr << ", skipped this on mei_type = " << mei_type << endl;
				}
			}
			
		// re-organize & clear under level
			dataOsPtr->ReOrganize();
			dataOsPtr->ClearUnderLevelMergeCells(); // to speed up
		// calculate LH
			for( MergeCellPtr merge_it = dataOsPtr->MergeData.begin(); merge_it != dataOsPtr->MergeData.end(); merge_it++ ) {
			  	rStats->SetRecordGL( merge_it );
			}
			rStats->AdjustUsedLoci( dataOsPtr );
			if ( focus_chr_str.compare("-1") != 0 ) { // single chr
				string split_vcf_name = vcf_name + "-" + focus_chr_str + "." + std::to_string(mei_type) + ".vcf";
				dataOsPtr->PrintGLasVcf( split_vcf_name, ptrMainOptions->ArgMap["Bam"], ptrMainOptions->ArgMap["GenomeFasta"], focus_chr_str );
			}
			else { // whole genome
				for( vector<string>::iterator current_chr = chr_list.begin(); current_chr != chr_list.end(); current_chr++ ) {
					if ( current_chr->length() > 5 && (!PSEUDO_CHR)) // skip pseudo chr
						continue;
					string split_vcf_name = vcf_name + "-" + (*current_chr) + "." + std::to_string(mei_type) + ".vcf";
					dataOsPtr->PrintGLasVcf( split_vcf_name, ptrMainOptions->ArgMap["Bam"], ptrMainOptions->ArgMap["GenomeFasta"], *current_chr );
				}
			}
		}
	}
// generate a whole vcf from all split vcf
	cout << "Merging & generating final vcf..." << endl;
	string final_name = work_dir + sample_name + "-final.vcf";
	ofstream final_vcf;
	final_vcf.open( final_name.c_str() );
	CheckOutFileStatus( final_vcf, final_name.c_str() );
	final_vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << endl;
	final_vcf.close();
	string final_cmd = "cat " + vcf_name + "* | awk '$6>=10' | grep PASS | sort -k1,1 -k2,2n >> " + final_name; 
	ExecuteCmd( final_cmd );
	
	cout << "Morphling Discover finished with no error reported. Check final output at: " << final_name << endl;
/* delete intermediates:
	ctrl dir: *.fastq, *-nsort.bam, *-remap-sort.bam, align-pe.sam*
	preprocess: *.sam
*/
}

 // set bool global parameters in Globals.h from OptMap in Options
void SetGlobalOptions( Options * ptrMainOptions )
{
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
void SetGlobalParameters( Options * ptrMainOptions )
{
	WIN = stoi(ptrMainOptions->ArgMap["Win"]);
	STEP = stoi(ptrMainOptions->ArgMap["Step"]);
	REF_CHR = ptrMainOptions->ArgMap["CtrlChr"];
	LEVEL = stoi( ptrMainOptions->ArgMap["Simplify"] );
	NON_OFFSET = LEVEL + 1;
}

// set bam-related parameters
void SetReadMapGlobals( Options * ptrMainOptions )
{
// avr read length & minimum read length
	int avr_read_len = stoi( ptrMainOptions->ArgMap["ReadLen"] );
	if ( avr_read_len < 0 ) {
		cout << "Estimating read length from bam: " << endl;
		avr_read_len = GetAvrReadLenFromBam( ptrMainOptions->ArgMap["Bam"] );
		cout << "    Average length = " << avr_read_len << " bp." << endl;
	}
	SetRMGReadLength( avr_read_len );

// avr ins size & minimum ins size	
	int avr_ins_size = stoi( ptrMainOptions->ArgMap["InsSize"] );
	if ( avr_ins_size < 0 ) {
		cout << "Estimating insert size from bam: " << endl;
		avr_ins_size = GetAvrInsSizeFromBam( ptrMainOptions->ArgMap["Bam"] );
		cout << "    Average insert size = " << avr_ins_size << " bp." << endl;
	}
	SetRMGinsSize( avr_ins_size );
	
	float dp = stof( ptrMainOptions->ArgMap["Depth"] );
	if ( dp < 0 ) {
		cout << "Estimating bam depth: " << endl;
		dp = EstimateBamDepth( ptrMainOptions->ArgMap["Bam"], avr_read_len, REF_CHR );
		cout << "    Rough depth = " << dp << "x." << endl;
	}
	SetRMGdepth( dp );
}







