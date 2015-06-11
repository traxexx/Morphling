#include "Sites.h"
#include "Utilities.h"
#include "Globals.h"
#include <fstream>
#include <iostream>
#include <algorithm> // all_of
#include <map>
#include <ctype.h> // isalpha/digit
#include <sstream>
#include <iomanip>  // std::setprecision

using std::stringstream;
using std::map;
using std::cout;
using std::cerr;
using std::endl;

// constructor: prepare for asembly
Sites::Sites( string & vcf_name, string & sample_list_name, string & out_vcf_name, string & me_list_name )
{
	InVcfName = vcf_name;
	// open files
	OutVcf.open( out_vcf_name.c_str() );
	CheckOutFileStatus( OutVcf, out_vcf_name.c_str() );
	SampleNames.resize( NSAMPLE);
// initialize	
	initializeSiteInfo();
	BamFileNames.resize(NSAMPLE);
	BamFiles.resize(NSAMPLE);
	BamFileHeaders.resize(NSAMPLE);
// vcf
	setSiteAndSamples();
// sample list
	setBamFilesFromSampleList( sample_list_name );
// consensus sequence
	loadMEsequence( me_list_name );
}

void Sites::initializeSiteInfo()
{
// count site number first
	string line;
	int site_size = 0;
	bool past_header = 0;
	ifstream in_vcf;
	in_vcf.open( InVcfName.c_str() );
	CheckInputFileStatus( in_vcf, InVcfName.c_str() );
	while( getline( in_vcf, line ) ) {
		if ( !past_header ) {
			if ( line[0] != '#' ) {
				past_header = 1;
				site_size++;
			}
			else {
				OutVcf << line << endl;
			}
		}
		else
			site_size++;
	}
	in_vcf.close();
	if ( site_size == 0 ) {
		cerr << "ERROR: [initializeSiteInfo] no available site in vcf: " << InVcfName << ". No need to do assembly!" << endl;
		exit(1);
	}
	else
		SiteInfo.resize( site_size );
}

void Sites::setSiteAndSamples()
{
// print added info to out vcf
	printAddedVcfHeaders();

// set site & gt list
	ifstream in_vcf;
	in_vcf.open( InVcfName.c_str() );
	CheckInputFileStatus( in_vcf, InVcfName.c_str() );
	string line;
	string last_line;
	bool past_header = 0;
	vector<PotSite>::iterator si_ptr = SiteInfo.begin();
	while( getline( in_vcf, line ) ) {
	// set header
		if ( !past_header ) {
			if ( line[0] == '#' ) {
				last_line = line;
				continue;
			}
			else { // use last line to add sample names
				stringstream ss;
				ss << last_line;
				string field;
				for( int i=0; i<9; i++ )
					getline( ss, field, '\t' );
				int idx = 0;
				while( getline( ss, field, '\t' ) ) {
					if ( idx >= NSAMPLE ) {
						cerr << "ERROR: vcf has more samples than sample list!" << endl;
						exit(1);
					}
					SampleNames[idx] = field;
					idx++;
				}
				last_line.clear();
				past_header = 1;
			}
		}
	// set site: pr(rr) <= 0.5
		stringstream ss;
		ss << line;
		string field;
		for( int i=0; i<9; i++ ) {
			getline( ss, field, '\t' );
			if ( i==1 )
				si_ptr->Position = stoi(field);	
			else if ( i==4 )
				si_ptr->MEtype = GetMEtypeFromAlt( field );
		}
		si_ptr->GtList = new bool[ NSAMPLE ];
		bool* bpt = si_ptr->GtList;
		while( getline( ss, field, '\t' ) ) {
			stringstream flss;
			flss << field;
			string subf;
			getline( flss, subf, ':' );
			getline( flss, subf, ':' );
			getline( flss, subf, ':' );
			if ( stoi(subf) > 3 )
				*bpt = 1;
			else
				*bpt = 0;
			bpt++;
		}
		si_ptr++;
	}
	in_vcf.close();

// print sample list line
	OutVcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for( int i=0; i<NSAMPLE; i++ )
		OutVcf << "\t" << SampleNames[i];
	OutVcf << endl;
}


void Sites::setBamFilesFromSampleList( string & sample_list_name )
{
// make a sample name map first
	map< string, string > nameMap;
	ifstream in_list;
	in_list.open( sample_list_name.c_str() );
	CheckInputFileStatus( in_list, sample_list_name.c_str() );
	string line;
	while( getline( in_list, line ) ) {
		string sample;
		stringstream ss;
		ss << line;
		getline( ss, sample, '\t' );
		string bam_name;
		getline( ss, bam_name, '\t' );
		getline( ss, bam_name, '\t' );
		nameMap[ sample ] = bam_name;
	}
	in_list.close();
	
// add to bam file name by the order of sample names
	for( int i=0; i<NSAMPLE; i++ ) {
		map< string, string >::iterator np = nameMap.find( SampleNames[i] );
		if ( np == nameMap.end() ) {
			cerr << "ERROR: [setBamFilesFromSampleList] vcf column " << SampleNames[i] << " does not exist in sample list!" << endl;
			exit(1);
		}
		BamFileNames[i] = np->second;
		nameMap.erase( SampleNames[i] );
	}
	if ( !nameMap.empty() ) {
		cerr << "ERROR: [setBamFilesFromSampleList] these samples do not exist in vcf: " << endl;
		cerr << "    ";
		for( map< string, string >::iterator it = nameMap.begin(); it != nameMap.end(); it++ )
			cerr << " " << it->first;
		cerr << endl;
		exit(1);
	}
	
// open bam files
	for( int i=0; i< NSAMPLE; i++ ) {
		BamFiles[i].OpenForRead( BamFileNames[i].c_str() );
		if ( !BamFiles[i].IsOpen() ) {
			cerr << "ERROR: Can't open " << BamFileNames[i] << endl;
			exit(1);
		}
		bool header_status = BamFiles[i].ReadHeader( BamFileHeaders[i] );
		if ( !header_status ) {
			cerr << "ERROR: Can't open header " << BamFileNames[i] << endl;
			exit(1);
		}
		string bai = BamFileNames[i] + ".bai";
		bool bai_status = BamFiles[i].ReadBamIndex( bai.c_str() );
		if ( !bai_status ) {
			cerr << "ERROR: Can't read bam index " << bai << endl;
			exit(1);
		}
	}
}

void Sites::printAddedVcfHeaders()
{
	OutVcf << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Estimated SV length\">" << endl;
	OutVcf << "##INFO=<ID=AVRDP,Number=1,Type=Float,Description=\"Average ALT depth in each 1/* sample\">" << endl;
	OutVcf << "##INFO=<ID=MPOS,Number=2,Type=Integer,Description=\"SV start and end on MEI consensus sequence\">" << endl;
	OutVcf << "##INFO=<ID=MISSING,Number=1,Type=Integer,Description=\"#Missing bases in assembled SV\">" << endl;
	OutVcf << "##INFO=<ID=STRAND,Number=1,Type=Char,Description=\"SV strand\">" << endl;
	OutVcf << "##INFO=<ID=ASBSAMPLES,Number=1,Type=Integer,Description=\"Total samples used in assembly\">" << endl;
}

void Sites::loadMEsequence( string & me_list_name )
{
// read me list
	ifstream me_list;
	me_list.open( me_list_name.c_str() );
	CheckInputFileStatus( me_list, me_list_name.c_str() );
	string line;
	vector<string> me_files;
	me_files.resize(3);
	while( getline( me_list, line ) ) {
		stringstream ss;
		ss << line;
		string me_name;
		getline( ss, me_name, '\t' );
		if ( me_name.compare("ALU") == 0 )
			getline(ss, me_files[0], '\t');
		else if ( me_name.compare("L1") == 0 )
			getline(ss, me_files[1], '\t');
		else if ( me_name.compare("SVA") == 0 )
			getline(ss, me_files[2], '\t');
		else {
			cerr << "ERROR: invalid ME name in MElist: " << me_name << endl;
			exit(1);
		}
	}
	me_list.close();
	if ( me_files[0].empty() || me_files[1].empty() || me_files[2].empty() ) {
		cerr << "ERROR: one or more ME fasta missing!" << endl;
		exit(1);
	}
// read me seqs & names
	MEseqs.resize(3);
	MEnames.resize(3);	
// set size first
	for( int i=0; i<3; i++ ) {
		ifstream infa;
		infa.open( me_files[i].c_str() );
		CheckInputFileStatus( infa, me_files[i].c_str() );
		string line;
		int rec_count = 0;
		while( getline( infa, line ) ) {
			if ( line[0] == '>' )
				rec_count++;
		}
		infa.close();
		if ( rec_count == 0 ) {
			cerr << "ERROR: fasta " << me_files[i] << " does not contain any fasta record!" << endl;
			exit(1);
		}
		MEseqs[i].resize( rec_count );
		MEnames[i].resize( rec_count);
	}
// load name & seq
	for( int mindex = 0; mindex<3; mindex++ ) {
		ifstream infa;
		infa.open( me_files[mindex].c_str() );
		CheckInputFileStatus( infa, me_files[mindex].c_str() );
		string line;
		int ri = 0;
		string seq; // sequence
		string name; // name
		while( getline( infa, line ) ) {
			if ( line.size() == 0 ) // skip empty lines
				continue;
			if ( line[0] == '>' ) { // name line
				if ( !seq.empty() ) {
					MEseqs[mindex][ri] = seq;
					MEnames[mindex][ri] = name;
					ri++;
				}
				seq.clear();
				int space_index = line.length();
				for( int i=1; i<(int)line.length(); i++ ) {
					if ( line[i] == ' ' || line[i] == '\t' ) {
						space_index = i;
						break;
					}
				}
				name = line.substr( 1, space_index - 1 );
			}
			else { // seq line
				seq += line;
			}
		}
		MEseqs[mindex][ri] = seq;
		MEnames[mindex][ri] = name;
		infa.close();
	}	
}


void Sites::AssemblySubtypes()
{
// open input vcf
	ifstream in_vcf;
	in_vcf.open( InVcfName.c_str() );
	CheckInputFileStatus( in_vcf, InVcfName.c_str() );
// loo through each site
	string line;
	while( getline( in_vcf, line ) ) {
		if ( line[0] != '#' )
			break;
	}
	bool first_line = 1;
	for( vector<PotSite>::iterator pinfo = SiteInfo.begin(); pinfo != SiteInfo.end(); pinfo++ ) {
	// check if position match
		if ( !first_line )
			getline( in_vcf, line );
		else
			first_line = 0;
		int pstart = GetTabLocation( 0, 1, line );
		int pend = GetTabLocation( pstart+1, 1, line );
		string pos_str = line.substr( pstart + 1, pend - pstart - 1 );
		if ( !std::all_of( pos_str.begin(), pos_str.end(), isdigit ) ) {
			cerr << "ERROR: POS contain non-digit chars at: " << line.substr(0, pend) << endl;
			exit(1);
		}
		if ( stoi( pos_str ) != pinfo->Position ) {
			cerr << "ERROR: Position does not match! vcf = " << line.substr(0, pend) << ", while cs = " << pinfo->Position << endl;
			exit(1);
		}
	// do assembly
		string chr = line.substr(0, pstart);
		AsbSite currentSite( chr, pinfo->Position, pinfo->GtList, BamFiles, BamFileHeaders, MEseqs[pinfo->MEtype] );
		currentSite.Assembly();
	// print
		printSingleRecord( pend, pinfo->MEtype, line, currentSite );
	}
	in_vcf.close();
	OutVcf.close();
}


void Sites::printSingleRecord( int pend, int mei_index, string & vline, AsbSite & cs )
{
	int alt_start = GetTabLocation( pend+1, 2, vline );
	int alt_end = GetTabLocation( alt_start+1, 1, vline);
	int info_end = GetTabLocation( alt_end+1, 3, vline );

// should set filter later
	if ( !cs.IsAssembled() ) { // for no-assembly site	
		OutVcf << vline.substr(0, info_end) << ";SVLEN=NA;SVCOV=NA;MISSING=NA;MPOS=NA,NA;STRAND=NA;ASBSAMPLES=0" << vline.substr(info_end) << endl;
	}
	else { // do print
		if ( cs.GetSubtype() >= (int)MEnames[mei_index].size() ) {
			cerr << "ERROR: at " << cs.GetPosition() << ", subtype = " << cs.GetSubtype() << ", but MEnames size = " << MEnames[mei_index].size() << endl;
			exit(1);
		}
		OutVcf << vline.substr(0, alt_start) << "\t<INS:ME:" << MEnames[mei_index][ cs.GetSubtype() ] << ">" << vline.substr(alt_end, info_end - alt_end) << ";SVLEN=" << cs.GetSVlength();
		OutVcf << ";SVCOV=" << std::setprecision(2) << std::fixed << cs.GetSVdepth();
		OutVcf << ";MISSING=" << cs.GetMissingBaseCount() << ";MPOS=" << cs.GetLeftMost() << "," << cs.GetRightMost();
		OutVcf << ";STRAND=";
		if( cs.GetStrand() )
			OutVcf << "+";
		else
			OutVcf << "-";
		OutVcf << ";ASBSAMPLES=" << cs.GetSampleCount();
		OutVcf << vline.substr(info_end) << endl;
	}
}


int GetMEtypeFromAlt( string & field)
{
	int index = field.length() - 2;
	int mtype;
	if ( field[index] == 'U' )
		mtype = 0;
	else if ( field[index] == '1' )
		mtype = 1;
	else if ( field[index] == 'A' )
		mtype =2;
	else {
		cerr << "ERROR: [GetMEtypeFromAlt] " << field << " is not a regular ALT from Morphling Genotype!" << endl;
		exit(1);
	}
	return mtype;
}


// search from search_start, end at 0-based position of the noccur(th)
int GetTabLocation( int search_start, int noccur, string & line )
{
	int n = 0;
	int loc = -1;
	for( int i=search_start; i<(int)line.length(); i++ ) {
		if ( line[i] == '\t' ) {
			n++;
			if ( n == noccur ) {
				loc = i;
				break;
			}	
		}
	}
	if ( loc < 0 ) {
		cerr << "ERROR: [GetTabLocation] Unable to find " << noccur << "th tab in " << line << endl;
		exit(1);
	}
	return loc;
}



