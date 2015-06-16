#include "Assembly.h"
#include "Sites.h"
#include "Globals.h"
#include "Utilities.h"
#include "bamUtility.h"
#include <iostream>
#include <utility>
#include <fstream>
#include <algorithm>    // std::count
#include "MeiSeqs.h" // load ME sequence
#include "SamFile.h" // read bam
#include "QC.h" // passQC
#include "Cluster.h" // pre-assemble class
#include <sys/stat.h> // stat() check file exists
#include <sstream>

using std::cout;
using std::endl;
using std::cerr;
using std::stringstream;

// pre-assemble
//		do re-mapping
//		output all score & align position on each MEI subtype
//		out format:
//			mei:key,i,i,..i,key,i....;...;... (4 semi-colon fields)
void PreAssemble( Options* ptrMainOptions )
{
// globals
	SetPreAsbGlobals( ptrMainOptions );

	vector<string> ChrNames;
	vector<int> PotLocation; // position
	vector<int> MeiType;
	vector<bool> Valid; // if GQ < 3, do not count read types
	
// load sites
	ifstream in_vcf;
	in_vcf.open( ptrMainOptions->ArgMap["Vcf"].c_str() );
	CheckInputFileStatus( in_vcf, ptrMainOptions->ArgMap["Vcf"].c_str() );
	string line;
	string sample_name= ptrMainOptions->ArgMap["Sample"];
	int index = 8;
	while( getline( in_vcf, line ) ) {
		if ( line.empty() ) {
			cerr << "Warning: skipped empty line in vcf." << endl;
			continue;
		}
		if ( line[0] == '#' ) {
			if ( line[1] == 'C' ) {
				stringstream ss;
				ss << line;
				string field;
				for( int i=0; i<9; i++ )
					getline( ss, field, '\t' );
				while( getline( ss, field, '\t' ) ) {
					index++;
					if ( sample_name.compare(field) == 0 )
						break;
				}
			}
			continue;
		}
		string field;
		stringstream ss;
		ss << line;
		getline( ss, field, '\t' );
		ChrNames.push_back( field );
		getline( ss, field, '\t' );
		if ( !std::all_of( field.begin(), field.end(), isdigit ) ) {
			cerr << "ERROR: vcf position field is not digit: " << field << endl;
			exit(1);
		}
		PotLocation.push_back( stoi(field) );
		getline( ss, field, '\t' );
		getline( ss, field, '\t' );
		getline( ss, field, '\t' );
		MeiType.push_back( GetMEtypeFromAlt(field) );
		int gt_start = GetTabLocation( 0, index, line );
		int gt_end = GetTabLocation( index + 1, 1, line );
		field = line.substr( gt_start+1, gt_end - gt_start - 1 );
		stringstream flss;
		flss << field;
		string subf;
		getline( flss, subf, ':' );
		getline( flss, subf, ':' );
		getline( flss, subf, ':' );
		if ( std::all_of( subf.begin(), subf.end(), isdigit ) ) {
			if ( stoi(subf) > 3 )
				Valid.push_back(1);
			else
				Valid.push_back(0);
		}
		else
			Valid.push_back(0);
	}
	if ( index < 9 ) {
		cerr << "ERROR: [PreAssemble] Can't find sample column: " << sample_name << " in vcf" << endl;
		exit(1);
	}
	in_vcf.close();

// load MEseq
	vector< vector<string> > MEseqs;
	vector< vector<string> > MEnames;
	loadMEsequence( ptrMainOptions->ArgMap["MElist"], MEseqs, MEnames );
		
	
// get info from bam
	SamFile bamIn;
	SamFileHeader bamHeader;
	OpenBamAndBai( bamIn, bamHeader, ptrMainOptions->ArgMap["Bam"] );
	ofstream out_info;
	ofstream seq_info;
	out_info.open( ptrMainOptions->ArgMap["Out"].c_str() );
	CheckOutFileStatus( out_info, ptrMainOptions->ArgMap["Out"].c_str() );
	string seq_info_name = ptrMainOptions->ArgMap["Out"] + ".seq";
	seq_info.open( seq_info_name.c_str() );
	CheckOutFileStatus( seq_info, seq_info_name.c_str() );

	for( int sp=0; sp<(int)PotLocation.size(); sp++ )	 {
		if ( !Valid[sp] ) { // skip invalid
			out_info << endl;
			seq_info << endl;
			continue;
		}
	// section
		int cluster_start = PotLocation[sp] - WIN/2;
		int cluster_end = PotLocation[sp] + WIN / 2;
		bool section_status = bamIn.SetReadSection( ChrNames[sp].c_str(), cluster_start, cluster_end );
		if ( !section_status )
			cerr << "Warning: [PreAssemble] unable set section: chr" << ChrNames[sp] << " " << cluster_start << "-" << cluster_end << ". Set as no pre-asseble info" << endl;
	// set read info
		Cluster scluster( MeiType[sp], MEseqs[ MeiType[sp] ] );
	  // add clip, record disc
		vector<DiscInfo> discVec;
		SamRecord sam_rec;
		while( bamIn.ReadRecord( bamHeader, sam_rec ) ) {
			if ( !PassQC(sam_rec) )
				continue;
			if ( sam_rec.getFlag() & 0x2 ) // proper
				scluster.AddProper( sam_rec );
			else { // save disc info
				if ( !AssemblyDiscPass(sam_rec) )
					continue;
				if ( sam_rec.getFlag() & 0x4 ) // do not use unmap as anchor
					continue;
				DiscInfo di;
				di.Chr = sam_rec.getMateReferenceName();
				di.Position = sam_rec.get1BasedMatePosition();
				di.MatePosition = sam_rec.get1BasedPosition();
				discVec.push_back( di );
			}
		}
	  // add from disc info
		for( vector<DiscInfo>::iterator pd = discVec.begin(); pd != discVec.end(); pd++ ) {
			bool section_status = bamIn.SetReadSection( pd->Chr.c_str(), pd->Position, pd->Position + 1 );
			if ( !section_status )
				cerr << "Warning: [add2ClusterFromSingleBam] Unable to set read section. Skip this DiscInfo record!" << endl;
//			bool got_mate = 0;
			while( bamIn.ReadRecord( bamHeader, sam_rec ) ) {
				if ( !PassQC(sam_rec) )
					continue;
				if ( !AssemblyDiscPass(sam_rec) )
					continue;
				if ( pd->Position == sam_rec.get1BasedPosition() &&
					pd->MatePosition == sam_rec.get1BasedMatePosition() &&
					ChrNames[sp].compare( sam_rec.getMateReferenceName() ) == 0 &&
					pd->Chr.compare( sam_rec.getReferenceName() ) == 0 ) {
					scluster.AddDisc( sam_rec );
//					got_mate = 1;
						break;
				}
			}
		}
		scluster.Print( out_info, seq_info );			
	}
	bamIn.Close();
	out_info.close();
	seq_info.close();
}

void Assembly( Options* ptrMainOptions )
{
// set globals
	SetAsbGlobals( ptrMainOptions );
	
	vector<string> PreAsbs;
	LoadAssemblySampleList( ptrMainOptions->ArgMap["Vcf"], ptrMainOptions->ArgMap["SampleList"], PreAsbs );
	
	CheckPreAssembles( PreAsbs );

// load sites & bam list
	Sites allSites( ptrMainOptions->ArgMap["Vcf"], PreAsbs, ptrMainOptions->ArgMap["Out"], ptrMainOptions->ArgMap["MElist"] );

// cluster reads in each sites & print
	allSites.AssemblySubtypes();

// report finish
	cout << "Morphling Assembly finished with no error reported. Check final vcf at: " << ptrMainOptions->ArgMap["Out"] << endl;
}


// load pre sample list
// check:
//	 whether these pres exists
//	 whether sample columns all match those in vcf

void LoadAssemblySampleList( string & vcf_name, string & sample_list_name, vector<string> & pres )
{
// get name from vcf
	ifstream vcf;
	vcf.open( vcf_name.c_str() );
	CheckInputFileStatus( vcf, vcf_name.c_str() );
	string line;
	map<string, bool> namemap;
	while( getline( vcf, line ) ) {
		if ( line.empty() )
			continue;
		if( line[0] == '#' ) {
			if ( line[1] == 'C' ) {
				stringstream ss;
				ss << line;
				string field;
				for( int i=0; i<9; i++ )
					getline(ss, field, '\t');
				while( getline(ss, field, '\t') )
					namemap[field];
			}
		}
		else
			break;
	}
	vcf.close();

// compare those from sample list
	ifstream slist;
	slist.open( sample_list_name.c_str() );
	CheckInputFileStatus( slist, sample_list_name.c_str() );
	while( getline( slist, line) ) {
		stringstream ss;
		ss << line;
		string sample;
		getline( ss, sample, '\t' );
		string pre;
		getline( ss, pre, '\t' );
		map<string, bool>::iterator pf = namemap.find( sample );
		if ( pf == namemap.end() ) {
			cerr << "ERROR: [LoadAssemblySampleList] can't find sample " << sample << " in vcf!" << endl;
			exit(1);
		}
		pres.push_back( pre );
		namemap.erase(sample);
	}
	slist.close();
	
// check if namemap empty
	if ( !namemap.empty() ) {
		cerr << "ERROR: [LoadAssemblySampleList] these samples do not have pre-assembly info: " << endl;
		int i=0;
		for( map<string, bool>::iterator pf = namemap.begin(); pf != namemap.end(); pf++,i++ ) {
			if (i > 10)
				break;
			cerr << " " << pf->first;
		}
		if (i>10)
			cerr << "...";
		cerr << endl;
		exit(1);
	}
}


// global sub function
void SetAsbGlobals( Options* ptrMainOptions )
{
	std::ifstream in_list( ptrMainOptions->ArgMap["SampleList"].c_str() );
	NSAMPLE = std::count(std::istreambuf_iterator<char>(in_list), std::istreambuf_iterator<char>(), '\n');
	in_list.close();
	if ( NSAMPLE <= 0 ) {
		cerr << "ERROR: empty sample list: " << ptrMainOptions->ArgMap["SampleList"] << endl;
		exit(1);
	}
}

void SetPreAsbGlobals( Options* ptrMainOptions )
{
	WIN = stoi( ptrMainOptions->ArgMap["Win"] );
	if ( WIN <= 0 ) {
		cerr << "ERROR: win size = " << WIN << ", please specify a valida win size use -Win option." << endl;
		exit(1);
	}
}

void CheckPreAssembles( vector<string> & pres)
{
// check if all pres exist
	for( vector<string>::iterator p = pres.begin(); p != pres.end(); p++ ) {
		struct stat buffer;   
		if (stat (p->c_str(), &buffer) != 0) {
			cerr << "ERROR: pre-assembly " << *p << " does not exist!" << endl;
			exit(1);
		} 
	}
}


bool AssemblyDiscPass( SamRecord & sam_rec )
{
// check length
	if ( sam_rec.getReadLength() < 30 )
		return 0;
// check insert size
//	const char * equal = sam_rec.getMateReferenceNameOrEqual();
//	if ( !equal ) {
//		cerr << "Warning: Unable to get mate reference name at: " << sam_rec.getReadName() << endl;
//		return 0;
//	}
//	if ( *equal == '=' ) {
//		if ( sam_rec.getInsertSize() < 100 )
//			return 0;
//	}			
// indicate print to discSam
	return 1;
}



















