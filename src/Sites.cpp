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
Sites::Sites( string & vcf_name, vector<string> & PreAsb, string & out_vcf_name, string & me_list_name )
{
// open files & set members
	InVcfName = vcf_name;
	OutVcf.open( out_vcf_name.c_str() );
	CheckOutFileStatus( OutVcf, out_vcf_name.c_str() );
	initializeSiteInfo();
	NpolyA.resize( Nsite, 0 );
	NpolyT.resize( Nsite, 0 );
	loadMEsequence( me_list_name, MEseqs, MEnames );
	loadPreAsb( PreAsb );	
	
// clean and load other seq	
	keepMaxSubtype();
	loadPreSeq( PreAsb );
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
			if ( line[0] != '#' )
				past_header = 1;
			else
				continue;
		}
		stringstream ss;
		ss << line;
		string field;
		getline( ss, field, '\t' );
		getline( ss, field, '\t' );
		getline( ss, field, '\t' );
		getline( ss, field, '\t' );
		getline( ss, field, '\t' );
		int mei = GetMEtypeFromAlt( field );
		MeiType.push_back(mei);
		site_size++;
	}
	in_vcf.close();
	if ( site_size == 0 ) {
		cerr << "ERROR: [initializeSiteInfo] no available site in vcf: " << InVcfName << ". No need to do assembly!" << endl;
		exit(1);
	}
	Nsite = site_size;
}


void Sites::loadPreAsb( vector<string> & PreAsb )
{
// initialize first
	Clusters.resize( Nsite );
	for( int i=0; i<Nsite; i++ )
		Clusters[i].resize(4);

// sample by sample
	int sp = 0;
	for( vector<string>::iterator p = PreAsb.begin(); p != PreAsb.end(); p++, sp++ )
		loadSinglePreAsb( *p, sp );
}

void Sites::loadSinglePreAsb( string & pre_name, int sp )
{
	ifstream pre;
	pre.open( pre_name.c_str() );
	CheckInputFileStatus( pre, pre_name.c_str() );
	string line;
	int site = -1;
	while( getline( pre, line ) ) {
		site++;
		if ( site >= Nsite ) {
			cerr << "ERROR: [loadSinglePreAsb] " << pre_name << " has more lines than #sites!" << ", pre name = " << pre_name << endl;
			exit(1);
		}
		if ( line.empty() ) // skip sample with no read info
			continue;
		int mei_index = line[0] - '0';
		if ( mei_index < 0 || mei_index > 2) {
			cerr << "ERROR: [loadSinglePreAsb] mei_index = " << mei_index << " at " << pre_name << ", line " << site+1 << ", pre name = " << pre_name << endl;
			exit(1);
		}
// add polyA & polyT
		int aend, tend;
		for( int i=2; i<(int)line.size(); i++) {
			if ( line[i] == ':' ) {
				aend = i;
				break;
			}
		}
		if ( aend >= (int)line.size() - 1 ) {
			cerr << "ERROR: [loadSinglePreAsb] no poly t record at line: " << site + 1 << ", pre_name=" << pre_name << endl;
			exit(1);
		}
		for( int i=aend+1; i<(int)line.size(); i++) {
			if ( line[i] == ':' ) {
				tend = i;
				break;
			}
		}
		string nb = line.substr(2, aend-2);
		if ( !std::all_of( nb.begin(), nb.end(), isdigit) ) {
			cerr << "ERROR: [loadSinglePreAsb] polyA field contains non-digit number " << nb << " at line: " << site+1 << ", pre_name=" << pre_name << endl;
			exit(1);
		}
		NpolyA[site] += stoi(nb);
		nb = line.substr(aend+1, tend-aend-1);
		if ( !std::all_of( nb.begin(), nb.end(), isdigit) ) {
			cerr << "ERROR: [loadSinglePreAsb] polyT field contains non-digit number " << nb << " at line: " << site+1 << ", pre_name=" << pre_name << endl;
			exit(1);
		}
		NpolyT[site] += stoi(nb);
		
		stringstream ss;
		ss << line.substr( tend+1 ); // remove mei_index:
		string subf;
		int cl = 0;
		while( getline( ss, subf, ':' ) ) {
			if ( cl > 3 ) {
				cerr << "ERROR: [loadSinglePreAsb] line has more than 4 fields in sample " << pre_name << ", line " << site+1 << ", pre name = " << pre_name << endl;
				exit(1);
			}
			if ( subf.find_first_not_of(';') == std::string::npos ) { // all ';', no read info in this cluster
				if ( subf.size() != MEnames[mei_index].size() - 1 ) {
					cerr << "ERROR: [loadSinglePreAsb] Empty line doesn't have ==subtype fields at line " << site+1 << ", pre name " << pre_name << endl;
					exit(1);
				}
				cl++;
				continue;
			}
			stringstream ss14;
			ss14 << subf;
			if ( Clusters[site][cl].empty() )
				Clusters[site][cl].resize( MEseqs[mei_index].size() );
			vector< subCluster >::iterator sub = Clusters[site][cl].begin();
			string field14;
			while( getline( ss14, field14, ';' ) ) {
				if ( sub == Clusters[site][cl].end() ) {
					cerr << "ERROR: [loadSinglePreAsb] #subfields > #subtype at " << pre_name << ", line " << site+1 << ", pre name = " << pre_name << endl;
					exit(1);
				}
				string fig;
				int idx = 0; // idx == 0 is map key
				int key;
				int evi[5];
				stringstream flss;
				flss << field14;
				while( getline( flss, fig, ',' ) ) {
					if ( !std::all_of( fig.begin(), fig.end(), isdigit ) ) {
						cerr << "ERROR: [loadSinglePreAsb] line contain non-digit chars in sample " << pre_name << ", str = " << fig << ", line " << site+1 << ", pre name = " << pre_name << endl;
						exit(1);
					}
					if ( idx == 0 ) {
						key = stoi(fig);
						idx++;
					}
					else if ( idx == 5 ) {
						evi[idx-1] = stoi(fig);
						EviInfo new_evi;
						new_evi.Boundary = evi[0];
						new_evi.LAlign = evi[1];
						new_evi.RAlign = evi[2];
						new_evi.Score = evi[3];
						new_evi.SeqKey = evi[4];
						new_evi.SampleKey = sp;
						(*sub)[key].push_back( new_evi );
						idx = 0;
					}
					else {
						evi[idx-1] = stoi(fig);
						idx++;	
					}
				}
				if ( idx != 0 ) {
					cerr << "ERROR: [loadSinglePreAsb] line doesn't have 5x fields in sample " << pre_name << ", line " << site+1 << ", field = " << field14 << ", pre name = " << pre_name << endl;
					exit(1);
				}
				sub++;
			}
			if (subf[ subf.size() - 1 ] == ';' ) // last field empty
				sub++;
			if ( sub != Clusters[site][cl].end() ) {
				cerr << "ERROR: [loadSinglePreAsb] line " << site+1 << ", cluster " << cl << " does not contain all subtype info! dist = " << Clusters[site][cl].end() - sub << ", pre name = "  << pre_name << endl;
				exit(1);
			}
			cl++;
		}
		if ( cl != 4 ) {
			cerr << "ERROR: [loadSinglePreAsb] line " << site+1 << " does not contain all 4 cluster info!" << endl;
			exit(1);
		}
	}
	pre.close();
	if ( site < Nsite-1 ) {
		cerr << "ERROR: [loadSinglePreAsb] " << pre_name << " has fewer lines than #sites!" << endl;
		exit(1);
	}
}

void Sites::keepMaxSubtype()
{
	Subtypes.resize( Nsite );
	Strands.resize( Nsite );
	for( int i=0; i<Nsite; i++ ) {
		setSinlgeSiteSubtypeAndStrand( Clusters[i], i );
	}
}


void Sites::setSinlgeSiteSubtypeAndStrand( vector< vector< subCluster > > & sclust, int sp )
{
	int subsize = MEnames[ MeiType[sp] ].size();
// calculate scores
	vector< vector<int> > SumScoreVec;
	SumScoreVec.resize(4);
	for( int cl=0; cl<4; cl++ ) {
		SumScoreVec[cl].resize( subsize, -1 );
		int sub_index = 0;
		for( vector< subCluster >::iterator psub = sclust[cl].begin(); psub != sclust[cl].end(); psub++, sub_index++ ) { // subtype
			if ( psub->empty() ) {
				SumScoreVec[cl][sub_index] = -1;
				continue;
			}
			// sum score
			int score_sum = 0;
			for( subCluster::iterator pread = psub->begin(); pread != psub->end(); pread++ ) { // each read
				for( int i=0; i<(int)pread->second.size(); i++ )
					score_sum += pread->second[i].Score;
			}
			SumScoreVec[cl][sub_index] = score_sum;
		}
	}

	int subtype;
	bool strand;	
// find most likely cluster by 2-end score
	vector<int> plus_strand;
	vector<int> minus_strand;
	plus_strand.resize( (int)SumScoreVec[0].size() );
	minus_strand.resize( (int)SumScoreVec[0].size() );
	for( int i=0; i<(int)SumScoreVec[0].size(); i++ ) {
		plus_strand[i] = SumScoreVec[0][i] + SumScoreVec[2][i];
		minus_strand[i] = SumScoreVec[1][i] + SumScoreVec[3][i];
	}
	vector<int>::iterator pmax_plus = std::max_element( plus_strand.begin(), plus_strand.end() );
	vector<int>::iterator pmax_minus = std::max_element( minus_strand.begin(), minus_strand.end() );
	if ( *pmax_plus == *pmax_minus ) {
		if ( *pmax_plus <= 0 ) {
			Subtypes[sp] = -1;
			return;
		}
	}
// set strand
	bool use_plus;
	if ( *pmax_plus > *pmax_minus )
		use_plus = 1;
	else if ( *pmax_plus < *pmax_minus )
		use_plus = 0;
	else {
		if ( NpolyA[sp] > NpolyT[sp] )
			use_plus = 1;
		else if ( NpolyA[sp] < NpolyT[sp] )
			use_plus = 0;
		else {
			cerr << "Warning: at site: " << sp << " plus==minus. sw=" << *pmax_plus << ", polyA=" << NpolyA[sp] << ". Use '+' strand!" << endl;
			use_plus = 1;
		}
	}
// clear other
	if ( use_plus ) {
		subtype = pmax_plus - plus_strand.begin();
		strand = 1;
		sclust[1].clear();
		sclust[3].clear();
		clearOtherSubcluster( sclust[0], subtype );
		clearOtherSubcluster( sclust[2], subtype );
	}
	else {
		subtype = pmax_minus - minus_strand.begin();
		strand = 0;
		sclust[0].clear();
		sclust[2].clear();
		clearOtherSubcluster( sclust[1], subtype );
		clearOtherSubcluster( sclust[3], subtype );
	}

// set
	Subtypes[sp] = subtype;
	Strands[sp] = strand;
}


void Sites::clearOtherSubcluster( vector< subCluster > & sc, int sub )
{
	for( int i=0; i<(int)sc.size(); i++ ) {
		if ( i == sub ) // skip the one kept
			continue;
		sc[i].clear();
	}
}


void Sites::loadPreSeq( vector<string> & PreAsb )
{
// initialize first
	Seqs.resize( Nsite );
	for( int i=0; i<Nsite; i++ )
		Seqs[i].resize(NSAMPLE);

// sample by sample
	int sp = 0;
	for( vector<string>::iterator p = PreAsb.begin(); p != PreAsb.end(); p++,sp++ ) {
		string seq_info_name = *p + ".seq";
		ifstream seq_info;
		seq_info.open( seq_info_name.c_str() );
		CheckInputFileStatus( seq_info, seq_info_name.c_str() );
		string line;
		int site = 0;
		while( getline(seq_info, line) ) {
			if ( site >= Nsite ) {
				cerr << "ERROR: [loadPreSeq] " << seq_info_name << " #lines > #sites!" << endl;
				exit(1);
			}
			if ( !line.empty() ) {
				stringstream ss;
				ss << line;
				string seq;
				while( getline( ss, seq, ',' ) )
					Seqs[site][sp].push_back( seq );
			}
			site++;
		}
		seq_info.close();
		if ( site != Nsite ) {
			cerr << "ERROR: [loadPreSeq] site = " << site << ", " << *p << " does not have " << Nsite << " lines!" << endl;
			exit(1);
		}
	}
}


/** assembly related functions ***/

void Sites::AssemblySubtypes()
{
// open input vcf
	ifstream in_vcf;
	in_vcf.open( InVcfName.c_str() );
	CheckInputFileStatus( in_vcf, InVcfName.c_str() );

// open vcf
	string line;
	while( getline( in_vcf, line ) ) {
		if ( line[0] == '#' ) {
			if ( line[1] == 'C' )
				printAddedVcfHeaders();
			OutVcf << line << endl;
		}
		else
			break;
	}
	
// loop through each site	
	for( int site=0; site<Nsite; site++ ) {
		if ( site > 0 )
			getline( in_vcf, line );
		int subtype = Subtypes[site];
		int meitype = MeiType[site];
		bool strand = Strands[site];
		int c1,c2;
		if ( strand ) {
			c1 = 0;
			c2 = 2;
		}
		else {
			c1 = 1;
			c2 = 3;
		}
		if ( subtype == -1 )
			printUnAssembledRecord( line, site );
		else {
		// check if it's a one-side hit
			if ( Clusters[site][c1].empty() ) { // no left info
				subCluster dummy;
				if ( Clusters[site][c2].empty() )
					cerr << "Warning: [AssemblySubtypes] no read info for both end at site " << site << ". Skip this site!" << endl;
				else {
					AsbSite cs( Seqs[site], dummy, Clusters[site][c2][subtype], strand, MEseqs[meitype][subtype] );
					printSingleRecord( line, site, cs );
				}
			}
			else { // left info exists
				if ( Clusters[site][c2].empty() ) { // no right info
					subCluster dummy;
					AsbSite cs( Seqs[site], Clusters[site][c1][subtype], dummy, strand, MEseqs[meitype][subtype] );
					printSingleRecord( line, site, cs );
				}
				else { // both end exists
					AsbSite cs( Seqs[site], Clusters[site][c1][subtype], Clusters[site][c2][subtype], strand, MEseqs[meitype][subtype] );
					printSingleRecord( line, site, cs );
				}
			}
		}
	}
	
// clear
	in_vcf.close();
	OutVcf.close();	
}


void Sites::printAddedVcfHeaders()
{
	OutVcf << "##INFO=<ID=SUB,Number=1,Type=Char,Description=\"MEI subtype\">" << endl;
	OutVcf << "##INFO=<ID=AVRDP,Number=1,Type=Float,Description=\"Average ALT depth in each 1/* sample\">" << endl;
	OutVcf << "##INFO=<ID=MPOS,Number=2,Type=Integer,Description=\"SV start and end on MEI consensus sequence\">" << endl;
	OutVcf << "##INFO=<ID=MISSING,Number=1,Type=Integer,Description=\"#Missing bases MPOS\">" << endl;
	OutVcf << "##INFO=<ID=STRAND,Number=1,Type=Char,Description=\"SV strand\">" << endl;
	OutVcf << "##INFO=<ID=ASBSAMPLES,Number=1,Type=Integer,Description=\"Total samples used in assembly\">" << endl;
//	OutVcf << "##INFO=<ID=VASBSAMPLES,Number=1,Type=Integer,Description=\"Total samples used with valid MEI reads in assembly\">" << endl;
	OutVcf << "##INFO=<ID=AVRPOLYA,Number=1,Type=FLOAT,Description=\"Average PolyA/T base count per sample\">" << endl;
}


void Sites::printUnAssembledRecord( string & vline, int site )
{
	int info_end = GetTabLocation( 0, 8, vline );
	int n = NpolyA[site] >= NpolyT[site] ? NpolyA[site] : NpolyT[site];
	OutVcf << vline.substr(0, info_end) << ";SUB=NA;SVLEN=NA;SVCOV=NA;MISSING=NA;MPOS=NA,NA;STRAND=NA;ASBSAMPLES=0;AVRPOLYA=" << n << vline.substr(info_end) << endl;
}

void Sites::printSingleRecord( string & vline, int site, AsbSite & cs )
{
	int info_end = GetTabLocation( 0, 8, vline );

// should set filter later
	if ( !cs.IsAssembled() ) { // for no-assembly site	
		OutVcf << vline.substr(0, info_end) << ";SUB=NA;SVLEN=NA;SVCOV=NA;MISSING=NA;MPOS=NA,NA;STRAND=NA;ASBSAMPLES=0" << vline.substr(info_end) << endl;
	}
	else { // do print
		OutVcf << vline.substr(0, info_end) << ";SUB=" << MEnames[ MeiType[site] ][ Subtypes[site] ]<< ";SVLEN=" << cs.GetSVlength();
		float svd = cs.GetSVdepth();
		if ( svd >= 0 )
			OutVcf << ";SVCOV=" << std::setprecision(2) << std::fixed << cs.GetSVdepth();
		else
			OutVcf << ";SVCOV=NA";
		OutVcf << ";MISSING=" << cs.GetMissingBaseCount() << ";MPOS=" << cs.GetLeftMost() << "," << cs.GetRightMost();
		OutVcf << ";STRAND=";
		if( Strands[site] )
			OutVcf << "+";
		else
			OutVcf << "-";
		int ns = cs.GetSampleCount();
		OutVcf << ";ASBSAMPLES=" << ns;
//		OutVcf << ";VASBSAMPLES=" << cs.GetValidSampleCount();
		if ( ns==0 )
			ns++;
		if ( Strands[site] )
			OutVcf << ";AVRPOLYA=" << std::setprecision(1) << (float)NpolyA[site] / ns;
		else
			OutVcf << ";AVRPOLYA=" << std::setprecision(1) << (float)NpolyT[site] / ns;
		OutVcf << vline.substr(info_end) << endl;
	}
}








