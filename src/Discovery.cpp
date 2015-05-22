#include "Discovery.h"
#include "Utilities.h"
#include "RefStats.h"

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h> // opendir

using std::string;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

// main function for doing discovery
void DiscoverMeiHits( const char* proper_prefix, const char* disc_prefix, const char* focus_chr,
	const char* ctrl_proper_prefix, const char* ctrl_disc_prefix, const char* out_prefix )
{
	string str_focus_chr = string(focus_chr);
  // let's first check if focus_chr data exist if not -1

  // loop through each mei_type	
	for(int mei_type = 0; mei_type <= 2; mei_type++) {
		RefStats * rsPtr = new RefStats( ctrl_proper_prefix, ctrl_disc_prefix, mei_type );
		OriginalStats * dataOsPtr = new OriginalStats();
		if ( focus_chr.compare("-1") == 0 ) { // whole genome: use opendir to get a list of all chrs
			addToDataViaDirectory( dataOsPtr, proper_prefix, disc_prefix );
		}
		else { // only discover on focus_chr
			string proper_chr_prefix = string(proper_prefix) + '-' + focus_chr;
			string disc_chr_prefix = string(disc_prefix) + '-' + focus_chr;
			dataOsPtr->Add( proper_chr_prefix.c_str(), disc_chr_prefix.c_str() );
		}
		dataOsPtr->ReOrganize();
		
	// calculate GL for each cell
		for( vector< MergeCell >::iterator merge_it = dataOsPtr->MergeData.begin(); merge_it != dataOsPtr->MergeData.end(); merge_it++ ) {
			rsPtr->SetRecordGL( merge_it->counts, merge_it->GL );
		}
		
		string vcf_name = string(out_prefix) + "-" + mei_type + ".vcf";
		PrintMeiHitsAsVcf( dataOsPtr->MergeData, vcf_name.c_str(), mei_type );
		
		delete rsPtr;
		delete dataOsPtr;
	}
}

// print a single type MEI vcf out
void PrintMeiHitsAsVcf( vector< MergeCell > & MergeData, const char* vcf_name, int mei_type )
{
	string mei_str = GetMeiTypeStr( mei_type );
	ofstream outVcf;
	outVcf.open( vcf_name );
	CheckOutFileStatus( outVcf, vcf_name );
	
	time_t raw_time;
	time(&raw_time);
	
	outVcf << "##fileformat=VCFv4.1" << endl;
	outVcf << "##filedate=" << ctime(&raw_time) << endl;
	outVcf << "##source=LHMEI" << endl;
	
	outVcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Most Likely Genotype\">" << endl;
	outVcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
	outVcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Call Quality\">" << endl;
	outVcf << "##FORMAT=<ID=PL,Number=.,Type=Integer,Description=\"Genotype Likelihoods for Genotypes in Phred Scale, for 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, ...\">" << endl;
	outVcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << SampleName << endl;

	int half_win = WIN / 2;
	for( map< string, vector< MergeCellPtr > >::iterator map_it = GenomeLocationMap.begin(); map_it != GenomeLocationMap.end(); map_it++ ) {
		for( vector< MergeCellPtr >::iterator data_ptr = map_it->second.begin(); data_ptr != map_it->second.end(); data_ptr++ ) {
			int location = data_ptr->win_index * STEP + half_win;
			outVcf << map_it->first << "\t" << location << ".\t.\t<" << MeiType << ">\t99\tPASS\tGT:DP:GQ:PL\t";
			int depth = GetVecSum( (*data_ptr)->ptr->counts );
			int gt_quality = GetMeiCountSum( (*data_ptr)->counts );
			Varint * variant = new Varint;
			SetGTandPL( variant, (data_ptr)->ptr->GL );
			outVcf << variant->genotype << ":" << depth << ":" << gt_quality << ":" << variant->pl[0] << "," << variant->pll[1] << "," << variant->pl[2] << endl;
		}
	}
	
	outVcf.close();
}


/********** inner function *************/
// add all proper & disc in their directory to data
// use in DiscoverMeiHits when focus_chr == -1
// assume proper & disc are in the same dir
void addToDataViaDirectory( OriginalStats * dataOsPtr, const char* proper_prefix, const char* disc_prefix )
{
	string proper_prefix_str = string(proper_prefix);
	int dir_end = proper_prefix_str.length();
	dir_end--;
	for( ; dir_end >= 0; dir_end-- ) {
		if ( proper_prefix[ dir_end ] == '/' )
			break;
	}
	string dir_name = dir_end > 0 ? proper_prefix_str.substr(0, dir_end) : string('.');
	struct dirent *pDirent;
	DIR *pDir;
	pDir = opendir( dir_name.c_str() );
  // sanity check
	if ( pDir == NULL ) {
		cerr << "ERROR: Cannot open directory: " << dir_name << endl;
		exit(1);
	}
  // loop through: get chr name list
  	vector< string > chr_name_list;
  	while( (pDirent = readdir(pDir)) != NULL ) {
  		string str_dname = string(pDirent->d_name);
  	// only read .3 since disc do not have .3
  		if ( pDirent->d_name[ str_dname.length() - 1 ] != '3' )
  			continue;
  		int dash_st = str_dname.length() - 3;
  		for( ;dash_st != dir_end; dash_st-- )
  			if ( str_dname[dash_st] == '-' )
  				break;
  		if ( dash_st == dir_end ) {
  			cerr << "ERROR: Can't find dash in proper name: " << str_dname << endl;
  			exit(1);
  		}
  		string chr_name = str_dname.substr( dash_st + 1, str_dname.length() - dasn_st - 4 );
  		chr_name_list.push_back( chr_name );
  	}
  	closedir( pDir );
  	
// add all to dataOs
	for( vector< string >::iterator it = chr_name_list.begin(); it != chr_name_list.end(); it++ ) {	
  		string proper_chr_prefix = string(proper_prefix) + '-' + *it;
		string disc_chr_prefix = string(disc_prefix) + '-' + *it;
		dataOsPtr->Add( proper_chr_prefix.c_str(), disc_chr_prefix.c_str() );
  	}
}



/**** Utility funcitons ******/

int GetVecSum( vector<int> & counts )
{
	int sum = 0;
	for( vector<int>::iterator it = counts.begin(); it != counts.end(); it++ )
		sum += (*it);
	return sum;
}

int GetMeiCountSum( vector<int> & counts )
{
	vector<int> as_mei = {4,5,8,9,12,13,16,17};
	int sum = 0;
	for( vector<int>::iterator it = asMei.begin(); it != asMei.end(); it++ )
		sum += counts[ *it ];
	return sum;
} 

// normalize GL to PL, and set genotype
void SetGTandPL( Varint * variant, vector<int> & GL )
{
	int type;
	if ( GL[0] >= GL[1] ) {
		type = GL[2] > GL[0] ? 2 : 0;
	}
	else { // 0 < 1
		type = GL[2] > GL[1] ? 2 : 1;
	}
// set genotype
	if (type == 0)
		varint->genotype = string("0/0");
	else if (type == 1)
		varint->genotype = string("0/1");
	else // 2
		variant->genotype = string("1/1");
// set pl
	for(int i=0; i<=2; i++) {
		if ( i == type ) { // largest
			variant->pl[i] = 0;
		}
		else { // others
			variant->pl[i] = GL[i] - GL[type];
		}
	}
}












