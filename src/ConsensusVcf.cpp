#include "ConsensusVcf.h"
#include "Utilities.h"
#include <iostream>
#include <sstream>
#include <iterator>
#include <time.h> // get file data
#include <GenomeSequence.h> // get ref allele

using std::to_string;
using std::cout;
using std::cerr;
using std::endl;
using std::stringstream;

ConsensusVcf::ConsensusVcf( int n_sample, int dist, string & chr_name ):
	chr( chr_name ),
	_nsample( n_sample ),
	_max_dist( dist )
{}

ConsensusVcf::~ConsensusVcf()
{}

void ConsensusVcf::SetSampleList( vector< vector<string> > sampleList )
{
	SampleNames.resize( _nsample );
	if ( _nsample != (int)sampleList.size() ) {
		cerr << "ERROR: total sample = " << _nsample << ", but sample list length = " << (int)sampleList.size() << "?" << endl;
		exit(1);
	}
	for( int i = 0; i < (int)sampleList.size(); i++ )
		SampleNames[i] = sampleList[i][0];
}

void ConsensusVcf::SetAltAlleleByMeiType( int mei_index )
{
	string alt = "<INS:ME:";
	if ( mei_index == 0 )
		alt += "ALU";
	else if ( mei_index == 1 )
		alt += "L1";
	else if ( mei_index == 2 )
		alt += "SVA";
	else {
		cerr << "ERROR: mei_index = " << mei_index << " does not match any MEI family1" << endl;
		exit(1);
	}
	alt += ">";
	this->alt_allele = alt;
}

void ConsensusVcf::SetFasta( string & fasta_name )
{
	_ref_fasta_name = fasta_name;
}

void ConsensusVcf::InitializeSdataFromSiteList( string & site_list_name )
{
	if ( !sdata.empty() )
		cerr << "Warning: sdata not empty. Replaced by new site list!" << endl;
		
	sdata.clear();
	ifstream site_list;
	site_list.open( site_list_name.c_str() );
	CheckInputFileStatus( site_list, site_list_name.c_str() );
	string line;
	while( getline( site_list, line) ) {
		stringstream ss;
		ss << line;
		string field;
		getline( ss, field, '\t' );
		int position = stoi( field );
		ConsensusVcfRecord * cv_ptr;
		cv_ptr = new ConsensusVcfRecord( position );
		sdata.push_back( cv_ptr );
	}
	site_list.close();
}

// don't need a valid header
// samples with same MEs must have same #lines
void ConsensusVcf::AddFromSingleVcf( int sample_index, string & vcf_name )
{
// read vcf & load to Data
	ifstream in_vcf;
	in_vcf.open( vcf_name.c_str() );
	CheckInputFileStatus( in_vcf, vcf_name.c_str() );
	string line;
	vector< ConsensusVcfRecord* >::iterator sptr = sdata.begin();
	while( getline( in_vcf, line ) ) {
		if ( sptr == sdata.end() ) {
			cerr << "ERROR: length of site list and sample: " << vcf_name << " does not match!" << endl;
			exit(1);
		}
// add to Data
		appendToCvcf( *sptr, sample_index, line );
		sptr++;
	}
// last sanity check
	if ( sptr != sdata.end() ) {
		cerr << "ERROR: length of site list > sample: " << vcf_name << " !" << endl;
		exit(1);	
	}
}

void ConsensusVcf::Polish()
{
	if ( this->_ref_fasta_name.empty() ) {
		cerr << "ERROR: in ConsensusVcf::Polish(),set ref fasta name before setting ref allele!" << endl;
		exit(1);
	}
	GenomeSequence* gs = new GenomeSequence( this->_ref_fasta_name );
	for( map<int, ConsensusVcfRecord* >::iterator mp = Data.begin(); mp != Data.end(); mp++ ) {
		mp->second->SetRefAllele( gs, this->chr );
		mp->second->EstimateAFviaEM();
		mp->second->EstimateFIC();
		mp->second->RefineInfoFields();
		mp->second->SetFilter();
	}
	delete gs;
}

// print Data as plain vcf
void ConsensusVcf::Print( string & out_name )
{
	ofstream out_vcf;
	out_vcf.open( out_name.c_str() );
	CheckOutFileStatus( out_vcf, out_name.c_str() );
	printFinalHeader( out_vcf ); // header + sample header
	for( map<int, ConsensusVcfRecord* >::iterator rec = Data.begin(); rec != Data.end(); rec++ ) {
		rec->second->Print( this->chr, out_vcf );
	}
	out_vcf.close();
}


void ConsensusVcf::appendToCvcf( ConsensusVcfRecord * cv, int & sp_index, string & line )
{
// parse line
		stringstream ss;
		ss << line;
		string field; // unnecessary fields
		getline( ss, field, '\t' ); // chr
		getline( ss, field, '\t' ); // position
		int location = stoi( field );
		getline( ss, field, '\t' ); // id
		getline( ss, field, '\t' );  // ref allele
		string alt_str;
		getline( ss, alt_str, '\t' );  // alt allele
		getline( ss, field, '\t' );  // qual
		getline( ss, field, '\t' );  // filter
		string info_str;
		getline( ss, info_str, '\t' );  // info
		getline( ss, field, '\t' );  // format
		string gl_str;
		getline( ss, gl_str, '\t' );  // gt & dp & gq & gl

	if( !cv->CheckPosition( location ) ) {
		cerr << "ERROR: site list & vcf does not match at: " << location << endl;
		exit(1);
	}
	cv->SetAltAllele( this->alt_allele );
	cv->SetInfoFields( sp_index, info_str );
	cv->SetGLFields( sp_index, gl_str );
}

// merge sdata to Data
void ConsensusVcf::MergeData()
{
// check if key exists in Data & add to Data
	for( vector< ConsensusVcfRecord* >::iterator sp = sdata.begin(); sp != sdata.end(); sp++ ) {
		int location = (*sp)->GetPosition();
		map<int, ConsensusVcfRecord* >::iterator data_ptr = Data.find( location );
		if ( data_ptr != Data.end() ) { // overlap key
			bool use_new = NewCvcfHigherRank( data_ptr->second, *sp );
			if ( use_new )
				data_ptr->second = *sp;
				
		}
		else { // new key
			Data[ location ] = *sp;
		}
	}

// clear sdata
	sdata.clear();
	
// check distance between keys
	map<int, ConsensusVcfRecord* >::iterator front_ptr = Data.begin();
	map<int, ConsensusVcfRecord* >::iterator rear_ptr = front_ptr++;
	if ( front_ptr == Data.end() ) { // only one record
		cerr << "Warning: only 1 vcf record exists in Data. Skip nearby merge!" << endl;
		return;
	}
	vector<int>  del_keys; // store keys to delete
	int front_index = 0;
	for( ; rear_ptr != Data.end(); rear_ptr++ ) {
		int dist = rear_ptr->first - front_ptr->first;
		if ( dist <= _max_dist ) {
			bool use_new = NewCvcfHigherRank( front_ptr->second, rear_ptr->second );
			if ( use_new )
				del_keys.push_back( front_index );
			else
				del_keys.push_back( front_index + 1 );
		}
		front_ptr++;
		front_index++;
	}
// delete keys
	for( vector<int>::iterator vp = del_keys.begin(); vp != del_keys.end(); vp++ ) {
		delete Data[*vp];
		Data.erase( *vp );
	}
}

// compare sum( variant_quality ) --> dosage sum
bool ConsensusVcf::NewCvcfHigherRank( ConsensusVcfRecord* old_cv, ConsensusVcfRecord* new_cv )
{
	int old_sum = old_cv->GetRawVariantQuality();
	int new_sum = new_cv->GetRawVariantQuality();
	if ( old_sum > new_sum )
		return 0;
	else if ( new_sum > old_sum )
		return 1;
		
	old_sum = old_cv->GetSumDosage();
	new_sum = new_cv->GetSumDosage();
	if ( old_sum > new_sum )
		return 0;
	else if ( new_sum > old_sum )
		return 1;
		
	cerr << "Warning: equal rank for position: " << old_cv->GetPosition() << " & " << new_cv->GetPosition() << endl;
	return 0;
}


void ConsensusVcf::printFinalHeader( ofstream & out_vcf )
{
	out_vcf << "##fileformat=VCFv4.1" << endl;
	time_t raw_time;
	time(&raw_time);
	out_vcf << "##fileDate=" << ctime(&raw_time) << endl;
	out_vcf << "##source=Morphling-v1.0" << endl;
	out_vcf << "##reference=file:"<< _ref_fasta_name << endl;
	out_vcf << "##contig=<ID=" << this->chr << ">" << endl;
	out_vcf << "##phasing=none" << endl;
	out_vcf << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << endl;
	out_vcf << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << endl;
	out_vcf << "##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">" << endl;
	out_vcf << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
	out_vcf << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
	out_vcf << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl;
	out_vcf << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl;
	out_vcf << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
	out_vcf << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Reference Allele Frequency\">" << endl;
	out_vcf << "##INFO=<ID=FIC,Number=A,Type=Float,Description=\"Genotype Likelihood based inbreeding Coefficient\">" << endl;
	out_vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	out_vcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
	out_vcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << endl;
	out_vcf << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Phread-scaled Genotype Likelihood Rounded To The Closest Integer\">" << endl;
	out_vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for( vector< string >::iterator ptr = SampleNames.begin(); ptr != SampleNames.end(); ptr++ )
		out_vcf << "\t" << *ptr;
	out_vcf << endl;
}


