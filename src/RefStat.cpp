#include "RefStat.h"
#include "morphError.h"
#include "genRef.h"
#include "bamUtility.h"
#include "Globals.h"
#include "qcCheck.h"
#include <fstream>
#include <iostream>
#include <math.h>

using std::endl;

/*

read region
extract reads from the region
remap

(class gen ref)
generate ref stats
print ref stats

*/


// the public one used to do all stuff
// need: remap output prefix
//		 ref coord bed prefix
//		 ME seq ref
void GenerateRefStats(
	SingleSampleInfo & msinfo, vector< vector<string> > & MEseq, vector< vector<string> > & MEname, 
	string & vcf_name, string & out_prefix, bool ref_exclusion )
{

// the input is (1~3):
//	chr (lift-overed)position subtype index

	string ref_dir = PATH + "refs/";
	string base_ctrl_bed = PATH + "refs/negative-base.bed";
	// get position from vcf if needed
	bool exclusion = ref_exclusion;
	string neg_ctrl_bed = out_prefix + ".neg.ctrl.bed";
	string exclude_bed = out_prefix + ".exclude.bed";
	if (exclusion) {
		string cmd = "cat " + vcf_name + " | cut -f1,2 | awk '$1==" + REF_CHR + "' | grep -v \"#\" | ";
		cmd += "awk '{print $1,$2-" + std::to_string(WIN/2) + ",$2+" + std::to_string(WIN/2);
		cmd += "}' | tr ' ' '\t' > " + exclude_bed;
		ExecuteCmd(cmd);
		if ( IsEmptyFile( exclude_bed ) ) // no record on REF_CHR. no exclusion
			exclusion = 0;
	}
	if (exclusion) { // generate negative-control bed by splitting base bed file and excluding the vcf positions
		int nex = GetFileLineCount( vcf_name );
		string split_out = out_prefix + ".split.bed";
		splitBaseBed( base_ctrl_bed, split_out, nex );
		string cmd = "bedtools intersect -a " + split_out + " -b " + exclude_bed + " -v > " + neg_ctrl_bed;
		ExecuteCmd(cmd);
	}
	else
		splitBaseBed( base_ctrl_bed, neg_ctrl_bed, 0);

	morphMessage("Starting remapping...");
	// remap
	map<int, int> ref_coord; // st -> end
	string unlift_name = PATH + "/refs/chr20-slice-MEI.bed";
	loadUnliftBed( ref_coord, unlift_name);
	string remap_bam_prefix = out_prefix + ".remap";
	makeRegionBam( ref_coord, msinfo, remap_bam_prefix );
	string remap_bam_name = remap_bam_prefix + ".bam";
	
	morphMessage("Remapping finished. Starting generating ctrl stat...");
	// generate all 3 ref
	SingleSampleInfo new_si = msinfo;
	new_si.bam_name = remap_bam_name;
	string ref_bed_prefix = ref_dir + "ref-MEI.";
	for(int m=0; m<NMEI; m++) {
		// hom
		string hom_out_name = out_prefix + ".hom-stat." + std::to_string(m);
		string bed_name = ref_bed_prefix + std::to_string(m);
		genRef gr( bed_name, new_si, MEseq[m], MEname[m], 0);
		gr.Print( hom_out_name );
		morphMessage( "Generated pos ctrl" );
		// neg
		genRef gr0( neg_ctrl_bed, msinfo, MEseq[m], MEname[m], 1);
		string neg_out_name = out_prefix + ".neg-stat." + std::to_string(m);
		gr0.Print( neg_out_name );
		morphMessage( "Generated neg ctrl" );
		// het
		vector< vector<float> > het_stat;
		gr.Export( het_stat );
		gr0.Export( het_stat );
		string het_out_name = out_prefix + ".het-stat." + std::to_string(m);
		PrintRefStats( het_stat, het_out_name );
		morphMessageNoTime( "Generated het ctrl" );
	}
	if (CLEAN) {
		string cmd = "rm -f " + out_prefix + ".exclude.bed " + out_prefix + ".neg.ctrl.bed ";
		cmd += remap_bam_name + " " + remap_bam_name + ".bai " + remap_bam_name + ".log ";
		cmd += remap_bam_name + ".qemp";
		ExecuteCmd(cmd);
	}
}

// print vector of vector float
// here empty element is allowed (didnot find ctrl mate)
void PrintRefStats( vector< vector<float> > & het_stat, string & out_name )
{
	std::ofstream outfile;
	outfile.open(out_name.c_str());
	for( int i=0; i<het_stat.size(); i++ ) {
		if (het_stat[i].empty()) // skip no-mate element
			continue;
		outfile << "1";
		if (het_stat[i].size() != NFIELD) {
			string str = "[PrintRefStats] at i=" + std::to_string(i) + ", het_stat doesn't have ";
			str += std::to_string(NFIELD) + " fields";
			morphError(str, 61);
		}
		for(int j=0; j<NFIELD; j++)
			outfile << "\t" << het_stat[i][j];
		outfile << endl;
	}
	outfile.close();
}

// st -> ed
void loadUnliftBed( std::map<int, int> & ref_coord, string & ref_bed_name)
{
	std::ifstream bedfile;
	bedfile.open( ref_bed_name.c_str() );
	if ( !bedfile.is_open() )
		morphErrorFile( ref_bed_name );
	string line;
	while( std::getline(bedfile, line) ) {
		std::stringstream ss;
		ss << line;
		string field;
		std::getline( ss, field, '\t' );
		std::getline( ss, field, '\t' );
		int pos = stoi(field);
		std::getline( ss, field, '\t' );
		ref_coord[pos] = stoi(field);
	}
	bedfile.close();
}


// split the bed file based on window size
// nex is #vcf-record to exclude
void splitBaseBed( string & base_ctrl_bed, string & split_out, int nex )
{
	// count total length first
	std::ifstream inbed;
	string line;
	inbed.open( base_ctrl_bed.c_str() );
	if (!inbed.is_open())
		morphErrorFile(base_ctrl_bed);
	int tlen = 0;
	while(std::getline(inbed, line)) {
		std::stringstream ss;
		ss << line;
		string field;
		std::getline(ss, field, '\t');
		std::getline(ss, field, '\t');
		int st = stoi(field);
		std::getline(ss, field, '\t');
		int ed = stoi(field);
		tlen += (ed - st);
	}
	inbed.close();
	// restrict window number to <6000 (with no vcf exclusion should be ~39133 on chr20)
	int neg_number = N_NEG_REF + nex;
	int padding = floor( (float)tlen / neg_number ); // use distance between windows to adjust

	std::ofstream out;
	inbed.open( base_ctrl_bed.c_str() );
	if (!inbed.is_open())
		morphErrorFile(base_ctrl_bed);
	out.open(split_out.c_str());
	if (!out.is_open())
		morphErrorFile(split_out);
	int nrec = 0;
	while(std::getline(inbed, line)) {
		std::stringstream ss;
		ss << line;
		string field;
		std::getline(ss, field, '\t');
		std::getline(ss, field, '\t');
		int st = stoi(field);
		std::getline(ss, field, '\t');
		int ed = stoi(field);
		int mid = (st + ed) / 2;
		if ( mid - WIN/2 - padding <= st ) {
			nrec++;
			out << REF_CHR << "\t" << mid-WIN/2 << "\t" << mid+WIN/2 << endl;
			continue;
		}
		int n = floor( (float)(mid - WIN/2 - st)/(WIN+padding));
		int beginning = mid - WIN/2 - n*(WIN+padding);
		for(int i=0; i<2*n+1; i++) {
			nrec++;
			out << REF_CHR << "\t" << beginning << "\t" << beginning + WIN << endl;
			beginning += WIN + padding;
		}
		if (nrec > neg_number) // restriction
			break;
	}
	out.close();
	inbed.close();
}


// extract
// bam2fastq
// remap sort clean
void makeRegionBam( map<int, int> & ref_coord, SingleSampleInfo & si, string & remap_bam_prefix)
{
	// extract
	SamFile bam;
	SamFileHeader bam_header;
	OpenBamAndBai( bam, bam_header, si.bam_name);

	string extract_bam_name = remap_bam_prefix + ".extract.bam";
	SamFile extract_bam;
	SamFileHeader extract_bam_header;
	extract_bam.OpenForWrite( extract_bam_name.c_str() );
	if (!extract_bam.IsOpen())
		morphErrorFile(extract_bam_name);
	extract_bam.WriteHeader(bam_header);
	SamRecord rec;

	/* regular region for bwa to estimate ins size: proper only
	string REF_CHR = "20";
	int REF_CHR_ST = 39600000;
	int REF_CHR_ED = 39830000;
	bool section_status = bam.SetReadSection( REF_CHR, REF_CHR_ST, REF_CHR_ED );
	if (!section_status) {
		cerr << "ERROR: cannot set read section: " << REF_CHR << REF_CHR_ST << "-" << REF_CHR_ED << endl;
		exit(32);
	}
	while( bam.ReadRecord(bam_header, rec) ) {
		if ( IsSupplementary(rec.getFlag()) )
			continue;
		if (!rec.getFlag() & 0x2)
			continue;
		if (rec.getReadLength() < avr_read_length/3)
			continue;
		extract_bam.WriteRecord(bam_header, rec);
	}
	*/

	// merge the region first
	vector<int> vst ( ref_coord.size(), 0 );
	vector<int> ved (ref_coord.size(), 0);
	int idx = 0;
	bool first = 1;
	for( map<int, int>::iterator t=ref_coord.begin(); t!=ref_coord.end(); t++ ) {
		int st = t->first - WIN*2;
		int ed = t->second + WIN*2;
		if ( first ) {
			vst[idx] = st;
			ved[idx]  = ed;
			idx++;
			first = 0;
			continue;
		}
		if (st <= ved[idx-1]) {
			if (st <= vst[idx-1])
				morphError("[makeRegionBam] something is wrong in merge");
			ved[idx-1] = ed;
			continue;
		}
		vst[idx] = st;
		ved[idx] = ed;
		idx++;
		continue;
	}
	if (idx == 0)
		morphError("[makeRegionBam] empty ref_coord");
	vst.resize(idx);
	ved.resize(idx);

	// output ref coord region to bam
	SamFile alt_sam;
	SamFileHeader alt_sam_header;
	OpenBamAndBai(alt_sam, alt_sam_header, si.bam_name);
	for( int i=0; i<vst.size(); i++ ) {
		bool section_status = bam.SetReadSection(REF_CHR.c_str(), vst[i], ved[i]);
		if (!section_status) // skip this section
			continue;
		while( bam.ReadRecord(bam_header, rec) ) {
			if ( IsSupplementary(rec.getFlag()) )
				continue;
			if (rec.getReadLength() < si.avr_read_length/3)
				continue;
			if ( OneEndUnmap(rec.getFlag()) )
				continue;
			if (!rec.getFlag() & 0x2) {
				if ( REF_CHR.compare(rec.getMateReferenceName()) == 0 &&  rec.get1BasedPosition() > rec.get1BasedMatePosition())
					continue;
				extract_bam.WriteRecord(bam_header, rec);
				SamRecord mate_rec;
				if( SetDiscMateRec( mate_rec, rec, alt_sam, alt_sam_header) )
					extract_bam.WriteRecord(bam_header, mate_rec);
			}
			else { // for proper reads, try not to output single-end reads
				if ( rec.get1BasedMatePosition() > ved[i] || rec.get1BasedMatePosition()+rec.getReadLength() < vst[i] )
					continue;
				extract_bam.WriteRecord(bam_header, rec);
			}
		}
	}
	alt_sam.Close();
	extract_bam.Close();
	bam.Close();

	// sort by name & bam2fastq
	string cmd = "samtools sort -n " + extract_bam_name + " " + remap_bam_prefix + ".extract.nsort";
	ExecuteCmd( cmd );
	cmd = "bam bam2FastQ --in " + remap_bam_prefix + ".extract.nsort.bam --readName --outBase " + remap_bam_prefix;
	ExecuteCmd(cmd);

	if (CLEAN) {
		cmd = "rm -f " + extract_bam_name + " " + remap_bam_prefix + ".extract.nsort.bam";
		ExecuteCmd(cmd);
	}

	// remap and clean
	remapAndClean( remap_bam_prefix );
}


void remapAndClean( string & remap_bam_prefix )
{
	// remap
	string slice_fasta = PATH + "refs/slice-chr20-hs37d5.fa";
	string cmd = "/net/mario/gotcloud/bin/bwa mem " + slice_fasta + " " + remap_bam_prefix + "_1.fastq ";
	cmd += remap_bam_prefix + "_2.fastq -t 4 > " + remap_bam_prefix + ".aln.sam";
	ExecuteCmd(cmd);
	// sort
	cmd = "samtools sort " + remap_bam_prefix + ".aln.sam " + remap_bam_prefix + ".sort";
	ExecuteCmd(cmd);
	// deup
	cmd = "bam dedup --recab --in " + remap_bam_prefix + ".sort.bam --out " + remap_bam_prefix;
	cmd += ".bam --refFile " + slice_fasta + " --force";
	ExecuteCmd(cmd);
	// index
	cmd = "samtools index " + remap_bam_prefix + ".bam";
	ExecuteCmd(cmd);

	if (!CLEAN)
		return;
	// clean intermediates
	cmd = "rm -f " + remap_bam_prefix + "_1.fastq " + remap_bam_prefix + "_2.fastq " + remap_bam_prefix + ".fastq ";
	cmd += remap_bam_prefix + ".aln.sam " + remap_bam_prefix + ".sort.bam ";
	ExecuteCmd(cmd);
}


// exit with error if command not executed
void ExecuteCmd( string & cmd )
{
	int status = system(cmd.c_str());
	if (status != 0) {
		string str = "fail to execute the following command:\n + cmd";
		morphError(str,1);
	}
}


// extract all "PASS" site by matching column 7
void generatePassVcf( string & raw_vcf_name, string & pass_vcf_name )
{
	string cmd = "cat " + raw_vcf_name + " | head -n 200 | grep \"#\" > " + pass_vcf_name;
	ExecuteCmd(cmd);
	cmd = "cat " + raw_vcf_name + " | awk '$7==\"PASS\"' >> " + pass_vcf_name;
	ExecuteCmd(cmd);
}



