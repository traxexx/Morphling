#include "Wrappers.h"
#include <iostream>
#include <unistd.h>

using std::cout;
using std::endl;

string GetExePath()
{
	string path;
    pid_t pid = getpid();
    char buf[20] = {0};
    sprintf(buf,"%d",pid);
    std::string _link = "/proc/";
    _link.append( buf );
    _link.append( "/exe");
    char proc[512];
    int ch = readlink(_link.c_str(),proc,512);
    if (ch != -1) {
        proc[ch] = 0;
        path = proc;
        std::string::size_type t = path.find_last_of("/");
        path = path.substr(0,t);
    }
    string fullFileName = path + string("/");
    return fullFileName;
}

void DisplayUsageInfo()
{
	cout << endl;
	cout << "Morphling -- A tool for detect Mobile Element Insertions (MEI) in whole genome sequencing data." << endl;
	cout << "(c) 2015 Sai Chen, Goncalo Abecasis" << endl;
	cout << endl;
	cout << "Please go to http://genome.sph.umich.edu/wiki/Morphling for the newest version." << endl;
	cout << endl;
	cout << "  Usage:" << endl;
	cout << "    Discover: Run discovery pipeline for single sample." << endl;
	cout << "    Genotype: Run genotyping pipeline on multi-samples." << endl;
}

void DisplayDiscoveryUsageInfo()
{
	cout << endl;
	cout << "Discover usage: [] is default. Please pecify those * options." << endl;
	cout << "    Use -h to show advanced options." << endl;
	cout << "    Test:  Run test commands. Override all single dash options." << endl;
	cout << endl;
	cout << " *  -Bam:  Raw bam files for LHMEI discovery. MUST be sorted by cooridinate." << endl;
	cout << " *  -WorkDir:  Working directory (if not exist, Morphling will create one)." << endl;
	cout << " *  -Mapper:  Mapper used to generate the raw bam. If not installed globally, please specify path. [/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa-mem]" << endl;
	cout << " *  -GenomeFasta: reference genome fasta. [/net/wonderland/home/saichen/reference/archive/hs37d5.fa]" << endl;
	cout << "    -Sample:  Sample name. [basename of -Bam]" << endl;
	cout << "    -Chr:  Chr to discover. -1 means whole genome. [-1]" << endl;
	cout << "    -Win:  Sliding window size in bp. [600]" << endl;
	cout << "    -Step:  Moving step of sliding window. [100]" << endl;
	cout << "    --passOnly: Only print pass variant. Toggle this when you don't need to run Genotype pipeline. [off]" << endl;
	cout << "    --verbose:  Print user defined options out. [off]" << endl;
	cout << endl;
}

void DisplayDiscoveryDetailedUsageInfo()
{
	DisplayDiscoveryUsageInfo();
	cout << endl;
	cout << "Parameters info. In general no need to change them." << endl;
	cout << "    -CtrlChr:  Slice & remapped chr. [20]" << endl;
	cout << "    -Simplify:  #Evidence reads required to exam a region. [1]" << endl;
	cout << "    -NonOffset: #More non-evidence reads required to exam a region. [1]" << endl;
	cout << "    -MElist:  List of MEI consensus sequences. [refs/MobileElement.list]" << endl;
	cout << "    -MEcoord:  List of MEI genome coordinate beds. [refs/MobileElement.coord]" << endl;
	cout << "    -SliceFA:  Genome fasta of sliced-MEI. [refs/slice-chr20-hs37d5.fa]" << endl;
	cout << "    -HetIndex: Het index file. [refs/hs37d5-chr20-MEI-slice.het-index]" << endl;
	cout << "    --keepIntermediates: Keep intermediate files. [off]" << endl;
	cout << "    --includeSingleAnchor:  Keep hits with single anchor in output vcf. [off]" << endl;
	cout << "    --pseudoChr:  Also run discovery on pseudoChr like GL000215.1. [off]" << endl;
	cout << "    --disableDPfilter: disable depth filter. [off]" << endl;
	cout << "    --noRefAllele: do not output reference allele in output vcf. [off]" << endl;
	cout << "    --noBreakPoint: do not refine break point by #spanning-reads. [off]" << endl;
	cout << endl;
}

void DisplayDiscoveryDebugUsageInfo()
{
	cout << endl;
	cout << "    -ReadLen:  specify read length. [-1]" << endl;
	cout << "    -InsSize:  specify insert size. [-1]" << endl;
	cout << "    -MeiType:  only discover the specific meitype. [-1]" << endl;
	cout << "    --noCtrlVcf:  do not output refLH.*.vcf for ctrl sliced bam. [Off]" << endl;
	cout << "    --debug:  Print detailed running info for debug. [off]" << endl;
	cout << "    --printNonVariant:  print all windows include GT=0/0. [off]" << endl;
	cout << "    --printRefStats:  print out refStats. [off]" << endl;
	cout << "    -refPrefix:  when --printRefStats toggled, prefix of refStats. [refStats]" << endl;
	cout << endl;
}

void DisplayGenotypeUsageInfo()
{
	cout << endl;
	cout << "Genotype usage: [] is default. Please pecify those * options." << endl;
	cout << "    Use -h to show advanced options." << endl;
	cout << "    Test:  Run test commands. Override all single dash options." << endl;
	cout << endl;
	cout << " *  -SampleList:  list of samples to genotype. Cols: sample, bam, vcf_prefix, ctrl_dir_prefix." << endl;
	cout << " *  -WorkDir:  Working directory (if not exist, Morphling will create one)." << endl;
	cout << "    -Chr:  Chr to genotype. -1 means whole genome. [-1]" << endl;
	cout << "    -Win:  Sliding window size in bp. [600]" << endl;
	cout << "    -MeiType:  only discover the specific meitype. -1 means all. [-1]" << endl;
	cout << "    --verbose:  Print user defined options out. [off]" << endl;
	cout << endl;
}

void DisplayGenotypeDetailedUsageInfo()
{

}