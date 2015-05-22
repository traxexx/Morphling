#include <string> 

#include "ComputeLHMEI.h"
#include "Debugs.h"
#include "Utilities.h"

#include "TakeChrOut.h" // generate single-chr ref

void ComputeLHMEI (Options * ptrMainOptions)
{
	int win_len = stoi(ptrMainOptions->ArgMap["Win"]);
	int step_len = stoi(ptrMainOptions->ArgMap["Step"]);
	std::string ref_chr = std::string("20");
	ReadMap * Rmap = new ReadMap( ptrMainOptions->ArgMap["Sample"], win_len, step_len, ref_chr, ptrMainOptions->ArgMap["GenomeFasta"], ptrMainOptions->ArgMap["MElist"]);
	
	std::string work_dir = ptrMainOptions->ArgMap["WorkDir"];
	std::string LHMEI_dir = ptrMainOptions->ArgMap["BinDir"];

/* delete after de-annotate below
	std::string temp_dir = work_dir + "/ctrl_tmps";
	std::string ctrl_bam = temp_dir + "/chr20-remap.sort.recal.bam";
	std::string cmd;
	int cmd_status;
	std::string bam_temp_dir = work_dir + "/bam_tmps";
*/

// set refs first
	std::string cmd = LHMEI_dir + "/bin/Adjust-refs.sh -i " + LHMEI_dir;
	int cmd_status = system(cmd.c_str());
	CheckCmdStatus(cmd, cmd_status);
	
	std::string temp_dir = work_dir + "/ctrl_tmps";
	if ( !ExistDoneFile(temp_dir, "Remap") ) { // do remap
		std::string bwa = std::string("/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa mem ");	
	// ctrl sam generate: mkdir + take out
		cmd = "mkdir -p " + temp_dir;
		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);	
	  // take out
	  	std::string chr20_sam = temp_dir + "/chr20.sam";
	  	bool takeout_status = TakeChrOut( ptrMainOptions->ArgMap["Bam"].c_str(), chr20_sam.c_str(), ref_chr.c_str() ); // 0 is normal
	  	if ( takeout_status ) {
	  		std::cerr << "ERROR: fail to take ref chr out from raw bam!" << std::endl; exit(1);
	  	}
  	// nsort
  		cmd = std::string("samtools sort -n ") + temp_dir + "/chr" + ref_chr + ".sam " + temp_dir + "/chr20.nsort";
  		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);
  	// bam to fastq
		cmd = std::string("bam bam2FastQ --in ") + temp_dir + "/chr20.nsort.bam --readName --outBase " + temp_dir + "/chr20";
		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);
 	 // re-map
  		std::string slice_ref = LHMEI_dir + "/refs/slice-chr20-hs37d5.fa ";
		cmd = bwa + slice_ref + temp_dir + "/chr20_1.fastq " + temp_dir + "/chr20_2.fastq > " + temp_dir + "/align-pe.sam";
		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);
	  // sort by coord
  		cmd = std::string("samtools sort ") + temp_dir + "/align-pe.sam " + temp_dir + "/chr20-remap.sort";
		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);	
	  // dedup
    	cmd = std::string("bam dedup --recab --in ") + temp_dir + "/chr20-remap.sort.bam --out " + temp_dir + "/chr20-remap.sort.recal.bam --force --refFile " + slice_ref + "--storeQualTag OQ --maxBaseQual 40";	
		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);	
		
// done file
		GenerateDoneFile(work_dir, "Remap");
	}
	
// calculate likelihood beased on ctrl
	std::string bam_temp_dir = work_dir + "/bam_tmps";
	cmd = std::string("mkdir -p ") + bam_temp_dir;
	cmd_status = system(cmd.c_str());
	CheckCmdStatus(cmd, cmd_status);
	
	Rmap->SetMapFromBam( ptrMainOptions->ArgMap["Bam"], ptrMainOptions->ArgMap["GenomeFasta"],
		bam_temp_dir, ptrMainOptions->ArgMap["MElist"], ptrMainOptions->ArgMap["MEcoord"]);

// set control (count read type, adjust)
	std::string ctrl_bam = temp_dir + "/chr20-remap.sort.recal.bam";	
	Rmap->SetControls( ctrl_bam, temp_dir );
/*	
	for( int mei_type = 0; mei_type < 2; mei_type++) {	
		Rmap->SetDataPL( mei_type, bam_temp_dir, temp_dir, ptrMainOptions->ArgMap["HetIndex"] );
		std::string outVcf = work_dir + "/Discovery." + getMeiTypeString(mei_type) + ".vcf";
		Rmap->PrintToVcf( outVcf, mei_type );
	}
*/
/*	
// concat vcf
	cmd = std::string("vcf-concat Discovery.Alu.vcf Discovery.L1.vcf Discovery.SVA.vcf > ") + work_dir + "/Unsort-merge.vcf";
	cmd_status = system(cmd.c_str());
	CheckCmdStatus(cmd, cmd_status);
// sort vcf
	cmd = std::string("vcf-sort ") + work_dir + "/Unsort-merge.vcf > " + work_dir + "/Discovery.All.vcf";
	cmd_status = system(cmd.c_str());
	CheckCmdStatus(cmd, cmd_status);
	
//	see if need to remove intermediate files;
	if (ptrMainOptions->OptMap.find("keepInterMediates") == ptrMainOptions->OptMap.end() || ptrMainOptions->OptMap["keepInterMediates"] == 1) {
		std::cout << "  --keepInterMediates not set. Clearing intermediate files..." << std::endl;
		cmd = std::string("rm -f -r ") + temp_dir + " " + work_dir + "/Unsort-merge.vcf " + bam_temp_dir;	
		cmd_status = system(cmd.c_str());
		CheckCmdStatus(cmd, cmd_status);
	}
//	PrintValidRawMapCounts(Rmap);
*/
}




















