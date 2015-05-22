#include "StringArray.h"
#include "StringHash.h"
#include "Parameters.h"
#include "Error.h"

#include "SamFile.h"
#include "SamValidation.h"

#include <string>
#include <vector>
#include <iostream>
#include <utility>


bool TakeChrOut( const char* inBam, const char* outSam, const char* ctrlChr );

int main(int argc, char * argv[])
{
	bool status = TakeChrOut("/net/wonderland/home/saichen/LHMEI_v4/usage_test/original.test.bam", "test/out.sam", "20");
	return status;
}

bool TakeChrOut( const char* inBam, const char* outSam, const char* ctrlChr )
{
	SamFile samIn, samOut;
	samIn.OpenForRead(inBam);
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Unable to open input bam: " << inBam << std::endl;
		exit(1);
	}

	samOut.OpenForWrite(outSam);
	if ( !samOut.IsOpen() ) {
		std::cerr << "ERROR: Unable to open output bam: " << outSam << std::endl;
		exit(1);
	}
	
	SamFileHeader samHeader;
	samIn.ReadHeader(samHeader);
	samOut.WriteHeader(samHeader);
	
	std::string bai = std::string(inBam) + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai.c_str() );
	if ( !bai_status ) {
		std::cerr << "ERROR: Unable to read bam index: " << outSam << ".bai, is it indexed?" << std::endl;
		exit(1);
	}

// list of non_ctrl chr
	std::string REF_CHR = std::string(ctrlChr);
	std::vector<std::string> chrNameVec;
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String val = samHeader.getReferenceLabel(i);
		std::string str_val = std::string(val.c_str());
		chrNameVec.push_back(str_val);
	}
	
// process bam
	for(int i=0; i< samHeader.getNumSQs(); i++) {
		std::string current_chr = chrNameVec[i];	
		bool section_status = samIn.SetReadSection(current_chr.c_str());
		if (!section_status) {
			std::cerr << "ERROR: Unable to set read section to chr " << current_chr << std::endl;
			exit(1);
		}
		
		int Counter = 0;
		SamRecord samRecord;
		
	// ref section
		if ( current_chr.compare(REF_CHR) == 0) {
			 while(samIn.ReadRecord(samHeader, samRecord)) {
			 	Counter++;
         		samOut.WriteRecord(samHeader, samRecord);
			 }
		}
		else { // non-ref section
 			while(samIn.ReadRecord(samHeader, samRecord)) {
         		Counter++;
         		int flag = samRecord.getFlag();
         		if (flag & 0x2)  continue;
         		if (flag & 0x8)  continue; // mate unmap
         		std::string mate_chr = samRecord.getMateReferenceName();
         		if (mate_chr.compare(REF_CHR) == 0)
	         		samOut.WriteRecord(samHeader, samRecord);
      		}		
		}
		if (current_chr.length() <= 5)
			std::cout << "Processed chr " << current_chr << ", total = " << Counter << " reads." << std::endl;
	}

	return 0;
}




