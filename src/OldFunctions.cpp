// only take single-chr out: use in single-chr LHMEI
void TakeChrOut( const char* inBam, const char* outSam, const char* ctrlChr );

void TakeChrOut( const char* inBam, const char* outSam, const char* ctrlChr )
{
	SamFile samIn, samOut;
	samIn.OpenForRead(inBam);
	samOut.OpenForWrite(outSam);
	
	std::string bai = std::string(inBam) + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai.c_str() );
	
	SanityCheckBams( samIn, samOut, bai_status );
	
	SamFileHeader samHeader;
	bool in_bam_header_status = samIn.ReadHeader(samHeader);
	if (!in_bam_header_status) {
		cerr << "ERROR: Fail to read header: " << inBam << endl;
		exit(1);
	}
	samOut.WriteHeader(samHeader);

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
			cerr << "ERROR: Unable to set read section to chr " << current_chr << endl;
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
			cout << "Processed chr " << current_chr << ", total = " << Counter << " reads." << endl;
	}

}
