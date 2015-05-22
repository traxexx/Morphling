#ifndef QC_H
#define QC_H

#include "SamRecord.h"

// contain read QC info fram a bam section or whole bam
class QC
{
  public:	
	QC();
	~QC();	
  // difference between this and below are whether you add them into QC class
  	bool PassQC( SamRecord & sam_rec );
  // print all private fields
	void PrintQCsummary();

  private:
	int Supplementary_Info;
	int PCR_Duplicates;
	int QC_Fail;
	int Secondary_Alignment;	
};

// check if this read pass general flag QC
bool PassQC( SamRecord & sam_rec );

// check if this read can be used in as disc-read
bool DiscSamPass( SamRecord & sam_rec );

// remove unnecessary fields from sam record before printing to disc sam in discovery stage
void simplifySamRec( SamRecord & sam_rec );

#endif
