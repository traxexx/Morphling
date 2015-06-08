#ifnedf ABSSITE_H
#define ABSSITE_H

// for storing mapping info of each read on each mei subtype
struct EviInfo
{
	int Boundary; // if minus, do not consider any mapping results on left; other wise right.
	int LAign;
	int RAlign
	int Score; // SW score
	string Seq; // sequence. if map to '-' strand, revert it
}

class AsbSite {
 public:
 	AsbSite( PotSite & ps, BamList & blist, vector<RefSeq*> & rs_vec );
 	void Assembly();
 	void Print( ofstream & out_vcf );
 	
 private:
 	int Meitype;
 	int Position;
 	bool* GtList;
 	vector<SamFile> * bFiles; // point to opened bam files
 	RefSeq * rSeq; // decide by mei type
 	vector< vector< map<int, EviInfo> > > Clusters; // left(+/-) right (+/-). strand based on MEI sequence (NOT genomic flanking)
 	string str_before_format;
 	string str_after_format;
 	int subtype;
 	int left_most;
 	int right_most;
	int missing_base;
	int basecount;
};

#endif

