#include "Sites.h"

void Sites::AssemblySubtype()
{
// something useful
	string out_vcf = ptrMainOptions->ArgMap["Out"];
	string mei_list = ptrMainOptions->ArgMap["MElist"];

// open out
	ofstream out_vcf;
	out_vcf.open( out_vcf_name.c_str() );
	CheckOutFileStatus( out_vcf, out_vcf_name.c_str() );
// ref mei fasta
	vector<RefSeq> rSeq;
	rSeq.resize(3);
	for( int i=0; i<3; i++ )
		rSeq;	
	
// header
	PrintHeader( out_vcf );
// loo through each site
	for( vector<PotSite>::iterator ps = SiteInfo.begin(); ps != SiteInfo.end(); ps++ ) {
		AsbSite currentSite( *ps, bamList, rSeq );
		currentSite.Assembly();
		currentSite.Print( out_vcf );
	}
	out_vcf.close();
}
