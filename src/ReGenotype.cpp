#include <map>
#include "ReGenotype.h"
#include "QC.h" // PassQC
#include "Globals.h" // WIN
#include "OriginalStats.h" // MergeCell, SetRecordGL
#include <utility>

using std::map;

void ResetVcfRecordFromBam( VcfRecord & vcf_rec, RefStats & rstats, vector<RefSeq*> & REF_SEQ, string & chr, int center, SamFile & samIn, SamFileHeader & samHeader )
{
	vector<int> raw_counts; // read-type count vector
	raw_counts.resize(18, 0);
	setReadCountInSection( raw_counts, chr, center, samIn, samHeader, REF_SEQ ); // set from bam

// set gl & if 1/1 or 1/0, set breakpoint
	vector<MergeCell> new_vec;
	new_vec.resize(1);
	MergeCellPtr merge_ptr = new_vec.begin();
	merge_ptr->dups = 1;
	merge_ptr->counts = raw_counts;
	raw_counts.clear();
	merge_ptr->GL.resize(3,0);
	rstats.SetRecordGL( merge_ptr );

	vcf_rec.SetChrName( chr );
	vcf_rec.SetPosition( center );
	vcf_rec.UpdateFromMergeCellPtr( merge_ptr );
	
// if need to refine break point && exist MEI, do refine
	if (REFINE_BREAK_POINT && vcf_rec.GetDosage() > 0)
		vcf_rec.SetBreakPointAndCIFromBam( samIn, samHeader );
}


void setReadCountInSection( vector<int> & raw_counts, string & chr, int center, SamFile & samIn, SamFileHeader & samHeader, vector<RefSeq*> & REF_SEQ )
{
	int st = center - WIN/2;
	int ed = center + WIN/2;
	bool section_status = samIn.SetReadSection( chr.c_str(), st, ed );
	if (!section_status) {
		std::cerr << "Warning: Unable to set read section: " << chr << ": " << st << "-" << ed << ". Set section depth = 0!" << std::endl;
		return;
	}

// proper reads	
	map<string, vector< DiscPair > > disc_recs; // record where the disc come from
	ProperDeck pDeck( REF_SEQ );
	SamRecord sam_rec;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
		bool pass_qc = PassQC( sam_rec );
		if ( !pass_qc )
			continue;
		if ( sam_rec.getFlag() & 0x2 ) {
			if ( sam_rec.getInsertSize() > 0 )
				pDeck.Add( sam_rec );
			else {
				RetrievedIndex rv = pDeck.RetrieveIndex( sam_rec );
				int index = getRetrievedIndexToRawCounts( rv );
				if (index >= 0) { // for read partially in window, only add clip
					if ( sam_rec.get1BasedPosition() < st || sam_rec.get1BasedAlignmentEnd() > ed ) {
						if ( index >= 2 )
							raw_counts[ index ]++;
					}
					else
						raw_counts[ index ]++;
				}
			}
		}
		else { // disc: rec info and wait to reset section later
			string mate_chr = sam_rec.getMateReferenceName();
		// check if this one is valid as anchor
			DiscPair new_pair( 1, 0, sam_rec, REF_SEQ );
			disc_recs[mate_chr].push_back( new_pair );
		}
	}
	
// disc reads
	for( map<string, vector< DiscPair > >::iterator chr_it = disc_recs.begin(); chr_it != disc_recs.end(); chr_it++ ) {
		for( vector< DiscPair >::iterator dp_it = chr_it->second.begin(); dp_it != chr_it->second.end(); dp_it++ ) {
			bool section_status = samIn.SetReadSection( chr_it->first.c_str(), dp_it->GetSecondAlignPosition(), dp_it->GetSecondAlignPosition() + WIN );
			if (!section_status) {
				std::cerr << "ERROR: Unable to set read section: " << chr << ": " << st << "-" << ed << std::endl;
				exit(1);
			}
			SamRecord sam_rec;
			while( samIn.ReadRecord(samHeader, sam_rec) ) {
				bool pass_qc = PassQC( sam_rec );
				if ( !pass_qc )
					continue;
				int position = sam_rec.get1BasedPosition();
				if ( position > dp_it->GetFirstAlignPosition() )
					break;
				if ( sam_rec.getFlag() & 0x2 )
					continue;
				bool same_pair = dp_it->IsSamePair( sam_rec );
				if ( !same_pair )
					continue;
			// now add ro raw stats: always use first as anchor
				dp_it->AddSecondToPair( sam_rec );
				Loci loci = dp_it->GetFirstLoci();
				if ( dp_it->FirstValid() ) {
					int index = getLociToRawCounts( loci );
					raw_counts[ index ]++;
				// clear & break
					break;		
				}
			}
		}
	}		
}

/*
  if invalid pair, return -1.
  Decide from rv.type
  Only mei_index sequence is included in REF_SEQ, so if position 3~5 contains 1, it is mei clip.
  Other rules follow in ReadMap
*/
int getRetrievedIndexToRawCounts( RetrievedIndex & rv )
{
	if ( !rv.valid )
		return -1;
		
	int real_type;
	if ( rv.type <= 2 ) // only proper
		real_type = (rv.type == 0) ? 0 : 1;
	else { // exist clip
		bool mei = (rv.type >> 2) > 0 ? 1 : 0;
		bool isBegin = rv.type & 4;
		real_type = ( (mei << 1) | isBegin ) & 3;
		if ( rv.mei_on_1st )
			real_type += 4;
	}
	return real_type;
}

/*
	-1 still may exist: 2nd too short, 2nd do not pass qc
	Only mei_index sequence is included in REF_SEQ, so if position 1~3 contains 1, it is mei disc.
	Other rules follow in ReadMap
*/
int getLociToRawCounts( Loci & loci)
{
	int real_type = loci.type > 0 ? 0 : 1;
	if (loci.sense_strand)
		real_type += 2;
	if (loci.mateMap)
		real_type += 4;
	return real_type;
}

