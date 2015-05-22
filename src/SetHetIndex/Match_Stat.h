/* Store mappability as bed graph
typedef struct {
	float Align;
    int CordLen;
} MappabilityList;
*/

#include <map>
#include <vector>
#include <string>

typedef struct {
	int chr;
	int win_index;
} CordElement;


typedef struct {
	float align;
	float gc;
	int type; // 1 = HOM, 0 = CTRL, -1 = HOM nearby region (skip), -2 = Non-formative (skip)
	int lift;
	bool used; // label whether it is used in SetHetTable. default = 0
} CtrlTableElement;

typedef struct {
	int chr;
	int hom_index;
	int hom_cord;
	int lift; // lift length
	std::string mei_type;
	CordElement * match_index;
} HetTableElement;

class CtrlTable {
	public:
		CtrlTable(const char * reference, const char * mappability, const char * slice_bed, int win, int step, int chr); // constructor
		void PrintHetIndex(const char * out);
		
	private:
		void ConstructCtrlTable(const char * reference);
		void ConstructMappability(const char * mappability);
		void ConstructGcFromFasta(const char * reference);
		void LabelCtrlTable(const char * slice_bed);
		void SetHetTable(int & chr);
		void setLiftLengthInCtrlTable( int & chr );
	// debug function
		void PrintMappabilityAndGc(const char * out);

		int win_size;
		int step_size;
		int het_element_size_;
		std::map<int, std::vector<CtrlTableElement> > ctrl_table_;
		std::vector<HetTableElement> Het_Table; // cord of neg corresponding to the certain hom
};
