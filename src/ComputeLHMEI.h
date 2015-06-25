#include "Options.h"
#include <string> 

using std::string;

void ComputeLHMEI (Options * ptrMainOptions);

void SetGlobalOptions( Options * ptrMainOptions );

void SetGlobalParameters( Options * ptrMainOptions );

void SetReadMapGlobals( Options * ptrMainOptions, string & qinfo_name );

string GetMeiNameFromIndex( int mt );

int SetQualThreshold( string & rec, string & qc_dir, string & bed_name, string & ex_name );

int getQualThredFromFile( string & qlog_name );


