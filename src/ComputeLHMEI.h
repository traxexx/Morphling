#include "Options.h"
#include <string> 

using std::string;

void ComputeLHMEI (Options * ptrMainOptions);

void SetGlobalOptions( Options * ptrMainOptions );

void SetGlobalParameters( Options * ptrMainOptions );

void SetReadMapGlobals( Options * ptrMainOptions, string & qinfo_name );

string GetMeiNameFromIndex( int mt );
