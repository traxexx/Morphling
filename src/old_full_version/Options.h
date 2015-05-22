#include <string>
#include <map>

/*
option format: type(bool: - or --) name(string) default(string) NO SPACE!!!
	option string: space for no default value. ArgMap only
		-step=100;-win=600;-bam= ;
	dummies: if --option exist, these options are invoked
		--read_bam:-bam,-mei;--read_merge:-merge;
		
Exception: when --verbose, display options
*/

class Options
{
  public:
	Options(int argc, char * argv[], std::string & arg_string, std::string & dummies);
	std::map<std::string, std::string> ArgMap;
	std::map<std::string, bool> OptMap;

  private:
  	void displayOptions();
  	void verifyOptions(); // if not pass, exit right now
};
