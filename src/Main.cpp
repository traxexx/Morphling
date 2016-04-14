#include <iostream>
#include "MultiSampleCalling.h"

void DisplayUsageInfo();

int main(int argc, char * argv[])
{
// if need to print help info only
	if (argc <= 1) {
		DisplayUsageInfo();
		return 0;
	}

	MultiSampleCalling( argc, argv );

	return 0; // successful execution
}

void DisplayUsageInfo()
{
	std::cout << "\nUsage:" << std::endl;
	std::cout << "    -SampleList    list of samples." << std::endl;
	std::cout << "    -OutVcf        output vcf name." << std::endl;
	std::cout << "[Optional]" << std::endl;
	std::cout << "    -Log           name of log file." << std::endl;
	std::cout << "    -SiteVcf       provide site vcf. Start directly from genotyping step." << std::endl;
	std::cout << "    -rStat         provide a list of ref stat. Skip generating ref stat step." << std::endl;
	std::cout << "    --siteOnly     only generate site vcf. No genotyping." << std::endl;
	std::cout << "    --statOnly     only generate ref stat." << std::endl;
	std::cout << "    --noClean      do not clean intermediate files." << std::endl;
//	std::cout << "    " << std::endl;
	std::cout << "\n" << std::endl;
}
