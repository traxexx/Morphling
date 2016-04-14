#ifndef SINGLESAMPLE_H
#define SINGLESAMPLE_H

#include <string>

typedef struct {
	std::string sample_name;
	std::string bam_name;
	float depth;
	int avr_read_length;
	int avr_ins_size;
	float var_avr_ins_size;
} SingleSampleInfo;

#endif