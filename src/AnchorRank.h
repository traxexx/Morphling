#ifndef ANCHORRANK_H
#define ANCHORRANK_H

#include <vector>

using std::vector;

int GetVariantPosterior( vector<float> & GL);

int GetSupportReadFraction( vector<int> & counts, int depth );

int GetProperReadFraction( vector<int> & counts, int depth );

int getSumSupportClips( vector<int> & counts );

int getSumSupportDiscs( vector<int> & counts);

int getSumSupportUnmaps( vector<int> & counts);

#endif