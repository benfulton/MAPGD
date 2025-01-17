#ifndef _ESTIMATE_INDIVIDUAL_H_
#define _ESTIMATE_INDIVIDUAL_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <ctime>

#include <algorithm>
#include <functional>

#include "interface.h" 
#include "pro-file.h" 
#include "individual-likelihood.h"
#include "likelihood.h"
#include "genotype.h"

int estimateInd(int, char **);
allele_stat estimate (Locus &site, models &model, std::vector<float_t> &gofs, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const size_t &MAXPITCH);

#endif 
