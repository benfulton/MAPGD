#ifndef POOLEDLIKELIHOOD_H_
#define POOLEDLIKELIHOOD_H_

#include "lnmultinomial.h"
#include "allele_stat.h"
#include <math.h>

void polymorphicmodel(allele_stat const &, float_t *);
void monomorphicmodel(allele_stat const &, float_t *);
void fixedmorphicmodel(allele_stat const &, float_t *);

#endif
