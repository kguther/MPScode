#ifndef TEST_FUNCTIONS_FOR_QN_CONSTRAINT
#define TEST_FUNCTIONS_FOR_QN_CONSTRAINT

#include "impBase.h"
#include <arcomp.h>

int checkQNConstraint(impBase &test);
int checkQNConstraint(impBase &test, int i);
int twositeCheck(impBase &test, arcomplex<double> *M);
#endif
