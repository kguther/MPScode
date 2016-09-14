#ifndef TEST_FUNCTIONS_FOR_QN_CONSTRAINT
#define TEST_FUNCTIONS_FOR_QN_CONSTRAINT

#include "mps.h"
#include <arcomp.h>

//int checkQNConstraint(impBase &test);
 
//This part is for iDMRG
int checkQNConstraint(mps &test, int i);
/*
int checkQNConstraint(impBase &test, int i);
int twositeCheck(impBase &test, mpsEntryType *M);
*/
#endif
