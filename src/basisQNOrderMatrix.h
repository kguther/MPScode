#ifndef MATRIX_TO_CREATE_QN_BLOCK_ORDERING
#define MATRIX_TO_CREATE_QN_BLOCK_ORDERING

#include "siteQNOrderMatrix.h"
#include "quantumNumber.h"
#include <vector>

//---------------------------------------------------------------------------------------------------//
// The basisQNOrderMatrix class (which has a somewhat strange name) contains the indices of all blocks
// of the MPS matrices for given good quantum numbers. The only necessary input are these quantum numbers
// and the bond dimensions of the MPS. Then, a bunch of vectors containing the indices can be generated.
// Index functions returning the global MPS index for given block indices are supplied.
//---------------------------------------------------------------------------------------------------//

class basisQNOrderMatrix{
 public:
  basisQNOrderMatrix(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin);
  basisQNOrderMatrix(dimensionTable &dimin, std::vector<pseudoQuantumNumber*> const &conservedQNsin);
  basisQNOrderMatrix(int iStart, int iStop, dimensionTable &dimin, std::vector<pseudoQuantumNumber*> const &conservedQNsin);
  basisQNOrderMatrix();
  void generateQNIndexTables();
  int numBlocksLP(int i)const {return localIndexTables[i].numBlocksLP();}
  int numBlocksRP(int i)const {return localIndexTables[i].numBlocksRP();}
  int lBlockSizeLP(int i, int iBlock)const {return localIndexTables[i].lBlockSizeLP(iBlock);}
  int rBlockSizeLP(int i, int iBlock)const {return localIndexTables[i].rBlockSizeLP(iBlock);}
  int lBlockSizeRP(int i, int iBlock)const {return localIndexTables[i].lBlockSizeRP(iBlock);}
  int rBlockSizeRP(int i, int iBlock)const {return localIndexTables[i].rBlockSizeRP(iBlock);}

  int nQNs() const{return localIndexTables[0].nQNs();}
  int validate()const;
  siteQNOrderMatrix const& getLocalIndexTable(int i)const;
 private:
  std::vector<siteQNOrderMatrix> localIndexTables;
  std::vector<pseudoQuantumNumber*> conservedQNs;
  dimensionTable dimInfo;
};

#endif
