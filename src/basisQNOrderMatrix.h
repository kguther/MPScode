#ifndef MATRIX_TO_CREATE_QN_BLOCK_ORDERING
#define MATRIX_TO_CREATE_QN_BLOCK_ORDERING

#include "mkl_complex_defined.h"
#include "quantumNumber.h"
#include <vector>


//---------------------------------------------------------------------------------------------------//
// The basisQNOrderMatrix class (which has a somewhat strange name) contains the indices of all blocks
// of the MPS matrices for given good quantum numbers. The only necessary input are these quantum numbers
// and the bond dimensions of the MPS. Then, a bunch of vectors containing the indices can be generated.
// Index functions returning the global MPS index for given block indices are supplied.
//---------------------------------------------------------------------------------------------------//


struct multInt{
  int aim;
  int si;
};


class basisQNOrderMatrix{
 public:
  basisQNOrderMatrix(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin);
  basisQNOrderMatrix();
  ~basisQNOrderMatrix();
  void initialize(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin);
  void generateQNIndexTables();
  int blockStructure(int const i, int const direction, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices);
  int aiBlockIndexLP(int const i, int const iBlock, int const j) const{return aiBlockIndicesLP[i][iBlock][j];}
  int siBlockIndexLP(int const i, int const iBlock, int const k) const{return siaimBlockIndicesLP[i][iBlock][k].si;}
  int aimBlockIndexLP(int const i, int const iBlock, int const k) const{return siaimBlockIndicesLP[i][iBlock][k].aim;}
  int aiBlockIndexRP(int const i, int const iBlock, int const k) const{return siaiBlockIndicesRP[i][iBlock][k].aim;}
  int siBlockIndexRP(int const i, int const iBlock, int const k) const{return siaiBlockIndicesRP[i][iBlock][k].si;}
  int aimBlockIndexRP(int const i, int const iBlock, int const j) const{return aimBlockIndicesRP[i][iBlock][j];}
  int aimBlockIndexSplit(int const i, int const iBlock, int const j) const{return aimBlockIndicesSplit[i][j];}
  int aiBlockIndexSplit(int const i, int const iBlock, int const k) const{return aiBlockIndicesSplit[i][k];}
  int siBlockIndexSplit(int const i, int const iBlock, int const j) const{return siBlockIndicesSplit[i][j];}
  int siBlockIndexSplitFixedaim(int const i, int const iBlock, int const k, int const j) const{return siBlockIndicesSplitFixedaim[i][k][j];}
  int lBlockSizeLP(int const i, int const iBlock){return siaimBlockIndicesLP[i][iBlock].size();}
  int rBlockSizeLP(int const i, int const iBlock){return aiBlockIndicesLP[i][iBlock].size();}
  int lBlockSizeRP(int const i, int const iBlock){return aimBlockIndicesRP[i][iBlock].size();}
  int rBlockSizeRP(int const i, int const iBlock){return siaiBlockIndicesRP[i][iBlock].size();}
  int siBlockSizeSplit(int const i, int const iBlock){return siBlockIndicesSplit[i].size();}
  int siBlockSizeSplitFixedaim(int const i, int const iBlock, int const k){return siBlockIndicesSplitFixedaim[i][k].size();}
  int aimBlockSizeSplit(int const i, int const iBlock){return aimBlockIndicesSplit[i].size();}
  int aiBlockSizeSplit(int const i, int const iBlock){return aiBlockIndicesSplit[i].size();}
  int numBlocksLP(int const i){return aiBlockIndicesLP[i].size();}
  int numBlocksRP(int const i){return aimBlockIndicesRP[i].size();}
 private:
  std::vector<quantumNumber> *conservedQNs;
  std::vector<std::vector<int> > *aiBlockIndicesLP;
  std::vector<std::vector<multInt> > *siaimBlockIndicesLP;
  std::vector<std::vector<int> > *aimBlockIndicesRP;
  std::vector<std::vector<multInt> > *siaiBlockIndicesRP;
  std::vector<int > *siBlockIndicesSplit;
  std::vector<int> *aimBlockIndicesSplit;
  std::vector<int> *aiBlockIndicesSplit;
  std::vector<std::vector<int> > *siBlockIndicesSplitFixedaim;
  dimensionTable dimInfo;
  void deleteTables();
  void splitIndexTables(int const i);
  std::complex<int> qnCriterium(int const iQN, int const i, int const aim, int const si, int const direction, int const pre);
};

#endif
