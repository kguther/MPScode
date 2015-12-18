#ifndef MATRIX_TO_CREATE_QN_BLOCK_ORDERING
#define MATRIX_TO_CREATE_QN_BLOCK_ORDERING

#include "mkl_complex_defined.h"
#include "quantumNumber.h"
#include <vector>

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
  int aimBlockIndexSplit(int const i, int const iBlock, int const j) const{return aimBlockIndicesSplit[i][iBlock][j];}
  int siBlockIndexSplit(int const i, int const iBlock, int const j) const{return siBlockIndicesSplit[i][iBlock][j];}
  int lBlockSizeLP(int const i, int const iBlock){return siaimBlockIndicesLP[i][iBlock].size();}
  int rBlockSizeLP(int const i, int const iBlock){return aiBlockIndicesLP[i][iBlock].size();}
  int lBlockSizeRP(int const i, int const iBlock){return aimBlockIndicesRP[i][iBlock].size();}
  int rBlockSizeRP(int const i, int const iBlock){return siaiBlockIndicesRP[i][iBlock].size();}
  int siBlockSizeSplit(int const i, int const iBlock){return siBlockIndicesSplit[i][iBlock].size();}
  int aimBlockSizeSplit(int const i, int const iBlock){return aimBlockIndicesSplit[i][iBlock].size();}
  int numBlocksLP(int const i){return aiBlockIndicesLP[i].size();}
  int numBlocksRP(int const i){return aimBlockIndicesRP[i].size();}
 private:
  std::vector<quantumNumber> *conservedQNs;
  std::vector<std::vector<int> > *aiBlockIndicesLP;
  std::vector<std::vector<multInt> > *siaimBlockIndicesLP;
  std::vector<std::vector<int> > *aimBlockIndicesRP;
  std::vector<std::vector<multInt> > *siaiBlockIndicesRP;
  std::vector<std::vector<int> > *siBlockIndicesSplit;
  std::vector<std::vector<int> > *aimBlockIndicesSplit;
  dimensionTable dimInfo;
  void deleteTables();
  void splitIndexTables(int const i);
};

#endif
