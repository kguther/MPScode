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
  int aiBlockIndexLP(int i, int iBlock, int j) const{return aiBlockIndicesLPAccess[reducedIndexFunction(i,iBlock,j)];}
  int siBlockIndexLP(int i, int iBlock, int k) const{return siaimBlockIndicesLPAccess[reducedIndexFunction(i,iBlock,k)].si;}
  int aimBlockIndexLP(int i, int iBlock, int k) const{return siaimBlockIndicesLPAccess[reducedIndexFunction(i,iBlock,k)].aim;}
  int aiBlockIndexRP(int i, int iBlock, int k) const{return siaiBlockIndicesRPAccess[reducedIndexFunction(i,iBlock,k)].aim;}
  int siBlockIndexRP(int i, int iBlock, int k) const{return siaiBlockIndicesRPAccess[reducedIndexFunction(i,iBlock,k)].si;}
  int aimBlockIndexRP(int i, int iBlock, int j) const{return aimBlockIndicesRPAccess[reducedIndexFunction(i,iBlock,j)];}
  int lBlockSizeLP(int const i, int const iBlock) const{return siaimBlockIndicesLP[i][iBlock].size();}
  int rBlockSizeLP(int const i, int const iBlock) const{return aiBlockIndicesLP[i][iBlock].size();}
  int lBlockSizeRP(int const i, int const iBlock) const{return aimBlockIndicesRP[i][iBlock].size();}
  int rBlockSizeRP(int const i, int const iBlock) const{return siaiBlockIndicesRP[i][iBlock].size();}
  int numBlocksLP(int const i) const{return aiBlockIndicesLP[i].size();}
  int numBlocksRP(int const i) const{return aimBlockIndicesRP[i].size();}
  int nQNs() const{return conservedQNs->size();}
 private:
  int maxNumBlocks, maxBlockSize;
  basisQNOrderMatrix(basisQNOrderMatrix const &source);
  basisQNOrderMatrix& operator=(basisQNOrderMatrix const &source);
  std::vector<quantumNumber> *conservedQNs;
  std::vector<int> aiBlockIndicesLPAccess, aimBlockIndicesRPAccess;
  std::vector<multInt> siaimBlockIndicesLPAccess, siaiBlockIndicesRPAccess;
  std::vector<std::vector<int> > *aiBlockIndicesLP;
  std::vector<std::vector<multInt> > *siaimBlockIndicesLP;
  std::vector<std::vector<int> > *aimBlockIndicesRP;
  std::vector<std::vector<multInt> > *siaiBlockIndicesRP;
  dimensionTable dimInfo;
  void generateAccessArrays();
  void deleteTables();
  std::complex<int> qnCriterium(int const iQN, int const i, int const aim, int const si, int const direction, int const pre);
  int reducedIndexFunction(int i, int iBlock, int k) const{return k+iBlock*maxBlockSize+i*maxBlockSize*maxNumBlocks;}
};

#endif
