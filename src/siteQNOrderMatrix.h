#ifndef SINGLE_SITE_INDEX_TABLE
#define SINGLE_SITE_INDEX_TABLE

#include "multInt.h"
#include "quantumNumber.h"

class siteQNOrderMatrix{
 public:
  siteQNOrderMatrix(){}
  siteQNOrderMatrix(int site, int lDL, int lDR, int ld, std::vector<quantumNumber> *conservedQNsin);
  siteQNOrderMatrix(int site, int lDL, int lDR, int ld, std::vector<pseudoQuantumNumber*> const &conservedQNsin);
  void blockStructure(int direction, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices);
  int aiBlockIndexLP(int iBlock, int j) const{return aiBlockIndicesLP[iBlock][j];}
  int siBlockIndexLP(int iBlock, int k) const{return siaimBlockIndicesLP[iBlock][k].si;}
  int aimBlockIndexLP(int iBlock, int k) const{return siaimBlockIndicesLP[iBlock][k].aim;}
  int aiBlockIndexRP(int iBlock, int k) const{return siaiBlockIndicesRP[iBlock][k].aim;}
  int siBlockIndexRP(int iBlock, int k) const{return siaiBlockIndicesRP[iBlock][k].si;}
  int aimBlockIndexRP(int iBlock, int j) const{return aimBlockIndicesRP[iBlock][j];}
  int lBlockSizeLP(int iBlock) const{return siaimBlockIndicesLP[iBlock].size();}
  int rBlockSizeLP(int iBlock) const{return aiBlockIndicesLP[iBlock].size();}
  int lBlockSizeRP(int iBlock) const{return aimBlockIndicesRP[iBlock].size();}
  int rBlockSizeRP(int iBlock) const{return siaiBlockIndicesRP[iBlock].size();}
  int numBlocksLP() const{return aiBlockIndicesLP.size();}
  int numBlocksRP() const{return aimBlockIndicesRP.size();}
  int nQNs() const{return conservedQNs.size();}
  int validate()const;
 private:
  int i, lDL, lDR, ld;
  std::vector<std::vector<int> > aiBlockIndicesLP;
  std::vector<std::vector<multInt> > siaimBlockIndicesLP;
  std::vector<std::vector<int> > aimBlockIndicesRP;
  std::vector<std::vector<multInt> > siaiBlockIndicesRP;
  std::vector<pseudoQuantumNumber*> conservedQNs;
  std::complex<int> qnCriterium(int iQN, int aim, int si, int direction, int pre);
  void loadConservedQNs(std::vector<quantumNumber> *conservedQNsin);
  void setUpTable();
};

#endif
