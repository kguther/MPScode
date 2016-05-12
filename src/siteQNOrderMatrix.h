#ifndef SINGLE_SITE_INDEX_TABLE
#define SINGLE_SITE_INDEX_TABLE

#include "multInt.h"
#include "quantumNumber.h"

//---------------------------------------------------------------------------------------------------//
// For the iDMRG, it proved to be useful to have as much stuff local as possible. Therefore, the index
// tables are now stored in separate objects for each site. There is no performance penalty, because
// in all applications, only the indices of a single site have to be considered.
//---------------------------------------------------------------------------------------------------//

class siteQNOrderMatrix{
 public:
  siteQNOrderMatrix(){}
  siteQNOrderMatrix(int site, int lDL, int lDR, int ld, std::vector<quantumNumber> *conservedQNsin);
  siteQNOrderMatrix(int site, int lDL, int lDR, int ld, std::vector<pseudoQuantumNumber*> const &conservedQNsin);
  void generateFull(siteQNOrderMatrix const &source);
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
  std::complex<int> qnLabelLP(int iBlock) const{return qnLabelsLP[iBlock];}
  std::complex<int> qnLabelRP(int iBlock) const{return qnLabelsRP[iBlock];}
  int nQNs() const{return conservedQNs.size();}
  int validate()const;
 private:
  int i, lDL, lDR, ld;
  std::vector<std::vector<int> > aiBlockIndicesLP;
  std::vector<std::vector<multInt> > siaimBlockIndicesLP;
  std::vector<std::vector<int> > aimBlockIndicesRP;
  std::vector<std::vector<multInt> > siaiBlockIndicesRP;
  std::vector<pseudoQuantumNumber*> conservedQNs;
  std::vector<std::complex<int> > qnLabelsLP, qnLabelsRP;
  std::complex<int> qnCriterium(int iQN, int aim, int si, int direction, int pre);
  void loadConservedQNs(std::vector<quantumNumber> *conservedQNsin);
  void setUpTable();
  //In contrast to setUpTable(), setUpTableFull() gets the labels by combining the pair indices, leading to an overhead and empty blocks. 
  //This is not useable in the optimization (will even fail) but needed in the enrichment
  void blockStructureFull(int direction, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices);
  void setUpTableFull();
};

void printIndexTable(siteQNOrderMatrix const &a);

#endif
