#ifndef TWO_SITE_QN_ORDER_MATRIX
#define TWO_SITE_QN_ORDER_MATRIX

#include "basisQNOrderMatrix.h"
#include "dimensionTable.h"
#include "pseudoQuantumNumber.h"
#include <vector>

//---------------------------------------------------------------------------------------------------//
// A class containing index tables for the entries of a twosite compound matrix that fullfill the qn constraint (i.e. the twosite constraint). It works much like the basisQNOrderMatrix, however, there is no left/right pairing since the pairing is always si/aim and sip/air in our case. The twosite approach is used in iDMRG. Since the table has to be generated before each step, it would be highly inefficient to make it global. Therefore, this table only describes the indices for two sites instead of the whole MPS.
//---------------------------------------------------------------------------------------------------//


class twositeQNOrderMatrix{
 public:
  twositeQNOrderMatrix();
  twositeQNOrderMatrix(int i, dimensionTable const &dimIn, pseudoQuantumNumber *conservedQNsin);
  int generateQNIndexTable();
  int aimBlockIndex(int iBlock, int k)const {return lBlockIndices[iBlock][k].aim;}
  int airBlockIndex(int iBlock, int j)const {return rBlockIndices[iBlock][j].aim;}
  int siBlockIndex(int iBlock, int k)const {return lBlockIndices[iBlock][k].si;}
  int sipBlockIndex(int iBlock, int j)const {return rBlockIndices[iBlock][j].si;}
  int lBlockSize(int iBlock)const {return lBlockIndices[iBlock].size();}
  int rBlockSize(int iBlock)const {return rBlockIndices[iBlock].size();}
  std::complex<int> blockQN(int iQN, int iBlock)const {return qnLabels[iQN][iBlock];}
  int numBlocks()const {return qnLabels[0].size();}
  //For now, only a single QN can be enforces here
  int nQNs()const {return 1;}
  int getSite()const {return site;}
 private:
  std::vector<std::vector<multInt> > lBlockIndices, rBlockIndices;
  std::vector<std::vector<std::complex<int> > > qnLabels;
  pseudoQuantumNumber *conservedQNs;
  dimensionTable dimInfo;
  int site;
  std::complex<int> qnCriterium(int iQN, int i, int ai, int si, int pre=1);
  void writeIndexTables(int i, int ld, int lD, std::vector<std::vector<multInt> > &target);
};

#endif
