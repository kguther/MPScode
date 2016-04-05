#infdef TWO_SITE_QN_ORDER_MATRIX
#define TWO_SITE_QN_ORDER_MATRIX

#include "basisQNOrderMatrix.h"
#include "dimensionTable.h"

class twositeQNOrderMatirx{
 public:
  twositeQNOrderMatrix(int i, dimensionTable &dimIn, std::vector<quantumNumber> *conservedQNsin);
  void generateQNIndexTable();
  int aimBlockIndex(int iBlock, int k)const {return lBlockIndices[iBlock][k].aim;}
  int airBlockIndex(int iBlock, int j)const {return rBlockIndices[iBlock][j].aim;}
  int siBlockIndex(int iBlock, int k)const {return lBlockIndices[iBlock][k].si;}
  int sipBlockIndex(int iBlock, int j)const {return rBlockIndices[iBlock][j].si;}
  int lBlockSize(int iBlock)const {return lBlockIndices[iBlock].size();}
  int rBlockSize(int iBlock)const {return rBlockIndices[iBlock].size();}
  int numBlocks()const {return lBlockIndices.size();}
  int nQNs()const {return conservedQNs->size();}
 private:
  std::vector<std::vector<multInt> > lBlockIndices, rBlockIndices;
  std::vector<quantumNumber> *conservedQNs;
  dimensionTable dimInfo;
  int site;
  int qnCriterium(int iQN, int aim, int air, int si, int sip);
}

#endif
