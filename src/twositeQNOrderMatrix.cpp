#include "twositeQNOrderMatrix.h"

twositeQNOrderMatrix::twositeQNOrderMatrix(int i, dimensionTable &dimIn, std::vector<quantumNumber> *conservedQNsin):
  conservedQNs(conservedQNsin),
  dimInfo(dimIn),
  site(i)
{
  generateQNIndexTable();
}

//---------------------------------------------------------------------------------------------------//

void twositeQNOrderMatrix::generateQNIndexTable(){
}

//---------------------------------------------------------------------------------------------------//

int twositeQNOrderMatrix::qnCriterium(int iQN, int aim, int air, int si, int sip){
  if((*conservedQNs)[iQN].groupOperation((*conservedQNs)[iQN].QNLabel(site-1,aim),(*conservedQNs)[iQN].QNLabel(si))==(*conservedQNs)[iQN].groupOperation((*conservedQNs)[iQN].QNLabel(site+1,air),(*conservedQNs)[iQN].QNLabel(sip),-1)){
    return 1;
  }
  return 0;
}
