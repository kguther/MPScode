#include "parameters.h"

problemParameters::problemParameters(localHSpaces const &din, int Lin, int Dwin, int nEigsin, int NumberQNs, std::complex<int> *QNconservedin, std::complex<int> * QNListin, double tReal, double tImag): d(din),L(Lin),Dw(Dwin),nEigs(nEigsin),nQNs(NumberQNs),t(std::complex<double>(tReal,tImag))
{
  filling.resize(nQNs);
  for(int iQN=0;iQN<nQNs;++iQN){
    filling[iQN]=QNconservedin[iQN].real()/static_cast<double>(L);
  }
  QNLocalList.resize(nQNs);
  for(int nQN=0;nQN<nQNs;++nQN){
    QNconserved.push_back(QNconservedin[nQN]);
    for(int si=0;si<d.maxd();++si){
      QNLocalList[nQN].push_back(QNListin[si+nQN*d.maxd()]);
    }
  }
}
