#include "timps.h"
#include <iostream>

timps::timps(int unitCellSizeIn, dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin):
  impBase(dimInfoIn),
  dimInfo(dimInfoIn),
  unitCellSize(unitCellSizeIn),
  position(unitCellSize/2-1)
{
  unitCell.resize(unitCellSize);
  dimInfo.setParameterL(unitCellSize);
  int lDL,lDR;
  int const d=dimInfo.d();
  std::vector<int> lDims(3,0);
  for(int i=0;i<unitCellSize;++i){
    lDL=dimInfo.locDimL(i);
    lDR=dimInfo.locDimR(i);
    lDims[0]=d;
    lDims[1]=lDR;
    lDims[2]=lDL;
    unitCell[i]=baseTensor<std::complex<double> >(lDims);
  }
  
  //setup pseudoQuantumNumbers
  conservedQNs.resize(conservedQNsin.size());
  for(int iQN=0;iQN<conservedQNsin.size();++iQN){
    conservedQNs[iQN]=manualQuantumNumber(dimInfo,conservedQNsin[iQN].QNValue(),conservedQNsin[iQN].localQNValue());
  }

  //Generate initial labels
  quantumNumber gqn;
  std::vector<std::complex<int> > QNlocal;
  std::complex<int> N;
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    QNlocal=conservedQNsin[iQN].localQNValue();
    N=conservedQNsin[iQN].QNValue();
    //Ensure that the taken indices are for a system of size unitCellSize
    gqn=quantumNumber(dimInfo,N,QNlocal);
    conservedQNs[iQN].indexLabelAccess()=gqn.indexLabelAccess();
  }
  setUpTables();

}

//---------------------------------------------------------------------------------------------------//

void timps::setUpTables(){
  //Initialize twosite index table for efficient optimization
  std::vector<pseudoQuantumNumber*> buf;
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    buf.push_back(&(conservedQNs[iQN]));
  }
  centralIndexTableVar=twositeQNOrderMatrix(position,dimInfo,buf);
  
  //Initialize index table for efficient contractions
  indexTableVar.initialize(dimInfo,buf);
  indexTableVar.generateQNIndexTables();
}

//---------------------------------------------------------------------------------------------------//

int timps::addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN, std::vector<std::complex<int> > const &source){
 
  int const deltaL=Lnew-dimInfo.L();
  if(deltaL==0){
    return 0;
  }

  std::cout<<"Old matrix dimensions: ";
  for(int i=0;i<unitCell.size();++i){
    std::cout<<unitCell[i].lDL()<<"x"<<unitCell[i].lDR()<<"\t";
  }
  std::cout<<std::endl;

  dimInfo.setParameterL(Lnew);
  
  int cpos;
  std::vector<int> locDims(3,0);
  baseTensor<std::complex<double> > dummy;
  for(int dI=0;dI<deltaL;++dI){
    cpos=position+dI+1;
    locDims[0]=dimInfo.d();
    locDims[1]=dimInfo.locDimR(i+dI);
    locDims[2]=dimInfo.locDimL(i+dI);
    dummy=baseTensor<std::complex<double> >(locDims);
    unitCell.insert(unitCell.begin()+cpos,dummy);
  }
  //Remove outdated tensors to save memory
  unitCell.pop_back();
  unitCell.erase(unitCell.begin());

  //Manually adjust the quantum numbers after the refinement
  int const targetSite=convertPosition(i)+Lnew-dimInfo.L();
  int const D=dimInfo.D();
  conservedQNs[0].indexLabelAccess().insert(conservedQNs[0].indexLabelAccess().begin()+targetSite,source.begin(),source.end());
  conservedQNs[0].indexLabelAccess().erase(conservedQNs[0].indexLabelAccess().begin(),conservedQNs[0].indexLabelAccess().begin()+D);
  conservedQNs[0].indexLabelAccess().erase(conservedQNs[0].indexLabelAccess().begin()+D*(unitCell.size()),conservedQNs[0].indexLabelAccess().end());

  //Update the index tables
  setUpTables();
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int timps::refineQN(int i, std::vector<std::complex<int> > const &newQN){
  int const D=dimInfo.D();
  int const targetSite=convertPosition(i);
  int const lDL=dimInfo.locDimL(i);
  for(int ai=0;ai<lDL;++ai){
    conservedQNs[0].indexLabelAccess()[ai+targetSite*D]=newQN[ai];
  }
  return 0;
}

void timps::subMatrixStart(lapack_complex_double *&pStart, int i, int si){
  int const targetSite=convertPosition(i);
  unitCell[targetSite].getPtr(pStart,si);
}
