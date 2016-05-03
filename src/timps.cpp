#include "timps.h"
#include <iostream>

timps::timps(int unitCellSizeIn, dimensionTable const &dimInfoBaseIn, std::vector<quantumNumber> const &conservedQNsin):
  impBase(dimInfoBaseIn),
  unitCellSize(unitCellSizeIn),
  position(0)
{
  unitCell.resize(unitCellSize);
  dimInfoBase.setParameterL(unitCellSize);
  int lDL,lDR;
  int const d=dimInfoBase.d();
  std::vector<int> lDims(3,1);
  //For the first step, no previous results are to be stored, use 1x1 dummy arrays instead
  for(int i=0;i<unitCellSize;++i){
    lDL=dimInfoBase.locDimL(i);
    lDR=dimInfoBase.locDimR(i);
    lDims[0]=d;
    lDims[1]=lDR;
    lDims[2]=lDL;
    unitCell[i]=baseTensor<std::complex<double> >(lDims);
  }
  
  //setup pseudoQuantumNumbers
  conservedQNs.resize(conservedQNsin.size());
  for(int iQN=0;iQN<conservedQNsin.size();++iQN){
    conservedQNs[iQN]=manualQuantumNumber(dimInfoBase,conservedQNsin[iQN].QNValue(),conservedQNsin[iQN].localQNValue());
  }

  //Generate initial labels
  quantumNumber gqn;
  std::vector<std::complex<int> > QNlocal;
  std::complex<int> N;
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    QNlocal=conservedQNsin[iQN].localQNValue();
    N=conservedQNsin[iQN].QNValue();
    //Ensure that the taken indices are for a system of size unitCellSize
    gqn=quantumNumber(dimInfoBase,N,QNlocal);
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
  centralIndexTableVar=twositeQNOrderMatrix(position,dimInfoBase,buf);
  
  //Initialize index table for efficient contractions
  //MAJOR PROBLEM: basisQNORderMatrix is global, BUT CAN ONLY USE LOCAL QNs AT THE UNIT CELL -> USE LOCAL siteQNORDER matrices
  indexTableVar.initialize(dimInfoBase,buf);
  indexTableVar.generateQNIndexTables();
}

//---------------------------------------------------------------------------------------------------//

int timps::addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN, std::vector<std::complex<int> > const &source){
 
  int const deltaL=Lnew-dimInfoBase.L();
  if(deltaL==0){
    return 0;
  }

  std::cout<<"Old matrix dimensions: ";
  for(int j=0;j<unitCell.size();++j){
    std::cout<<unitCell[j].lDL()<<"x"<<unitCell[j].lDR()<<"\t";
  }
  std::cout<<std::endl;

  dimInfoBase.setParameterL(Lnew);
  
  int cpos;
  std::vector<int> locDims(3,0);
  baseTensor<std::complex<double> > dummy;
  for(int dI=0;dI<deltaL;++dI){
    cpos=position+dI+1;
    locDims[0]=dimInfoBase.d();
    locDims[1]=dimInfoBase.locDimR(i+dI);
    locDims[2]=dimInfoBase.locDimL(i+dI);
    dummy=baseTensor<std::complex<double> >(locDims);
    unitCell.insert(unitCell.begin()+cpos,dummy);
  }

  //Remove outdated tensors to save memory
  if(Lnew>unitCellSize){
    unitCell.pop_back();
    unitCell.erase(unitCell.begin());
  }

  //Manually adjust the quantum numbers after the refinement
  int const targetSite=convertPosition(i);
  int const D=dimInfoBase.D();
  conservedQNs[0].indexLabelAccess().insert(conservedQNs[0].indexLabelAccess().begin()+targetSite,source.begin(),source.end());
  if(Lnew>unitCellSize+2){
    conservedQNs[0].indexLabelAccess().erase(conservedQNs[0].indexLabelAccess().begin(),conservedQNs[0].indexLabelAccess().begin()+D);
    conservedQNs[0].indexLabelAccess().erase(conservedQNs[0].indexLabelAccess().begin()+D*(unitCell.size()),conservedQNs[0].indexLabelAccess().end());
  }
  //Shift the U(1)-QN of the right half-system
  int deltaN=targetQN[0].real()-conservedQNs[0].QNValue().real();
  int const rightEnd=targetSite+source.size();
  for(int m=rightEnd;m<(3+unitCellSize)*D;++m){
    conservedQNs[0].indexLabelAccess()[m]+=deltaN;
  }
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    conservedQNs[iQN].setTargetQN(targetQN[iQN]);
  }
  

  //Update the index tables
  setUpTables();
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int timps::refineQN(int i, std::vector<std::complex<int> > const &newQN){
  int const D=dimInfoBase.D();
  int const targetSite=convertPosition(i);
  int const lDL=dimInfoBase.locDimL(i);
  for(int ai=0;ai<lDL;++ai){
    conservedQNs[0].indexLabelAccess()[ai+targetSite*D]=newQN[ai];
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void timps::subMatrixStart(std::complex<double> *&pStart, int i, int si){
  int const targetSite=convertPosition(i);
  unitCell[targetSite].getPtr(pStart,si);
}
