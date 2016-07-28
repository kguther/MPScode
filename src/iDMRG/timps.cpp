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
  reducedIndexLabelQNs.resize(conservedQNs.size());
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    QNlocal=conservedQNsin[iQN].localQNValue();
    N=conservedQNsin[iQN].QNValue();
    //Ensure that the taken indices are for a system of size unitCellSize
    gqn=quantumNumber(dimInfoBase,N,QNlocal);
    conservedQNs[iQN].indexLabelAccess()=gqn.indexLabelAccess();
    reducedIndexLabelQNs[iQN]=conservedQNs[iQN];
  }
  setUpTwositeTable();
}

//---------------------------------------------------------------------------------------------------//

void timps::setUpSingleSiteTables(){  
  std::vector<pseudoQuantumNumber*> buf;
  int const D=dimInfoBase.D();
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    buf.push_back(&(conservedQNs[iQN]));
  }
 
  //Initialize index table for efficient contractions
  //Generate index tables for the unit cell
  int const lDL=dimInfoBase.locDimL(currentSite());
  int const lDR=dimInfoBase.locDimR(currentSite());
  int const ld=dimInfoBase.d();
  int const lDLp=dimInfoBase.locDimL(currentSite()+unitCellSize-1);
  int const lDRp=dimInfoBase.locDimR(currentSite()+unitCellSize-1);
  leftTable=siteQNOrderMatrix(position,lDL,lDR,ld,buf);
  //The labels have to be taken after refinement. Then, left/right labels have been separated, leading to the +1 instead of -1
  rightTable=siteQNOrderMatrix(position+unitCellSize+1,lDLp,lDRp,ld,buf);
}

//---------------------------------------------------------------------------------------------------//

void timps::setUpTwositeTable(){
  //The centralIndexTable has to be set at another point than the left/right Tables (because it is used for something else)
  std::vector<pseudoQuantumNumber*> buf;
  int const D=dimInfoBase.D();
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    //This can only be done after the first step
    if(dimInfoBase.L()>unitCellSize){
    //Copy only the central index labels to the reduced QNs. This always includes only the left/right selection made in the last optimization
      reducedIndexLabelQNs[iQN].indexLabelAccess()=std::vector<std::complex<int> >(conservedQNs[iQN].indexLabelAccess().begin()+D,conservedQNs[iQN].indexLabelAccess().begin()+(2+unitCellSize)*D);
    }
    buf.push_back(&(reducedIndexLabelQNs[iQN]));
  }
  
  //Initialize twosite index table for efficient optimization
  centralIndexTableVar=twositeQNOrderMatrix(position,dimInfoBase,buf);
}

//---------------------------------------------------------------------------------------------------//

int timps::addSite(int Lnew, int i){

  //This comes last in the setup for the next step, here, the previous results are deleted and replaced with empty tensors of the next dimensions. 
  //If the previous results are still required, they have to be backuped prior to this step
 
  int const deltaL=Lnew-dimInfoBase.L();
  if(deltaL!=2){
    return 1;
    //ADD EXCEPTION
  }

  std::cout<<"Old matrix dimensions:\n";
  for(int j=0;j<unitCell.size();++j){
    std::cout<<unitCell[j].lDL()<<"x"<<unitCell[j].lDR()<<"\t";
  }
  std::cout<<std::endl;

  dimInfoBase.setParameterL(Lnew);

  setUpTwositeTable();
  
  std::vector<int> locDims(3,0);
  for(int dI=0;dI<unitCellSize;++dI){
    locDims[0]=dimInfoBase.d();
    //i is a coordinate in the new system
    locDims[1]=dimInfoBase.locDimR(i+dI);
    locDims[2]=dimInfoBase.locDimL(i+dI);
    unitCell[dI]=baseTensor<std::complex<double> >(locDims);
  }

  return 0;
}

//---------------------------------------------------------------------------------------------------//

int timps::refineQN(int i, std::vector<std::complex<int> > const &leftSideLabels, std::vector<std::complex<int> > const &rightSideLabels, std::vector<std::complex<int> > const &targetQN){
  int const D=dimInfoBase.D();
  int const targetSite=position+1;
  //Manually adjust the quantum number labels
  std::vector<std::complex<int> > &indexLabel=conservedQNs[0].indexLabelAccess();

  /*
  std::cout<<std::endl;
  std::cout<<"Pre-refinement:\n";
  for(int i=0;i<unitCellSize+3;++i){
    for(int ai=0;ai<D;++ai){
      std::cout<<indexLabel[ai+D*i]<<"\t";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<"Refining with:\n";
  int const lDR=dimInfoBase.locDimL(i);
  for(int m=0;m<lDR;++m){
    std::cout<<leftSideLabels[m]<<"\t"<<rightSideLabels[m]<<std::endl;
  }
  std::cout<<std::endl;
  */

  //ADD EXCEPTION
  //both input labels have to be of size D

  //Cut the outdated labels
  if(dimInfoBase.L()!=unitCellSize){
    indexLabel.erase(indexLabel.begin(),indexLabel.begin()+D);
    indexLabel.erase(indexLabel.begin()+(indexLabel.size()-D),indexLabel.end());
  }

  indexLabel.insert(indexLabel.begin()+targetSite*D,leftSideLabels.begin(),leftSideLabels.end());
  indexLabel.insert(indexLabel.begin()+(targetSite+unitCellSize)*D,rightSideLabels.begin(),rightSideLabels.end());
  //The middle labels are just left where they are

  //Shift the U(1)-QN of the right half-system
  int const deltaN=targetQN[0].real()-conservedQNs[0].QNValue().real();
  int const rightEnd=(targetSite+unitCellSize)*D;
  
  if(rightEnd>=indexLabel.size())
    std::cout<<"CRITICAL ERROR: INVALID INPUT\n";
 
  
  for(int m=rightEnd;m<indexLabel.size();++m){
    indexLabel[m]+=deltaN;
  }
  for(int iQN=0;iQN<conservedQNs.size();++iQN){
    conservedQNs[iQN].setTargetQN(targetQN[iQN]);
  }

  /*
  std::cout<<std::endl;
  std::cout<<std::endl;
  for(int i=0;i<unitCellSize+3;++i){
    for(int ai=0;ai<D;++ai){
      std::cout<<indexLabel[ai+D*i]<<"\t";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<std::endl;
  */

  //Update the index tables
  setUpSingleSiteTables();
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void timps::subMatrixStart(std::complex<double> *&pStart, int i, int si){
  int const targetSite=convertPosition(i);
  unitCell[targetSite].getPtr(pStart,si);
}
