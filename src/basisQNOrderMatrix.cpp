#include "basisQNOrderMatrix.h"
#include <iostream>

//---------------------------------------------------------------------------------------------------//
// There is not really much to see here, a bunch of vectors is supplied which store the indices of 
// the blocks of the MPS matrices.
//---------------------------------------------------------------------------------------------------//

basisQNOrderMatrix::basisQNOrderMatrix():
  nQN(0)
{}

//---------------------------------------------------------------------------------------------------//
// The basisQNOrderMatrix only relies on the functionality of pseudoQuantumNumbers, therefore, we store pointers to those internally. There still has to be a constructor and an initialize() function for quantumNumbers though, since these are called by the mps (could be adjusted, but has no priority). 
//---------------------------------------------------------------------------------------------------//

basisQNOrderMatrix::basisQNOrderMatrix(dimensionTable const &dimin, std::vector<quantumNumber> *conservedQNsin):
  dimInfo(dimin),
  nQN(conservedQNsin->size())
{
  conservedQNs.resize(nQN);
  for(int m=0;m<nQN;++m){
    conservedQNs[m]=&((*conservedQNsin)[m]);
  }
  generateQNIndexTables();
}

//---------------------------------------------------------------------------------------------------//

basisQNOrderMatrix::basisQNOrderMatrix(dimensionTable const &dimin, std::vector<pseudoQuantumNumber*> const &conservedQNsin):
  dimInfo(dimin),
  conservedQNs(conservedQNsin),
  nQN(conservedQNsin.size())
{
  generateQNIndexTables();
}

//---------------------------------------------------------------------------------------------------//

basisQNOrderMatrix::basisQNOrderMatrix(int iStart, int iStop, dimensionTable const &dimin, std::vector<pseudoQuantumNumber*> const &conservedQNsin):
  dimInfo(dimin),
  conservedQNs(conservedQNsin),
  nQN(conservedQNsin.size())
{
  int cI;
  localIndexTables.resize(iStop-iStart+1);
  for(int i=0;i<localIndexTables.size();++i){
    cI=i+iStart;
    //Write only a few local index tables, for sites in range from iStart to iStop (including both)
    localIndexTables[i]=siteQNOrderMatrix(cI,dimInfo.locDimL(cI),dimInfo.locDimR(cI),dimInfo.locd(cI),conservedQNs);
    //ADD EXCEPTION
  }
}

//---------------------------------------------------------------------------------------------------//

siteQNOrderMatrix const& basisQNOrderMatrix::getLocalIndexTable(int i)const {
  //ADD EXCEPTION
  return localIndexTables[i];
}

//---------------------------------------------------------------------------------------------------//

siteQNOrderMatrix& basisQNOrderMatrix::getLocalIndexTable(int i) {
  //ADD EXCEPTION
  return localIndexTables[i];
}


//---------------------------------------------------------------------------------------------------//
// Generates the tables containing the virtual bond indices for the block indices. Returns 0 at success, 1 if the quantum number targeted is invalid and 2 if the labeling scheme itself is invalid (i.e. contains non-normalizable blocks).
//---------------------------------------------------------------------------------------------------//

void basisQNOrderMatrix::generateQNIndexTables(){
  int L=dimInfo.L();
  int cumulativeBlockSize;
  localIndexTables.resize(L);
  for(int i=0;i<L;++i){
    localIndexTables[i]=siteQNOrderMatrix(i,dimInfo.locDimL(i),dimInfo.locDimR(i),dimInfo.locd(i),conservedQNs);
  }
  
  //Check if labeling scheme is valid
  int info=validate();
  if(info && 0){
    std::cout<<"CRITICAL ERROR: Invalid QN labeling scheme at site "<<abs(info)-1<<"\n";
    for(int i=0;i<dimInfo.L();++i){
      if(i+1==abs(info)){
	siteQNOrderMatrix output=localIndexTables[i];
	// This part is used to test QN labeling schemes for their useability. It prints out the block indices and their QN labels.
	std::cout<<"Right labels:\n";
	for(int aim=0;aim<dimInfo.locDimL(i+1);++aim){
	  std::cout<<aim<<" with label "<<(conservedQNs[0])->QNLabel(i,aim)<<std::endl;
	}
	std::cout<<"Left labels:\n";
	for(int aim=0;aim<dimInfo.locDimL(i);++aim){
	  std::cout<<aim<<" with label "<<(conservedQNs[0])->QNLabel(i-1,aim)<<std::endl;
	}
	if(info>0){
	  for(int iBlock=0;iBlock<numBlocksRP(i);++iBlock){
	    std::cout<<"Right indices: "<<std::endl;
	    for(int j=0;j<rBlockSizeRP(i,iBlock);++j){
	      std::cout<<output.aiBlockIndexRP(iBlock,j)<<" with label "<<(conservedQNs[0])->QNLabel(i,output.aiBlockIndexRP(iBlock,j))<<"\t"<<output.siBlockIndexRP(iBlock,j)<<" with label "<<(conservedQNs[0])->QNLabel(output.siBlockIndexRP(iBlock,j))<<std::endl;
	    }
	    std::cout<<"Left indices: \n";
	    for(int j=0;j<lBlockSizeRP(i,iBlock);++j){
	      std::cout<<output.aimBlockIndexRP(iBlock,j)<<" with label "<<(conservedQNs[0])->QNLabel(i-1,output.aimBlockIndexRP(iBlock,j))<<std::endl;
	    }
	  }
	}
	else{
	  for(int iBlock=0;iBlock<numBlocksLP(i);++iBlock){
	    std::cout<<"Left indices: "<<std::endl;
	    for(int j=0;j<lBlockSizeLP(i,iBlock);++j){
	      std::cout<<output.aimBlockIndexLP(iBlock,j)<<" with label "<<(conservedQNs[0])->QNLabel(i-1,output.aimBlockIndexLP(iBlock,j))<<"\t"<<output.siBlockIndexLP(iBlock,j)<<" with label "<<(conservedQNs[0])->QNLabel(output.siBlockIndexLP(iBlock,j))<<std::endl;
	    }
	    std::cout<<"Right indices: \n";
	    for(int j=0;j<rBlockSizeLP(i,iBlock);++j){
	      std::cout<<output.aiBlockIndexLP(iBlock,j)<<" with label "<<(conservedQNs[0])->QNLabel(i,output.aiBlockIndexLP(iBlock,j))<<std::endl;
	    }
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// The validate() function returns -1 if there are non-normalizable blocks and 0 else.
//---------------------------------------------------------------------------------------------------//

int basisQNOrderMatrix::validate()const {
  int info=0;
  for(int i=0;i<dimInfo.L();++i){
    info+=localIndexTables[i].validate();
  }
  return info;
}
