#include "mps.h"
#include "arrayprocessing.h"
#include "mkl_complex_defined.h"
#include "exceptionClasses.h"
#include <cmath>
#include <memory>
#include <iostream>

void translateBlocks(baseTensor<std::complex<double> > const &source, siteQNOrderMatrix const &sourceTable, baseTensor<std::complex<double> > &target, siteQNOrderMatrix const &targetTable);

mps::mps():stateArray()
{}

//---------------------------------------------------------------------------------------------------//

mps::mps(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin){
  stateArray::initialize(dimInfoIn);
  setUpQNs(conservedQNsin);
  createInitialState();
}

//---------------------------------------------------------------------------------------------------//

mps::mps(mps const &source):
  stateArray(source),
  conservedQNs(source.conservedQNs)
{
  loadIndexTables();
}

//---------------------------------------------------------------------------------------------------//

mps& mps::operator=(mps const &source){
  stateArray::operator=(source);
  conservedQNs=source.conservedQNs;
  loadIndexTables();
  return *this;
}

//---------------------------------------------------------------------------------------------------//

void mps::generate(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin){
  stateArray::initialize(dimInfoIn);
  setUpQNs(conservedQNsin);
  createInitialState();
}

//---------------------------------------------------------------------------------------------------//

void mps::setUpQNs(std::vector<quantumNumber> const &conservedQNsin){
  //Often-used function for setting the internal quantum numbers of an mps
  conservedQNs=conservedQNsin;
  loadIndexTables();
}

//---------------------------------------------------------------------------------------------------//

void mps::loadIndexTables(){
  int info=0;
  if(conservedQNs.size()){
    nQNs=conservedQNs.size();
  }
  else{
    nQNs=0;
    info=-1;
  }
  if(nQNs){
    indexTableVar=basisQNOrderMatrix(dimInfo,&conservedQNs);
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::loadIndexTablesNoexcept(){
  if(conservedQNs.size()){
    nQNs=conservedQNs.size();
  }
  else{
    nQNs=0;
  }
  if(nQNs){
    for(int i=0;i<dimInfo.L();++i){
      try{
	indexTableVar.getLocalIndexTable(i)=siteQNOrderMatrix(i,dimInfo.locDimL(i),dimInfo.locDimR(i),dimInfo.locd(i),&conservedQNs);
      }
      catch(empty_table &err){
	//if a local table is empty, the labels on the next site have to be adapted to create non-empty blocks
	for(int iQN=0;iQN<nQNs;++iQN){
	  if(i!=(dimInfo.L()-1)){
	    //this actually can not happen on the last site except for very small L or D since there is no freedom in picking the labels for the last site for reasonably large systems
	    conservedQNs[iQN].adaptLabels(i,1);
	  }
	}
	--i;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::refineQNLabels(int i, int iQN, std::vector<std::complex<int> > const &source){
  conservedQNs[iQN].refine(i,source);
  loadIndexTables();
}

//---------------------------------------------------------------------------------------------------//

void mps::adaptLabels(int i, int direction){
  for(int iQN=0;iQN<nQNs;++iQN){
    conservedQNs[iQN].adaptLabels(i,direction);
  }
  loadIndexTablesNoexcept();
}

//---------------------------------------------------------------------------------------------------//

void mps::createInitialState(){
  int lDL, lDR, ld;
  if(!nQNs){
    for(int i=0;i<L;++i){
      lDL=locDimL(i);
      lDR=locDimR(i);
      lDL=(lDR<lDL)?lDR:lDL;
      ld=locd(i);
      for(int si=0;si<ld;++si){
	for(int aim=0;aim<lDL;++aim){
	  global_access(i,si,aim,aim)=1;
	}
      }
    }
  }
  else{
    //This sets all entries with indices fullfilling the QN constraint to 1.
    int numBlocks, lBlockSize, rBlockSize;
    siteQNOrderMatrix tmp;
    for(int i=0;i<L;++i){
      tmp=indexTableVar.getLocalIndexTable(i);
      numBlocks=tmp.numBlocksLP();
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	rBlockSize=tmp.rBlockSizeLP(iBlock);
	lBlockSize=tmp.lBlockSizeLP(iBlock);
	for(int j=0;j<rBlockSize;++j){
	  for(int k=0;k<lBlockSize;++k){
	    global_access(i,tmp.siBlockIndexLP(iBlock,k),tmp.aiBlockIndexLP(iBlock,j),tmp.aimBlockIndexLP(iBlock,k))=1;
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

int mps::setParameterD(int Dnew){
  //this is insufficient as the old labels and the new labels have different index tables
  //-> entries belonging to some block are now shuffled arbitrarily
  dimInfo.setParameterD(Dnew);
  for(int iQN=0;iQN<nQNs;++iQN){
    conservedQNs[iQN].setParameterD(Dnew);
  }
  if(nQNs){
    basisQNOrderMatrix const backupIndexTableVar=indexTableVar;
    indexTableVar=basisQNOrderMatrix(dimInfo,&conservedQNs);
    std::vector<int> siteDims(3);
    for(int i=0;i<dimInfo.L();++i){
      siteDims[2]=dimInfo.locDimL(i);
      siteDims[1]=dimInfo.locDimR(i);
      siteDims[0]=dimInfo.locd(i);
      baseTensor<std::complex<double> > newSiteMatrix(siteDims);
      translateBlocks(getStateArrayEntry(i),backupIndexTableVar.getLocalIndexTable(i),newSiteMatrix,indexTableVar.getLocalIndexTable(i));
      setStateArrayEntry(i,newSiteMatrix);
    }
    return 0;
  }
  else{
  return stateArray::setParameterD(Dnew);
  }
}

//---------------------------------------------------------------------------------------------------//

int mps::setParameterL(int Lnew){
  dimInfo.setParameterL(Lnew);
  for(int iQN=0;iQN<nQNs;++iQN){
    conservedQNs[iQN].setParameterL(Lnew);
  }
  if(nQNs){
    indexTableVar=basisQNOrderMatrix(dimInfo,&conservedQNs);
  }
  return stateArray::setParameterL(Lnew);
}

//---------------------------------------------------------------------------------------------------//
// The following functions are left/right normalizing the matrices of a site after 
// optimization and multiplying the remainder to the matrices of the site to the left/right.
//---------------------------------------------------------------------------------------------------//

void mps::leftNormalizeState(int i){
  if(nQNs){
    leftNormalizeStateBlockwise(i);
  }
  else{
    leftNormalizePrimitive(i);
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::rightNormalizeState(int i){
  if(nQNs){
    return rightNormalizeStateBlockwise(i);
  }
  else{
    rightNormalizePrimitive(i);
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::normalizeFinal(int i){
  lapack_complex_double normalization;
  int site, lcD,ld;
  if(i){
    site=0;
    lcD=locDimR(0);
  }
  else{
    site=L-1;
    lcD=locDimL(L-1);
  }
  ld=locd(site);
  //Normalize last matrices to maintain normalization of state
  lapack_complex_double *finalArray;
  subMatrixStart(finalArray,site);
  normalization=cblas_dznrm2(ld*lcD,finalArray,1);
  normalization=1.0/normalization;
  cblas_zscal(ld*lcD,&normalization,finalArray,1);
}

//---------------------------------------------------------------------------------------------------//
// Now, these two functions also left-/rightnormalize a site matrix, but they use the block structure
// of the MPS matrices when using conserved QNs. In particular, they only make sense when using conserved
// QNs. There, the normalization is executed on each block individually, preserving the QN constraint. 
// IMPORTANT: IN GENERAL, THE BLOCKS ARE NOT SQUARE. THE QN LABELING SCHEME MUST BE DESIGNED SUCH
// THAT EACH BLOCK CAN BE BROUGHT INTO CANONICAL FORM. THIS IS A NONTRIVIAL CONSTRAINT.
//---------------------------------------------------------------------------------------------------//

void mps::leftNormalizeStateBlockwise(int i){
  int ld, lDR, lDL, lDRR;
  int lBlockSize,rBlockSize;
  int aiCurrent, siCurrent, aimCurrent, aipCurrent;
  lapack_complex_double *M, *R;
  lapack_complex_double *Rcontainer, *Qcontainer;
  lapack_complex_double *inputA;
  lapack_complex_double zone=1.0, zzero=0.0;
  lapack_int info;
  lDL=locDimL(i);
  lDR=locDimR(i);
  lDRR=locDimR(i+1);
  ld=locd(i);
  std::unique_ptr<lapack_complex_double[]> RP(new lapack_complex_double[lDR*lDR]);
  std::unique_ptr<lapack_complex_double[]> MP, QP, RcP;
  R=RP.get();
  siteQNOrderMatrix const localIndexTable=indexTableVar.getLocalIndexTable(i);
  for(int iBlock=0;iBlock<localIndexTable.numBlocksLP();++iBlock){
    rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
    lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
    //rBlockSize is required to be smaller than or equal to lBlockSize
    if(rBlockSize!=0 && lBlockSize!=0){
      //We do not normalize empty blocks
      MP.reset(new lapack_complex_double[lBlockSize*rBlockSize]);
      M=MP.get();
      //First, copy the content of the block to some dummy array
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  convertIndicesLP(localIndexTable,j,k,iBlock,siCurrent,aiCurrent,aimCurrent);
	  M[k+j*lBlockSize]=global_access(i,siCurrent,aiCurrent,aimCurrent);
	}
      }
      RcP.reset(new lapack_complex_double[rBlockSize*rBlockSize]);
      QP.reset(new lapack_complex_double[rBlockSize]);
      Rcontainer=RcP.get();
      Qcontainer=QP.get();
      //Make a QR decomposition of that dummy array
      info=LAPACKE_zgeqrf(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zgeqrf: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
      }
      auxiliary::lowerdiag(rBlockSize,rBlockSize,M,Rcontainer,lBlockSize);
      info=LAPACKE_zungqr(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zungqr: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
      }
      //Copy the resulting Q-matrix back into the block of the MPS matrix
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  convertIndicesLP(localIndexTable,j,k,iBlock,siCurrent,aiCurrent,aimCurrent);
	  global_access(i,siCurrent,aiCurrent,aimCurrent)=M[k+j*lBlockSize];
	}
	//And the resulting R-matrix into the corresponding elements of a uncompressed matrix R
	for(int l=0;l<rBlockSize;++l){
	  aiCurrent=localIndexTable.aiBlockIndexLP(iBlock,j);
	  aipCurrent=localIndexTable.aiBlockIndexLP(iBlock,l);
	  R[aipCurrent+aiCurrent*lDR]=Rcontainer[l+j*rBlockSize];
	}
      }
    }
  }
  //Finally, multiply the uncompressed matrix R into the next MPS matrix. Block structure does not need to be adressed explicitly, since they both have the required structure.
  std::unique_ptr<lapack_complex_double[]> inputAP(new lapack_complex_double[lDR*lDRR]);
  inputA=inputAP.get();
  lapack_complex_double *localMatrix;
  for(int si=0;si<ld;++si){
    subMatrixStart(localMatrix,i+1,si);
    auxiliary::arraycpy(lDRR,lDR,localMatrix,inputA);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,lDR,&zone,R,lDR,inputA,lDR,&zzero,localMatrix,lDR);
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::rightNormalizeStateBlockwise(int i){
  int ld, lDR, lDL, lDLL;
  int lBlockSize,rBlockSize;
  int aiCurrent, siCurrent, aimCurrent, aimpCurrent;
  lapack_complex_double *M, *R;
  lapack_complex_double *Rcontainer, *Qcontainer;
  lapack_complex_double *inputA;
  lapack_complex_double zone=1.0, zzero=0.0;
  lapack_int info;
  lDL=locDimL(i);
  lDR=locDimR(i);
  lDLL=locDimL(i-1);
  ld=locd(i);
  std::unique_ptr<lapack_complex_double[]> RP(new lapack_complex_double[lDL*lDL]);
  std::unique_ptr<lapack_complex_double[]> MP, RcP, QP;
  R=RP.get();
  siteQNOrderMatrix const localIndexTable=indexTableVar.getLocalIndexTable(i);
  for(int iBlock=0;iBlock<indexTableVar.numBlocksRP(i);++iBlock){
    lBlockSize=localIndexTable.lBlockSizeRP(iBlock);
    rBlockSize=localIndexTable.rBlockSizeRP(iBlock);
    //lBlockSize is required to be smaller than or equal to rBlockSize
    if(rBlockSize!=0 && lBlockSize!=0){
      MP.reset(new lapack_complex_double[lBlockSize*rBlockSize]);
      M=MP.get();
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  convertIndicesRP(localIndexTable,j,k,iBlock,siCurrent,aiCurrent,aimCurrent);
	  M[k+j*lBlockSize]=global_access(i,siCurrent,aiCurrent,aimCurrent);
	}
      }
      RcP.reset(new lapack_complex_double[lBlockSize*lBlockSize]);
      QP.reset(new lapack_complex_double[rBlockSize]);
      Rcontainer=RcP.get();
      Qcontainer=QP.get();
      info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zgerqf: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
      }
      auxiliary::lowerdiag(lBlockSize,lBlockSize,M+(rBlockSize-lBlockSize)*lBlockSize,Rcontainer);
      info=LAPACKE_zungrq(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,lBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zungrq: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
      }
      for(int j=0;j<lBlockSize;++j){
	for(int k=0;k<rBlockSize;++k){
	  convertIndicesRP(localIndexTable,k,j,iBlock,siCurrent,aiCurrent,aimCurrent);
	  global_access(i,siCurrent,aiCurrent,aimCurrent)=M[j+k*lBlockSize];
	}
	for(int l=0;l<lBlockSize;++l){
	  aimCurrent=localIndexTable.aimBlockIndexRP(iBlock,j);
	  aimpCurrent=localIndexTable.aimBlockIndexRP(iBlock,l);
	  R[aimpCurrent+aimCurrent*lDL]=Rcontainer[l+j*lBlockSize];
	}
      }
    }
  }
  std::unique_ptr<lapack_complex_double[]> inputAP(new lapack_complex_double[lDL*lDLL]);
  lapack_complex_double *localMatrix;
  inputA=inputAP.get();
  for(int si=0;si<ld;++si){
    subMatrixStart(localMatrix,i-1,si);
    auxiliary::arraycpy(lDLL*lDL,localMatrix,inputA);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,lDL,&zone,inputA,lDLL,R,lDL,&zzero,localMatrix,lDLL);
  }
}

//---------------------------------------------------------------------------------------------------//
// Simple, non-blockwise version of left/right normalization. Does not require a valid QN labeling 
// scheme
//---------------------------------------------------------------------------------------------------//

void mps::leftNormalizePrimitive(int i){
  lapack_int info;
  int D1, D2, D3, ld;
  //Yep, the local dimensions are just named D1, D2, D3 - the decomposition is of a d*D1xD2 matrix
  D1=locDimL(i);
  D2=locDimR(i);
  D3=locDimR(i+1);
  ld=locd(i);
  lapack_complex_double *Rcontainer, *Qcontainer, *localMatrix;
  const lapack_complex_double zone=1.0;
  std::unique_ptr<lapack_complex_double[]> Qp(new lapack_complex_double[D2]);
  Qcontainer=Qp.get();//Used for storage of lapack-internal matrices
  std::unique_ptr<lapack_complex_double[]> Rp(new lapack_complex_double[D2*D2]);//Used for storage of R from RQ decomposition
  Rcontainer=Rp.get();
  //Enable use of LAPACK_ROW_MAJOR which is necessary here due to the applied storage scheme
  for(int si=0;si<ld;++si){
    subMatrixStart(localMatrix,i,si);
    auxiliary::transp(D2,D1,localMatrix);
  }
  //Use thin QR decomposition
  subMatrixStart(localMatrix,i);
  info=LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,ld*D1,D2,localMatrix,D2,Qcontainer);
  auxiliary::upperdiag(D2,D2,localMatrix,Rcontainer);
  //Only first D2 columns are used -> thin QR (below icrit, this is equivalent to a full QR)
  info=LAPACKE_zungqr(LAPACK_ROW_MAJOR,ld*D1,D2,D2,localMatrix,D2,Qcontainer);
  if(info){
    std::cout<<"Error in LAPACKE_zungqr: "<<info<<" At site: "<<i<<" With dimensions: "<<D2<<"x"<<D1<<" and local Hilbert space dimension: "<<ld<<std::endl;
  }
  auxiliary::transp(D2,D2,Rcontainer);
  for(int si=0;si<ld;++si){
    subMatrixStart(localMatrix,i,si);
    auxiliary::transp(D1,D2,localMatrix);
    subMatrixStart(localMatrix,i+1,si);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,localMatrix,D2);
    //here, R is packed into the matrices of the next site
  }                                             
}

//---------------------------------------------------------------------------------------------------//

void mps::rightNormalizePrimitive(int i){
  lapack_int info;
  //This time, we decompose a D2xd*D1 matrix and multiply the D2xD2 matrix R with a D3xD2 matrix
  int D1,D2,D3,ld;
  D1=locDimR(i);
  D2=locDimL(i);
  D3=locDimL(i-1);
  ld=locd(i);
  lapack_complex_double *Rcontainer, *Qcontainer, *localMatrix;
  const lapack_complex_double zone=1.0;
  std::unique_ptr<lapack_complex_double[]> QP(new lapack_complex_double[ld*D1]);
  std::unique_ptr<lapack_complex_double[]> RP(new lapack_complex_double[D2*D2]);
  Qcontainer=QP.get();
  Rcontainer=RP.get();
  //Thats how zgerqf works: the last D2 columns contain the upper trigonal matrix R, to adress them, move D2 from the end
  subMatrixStart(localMatrix,i);
  info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,D2,ld*D1,localMatrix,D2,Qcontainer);
  //lowerdiag does get an upper trigonal matrix in column major ordering, dont get confused
  auxiliary::lowerdiag(D2,D2,localMatrix+D2*(ld*D1-D2),Rcontainer);
  info=LAPACKE_zungrq(LAPACK_COL_MAJOR,D2,ld*D1,D2,localMatrix,D2,Qcontainer);
  if(info){
    std::cout<<"ERROR IN LAPACKE_zungrq:"<<info<<" At site: "<<i<<" With dimensions: "<<D2<<"x"<<D1<<" and local Hilbert space dimension: "<<ld<<std::endl;
  }
  for(int si=0;si<ld;++si){
    subMatrixStart(localMatrix,i-1,si);
    cblas_ztrmm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,D3,D2,&zone,Rcontainer,D2,localMatrix,D3);
  }  
}

//---------------------------------------------------------------------------------------------------//
// These functions convert block-internal indices to normal MPS bond indices.
//---------------------------------------------------------------------------------------------------//

void mps::convertIndicesLP(siteQNOrderMatrix const &localIndexTable, int j, int k, int iBlock, int &si, int &ai, int &aim){
  ai=localIndexTable.aiBlockIndexLP(iBlock,j);
  aim=localIndexTable.aimBlockIndexLP(iBlock,k);
  si=localIndexTable.siBlockIndexLP(iBlock,k);
}

//---------------------------------------------------------------------------------------------------//

void mps::convertIndicesRP(siteQNOrderMatrix const &localIndexTable, int j, int k, int iBlock, int &si, int &ai, int &aim){
  aim=localIndexTable.aimBlockIndexRP(iBlock,k);
  ai=localIndexTable.aiBlockIndexRP(iBlock,j);
  si=localIndexTable.siBlockIndexRP(iBlock,j);
}

//---------------------------------------------------------------------------------------------------//
// Outdated function to enforce the QN constraint on some MPS.
//---------------------------------------------------------------------------------------------------//

void mps::restoreQN(int i){
  int lDL, lDR, ld;
  for(int iQN=0;iQN<nQNs;++iQN){
    ld=locd(i);
    lDR=locDimR(i);
    lDL=locDimL(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  if(real((conservedQNs[iQN].QNLabel(i,ai)-conservedQNs[iQN].QNLabel(i-1,aim)-conservedQNs[iQN].QNLabel(si))) && abs(global_access(i,si,ai,aim))>0.000001){
	    global_access(i,si,ai,aim)=0;
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// Functions for obtaining the entanglement spectrum and entropy of the MPS
//---------------------------------------------------------------------------------------------------//

void mps::getEntanglementSpectrumOC(int i, double &S, std::vector<double> &spectrum){
  //compute the entanglement spectrum at some site - this requires the site to be the orthogonality center
  int ld, lDR, lDL;
  ld=locd(i);
  lDL=locDimL(i);
  lDR=locDimR(i);
  lapack_complex_double *currentM;
  std::unique_ptr<double> diagsP(new double[lDR]);
  double *diags=diagsP.get();
  std::unique_ptr<lapack_complex_double[]> AnewP(new lapack_complex_double[lDL*lDR*ld]);
  lapack_complex_double *Anew=AnewP.get();
  subMatrixStart(currentM,i);
  auxiliary::arraycpy(ld*lDL*lDR,currentM,Anew);
  LAPACKE_zgesdd(LAPACK_COL_MAJOR,'N',ld*lDL,lDR,Anew,ld*lDL,diags,0,1,0,1);
  spectrum.clear();
  for(int m=0;m<lDR;++m){
    spectrum.push_back(diags[m]);
  }
  S=0;
  for(int m=0;m<spectrum.size();++m){
    if(std::abs(spectrum[m])>1e-12){
      S+=spectrum[m]*spectrum[m]*log(spectrum[m]*spectrum[m]);
    }
  }
  S*=-1;
}

//---------------------------------------------------------------------------------------------------//

void mps::getEntanglementEntropy(std::vector<double> &S, std::vector<std::vector<double> > &spectra){
  //compute the entanglement spectrum and entropy at all sites
  S.resize(L);
  spectra.resize(L);
  for(int i=L-1;i>0;--i){
    rightNormalizePrimitive(i);
  }
  normalizeFinal(1);
  for(int i=0;i<L;++i){
    getEntanglementSpectrumOC(i,S[i],spectra[i]);
    if(i<L-1)
      leftNormalizePrimitive(i);
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::getEntanglementSpectrum(int i, double &S, std::vector<double> &spectra){
  for(int j=L-1;j>0;--j){
    rightNormalizeState(j);
  }
  normalizeFinal(1);
  for(int j=0;j<i;++j){
    leftNormalizeState(j);
  }
  getEntanglementSpectrumOC(i,S,spectra);
}

//---------------------------------------------------------------------------------------------------//

//For insight into the algorithm
void printQNLabels(mps const &test){
  quantumNumber gqn=test.getConservedQNs()[0];
  dimensionTable dimInfo=test.getDimInfo();
  for(int i=0;i<dimInfo.L()+1;++i){
    std::cout<<"Labels at bond "<<i<<": "<<std::endl;
    for(int aim=0;aim<dimInfo.locDimL(i);++aim){
      std::cout<<gqn.QNLabel(i-1,aim)<<"\t";
    }
    std::cout<<std::endl;
  }
}

//For adaption of D
void translateBlocks(baseTensor<std::complex<double> > const &source, siteQNOrderMatrix const &sourceTable, baseTensor<std::complex<double> > &target, siteQNOrderMatrix const &targetTable){
  //works only if sourceTable is contained in targetTable
  int const numBlocks=sourceTable.numBlocksLP();
  int lBlockSize, rBlockSize;
  std::vector<int> sourceIndices(3);
  std::vector<int> targetIndices(3);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=sourceTable.lBlockSizeLP(iBlock);
    rBlockSize=sourceTable.rBlockSizeLP(iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	sourceIndices[0]=sourceTable.siBlockIndexLP(iBlock,k);
	sourceIndices[1]=sourceTable.aiBlockIndexLP(iBlock,j);
	sourceIndices[2]=sourceTable.aimBlockIndexLP(iBlock,k);
	targetIndices[0]=targetTable.siBlockIndexLP(iBlock,k);
	targetIndices[1]=targetTable.aiBlockIndexLP(iBlock,j);
	targetIndices[2]=targetTable.aimBlockIndexLP(iBlock,k);
	target(targetIndices)=source(sourceIndices);
      }
    }
  }
}