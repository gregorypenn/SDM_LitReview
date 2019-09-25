#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#include "maximize.oxh"
#include "ReportMCMC.ox"
static decl Gall,GT1,Gt,M,N,NpredPA,NpredPO;

//PA: likelihood for LGCP
likelihood0(const p,const vf,const mC,const mXG,const vbeta,const vtheta,const vArea)
{
 decl loglike;
 loglike=-exp(vf+mXG*vbeta)'*vArea+mC'*(vf+mXG*vbeta); 	//normalized
 return {loglike};
}

NNGP_completion(const N0,const mC,const vInd)
{
 decl i;
 decl vF=zeros(N0,1);
 decl mB=zeros(N0,M);
 decl mCN,mICNN;
 for(i=0;i<N0;i++){
   if(i==0){
     vF[0]=mC[0][0];
   }else if(i<M){
     mCN=mC[i][:i-1];
     mICNN=invert(mC[:i-1][:i-1]);
     mB[i][:i-1]=mCN*mICNN;
     vF[i]=mC[i][i]-mCN*mICNN*mCN';
   }else{
     mCN=mC[i][vInd[i][]];
     mICNN=invert(mC[vInd[i][]'][vInd[i][]]);
     mB[i][]=mCN*mICNN;
     vF[i]=mC[i][i]-mCN*mICNN*mCN';
   }
 }
 return {mB,vF};
}

NNGP_completion_Inv(const N0,const mB,const vF,const vInd)
{
 decl i;
 decl mIA=unit(N0,N0);
 decl mDiag=zeros(N0,N0);
 for(i=0;i<N0;i++){
   mDiag[i][i]=1/vF[i];
   if(i==0){
   }else if(i<M){
	 mIA[i][:i-1]=-mB[i][:i-1];
   }else{
	 mIA[i][vInd[i][]]=-mB[i][];
   }
 }
 decl mInv=mIA'*mDiag*mIA;
 return {mInv};
}

//Sampling NNGP 
Sampling_W_NNGP(const N0,const vZ,const mX,const vWN,const vfPA,const vfPO,
const mB,const vF,const vInd,const vIndK,const vOrder,const vCt,const vgamma,const vtheta)
{
 decl i,t,j;
 decl vmu=zeros(N0,1);
 decl vvar=zeros(N0,1);
 decl vWNpos=vWN;
 decl sumF=zeros(N0,1);
 decl sumB=zeros(N0,1);
 decl vFp=vtheta[7]*vfPA+vtheta[8]*vfPO;
 decl vWN0;
 for(i=0;i<N0;i++){
   for(j=0;j<vCt[i];j++){
	   sumF[i]+=(mB[vIndK[i][j]][vOrder[i][j]]^2)/vF[vIndK[i][j]];
	   vWN0=vWNpos[vInd[vIndK[i][j]][]'];
	   sumB[i]+=(mB[vIndK[i][j]][vOrder[i][j]]*(vWNpos[vIndK[i][j]]-(mB[vIndK[i][j]][]*vWN0[]-mB[vIndK[i][j]][vOrder[i][j]]*vWN0[vOrder[i][j]])))/vF[vIndK[i][j]];
   }
       vvar[i]=1/(1/vtheta[6]+1/vF[i]+sumF[i]);
       vmu[i]=vvar[i]*((vZ[i]-mX[i][]*vgamma-vFp[i])/vtheta[6]+mB[i][]*vWNpos[vInd[i][]']/vF[i]+sumB[i]);
	   vWNpos[i]=vmu[i]+sqrt(vvar[i])*rann(1,1);
 }
 return {vWNpos};
}

Sampling_F_NNGP(const PA,const N0,const mCt,const mXG,const vfpri,const vfex,const mB,const vF,
const vInd,const vIndK,const vOrder,const vCt,const vZ,const vw,const mX,const vIndInv,
const vCInv,const vbeta,const vgamma,const vtheta,const vAreat)
{
 decl i,t,j;
 decl vmu=zeros(N0,1);
 decl vvar=zeros(N0,1);
 decl sumF=zeros(N0,1);
 decl sumB=zeros(N0,1);
 decl vfpos=zeros(N0,1);
 decl vlambdapri=mXG*vbeta+vfpri;
 decl vlambdapos=vlambdapri;
 decl vlambdapro=vlambdapri;
 decl like0=0;
 decl like1=0;
 decl d=0;
 decl vlam0,vlam1,alpha,u;
 for(i=0;i<N0;i++){
   if(i==0){
     vlambdapro[i]=mXG[i][]*vbeta+sqrt(vF[i])*rann(1,1);
   }else if(i<M){
     vlambdapro[i]=mXG[i][]*vbeta+mB[i][:i-1]*(vlambdapos[:i-1]-mXG[:i-1][]*vbeta)+sqrt(vF[i])*rann(1,1);
   }else{
     vlambdapro[i]=mXG[i][]*vbeta+mB[i][]*(vlambdapos[vInd[i][]']-mXG[vInd[i][]'][]*vbeta)+sqrt(vF[i])*rann(1,1);
   }
   like0=-exp(vlambdapos[i])*vAreat[i]+mCt[i]*vlambdapos[i];
   like1=-exp(vlambdapro[i])*vAreat[i]+mCt[i]*vlambdapro[i];
   if(PA==1){
     for(j=0;j<vCInv[i];j++){
	   like0+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[7]*(vlambdapos[i]-mXG[i][]*vbeta)-vtheta[8]*vfex[vIndInv[i][j]])^2;
	   like1+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[7]*(vlambdapro[i]-mXG[i][]*vbeta)-vtheta[8]*vfex[vIndInv[i][j]])^2;
     }
   }else{
     for(j=0;j<vCInv[i];j++){
	   like0+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[8]*(vlambdapos[i]-mXG[i][]*vbeta)-vtheta[7]*vfex[vIndInv[i][j]])^2;
	   like1+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[8]*(vlambdapro[i]-mXG[i][]*vbeta)-vtheta[7]*vfex[vIndInv[i][j]])^2;
     }
   }
   for(j=0;j<vCt[i];j++){
	 vlam0=vlambdapos[vInd[vIndK[i][j]][]']-mXG[vInd[vIndK[i][j]][]'][]*vbeta;
	 vlam1=vlambdapro[vInd[vIndK[i][j]][]']-mXG[vInd[vIndK[i][j]][]'][]*vbeta;
	 like0+=-((vlambdapos[vIndK[i][j]]-mXG[vIndK[i][j]][]*vbeta-mB[vIndK[i][j]][]*vlam0)^2)/(2*vF[vIndK[i][j]]);
     like1+=-((vlambdapro[vIndK[i][j]]-mXG[vIndK[i][j]][]*vbeta-mB[vIndK[i][j]][]*vlam1)^2)/(2*vF[vIndK[i][j]]);
   }
   alpha=min(1,exp(like1-like0));
   u=ranu(1,1);
   if(u<alpha){
	 vlambdapos[i]=vlambdapro[i];
	 d+=1;
   }else{
	 vlambdapro[i]=vlambdapri[i];
   }
   vfpos[i]=vlambdapos[i]-mXG[i][]*vbeta;
 }
 print("d=",d);
 return {vfpos};
}

Sampling_F_NNGP_PO(const PA,const N0,const vIndgrid,const mCt,const mXG,const vfpri,const vfex,const mB,const vF,
const vInd,const vIndK,const vOrder,const vCt,const vZ,const vw,const mX,const vIndInv,
const vCInv,const vbeta,const vgamma,const vtheta,const vAreat)
{
 decl i,t,j;
 decl vmu=zeros(N0,1);
 decl vvar=zeros(N0,1);
 decl sumF=zeros(N0,1);
 decl sumB=zeros(N0,1);
 decl vfpos=zeros(N0,1);
 decl vlambdapri=mXG*vbeta+vfpri;
 decl vlambdapos=vlambdapri;
 decl vlambdapro=vlambdapri;
 decl like0=0;
 decl like1=0;
 decl d=0;
 decl k=0;
 decl vlam0,vlam1,alpha,u;
 for(i=0;i<N0;i++){
   if(i==0){
     vlambdapro[i]=mXG[i][]*vbeta+sqrt(vF[i])*rann(1,1);
   }else if(i<M){
     vlambdapro[i]=mXG[i][]*vbeta+mB[i][:i-1]*(vlambdapos[:i-1]-mXG[:i-1][]*vbeta)+sqrt(vF[i])*rann(1,1);
   }else{
     vlambdapro[i]=mXG[i][]*vbeta+mB[i][]*(vlambdapos[vInd[i][]']-mXG[vInd[i][]'][]*vbeta)+sqrt(vF[i])*rann(1,1);
   }
   if(vIndgrid[i]==1){
     like0=-exp(vlambdapos[i])*vAreat[i]+mCt[i]*vlambdapos[i];
     like1=-exp(vlambdapro[i])*vAreat[i]+mCt[i]*vlambdapro[i];
	 k+=1;
   }else{
     like0=0;
     like1=0;
   }
   if(PA==1){
     for(j=0;j<vCInv[i];j++){
	   like0+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[7]*(vlambdapos[i]-mXG[i][]*vbeta)-vtheta[8]*vfex[vIndInv[i][j]])^2;
	   like1+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[7]*(vlambdapro[i]-mXG[i][]*vbeta)-vtheta[8]*vfex[vIndInv[i][j]])^2;
     }
   }else{
     for(j=0;j<vCInv[i];j++){
	   like0+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[8]*(vlambdapos[i]-mXG[i][]*vbeta)-vtheta[7]*vfex[vIndInv[i][j]])^2;
	   like1+=-(1/2)*(vZ[vIndInv[i][j]]-mX[vIndInv[i][j]][]*vgamma-vw[vIndInv[i][j]]-vtheta[8]*(vlambdapro[i]-mXG[i][]*vbeta)-vtheta[7]*vfex[vIndInv[i][j]])^2;
     }
   }
   for(j=0;j<vCt[i];j++){
	 vlam0=vlambdapos[vInd[vIndK[i][j]][]']-mXG[vInd[vIndK[i][j]][]'][]*vbeta;
	 vlam1=vlambdapro[vInd[vIndK[i][j]][]']-mXG[vInd[vIndK[i][j]][]'][]*vbeta;
	 like0+=-((vlambdapos[vIndK[i][j]]-mXG[vIndK[i][j]][]*vbeta-mB[vIndK[i][j]][]*vlam0)^2)/(2*vF[vIndK[i][j]]);
     like1+=-((vlambdapro[vIndK[i][j]]-mXG[vIndK[i][j]][]*vbeta-mB[vIndK[i][j]][]*vlam1)^2)/(2*vF[vIndK[i][j]]);
   }
   alpha=min(1,exp(like1-like0));
   u=ranu(1,1);
   if(u<alpha){
	 vlambdapos[i]=vlambdapro[i];
	 d+=1;
   }else{
	 vlambdapro[i]=vlambdapri[i];
   }
   vfpos[i]=vlambdapos[i]-mXG[i][]*vbeta;
 }
 print("d=",d);
 return {vfpos};
}


////Sampling hyperparameters 
Sampling_theta_NNGP(const ID,const N0,const mD,const mBpri,const vFpri,const vw,const vInd,const vthetapri)
{
  decl i,j;
  decl vthetapro=vthetapri;
  decl vthetapos=vthetapri;
  decl mBpro=zeros(N0,M);
  decl vFpro=zeros(N0,1);
  decl vMupro=zeros(N0,1);
  decl vMupri=zeros(N0,1);
  do{
  vthetapro[2*ID:2*ID+1]=vthetapri[2*ID:2*ID+1]+0.1*rann(2,1);
  }while(vthetapro[2*ID]<0 || vthetapro[2*ID+1]<0);
  decl mCovpro=vthetapro[2*ID]*exp(-vthetapro[2*ID+1]*mD);
  [mBpro,vFpro]=NNGP_completion(N0,mCovpro,vInd);
  for(i=0;i<N0;i++){
    if(i==0){
    }else if(i<M){
      vMupro[i]=mBpro[i][:i-1]*vw[:i-1];
      vMupri[i]=mBpri[i][:i-1]*vw[:i-1];
    }else{
      vMupro[i]=mBpro[i][]*vw[vInd[i][]'];
      vMupri[i]=mBpri[i][]*vw[vInd[i][]'];
    }
  }
  //NNGP completion
  decl loglike0=-(1/2)*sumc(log(vFpri))-sumc(((vw-vMupri).^2)./(2*vFpri));
  decl loglike1=-(1/2)*sumc(log(vFpro))-sumc(((vw-vMupro).^2)./(2*vFpro));
  decl u=ranu(1,1);
  decl alpha=min(1,exp(loglike1+log(densgamma(1/vthetapro[2*ID],2,0.1))-loglike0-log(densgamma(1/vthetapri[2*ID],2,0.1))));
  decl mBpos,vFpos,mLpos,mHpos,mInvCpos;
  if(u<alpha){
   vthetapos=vthetapro;
   mBpos=mBpro;
   vFpos=vFpro;
  }else{
   vthetapos=vthetapri;
   mBpos=mBpri;
   vFpos=vFpri;
  }
  return {vthetapos,mBpos,vFpos};
}

									 
Sampling_deltaPA(const vZ,const mX,const vgamma,const vw,const vIntPA,const vIntPO,const vtheta)
{
  decl i,l;
  decl var=1/(vIntPA'*vIntPA/vtheta[6]+1/100);
  decl sumMU=vIntPA'*(vZ-mX*vgamma-vw-vIntPO*vtheta[8]);
  decl vmu=var*((1/vtheta[6])*sumMU);
  decl deltaPApos;
    deltaPApos=vmu+sqrt(var)*rann(1,1);

  return {deltaPApos};
}

Sampling_deltaPO(const vZ,const mX,const vgamma,const vw,const vIntPA,const vIntPO,const vtheta)
{
  decl i,l;
  decl var=1/(vIntPO'*vIntPO/vtheta[6]+1/100);
  decl sumMU=vIntPO'*(vZ-mX*vgamma-vw-vIntPA*vtheta[7]);
  decl vmu=var*((1/vtheta[6])*sumMU);
  decl deltaPOpos;
  decl a=(0-vmu)/sqrt(var);
  decl u=probn(a)+(1-probn(a))*ranu(1,1);
  decl x=quann(u);
  deltaPOpos=vmu+sqrt(var)*x;
  return {deltaPOpos};
}

Sampling_latent(const vY,const mX,const vgamma,const vw,const vIntPA,const vIntPO,const vtheta)
{
  decl i,l;
  decl vZ=zeros(N,1);
  decl a,x,u;
  decl mCov=mX*vgamma+vw+vtheta[7]*vIntPA+vtheta[8]*vIntPO;
  for(i=0;i<N;i++){
    	a=(0-mCov[i])/sqrt(vtheta[6]);
	  if(vY[i]==1){
		u=probn(a)+(1-probn(a))*ranu(1,1);
	  }else{
		u=probn(a)*ranu(1,1);
	  }
		x=quann(u);
		vZ[i]=mCov[i]+sqrt(vtheta[6])*x;
  }
  return {vZ};
}

Sampling_gammaPA(const vZ,const mX,const vw,const vIntPA,const vIntPO,const vtheta)
{
  decl i,l;
  decl p=columns(mX);
  decl sumX=mX'*mX;
  decl sumMU=mX'*(vZ-vw-vtheta[7]*vIntPA-vtheta[8]*vIntPO);
  decl mSig=invert(sumX/vtheta[6]+(1/100)*unit(p,p));
  decl vmu=mSig*(1/vtheta[6])*sumMU;
  decl Chol=choleski(mSig);
  decl vgamma=vmu+Chol*rann(p,1);
  return {vgamma};
}

//Sampling covariates coefficients
Sampling_beta(const p,const vf,const mCt,const mXGt,const vbetapri,const vtheta,const vArea)
{
  decl i,j;
  decl vbetapro=vbetapri;
  decl vbetapos=vbetapri;
  decl likepri,likepro;
  decl u=ranu(1,1);
  decl alpha=0;
  //generating vars via Gibbs
  for(i=0;i<(p+1);i++){
  vbetapro[i]=vbetapri[i]+0.1*rann(1,1);
  [likepri]=likelihood0(p,vf,mCt,mXGt,vbetapos,vtheta,vArea);
  [likepro]=likelihood0(p,vf,mCt,mXGt,vbetapro,vtheta,vArea);
  u=ranu(1,1);
  alpha=min(1,exp(likepro-((vbetapro[i]^2)/200)-likepri+((vbetapri[i]^2)/200)));
  if(u<alpha){
   vbetapos[i]=vbetapro[i];
  }else{
   vbetapos[i]=vbetapri[i];
   vbetapro[i]=vbetapri[i];
  }
  }
  return {vbetapos};
}

Predict(const Imodel,const p,const mXpredPA,const vw,const vIntPA,const vIntPO,const vIndpred,const vFpredPA,const vetapred,const vgamma,const vtheta)
{
  decl i;
  decl vZpred=zeros(NpredPA,1);
  //Model(c)
  if(Imodel==0){
    vZpred=mXpredPA*vgamma;
  }
  if(Imodel==1){
    for(i=0;i<NpredPA;i++){
      vZpred[i]=mXpredPA[i][]*vgamma+vFpredPA[i][]*vw[vIndpred[i][]']+vetapred[i];
    }
  }	
  if(Imodel==2){
    for(i=0;i<NpredPA;i++){
      vZpred[i]=mXpredPA[i][]*vgamma+vFpredPA[i][]*vw[vIndpred[i][]']+vetapred[i]+vIntPA[i]*vtheta[7];
    }
  }	
  if(Imodel==3){
    vZpred=mXpredPA*vgamma+vIntPA*vtheta[7];
  }
  if(Imodel==4){
    for(i=0;i<NpredPA;i++){
      vZpred[i]=mXpredPA[i][]*vgamma+vFpredPA[i][]*vw[vIndpred[i][]']+vetapred[i]+vIntPO[i]*vtheta[8];
    }
  }	
  if(Imodel==5){
    vZpred=mXpredPA*vgamma+vIntPO*vtheta[8];
  }
  if(Imodel==6){
    vZpred=mXpredPA*vgamma+vIntPA*vtheta[7]+vIntPO*vtheta[8];
  }
  if(Imodel==7){
    for(i=0;i<NpredPA;i++){
      vZpred[i]=mXpredPA[i][]*vgamma+vFpredPA[i][]*vw[vIndpred[i][]']+vetapred[i]+vIntPA[i]*vtheta[7]+vIntPO[i]*vtheta[8];
    }
  }	
  return {vZpred};
}


main()
{
//Real Data
 decl i,j,l;
 decl file1=loadmat("Covariate_POGrid.csv");	//covariate surface
 decl file2=loadmat("NewEnglandPA.csv");		//
 decl file3=loadmat("GridNE80.csv");
 decl file4=loadmat("GridNE320.csv");			//Grid for 
 decl vIDall=file4[][13];   //All
 decl vIDT1=file4[][14];	//Grid where T(s)=1
 decl Ipred=0; //0:inference, 1:prediction
 decl dprob=1;	 //Proportion of training samples
 decl Imodel=0; //0:(a), 1:(b), 2:(c), 3:(d), 4:(c'), 5:(d'), 6:(e), 7:(f)
 print("Imodel=",Imodel);
 Gt=rows(file3);
 Gall=sumc(vIDall);
 GT1=sumc(vIDT1);
 print("Gall=",Gall);
 print("GT1=",GT1);
 decl dGC=rows(file1);//number of total grids for covariate surface
 decl scale=10000;	 //scale for UTG
 decl p=7;			 //dimension of covariates
 decl dsp=0;		 //Index for species: 0:MR, 1:OB, 2:JB, 3:GB, 4:AO, 5:BB, 6:GM 
 decl mCPA=file3[][4];
 decl mCPO=zeros(Gall,1);
 decl vAreaPO=zeros(Gall,1);
 decl vAreaPA=file3[][3]/sumc(file3[][3]);
 decl mGT1=zeros(GT1,2);
 decl mGall=zeros(Gall,2);
 decl vIPA=ones(Gt,1);
 decl vIPO=ones(Gall,1);
 decl l1=0;
 decl l2=0;
 for(i=0;i<rows(file4);i++){
   if(vIDT1[i]==1){
	 mGT1[l1][]=file4[i][:1]./scale;
	 l1+=1;
   }	 
   if(vIDall[i]==1){
	 mGall[l2][]=file4[i][:1]./scale;
	 mCPO[l2]=file4[i][6+dsp];
	 vAreaPO[l2]=file4[i][3];
	 if(vIDT1[i]==0){
	   vIPO[l2]=0;
	   mCPO[l2]=0;
	   vAreaPO[l2]=0;
	 }
	 l2+=1;
   }	 
 }
 vAreaPO=vAreaPO/sumc(vAreaPO);
 print("vAreaPO=",vAreaPO);
 print("CPO=",sumc(mCPO));
 print("CPA=",sumc(mCPA));
 N=rows(file2);
 print("dsp=",dsp);
 print("N=",N);
 ///////////////
 //PA analysis//
 ///////////////
 decl mean1,var1;
 for(i=0;i<p;i++){
   mean1=meanc(file1[][4+i]);
   var1=varc(file1[][4+i]);
   file1[][4+i]=(file1[][4+i]-ones(dGC,1)*mean1)/sqrt(var1);
   file2[][4+7+i]=(file2[][4+7+i]-ones(N,1)*mean1)/sqrt(var1);
 }
 print("min2=",minc(file1[][4+1]));
 print("max5=",maxc(file1[][4+4]));
 decl d=0;
 decl mGs=zeros(N,2);	//locations for PA
 decl mGspred=zeros(N,2);	//locations for PA
 decl mG0=file1[][2:3]./scale;	//covariate grids
 decl mGt=file3[][:1]./scale;	//covariate grids
 decl vY=zeros(N,1);	
 decl vYpred=zeros(N,1);
 decl mX=ones(N,p+1);
 decl mXGt=ones(Gt,p+1);
 decl mXGall=ones(Gall,p+1);
 decl mXpred=ones(N,p+1);
 decl vIndreserve=zeros(N,1);
 decl u,I;
 d=0;
 NpredPA=0;
 ranseed(10);
 for(i=0;i<N;i++){
   u=ranu(1,1);
   if(u<dprob){
     vY[i-NpredPA]=file2[i][4+dsp];
	 mGs[i-NpredPA][]=file2[i][2:3]./scale;	
     mX[i-NpredPA][1:]=file2[i][4+7:4+7+p-1];
	 if(vY[i-NpredPA]==1){
     }else{
	   d+=1;
     }
   }else{
     if(Ipred==1){
       vYpred[NpredPA]=file2[i][4+dsp];
	   mGspred[NpredPA][]=file2[i][2:3]./scale;	
	   mXpred[NpredPA][1:]=file2[i][4+7:4+7+p-1];
	   vIndreserve[NpredPA]=i;
	   NpredPA+=1;
	 }
   }
 }
// print("vYpred=",vYpred[:NpredPA-1]);
// print("mGspred=",mGspred[:NpredPA-1][]);
 N=N-NpredPA;
 print("N=",N);
 mGs=mGs[:N-1][];
 vY=vY[:N-1][];
 mX=mX[:N-1][];
 decl vZ=zeros(N,1);
 decl mDist,mDists,mDistp;
 ranseed(10);
 for(i=0;i<Gt;i++){
   mDist=sqrt(fabs(mG0[][0]-ones(dGC,1)*mGt[i][0]).^2+fabs(mG0[][1]-ones(dGC,1)*mGt[i][1]).^2);
   I=sortcindex(mDist);
   mXGt[i][1:]=file1[I[0]][4:4+p-1];
 }
 for(i=0;i<Gall;i++){
   mDist=sqrt(fabs(mG0[][0]-ones(dGC,1)*mGall[i][0]).^2+fabs(mG0[][1]-ones(dGC,1)*mGall[i][1]).^2);
   I=sortcindex(mDist);
   mXGall[i][1:]=file1[I[0]][4:4+p-1];
 }
 print("min2=",minc(mXGt[][1+1]));
 print("max5=",maxc(mXGt[][1+4]));
 decl vInd=zeros(N,1);
 decl vIndPO=zeros(N,1);
 decl vIndG=zeros(N,1);
 decl minD;
 decl d0=0;
 decl d1=0;
 for(i=0;i<N;i++){
   mDist=sqrt(fabs(mGt[][0]-ones(Gt,1)*mGs[i][0]).^2+fabs(mGt[][1]-ones(Gt,1)*mGs[i][1]).^2);
   I=sortcindex(mDist);
   vInd[i]=I[0];
   mDist=sqrt(fabs(mGall[][0]-ones(Gall,1)*mGs[i][0]).^2+fabs(mGall[][1]-ones(Gall,1)*mGs[i][1]).^2);
   I=sortcindex(mDist);
   vIndPO[i]=I[0];
 }
 decl vIndInv=zeros(Gt,300);
 decl vCInv=zeros(Gt,1);
 for(i=0;i<N;i++){
   vIndInv[vInd[i]][vCInv[vInd[i]]]=i;
   vCInv[vInd[i]]+=1;
 }
 decl vIndInvPO=zeros(Gall,300);
 decl vCInvPO=zeros(Gall,1);
 for(i=0;i<N;i++){
   vIndInvPO[vIndPO[i]][vCInvPO[vIndPO[i]]]=i;
   vCInvPO[vIndPO[i]]+=1;
 }

 
 decl mCorr=zeros(p,p);
 decl mCov=0;
 decl mD=zeros(Gt,2);
 for(i=0;i<p;i++){
   for(j=0;j<i;j++){
     mD[][0]=mXGt[][1+i];
     mD[][1]=mXGt[][1+j];
     mCov=variance(mD);
	 mCorr[i][j]=mCov[0][1]/(sqrt(varc(mXGt[][1+i]))*sqrt(varc(mXGt[][1+j])));
   }
 }
 print("mCorr=",mCorr);
 decl q=p;
 decl q2=p;

 decl mDft=sqrt(fabs(mGt[][0]-mGt[][0]').^2+fabs(mGt[][1]-mGt[][1]').^2);
 decl mDall=sqrt(fabs(mGall[][0]-mGall[][0]').^2+fabs(mGall[][1]-mGall[][1]').^2);
 decl mDw=sqrt(fabs(mGs[][0]-mGs[][0]').^2+fabs(mGs[][1]-mGs[][1]').^2);
 print("maxDft=",maxr(maxc(mDft)));
 print("phiu=",-log(0.01)/minc(minr(mDft)));
 print("phil=",-log(0.05)/maxc(maxr(mDft)));
 print("maxDw=",maxr(maxc(mDw)));
 print("phiu=",-log(0.01)/minc(minr(mDw)));
 print("phil=",-log(0.05)/maxc(maxr(mDw)));

 ////////////////
 //NNGP indices//
 ////////////////
 M=15;
 decl vIndw=zeros(N,M);	   //Neighbor indices
 decl vIndf=zeros(Gt,M);	   //Neighbor indices
 decl vIndfPO=zeros(Gall,M);	   //Neighbor indices
 decl vI;
 for(i=0;i<N;i++){
   if(i==0){
	 vIndw[i][0]=0;
   }else if(i<M){
	 for(j=0;j<i;j++){
	   vIndw[i][j]=j;
	 }
   }else{
     vI=sortcindex(mDw[i][:i-1]');
	 vIndw[i][]=vI[:M-1]';
   }
 }
 for(i=0;i<Gt;i++){
   if(i==0){
	 vIndf[i][0]=0;
   }else if(i<M){
	 for(j=0;j<i;j++){
	   vIndf[i][j]=j;
	 }
   }else{
     vI=sortcindex(mDft[i][:i-1]');
	 vIndf[i][]=vI[:M-1]';
   }
 }
 for(i=0;i<Gall;i++){
   if(i==0){
	 vIndfPO[i][0]=0;
   }else if(i<M){
	 for(j=0;j<i;j++){
	   vIndfPO[i][j]=j;
	 }
   }else{
     vI=sortcindex(mDall[i][:i-1]');
	 vIndfPO[i][]=vI[:M-1]';
   }
 }
 decl vIndKw=zeros(N,10*M);
 decl vCtw=zeros(N,1);
 decl vOrderw=zeros(N,10*M);
 decl vIndKf=zeros(Gt,10*M);
 decl vCtf=zeros(Gt,1);
 decl vOrderf=zeros(Gt,10*M);
 decl vIndKfPO=zeros(Gall,10*M);
 decl vCtfPO=zeros(Gall,1);
 decl vOrderfPO=zeros(Gall,10*M);
 for(i=0;i<N;i++){
   if(i==0){
   }else if(i<M){
   for(j=0;j<i;j++){
	 vIndKw[vIndw[i][j]][vCtw[vIndw[i][j]]]=i;
	 vOrderw[vIndw[i][j]][vCtw[vIndw[i][j]]]=j;
	 vCtw[vIndw[i][j]]+=1;
   }
   }else{
   for(j=0;j<M;j++){
	 vIndKw[vIndw[i][j]][vCtw[vIndw[i][j]]]=i;
	 vOrderw[vIndw[i][j]][vCtw[vIndw[i][j]]]=j;
	 vCtw[vIndw[i][j]]+=1;
   }
   }
 }
 for(i=0;i<Gt;i++){
   if(i==0){
   }else if(i<M){
   for(j=0;j<i;j++){
	 vIndKf[vIndf[i][j]][vCtf[vIndf[i][j]]]=i;
	 vOrderf[vIndf[i][j]][vCtf[vIndf[i][j]]]=j;
	 vCtf[vIndf[i][j]]+=1;
   }
   }else{
   for(j=0;j<M;j++){
	 vIndKf[vIndf[i][j]][vCtf[vIndf[i][j]]]=i;
	 vOrderf[vIndf[i][j]][vCtf[vIndf[i][j]]]=j;
	 vCtf[vIndf[i][j]]+=1;
   }
   }
 }
 for(i=0;i<Gall;i++){
   if(i==0){
   }else if(i<M){
   for(j=0;j<i;j++){
	 vIndKfPO[vIndfPO[i][j]][vCtfPO[vIndfPO[i][j]]]=i;
	 vOrderfPO[vIndfPO[i][j]][vCtfPO[vIndfPO[i][j]]]=j;
	 vCtfPO[vIndfPO[i][j]]+=1;
   }
   }else{
   for(j=0;j<M;j++){
	 vIndKfPO[vIndfPO[i][j]][vCtfPO[vIndfPO[i][j]]]=i;
	 vOrderfPO[vIndfPO[i][j]][vCtfPO[vIndfPO[i][j]]]=j;
	 vCtfPO[vIndfPO[i][j]]+=1;
   }
   }
 }
 
 ///////////////
 //MCMC tuning//
 ///////////////
 decl istart=20000; //number of burn-in period
 decl iLoop=40000; //total number of MCMC iteration
 decl thin=(iLoop-istart)/100;   //thinning for posterior samples
 decl betaPA_result=new matrix[(iLoop-istart)][q2+1];	
 decl betaPO_result=new matrix[(iLoop-istart)][q2+1];	
 decl gamma_result=new matrix[(iLoop-istart)][q+1];	
 decl theta_result=new matrix[(iLoop-istart)][9];
 decl theta2_result=new matrix[(iLoop-istart)][3];
 decl vZ_result=new matrix[(iLoop-istart)/thin][N];
 decl vW_result=new matrix[(iLoop-istart)/thin][N];
 decl vfPA_result=new matrix[(iLoop-istart)/thin][Gt];
 decl vIntPA_result=new matrix[(iLoop-istart)/thin][Gt];
 decl vfPO_result=new matrix[(iLoop-istart)/thin][Gt];
 decl vIntPO_result=new matrix[(iLoop-istart)/thin][Gt];
 decl vZpre_result=new matrix[(iLoop-istart)/thin][NpredPA];
 decl icount1=0;			  
 decl vbetapos=zeros(q+1,2);
 vbetapos[0][0]=log(sumc(mCPA)/sumc(vAreaPA))/1.5;
 vbetapos[0][1]=log(sumc(mCPO)/sumc(vAreaPO))/1.5;
 print("beta0=",vbetapos[0][]);
 decl vgammapos=zeros(q+1,1);
 decl vthetapos=ones(9,1);
 vthetapos[0]=3;   //thetaPA
 vthetapos[1]=2;   //thetaPA
 vthetapos[2]=3;   //thetaPO
 vthetapos[3]=2;   //thetaPO
 vthetapos[4]=1;   //thetaw
 vthetapos[5]=1;   //thetaw
 vthetapos[6]=1;   //tau, fixed
 vthetapos[7]=0;   //deltaPA
 vthetapos[8]=0;   //deltaPO
 print("phi=",vthetapos[1]);
 decl deltapos=0;
 print("deltapos=",deltapos);
 decl time = timer();
 decl nuwpos=zeros(N,1);
 if(Imodel==1 || Imodel==2 || Imodel==4 || Imodel==7){
   nuwpos=rann(N,1);
 }else{
   nuwpos=zeros(N,1);
 }
 print("nuwpos=",nuwpos[:5]);
 decl mCw=vthetapos[4]*exp(-vthetapos[5]*mDw);
 decl mLwpos=choleski(mCw);
 decl vwpos=mLwpos*nuwpos;
 decl phiU=-log(0.01)/minc(minr(mDw));
 decl phiL=-log(0.05)/maxc(maxr(mDw));
 decl nutPApos=rann(Gt,1);
 decl nutPOpos=rann(Gall,1);
 decl mCfPA=vthetapos[0]*exp(-vthetapos[1]*mDft);
 decl mCfPO=vthetapos[2]*exp(-vthetapos[3]*mDall);
 decl mLtPApos=choleski(mCfPA);
 decl mLtPOpos=choleski(mCfPO);
 decl vftPApos=mLtPApos*nutPApos;
 decl vftPOpos=mLtPOpos*nutPOpos;
 decl vFwpos,mBwpos,vFfPApos,mBfPApos,vFfPOpos,mBfPOpos;
 //NNGP completion
 [mBwpos,vFwpos]=NNGP_completion(N,mCw,vIndw);
 [mBfPApos,vFfPApos]=NNGP_completion(Gt,mCfPA,vIndf);
 [mBfPOpos,vFfPOpos]=NNGP_completion(Gall,mCfPO,vIndfPO);
 
 decl icount2=0;
 decl likepos;
 decl likelogit;

 //////////////
 //Prediction//
 //////////////
 decl mBpred,vFpred,mCovpred,vetapred;
 decl mICwpos,mLpredPA,vFpredPA;
 decl vZpre;
 decl vIndpredPA=zeros(NpredPA,1);
 decl vIndpredPO=zeros(NpredPA,1);
 decl vIndpred=zeros(NpredPA,M);	   //Neighbor indices
 if(Ipred==1){
   vYpred=vYpred[:NpredPA-1][];
   mGspred=mGspred[:NpredPA-1][];
   mXpred=mXpred[:NpredPA-1][];
   for(i=0;i<NpredPA;i++){
     mDistp=sqrt(fabs(mGt[][0]-ones(Gt,1)*mGspred[i][0]).^2+fabs(mGt[][1]-ones(Gt,1)*mGspred[i][1]).^2);
     I=sortcindex(mDistp);
     vIndpredPA[i]=I[0];
     mDistp=sqrt(fabs(mGall[][0]-ones(Gall,1)*mGspred[i][0]).^2+fabs(mGall[][1]-ones(Gall,1)*mGspred[i][1]).^2);
     I=sortcindex(mDistp);
     vIndpredPO[i]=I[0];
   }
   decl mDspred=sqrt(fabs(mGspred[][0]-mGs[][0]').^2+fabs(mGspred[][1]-mGs[][1]').^2);
   for(i=0;i<NpredPA;i++){
     vI=sortcindex(mDspred[i][]');
     vIndpred[i][]=vI[:M-1]';
   }
 }
 decl mDs1=sqrt(fabs(mGspred[][0]-mGspred[][0]').^2+fabs(mGspred[][1]-mGspred[][1]').^2);
 decl mDs2=sqrt(fabs(mGspred[][0]-mGs[][0]').^2+fabs(mGspred[][1]-mGs[][1]').^2);
 decl vZpremean=zeros(NpredPA,1);

 /////////////////////////////////
 //MCMC iteration for fitting LGCP for origin point patterns
 ranseed(10);
 time = timer();
 d=0;
 d1=0;
 for(i=0,icount1=0;i<iLoop;i++)				   
 {
  //point process models
  if(Imodel==2 || Imodel==3 || Imodel==6 || Imodel==7){
    [vftPApos]=Sampling_F_NNGP(1,Gt,mCPA,mXGt,vftPApos,vftPOpos[vIndPO],mBfPApos,vFfPApos,vIndf,vIndKf,vOrderf,vCtf,vZ,vwpos,mX,vIndInv,vCInv,vbetapos[][0],vgammapos,vthetapos,vAreaPA);
    [vbetapos[][0]]=Sampling_beta(p,vftPApos,mCPA,mXGt,vbetapos[][0],vthetapos,vAreaPA);
    [vthetapos[7]]=Sampling_deltaPA(vZ,mX,vgammapos,vwpos,vftPApos[vInd],vftPOpos[vIndPO],vthetapos);
    [vthetapos,mBfPApos,vFfPApos]=Sampling_theta_NNGP(0,Gt,mDft,mBfPApos,vFfPApos,vftPApos,vIndf,vthetapos);
  }
  if(Imodel==4 || Imodel==5 || Imodel==6 || Imodel==7){
    [vftPOpos]=Sampling_F_NNGP_PO(0,Gall,vIPO,mCPO,mXGall,vftPOpos,vftPApos[vInd],mBfPOpos,vFfPOpos,vIndfPO,vIndKfPO,vOrderfPO,vCtfPO,vZ,vwpos,mX,vIndInvPO,vCInvPO,vbetapos[][1],vgammapos,vthetapos,vAreaPO);
    [vbetapos[][1]]=Sampling_beta(p,vftPOpos,mCPO,mXGall,vbetapos[][1],vthetapos,vAreaPO);
    [vthetapos[8]]=Sampling_deltaPO(vZ,mX,vgammapos,vwpos,vftPApos[vInd],vftPOpos[vIndPO],vthetapos);
    [vthetapos,mBfPOpos,vFfPOpos]=Sampling_theta_NNGP(1,Gall,mDall,mBfPOpos,vFfPOpos,vftPOpos,vIndfPO,vthetapos);
  }
  //PA sampling
  [vZ]=Sampling_latent(vY,mX,vgammapos,vwpos,vftPApos[vInd],vftPOpos[vIndPO],vthetapos);
  [vgammapos]=Sampling_gammaPA(vZ,mX,vwpos,vftPApos[vInd],vftPOpos[vIndPO],vthetapos);
  if(Imodel==1 || Imodel==2 || Imodel==4 || Imodel==7){
    [vwpos]=Sampling_W_NNGP(N,vZ,mX,vwpos,vftPApos[vInd],vftPOpos[vIndPO],mBwpos,vFwpos,vIndw,vIndKw,vOrderw,vCtw,vgammapos,vthetapos);
    [vthetapos,mBwpos,vFwpos]=Sampling_theta_NNGP(2,N,mDw,mBwpos,vFwpos,vwpos,vIndw,vthetapos);
  }
  print("i=",i);
  print("vbeta=",vbetapos');
  print("vgamma=",vgammapos');
  print("vtheta=",vthetapos');
  if(i>=istart){
    betaPA_result[i-istart][:p]=vbetapos[][0]';
    betaPO_result[i-istart][:p]=vbetapos[][1]';
    gamma_result[i-istart][:p]=vgammapos';
    theta_result[i-istart][:8]=vthetapos';
    theta2_result[i-istart][0]=vthetapos[0]*vthetapos[1];
    theta2_result[i-istart][1]=vthetapos[2]*vthetapos[3];
    theta2_result[i-istart][2]=vthetapos[4]*vthetapos[5];
  }
  if(Ipred==1){
    if(i>=istart && fmod(i-istart,thin)==0){
      vZ_result[d][]=vZ';
      vW_result[d][]=vwpos';
      vfPA_result[d][]=vftPApos';
      vIntPA_result[d][]=(vftPApos+mXGt*vbetapos[][0])';
      vfPO_result[d][]=vftPOpos';
      vIntPO_result[d][]=(vftPOpos+mXGall*vbetapos[][1])';
      vFpredPA=zeros(NpredPA,M);
      vetapred=zeros(NpredPA,1);
      for(j=0;j<NpredPA;j++){
	    mICwpos=invert(vthetapos[4]*exp(-vthetapos[5]*mDw[vIndpred[j][]'][vIndpred[j][]]));
        vFpredPA[j][]=vthetapos[4]*exp(-vthetapos[5]*mDs2[j][vIndpred[j][]])*mICwpos;
        mCovpred=vthetapos[4]-vFpredPA[j][]*(vthetapos[4]*exp(-vthetapos[5]*mDs2[j][vIndpred[j][]])');
	    vetapred[j]=sqrt(mCovpred)*rann(1,1);
      }
      [vZpre]=Predict(Imodel,q,mXpred,vwpos,vftPApos[vIndpredPA],vftPOpos[vIndpredPO],vIndpred,vFpredPA,vetapred,vgammapos,vthetapos);
      vZpre_result[d][]=vZpre';
      d+=1;
    }
  }
 }
 println( "Execution time2:   ", "%62s", timespan( time ) );
 println( "Program finished2: ", "%62s", date() );
 if(Ipred==1){
   savemat("vZ_Fusion3d.csv",vZ_result');
   savemat("vW_Fusion0f.csv",vW_result');
   savemat("vfPA_Fusion0f.csv",vfPA_result');
   savemat("vIntPA_Fusion0f.csv",vIntPA_result');
   savemat("vfPO_Fusion0f.csv",vfPO_result');
   savemat("vIntPO_Fusion0f.csv",vIntPO_result');
   savemat("NEFusion2PredmY6f.csv",vZpre_result');
 }


 decl fileC;
 fileC=fopen("PAPO","a");
 fclose(fileC);

 decl outC=new ReportMCMC(gamma_result[][]);
 outC.SetBandwidth((iLoop-istart)/10);	  // change the bandwidth.   Default= number of sample x 0.1
 outC.SetFormat(500);	      // change the output format. Default="%#13.5g" // change the figure type.   Default= eps
 outC.SetOutfileName("PAPO"); // change outputfile name. default=Output
 outC.SetVarNames({"gamma0","gamma1","gamma2","gamma3","gamma4","gamma5","gamma6","gamma7"}); // Set variable names. Default=Param1,..
 outC.Report(); 	


}