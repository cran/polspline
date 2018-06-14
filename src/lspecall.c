/*
*  Copyright (C) 1993--2018  Charles Kooperberg
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  The text of the GNU General Public License, version 2, is available
*  as http://www.gnu.org/copyleft or by writing to the Free Software
*  Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
#include <math.h>
#include <stdio.h>
#include "R.h"
#define Salloc(n, t)  (t *)R_alloc((long)(n), (int)sizeof(t))
#define NBMAX 70
#define PIL 3.141592653589793116
#define TINY 1.0e-20
#define DIM5 NBMAX+5

void F77_NAME(xdsifa)(double[][DIM5], int *, int *, int *, int *);
void F77_NAME(xdsisl)(double[][DIM5], int *, int *, int *, double *);
void F77_NAME(xdsidi)(double[][DIM5], int *, int *, int *, double *, int *, double *, int *);

static void tslusolve(),tsluinverse(),tsintsum(),tsallocer(),tspsps2(),tsbasis();
static void tsb1(),tsb2(),tsb3(),tsb4(),tsb5();
static double *tssdvec(),**tssdmat(),***tssdtri(),**uumm,*uuaa,*uuww,*uuvv1,*uuvv2;
static double *uubetan,**xxcumul,tsraod(),tsraoc(),tsnew(),tslogall();
static int *tssivec(),*uuika,silent,tsadd(),tsrem();
/******************************************************************************/
/* This function controls the updown movements */
/******************************************************************************/

void tspspsx(dims)
int *dims;
{
      dims[0]=NBMAX;
      return;
}

void tspsps(dims,data,knots,atoms,alpha,logs,theta,ad,mass)
	int *dims,*atoms,*ad;
	double *data,*knots,*alpha,*logs,*theta,*mass;

/* for most variables see tspsps2 below */
{

/* dims:
   0 - nx (SAMPLE SIZE) 
   1 - maxdim 
   2 - dimatt (should we attain maxdim?) 
   3 - maxknots 
   4 - ktsatt 
   5 - nknots 
   6 - maxatoms 
   7 - spkatt 
   8 - natoms
   9 - odd
   10 - repeat
   11 - error 
   12 - mind
   yy - winning log-likelihood, last time
   i,j - counter
   xspk - original atoms
   er2 - number of problems
   cdims - copy of dims
   xkts - original kts
   xx - x-coordinates
   zz - log(data)
   mass - minimum mass in a atom
   zzz - utility */

   int i,j,er2=0,*xspk,*nothere;
   int cdims[13];
   double yy=0.,*xx,*zz,*xkts,zzz=0.;

/* we only want the length of NBMAX */
   if(dims[0]<0){
      dims[0]=NBMAX;
      return;
   }
   tsallocer();

/* allocation */
   xx=tssdvec(dims[0]+1);
   zz=tssdvec(dims[0]);
   nothere=tssivec(dims[0]);
   xkts=tssdvec(dims[0]);
   xspk=tssivec(dims[0]);
   silent=1-dims[11];
   dims[11]=0;

/* at least one repeat */
   if(dims[10]<1)dims[10]=1;

/* initialize the x-values - this makes a difference whether the original
   series was odd or even */
   if(dims[9]==0){
      for(i=0;i<dims[0];i++)xx[i]=PIL*(double)i/(double)(dims[0]-1);
   }
   else {
      for(i=0;i<dims[0];i++)xx[i]=PIL*(double)i/(double)(dims[0]-0.5);
   }
   for(i=0;i<dims[0];i++)nothere[i]=0;

/* the log of the data */
   dims[12]=1;
   for(i=0;i<dims[0];i++)zz[i]=log(data[i]);

/* and up and down and up and down and up and down */
   for(i=0;i<13;i++)cdims[i]=dims[i];
   xxcumul=tssdmat(11+NBMAX*4,dims[0]+1);
   for(j=0;j<=dims[10];j++){
      

/* save what we might need later */
       
      for(i=0;i<dims[0];i++)xkts[i]=knots[i];
      for(i=0;i<dims[0];i++)xspk[i]=atoms[i];
      for(i=2;i<13;i++)cdims[i]=dims[i];
      do{
/* recover */
         for(i=0;i<dims[0];i++)knots[i]=xkts[i];
         for(i=0;i<dims[0];i++)atoms[i]=xspk[i];
         for(i=0;i<13;i++)dims[i]=cdims[i];
         if(j==0){
            if(dims[10]>1 && dims[7]==0){
               dims[6]=0;           
               dims[7]=1;
            }   
            else j=1;
         }
         tspsps2(dims,data,knots,*alpha,logs,theta,ad,xx,zz,atoms,nothere,
                                                          *mass);

/* if we did not converge, initialize again */
         if(dims[11]!=0){
            er2++;
            if(er2==1){
               for(i=0;i<dims[0];i++)zzz+=data[i];
               for(i=0;i<dims[0];i++)data[i]+=0.000001*zzz/((double)dims[0]);
               for(i=0;i<dims[0];i++)zz[i]=log(data[i]);
            }
            else {
               dims[12]+=1;
               cdims[12]+=1;
            }
         }
         if(j>0 && dims[5]+dims[8]>cdims[1]-5 && dims[2]==0){
            cdims[1]=dims[5]+dims[8]+5;
            if(cdims[1]>NBMAX-5)cdims[1]=NBMAX-5;
         } 
         if(j==0 && dims[5]+dims[8]>cdims[1]-1 && dims[2]==0){
            cdims[1]=dims[5]+dims[8]+1;
            if(cdims[1]>NBMAX-5)cdims[1]=NBMAX-1;
         } 

/* we really didn't converge */
      }while(dims[11]!=0 && dims[12]<6);
      if(j==0){
         dims[7]= 0; 
         dims[6]= -1; 
      }

/* return on error or on only repeat */
      if(dims[10]==1)dims[10]=0;
      if(dims[10]==0 || dims[11]!=0) return;
      for(i=3;i<13;i++)cdims[i]=dims[i];

/* there was no improvement, or the winner was fitted during addition */
      if(j>0){
         if(fabs((double)(yy-logs[dims[5]+dims[8]]))<TINY 
             ||ad[dims[5]+dims[8]]==0)
         {
            if(dims[5]+dims[8]!=dims[1]){
               dims[10]= -j;
            }
            else{
               yy=logs[dims[5]+dims[8]];
            }
         }
         else{
            yy=logs[dims[5]+dims[8]];
         }
      }
   }

/* go once more down */
   if(dims[10]<0 && dims[5]+dims[8]>1){
/* maximum number is present number */
      dims[1]=dims[5]+dims[8];
/* attain nothing */
      dims[2]=0;
      dims[4]=0;
      dims[7]=0;
      tspsps2(dims,data,knots,*alpha,logs,theta,ad,xx,zz,atoms,nothere,*mass);
   }

   dims[10]= -dims[10];
}
/******************************************************************************/
/* does the work */
/******************************************************************************/
static void tspsps2(dims,data,knots,alpha,logs,theta,ad,xx,zz,atoms,nothere,mass)
int *dims,*atoms,*ad,*nothere;
double *data,*knots,alpha,*logs,*xx,*zz,*theta,mass;
{
/* dims - various integer parameters
   data - periodogram
   knots - starting knots/best knots
   alpha - penalty parameter (bic)
   logs - log-likelihood of fitted models
   theta - coefficients (in powerbasis format)
   ad - was a model fit during the addition (0), or deletion (1) stage
   zz - log(xx)
   xx - frequencies of the periodogram   
   atoms - starting atoms/best atoms
   nothere - at which indices is no atom allowed
   mass - minimum mass at an index
   mind - minimum distance between knots */

   int nx=dims[0],maxd=dims[1],atd=dims[2],maxk=dims[3],atk=dims[4],nk=dims[5];
   int maxs=dims[6],ats=dims[7],ns=dims[8],er=dims[11],mind=dims[12],*fl;
   int nd=dims[5]+dims[8],add=0,i,j,*spk,cank,cans,id=dims[5];
   double **info,*score,**basis,*beta,logl,*kts,aic,aicmn,**coef,*bb;
   double ***coef2,**cumul;

/* nx - sample-size
   maxd - maximum number of dimensions
   atd - should maxd be attained (1=yes, 0=no)
   maxk - maximum number of knots
   atk - should maxk be attained (1=yes, 0=no)
   nk - number of (starting) knots
   maxs - maximum number of atoms
   ats - should maxs be attained (1=yes, 0=no)
   ns - number of (starting) atoms
   er - error criterion
   mind - minimum distance between knots
   nd - number of dimensions
   add - addition stage (0), deletion stage (1) or finished (>1)
   i,j - counter
   cank - can we add a knot?
   cans - can we add a atom?
   tsrem() - removes a knot
   tsadd() - adds a knot
   spk - present atoms
   info - hessian
   score - score vector
   basis - matrix of basis functions (smooth part) (nx x nk)
   beta - coefficients of basis functions
   logl - loglikelihood present model
   kts - present knots
   aic - present aic
   aicmn - minimum aic
   tsnew() - newton raphson
   coef - translates powerbasis and B-spline basis
   coef2 - translates powerbasis and B-spline basis
   bb - misfit - residual ratios
   cumul - cumulative sums: xx to power * (1,misfit) * (1,basis)
   tsbasis() - computes the basis */

/* allocate storage */
   score=tssdvec(NBMAX);
   beta=tssdvec(NBMAX);
   info=tssdmat(NBMAX,NBMAX);
   coef=tssdmat(NBMAX,NBMAX+4);
   coef2=tssdtri(4,NBMAX,NBMAX+1);
   basis=tssdmat(nx,NBMAX);
   fl=tssivec(2*NBMAX);
   kts=tssdvec(NBMAX);
   spk=tssivec(NBMAX);
   bb=tssdvec(nx);
   cumul=xxcumul;

/* define some cumul elements */
   for(i=0;i<11+NBMAX*4;i++){
      for(j=0;j<=nx;j++){
         cumul[i][j]=0.;
      }
   }
   for(i=nx-1;i>0;i--){ 
      cumul[7][i-1]=cumul[7][i]+1.;
      cumul[8][i-1]=cumul[8][i]+xx[i];
      cumul[9][i-1]=cumul[9][i]+xx[i]*xx[i];
      cumul[10][i-1]=cumul[10][i]+xx[i]*xx[i]*xx[i];
   }

/* initialize ad and logs and aicmn */
   for(i=0;i<NBMAX;i++){
      ad[i]=2;
      logs[i]=0.;
      theta[i]=0.;
   }
   aicmn=pow((double)10.,(double)300.);

/* copy knots into kts; initialize aicmn */
   if(nk>0)add=7;
   for(i=0;i<ns;i++)spk[i]=atoms[i];
   for(i=0;i<nk;i++)kts[i]=knots[i];

/* start the loop */
   do{

/* remove a knot if add==1 */
      if(add==1) nd=tsrem(info,&nk,beta,kts,coef,&ns,spk,nd);

/* add a knot */
      if(add==0){

         if(nk==1 && id==0 && kts[0]<0.000000001){
/* replace the knot at 0 */
            nd=tsadd(basis,info,nd,nx,kts,mind,xx,bb,cumul,spk,&nk,&ns,2,0,
               nothere);
         }
         cank=1;
         cans=1;
         if(nk==0) cans=0;
         if(nk==maxk && atk==1) cank=0;
         if(ns==maxs && ats==1) cans=0;
         if(nk==maxd-maxs && ats==1) cans=0;
         if(ns==maxd-maxk && atk==1) cank=0;
         if(nothere[0]>=10) cans=0;
/* what can we add? */
         i=tsadd(basis,info,nd,nx,kts,mind,xx,bb,cumul,spk,&nk,&ns,cank,cans,
            nothere);

/* can we keep on adding?*/
         if(i> -10)nd=i;
         else add=1;
      }

/* this is the first time, and we specified knots */
      if(add==7)add=0;

/* compute the continuous basis functions */
      tsbasis(basis,kts,nx,nk,coef,xx,coef2,fl);

/* fit using newton raphson */
      do{
         logl=tsnew(data,nx,beta,&er,score,info,nk,zz,bb,coef2,xx,kts,cumul,
            nd,spk,basis,ns,nothere,mass,fl);
         if(silent==1)(void)Rprintf("==> %.2f (%d)\n",logl,nd);
         if(er<0){
            ns--;
            nd--;
         }
      }while(er<0);

/* if er==0 there was no error */
      if(er==0){

/* store the loglikelihood, if this is the best model of this dimension */
         if(ad[nd-1]==2||logl>logs[nd-1]){
            ad[nd-1]=add;
            if(nd==1)ad[nd-1]=0;
            logs[nd-1]=logl;
         }

/* compute aic */
         aic=nd*alpha-2*logl;

/* is this the best model up to now: store beta as theta, knots, aic and dims */
         if(aic<aicmn){
            aicmn=aic;
            for(i=0;i<nk;i++)knots[i]=kts[i];
            for(i=0;i<ns;i++)atoms[i]=spk[i];
            for(i=0;i<nk+4;i++){
               theta[i]=0.;
               for(j=0;j<nk;j++) theta[i]+=beta[j]*coef[j][i];
            }
            for(i=nk;i<nd;i++)theta[i+4]=beta[i];
            dims[5]=nk;
            dims[8]=ns;
         }

/* if we are deleting, and we won't improve, call it quits */
         if(add==1 && (aic-aicmn>=alpha*(nd-1) || nd==1))add=3;

/* if we are adding, and we have the maximum number of knots, start deleting */
         if(add==0){
            if(nd==maxd)add=1;
            if(nk==maxk && ns==maxs && atk==1 && ats==1)add=1;
         }

/* if we are adding, and we have make no improvement, start deleting */
         if(add==0 && (atd==0 || atk*ats==0)){
            for(i=2;i<nd-2;i++){
               if(logs[nd-1]-logs[i-1]<((nd-i)/2.-0.5) && ad[i-1]!=2){
                  add=1;
                  dims[1]=nd;
               }
            }
         }
      }

/* if *er is not 0, we quit */
      else {
         add=3;
         dims[11]=er;
      }
   }while(add==0||add==1);

   return;
}

/******************************************************************************/
/* adds a basis function */
/******************************************************************************/
static  int tsadd(basis,info,nd,n,kts,mind,xx,bb,cumul,spk,nk,ns,cank,cans,nothere)
double **basis,**info,*kts,*xx,*bb,**cumul;
int n,*nk,*ns,mind,nd,*spk,cank,cans,*nothere;
{
   int i,j,k,bestloc= -1,*ika,bestlod= -1,ix;
   double bestraoc= -1.,nowrao,**mm,r1,*aa,*ww,bestraod= -1.;

/* i,j - counter
   r1 - utility
   bestloc - index of best knot
   bestlod - index of best atom
   ika - index of present knots
   bb - present fit
   bestraoc - rao corresponding to bestloc
   bestraod - rao corresponding to bestlod
   nowrao - rao underconsideration
   tsraoc() - computes rao for a continuous knot
   tsraod() - computes rao for a atom
   aa,ww - initialized for use in tsraoc() */

/* if this is the first knot, add it at 0 */
   if((*nk)==0){
      (*nk)++;
      kts[0]=0.;
      return nd+1;
   }

/* allocate storage */
   ika=uuika;
   mm=uumm;
   aa=uuaa;
   ww=uuww;

/* the even a-elements (later coefficients of the extra basis function) */
   aa[2]=3*kts[0];
   aa[4]=aa[2]*kts[0];
   aa[6]=aa[4]*kts[0]/3.;

/* invert the hessian */
   for(i=0;i<nd;i++) for(j=0;j<nd;j++) mm[i][j]=info[i][j];
   tsluinverse(mm,nd);

/* compute some cumul elements */
   if(cank>0){
      for(i=n-1;i>0;i--){
         for(j=0;j<(*nk);j++){
            k=11+4*j;
            r1=basis[i][j]*bb[i];
            cumul[k][i-1]=cumul[k][i]+r1;
            r1=r1*xx[i];
            cumul[k+1][i-1]=cumul[k+1][i]+r1;
            r1=r1*xx[i];
            cumul[k+2][i-1]=cumul[k+2][i]+r1;
            r1=r1*xx[i];
            cumul[k+3][i-1]=cumul[k+3][i]+r1;
         }
      }

/* find the indices of the present knots */
      j=0;
      for(i=0;(i<(n-1)&&j<(*nk));i++){
         if(kts[j]<(xx[i]+xx[i+1])/2.){
            ika[j]=i;
            j++;
         }
      }
      if(j<(*nk))ika[j]=n-1;
      ika[*nk]=n+mind;
      if(ika[0]==1)ika[0]=2;

/* some stuff useful for tsraoc() */
      ww[NBMAX]=cumul[10][ika[0]]-cumul[3][ika[0]];
      ww[NBMAX]+=aa[2]*(cumul[2][ika[0]]-cumul[9][ika[0]]);
      ww[NBMAX]-=aa[4]*(cumul[1][ika[0]]-cumul[8][ika[0]]);
      ww[NBMAX]+=aa[6]*(cumul[0][ika[0]]-cumul[7][ika[0]]);
      for(j=0;j<(*nk);j++){
         i=13+j*4;
         ww[j]=cumul[1+i][ika[0]]-aa[2]*cumul[i][ika[0]]
            +aa[4]*cumul[i-1][ika[0]]-aa[6]*cumul[i-2][ika[0]];
      }

/* search in between pairs of knots */
      bestraoc= -1.;
      j=0;
      i=1;
      if(ika[0]<=mind){
         i=ika[0]+mind+1;
         if(ika[0]==0)i++;
         j=1;
      }
      ix=i;
      for(i=ix;i<n;i++){
         if(i<ika[j]-mind){  
            nowrao=
               tsraoc(i,mm,n,(*nk),xx,kts,cumul,ika[0],ww,aa,bb,basis,nd,spk);
            if(nowrao>bestraoc){
               bestraoc=nowrao;
               bestloc=i;
            }
         }
         else{
            i=ika[j]+mind;
            j++;
         }
      }
   }
 
/* discrete search */ 
   if(cans>0){
      j=0;
      spk[(*ns)]=n+1;
      for(i=1;i<n;i++){
         if(i==spk[j]){
            j++;
        /*  i++;  remove this line if atoms are next to each other */
         }
         else {
            if(bb[i]>1. && nothere[i]==0){
               nowrao=tsraod(i,mm,(*nk),nd,xx,bb,basis,n);
               if(nowrao>bestraod){
                  bestraod=nowrao;
                  bestlod=i;
               }
            }
         }
      }
   }
/* is there anything? */  
   if(bestloc<0 && bestlod<0)return -100;

/* record the knot, sort the knots */
   if(bestraoc>bestraod){
      if(bestloc<n)kts[(*nk)]=xx[bestloc];
      else kts[(*nk)]=PIL;
     if(silent==1)(void)Rprintf("add knot at %.3f (%.3f)  ",kts[(*nk)],bestraoc);
      for(i=(*nk);i>0;i--)if(kts[i]<kts[i-1]){
         bestraoc=kts[i-1];
         kts[i-1]=kts[i];
         kts[i]=bestraoc;
      }
      (*nk)++;
      if(cank==2){
         kts[0]=kts[1];
         (*nk)=1;
         nd--;
      }
   }
   else{
      spk[(*ns)]=bestlod;
     if(silent==1)(void)Rprintf("add atom at %.3d (%.3f)  ",spk[(*ns)],bestraod);
      for(i=(*ns);i>0;i--)if(spk[i]<spk[i-1]){
         bestlod=spk[i-1];
         spk[i-1]=spk[i];
         spk[i]=bestlod;
      }
      (*ns)++;
   }

   return nd+1;
}

/******************************************************************************/
/* computes the rao statistic atom */
/******************************************************************************/
static double tsraod(loc,mm,nk,nd,xx,bb,basis,n)
double *xx,**mm,**basis,*bb;
int nd,nk,loc,n;
{
   int i,j;
   double *vv,bbn,r1;
   vv=uuvv1;
/* the completely discrete elements */
   bbn=1.-bb[loc];
   vv[nd]= -bb[loc];
   for(i=0;i<nk;i++) vv[i]= -bb[loc]*basis[loc][i];
   for(i=nk;i<nd;i++) vv[i]= 0.;

/* if there is a atom at pi */
   if(loc==n-1 && xx[n-1]>=3.1415926){
      bbn=0.5*bbn;
      vv[nd]=0.5*vv[nd];
      vv[0]=0.5*vv[0];
      if(nk>1){
         vv[1]=0.5*vv[1];
         if(nk>3){
            vv[3]=0.5*vv[3];
         }
      }
   }

/* compute */
   r1=0.;
   for(i=0;i<nd;i++) for(j=0;j<i;j++) r1-= vv[i]*vv[j]*mm[i][j];
   r1=r1*2.;
   for(i=0;i<nd;i++) r1-= vv[i]*vv[i]*mm[i][i];
   r1+=vv[nd];
   bbn= -bbn*bbn/r1;

   return bbn;
}

/******************************************************************************/
/* computes the rao statistic continuous knot */
/******************************************************************************/
static double tsraoc(loc,mm,n,nk,xx,kts,cumul,ika,ww,aa,bb,basis,nd,spk)
int n,nd,loc,ika,nk,*spk;
double *xx,*kts,**mm,**cumul,*aa,*ww,*bb,**basis;
{
   int i,j;
   double *vv,bbn,xlc,r1,b1,b2,aax[8],bsx,bbpi,b1pi,b3pi,x,y;

/* i,j - counter
   vv - extra column hessian, work copy
   bbn - extra element score
   xlc - new knot 
   b1,b2,aax - parameters 
   r1 - utility */

/* allocate storage, initialize */
   bbpi=bb[n-1];
   b1pi=basis[n-1][1];
   b3pi=basis[n-1][3];
   vv=uuvv2;
   for(i=0;i<nk;i++) vv[i]=ww[i];
   vv[nd]=0.;


/* location */
   if(loc!=n) xlc=xx[loc];
   else xlc=PIL;
   if(loc==n)loc=n-1;

/* parameters */
   aa[1]= (3.*((PIL-kts[0])*(PIL-kts[0])-(PIL-xlc)*(PIL-xlc)))/(2*PIL);
   aa[3]=3*xlc; 
   aa[5]=aa[3]*xlc; 
   aa[7]=aa[5]*xlc/3.;

/* extra element score */
   bbn=ww[NBMAX]+aa[1]*(cumul[2][0]-cumul[9][0])+cumul[3][loc]-cumul[10][loc]
       -aa[3]*(cumul[2][loc]-cumul[9][loc])+aa[5]*(cumul[1][loc]-cumul[8][loc])
       -aa[7]*(cumul[0][loc]-cumul[7][loc]);

/* extra column hessian */
   for(j=0;j<nk;j++){
      i=13+j*4;
      vv[j]+= -aa[1]*cumul[i][0]-cumul[1+i][loc]+aa[3]*cumul[i][loc]
             -aa[5]*cumul[i-1][loc]+aa[7]*cumul[i-2][loc];
   }
   if(xlc>kts[0]){
      vv[nd]=(cumul[4][ika]-cumul[4][0])*aa[1]*aa[1];

/* sequence important */
      aax[3]=aa[1]+aa[2]-aa[3]; 
      aax[5]=aa[5]-aa[4]; 
      aax[7]=aa[6]-aa[7]; 
      aax[2]=aa[1]+aa[2]; 
      for(i=0;i<7;i++){
         if(i>3){
            if(i==6){
               b1=1.;
            }
            else{
               if(i==5){
                  b1= -2.*aax[2];
               }
               else{
                  b1=2.*aa[4]+aax[2]*aax[2]; 
                  b2=aax[3]*aax[3];
               }
            }
         }
         else{
            if(i>1){
               if(i==3){
                  b1= -2.*aa[6]-2.*aax[2]*aa[4]; 
                  b2=2.*aax[3]*aax[5];
               }
               else{
                  b1=aa[4]*aa[4]+2.*aa[6]*aax[2]; 
                  b2=aax[5]*aax[5]+2.*aax[7]*aax[3];
               }
            }
            else{
               if(i==1){
                  b1= -2.*aa[4]*aa[6]; 
                  b2=2.*aax[5]*aax[7];
               }
               else{
                  b1=aa[6]*aa[6]; 
                  b2=aax[7]*aax[7];
               }
            }
         }
         vv[nd]+=b1*(cumul[i][loc]-cumul[i][ika]);
         if(i<5)vv[nd]-=b2*cumul[i][loc];
      } 
   }
   else{
      vv[nd]=(cumul[4][loc]-cumul[4][0])*aa[1]*aa[1];

/* sequence important */
      aax[2]=aa[1]+aa[2]-aa[3]; 
      aax[4]=aa[5]-aa[4]; 
      aax[6]=aa[6]-aa[7]; 
      aax[3]=aa[1]-aa[3]; 
      for(i=0;i<7;i++){
         if(i>3){
            if(i==6){
               b1=1.;
            }
            else{
               if(i==5){
                  b1=2.*aax[3];
               }
               else{
                  b1=2.*aa[5]+aax[3]*aax[3]; 
                  b2=aax[2]*aax[2];
               }
            }
         }
         else{
            if(i>1){
               if(i==3){
                  b1= -2.*aa[7]+2.*aax[3]*aa[5]; 
                  b2=2.*aax[2]*aax[4];
               }
               else{
                  b1=aa[5]*aa[5]-2.*aa[7]*aax[3]; 
                  b2=aax[4]*aax[4]+2.*aax[6]*aax[2];
               }
            }
            else{
               if(i==1){
                  b1= -2.*aa[5]*aa[7]; 
                  b2=2.*aax[4]*aax[6];
               }
               else{
                  b1=aa[7]*aa[7]; 
                  b2=aax[6]*aax[6];
               }
            }
         }
         vv[nd]+=b1*(cumul[i][ika]-cumul[i][loc]);
         if(i<5)vv[nd]-=b2*cumul[i][ika];
      }
   }

/* subtract half times the last element */
   if(xx[n-1]>3.1415926){
      bsx=PIL-kts[0];
      bsx=(aa[1]*PIL*PIL+(PIL-xlc)*(PIL-xlc)*(PIL-xlc)-bsx*bsx*bsx)/2.;
      bbn+=bsx*(1.-bbpi);
      bbpi=bsx*bbpi;
      vv[nd]+=bsx*bbpi*2.;
      vv[0]+=bbpi;
      if(nk>1){
         vv[1]+=bbpi*b1pi;
         if(nk>3)vv[3]+=bbpi*b3pi;
      }
   }

   for(i=nk;i<nd;i++){
      x=xx[spk[i-nk]];
      y=aa[1]*x*x;
      if(x>xlc)y+=(x-xlc)*(x-xlc)*(x-xlc);
      if(x>kts[0])y-=(x-kts[0])*(x-kts[0])*(x-kts[0]);
      vv[i]= -bb[spk[i-nk]]*y;
      if(x>3.1415926)vv[i]=0.5*vv[i];
   }
   if(xlc>xx[n-1])loc=n;
   

/* compute*/
   r1=0.;
   for(i=0;i<nd;i++) for(j=0;j<i;j++) r1-= vv[i]*vv[j]*mm[i][j];
   r1=r1*2.;
   for(i=0;i<nd;i++) r1-= vv[i]*vv[i]*mm[i][i];
   r1+=vv[nd];
   bbn= -bbn*bbn/r1;

   return bbn;
}

/******************************************************************************/
/* removes a knot */
/******************************************************************************/
static int tsrem(info,nk,beta,kts,coef,ns,spk,nd)
double **coef,*beta,**info,*kts;
int *nk,*ns,*spk,nd;
{
   int irmax=4,i,j,k;
   double se,phi,ratmax=0;

/* irmax - index of function to be removed
   i,j,k - counters
   se - standard error of phi
   phi - one power basis coefficient
   ratmax - maximum ratio |phi|/se */

/* the last one, who cares? */
   if(nd==2){
      (*nk)=1;
      (*ns)=0;
      return 1;
   }

/* invert */
   tsluinverse(info,nd);

/* for each knot */
   if((*nk)>1){
      for(i=0;i<(*nk);i++){
         phi = 0.;
         se = 0.;

/* compute phi and se^2 */
         for(j=0;j<(*nk);j++){
            phi+=beta[j]*coef[j][i+4];
            for(k=0;k<(*nk);k++) se-=coef[j][i+4]*coef[k][i+4]*info[j][k];
         }

/* what is the ratio? */
         phi=fabs(phi);
         if(se>0) se=sqrt(se);
         else se=0.;
         if(se>phi*ratmax){
            ratmax=se/phi;
            irmax=i;
         }
      }
   }

/* for each atom */
   for(i=(*nk);i<nd;i++){
      phi=fabs(beta[i]);
      se= -info[i][i];
      if(se>0) se=sqrt(se);
      else se=0.;
      if(se>phi*ratmax){
         ratmax=se/phi;
         irmax=i;
      }
   }

/* remove the loser  - if it is a knot */
   if(irmax<(*nk)){
      if(silent==1)(void)Rprintf("del knot at %.3f (",kts[irmax]);
      for(j=irmax;j<((*nk)-1);j++)kts[j]=kts[j+1];
      (*nk)-=1;
   }

/* remove the loser  - if it is a atom */
   else{
      irmax-=(*nk);
      if(silent==1)(void)Rprintf("del atom at %.3d (",spk[irmax]);
      for(j=irmax;j<((*ns)-1);j++)spk[j]=spk[j+1];
      (*ns)-=1;
   }
   if(silent==1)Rprintf("%.3f)  ",1./ratmax);
   return nd-1;
}

/******************************************************************************/
/* does the newton raphson iterations */
/******************************************************************************/
static double tsnew(data,n,beta,er,score,info,nk,zz,ff,coef2,xx,kk,cumul,
		           nd,spk,basis,ns,nothere,mass,fl)
double *data,*beta,*score,**info,*zz,*ff,***coef2,*xx,*kk,**cumul,**basis,mass;
int n,*er,nk,ns,nd,*spk,*nothere,*fl;
{
   int i,j,k,i1,jx;
   double logold,lognew,*betan,r,zerror;
   double uu[7],vv[7],xz;
   int k1,k2,k3,l,l2;

/* i,j,k - counter
   tslogall() computes score, hessian and loglikelihood
   logold - old loglikelihood
   lognew - new loglikelihood
   betan - new beta
   r - utility
   zerror - convergence criterion */

/* allocate, initialize */
   betan=uubetan;
   *er=0;
   for(i=0;i<nd;i++){
      beta[i]=0.;
      for(j=0;j<=i;j++)info[i][j]=0.;
   }

/* starting values by least squares */
   k2=0;
   k1=0;
   k3=0;
   for(j=0;j<4;j++)vv[j]=0;
   if(fabs(xx[0]-kk[0])<TINY)k1=1;
/* continuous part */
   for(k=1;k<n;k++){
      vv[0]+=zz[k];
      vv[1]+=zz[k]*xx[k];
      vv[2]+=zz[k]*xx[k]*xx[k];
      vv[3]+=zz[k]*xx[k]*xx[k]*xx[k];
      if(k==(n-1)||(xx[k+1]>=kk[k1] && k1<nk))k3=k;
      if(k3>0){
         for(i=0;i<nk;i++)for(j=0;j<4;j++)beta[i]+=coef2[j][i][k1]*vv[j];
         tsintsum(uu,k2,k3,xx[1]);
         for(i=0;i<nk;i++)for(j=0;j<=i;j++)for(l=0;l<4;l++)for(l2=0;l2<4;l2++)
            info[i][j]+=coef2[l][i][k1]*coef2[l2][j][k1]*uu[l+l2];
         k1++;
         for(j=0;j<4;j++)vv[j]=0;
         k3=0;
         k2=k;
      }
   }
/* discrete part */
   for(i=nk;i<nd;i++){
      beta[i]=zz[spk[i-nk]];
      info[i][i]=1.;
      for(j=0;j<nk;j++){
         info[i][j]=basis[spk[i-nk]][j];
         info[j][i]=basis[spk[i-nk]][j];
      }
   }
   for(i=0;i<nd;i++)for(j=i+1;j<nd;j++)info[i][j]=info[j][i];
   tslusolve(info,nd,beta);

/* the fun starts */
   for(i=0;i<100;i++){

/* compute it all */
      logold=tslogall(ff,beta,score,info,n,nk,data,1,xx,cumul,basis,nd,spk,ns,
              fl);

/* solve the system */
      tslusolve(info,nd,score);

/* what is the new beta? */
      for(j=0;j<nd;j++)betan[j]=beta[j]-score[j];

/* what is its loglikelihood, what is the zerror */
      lognew=tslogall(ff,betan,score,info,n,nk,data,0,xx,cumul,basis,nd,spk,ns,
              fl);
      zerror=lognew-logold;

/* do we need to step-half? */
      if(lognew<logold-1.0E-7){
         zerror=1.;
         r=1.;
         do{

/* step halving problems */
            if(r<0.000001){
               (*er)=1;
               return logold;
            }

/* step halving */
            r=r/2.;
            for(j=0;j<nd;j++)betan[j]=beta[j]-score[j]*r;
            lognew=tslogall(ff,betan,score,info,n,nk,data,0,
               xx,cumul,basis,nd,spk,ns,fl);
         }while(lognew<logold-1.0E-7);
      }

/* record the new beta */
      for(j=0;j<nd;j++)beta[j]=betan[j];

/* did we converge? */
      if(zerror<0.000001){

/* compute the good stuff */
         lognew=tslogall(ff,beta,score,info,n,nk,data,2,
            xx,cumul,basis,nd,spk,ns,fl);
         for(j=nk;j<nd;j++){
            xz=0.;
            for(i1=0;i1<nk;i1++)xz+=beta[i1]*basis[spk[j-nk]][i1];
            beta[j]=exp(xz+beta[j])-exp(xz);
         }
/* too small atoms */
         for(j=nk;j<nd;j++){
            if(beta[j]<mass){
               nothere[spk[j-nk]]=1;
               nothere[0]++;
               jx=j;
               for(j=jx;j<nd-1;j++)spk[j-nk]=spk[j+1-nk];
               (*er)= -1;
               return 0.;
            }
         }
         return lognew;
      }
   }

/* no convergence */
   (*er)=2;
   return logold;
}

/******************************************************************************/
/* compute loglikelihood, hessian and score */
/******************************************************************************/
static double tslogall(bb,beta,score,info,n,nk,data,what,xx,cumul,basis,nd,spk,ns,fl)
double *beta,*score,**info,*data,*bb,*xx,**cumul,**basis;
int n,nk,what,nd,ns,*spk,*fl;

/* what: 0 loglikelihood only, 1 also score and hessian */
{
   int i,j,k,uu,kx;
   double logl,r1,b1,b3;

/* i,j,k - counter
   logl - loglikelihood
   bb - fitted value */

/* initializations */
   b1=basis[n-1][1];
   b3=basis[n-1][3];
   logl=0.;
   if(what>0){
      for(i=0;i<nd;i++){
         score[i]=0.;
         for(j=0;j<nd;j++)info[i][j]=0.;
      }
   }

/* get the fitted values continuous part */
   for(k=1;k<n;k++){
      bb[k]=0;
      for(i=0;i<nk;i++)bb[k]+=beta[i]*basis[k][i];
   }

/* get the fitted values discrete part */
   for(k=nk;k<nd;k++) bb[spk[k-nk]]+=beta[k];

/* get the log-likelihood */
   for(k=1;k<n;k++){
      logl-=bb[k];
      bb[k]=data[k]*exp(-bb[k]);
      logl-=bb[k];
   }

/* correct if the last observation is at PIL */
   if(xx[n-1]>=3.1415926) logl+=0.5*(bb[n-1]-log(bb[n-1]/data[n-1]));

/* update score and hessian - the continuous elements */
   if(what==2){
      for(k=n-1;k>0;k--){
         r1=bb[k];
         cumul[0][k-1]=cumul[0][k]+bb[k];
         for(j=1;j<7;j++){
            r1=r1*xx[k];
            cumul[j][k-1]=cumul[j][k]+r1;
         }
      }
   }
   if(what>0){
      for(i=0;i<nk;i++){
         for(k=fl[2*i];k<fl[2*i+1];k++)score[i]-=basis[k][i]*(1.-bb[k]);
         for(j=0;j<=i;j++){
            k=fl[2*i];
            if(fl[2*j]>k)k=fl[2*j];
            uu=fl[2*i+1];
            if(fl[2*j+1]>uu)uu=fl[2*j+1];
            kx=k;
            for(k=kx;k<uu;k++)info[i][j]-=basis[k][i]*bb[k]*basis[k][j];
         }
      }

/* correct if the last observation is at PIL - only basis 0, 1 and 3 */
      if(xx[n-1]>=3.1415926){
         bb[n-1]=0.5*bb[n-1];
         score[0]+=(0.5-bb[n-1]);
         info[0][0]+=bb[n-1];
         if(nk>1){
            score[1]+=b1*(0.5-bb[n-1]);
            info[1][0]+=b1*bb[n-1];
            info[1][1]+=b1*b1*bb[n-1];
            if(nk>3){
               score[3]+=b3*(0.5-bb[n-1]);
               info[3][0]+=b3*bb[n-1];
               info[3][1]+=b3*b1*bb[n-1];
               info[3][3]+=b3*b3*bb[n-1];
            }
         }
         bb[n-1]=2.*bb[n-1];
      }
      for(i=0;i<nk;i++)for(j=i+1;j<nk;j++)info[i][j]=info[j][i];

/* the completely discrete elements */
      for(i=nk;i<nd;i++){
         score[i]=1.-bb[spk[i-nk]];
         info[i][i]= -bb[spk[i-nk]];
      }

/* the mixed info elements */
      for(i=nk;i<nd;i++){
         for(j=0;j<nk;j++){
            info[i][j]= -bb[spk[i-nk]]*basis[spk[i-nk]][j];
            info[j][i]= info[i][j];
         }
      }

/* if there is a atom at pi */
      if(ns>0){
         if(spk[ns-1]==n-1 && xx[n-1]>=3.1415926){
            score[nd-1]=0.5*score[nd-1];
            info[nd-1][nd-1]=0.5*info[nd-1][nd-1];
            info[0][nd-1]=0.5*info[0][nd-1];
            info[nd-1][0]=info[0][nd-1];
            if(nk>1){
               info[1][nd-1]=0.5*info[1][nd-1];
               info[nd-1][1]=info[1][nd-1];
               if(nk>3){
                  info[3][nd-1]=0.5*info[3][nd-1];
                  info[nd-1][3]=info[3][nd-1];
               }
            }
         }
      }
   }

   return logl;
}

/******************************************************************************/
/* compute the basis functions */
/******************************************************************************/
static void tsbasis(basis,kk,n,nb,coef,xx,coef2,fl)
double **basis,*kk,**coef,*xx,***coef2;
int n,nb,*fl;
{
   int i,j,k;
   double r1;

/* i,j,k - counters
   tsb1() -get coef for basis function 1
   tsb2() -get coef for basis function 2
   tsb3() -get coef for basis function 3
   tsb4() -get coef for basis function 4
   tsb5() -get coef for basis function 5 and later */

/* initialize */
   for(i=0;i<NBMAX;i++) for(j=0;j<NBMAX+4;j++) coef[i][j]=0.;

/* get the coefficients */
   tsb1(coef,0,0);
   if(nb>1)tsb2(coef,kk,nb+2,1,nb-2);
   if(nb>2)tsb3(coef,kk,4,2,0);
   if(nb>3)tsb4(coef,kk,nb,3,nb-4);
   for(i=5;i<=nb;i++)tsb5(coef,kk,i-1,i-1,i-5);

/* compute the values */
   for(i=0;i<n;i++){
      basis[i][0]=1.;
      for(j=1;j<nb;j++) basis[i][j]=coef[j][0]+xx[i]*xx[i]*coef[j][2];
   }
   for(k=4;k<nb+4;k++) {
      for(i=n-1;i>=0;i--){
         if(xx[i]>kk[k-4]){
            r1=(xx[i]-kk[k-4])*(xx[i]-kk[k-4])*(xx[i]-kk[k-4]);
            for(j=1;j<nb;j++) basis[i][j]+= coef[j][k]*r1;
         }
         else i= -1;
      }
   } 
   for(i=0;i<nb;i++){
      fl[2*i]=1;
      fl[2*i+1]=n;
      if(i>=4) {
         fl[2*i]=(int)(kk[i-4]*(double)n/PIL)-1;
         if(fl[2*i]<1) fl[2*i]=1;
         fl[2*i+1]=(int)(kk[i]*(double)n/PIL)+1;
         fl[2*i+1]=n;
      }
   }

   for(i=0;i<nb;i++)for(j=0;j<nb+1;j++)for(k=0;k<4;k++)coef2[k][i][j]=0;

/* basisfunction */
   for(i=0;i<nb;i++){

/* interval */
      for(j=0;j<=nb;j++){
         coef2[0][i][j]+=coef[i][0];
         coef2[2][i][j]+=coef[i][2];

/* knot */
         for(k=0;k<j;k++){
            coef2[0][i][j]+= -coef[i][k+4]*kk[k]*kk[k]*kk[k];
            coef2[1][i][j]+= 3.*coef[i][k+4]*kk[k]*kk[k];
            coef2[2][i][j]+= -3.*coef[i][k+4]*kk[k];
            coef2[3][i][j]+= coef[i][k+4];
         }
      }
   }
}

/******************************************************************************/
/* basis functions 5 and beyond - B-splines */
/******************************************************************************/
static void tsb5(coef,kk,col,row,i)
int col,row,i;
double **coef,*kk;
{
   double uu[5],vv[5],r0,r1;
   int j;
   uu[0]=1.;
   uu[4]=0.;
   uu[2]= (kk[i+1]-kk[i+0])*(kk[i+3]-kk[i+0])/
                           ((kk[i+2]-kk[i+1])*(kk[i+3]-kk[i+2]));
   uu[3]= (kk[i+1]-kk[i+0])*(kk[i+0]-kk[i+2])/
                           ((kk[i+3]-kk[i+1])*(kk[i+3]-kk[i+2]));
   uu[1]= -1.-uu[3]-uu[2];
   r0=0.;
   for(j=0;j<4;j++)
      r0+=uu[j]*(kk[i+4]-kk[i+j])*(kk[i+4]-kk[i+j])*(kk[i+4]-kk[i+j]);
   vv[1]=1.;
   vv[0]=0.;
   vv[3]= (kk[i+2]-kk[i+1])*(kk[i+4]-kk[i+1])/
                             ((kk[i+3]-kk[i+2])*(kk[i+4]-kk[i+3]));
   vv[4]= (kk[i+2]-kk[i+1])*(kk[i+1]-kk[i+3])/
                             ((kk[i+4]-kk[i+2])*(kk[i+4]-kk[i+3]));
   vv[2]= -1.-vv[3]-vv[4];
   r1=0.;
   for(j=1;j<5;j++)
      r1+=vv[j]*(kk[i+4]-kk[i+j])*(kk[i+4]-kk[i+j])*(kk[i+4]-kk[i+j]);
   for(j=0;j<5;j++)coef[row][col+j]=uu[j]-r0/r1*vv[j];
}

/******************************************************************************/
/* basis function 1 - constant */
/******************************************************************************/
static void tsb1(coef,col,row)
int col,row;
double **coef;
{
   coef[row][col]=1.;
}

/******************************************************************************/
/* basis function  2 squarish */
/******************************************************************************/
static void tsb2(coef,kk,col,row,i)
int col,row,i;
double **coef,*kk;
{
   coef[row][2]=1.;
   coef[row][col]=2*PIL/(3.*((PIL-kk[i+1])*(PIL-kk[i+1])-(PIL-kk[i])*(PIL-kk[i])));
   coef[row][col+1]= -coef[row][col];
}

/******************************************************************************/
/* basis function 3 squarish */
/******************************************************************************/
static void tsb3(coef,kk,col,row,i)
int col,row,i;
double **coef,*kk;
{
   double r;
   coef[row][0]=1.;
   coef[row][col+1]=PIL/(3.*(kk[i+1]-kk[i]));
   coef[row][col]= -coef[row][col+1];
   coef[row][2]=coef[row][col+1]*(PIL-kk[i+1])*(PIL-kk[i+1]);
   coef[row][2]+=coef[row][col]*(PIL-kk[i])*(PIL-kk[i]);
   coef[row][2]= -(coef[row][2]*3.)/(2.*PIL);
   r=1+PIL*PIL*coef[row][2]+coef[row][col]*(PIL-kk[i])*(PIL-kk[i])*(PIL-kk[i])+
     coef[row][col+1]*(PIL-kk[i+1])*(PIL-kk[i+1])*(PIL-kk[i+1]);
   coef[row][col+1]=coef[row][col+1]/(1.-r);
   coef[row][col]=coef[row][col]/(1.-r);
   coef[row][2]=coef[row][2]/(1.-r);
}

/******************************************************************************/
/* basis function 4 S-ish */
/******************************************************************************/
static void tsb4(coef,kk,col,row,i)
int col,row,i;
double **coef,*kk;
{
   coef[row][col]=1.;
   coef[row][col+2]= (kk[i+1]-kk[i])*(kk[i+3]-kk[i])/
                  ((kk[i+2]-kk[i+1])*(kk[i+3]-kk[i+2]));
   coef[row][col+3]= (kk[i+1]-kk[i])*(kk[i]-kk[i+2])/
                  ((kk[i+3]-kk[i+1])*(kk[i+3]-kk[i+2]));
   coef[row][col+1]= -1.-coef[row][col+3]-coef[row][col+2];
}

/******************************************************************************/
/* rest is adapted numerical recipes and cmlib stuff */
/******************************************************************************/
/* static void tsrerror(error_text)                                           */
/* char error_text[];                                                         */
/* {                                                                          */
/*    void exit();                                                            */
/*    error("%s\n this is serious.....\n",error_text);                        */
/*    exit(1);                                                                */
/* }                                                                          */
/*                                                                            */
/******************************************************************************/
/* i-vector */
/******************************************************************************/
static int *tssivec(nh)
int nh;
{
   int *v,i;
   v=(int *)Salloc(nh+1,int);
   for(i=0;i<=nh;i++)v[i]=0;
   return v;
}


/******************************************************************************/
/* d-vector */
/******************************************************************************/
static double *tssdvec(nh)
int nh;
{
   double *v;
   int i;
   v=(double *)Salloc(nh+1,double);
   for(i=0;i<=nh;i++)v[i]=0.;
   return v;
}

/******************************************************************************/
/* d-3d array */
/******************************************************************************/
static double ***tssdtri(r,c,s)
int r,c,s;
{
   int i,j,k;
   double ***m;
   m=(double ***) Salloc(r+1,double**); 
   for(i=0;i<=r;i++) m[i]=(double **) Salloc(c+1,double*);
   for(i=0;i<=r;i++) for(j=0;j<=c;j++){
       m[i][j]=(double *) Salloc(s+1,double); 
       for(k=0;k<=s;k++)m[i][j][k]=0.;
   }
   return m;
}

/******************************************************************************/
/* d-matrix */
/******************************************************************************/

static double **tssdmat(nrh,nch)
int nrh,nch;
{
   int i,j;
   double **m;
   m=(double **) Salloc(nrh+1,double*); 
   for(i=0;i<=nrh;i++){
      m[i]=(double *) Salloc(nch+1,double); 
      for(j=0;j<=nch;j++)m[i][j]=0.;
   }
   return m;
}
/******************************************************************************/
/* cumulative sums */
/******************************************************************************/
static void tsintsum(rr,k0,k1,f)
int k0,k1;
double f,rr[];
{
   double ff,l0,l1,m0,m1;
   l0=(double)k0*(k0+1.);
   l1=(double)k1*(k1+1.);
   m0=(double)2.*k0+1.;
   m1=(double)2.*k1+1.;
   rr[0]=(double)(k1-k0);
   rr[1]=f*(l1-l0)/2.;
   ff=f*f;
   rr[2]=ff*(l1*m1-l0*m0)/6.;
   ff=ff*f;
   rr[3]=ff*(l1*l1-l0*l0)/4.;
   ff=ff*f;
   rr[4]=ff*(l1*m1*(3*l1-1)-l0*m0*(3*l0-1))/30.;
   ff=ff*f;
   rr[5]=ff*(l1*l1*(2*l1-1)-l0*l0*(2*l0-1))/12.;
   rr[6]=ff*f*(l1*m1*(1.-3.*k1*(1.-(double)k1*k1*(double)(2.+k1)))
                     -l0*m0*(1.-3.*k0*(1.-(double)k0*k0*(double)(2.+k0))))/42;
}
/******************************************************************************/
static void tslusolve(a,n,b)
int n;
double **a,*b;
{
   static double aa[DIM5][DIM5],bb[DIM5];
   static int kpvt[DIM5];
   int i,j,info;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
      bb[i]=b[i];
   }
   i=DIM5;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
/* if(info|=0)tsrerror("oops in lusolve"); */
   F77_CALL(xdsisl)(aa,&i,&n,kpvt,bb);
   for(i=0;i<n;i++)b[i]=bb[i];
}
/******************************************************************************/
static void tsluinverse(a,n)
int n;
double **a;
{
   static double aa[DIM5][DIM5],bb[DIM5];
   double det[2];
   static int kpvt[DIM5];
   int i,j,info,inert[3];
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
   }
   i=DIM5;
   j=1;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
/* if(info|=0)tsrerror("oops in luinverse"); */
   F77_CALL(xdsidi)(aa,&i,&n,kpvt,det,inert,bb,&j);
   for(i=0;i<n;i++){
      for(j=0;j<i;j++)a[i][j]=aa[i][j];
      for(j=i;j<n;j++)a[i][j]=aa[j][i];
   }
}
/******************************************************************************/
static void tsallocer()
{
   uuika=tssivec(NBMAX+1);

   uumm=tssdmat(NBMAX,NBMAX);

   uuaa=tssdvec(10);
   uuww=tssdvec(NBMAX+1);
   uuvv1=tssdvec(NBMAX);
   uuvv2=tssdvec(NBMAX);
   uubetan=tssdvec(NBMAX+10);
}
