/*
*  Copyright (C) 1995--2002  Charles Kooperberg
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
/* this function describes the basic structures */

#include <math.h>
#include <stdio.h>
#include "R.h"
#define Salloc(n, t)  (t *)R_alloc((long)(n), (int)sizeof(t))


/* we want to be able to use those everywhere */

#define MAXSPACE 250
#define MAXKNOTS 10
#define DIM5 MAXSPACE+5

void F77_NAME(xdsifa)(double[][DIM5], int *, int *, int *, int *);
void F77_NAME(xdsisl)(double[][DIM5], int *, int *, int *, double *);
void F77_NAME(xdsidi)(double[][DIM5], int *, int *, int *, double *, int *, double *, int *);
void F77_NAME(xdgefa)(double[][DIM5], int *, int *, int *, int *);
void F77_NAME(xdgedi)(double[][DIM5], int *, int *, int *, double *, double *, int *);

/* MAXSPACE - maximum dimensionality of the model
   MAXKNOTS - maximum number of knots for one covariate */

struct datastruct {
   int ndata,ncov,*bincov,nclass,*yy,*icov;
   double **work,**work2,*wgt,wgtsum;
};

/* datastruct is a structure containing all information about the data. At any 
   time there is only one datastruct, which is typically called data.

   ndata  - number of datapoints 
   ncov   - number of covariates 
   nclass - how many classes are there
   bincov - are the covariates binary? 0=no, 1=yes binary cov should be 0-1 
   yy     - response
   cov    - covariates cov[i][j] is covariate j for observation i 
   work   - also to keep exp(theta-c)
   work2  - also to keep theta-bar(theta)
   wgt    - case weights */

struct space {
   int ndim,nbas;
   double aic,**info,*score,**infox,epsilon,logl;
   struct basisfunct *basis; 
   struct subdim **sub; 
};

/* space is the basic structure containing a model. The main ingredients are a
   (sort of double) representation of the basisfunctions: by means of basis on
   a basisfunction by basisfunction scale and by means of sub on a subdimension
   scale

   ndim    - the dimensionality of the space
   nbas    - number of basis functions ndim=nbas*(nclass-1)
   aic     - the aic value of the present model - only accurate after the model
             has been fitted.
   info    - the hessian of the model
   infox   - epsilon term for the hessian
   epsilon - epsilon for the penalty term
   score   - the score vector of the model
   basis   - the array of basisfunctions
   sub     - the matrix of subdimensions element:
             [i][data.ncov]  (0<=i<data.ncov) the subdimension belonging to a
                                              covariate by itself.
             [i][j]  (0<=i<j<data.ncov)       the subdimension belonging to the
                                              tensor-product of two covariates
             others                           not used. */

struct basisfunct {
   int b1,b2,t1,t2,*link1,*link2,ib,j1,j2;
   double *beta;
};

/* a structure describing one basisfunction (surprising isn't it?)
   b1 and b2 - indicate the subdimension to which this basisfunction belongs
   t1 and t2 - the rank of the knots in subdimensions b1 and b2, if applicable
               if these are -1, while b1 and b2 would indicate a covariate, the
               basisfunction is linear in this covariate.
   beta      - betas of this basisfunction
   ib        - number of independent basis functions
   link1     - which class referes to which basis function
               if link1[i] of basisfunct[k] is j, the i-th class belongs to
               group j for basis function k
   link2     - which class referes to which basis function
               if link2[l] of basisfunct[k] is m, the beta for group l for
               basis function k is the m-th element of the beta-vector
               thus if link2[link1[i]] of basisfunct[k] is m, the beta for
               class for basis function k is the m-th element of the beta-vector
   */

struct subdim {
   int dim1,**kts1;
   double *ktsc;
};

/* a structurte describing a subdimension - it can depend on either one
   or two covariates or one covariate and time.
   dim1   - dimensionality of component 
   kts1   - binary - if =0 indicates that basisfunction is not in the model,
            otherwise number of basisfunction.
   ktsc   - if this is a subdimension depending on 1 covariate, the knots */

/* workspaces */

static int padddim();
static double padders1();
static double padders2();
static double pearch();
static double pestbasis();
static double cripswap();
static double prao();
static int pindyl();
static int pindyr();
static int pindl();
static int pindr();
static int pindx();
static int pindm();
static int plocation();
static void paddbasis();
static void soutspace();
static double petvector2();
static void petvector();
static struct datastruct *pdefinedata();
static struct space *pdefinespace();
static struct basisfunct *pdefinebasis();
static struct subdim **pdefinedim();
static int *ispvector();
static double *dspvector();
static int **ipmatrix();
static double **dpmatrix();
static double aiccv();
static void aicb2();
static int aicbest();
static void predefinespace();
static void proj();
static double pnewton();
static double pcompall();
static int dlink();
static void getinfox();
static double pylog();
static double pcomp2();
static int prembas();
static void premdim();
static void puuu();
static void poly();
static void pconstant();
static void pswapspace();
static void computeloss();
static double getcrit();
static int plumbertester();
static void Ppsort();
static void xpsort();
static int lusolinv();

static double **w1,**w2,**w3,*v1,*v2,*v3,*v4,*v5,*v6,*v7,*v8;
static int **iw1,*iv1,*iv2;
static float *trcov,*tecov;

/* misc stuff */

static int maxdim;
/******************************************************************************/

/* this routine pearches all dimensions for something to add */

static int padddim(current,new,newt,data,mind,exclude,silent,lins)
struct space *current,*new,*newt;
struct datastruct *data;
int mind,**exclude,silent,*lins;

/* current - current space
   new     - copy of current space to play with
   newt    - best addition space up to now
   data    - data
   mind    - minimum distance (order statistics) between knots 
   exclude - if exclude = 1 this term is never included 
   silent  - should diagnostic output be printed 
   lins    - which dimensions should only be added linear */

{
   int i,j,ncov=(*data).ncov;
   double criterion;

/* i,j        - counters
   ncov       - save typing
   criterion  - criterion
   padders1   - does the work - something linear
   padders2   - does the work - something with a knot
   pswapspace - copies one space into another 
   puuu       - print some diagnostic statistics */

/* initialization */
   criterion = -pow((double)10.,(double)20.);

/* get a space to play with */
   pswapspace(newt,current,data);

/* ready for inverse: w3 is hessian inverse */
   for(i=0;i<((*current).nbas)*((*data).nclass);i++){
      for(j=0;j<((*current).nbas)*((*data).nclass);j++) w1[i][j]= -w3[i][j];
   }
   
/* pearches all the dimensions: knots*/
   for(i=0;i<ncov;i++) for(j=i+1;j<=ncov;j++){
/* only select existing dimensions */
      if(i==ncov || i!=j) if(j==ncov || exclude[i][j]==0){
         criterion=padders1(i,j,new,newt,criterion,data);
      }
   }

/* pearches all the dimensions: linears */
   for(i=0;i<ncov;i++) for(j=i+1;j<=ncov;j++){
/* only select existing dimensions */
      if(i==ncov || i!=j) if(j==ncov || exclude[i][j]==0){
         criterion=padders2(i,j,current,new,newt,criterion,data,mind,lins);
      }
   }

/* copy the result if there is success */
   if(criterion>0.){
      pswapspace(current,new,data);
/* announce the result */
      if(silent!=1){
         i=(*current).nbas-1;
         puuu(current,(*current).basis[i].b1,(*current).basis[i].b2,
            (*current).basis[i].t1,(*current).basis[i].t2,ncov,0);
         (void)Rprintf("(rao=%.2f)\n",criterion);
      }
      return 1;
   }
/* failure */
   else return 0;
}

/******************************************************************************/

/* this routine pearches a subdimension for a supspace to add */

static double padders1(i0,j0,new,newt,crit,data)
int i0,j0;
struct space *new,*newt;
struct datastruct *data;
double crit;

/* i0,j0     - which subspace (see pstruct)
   new       - will be the best space with additions up to now.
   newt      - actually a copy of current, we play with it until we are done
   current   - sometimes we need two of them (see newt)
   crit      - the best rao statistic (chi-square p-value) up to now
   data      - structure containing the data */

{
   int i,j,d1,d2,d3;

/* pswapspace- copies one space into another
   i,j       - counter
   pestbasis - does the work for 2d dimensions and 1d where no pearch needed
   critx     - possibly optimal criterion 
   d1,d2,d3  - save typing */


   d1=(*newt).sub[i0][j0].dim1;
   d2=(*newt).sub[i0][(*data).ncov].dim1;
   d3=(*newt).sub[j0][(*data).ncov].dim1;
/* a 1-d space */ 
   if(j0==(*data).ncov){
/* a covariate that has not yet been entered */
      if(d1==0) crit=pestbasis(new,newt,crit,data,i0,j0,0,-1,(double)0);

/* a covariate that has been entered before */
      else return crit;
   }

/* a 2-d space */
   else{

/* linear x linear */
      if(d2>0 && d3>0 ){
         if(d1==0){
            crit=pestbasis(new,newt,crit,data,i0,j0,-1,-1,(double)0);
         }
         else{
             
            for(i=0;i<d2-1;i++){
               if((*newt).sub[i0][j0].kts1[i+1][0]>0){
/* knot x knot */
                  for(j=0;j<d3-1;j++){
                     if((*newt).sub[i0][j0].kts1[i+1][j+1]==0 &&
                        (*newt).sub[i0][j0].kts1[0][j+1]>0){
                        crit=pestbasis(
                           new,newt,crit,data,i0,j0,i,j,(double)0);
                     }
                  }
               }
               else{
/* knot x linear */
                  crit=pestbasis(new,newt,crit,data,i0,j0,i,-1,(double)0);
               }
            }
            for(j=0;j<d3-1;j++){
               if((*newt).sub[i0][j0].kts1[0][j+1]==0){
/* linear x knot */
                  crit=pestbasis(new,newt,crit,data,i0,j0,-1,j,(double)0);
               }
            }
         }
      }
   }
   return crit;
}
/******************************************************************************/

/* this routine pearches a subdimension for a supspace to add */

static double padders2(i0,j0,current,new,newt,crit,data,mind,lins)
int i0,j0,mind,*lins;
struct space *new,*newt,*current;
struct datastruct *data;
double crit;

/* i0,j0     - which subspace (see pstruct)
   new       - will be the best space with additions up to now.
   newt      - actually a copy of current, we play with it until we are done
   current   - sometimes we need two of them (see newt)
   crit      - the best rao statistic (chi-square p-value) up to now
   data      - structure containing the data 
   mind      - minimum distance (in order statistics) between knots 
   lins      - which dimensions can only be added linear */

{
   int d1;
   double critx;

/* pswapspace- copies one space into another
   pearch    - does the work for 1d dimensions
   critx     - possibly optimal criterion 
   d1        - save typing */

   d1=(*newt).sub[i0][j0].dim1;
/* a 1-d space */ 
   if(j0==(*data).ncov){
/* a covariate that has not yet been entered */
      if(d1==0) return crit;

/* a covariate that has been entered before */
      else if((*data).bincov[i0]==0 && d1<MAXKNOTS){
         if(lins[i0]==0) critx=pearch(current,newt,data,i0,mind,crit);
         else critx= -100.;
         if(critx>crit){
            crit=critx;
            pswapspace(new,current,data);
         }
      }
   }
   return crit;
}
/******************************************************************************/

/* if a new knot is to be added in a one-covariate dimension or in time, we 
   have to pearch, and that is what we do in this routine */

static double pearch(new,newt,data,i0,mind,crito)
struct space *newt,*new;
double crito;
struct datastruct *data;
int i0,mind;

/* new   - the best added space up to now
   newt  - a space to which we can add
   data  - data
   i0    - first coordinate of the subdimension (second is data.ncov)
   mind  - minimum distance (in order statistics) between knots 
   crito - old criterion */

{
   double *sorted,critnew,crit,crit2,*kts;
   int i,lgth,iloc,lloc,bloc,uloc,iloc2,ll,uu,nx,l;

/* sorted  - sorted data or covariate
   critnew - new criterion
   crit    - best criterion up to now
   crit2   - alternate new criterion
   pestbasis - compute criterion for a basis
   kts     - already used knots
   i       - counter
   lgth    - number of already used knots
   iloc    - present location under study
   lloc    - lower bound to best location
   bloc    - best location up to now
   uloc    - upper bound to best location
   iloc2   - other location under study
   ll      - candidate for lloc
   uu      - candidate for uloc 
   nx      - (*data).ndata 
   l       - emergency break
   Ppsort  - sorting routine 
   pind..  - pind location for new knot under various circumstances */
   
/* initialization */
   bloc = -1;
   crit = -pow((double)10.,(double)20.);
   sorted = v4;

/* find lgth, create kts: already used knots */
   lgth = (*newt).sub[i0][(*data).ncov].dim1-1;
   kts = v3;
   for(i=0;i<lgth;i++) kts[i]=(*newt).sub[i0][(*data).ncov].ktsc[i];
   nx = (*data).ndata;
   for(i=0;i<nx;i++) sorted[i]=trcov[i0+(*data).icov[i]-1];
   Ppsort(sorted,nx);

/* find the interval */
   for(i= -2; i<=lgth;i++){
      if(lgth>0  && i== -2)i=0;
/* before first knot */
      if(i== 0   && lgth>0) iloc=pindl(&ll,&uu,mind,sorted,nx,kts[0]);
/* after last knot */
      if(i== lgth&& lgth>0) iloc=pindr(&ll,&uu,mind,sorted,nx,kts[lgth-1]);
/* first knot */
      if(i== 0   && lgth==0)iloc=pindx(&ll,&uu,nx,0,mind);
      if(i== -1  && lgth==0)iloc=pindx(&ll,&uu,nx,1,mind);
      if(i== -2  && lgth==0)iloc=pindx(&ll,&uu,nx,2,mind);
/* in between knots */
      if(i>0    && i<lgth)  iloc=pindm(&ll,&uu,mind,sorted,nx,kts[i-1],kts[i]);
/* possible location */
      if(iloc>=0){
         critnew=pestbasis(new,newt,crit,data,i0,(*data).ncov,0,0,sorted[iloc]);
/* improvement */
         if(critnew>crit){
            lloc=ll;
            uloc=uu;
            bloc=iloc;
            crit=critnew;
         }
      }
   }
   if(bloc<0)return -1;
   if(crit<crito-25 && lgth==0)return -1;

/* as long as the locations are different, do interval halving */
   l= -1;
   do{
      l++;
      if(l>=1 && crit<crito-25)return -1;
      if(l>=3 && crit<crito-10)return -1;
      if(sorted[uloc]>sorted[lloc]){
         iloc2=pindyr(uloc,bloc,sorted);
/* two pearch points, the upper one */
         if(iloc2>=0){
            crit2=pestbasis(new,newt,crit,data,i0,(*data).ncov,0,
                                                    0,sorted[iloc2]);
         }
         else crit2=crit;

/* two pearch points, the lower one */
         iloc=pindyl(bloc,lloc,sorted);
         if(iloc>=0){
            critnew=pestbasis(new,newt,crit2,data,i0,(*data).ncov,0,
                                                  0,sorted[iloc]);
         }
         else critnew=crit;
/* the middle one is the best, we call it quits */
         if(crit>=critnew && crit>=crit2){
            lloc=uloc;
         }
         else{
/* the lower pearch point is the best */
            if(critnew>crit2){
               uloc=bloc;
               bloc=iloc;
               crit=critnew;
            }
            else{
/* the upper pearch point is the best */
               lloc=bloc;
               bloc=iloc2;
               crit=crit2;
            }
         }
      }
   }while(sorted[uloc]>sorted[lloc]);
   return crit;
}
/******************************************************************************/

/* after another routine has decided to check the rao-criterion for a model
   with an added basis, this routine first adds the basis (addbasis), then
   it checks the criterion (cripswap) - there are lots of possibilities to
   check. */

static double pestbasis(new,newt,criterion,data,i0,j0,ki,kj,ti)
double ti,criterion;
int i0,j0,ki,kj;
struct datastruct *data;
struct space *new,*newt;

/* new       - best space with added dimensions
   newt      - space to which dimensions are added 
   data      - data
   criterion - best rao statistic up to now
   i0,j0     - indicate which subdimension is going to be changed
   ki,kj     - ranknumber of knots to be added
   ti        - some sort of knot to be added */

{

   double arg[4];
   int ncov=(*data).ncov;

/* cripswap  - computes rao and if there is improvement swaps the space
   ncov      - save typing */

/* most common occurences - preset for linear in covariates */
   arg[0]= -1.;
   arg[1]= -1.;
   arg[2]= -1.;
   arg[3]= -1.;

/* 1 covariate subdimension */
/* this is not the first (i.e. linear) space */
   if(j0==ncov) if((*newt).sub[i0][j0].dim1>0){
/* what is the knot to be added */
      arg[0]=ti;
      arg[2]=(*newt).sub[i0][ncov].dim1-1;
      (*newt).sub[i0][ncov].ktsc[(*newt).sub[i0][ncov].dim1-1] =ti;
   }

/* a crossproduct subdimension */
   if(j0<ncov){
/* if it is the first one, it is linear*time-knot[ki] */
      if(ki>=0){
         arg[2]=ki;
         arg[0]=(*newt).sub[j0][ncov].ktsc[(int)arg[2]];
      }
      if(kj>=0){
         arg[3]=kj;
         arg[1]=(*newt).sub[j0][ncov].ktsc[(int)arg[3]];
      }
      (*newt).sub[i0][j0].kts1[ki+1][kj+1]=1;
   }

   paddbasis(i0,j0,arg,data,&((*newt).basis[(*newt).nbas]));
/* compute rao. possibly swap */
   criterion=cripswap(newt,data,new,criterion,i0,j0);
   if(j0<ncov)(*newt).sub[i0][j0].kts1[ki+1][kj+1]=0;

/* done */
   return criterion;
}
/******************************************************************************/

/* after the space newt has been updated with extra new basisfunctions, it
   computes the criterion. If this is an improvement it copies the basis into
   new, then it restores newt. */

static double cripswap(newt,data,new,criterion,i0,j0)
double criterion;
struct datastruct *data;
int i0,j0;
struct space *new,*newt;

/* criterion - best rao p-value up to now
   data      - the data
   i0,j0     - which subdimension is being altered
   new       - best model with additions
   newt      - model tested whether it is better */

{
   double crit;

/* crit      - present value of criterion
   prao     - computes rao-criterion
   pswapspace - copies one space into another */

/* update the dimension parameters */

   ((*newt).ndim)+=(*data).nclass;
   ((*newt).nbas)+=1;
   ((*newt).sub[i0][j0].dim1)+=1;

   crit=prao(newt,data);

/* if there is improvement, copy the space */
   if(crit>criterion){
      pswapspace(new,newt,data);
      criterion=crit;
   }

/* change back the dimensions */

   ((*newt).ndim) -= (*data).nclass;
   ((*newt).nbas) -= 1;
   ((*newt).sub[i0][j0].dim1) -= 1;

   return criterion;  
}
/******************************************************************************/

/* this routine computes the extra elements of hessian and score then it
   computes rao - the routine is very much like the routine complog, which
   is part of Newton, except that it does not compute the log-likelihood 
   and it makes use of the fact that part of b0, b1 and b2 might be known and
   completely at the end, it computes the rao statistic.                      */

static double prao(spc,data)
struct space *spc;
struct datastruct *data;

/* spc   - the present model 
   data  - the data */
{
   double raoc=0,**hhh,*ss,*ss2,rtemp,xx,epsi=(*spc).epsilon,yy,yy4,d0,d1,**hh2;
   int i,j,k,k2,itemp,nclass=(*data).nclass,nbas=(*spc).nbas,k3;
   int ndim=(*spc).ndim,extra=nbas-1,alhere=(nbas-1)*nclass,ncl;
   double *wh1,*wh2,*wh3;

/* i,j,k,k2   - counters
   ss,ss2  - score
   raoc    - rao
   hhh,hh2     - hessian 
   nclass, nbas, ndim, extra, alhere - save typing
   rtmep, itemp - frequently used 
   xx,epsi - save typing
   yy, yy4, d0, d1, k3 - half products
   wh1, wh2, wh3  - half products */

/* allocation  and initlization equal to 0 */
   xx=epsi*2./(nclass+1.);
   ncl=nclass+1;
   hhh=w1;
   hh2=w2;
   ss=v1;
   ss2=v2;

/* initialization equal to zero */
   for(i=alhere;i<=ndim;i++){
      ss[i-alhere]=0.;
      ss2[i-alhere]=0.;
      for(j=0;j<=i;j++){
         hhh[i][j]=0.;
         hhh[j][i]=0.;
      }
   }

/* now circle the datapoints */
   for(i=0;i<(*data).ndata;i++){
      wh1=(*data).work[i];
      wh2=(*data).work2[i];
      petvector(spc,data,v7,v8,i);
      yy= -xx*(*data).wgt[i];
      for(k=0;k<nclass;k++){
         rtemp=v8[extra]*wh1[k];
         ss[k] -= rtemp;
         itemp=alhere+k;
         wh3=hhh[itemp];
         yy4=yy*v7[extra];
         for(k2=0;k2<nclass;k2++) wh3[k2] += yy4-rtemp*wh1[k2];
         wh3[k]+=rtemp-ncl*yy4;
         for(j=1;j<nbas;j++){
            k3=j*nclass;
            if(v7[j]!=0.0){
               d0=rtemp*v7[j];
               d1=yy4*v7[j];
               for(k2=0;k2<nclass;k2++) wh3[k2+k3] += d1-d0*wh1[k2];
               wh3[k+k3]+=v7[j]*rtemp-ncl*yy4*v7[j];
            }
         }
         for(k2=0;k2<=nclass;k2++) ss[k] += wh2[k2]*v7[extra];
         ss[k] -= wh2[k]*v7[extra]*ncl;
      }
      if((*data).yy[i]!=nclass) ss[(*data).yy[i]]+= v8[extra];
   }
   
/* compute rao */
   for(k=0;k<nclass;k++){
      itemp=alhere+k;
      for(j=0;j<nbas-1;j++){
         for(k2=0;k2<nclass;k2++){
            hhh[k2+j*nclass][itemp]=0.;
            for(i=0;i<nbas-1;i++)for(k3=0;k3<nclass;k3++)
               hhh[k2+j*nclass][itemp]+=
               hhh[k2+j*nclass][k3+i*nclass]*hhh[itemp][k3+i*nclass];
         }
      }
   }
   raoc=0;
   for(k=0;k<nclass;k++){
      ss2[k]=ss[k];
      itemp=alhere+k;
      for(k2=0;k2<nclass;k2++){
         j=alhere+k2;
         for(i=0;i<nbas-1;i++)for(k3=0;k3<nclass;k3++)hhh[j][itemp]-=
            hhh[j][k3+i*nclass]*hhh[k3+i*nclass][itemp];
         hh2[k2][k]=hhh[j][itemp];
      }  
   }
   j=lusolinv(hh2,nclass,ss2,2);
   if(j>0) for(j=0;j<nclass;j++) raoc+=ss2[j]*ss[j];
   if(raoc>1000 && nbas>3)raoc=0.;
   return raoc;
}
/******************************************************************************/
/* pinds a new location in an interval (l,b) - that is the lower end might not
   have been tested yet */
static int pindyl(u,l,x)
int l,u;
double *x;
{
   int i;
   if(x[l]==x[u])return -1;
   i=(u+l-1)/2;
   if(x[i]!=x[u])return i;
   i=(i+l)/2;
   if(x[i]!=x[u])return i;
   return l;
}
/******************************************************************************/
/* pinds a new location in an interval (l,u) - that is the upper end might not
   have been tested yet */
static int pindyr(u,l,x)
int l,u;
double *x;
{
   int i;
   if(x[l]==x[u])return -1;
   i=(u+l+1)/2;
   if(x[i]!=x[l])return i;
   i=(i+u)/2;
   if(x[i]!=x[l])return i;
   return u;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (0,knot1) */
static int pindl(ll,uu,mind,x,nx,knt)
double *x,knt;
int nx,*ll,*uu,mind;
/* ll - lowest number we can pearch on in the future
   uu - highest number we can pearch on in the future
   mind minimum distance between knots
   x  - data
   nx - length of data
   knt- knot */
{

/* i  - utility
   plocation - pinds uu */

   int i;

   (*uu)=plocation(0,x,nx,knt);
   if((*uu)<2*mind)return -1;
   i=((*uu)-1)/2;
   if((*uu)-i<mind+1)i=(*uu)-mind-1;
   *ll=mind;
   *uu=(*uu)-mind-1;
   return i;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (knot-last,nx-1) */

static int pindr(ll,uu,mind,x,nx,knt)
double *x,knt;
int nx,*ll,*uu,mind;
/* ll - lowest number we can pearch on in the future
   uu - highest number we can pearch on in the future
   mind minimum distance between knots
   x  - data
   nx - length of data
   knt- knot */
{
   int i;
/* i  - utility
   plocation - pinds ll */

   (*ll)=plocation(1,x,nx,knt);
   if(nx-1-(*ll)<2*mind)return -1;
   i=(nx+(*ll))/2;
   if(i-(*ll)<mind+1)i=(*ll)+mind+1;
   *uu=nx-1-mind;
   *ll=(*ll)+mind+1;
   return i;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (0,nx-1) */

static int pindx(ll,uu,nx,i,mind)
int nx,*ll,*uu,mind,i;
/* ll - lowest number we can pearch on in the future
   uu - highest number we can pearch on in the future
   nx - length of data */
{
   if(i==0){
      *ll=mind;
      *uu=nx/2;
      if((*uu)>nx-mind-1)(*uu)=nx-mind-1;
   }
   if(i==1){
      *ll=nx/4;
      if((*ll)>mind)(*ll)=mind;
      *uu=3*nx/4;
      if((*uu)>nx-mind-1)(*uu)=nx-mind-1;
   }
   if(i==2){
      *ll=nx/2;
      if((*ll)>mind)(*ll)=mind;
      *uu=nx-1-mind;
   }
   if((*ll)>(*uu))return -1;
   return (int)((*ll)+(*uu))/2;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (k0,k1) */

static int pindm(ll,uu,mind,x,nx,k0,k1)
double *x,k0,k1;
int nx,*ll,*uu,mind;
/* ll - lowest number we can pearch on in the future
   uu - highest number we can pearch on in the future
   mind minimum distance between knots
   x  - data
   nx - length of data
   k0 - knot
   k1 - knot */

{
/* plocation - pinds ll */


   (*ll)=plocation(1,x,nx,k0);
   (*uu)=plocation(0,x,nx,k1);
   if((*uu)-(*ll)<2*mind+1)return -1;
   *uu=(*uu)-mind-1;
   *ll=(*ll)+mind+1;
   return ((*uu)+(*ll))/2;
}
/******************************************************************************/
/* finds the lowest (if what = 0) or the highest (if what = 1) index of x for
   which x==k */

static int plocation(what,x,nx,k)
int nx,what;
double k,*x;

/* what - see above
   x    - data
   nx   - length data
   k    - see above */

{
   int i;
   if(what==1){
      if(x[0]>k)return 0;
      if(x[nx-1]<=k)return nx-1;
      for(i=0;i<nx-1;i++){
         if(x[i+1]>k && x[i]<=k) return i;
      }
   }
   if(x[nx-1]<k)return nx-1;
   if(x[0]>=k)return 0;
   for(i=1;i<nx;i++){
      if(x[i]>=k && x[i-1]<k)return i;
   }
   return nx;
}
/******************************************************************************/

/* after another routine has decided to add a basis function, this routine
   actually adds the basis functions */

static void paddbasis(i0,j0,arg,data,basis)
double *arg;
struct datastruct *data;
struct basisfunct *basis;
int i0,j0;

/* i0,j0 - which subdimension does this new basisfunction belong to
   arg   - elements 2 and 3: ranknumber of the knot
           elements 0 and 1: the knot itself.
   data  - the data
   basis - the basisfunction which is added */

{
   int i;

/* i - counter
   pct1,pct2 - functions to compute version */

/* update the scalar part of the basis function */
   for(i=0;i<(*data).nclass;i++)(*basis).beta[i]=0.;
   (*basis).b1=i0;
   (*basis).b2=j0;
   (*basis).t1=arg[2];
   (*basis).t2=arg[3];
}
/******************************************************************************/

/* the S-I/O routine */
void spolyx(intpars)
int *intpars;
{
   intpars[1]=MAXSPACE;
   intpars[0]=MAXKNOTS;
   return;
}

void spoly(intpars,cls,t1cov,iloss,penalty,bbtt,cckk,vexclude,lins,logs,wgt,
           tcls,t2cov,twgt,bbb,xxx)
double *penalty,*iloss,*logs,*bbtt,*cckk,*wgt,*twgt,*bbb;
int *lins,*vexclude,*intpars,*cls,*tcls;
float *xxx,*t1cov,*t2cov;

/* intpars- all sorts of integer parameters
   cls    - the classes
   trcov   - the covariates
   iloss   - loss function 
   penalty- penalty 
   bbtt   - summarizes the model (basisfunctions)
   cckk   - summarizes the model (1d subspaces) 
   vexclude - exclude in vector form
   lins   - which terms should be linear
   logs/ad- record information about the models
   wgt    - case-weights
   tcls,tecov,twgt - as cls, trcov and wgt, but for the test-data 
   aics - losses for models */
{
   int i,j,k,**exclude,mindis,ndmax,silent,fitter,strt,it,il;
   int cv,cvx,j0,j1,bb1,bb2,*numbers,naction,sngle,*ad;
   double **ranges,**losses,*aics,**meas;
   struct space *best,*s1,*s2,*s3,*s4;
   struct datastruct *data,*tdata,*zdata;
   double **loss;

/* i,k,l         - counters 
   exclude       - which terms to exclude
   mindis        - minimum distance between knots
   ndmax         - maximum number of dimensions
   naction       - number of possible actions
   silent        - should diagnostic output be printed? (0=yes, 1=no) 
   fitter        - was a starting model provided?
   strt          - start with constant (0) ot linear (1) 
   zdata         - data storage for cv
   it            - 1: testset, 0: AIC, 2: cv
   cv            - number of cross validations
   cvx           - present cross validations
   j0            - counter of cv-cases
   j1            - counter of cv-test cases
   bb1           - lowest number that is a cv-test case
   bb2           - highest number +1 that is a cv-test case
   numbers       - numbers[i] has the nr of different models for CV nr i
   aicbest       - finds out for which aic-values a particular model is best
   ranges        - ranges[i][j] for CV nr i the j-th aic break point
   losses        - losses[i][j] the loss between ranges[i][j] & ranges[i][j+1]
   aiccv         - finds the aic cv-values
   best          - the final fitted model
   s1,s2,s3      - storage for spaces
   pdefinespace- allocates storage for a space
   predefinespace- initializes a space
   data          - structure containing all information about the data
   tdata         - structure containing all information about the testset
   pdefinedata - allocates storage for data
   poly        - that is where it happens
   soutspace   - arranges the output 
   freespace   - cleans up the bits
   initbasis   - initializes one basis function */

/* intpars 0 = ndata/1 = nclass/2 = ncov/3 = mindis/4 = ndmax/5 = silent/
   /6 = fitter/7 = cv/8 = it/9 = ndata (testset)/10 = naction/11 = il */

/* if intpars[0]<0 we only want the parameters */
   trcov=t1cov;
   tecov=t2cov;
   if(intpars[0]<0){
      intpars[1]=MAXSPACE;
      intpars[0]=MAXKNOTS;
      return;
   }
/* allocations of workspace */
   maxdim=intpars[4];
   if(maxdim<0)maxdim= -maxdim;
   iw1=ipmatrix(maxdim+1,maxdim+1);
   w1=dpmatrix(maxdim+1,maxdim+1);
   w2=dpmatrix(maxdim+1,maxdim+1);
   w3=dpmatrix(maxdim+1,maxdim+1);
   meas=dpmatrix(7,maxdim+1);
   iv1=ispvector(maxdim+1);
   iv2=ispvector(maxdim+1);
   v1=dspvector(maxdim+1);
   v2=dspvector(maxdim+1);
   v3=dspvector(maxdim+1);
   v4=dspvector(intpars[9]+intpars[0]+3);
   v5=dspvector(maxdim+1);
   v6=dspvector(maxdim+1);
   v7=dspvector(maxdim+1);
   v8=dspvector(maxdim+1);
   aics=dspvector(maxdim+1);
   ad=ispvector(maxdim+1);
   for(i=0;i<=maxdim;i++)ad[i]= -1;

/* get some parameters */
   mindis=intpars[3];
   ndmax=intpars[4];
   silent=intpars[5];
   fitter=intpars[6];
   cv=intpars[7];
   it=intpars[8];
   naction=intpars[10];
   il=intpars[11];
   sngle=intpars[12];
   losses=dpmatrix(cv+2,maxdim+1);

/* arrange the data  the test data and the loss matrix */
   data=pdefinedata(intpars[0],intpars[2],intpars[1],intpars[0],cls,wgt,1);
   if(it==1){
/* a test set and a loss matrix */
      tdata=pdefinedata(intpars[9],intpars[2],intpars[1],intpars[9],
            tcls,twgt,-1);
      loss=dpmatrix(intpars[10],intpars[1]);
      for(i=0;i<intpars[10];i++){
         for(j=0;j<intpars[1];j++) loss[i][j]=iloss[i+intpars[10]*j];
      }
   }
   if(it==0) {
/* the test set is fake */
      tdata=pdefinedata(2,intpars[2],intpars[1],0,tcls,twgt,-1);
      loss=dpmatrix(intpars[10],intpars[1]);
      for(i=0;i<intpars[10];i++){
         for(j=0;j<intpars[1];j++) loss[i][j]=iloss[i+intpars[10]*j];
      }
   }
   if(it==2) {
/* cv - zdata keeps the data safe - also, a loss function */
      tdata=pdefinedata(intpars[0],intpars[2],intpars[1],intpars[0],
            cls,wgt,1);
      zdata=pdefinedata(intpars[0],intpars[2],intpars[1],intpars[0],
            cls,wgt,1);
      loss=dpmatrix(intpars[10],intpars[1]);
      for(i=0;i<intpars[10];i++){
         for(j=0;j<intpars[1];j++) loss[i][j]=iloss[i+intpars[10]*j];
      }
      numbers=ispvector(cv);
      ranges=dpmatrix(cv,maxdim+1);
   }

/* check whether things are binary */
   for(k=0;k<(*data).ncov;k++){
      (*data).bincov[k]=1;
      for(i=0;i<(*data).ndata;i++){
         if(trcov[(*data).icov[i]+k-1]!=0. && trcov[(*data).icov[i]+k-1]!=1.){
            if(it==2) (*zdata).bincov[k]=0;
            (*data).bincov[k]=0;
            i=(*data).ndata;
         }
      }
   }
   if(it==2)for(k=0;k<(*data).ncov;k++)(*zdata).bincov[k]=(*data).bincov[k];

/* initialize exclude */
   exclude=ipmatrix((*data).ncov+1,(*data).ncov+1);
   k=(int)fabs((double)vexclude[0]);
   for(i=0;i<k;i++){
      if(vexclude[2*i+1]==0)vexclude[2*i+1]=(*data).ncov+1;
      if(vexclude[2*i+2]==0)vexclude[2*i+2]=(*data).ncov+1;
      exclude[vexclude[2*i+1]-1][vexclude[2*i+2]-1]=1;
      exclude[vexclude[2*i+2]-1][vexclude[2*i+1]-1]=1;
   }

/* if the 1st element of vexclude was < 0 it was vinclude */
   if(vexclude[0]<0){
      for(k=0;k<=(*data).ncov;k++){
         for(i=0;i<=(*data).ncov;i++)exclude[i][k]=1-exclude[i][k];
      }
   }

/* make a space */
   best=pdefinespace(data);
   s1=pdefinespace(data);
   s2=pdefinespace(data);
   s3=pdefinespace(data);
   s4=pdefinespace(data);

   for(cvx=0;cvx<=cv;cvx++){
      if(cv!=0 && silent!=1){
         if(cv==cvx)(void)Rprintf("Final run starts now.\n");
         else (void)Rprintf("CV-run %d starts now.\n",cvx+1);
         (void)fflush(stdout);
      }
/* get ad ready */
      if(cvx!=0){
         predefinespace(zdata,best);
         predefinespace(zdata,s1);
         predefinespace(zdata,s2);
         predefinespace(zdata,s3);
         predefinespace(zdata,s4);
      }
      for(i=0;i<maxdim;i++)ad[i]= -1;
      if(it==2){
/* is this the final cv-run ? */
         if(cvx==cv){
/* find the optimal cv-value */
            (*penalty)=aiccv(ranges,losses,numbers,cv,xxx,il);
            if(silent!=1){ 
               (void)Rprintf("CV-Penalty range is %5.2f",xxx[0]);
               if(xxx[1]<0)(void)Rprintf(" - Inf");
               else (void)Rprintf("- %5.2f",xxx[1]);
               (void)Rprintf(" with loss %.2f.\n",xxx[2]);
               (void)Rprintf("Penalty value to be used is %5.2f\n",xxx[3]);
               (void)fflush(stdout);
            }
            data=zdata;
            it=0;
         }
         else{
/* separate test set and training set */
            bb1=((*zdata).ndata*cvx)/cv;
            bb2=((*zdata).ndata*(cvx+1))/cv;
            (*data).ndata=(*zdata).ndata+bb1-bb2;
            (*tdata).ndata=bb2-bb1;
            j0=0;
            j1=0;
            (*tdata).wgtsum=0.;
            (*data).wgtsum=0.;
            for(i=0;i<(*zdata).ndata;i++){
               if(i>=bb1 && i<bb2){
                  (*tdata).wgt[j1]=(*zdata).wgt[i];
                  (*tdata).yy[j1]=(*zdata).yy[i];
                  (*tdata).wgtsum+=(*tdata).wgt[j1];
                  (*tdata).icov[j1]=(*zdata).icov[i];
                  j1++;
               }
               else{
                  (*data).wgt[j0]=(*zdata).wgt[i];
                  (*data).yy[j0]=(*zdata).yy[i];
                  (*data).wgtsum+=(*data).wgt[j0];
                  (*data).icov[j0]=(*zdata).icov[i];
                  j0++;
               }
            }
         }
      }
/* we should initialize the space if we started with a fit */
      strt=0;
      if(fitter>0){
         (*best).nbas=fitter/(*data).nclass;
         (*best).ndim=fitter;
         j=0;
/* record the knots from cckk */
         for(i=0;i<(*data).ncov;i++){
            j++;
            (*best).sub[i][(*data).ncov].dim1=cckk[j];
            for(k=1;k<cckk[0];k++){
               j++;
               (*best).sub[i][(*data).ncov].ktsc[k-1]=cckk[j];
            }
         }
/* record the basisfunctions */
         j=0;
         for(i=0;i<(*best).nbas;i++){
            (*best).basis[i].b1=bbtt[j]-1; 
            if((*best).basis[i].b1== -2)(*best).basis[i].b1= (*data).ncov;
            j++;
            (*best).basis[i].t1=bbtt[j]-1; 
            if((*best).basis[i].t1== -2)(*best).basis[i].t1= -1;
            j++;
            (*best).basis[i].b2=bbtt[j]-1; 
            if((*best).basis[i].b2== -2)(*best).basis[i].b2= (*data).ncov;
            j++;
            (*best).basis[i].t2=bbtt[j]-1; 
            if((*best).basis[i].t2== -2)(*best).basis[i].t2= -1;
            j++;
            for(k=0;k<(*data).nclass;k++){
               (*best).basis[i].beta[k]=bbtt[j];
               j++;
            }
            (*best).basis[i].ib=(*data).nclass;
            for(k=0;k<(*data).nclass;k++)(*best).basis[i].link1[k]=k;
            if((*best).basis[i].b2!=(*data).ncov){
               ((*best).sub[(*best).basis[i].b1][(*best).basis[i].b2].dim1)++;
               (*best).sub[(*best).basis[i].b1][(*best).basis[i].b2].
                  kts1[(*best).basis[i].t1+1][(*best).basis[i].t2+1]=1;
            }
            else{
               if(i>0 && (*best).basis[i].t1== -1){
                  ((*best).sub
                     [(*best).basis[i].b1][(*best).basis[i].b2].dim1)++;
               }
            }
         }
         strt = -1;
      }
   
/* do the work */
      poly(best,data,loss,-(*penalty),ndmax,mindis,exclude,strt,silent,meas,ad,
         lins,tdata,it,losses[cvx],s1,s2,s3,naction,il,sngle,s4);
/* organize for cv */
      if(it==2) numbers[cvx]=aicbest(ad,ranges[cvx],losses[cvx],meas[0]);
      else for(i=0;i<maxdim;i++)aics[i]=losses[cvx][i];
   }
   aicb2(ad,aics,meas,logs);
/* admire the results */
   proj(best,data,bbb,iloss);
   soutspace(best,data,bbtt,cckk);
   intpars[0]=(*best).ndim; 
   intpars[1]=(*data).nclass; 
   intpars[2]=(*best).nbas; 
   return;
}
/******************************************************************************/

/* this is an output routine, it writes the matrices bbtt and cckk - which are
   given as output to S */

static void soutspace(spc,data,bbtt,cckk)
struct space *spc;
struct datastruct *data;
double *cckk,*bbtt;

/* spc   - structure describing the model
   data  - data 
   bbtt  - matrix describing the basisfunctions
   cckk  - matrix describing the 1d subspaces  */

{
   int i,j,k,l;

   l=MAXKNOTS+1;

/* cckk lines - covariate number, dimension of space (#kts+1), kts */
   for(j=0;j<(*data).ncov;j++){
      cckk[j*l]=(*spc).sub[j][(*data).ncov].dim1-1;
      if((*spc).sub[j][(*data).ncov].dim1==0)cckk[j*l]=0;
      for(i=0;i<(*spc).sub[j][(*data).ncov].dim1-1;i++){
         cckk[j*l+i+1]=(*spc).sub[j][(*data).ncov].ktsc[i];
      }
      for(k=(*spc).sub[j][(*data).ncov].dim1;k<=MAXKNOTS;k++)cckk[(j+1)*l+k]=0.;
   }

/* rest: basisfunctions  (bbtt) */
   for(j=0;j<(*spc).nbas;j++){
      if((*spc).basis[j].b1>=0){
         (*spc).basis[j].b1+=1;
      }
      if((*spc).basis[j].b2>=0){
         (*spc).basis[j].b2+=1;
         if((*spc).basis[j].b2>(*data).ncov) (*spc).basis[j].b2= -1;
      }
      (*spc).basis[j].t1+=1;
      (*spc).basis[j].t2+=1;

/* which variable, which knot, which variable#2, which knot #2, beta 
   variable=0: time, knot=0: constant, otherwise: knot number
   variable>0: covariate, knot=0: linear, otherwise: knot number 
   variable#2 = -1: 1d basisfunction */
      bbtt[j*(4+(*data).nclass)+0]=(*spc).basis[j].b1;
      bbtt[j*(4+(*data).nclass)+1]=(*spc).basis[j].t1;
      bbtt[j*(4+(*data).nclass)+2]=(*spc).basis[j].b2;
      bbtt[j*(4+(*data).nclass)+3]=(*spc).basis[j].t2;
      for(k=0;k<(*data).nclass;k++)
         bbtt[j*(4+(*data).nclass)+4+k]=(*spc).basis[j].beta[k];
   }
   return;
}

/******************************************************************************/
/* gets the basisfunctions */
static double petvector2(best,data,i,k)
struct space *best;
struct datastruct *data;
int i,k;

/* best  - the model
   data - the data */
{
   int j,b1,b2,t1,t2,ncov=(*data).ncov;
   double xx,val;
   float *cc;

/* j   - counter
   b1,b2,t1,t2 - b1,b2,t1,t2 for the present basisfunction
   xx    - the second half 
   ndata - number of datapoints
   ncov  - number of covariates
   cov   - covariates */

   b1=(*best).basis[i].b1;
   t1=(*best).basis[i].t1;
   b2=(*best).basis[i].b2;
   t2=(*best).basis[i].t2;
/* circle the basisfunctions */
      j=k;
      k=(*data).icov[j];
      if(k>0)cc= &(trcov[k-1]);
      else cc= &(tecov[-k-1]);
/* if it is time only it is easy */
      val=0.;
      if(b1==ncov){
         val=1.;
      }
      else{

/* then first take the first component of the basisfunction */
         val=cc[b1];

/* -1 means linear, otherwise it is piecewise linear */
         if(t1> -1){
            val-=(*best).sub[b1][ncov].ktsc[t1];
            if(val<0.) val=0.;
         }
      }
/* and then the second component of the basisfunction */
      if(b2!=ncov && b2!= -1){
         xx=cc[b2];

/* -1 means linear, otherwise it is piecewise linear */
         if(t2> -1){
            xx-=(*best).sub[b2][ncov].ktsc[t2];
            if(xx<0.) xx=0.;
         }
         val=val*xx;
      }
   return val;
}
/******************************************************************************/
/* gets the basisfunctions */
static void petvector(best,data,val,wal,j)
struct space *best;
struct datastruct *data;
double *val,*wal;
int j;

/* best  - the model
   data - the data */
{
   int i,k,b1,b2,t1,t2,ncov=(*data).ncov;
   double xx;
   float *cc;

/* i   - counter
   b1,b2,t1,t2 - b1,b2,t1,t2 for the present basisfunction
   xx    - the second half 
   ndata - number of datapoints
   ncov  - number of covariates
   cov   - covariates */

/* circle the basisfunctions */
   k=(*data).icov[j];
   if(k>0)cc= &(trcov[k-1]);
   else cc= &(tecov[-k-1]);
   for(i=0;i<(*best).nbas;i++){

/* if it is time only it is easy */
      val[i]=0.;
      b1=(*best).basis[i].b1;
      if(b1==ncov){
         val[i]=1.;
      }
      else{
         t1=(*best).basis[i].t1;

/* then first take the first component of the basisfunction */
         val[i]=cc[b1];

/* -1 means linear, otherwise it is piecewise linear */
         if(t1> -1){
            val[i]-=(*best).sub[b1][ncov].ktsc[t1];
            if(val[i]<0.) val[i]=0.;
         }
      }
/* and then the second component of the basisfunction */
      b2=(*best).basis[i].b2;
      if(b2!=ncov && b2!= -1){
         xx=cc[b2];
         t2=(*best).basis[i].t2;

/* -1 means linear, otherwise it is piecewise linear */
         if(t2> -1){
            xx-=(*best).sub[b2][ncov].ktsc[t2];
            if(xx<0.) xx=0.;
         }
         val[i]=val[i]*xx;
      }
      wal[i]=val[i]*((*data).wgt[j]);
   }
   return;
}
/******************************************************************************/
/* this function allocates storage for a data structure */

static struct datastruct *pdefinedata(ndata,ncov,nclass,xndata,cls,wgt,icov)
int ncov,ndata,nclass,xndata,icov;
int *cls;
double *wgt;
{
   struct datastruct *newdata;
   int i;
   newdata=(struct datastruct *)Salloc(1,struct datastruct);
   (*newdata).work=dpmatrix(ndata,nclass+1);
   (*newdata).work2=dpmatrix(ndata,nclass+2);
   (*newdata).bincov=ispvector(ncov);
   (*newdata).wgt=dspvector(ndata);
   (*newdata).yy=ispvector(ndata);
   (*newdata).icov=ispvector(ndata);
   (*newdata).ndata=xndata;
   (*newdata).nclass=nclass-1;
   (*newdata).ncov=ncov;
   (*newdata).wgtsum=0.;
   for(i=0;i<xndata;i++){
      (*newdata).yy[i]=cls[i];
      (*newdata).wgt[i]=wgt[i];
      (*newdata).wgtsum+=wgt[i];
      (*newdata).icov[i]=icov*i*ncov+icov;
   }
   return newdata;
}
/******************************************************************************/
/* this function allocates storage for a space */

static struct space *pdefinespace(data)
struct datastruct *data;

{
   int ncov=(*data).ncov,nclass=(*data).nclass;
   struct space *newspace;
   int i,j;

/* this routine is mainly a copy from the numerical recipes routines
   ncov,ndata,nclass - save typing
   i,j         - counters
   newspace    - thing to be initialized
   pdefinebasis - initializes an array of basisfunctions
   pdefinedim   - initializes a matrix of subdimensions */

/* basic allocation */
   newspace=(struct space *)Salloc(1,struct space);

/* the simple elements */
   (*newspace).info=dpmatrix(maxdim,maxdim);
   (*newspace).infox=dpmatrix(maxdim,maxdim);
   (*newspace).score=dspvector(maxdim);
   (*newspace).epsilon=.00001;
   (*newspace).ndim=0;
   (*newspace).nbas=0;
   (*newspace).aic=0.;
   (*newspace).epsilon=0.;
   (*newspace).logl=0.;

/* defines the basisfunctions, initializes */
   (*newspace).basis=pdefinebasis();
   for(i=0;i<maxdim;i++) {
      (*newspace).basis[i].b1= -1;
      (*newspace).basis[i].b2= -1;
      (*newspace).basis[i].t1= -1;
      (*newspace).basis[i].t2= -1;
      (*newspace).basis[i].ib= nclass;
      (*newspace).basis[i].beta=dspvector(nclass+1);
      (*newspace).basis[i].link1=ispvector(nclass);
      (*newspace).basis[i].link2=ispvector(nclass);
      for(j=0;j<nclass;j++)(*newspace).basis[i].link1[j]=j;
   }

/* defines the subdimensions */
   (*newspace).sub=pdefinedim(ncov+1);

/* for the 2-covariate subdimensions */
   for(i=0;i<ncov;i++) for(j=i+1;j<ncov;j++){
      (*newspace).sub[i][j].kts1=ipmatrix(MAXKNOTS+1,MAXKNOTS+1);
      (*newspace).sub[i][j].dim1=0;
   }

/* for the 1-covariate subdimensions */
   for(j=0;j<ncov;j++){
      (*newspace).sub[j][ncov].ktsc=dspvector(MAXKNOTS);
      (*newspace).sub[j][ncov].dim1=0;
   }

   return newspace;
}

/******************************************************************************/
/* this function allocates storage for an array of basisfunctions */
static struct basisfunct *pdefinebasis()
{
   struct basisfunct *nb;

   nb=(struct basisfunct *)Salloc((maxdim),struct basisfunct);
   return nb;
}

/******************************************************************************/
/* this function allocates storage for a matrix of subdimensions */
static struct subdim **pdefinedim(ncov)
int ncov;
{
   struct subdim **newdim;
   int i;

   newdim=(struct subdim **)Salloc(ncov,struct subdim *);
   for(i=0;i<ncov;i++)
      newdim[i]=(struct subdim *)Salloc(ncov,struct subdim);
   return newdim;
}

/******************************************************************************/
static int *ispvector(l)
int l;
/* allocate an int vector with subscript range v[0...l] */
{
   int *v,i;
   v=(int *)Salloc((l+1),int);
   for(i=0;i<=l;i++)v[i]=0;
   return v;
}
/******************************************************************************/
static double *dspvector(l)
int l;
/* allocate a double vector with subscript range v[0...l] */
{
   double *v;
   int i;
   v=(double *)Salloc((l+1),double);
   for(i=0;i<=l;i++) v[i]=0.;
   return v;
}
/******************************************************************************/
static int **ipmatrix(r,c)
int r,c;
/* allocate an int matrix with subscript range m[0..r][0..c] */
{
   int i,j,**m;
   m=(int **) Salloc((r+1),int*);
   for(i=0;i<=r;i++){
       m[i]=(int *) Salloc( (c+1),int);
       for(j=0;j<=c;j++)m[i][j]=0;
   }
   return m;
}
/******************************************************************************/
static double **dpmatrix(r,c)
int r,c;
/* allocate a double matrix with subscript range m[0..r][0..c] */
{
   int i,j;
   double **m;
   m=(double **) Salloc((r+1),double*);
   for(i=0;i<=r;i++){
       m[i]=(double *) Salloc((c+1),double);
       for(j=0;j<=c;j++)m[i][j]=0.;
   }
   return m;
}
/******************************************************************************/
static double aiccv(ranges,losses,numbers,k,xio,il)
int k,*numbers,il;
double **ranges,**losses;
float *xio;
{
   int *where,i,j,l,m;
   double bestl,bestu,bestx=0.,nowl,nowu,nowx,maxu;
   where=iv1;
   for(i=0;i<=maxdim;i++)where[i]=0;
   maxu=0;
   for(i=0;i<k;i++){
      if(ranges[i][numbers[i]-1]>maxu)maxu=ranges[i][numbers[i]-1];
   }
   for(i=0;i<k;i++) ranges[i][numbers[i]]=maxu*2.;
   nowu=0.;
   j=0;
   m=5;
   do{
      nowl=nowu;
      nowu=1.5*maxu;
      l= -1;
      for(i=0;i<k;i++) if(ranges[i][where[i]+1]<nowu){
         l=i;
         nowu=ranges[i][where[i]+1];
      }
      nowx=0.;
      if(il!=0)for(i=0;i<k;i++)nowx+=losses[i][where[i]];
      else for(i=0;i<k;i++)nowx-=losses[i][where[i]];
      if(m<997){
         xio[m]=nowl;
         xio[m+1]=nowu;
         xio[m+2]=nowx;
         m=m+3;
      }
      if(j==0||nowx<=bestx){
         j=1;
         bestx=nowx;
         bestu=nowu;
         bestl=nowl;
      }
      if(l!= -1)where[l]++;
   }while(l != -1);
   xio[0]=bestl;
   xio[1]=bestu;
   xio[2]=bestx;
   if(il==0)xio[2]= -bestx;
   xio[4]=m;
   if(bestu>maxu){
      bestu=maxu;
      bestu=1.0e+30;
      xio[1]= -1.;
   }
   if(bestl<=0)bestl=0.; 
   xio[3]=sqrt(bestl*bestu);
   return xio[3];
}
/******************************************************************************/
static void aicb2(ads,aics,meas,logls)
int *ads;
double *aics,*logls,**meas;
{
   int i,j,k,*k1,*k2;
   double *d1,*d3,*d4,*d2;
   k1=iv1;
   k2=iv2;
   d1=v1;
   d2=v2;
   d3=v3;
   d4=v6;
   j=0;
   for(i=0;i<maxdim;i++) if(ads[i]>=0){
      k1[j]=i+1;
      k2[j]=ads[i];
      d1[j]=meas[0][i];
      d2[j]=aics[i];
      d3[j]= -2.;
      d4[j]= -1.;
      j++;
   }
   if(j>1){
      for(i=0;i<j-1;i++){
         d3[i]=(d1[i+1]-d1[i])/(k1[i+1]-k1[i]);
         for(k=i+2;k<j;k++){
             if((d1[k]-d1[i])/(k1[k]-k1[i])>d3[i])
                d3[i]=(d1[k]-d1[i])/(k1[k]-k1[i]);
         }
         d4[i+1]=(d1[0]-d1[i+1])/(k1[0]-k1[i+1]);
         for(k=1;k<=i;k++){
             if((d1[k]-d1[i+1])/(k1[k]-k1[i+1])<d4[i+1])
                d4[i+1]=(d1[k]-d1[i+1])/(k1[k]-k1[i+1]);
         }
      }
   }
   logls[0]=j;
   for(i=0;i<j;i++){
      logls[i*11+1]=k1[i];
      logls[i*11+3]=d1[i];
      logls[i*11+2]=d2[i];
      logls[i*11+4]=meas[1][k1[i]-1];
      logls[i*11+5]=meas[2][k1[i]-1];
      logls[i*11+6]=meas[3][k1[i]-1];
      logls[i*11+7]=meas[4][k1[i]-1];
      logls[i*11+8]=meas[5][k1[i]-1];
      logls[i*11+9]=k2[i];
      if(i==0){
         logls[i*11+10]=2.*d3[i];
         logls[i*11+11]=4.*d3[i];
      }
      if(i==j-1){
         logls[i*11+10]=0.;
         logls[i*11+11]=2.*d4[i];
      }
      if(i>0 && i<j-1){
         if(d4[i]>=d3[i]){
            logls[i*11+10]=2.*d3[i];
            logls[i*11+11]=2.*d4[i];
         }
         else{
            logls[i*11+10]= -1;
            logls[i*11+11]= -1;
         }
      }
   }
}
/******************************************************************************/
static int aicbest(ads,ranges,losses,logls)
int *ads;
double *ranges,*losses,*logls;
{
   int i,j,k,l,*k1,*k2;
   double *d1,*d3,*d4,*d2;
   k1=iv1;
   k2=iv2;
   d1=v1;
   d2=v2;
   d3=v3;
   d4=v6;
   j=0;
   for(i=0;i<maxdim;i++) if(ads[i]>=0){
      k1[j]=i+1;
      d1[j]=logls[i];
      d2[j]=losses[i];
      d3[j]= -2.;
      d4[j]= -1.;
      j++;
   }
   if(j>1){
      for(i=0;i<j-1;i++){
         d3[i]=(d1[i+1]-d1[i])/(k1[i+1]-k1[i]);
         for(k=i+2;k<j;k++){
             if((d1[k]-d1[i])/(k1[k]-k1[i])>d3[i])
                d3[i]=(d1[k]-d1[i])/(k1[k]-k1[i]);
         }
         d4[i+1]=(d1[0]-d1[i+1])/(k1[0]-k1[i+1]);
         for(k=1;k<=i;k++){
             if((d1[k]-d1[i+1])/(k1[k]-k1[i+1])<d4[i+1])
                d4[i+1]=(d1[k]-d1[i+1])/(k1[k]-k1[i+1]);
         }
      }
   }
   for(i=1;i<j-1;i++)if(d4[i]<d3[i]){
      d4[i]= -3.;
      d3[i]= -3.;
   }
   l=0;
   for(i=0;i<j;i++) if(i==0 || i==(j-1) || d4[i]>0){
      if(i==j-1)d3[l]=0.;
      else d3[l]=2.*d3[i];
      k2[l]=k1[i];
      d2[l]=d2[i];
      l++;
   }
   for(i=0;i<l;i++){
      ranges[l-1-i]=d3[i];
      losses[l-1-i]=d2[i];
   }
   return l;
}
/******************************************************************************/
/* this function re-initializes storage for a space */

static void predefinespace(data,spc)
struct datastruct *data;
struct space *spc;

{
   int ncov=(*data).ncov,nclass=(*data).nclass,i,j,k,l;

/* ncov,ndata,nclass - save typing
   i,j,k,l     - counters
   newspace    - thing to be initialized */

/* defines the basisfunctions, initializes */
   for(i=0;i<maxdim;i++) {
      (*spc).basis[i].b1= -1;
      (*spc).basis[i].b2= -1;
      (*spc).basis[i].t1= -1;
      (*spc).basis[i].t2= -1;
      (*spc).basis[i].ib= nclass;
      for(j=0;j<nclass;j++)(*spc).basis[i].beta[j]=0.;
      for(j=0;j<nclass;j++)(*spc).basis[i].link1[j]=j;
      for(j=0;j<=nclass;j++)(*spc).basis[i].link2[j]=0;
   }

/* for the 2-covariate subdimensions */
   for(i=0;i<ncov;i++) for(j=i+1;j<ncov;j++){
      for(k=0;k<MAXKNOTS+1;k++) for(l=0;l<MAXKNOTS+1;l++)
          (*spc).sub[i][j].kts1[k][l]=0;
      (*spc).sub[i][j].dim1=0;
   }

/* for the 1-covariate subdimensions */
   for(j=0;j<ncov;j++){
      for(k=0;k<MAXKNOTS;k++) (*spc).sub[j][ncov].ktsc[k]=0.;
      (*spc).sub[j][ncov].dim1=0;
   }
}

/******************************************************************************/
static void proj(spc,data,bbb,anova)
struct space *spc;
struct datastruct *data;
double *bbb,*anova;
{
   double x0,x1,x2,x3,x4,vaa,vbb;
   int i,j,k,k0,k1,k2,k3,k4,k5,k6,ncov=(*data).ncov,nbas=(*spc).nbas;
   int nclass=(*data).nclass;
   for(i=0;i<nbas;i++) for(j=0;j<nbas;j++)w1[i][j]=0.;
   w1[0][0]=1.;
   for(j=0;j<nbas;j++){
      iw1[j][0]=(*spc).basis[j].b1;
      iw1[j][1]=(*spc).basis[j].b2;
      iw1[j][2]=(*spc).basis[j].t1;
      iw1[j][3]=(*spc).basis[j].t2;
   }
   for(j=0;j<nbas;j++){
      (*spc).basis[j].j1=0;
      (*spc).basis[j].j2=0;
      if(iw1[j][0]!=ncov && iw1[j][1]==ncov && iw1[j][2]<0 ){
         x0=0;
         for(k=0;k<(*data).ndata;k++){
            vaa=petvector2(spc,data,j,k);
            x0+=vaa*(*data).wgt[k];
         }
         x0=x0/(*data).wgtsum;
         w1[j][j]=1.;
         w1[j][0]= -x0;
      }
   }
   for(j=0;j<nbas;j++){
      if(iw1[j][0]!=ncov && iw1[j][1]==ncov && iw1[j][2]>=0 ){
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==ncov && iw1[i][2]== -1){
               k=i;
               i=nbas;
            }
         }
         x0=0.;x1=0.;x2=0.;x3=0.;
         i=k;
         for(k=0;k<(*data).ndata;k++){
            vaa=petvector2(spc,data,i,k);
            vbb=petvector2(spc,data,j,k);
            x4=(*data).wgt[k]*vaa;
            x0+=x4;
            x1+=(*data).wgt[k]*vbb;
            x2+=x4*vaa;
            x3+=x4*vbb;
         }
         x3=(*data).wgtsum*x3-x0*x1;
         x3=x3/((*data).wgtsum*x2-x0*x0);
         x2=(x1-x3*x0)/(*data).wgtsum;
         w1[j][j]=1.;
         w1[j][i]= -x3;
         w1[j][0]= -x2;
      }
   }
   for(j=0;j<nbas;j++){
      if(iw1[j][0]!=ncov && iw1[j][1]!=ncov && iw1[j][2]<0 && iw1[j][3]<0){
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==ncov && iw1[i][2]== -1){
               k1=i;
               i=nbas;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][1]&&iw1[i][1]==ncov && iw1[i][2]== -1){
               k2=i;
               i=nbas;
            }
         }
         w1[j][j]=1.;
         w1[j][k1]=w1[k2][0];
         w1[j][k2]=w1[k1][0];
         w1[j][0]=w1[k2][0]*w1[k1][0];
         (*spc).basis[j].j1=k1;
         (*spc).basis[j].j2=k2;
      }
   }
   for(j=0;j<nbas;j++){
      if(iw1[j][0]!=ncov && iw1[j][1]!=ncov && iw1[j][2]<0 && iw1[j][3]>=0){
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==ncov && iw1[i][2]== -1){
               k1=i;
               i=nbas;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][1]&&iw1[i][1]==ncov){
               if(iw1[i][2]== -1)k2=i;
               if(iw1[i][2]==iw1[j][3])k3=i;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==iw1[j][1] && iw1[i][2]== -1 
               && iw1[i][3]== -1){
               k4=i;
               i=nbas;
            }
         }
         w1[j][ 0]=w1[k3][ 0]*w1[k1][ 0];
         w1[j][k2]=w1[k3][k2]*w1[k1][ 0];
         w1[j][k3]=w1[k3][k3]*w1[k1][ 0];
         w1[j][k1]=w1[k3][ 0]*w1[k1][k1];
         w1[j][k4]=w1[k3][k2]*w1[k1][k1];
         w1[j][ j]=w1[k3][k3]*w1[k1][k1];
         (*spc).basis[j].j1=k1;
         (*spc).basis[j].j2=k3;
      }
   }
   for(j=0;j<nbas;j++){
      if(iw1[j][0]!=ncov && iw1[j][1]!=ncov && iw1[j][2]>=0 && iw1[j][3]<0){
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][1]&&iw1[i][1]==ncov && iw1[i][2]== -1){
               k1=i;
               i=nbas;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==ncov){
               if(iw1[i][2]== -1)k2=i;
               if(iw1[i][2]==iw1[j][2])k3=i;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==iw1[j][1] && iw1[i][2]== -1 
               && iw1[i][3]== -1){
               k4=i;
               i=nbas;
            }
         }
         w1[j][ 0]=w1[k3][ 0]*w1[k1][ 0];
         w1[j][k2]=w1[k3][k2]*w1[k1][ 0];
         w1[j][k3]=w1[k3][k3]*w1[k1][ 0];
         w1[j][k1]=w1[k3][ 0]*w1[k1][k1];
         w1[j][k4]=w1[k3][k2]*w1[k1][k1];
         w1[j][ j]=w1[k3][k3]*w1[k1][k1];
         (*spc).basis[j].j1=k3;
         (*spc).basis[j].j2=k1;
      }
   }
   for(j=0;j<nbas;j++){
      if(iw1[j][0]!=ncov && iw1[j][1]!=ncov && iw1[j][2]>=0 && iw1[j][3]>=0){
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==ncov){
               if(iw1[i][2]== -1)k0=i;
               if(iw1[i][2]==iw1[j][2])k1=i;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][1]&&iw1[i][1]==ncov){
               if(iw1[i][2]== -1)k2=i;
               if(iw1[i][2]==iw1[j][3])k3=i;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][0]&&iw1[i][1]==iw1[j][1] && iw1[i][2]== -1 
               && iw1[i][3]== -1){
               k4=i;
               i=nbas;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][1]&&iw1[i][1]==iw1[j][1] && iw1[i][3]== -1 
               && iw1[i][2]==iw1[j][2] ){
               k5=i;
               i=nbas;
            }
         }
         for(i=1;i<nbas;i++){
            if(iw1[i][0]==iw1[j][1]&&iw1[i][1]==iw1[j][1] && iw1[i][2]== -1 
               && iw1[i][3]==iw1[j][3] ){
               k6=i;
               i=nbas;
            }
         }
         w1[j][ 0]=w1[k1][ 0]*w1[k3][ 0];
         w1[j][k0]=w1[k1][k0]*w1[k3][ 0];
         w1[j][k1]=w1[k1][k1]*w1[k3][ 0];
         w1[j][k2]=w1[k1][ 0]*w1[k3][k2];
         w1[j][k3]=w1[k1][ 0]*w1[k3][k3];
         w1[j][k4]=w1[k1][k0]*w1[k3][k2];
         w1[j][k5]=w1[k1][k1]*w1[k3][k2];
         w1[j][k6]=w1[k1][k0]*w1[k3][k3];
         w1[j][ j]=w1[k1][k1]*w1[k3][k3];
         (*spc).basis[j].j1=k1;
         (*spc).basis[j].j2=k3;
      }
   }
   for(i=0;i<nbas;i++) for(j=0;j<nbas;j++)w2[i][j]=0.;
   for(i=0;i<(*data).ndata;i++){
      petvector(spc,data,v4,v5,i);
      for(j=0;j<nbas;j++)v5[j]=0.;
      for(j=0;j<nbas;j++)if((*spc).basis[j].j1==0) for(k=0;k<nbas;k++) 
         v5[j]+=v4[k]*w1[j][k];
      for(j=0;j<nbas;j++){
         x0=v5[j]*(*data).wgt[i];
         if((*spc).basis[j].j1==0)for(k=0;k<nbas;k++)w2[k][j]+=v5[k]*x0;
      }
   }
   for(j=0;j<nbas;j++) for(k=0;k<nbas;k++)w2[k][j]=w2[k][j]/(*data).wgtsum;
   (void)lusolinv(w1,nbas,v1,0);
   for(i=0;i<nbas*nclass+1;i++)bbb[i]=0.;
   for(i=0;i<nbas;i++){
      for(j=0;j<nclass;j++){
         k0=j+i*(nclass+1);
         for(k1=0;k1<nbas;k1++){
            bbb[k0]+=w1[k1][i]*(*spc).basis[k1].beta[j];
         }
      }
   }
   for(i=0;i<nbas;i++){
      x0=0;
      for(j=0;j<nclass;j++) x0+=bbb[j+i*(nclass+1)];
      x0=x0/(double)(nclass+1);
      for(j=0;j<nclass+1;j++) bbb[j+i*(nclass+1)]-=x0;
   }
   for(j=0;j<nbas*3+1;j++)anova[j]=0.;
   anova[0]=3;
   anova[1]= -1;
   anova[2]= -1;
   for(j=0;j<nclass+1;j++)anova[3]+=bbb[j]*bbb[j];
   for(i=0;i<ncov;i++){
      k0= -1;
      k1=anova[0];
      anova[k1+1]=i+1;
      anova[k1+2]= -1;
      for(j=1;j<nbas;j++)if(iw1[j][0]==i && iw1[j][1]== ncov){
         for(k=1;k<nbas;k++)if(iw1[k][0]==i && iw1[k][1]== ncov){
            if(k0== -1){ 
               anova[0]+=3.;
               k0=0;
            }
            for(k2=0;k2<nclass+1;k2++)
               anova[k1+3]+=w2[j][k]*bbb[k2+j*(nclass+1)]*bbb[k2+k*(nclass+1)];
         }
      }
   }
   for(i=0;i<ncov;i++)for(k3=0;k3<ncov;k3++){
      k0= -1;
      k1=anova[0];
      anova[k1+1]=i+1;
      anova[k1+2]=k3+1;
      for(j=1;j<nbas;j++)if(iw1[j][0]==i && iw1[j][1]== k3){
         for(k=1;k<nbas;k++)if(iw1[k][0]==i && iw1[k][1]== k3){
            if(k0== -1){ 
               anova[0]+=3.;
               k0=0;
            }
            for(k2=0;k2<nclass+1;k2++)
               anova[k1+3]+=bbb[k2+j*(nclass+1)]*bbb[k2+k*(nclass+1)]*
                  w2[(*spc).basis[j].j1][(*spc).basis[k].j1]*
                  w2[(*spc).basis[j].j2][(*spc).basis[k].j2];
         }
      }
   }
   x0=0;
   i=anova[0]/3.;
   for(j=1;j<=i;j++)x0+=anova[j*3];
   for(j=1;j<=i;j++)anova[j*3]=anova[j*3]/x0;
   for(i=0;i<nbas;i++){
      if((*spc).basis[i].j1==0)x0=w2[i][i];
      else x0=w2[(*spc).basis[i].j1][(*spc).basis[i].j1]
             *w2[(*spc).basis[i].j2][(*spc).basis[i].j2];
      x0=sqrt(x0);
      for(j=0;j<nclass+1;j++)bbb[j+i*(nclass+1)]=bbb[j+i*(nclass+1)]*x0;
   }
}

/******************************************************************************/
/* this file contains the Newton-Raphson iteration loop subroutine, and the
   functions which it needs */

/* this routine does the Newton-Raphson iteration loop */

static double pnewton(spc,data)
struct datastruct *data; 
struct space *spc;

/* spc       - the present model
   data      - the data  */

{
   double zerror=0.0001,lnew=1,logl=0;
   int i,iter,ihalf,maxiter=300,j,k,l,mm;

/* zerror  - convergence criterion
   lnew   - new loglikelihood
   pcompall- routine to compute hessian, score and likelihood
   logl   - log-likelihood 
   iter   - iteration counter
   ihalf  - how many times has there been step-halving
   maxiter- maximum number of iterations 
   i,j,k,l  - counter 
   getinfox - gets the epsilon term of the hessian
   dlink - linkercombination */

/* initialize link2 */
   l= -1;
   for(j=0;j<(*spc).nbas;j++) for(k=0;k<(*spc).basis[j].ib;k++){
      l++;
      (*spc).basis[j].link2[k]=l;
   }
   if(l!=(*spc).ndim-1)error("serious linky-dim problem");

/* allocation/initialization/start iteration */
   getinfox(spc,data);
   mm=0;

   for(iter=0;iter<maxiter;iter++){
/* compute score, info and old loglikelihood */
      if(iter==0 || mm==3 || iter >7){
         logl=pcompall(spc,data,1);
         if(logl<(*data).ndata*log((double)(1./((*data).nclass+1.)))){
            for(j=0;j<(*spc).nbas;j++) for(k=0;k<(*data).nclass;k++) 
               (*spc).basis[j].beta[k]=0.;
            logl=pcompall(spc,data,1);
         }
         i=lusolinv((*spc).info,(*spc).ndim,v2,1);
         mm=2;
      }   
      else{
         if(mm==7){
            iter=1;
            mm=1;
            for(j=0;j<(*spc).nbas;j++) for(k=0;k<(*data).nclass;k++) 
               (*spc).basis[j].beta[k]=0.;
            logl=pcompall(spc,data,1);
         }
         else{
            for(j=0;j<(*spc).nbas;j++) for(k=0;k<(*data).nclass;k++){
               i=dlink(spc,j,k);
               if(i>=0) v5[i]=(*spc).basis[j].beta[k];
            }
            logl=pcomp2(spc,data);
            mm=1;
         }
      }

/* Rprintf("%d %d %f\n",iter,mm,logl);fflush(stdout); */
        
/* solve system */
      for(j=0;j<(*spc).ndim;j++)v3[j]=(*spc).score[j];
      for(j=0;j<(*spc).nbas;j++){
         for(k=0;k<(*data).nclass;k++){
            i=dlink(spc,j,k);
            if(i>=0) v6[i]=(*spc).basis[j].beta[k];
         }
      }
      for(i=0;i<(*spc).ndim;i++){
         (*spc).score[i]=0.;
         for(j=0;j<(*spc).ndim;j++)(*spc).score[i]+=v3[j]*(*spc).info[i][j];
      }
      ihalf=1;
      do{

/* compute new loglikelihood (cheat on beta) */
         for(j=0;j<(*spc).nbas;j++){
            for(k=0;k<(*data).nclass;k++){
               i=dlink(spc,j,k);
               if(i>=0) (*spc).basis[j].beta[k] -= (*spc).score[i];
            }
         }
         lnew=pcompall(spc,data,0);
/* Rprintf("%f ",lnew);fflush(stdout); */
/* get beta back again */
         for(j=0;j<(*spc).nbas;j++){
            for(k=0;k<(*data).nclass;k++){
               i=dlink(spc,j,k);
               if(i>=0) (*spc).basis[j].beta[k] += (*spc).score[i];
            }
         }

/* step halving if required */
         if(lnew<logl-zerror){ 

/* oops, too much step halving */
            if(ihalf>8192 && mm==2 && iter>0)return 200.;
            if(ihalf>8192 && mm==2)mm=7;
            if(ihalf>8192 && mm!=7){
               mm=3;
               (void)Rprintf("step half ouch...\n");(void)fflush(stdout);
            }

/* the actual halving */
            ihalf=ihalf*2;
            for(j=0;j<(*spc).ndim;j++) (*spc).score[j]=(*spc).score[j]/2;
         }
      }while(lnew<logl-zerror && mm!=3 && mm!=7);

/* record the new solution: add the step */
      if(mm!=3 && mm!=7)for(j=0;j<(*spc).nbas;j++){
         for(k=0;k<(*data).nclass;k++){
            i=dlink(spc,j,k);
            if(i>=0) (*spc).basis[j].beta[k] -= (*spc).score[i];
         }
      }

/* did we converge */
      if(plumbertester(logl)+plumbertester(lnew)!=6 &&mm!=3 &&mm!=7)
         return 200.;
      if(lnew-logl<zerror &&mm!=3 &&mm!=7) iter=maxiter+1000; 
   }
/* did we finish because we converged? */
   if(iter<maxiter+500) return 200.;

/* get a fresh copy of info */
   for(j=0;j<(*spc).nbas;j++) for(k=0;k<(*data).nclass;k++){
      i=dlink(spc,j,k);
      if(i>=0) v5[i]=(*spc).basis[j].beta[k];
   }
   /* lnew=pcomp2(spc,data); */
   lnew=pcompall(spc,data,2); 
   for(i=0;i<(*spc).ndim;i++)for(j=0;j<(*spc).ndim;j++)
      w3[i][j]=(*spc).info[i][j];
   i=lusolinv(w3,(*spc).ndim,v2,1);
   return lnew;
}

/******************************************************************************/
/* this routine computes hessian score and log-likelihood */

static double pcompall(spc,data,what)
struct space *spc;
struct datastruct *data;
int what;

/* spc   - the present model
   data  - the data 
   what  - loglikelihood (0) or also score and hessian? */

{
   int i,j,k,j2,k2,i2,i1,itemp,i3,j4;
   int nclass=(*data).nclass,nbas=(*spc).nbas,j3,ncl;
   double logl,rtemp,**xinfo,*xscore,**dwk,**dwk2,epsilon=(*spc).epsilon,xx;
   double d1,*dwl,*dwl2;
/* i      - typically counter data
   j      - typically counter basisfunctions
   k      - typically counter classes
   j2     - typically counter basisfunctions
   k2     - typically counter classes
   logl   - loglikelihood 
   itemp,rtemp,i2 - save typing
   nclass,nbas,dwk - save typing
   xinfo and xscore - info and score were all basis functions independent
   dlink  - link combination */

/* the computations are first all carried out as if all elements are independent
   lumping is done afterwards */

/* initializations  and allocations */
   ncl=nclass+1;
   dwk=(*data).work;
   dwk2=(*data).work2;
   xinfo=w1;
   xscore=v1;
   logl=0.;
   if(what!=0){
      for(i=0;i<nclass*nbas;i++){
         xscore[i]=0.;
         (*spc).score[i]=0.;
         if(what<=2)for(j=0;j<nclass*nbas;j++){
             xinfo[i][j]= -(*spc).infox[i][j];
             (*spc).info[i][j]=0.;
         }
      }
   }

/* now circle the datapoints */
   for(i=0;i<(*data).ndata;i++){
      dwl=dwk[i];
      dwl2=dwk2[i];
      petvector(spc,data,v7,v8,i);
      xx=epsilon*(*data).wgt[i];
/* compute theta */
      dwl[nclass]=0.;
      dwl2[nclass]=0.;
      dwl2[ncl]=0.;
      for(k=0;k<nclass;k++)dwl[k]= (*spc).basis[0].beta[k];
      for(j=1;j<nbas;j++){
         if(v7[j]!=0.){
            for(k=0;k<nclass;k++) { 
               dwl[k]+= v7[j]*(*spc).basis[j].beta[k];
            }
         }
      }
      for(k=0;k<nclass;k++) { 
         dwl2[k]=dwl[k];
         dwl2[ncl]+=dwl[k];
         if(dwl[k]<600)dwl[k]=exp(dwl[k]);
         else dwl[k]=exp((double)600);
      }
      for(k=0;k<nclass;k++){
         dwl[nclass]+=dwl[k];
         dwl2[k]=dwl2[k]-dwl2[ncl]/((double)(ncl));
      }
/* compute -c */
      dwl2[nclass]=dwl2[nclass]-dwl2[ncl]/((double)(ncl));
      dwl[nclass]=1./(1.+dwl[nclass]);
/* make them theta - c */
      for(k=0;k<nclass;k++)  dwl[k]=dwl[k]*dwl[nclass]; 
/* update the log-likelihhod */
      logl+=pylog(dwl[(*data).yy[i]])*((*data).wgt[i]);  
/* epsilon part */
      if(what!=2)for(k=0;k<=nclass;k++) logl-= xx*dwl2[k]*dwl2[k];
/* exponentiate theta-c */
      if(what!=0){
         xx=xx*2./(double)(ncl);
/* preparation */
         for(k=0;k<nclass;k++){
            dwl2[k]=xx*dwl2[k];
            w2[0][k]=dwl[k];
            w3[0][k]=dwl2[k];
         }
         dwl2[nclass]=xx*dwl2[nclass];
         w3[0][nclass]=dwl2[nclass];
         for(j=1;j<nbas;j++){
            if(v7[j]!=0.){
               for(k=0;k<nclass;k++){
                  w3[j][k]=v7[j]*dwl2[k];
                  w2[j][k]=v7[j]*dwl[k];
               }
               w3[j][nclass]=v7[j]*dwl2[nclass];
            }
         }
/* score */
         for(j=0;j<nbas;j++){
            if(v7[j]!=0.){
               d1=0;
               j3=j*nclass;
               for(k=0;k<=nclass;k++) d1 += w3[j][k];
               for(k=0;k<nclass;k++){
                  rtemp=v8[j]*dwl[k];
                  itemp=k+j3;
                  xscore[itemp]+= d1-dwl2[k]*v7[j]*(ncl)- rtemp; 
/* hessian */
                  for(j2=j;j2<nbas;j2++) { 
                     if(v7[j2]!=0.){
                        j4=j2*nclass;
                        xinfo[itemp][k+j4]-= rtemp*v7[j2];
                        for(k2=0;k2<nclass;k2++){
                           xinfo[itemp][k2+j4]+=rtemp*w2[j2][k2];
                        }  
                     }
                  }
               }
/* delta parts */
               if(((*data).yy[i])!=nclass) xscore[(*data).yy[i]+j3]+= v8[j];
            }
         }
      }
   }

   if(what!=0){
/* symmaterize the hessian */
      for(k=0;k<nbas*nclass-1;k++) for(j=k+1;j<nbas*nclass;j++){
         xinfo[j][k]=xinfo[k][j];
      }
/* rearrange score and info according to linky */
      for(j=0;j<nbas;j++) for(k=0;k<nclass;k++){
         i1=dlink(spc,j,k);
         if(i1>=0){
            i2=j*nclass+k;
            (*spc).score[i1] += xscore[i2];
            if(what<=2)for(j2=0;j2<nbas;j2++) for(k2=0;k2<nclass;k2++){
               i3=dlink(spc,j2,k2);
               if(i3>=0) (*spc).info[i3][i1]+=xinfo[j2*nclass+k2][i2];
            }
         }
      }
   }
   (*spc).logl=logl;
   return logl;
}

/******************************************************************************/
/* get one link element */
static int dlink(spc,j,k)
struct space *spc;
int j,k;
{
   if((*spc).basis[j].link1[k]< 0)return -1;
   return (*spc).basis[j].link2[(*spc).basis[j].link1[k]];
}
/******************************************************************************/
/* get the info-x term */
static void getinfox(spc,data)
struct space *spc;
struct datastruct *data;
{
   int ndata=(*data).ndata,nclass=(*data).nclass,i,j1,k1,j2,k2,l1,j3;
   int nbas=(*spc).nbas,kk=nclass+1;
   double x,y,epsilon=2.*(*spc).epsilon,xkk;
   epsilon=epsilon/(double)(nclass+1);
   for(j1=0;j1<nbas;j1++) for(k1=0;k1<nclass;k1++){
      l1=k1+j1*nclass;
      for(j2=j1;j2<nbas;j2++) for(k2=0;k2<nclass;k2++){
         (*spc).infox[l1][k2+j2*nclass]=0.;
      }
   }
   for(i=0;i<ndata;i++){
      petvector(spc,data,v7,v8,i);
      y=(*data).wgt[i]*epsilon;
      for(j1=0;j1<nbas;j1++)for(j2=j1;j2<nbas;j2++){
         x=y*v7[j1]*v7[j2];
         xkk = x*kk;
         j3=j1*nclass;
         for(k1=0;k1<nclass;k1++){
            l1=k1+j3;
            for(k2=0;k2<nclass;k2++) (*spc).infox[l1][k2+j2*nclass] -= x;
            (*spc).infox[l1][k1+j2*nclass] += xkk;
         }
      }
   }
   for(k1=0;k1<nbas*nclass-1;k1++) for(j1=k1+1;j1<nbas*nclass;j1++)
      (*spc).infox[j1][k1]=(*spc).infox[k1][j1];
}
/*************/
static double pylog(x)
double x;
{
if(x < 10.e-250)return (double)(-576.);
else return log(x);
}
/*************/
/* this routine computes hessian score and log-likelihood  quasi newton*/

static double pcomp2(spc,data)
struct space *spc;
struct datastruct *data;

/* spc   - the present model
   data  - the data */

{
   int i,j,k,i2,i1,itemp;
   int nclass=(*data).nclass,nbas=(*spc).nbas,j3,ncl;
   double logl,rtemp,*xscore,**dwk,**dwk2,epsilon=(*spc).epsilon,xx;
   double d1,*dwl,*dwl2,d0;
/* i      - typically counter data
   j      - typically counter basisfunctions
   k      - typically counter classes
   logl   - loglikelihood 
   itemp,rtemp,i2 - save typing
   nclass,nbas,dwk - save typing
   xinfo and xscore - info and score were all basis functions independent
   dlink  - link combination */

/* the computations are first all carried out as if all elements are independent
   lumping is done afterwards */

/* initializations  and allocations */
   ncl=nclass+1;
   dwk=(*data).work;
   dwk2=(*data).work2;
   xscore=v1;
   logl=0.;
   for(i=0;i<nclass*nbas;i++){
      xscore[i]=0.;
      (*spc).score[i]=0.;
   }

/* now circle the datapoints */
   for(i=0;i<(*data).ndata;i++){
      dwl=dwk[i];
      dwl2=dwk2[i];
      petvector(spc,data,v7,v8,i);
      xx=epsilon*(*data).wgt[i];
/* compute theta */
      dwl[nclass]=0.;
      dwl2[nclass]=0.;
      dwl2[ncl]=0.;
      for(k=0;k<nclass;k++)dwl[k]= (*spc).basis[0].beta[k];
      for(j=1;j<nbas;j++){
         if(v7[j]!=0.){
            for(k=0;k<nclass;k++) { 
               dwl[k]+= v7[j]*(*spc).basis[j].beta[k];
            }
         }
      }
      for(k=0;k<nclass;k++) { 
         dwl2[k]=dwl[k];
         dwl2[ncl]+=dwl[k];
         if(dwl[k]<600)dwl[k]=exp(dwl[k]);
         else dwl[k]=exp((double)600);
      }
      for(k=0;k<nclass;k++){
         dwl[nclass]+=dwl[k];
         dwl2[k]=dwl2[k]-dwl2[ncl]/((double)(ncl));
      }
/* compute -c */
      dwl2[nclass]=dwl2[nclass]-dwl2[ncl]/((double)(ncl));
      dwl[nclass]=1./(1.+dwl[nclass]);
/* make them theta - c */
      for(k=0;k<nclass;k++)  dwl[k]=dwl[k]*dwl[nclass]; 
/* update the log-likelihhod */
      logl+=pylog(dwl[(*data).yy[i]])*((*data).wgt[i]);  
/* epsilon part */
      for(k=0;k<=nclass;k++) logl-= xx*dwl2[k]*dwl2[k];
/* exponentiate theta-c */
      xx=xx*2./(double)(ncl);
/* preparation */
      for(k=0;k<nclass;k++){
         dwl2[k]=xx*dwl2[k];
         w2[0][k]=dwl[k];
         w3[0][k]=dwl2[k];
      }
      dwl2[nclass]=xx*dwl2[nclass];
      w3[0][nclass]=dwl2[nclass];
      for(j=1;j<nbas;j++){
         if(v7[j]!=0.){
            for(k=0;k<nclass;k++){
               w3[j][k]=v7[j]*dwl2[k];
               w2[j][k]=v7[j]*dwl[k];
            }
            w3[j][nclass]=v7[j]*dwl2[nclass];
         }
      }
/* score */
      for(j=0;j<nbas;j++){
         if(v7[j]!=0.){
            d1=0;
            j3=j*nclass;
            for(k=0;k<=nclass;k++) d1 += w3[j][k];
            for(k=0;k<nclass;k++){
               rtemp=v8[j]*dwl[k];
               itemp=k+j3;
               xscore[itemp]+= d1-dwl2[k]*v7[j]*(ncl)- rtemp; 
            }
            if(((*data).yy[i])!=nclass) xscore[(*data).yy[i]+j3]+= v8[j];
         }
      }
   }

/* rearrange score and info according to linky */
   for(j=0;j<nbas;j++) for(k=0;k<nclass;k++){
      i1=dlink(spc,j,k);
      if(i1>=0){
         i2=j*nclass+k;
         (*spc).score[i1] += xscore[i2];
      }
   }
   for(j=0;j<(*spc).ndim;j++){
      v5[j]-=v6[j];
      v6[j]=(*spc).score[j]-v3[j];
   }
   for(j=0;j<(*spc).ndim;j++){
      v3[j]=0.;
      for(i=0;i<(*spc).ndim;i++) v3[j]+=v6[i]*(*spc).info[i][j];
   }
   d0=0.;
   d1=0.;
   for(j=0;j<(*spc).ndim;j++){
      d0+=v6[j]*v5[j];
      d1+=v3[j]*v6[j];
   }
   for(j=0;j<(*spc).ndim;j++) for(i=0;i<(*spc).ndim;i++)
      (*spc).info[i][j]+=v5[i]*v5[j]/d0-v3[i]*v3[j]/d1;
   (*spc).logl=logl;
   return logl;
}

/******************************************************************************/

/* this routine pearches all dimensions for a basis function to remove */

static int prembas(spc,data,silent)
struct space *spc;
struct datastruct *data;
int silent;

/* spc   - the model from which to remove something
   data  - data
   silent- should diagnostic output be printed?  */

{
   int nclass=(*data).nclass,ncov=(*data).ncov;
   int nb1,nt1,nb2,nt2,j,k,j2,eligible,i,bbi,l,i1;
   int nbas=(*spc).nbas,bb1,bb2,bt1,bt2,**tlink;
   double criterion,wald,**tinfo;
/* nclass,ncov,nbas- save typing
   nb1,nt1,nb2,nt2 - b and t attributes of the basis function being studied
   aj,ak           - a beta which corresponds to them
   sj,sk           - their score components
   j,k             - the independent components that we try to equate
   j2              - possible conflicting basis function
   bbi,bbj,bbk     - best i, j and k
   n               - if 0, a potential conflict
   eligilble       - conflicting basisfunction
   i               - the basis function that we study
   l               - for initialization of links 
   k1,k2,k3,k4     - used for old-new link relations
   tlink           - present link relations
   tinfo           - used when sweeping info
   baj,bak,bsj,bsk - best aj,ak,sj,sk
   ii,i1,j1        - counters 
   bb1,bb2,bt1     - b's and t's for removal 
   wald            - wald statistic
   criterion       - best wald up to now
   baj,bak         - best aj and ak 
   dlink           - find a double link */

   tlink=iw1;
   for(i=0;i<nbas;i++)for(j=0;j<nclass;j++)tlink[i][j]=dlink(spc,i,j);

/* initialization */
   criterion=pow((double)10.,(double)100.);
   bbi = -1;
   tinfo=w1;
   for(i=0;i<(*spc).ndim;i++){
      for(j=0;j<(*spc).ndim;j++){
         tinfo[i][j]=(*spc).info[i][j];
         (*spc).info[i][j]=w3[i][j];
      }
   }

/* circle all the basisfunctions except for the first one */
   for(i=0;i<nbas;i++){
/* what are the b and t attributes ? */
      nb1=(*spc).basis[i].b1;
      nt1=(*spc).basis[i].t1;
      nb2=(*spc).basis[i].b2;
      nt2=(*spc).basis[i].t2;

/* we propose to equate classes j and k - find an involved beta (aj, ak) and
   the involved scores/info (ak, sk) first for j */
      eligible=1;
/* see whether this is an allowable removal  - run the other basisfunctions */
      for(j2=1;j2<nbas;j2++) if(j2!=i){
/* constant - all possible confilcts */
         if(i==0)eligible=0;
         if(i!=0){
            if(nb2==ncov){
/* potential conflicts of a 1 d selected knot with a 2 d existing  knot */
               if((*spc).basis[j2].b1==nb1&&(*spc).basis[j2].t1==nt1)eligible=0;
               if((*spc).basis[j2].b2==nb1&&(*spc).basis[j2].t2==nt1)eligible=0;
/* potential conflicts of a 1 d selected linear with a 1/2d existing anything */
               if(nt1== -1 && (*spc).basis[j2].b1==nb1)eligible=0;
               if(nt1== -1 && (*spc).basis[j2].b2==nb1)eligible=0;
            }
            else{
/* potential conflicts between 2 d linears and 2 d knots */
               if((*spc).basis[j2].b1==nb1 && (*spc).basis[j2].b2==nb2) {
                  if(nt2== -1 && nt1== -1) eligible=0;
                  if(nt1== -1 && nt2>= 0 && (*spc).basis[j2].t2==nt2)eligible=0;
                  if(nt1>= 0 && nt2== -1 && (*spc).basis[j2].t1==nt1)eligible=0;
               }
            }
         }
/* are they already the same on a higher level? */
      }

/* if we are eligible, compute wald */
      if(eligible==1){
         for(j=0;j<(*data).nclass;j++){
            for(k=0;k<(*data).nclass;k++){
               w3[j][k]=(*spc).info[dlink(spc,i,j)][dlink(spc,i,k)];
            }
         }
         for(j=0;j<(*data).nclass;j++)v1[j]=(*spc).basis[i].beta[j];
         (void)lusolinv(w3,(*data).nclass,v1,2);
         wald=0.;  
         for(j=0;j<(*data).nclass;j++) wald+=(*spc).basis[i].beta[j]*v1[j];
         wald=fabs(wald);
/* did we improve ? */
         if(plumbertester(wald)!=2 && wald<criterion){
            bbi=i;
            criterion=wald;
         }
      }
   }
   if(bbi!= -1){
              
/* get the beta shift */
   j2=bbi*(*data).nclass;
   for(i=0;i<(*spc).nbas;i++){
      if(i!=bbi) for(j=0;j<(*data).nclass;j++){
         k=i*(*data).nclass+j;        
         v3[k]=0.;
         for(i1=0;i1<(*data).nclass;i1++){
            v3[k]+=(*spc).basis[bbi].beta[i1]*tinfo[j2+i1][k];
         }
      }
   }
   i1=((*spc).nbas-1)*(*data).nclass;
   for(i=0;i<(*spc).ndim;i++){
      for(j=0;j<(*data).nclass;j++){
         tinfo[i][j2+j]=tinfo[i][i1+j];
      }
   }
   for(i=0;i<(*spc).ndim;i++){
      for(j=0;j<(*data).nclass;j++){
         tinfo[j2+j][i]=tinfo[i1+j][i];
      }
   }
   for(j=0;j<(*data).nclass;j++){
      v3[j2+j]=v3[i1+j];
   }
/* reduce dim */
   (*spc).ndim-=(*data).nclass;
   (void)lusolinv(tinfo,(*spc).ndim,v3,2);

/* now remove the worst dimension */
   (*spc).nbas-=1;
   nbas-=1;
   bb1=(*spc).basis[bbi].b1;
   bb2=(*spc).basis[bbi].b2;
   bt1=(*spc).basis[bbi].t1;
   bt2=(*spc).basis[bbi].t2;

/* announce the results */
   if(silent!=1){
      puuu(spc,bb1,bb2,bt1,bt2,(*data).ncov,1);
      (void)Rprintf("(wald=%.2f)\n",criterion);
   }


/* move the last basisfunction to the one that is removed */
   (*spc).basis[bbi].b1=(*spc).basis[nbas].b1;
   (*spc).basis[bbi].b2=(*spc).basis[nbas].b2;
   (*spc).basis[bbi].t1=(*spc).basis[nbas].t1;
   (*spc).basis[bbi].t2=(*spc).basis[nbas].t2;
   (*spc).basis[bbi].ib=(*spc).basis[nbas].ib;
   for(i=0;i<nclass+1;i++){
      (*spc).basis[bbi].beta[i]=(*spc).basis[nbas].beta[i];
      (*spc).basis[bbi].link1[i]=(*spc).basis[nbas].link1[i];
      tlink[bbi][i]=tlink[nbas][i];
   }

/* it is a 1d basisfunction */
   if(bb2==ncov){
/* change the number of knots */
      (*spc).sub[bb1][ncov].dim1-=1;

/* change all the other pointers */
      if((*spc).sub[bb1][ncov].dim1>0){

/* first in the basisfunctions */
         for(j=0;j<nbas;j++){
            if((*spc).basis[j].b1==bb1)
               if((*spc).basis[j].t1>bt1)(*spc).basis[j].t1-=1;
            if((*spc).basis[j].b2==bb1)
               if((*spc).basis[j].t2>bt1)(*spc).basis[j].t2-=1;
         }
      }

/* in the knots themselves */
      for(j=bt1;j>(-1)&&j<(*spc).sub[bb1][ncov].dim1;j++){
         (*spc).sub[bb1][ncov].ktsc[j]=(*spc).sub[bb1][ncov].ktsc[j+1];
      }
   }
   else{
/* if it is a two variable dimension */
      (*spc).sub[bb1][bb2].dim1-=1;
   }
/* shift the beta */
   for(j=0;j<(*spc).nbas;j++){
      for(i=0;i<(*data).nclass;i++){
         (*spc).basis[j].beta[i]+=v3[j*(*data).nclass+i];
      }
   }

/* initialize link2 */
   l= -1;
   for(j=0;j<(*spc).nbas;j++){
      for(k=0;k<(*spc).basis[j].ib;k++){
         l++;
         (*spc).basis[j].link2[k]=l;
      }
   }
   }
   return bbi;
}
/******************************************************************************/

/* this routine pearches all dimensions for something to remove - it is a mess*/

static void premdim(spc,data,silent,dwald,iwald)
struct space *spc;
struct datastruct *data;
int silent,*iwald;
double *dwald;

/* spc   - the model from which to remove something
   data  - data
   silent- should diagnostic output be printed?  */

{
   int nclass=(*data).nclass,ncov=(*data).ncov;
   int nb1,nt1,nb2,nt2,aj,ak,sj,sk,j,k,j2,n,eligible,i,bbi,bbj,bbk,l;
   int nbas=(*spc).nbas,bb1,bb2,bt1,ii,baj,bak,bsj,bsk,**tlink;
   int k1,k2,k3,k4,i1,j1;
   double criterion,wald,xx,**tinfo;
/* nclass,ncov,ndata,nbas - save typing
   nb1,nt1,nb2,nt2 - b and t attributes of the basis function being studied
   aj,ak           - a beta which corresponds to them
   sj,sk           - their score components
   j,k             - the independent components that we try to equate
   j2              - possible conflicting basis function
   bbi,bbj,bbk     - best i, j and k
   n               - if 0, a potential conflict
   eligilble       - conflicting basisfunction
   i               - the basis function that we study
   l               - for initialization of links 
   k1,k2,k3,k4     - used for old-new link relations
   tlink           - present link relations
   tinfo           - used when sweeping info
   baj,bak,bsj,bsk - best aj,ak,sj,sk
   ii,i1,j1        - counters 
   bb1,bb2,bt1     - b's and t's for removal 
   wald            - wald statistic
   criterion       - best wald up to now
   baj,bak         - best aj and ak 
   dlink           - find a double link */

   tlink=iw1;
   for(i=0;i<nbas;i++)for(j=0;j<nclass;j++)tlink[i][j]=dlink(spc,i,j);

/* initialization */
   criterion=pow((double)10.,(double)100.);
   if((*iwald)==0)(void)lusolinv((*spc).info,(*spc).ndim,v1,1); 

/* circle all the basisfunctions except for the first one */
   for(i=0;i<nbas;i++){
/* what are the b and t attributes ? */
      nb1=(*spc).basis[i].b1;
      nt1=(*spc).basis[i].t1;
      nb2=(*spc).basis[i].b2;
      nt2=(*spc).basis[i].t2;

/* we propose to equate classes j and k - find an involved beta (aj, ak) and
   the involved scores/info (ak, sk) first for j */
      for(j=0;j<(*spc).basis[i].ib;j++){
         for(k=j;k<(*spc).basis[i].ib;k++){
            eligible=1;
            for(ii=0;ii<nclass;ii++){
               if((*spc).basis[i].link1[ii]==j){
                  aj=ii;
                  sj=(*spc).basis[i].link2[j];
                  ii=nclass+1;
               }
            }
/* then for k, if j==k we try to make this class equal to 0 */
            if(j!=k){
                  for(ii=0;ii<nclass;ii++){
                  if((*spc).basis[i].link1[ii]==k){
                     ak=ii;
                     sk=(*spc).basis[i].link2[k];
                     ii=nclass+1;
                  }
               }
            }
/* see whether this is an allowable removal  - run the other basisfunctions */
            for(j2=1;j2<nbas;j2++) if(j2!=i){
               n=1;
/* constant - all possible confilcts */
               if(i==0)n=0;
               if(i!=0){
                  if(nb2==ncov){
/* potential conflicts of a 1 d selected knot with a 2 d existing  knot */
                     if((*spc).basis[j2].b1==nb1&&(*spc).basis[j2].t1==nt1)n=0;
                     if((*spc).basis[j2].b2==nb1&&(*spc).basis[j2].t2==nt1)n=0;
/* potential conflicts of a 1 d selected linear with a 1/2d existing anything */
                     if(nt1== -1 && (*spc).basis[j2].b1==nb1)n=0;
                     if(nt1== -1 && (*spc).basis[j2].b2==nb1)n=0;
                  }
                  else{
/* potential conflicts between 2 d linears and 2 d knots */
                     if((*spc).basis[j2].b1==nb1 && (*spc).basis[j2].b2==nb2) {
                        if(nt2== -1 && nt1== -1) n=0;
                        if(nt1== -1 && nt2>= 0 && (*spc).basis[j2].t2==nt2)n=0;
                        if(nt1>= 0 && nt2== -1 && (*spc).basis[j2].t1==nt1)n=0;
                     }
                  }
               }
/* are they already the same on a higher level? */
               if(n==0){
                  if(j==k && (*spc).basis[j2].link1[aj]!= -1)eligible=0;
                  if(j!=k && (*spc).basis[j2].link1[aj]
                           !=(*spc).basis[j2].link1[ak])eligible=0;
               }
               if(eligible==0)j2=nbas;
            }

/* if we are eligible, compute wald */
            if(eligible==1){
/* easy if a beta has to be put equal to 0 */
               if(j==k){
                  wald=fabs((*spc).basis[i].beta[aj]/
                       sqrt(fabs((*spc).info[sj][sj])));
                  sk= -1;
                  ak= -1;
               }
               else{
/* otherwise select a 2x2 part of score and beta */
                  wald=fabs((*spc).basis[i].beta[aj]-(*spc).basis[i].beta[ak])/
                       sqrt(fabs(-(*spc).info[sj][sj]+2*(*spc).info[sj][sk]
                                                 -(*spc).info[sk][sk]));
               }
/* did we improve ? */
               if(wald<criterion){
                  bbi=i;
                  bbj=j;
                  bbk=k;
                  baj=aj;
                  bak=ak;
                  bsj=sj;
                  bsk=sk;
                  criterion=wald;
               } 
            }
         }
      }
   }
              
/* reduce dim */

   (*spc).ndim-=1;

/* sweep info */
   if(bsk!= -1){
      for(i=0;i<=(*spc).ndim;i++){
         (*spc).info[i][bsj]=(*spc).info[i][bsj]-(*spc).info[i][bsk];
      }
      for(i=0;i<=(*spc).ndim;i++){
         (*spc).info[bsj][i]=(*spc).info[bsj][i]-(*spc).info[bsk][i];
      }
   }
   for(i=0;i<=(*spc).ndim;i++) for(j=0;j<=(*spc).ndim;j++)
      w2[i][j]=(*spc).info[i][j];
   (void)lusolinv(w2,(*spc).ndim+1,v1,1);
   for(i=0;i<=(*spc).ndim;i++)v5[i]=w2[i][bsj];
   v5[bsj]=v5[(*spc).ndim];
   tinfo=w1;
   for(i=0;i<=(*spc).ndim;i++) if(i!=bsj){
      for(j=0;j<=(*spc).ndim;j++) if(j!=bsj)tinfo[i][j]=(*spc).info[i][j]-
            (*spc).info[i][bsj]*(*spc).info[j][bsj]/(*spc).info[bsj][bsj];
   }
   for(i=0;i<=(*spc).ndim;i++){
      tinfo[i][bsj]=0.;
      tinfo[bsj][i]=0.;
   }
/* get the beta shift */
   for(i=0;i<=(*spc).ndim;i++)
      for(j=0;j<=(*spc).ndim;j++)w2[i][j]=tinfo[i][j];
   for(i=0;i<=(*spc).ndim;i++) w2[i][bsj]=w2[i][(*spc).ndim];
   for(i=0;i<=(*spc).ndim;i++) w2[bsj][i]=w2[(*spc).ndim][i];
   for(i=0;i<(*spc).ndim;i++){
      v6[i]=0.;
      for(j=0;j<(*spc).ndim;j++)v6[i]-=w2[i][j]*v5[j];
   }
   if(bsk== -1)
      for(i=0;i<(*spc).ndim;i++)v6[i]=v6[i]*(*spc).basis[bbi].beta[baj];
   else for(i=0;i<(*spc).ndim;i++)
      v6[i]=v6[i]*((*spc).basis[bbi].beta[baj]-(*spc).basis[bbi].beta[bak]);
   for(i=0;i<nbas;i++)for(j=0;j<nclass;j++){
      k=dlink(spc,i,j);
      if(k>=0){
         if(k==(*spc).ndim)k=bsj;
         (*spc).basis[i].beta[j]-=v6[k];
      }
   }

/* announce the results */
   if(silent!=1){
      if((*spc).basis[bbi].ib>1) puuu(spc,(*spc).basis[bbi].b1,
           (*spc).basis[bbi].b2,(*spc).basis[bbi].t1,(*spc).basis[bbi].t2,
           (*data).ncov,2);
      else puuu(spc,(*spc).basis[bbi].b1,
           (*spc).basis[bbi].b2,(*spc).basis[bbi].t1,(*spc).basis[bbi].t2,
           (*data).ncov,1);
      (void)Rprintf("(wald=%.2f)\n",criterion*criterion);
   }

/* now remove the worst dimension; if ib==1 this is everything */
   if((*spc).basis[bbi].ib==1){
      (*spc).nbas-=1;
      nbas-=1;
      bb1=(*spc).basis[bbi].b1;
      bb2=(*spc).basis[bbi].b2;
      bt1=(*spc).basis[bbi].t1;

/* move the last basisfunction to the one that is removed */
      (*spc).basis[bbi].b1=(*spc).basis[nbas].b1;
      (*spc).basis[bbi].b2=(*spc).basis[nbas].b2;
      (*spc).basis[bbi].t1=(*spc).basis[nbas].t1;
      (*spc).basis[bbi].t2=(*spc).basis[nbas].t2;
      (*spc).basis[bbi].ib=(*spc).basis[nbas].ib;
      for(i=0;i<nclass+1;i++){
         (*spc).basis[bbi].beta[i]=(*spc).basis[nbas].beta[i];
         (*spc).basis[bbi].link1[i]=(*spc).basis[nbas].link1[i];
         tlink[bbi][i]=tlink[nbas][i];
      }

/* it is a 1d basisfunction */
      if(bb2==ncov){
/* change the number of knots */
         (*spc).sub[bb1][ncov].dim1-=1;

/* change all the other pointers */
         if((*spc).sub[bb1][ncov].dim1>0){

/* first in the basisfunctions */
            for(j=0;j<nbas;j++){
               if((*spc).basis[j].b1==bb1)
                  if((*spc).basis[j].t1>bt1)(*spc).basis[j].t1-=1;
               if((*spc).basis[j].b2==bb1)
                  if((*spc).basis[j].t2>bt1)(*spc).basis[j].t2-=1;
            }
         }

/* in the knots themselves */
         for(j=bt1;j>(-1)&&j<(*spc).sub[bb1][ncov].dim1;j++){
            (*spc).sub[bb1][ncov].ktsc[j]=(*spc).sub[bb1][ncov].ktsc[j+1];
         }
      }
      else{
/* if it is a two variable dimension */
         (*spc).sub[bb1][bb2].dim1-=1;
      }
   }

/* or merge two links */
   else {
      i=(*spc).basis[bbi].link1[baj];
      if(bbj==bbk){
/* a merger with emptiness */
         for(j=0;j<nclass;j++){
            if((*spc).basis[bbi].link1[j]==i){
               (*spc).basis[bbi].link1[j]= -1;
               (*spc).basis[bbi].beta[j]=0.;
            }
            if((*spc).basis[bbi].link1[j]>i) (*spc).basis[bbi].link1[j]-= 1;
         }
      }
      else{
/* a true merger */
         k=(*spc).basis[bbi].link1[bak];
         if(k>i){
             j=k;
             k=i;
             i=j;
         }
         xx=(*spc).basis[bbi].beta[bak];
         for(j=0;j<nclass;j++){
            if((*spc).basis[bbi].link1[j]==k) (*spc).basis[bbi].beta[j]=xx;
            if((*spc).basis[bbi].link1[j]==i){
               (*spc).basis[bbi].beta[j]=xx;
               (*spc).basis[bbi].link1[j]=k;
            }
            if((*spc).basis[bbi].link1[j]>i) (*spc).basis[bbi].link1[j]-= 1;
         }
      }
      (*spc).basis[bbi].ib -= 1;
   }
   (*iwald)=(*iwald)+1;
   (*dwald)=(*dwald)+criterion*criterion;

/* initialize link2 */
   l= -1;
   for(j=0;j<(*spc).nbas;j++){
      for(k=0;k<(*spc).basis[j].ib;k++){
         l++;
         (*spc).basis[j].link2[k]=l;
      }
   }
   for(i=0;i<nbas;i++)for(j=0;j<nclass;j++){
      k1=dlink(spc,i,j);
      if(k1>=0){
         k2=tlink[i][j];
         if(bsk!= -1&&k2==bsj)k2=bsk;
         for(i1=0;i1<nbas;i1++)for(j1=0;j1<nclass;j1++){
            k3=dlink(spc,i1,j1);
            if(k3>=0){
               k4=tlink[i1][j1];
               if(bsk!= -1&&k4==bsj)k4=bsk;
               (*spc).info[k1][k3]=tinfo[k2][k4];
            }
         }
      }
   }
   return;
}
/******************************************************************************/
static void puuu(spc,b1,b2,t1,t2,ncov,ii)
struct space *spc;
int b1,b2,t1,t2,ncov,ii;
{
   if(ii==0)(void)Rprintf("   add: ");
   if(ii==1)(void)Rprintf("remove: ");
   if(ii==2)(void)Rprintf(" merge: ");
   if(b1!=ncov){   
      (void)Rprintf("cov(%d",b1+1);
      if(b2==ncov)(void)Rprintf(")=(");
      else (void)Rprintf(",%d)=(",b2+1);
      if(t1!= -1)(void)Rprintf("knot=%.2f",(*spc).sub[b1][ncov].ktsc[t1]);
      else (void)Rprintf("linear");
      if(b2==ncov)(void)Rprintf(") ");
      else {
         if(t2!= -1)(void)Rprintf(",%.2f)",(*spc).sub[b2][ncov].ktsc[t2]);
         else (void)Rprintf(",linear) ");
      }
   }
   else (void)Rprintf("constant ");
}
/******************************************************************************/

/* This program does the main control */

static void poly(best,data,loss,pen,ndmax,mind,exclude,strt,silent,logs,ad,lins,
    tdata,it,aics,current,new,trynew,naction,il,xsingle,newx)
struct space *best,*current,*new,*trynew,*newx;
struct datastruct *data,*tdata;
double **loss,**logs,*aics,pen;
int ndmax,mind,**exclude,strt,silent,*ad,*lins,it,naction,il,xsingle;

/* best        - the best model up to now
   data        - the data 
   loss        - the loss function. 
   pen         - penalty 
   ndmax       - maximum number of basisfunctions
   mind        - minimum distance (in order statistics) between knots 
   naction     - number of possible different actions
   exclude     - which terms should be excluded from the model 
   strt        - 0: start with constant, 1 start with linear, -1 start with fit
   silent      - should diagnostic output be printed? 
   ad          - is the best model during addition or deletion 
   lins        - dimensions that can only be added linear 
   tdata       - testset data
   it          - are we using a testset or aic? 
   aics        - all the losses */
{
   double dwald;
   int add=1,i,ndm2,iwald,okd;

/* getcrit       - computes the criterion (AIC or loss)
   pnewton     - fits a model using NR
   current       - the present model
   new,trynew    - storage for a space, used by padddim, prembas and premdim
   pdefinespace-allocates storage for a space
   add           - are we still adding?
   padddim     - adds dimensions to a space
   i             - counter
   ndm2          - copy of ndmax on entrance
   ndim          - save typing
   pconstant   - initializes a constant hazard space
   pswapspace  - copies one space into another 
   prembas     - remove basisfunctions from a space  
   premdim     - remove dimensions from a space  */

/* swap spaces */
   pswapspace(current,best,data);

/* initialization */
   if(silent==0 && it==0)(void)Rprintf
      ("dim        AIC log-likeli log-like/n resub-ls/n sq-error/n\n");
   if(silent==0 && it!=0){
      (void)Rprintf("dim  measure <== training set / n ==>");
      (void)Rprintf(" <==== test set / n ====>\n");
      if(il==0)(void)Rprintf("    log-like ");
      if(il==1)(void)Rprintf("        loss ");
      if(il==2)(void)Rprintf("    sq-error ");
      (void)Rprintf("log-like    loss  sq-err ");
      (void)Rprintf("log-like    loss  sq-err\n");
   }
   ndm2=ndmax;
   if(ndmax<0)ndmax= -ndmax;
   (*best).aic=pow((double)10.,(double)150.);
   for(i=0;i<maxdim;i++)logs[0][i]= -pow((double)10.,(double)150.);

/* specifies constant model */
   if(strt>=0)pconstant(current,data);

/* we start in adding mode */
   do{

/* fits the model */
      (*current).aic=pnewton(current,data);
      if((*current).aic>190){
         pswapspace(current,newx,data);
         add=(int)pnewton(current,data);
         ndmax=(*current).nbas;
         add=0;
         Rprintf("warning - model size was reduced\n");fflush(stdout);
      }
      else{

/* compute aic */
         (*current).aic=getcrit(current,tdata,data,it,loss,silent,logs,ad,1,
              aics,&pen,naction,il);
         if((*current).ndim>=ndmax-(*data).nclass+1)add=0;

/* did we improve */
         if((*current).aic<=(*best).aic+0.00000001){
            pswapspace(best,current,data);
            if(silent==0)(void)Rprintf(" best up to now!");
         }
         if(silent==0)(void)Rprintf("\n");
         if(silent==0)(void)fflush(stdout);

/* adds dimensions, computes new starting values */
         if(add==1 && ndm2<0){
            for(i=2;i<(*current).nbas-2;i++){
               if(logs[0][(*current).ndim-1]-
                  logs[0][i-1]<((*current).nbas-i)/2.-0.5){
                  add=0;
                  ndmax=(*current).nbas;
               }
            }
         }
         if(add==1){
            pswapspace(newx,current,data);
            add=padddim(current,new,trynew,data,mind,exclude,silent,lins);
            if(add!=1) ndmax=(*current).nbas;
         }
         if(silent==0)(void)fflush(stdout);
     }

/* keep on adding? */
   }while(add==1);

/* start deleting */
   if(xsingle!=0)do{

/* removes dimensions, computes new starting values */
      if(xsingle==1){
         dwald=0.;
         iwald=0;
         do{
            if(ndmax>1) premdim(current,data,silent,&dwald,&iwald);
         }while(iwald< 10 && dwald<2.5 && (*current).ndim>25);
      }
      else{
         okd = prembas(current,data,silent);
      }
      if(okd!=-1){
      if(silent==0)(void)fflush(stdout);
      (*current).aic=pnewton(current,data);
      if((*current).aic > 190.)add=17;
      else{

/* compute aic */
         (*current).aic=getcrit(current,tdata,data,it,loss,silent,logs,ad,0,
            aics,&pen,naction,il);

/* did we improve */
         if((*current).aic<=(*best).aic+0.00000001){
            pswapspace(best,current,data);
            if(silent==0)(void)Rprintf(" best up to now!");
         }
         if(silent==0)(void)Rprintf("\n");
         if(silent==0)(void)fflush(stdout);
      }
      }
/* does further deleting make sense */
   }while(okd!=-1 && (
    ((*current).aic-(*best).aic< -pen*((*current).ndim-(*data).nclass)&&it==0)||
     ((*current).ndim>(*data).nclass && it!=0 && add!=17)));
}

/******************************************************************************/

/* this function initializes a constant hazard space */

static void pconstant(spc,data)
struct space *spc;
struct datastruct *data;

/* spc   - space to be initialized
   data  - the data  */
{
   int i,ncov=(*data).ncov,nclass=(*data).nclass;

/* i - counter
   ncov,nclass,ndata - save typing */

   (*spc).ndim=(*data).nclass;
   (*spc).nbas=1;
   (*spc).basis[0].b1=ncov;
   (*spc).basis[0].b2=ncov;

/* initialize the values and the starting beta */
   for(i=0;i<nclass;i++) (*spc).basis[0].beta[i]= 0.;
}

/******************************************************************************/

/* this function copies one space into another - just element by element */

static void pswapspace(spout,spin,data)
struct space *spin,*spout;
struct datastruct *data;

/* spin  - input space
   spout - output space
   data */

{
   int i,j,k,l,m;
   int ncov=(*data).ncov,nclass=(*data).nclass;

/* i,j,k,l - counters 
   ndata - number of datapoints
   nclass- number of classes
   ncov  - number of covariates */

/* the space-core part */
   (*spout).ndim=(*spin).ndim;
   (*spout).nbas=(*spin).nbas;
   (*spout).aic=(*spin).aic;
   for(i=0;i<(*spin).ndim;i++){
      (*spout).score[i]=(*spin).score[i];
      for(j=0;j<(*spin).ndim;j++) (*spout).info[i][j]=(*spin).info[i][j];
   }

/* the basis-functions part */
   for(i=0;i<(*spin).nbas;i++){
      for(j=0;j<nclass;j++)(*spout).basis[i].beta[j]=(*spin).basis[i].beta[j];
      (*spout).basis[i].b1=(*spin).basis[i].b1;
      (*spout).basis[i].b2=(*spin).basis[i].b2;
      (*spout).basis[i].t1=(*spin).basis[i].t1;
      (*spout).basis[i].t2=(*spin).basis[i].t2;
      (*spout).basis[i].ib=(*spin).basis[i].ib;
      for(j=0;j<nclass;j++){
         (*spout).basis[i].link1[j]=(*spin).basis[i].link1[j];
         (*spout).basis[i].link2[j]=(*spin).basis[i].link2[j];
      }
   }

/* the subdimensions part, first the 2-covariate ones */
   m=MAXKNOTS+1;
   if((*spin).nbas<m)m=(*spin).nbas;
   for(i=0;i<ncov;i++){
      for(j=i+1;j<ncov;j++){
         (*spout).sub[i][j].dim1=(*spin).sub[i][j].dim1;
         if((*spout).sub[i][j].dim1>0){
            for(k=0;k<m;k++){
               for(l=0;l<m;l++){
                  (*spout).sub[i][j].kts1[k][l]=(*spin).sub[i][j].kts1[k][l];
               }
            }
         }
      }
   }

/* the subdimensions part, the 1-covariate ones */
   for(j=0;j<ncov;j++){
      (*spout).sub[j][ncov].dim1=(*spin).sub[j][ncov].dim1;
      for(k=0;k<(*spout).sub[j][ncov].dim1-1;k++){
         (*spout).sub[j][ncov].ktsc[k]=(*spin).sub[j][ncov].ktsc[k];
      }
   }
}
/******************************************************************************/
/* computes the loss for the testset case */
static void computeloss(spc,data,loss,naction,res)
struct datastruct *data;
struct space *spc;
double **loss,res[3];
int naction;

/* spc  - the model
   data - the data
   naction - number of different actions
   loss - the lossfunction */

{
   int i,j,k,i1,i2,ncov=(*data).ncov,nclass=(*data).nclass,ll;
   double *chances,value2,value;
   float *cc;
   
/* i,j,k        - counters
   ncov, nclass, i1, i2 - save typing
   loser        - the loss
   chances      - log(probability)
   value2       - piecewise constant #2 */

   for(i=0;i<3;i++)res[i]=0;
   chances=v1;
   ll=1;
   if(naction==nclass) for(i=0;i<nclass;i++) for(j=0;j<nclass;j++){
      if(i==j && fabs(loss[i][j])>0.0000001)ll=0;
      if(i!=j && fabs(loss[i][j]-1)>0.0000001)ll=0;
      if(ll==0){i=nclass;j=nclass;}
   }

/* circle the data points */
   for(i=0;i<(*data).ndata;i++){
      j=(*data).icov[i];
      if(j>0)cc= &(trcov[j-1]);
      else cc= &(tecov[-j-1]);

/* initialize */
      for(j=0;j<=nclass;j++) chances[j]=0.;

/* circle the basis functions */
      for(k=0;k<(*spc).nbas;k++){
/* compute the basis function: the first component */
         if(k>0){
            value=1;
            i1=(*spc).basis[k].t1;
            i2=(*spc).basis[k].b1;
            if(i1== -1)value=cc[i2];
            else{
               value=cc[i2]-(*spc).sub[i2][ncov].ktsc[i1];
               if(value<0) value=0.;
            }
/* the second component */
            i2=(*spc).basis[k].b2;
            if(i2!=ncov && value!=0.){
               value2=1.;
               i1=(*spc).basis[k].t2;
               if(i1== -1)value2=cc[i2];
               else{
                  value2=cc[i2]-(*spc).sub[i2][ncov].ktsc[i1];
                  if(value2<0)value2=0.;
               }
               value= value*value2;
            }
            for(j=0;j<nclass;j++) chances[j]+=(*spc).basis[k].beta[j]*value;
         }
         else for(j=0;j<nclass;j++) chances[j]+=(*spc).basis[k].beta[j];
      }
/* make it probabilities */
      for(j=0;j<=nclass;j++){
         if(chances[j]<600.)chances[j]= exp(chances[j]);
         else chances[j]=exp((double)600.);
      }
      value=0.;
      for(j=0;j<=nclass;j++)value+=chances[j];
      for(j=0;j<=nclass;j++)chances[j]=chances[j]/value;
/* who wins */
      res[0]+=log(chances[(*data).yy[i]])*(*data).wgt[i];
      k=0;
      value=0;
      for(j=0;j<naction;j++){
         value2=0;
         for(i1=0;i1<=nclass;i1++) value2+=chances[i1]*loss[j][i1];
         if(value2<value || j==0){
            k=j;   
            value=value2;
         }
      }
      res[1]+=loss[k][(*data).yy[i]]*(*data).wgt[i];
      res[2]+=(1.-chances[(*data).yy[i]])*(1.-chances[(*data).yy[i]])
                                                            *(*data).wgt[i];
   }
}
/******************************************************************************/
/* this function computes the criterion and does some other bookkeeping */
static double getcrit(spc,tdata,data,it,loss,silent,logs,ad,ii,aics,pen,naction,il)
int it,silent,*ad,ii,naction,il;
double **loss,**logs,*aics,*pen;
struct space *spc;
struct datastruct *tdata,*data;

/* spc   - the space
   tdata - test data
   it    - aic(0) or testset 
   loss  - loss matrix
   pen   - penalty
   silent- print diagnostic output?
   logs  - log-likelihoods
   ad    - was model fit during addition or deletion?
   naction - number of possible actions
   ii    - 0 (deletion stage) or 1 (addition stage) */
{
   double crit,res[3],tes[3];
   int i;
/* crit  - the criterion
   computeloiss - computes it in the case of a testset */
   
   computeloss(spc,data,loss,naction,res);
   if(it!=0)computeloss(spc,tdata,loss,naction,tes);
   if(it==0)crit=(*spc).aic;
   else crit=tes[il];
   if(it==0)crit=2*crit+(*pen)*(*spc).ndim;
   if((crit<aics[(*spc).ndim-1]&&it!=0)||
      (crit>aics[(*spc).ndim-1]&&it==0) || ad[(*spc).ndim-1]== -1 || ii==1){
      ad[(*spc).ndim-1]=ii;
      aics[(*spc).ndim-1]=crit;
      for(i=0;i<3;i++)logs[i][(*spc).ndim-1]=res[i];
      if(it==0)for(i=3;i<6;i++)logs[i][(*spc).ndim-1]=0;
      else for(i=3;i<6;i++)logs[i][(*spc).ndim-1]=tes[i-3];
   }
   if(silent==0){
      if(it==0)
         (void)Rprintf("%3d %10.4f %10.4f %10.4f %10.4f %10.4f ",
           (*spc).ndim,-crit,res[0],res[0]/(*data).wgtsum,res[1]/(*data).wgtsum,
            res[2]/(*data).wgtsum);
      else{
         (void)Rprintf("%3d %8.3f %8.3f %7.3f %7.3f ",(*spc).ndim,crit,
            res[0]/(*data).wgtsum,res[1]/(*data).wgtsum,res[2]/(*data).wgtsum);
         (void)Rprintf("%8.3f %7.3f %7.3f ",
          tes[0]/(*tdata).wgtsum,tes[1]/(*tdata).wgtsum,tes[2]/(*tdata).wgtsum);
      }
   }
   if(it!=0 && il==0)crit= -crit;
   if(it==0)crit= -crit;
   return crit;
}
/******************************************************************************/
/* checks a number */
static int plumbertester(aa)
double aa;
/* if aa = -Inf: 0
      aa = +Inf: 1
      aa =  NaN: 2
      otherwise: 3 */
{
   int i1=0,i2=0,i3=0,i4=0;
   if(aa< 2.)i1=1;
   if(aa> 0.)i2=1;
   if(aa< pow(10.,200.))i3=1;
   if(aa> -pow(10.,200.))i4=1;
   if(i1+i2+i3+i4>=3)return 3;
   if(i2==1 && i4==1)return 1;
   if(i1==1 && i3==1)return 0;
   return 2;
}
/******************************************************************************/
/* sort, put result in rb */
static void Ppsort(ra,n)
int n;
double *ra;
{
   xpsort(ra-1,n);
}
/******************************************************************************/
/* sort */
static void xpsort(ra,n)
int n;
double *ra;
{
   int l,j,ir,i;
   double rra;

   l=(n >> 1)+1;
   ir=n;
   for (;;) {
      if (l > 1) rra=ra[--l];
      else {
          rra=ra[ir];
         ra[ir]=ra[1];
         if (--ir == 1) {
            ra[1]=rra;
            return;
         }
      }
      i=l;
      j=l << 1;
      while (j <= ir) {
         if (j < ir && ra[j] < ra[j+1]) ++j;
         if (rra < ra[j]) {
            ra[i]=ra[j];
            j += (i=j);
         }
         else j=ir+1;
      }
      ra[i]=rra;
   }
}
/******************************************************************************/
static int lusolinv(a,n,b,k)
int n,k;
double **a,*b;
/* various lu things
   k=0 inverse, non symmetric
   k=1 inverse, symmetric
   k=2 solve, symmetric */
{
   double aa[DIM5][DIM5],bb[DIM5],det[2];
   int kpvt[DIM5],info,i,j,inert[3];
   if(k<2) for(i=0;i<n;i++) for(j=0;j<n;j++)aa[i][j]=a[j][i];
   else for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
      bb[i]=b[i];
   }
   i=DIM5;
   j=1;
   if(k==0){
      F77_CALL(xdgefa)(aa,&i,&n,kpvt,&info);
      F77_CALL(xdgedi)(aa,&i,&n,kpvt,det,bb,&j);
      for(i=0;i<n;i++) for(j=0;j<n;j++)a[i][j]=aa[j][i];
   }
   if(k==1){
      F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
      F77_CALL(xdsidi)(aa,&i,&n,kpvt,det,inert,bb,&j);
      for(i=0;i<n;i++) for(j=i;j<n;j++)a[i][j]=aa[j][i];
      for(i=0;i<n;i++) for(j=0;j<i;j++)a[i][j]=aa[i][j];
   }
   if(k==2){
      F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
      if(info!=0)return 0;
      F77_CALL(xdsisl)(aa,&i,&n,kpvt,bb);
      for(i=0;i<n;i++)b[i]=bb[i];
   }
   return 1;
}
/******************************************************************************/
