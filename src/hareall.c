/*
*  Copyright (C) 1993--2002  Charles Kooperberg
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
/* this file contains the main body of the program, a few small routines on 
   which it depends and two routines that deal with space-structures */

/* this function describes the basis data structure */

#include <math.h>
#include <stdio.h>
#include "R.h"
#define Salloc(n, t)  (t *)R_alloc((long)(n), (int)sizeof(t))

/* we want to be able to use those everywhere */

#define MAXSPACE 53
#define MAXKNOTS 1
#define DIM5 MAXSPACE+5

/* MAXSPACE - maximum dimensionality of the model
   MAXKNOTS - maximum number of knots for one covariate */

void F77_NAME(xdsifa)(double[][DIM5], int *, int *, int *, int *);
void F77_NAME(xdsisl)(double[][DIM5], int *, int *, int *, double *);
void F77_NAME(xdsidi)(double[][DIM5], int *, int *, int *, double *, int *, double *, int *);
void F77_NAME(xssort)(double *, double *, int *, int *);
void F77_NAME(xdsico)(double[][DIM5], int *, int *, int *,  double *, double *);


struct datastruct {
   int ndata,ncov,*delta,*bincov,*same;
   double *times,**cov;
};

/* datastruct is a structure containing all information about the data. At any 
   time there is only one datastruct, which is called data.

   ndata  - number of datapoints 
   ncov   - number of covariates 
   delta  - delta (censor indiactor) 0=censored 1=uncensored 
   bincov - are the covariates binary? 0=no, 1=yes binary cov should be 0-1 
   same   - if the covariate structure of case i is the same as for case i-1 
   times  - observation/censoring times 
   cov    - covariates cov[i][j] is covariate j for observation i */

struct space {
   int ndim,nknots;
   double *knots,aic,**info,*score,**b0,**b1,*b2,**xtx;
   struct basisfunct *basis; 
   struct subdim **sub; 
};

/* space is the basic structure containing a model. The main ingredients are a
   (sort of double) representation of the basisfunctions: by means of basis on
   a basisfunction by basisfunction scale and by means of sub on a subdimension
   scale

   ndim    - the dimensionality of the space
   nknots  - the number of time-knots
   knots   - the time-knots
   aic     - the aic value of the present model - only accurate after the model
             has been fitted.
   info    - the hessian of the model
   score   - the score vector of the model
   b0      - first element constant term of lambda[between-2-knots,datapoint]
             others constant term basisfunction(j)[between-2-knots,datapoint]
   b1      - b1, but linear term
   b2      - first element, lambda in a datapoint others, basisfunction(j) in a
             datapoint
   basis   - the array of basisfunctions
   sub     - the matrix of subdimensions element:
             [i][data.ncov]  (0<=i<data.ncov) the subdimension belonging to a
                                              covariate by itself.
             [data.ncov][data.ncov]           the subdimension belonging to time
             [data.ncov][i]  (0<=i<data.ncov) the subdimension belonging to the
                                              tensor-product between time and
                                              a covariate
             [i][j]  (0<=i<j<data.ncov)       the subdimension belonging to the
                                              tensor-product of two covariates
             others                           not used. */

struct basisfunct {
   int b1,b2,t1,t2,iknots;
   double beta,*values,*values2,se;
};

/* a structure describing one basisfunction (surprising isn't it?)
   b1 and b2 - indicate the subdimension to which this basisfunction belongs
   t1 and t2 - the rank of the knots in subdimensions b1 and b2, if applicable
               if these are -1, while b1 and b2 would indicate a covariate, the
               basisfunction is linear in this covariate.
   iknots    - as t1 or t2, but for the time-knot
   beta      - beta of this basisfunction
   values    - value of the not-time part of the basis-function in datapoints 
   values2   - value of the complete basis-function in datapoints */

struct subdim {
   short int dim1,**kts1;
   float *ktsc;
};

/* a structurte describing a subdimension - it can depend on either one
   or two covariates or one covariate and time.
   dim1   - dimensionality of component 
   kts1   - binary - if =0 indicates that basisfunction is not in the model,
            otherwise number of basisfunction.
   ktsc   - if this is a subdimension depending on 1 covariate, the knots */

/* allocation stuff used a lot */

static int glusolve2(),humbertester();
static void nrerror(),glusolve();
static int *newtonwhere;
static double *searchsorted,*searchkts,*searchsorted2,*remdimy,*raoss;
static double *raoscorecopy,*newtonscp,*compallss,*complogbasis0,*complogbasis1;
static double *remdimxty,**raohhh,**getsexx,**compallhhh,**remdimxtx;
static double *dgvector(),**dgmatrix(),testbasis();
static float *ddgvector();
static short **iigmatrix(),*iigvector();
static int *igvector(),**igmatrix(),zlocation();
static double newton(),adders(),search(),eint(),xeint();
static struct space *definegspace(),*hdefinegspace();
static int gadddim(),gindl(),gindr(),gindx(),gindm(),gindyl(),gindyr();
static void constant(),swapgspace(),gremdim(),gluinverse(),getse();
static void soutgspace(),houtgspace(),poutgspace(),uuu(),sort(),dsort();
static double critswap();
static void addbasis(),tswapout(),upbasis(),veint(),upbasis2();
static int tswapin();
static double fct1(),fct2(),grao(),condition(),compall(),complog(),hcomplog();
static void cleanupt(),cleanup1(),basisswap();
static struct datastruct *definedata();
static void getvector(),hareallocer();
static struct basisfunct *definebasis();
static struct subdim **definedim();
static struct basisfunct *hdefinebasis();
static struct subdim **hdefinedim();
static void getvectors(),getthosep(),upbasis3();

/******************************************************************************/

/* This program does the main control */
static void hare(best,data,alpha,ndmax,mind,exclude,strt,silent,logs,ad,lins)
struct space *best;
struct datastruct *data;
double alpha,*logs;
int ndmax,mind,**exclude,strt,silent,*ad,*lins;

/* best        - the best model up to now
   data        - the data 
   alpha       - alpha for bic
   ndmax       - maximum number of basisfunctions
   mind        - minimum distance (in order statistics) between knots 
   exclude     - which terms should be excluded from the model 
   strt        - 0: start with constant, 1 start with linear, -1 start with fit
   silent      - should diagnostic output be printed? 
   ad          - is the best model during addition or deletion 
   lins        - dimensions that can only be added linear */
{
   struct space *current,*new,*trynew;
   int add=1,i,oops=0,ndm2;

/* newton      - fits a model using NR
   ndm2        - copy of ndmax on entrance
   add         - are we still adding?
   gadddim      - adds dimensions to a space
   current     - the present model
   new,trynew  - storage for a space, used by gadddim and remdim
   definegspace-allocates storage for a space
   constant    - initializes a constant hazard space
   swapgspace   - copies one space into another 
   gremdim      - remove dimensions from a space 
   getse       - get ses of the coefficients */

/* allocates storage for spaces */
   new=definegspace((*data).ncov,(*data).ndata);
   trynew=definegspace((*data).ncov,(*data).ndata);
   current=definegspace((*data).ncov,(*data).ndata);
   if(strt<0)swapgspace(current,best,(*data).ndata,(*data).ncov);


/* initialization */
   ndm2=ndmax;
   if(ndmax<0)ndmax=-ndmax;
   (*best).aic=pow((double)10.,(double)150.);
   for(i=0;i<MAXSPACE;i++)logs[i]=-pow((double)10.,(double)150.);

/* specifies constant model or linear, plus starting values */
   if(strt>=0)constant(current,data,strt);

/* we start in adding mode */
   do{

/* fits the model */
      (*current).aic=newton(current,data,0,silent,&oops);

/* compute aic */
      logs[(*current).ndim-1]=(*current).aic;
      ad[(*current).ndim-1]=1;
      (*current).aic=(*current).ndim*alpha-2*(*current).aic;
      if((*current).ndim==ndmax)add=0;

/* did we improve */
      if((*current).aic<(*best).aic){
         getse(current);
         swapgspace(best,current,(*data).ndata,(*data).ncov);
      }

/* adds dimensions, computes new starting values */
      if(add==1 && ndm2<0){
         for(i=2;i<(*current).ndim-2;i++){
            if(logs[(*current).ndim-1]-logs[i-1]<((*current).ndim-i)/2.-0.5){
               add=0;
               ndmax=(*current).ndim;
            }
         }
      }
      if(add==1){
         add=gadddim(current,new,trynew,data,mind,exclude,silent,lins);
         if(add!=1) ndmax=(*current).ndim;
      }

/* keep on adding? */
   }while(add==1);

/* the last addition space is the first best space */
   (*current).aic=newton(current,data,1,silent,&oops);
   logs[(*current).ndim-1]=(*current).aic;
   (*current).aic=(*current).ndim*alpha-2*(*current).aic;

/* start deleting */
   do{

/* removes dimensions, computes new starting values */
      if(ndmax>1)gremdim(current,(*data).ncov,(*data).ndata,silent);
      if((*best).ndim==(*current).ndim+1 && ad[(*best).ndim-1]==0){
         for(i=0;i<(*best).ndim;i++){
            (*best).basis[i].se=(*current).basis[i].se;
         }
      }
      (*current).aic=newton(current,data,1,silent,&oops);
      if(oops==1)(*current).aic=newton(current,data,2,silent,&oops);

/* compute aic */
      if((*current).aic>logs[(*current).ndim-1]){
         logs[(*current).ndim-1]=(*current).aic;
         ad[(*current).ndim-1]=0;
      }
      (*current).aic=(*current).ndim*alpha-2*(*current).aic;

/* did we improve */
      if((*current).aic<(*best).aic){
         swapgspace(best,current,(*data).ndata,(*data).ncov);
      }

/* does further deleting make sense */
   }while((*current).aic-(*best).aic<alpha*((*current).ndim-1));
   if((*best).ndim==(*current).ndim && ad[(*best).ndim-1]==0){
      gluinverse((*current).info,(*best).ndim);
      for(i=0;i<(*best).ndim;i++)
                            (*best).basis[i].se=sqrt(-(*current).info[i][i]);
   }
}

/******************************************************************************/

/* this function initializes a constant hazard space */

static void constant(spc,data,strt)
struct space *spc;
struct datastruct *data;
int strt;

/* spc   - space to be initialized
   data  - the data 
   strt  - start with a constant (0) or linear (1) */
{
   int i,j=0,k;
   double r=0.;

/* i,k - counter
   j - counts the number of deltapoints
   r - counts the data */

   (*spc).ndim=1;
   (*spc).nknots=0;
   (*spc).basis[0].b1=(*data).ncov;
   (*spc).basis[0].b2=(*data).ncov;

/* initialize the values and compute the starting beta */
   for(i=0;i<(*data).ndata;i++){
      (*spc).basis[0].values[i]=1.;
      (*spc).basis[0].values2[i]=1.;
      (*spc).xtx[0][0]=(double)(*data).ndata;
      r+=(*data).times[i];
      j+=(*data).delta[i];
   }
   (*spc).basis[0].beta= -log(r/(double)j);

/* add the linear spaces */
   if(strt!=0){
      for(i=0;i<(*data).ncov;i++){
         (*spc).ndim++;
         (*spc).basis[i+1].b1=i;
         (*spc).basis[i+1].b2=(*data).ncov;
         (*spc).basis[i+1].beta=0.;
         (*spc).basis[i+1].se=0.;
         (*spc).sub[i][(*data).ncov].dim1=1;
         for(j=0;j<(*data).ndata;j++){
            (*spc).basis[i+1].values[j]=(*data).cov[i][j];
            (*spc).basis[i+1].values2[j]=(*data).cov[i][j];
         } 
         (*spc).xtx[0][i+1]=0.;
         (*spc).xtx[i+1][0]=0.;
         for(j=0;j<(*data).ndata;j++){
            (*spc).xtx[0][i+1]+=(*spc).basis[i+1].values2[j];
         }
         (*spc).xtx[0][i+1]=(*spc).xtx[i+1][0];
         for(j=0;j<(*data).ncov;j++){
            (*spc).xtx[i+1][j+1]=0.;
            for(k=0;k<(*data).ndata;k++){
               (*spc).xtx[0][i+1]+=(*spc).basis[i+1].values2[k]*
                                   (*spc).basis[j+1].values2[k];
            }
         }
      }
   }
}

/******************************************************************************/

/* this function copies one space into another - just element by element */

static void swapgspace(spout,spin,ndata,ncov)
struct space *spin,*spout;
int ndata,ncov;

/* spin  - input space
   spout - output space
   ndata - number of datapoints
   ncov  - number of covariates */

{
   int i,j,k,l,m;

/* i,j,k,l - counters */

/* the space-core part */
   (*spout).ndim=(*spin).ndim;
   (*spout).nknots=(*spin).nknots;
   (*spout).aic=(*spin).aic;
   for(i=0;i<(*spin).nknots;i++)  (*spout).knots[i]=(*spin).knots[i];
   for(i=0;i<=(*spin).ndim;i++){
      (*spout).b2[i]=(*spin).b2[i];
      for(j=0;j<=(*spin).nknots+1;j++){
         (*spout).b0[j][i]=(*spin).b0[j][i];
         (*spout).b1[j][i]=(*spin).b1[j][i];
      }
   }
   for(i=0;i<(*spin).ndim;i++){
      (*spout).score[i]=(*spin).score[i];
      for(j=0;j<(*spin).ndim;j++) (*spout).info[i][j]=(*spin).info[i][j];
   }

/* xtx */
   for(i=0;i<(*spin).ndim;i++)
      for(j=0;j<(*spin).ndim;j++)
         (*spout).xtx[i][j]=(*spin).xtx[i][j];

/* the basis-functions part */
   for(i=0;i<(*spin).ndim;i++){
      (*spout).basis[i].beta=(*spin).basis[i].beta;
      (*spout).basis[i].se=(*spin).basis[i].se;
      (*spout).basis[i].iknots=(*spin).basis[i].iknots;
      (*spout).basis[i].b1=(*spin).basis[i].b1;
      (*spout).basis[i].b2=(*spin).basis[i].b2;
      (*spout).basis[i].t1=(*spin).basis[i].t1;
      (*spout).basis[i].t2=(*spin).basis[i].t2;
      for(j=0;j<ndata;j++){
         (*spout).basis[i].values[j]=(*spin).basis[i].values[j];
         (*spout).basis[i].values2[j]=(*spin).basis[i].values2[j];
      }
   }

/* the subdimensions part, first the 2-covariate ones */
   for(i=0;i<ncov;i++){
      for(j=i+1;j<ncov;j++){
         (*spout).sub[i][j].dim1=(*spin).sub[i][j].dim1;
         if((*spout).sub[i][j].dim1>0){
            for(k=0;k<(*spout).sub[i][ncov].dim1+1;k++){
               for(l=0;l<(*spout).sub[j][ncov].dim1+1;l++){
                  (*spout).sub[i][j].kts1[k][l]=(*spin).sub[i][j].kts1[k][l];
               }
            }
         }
      }
   }

/* the subdimensions part, the 1-covariate and time ones */
   m=MAXKNOTS+1;
   if((*spin).ndim<m)m=(*spin).ndim;
   for(j=0;j<ncov;j++){
      (*spout).sub[ncov][j].dim1=(*spin).sub[ncov][j].dim1;
      if((*spout).sub[ncov][j].dim1>0){
         for(k=0;k<m;k++){
            for(l=0;l<m;l++){
               (*spout).sub[ncov][j].kts1[k][l]=(*spin).sub[ncov][j].kts1[k][l];
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

/* this routine searches all dimensions for something to add */

static int gadddim(current,new,newt,data,mind,exclude,silent,lins)
struct space *current,*new,*newt;
struct datastruct *data;
int mind,**exclude,silent,*lins;

/* current - current space
   new     - copy of current space to play with
   newt    - best addition space up to now
   data    - data
   mind    - minimum distance (order statistics) between knots 
   exclude - if exclude = 1 this term is never included 
   lins    - which dimensions should onle be added linear
   silent  - should diagnostic output be printed */

{
   int i,j;
   double criterion;

/* i,j       - counters
   criterion - criterion
   adders    - does the work
   swapspace - copies one space into another 
   uuu       - print some diagnostic statistics */

/* initialization */
   criterion = -pow((double)10.,(double)20.);

/* get a space to play with */
   swapgspace(newt,current,(*data).ndata,(*data).ncov);
   
/* searches all the dimensions */
   criterion=adders((*data).ncov,(*data).ncov,current,new,newt,criterion,
                                                        data,mind,lins);
   for(i=0;i<(*data).ncov;i++){
      for(j=i+1;j<=(*data).ncov;j++){
         if(i==(*data).ncov || i!=j){
            if(j==(*data).ncov || exclude[i][j]==0){
               criterion=adders(i,j,current,new,newt,criterion,data,mind,lins);
            }
         }
      }
      if(exclude[(*data).ncov][i]==0)
         criterion=adders((*data).ncov,i,current,new,newt,criterion,data,mind,
                                                                    lins);
   }

/* copy the result success */
   if(criterion>0.){
      swapgspace(current,new,(*data).ndata,(*data).ncov);
      i=(*current).ndim-1;
      if(silent!=1)uuu(current,(*current).basis[i].b1,(*current).basis[i].b2,
         (*current).basis[i].t1,(*current).basis[i].t2,(*data).ncov,0);
      if(silent!=1)(void)Rprintf("(rao= %.2f) ",criterion);
      return 1;
   }
/* failure */
   else{
      return 0;
   }
}

/******************************************************************************/

/* this routine searches a subdimension for a supspace to add */

static double adders(i0,j0,current,new,newt,criterion,data,mind,lins)
int i0,j0,mind,*lins;
struct space *new,*newt,*current;
struct datastruct *data;
double criterion;

/* i0,j0     - which subspace (see hstruct)
   current   - sometimes we need two of them (see newt)
   new       - will be the best space with additions up to now.
   newt      - actually a copy of current, we play with it until we are done
   criterion - the best rao statistic (chi-square p-value) up to now
   data      - structure containing the data 
   mind      - minimum distance (in order statistics) between knots */

{
   int i,j;
   double  crit1;

/* swapspace- copies one space into another
   i,j      - counter
   testbasis- does the work for 2d dimensions
   search   - does the work for 1d dimensions
   crit1    - possibly optimal criterion */

/* a 1-d space */
   if(j0==(*data).ncov){
      if(i0==(*data).ncov){
/* search for a t-knot */
         if((*newt).nknots<MAXKNOTS && lins[i0]==0){
            criterion=search(new,newt,data,i0,mind);
         }
      }

/* a covariate that has not yet been entered */
      else{
         if((*newt).sub[i0][j0].dim1==0){
           criterion=testbasis(new,newt,criterion,data,i0,j0,0,-1,(double)0);
         }

/* a covariate that has been entered before */
         else{
            if((*data).bincov[i0]==0){
               if((*newt).sub[i0][j0].dim1<MAXKNOTS){
                  if(lins[i0]==0)crit1=search(current,newt,data,i0,mind);
                  if(crit1>criterion && lins[i0]==0){
                     criterion=crit1;
                     swapgspace(new,current,(*data).ndata,(*data).ncov);
                  }
               }
            }
         }
      }
   }

/* a 2-d space */
   else{

/* covariate x time */
      if(i0==(*data).ncov && (*newt).sub[j0][(*data).ncov].dim1>0){
         for(i=0;i<(*newt).nknots;i++){
            if((*newt).sub[i0][j0].kts1[i+1][0]>0){
               for(j=0;j<(*newt).sub[j0][(*data).ncov].dim1-1;j++){
                  if((*newt).sub[i0][j0].kts1[i+1][j+1]==0){
/* knot x knot */
                     criterion=testbasis(new,newt,criterion,data,i0,j0,i,j,
                                                                    (double)0);
                  }
               }
            }
/* knot x linear */
            else{
               criterion=testbasis(new,newt,criterion,data,i0,j0,i,-1,
                                                                    (double)0);
            }
         }
      }
      if(i0!=(*data).ncov){
/* linear x linear */
         if((*newt).sub[i0][j0].dim1==0){
            if((*newt).sub[i0][(*data).ncov].dim1>0 &&
               (*newt).sub[j0][(*data).ncov].dim1>0 ){
               criterion=
                  testbasis(new,newt,criterion,data,i0,j0,-1,-1,(double)0);
            }
         }
         else{
            for(i=0;i<(*newt).sub[i0][(*data).ncov].dim1-1;i++){
               if((*newt).sub[i0][j0].kts1[i+1][0]>0){
                  for(j=0;j<(*newt).sub[j0][(*data).ncov].dim1-1;j++){
                     if((*newt).sub[i0][j0].kts1[i+1][j+1]==0 &&
                        (*newt).sub[i0][j0].kts1[0][j+1]>0){
/* knot x knot */
                        criterion=testbasis(new,newt,criterion,data,i0,j0,i,j,
                                                                    (double)0);
                     }
                  }
               }
               else{
/* knot x linear */
                  criterion=testbasis(new,newt,criterion,data,i0,j0,i,-1,
                                                                    (double)0);
               }
            }
            for(j=0;j<(*newt).sub[j0][(*data).ncov].dim1-1;j++){
               if((*newt).sub[i0][j0].kts1[0][j+1]==0){
/* linear x knot */
                  criterion=testbasis(new,newt,criterion,data,i0,j0,-1,j,
                                                                    (double)0);
               }
            }
         }
      }
   }
   return criterion;
}
/******************************************************************************/

/* if a new knot is to be added in a one-covariate dimension or in time, we 
   have to search, and that is what we do in this routine */

static double search(new,newt,data,i0,mind)
struct space *newt,*new;
struct datastruct *data;
int i0,mind;

/* new   - the best added space up to now
   newt  - a space to which we can add
   data  - data
   i0    - first coordinate of the subdimension (second is data.ncov)
   mind  - minimum distance (in order statistics) between knots */

{
   double *sorted,*sorted2,critnew,crit,crit2,*kts;
   int i,lgth,iloc,lloc,bloc,uloc,iloc2,ll,uu,nx;

/* sorted  - sorted data or covariate
   sorted2 - uncensored data
   critnew - new criterion
   crit    - best criterion up to now
   crit2   - alternate new criterion
   testbasis - compute criterion for a basis
   i       - counter
   kts     - already used knots
   lgth    - number of already used knots
   iloc    - present location under study
   lloc    - lower bound to best location
   bloc    - best location up to now
   uloc    - upper bound to best location
   iloc2   - other location under study
   sort    - sorting routine 
   ll      - candidate for lloc
   uu      - candidate for uloc 
   nx      - (*data).ndata 
   find..  - find location for new knot under various circumstances */
   
/* initialization */
   bloc=-1;
   crit=-pow((double)10.,(double)20.);
   sorted=searchsorted;
   sorted2=searchsorted2;

/* find lgth, create kts: already used knots */
   if(i0==(*data).ncov){
      lgth=(*newt).nknots;
      kts=searchkts;
      for(i=0;i<lgth;i++) kts[i]=(*newt).knots[i];
      nx=0;
      for(i=0;i<(*data).ndata;i++){
         if((*data).delta[i]==1){
            sorted2[nx]=(*data).times[i];
            nx++;
         }
      }
      sort(sorted,sorted2,nx);
   }
   else{
      lgth=(*newt).sub[i0][(*data).ncov].dim1-1;
      kts=searchkts;
      for(i=0;i<lgth;i++) kts[i]=(*newt).sub[i0][(*data).ncov].ktsc[i];
      nx = (*data).ndata;
      sort(sorted,(*data).cov[i0],nx);
   }

/* find the interval */
   for(i=-2;i<=lgth;i++){
      if(lgth>0 && i==-2)i=0;
/* before first knot */
      if(i==0   && lgth>0) iloc=gindl(&ll,&uu,mind,sorted,nx,kts[0]);
/* after last knot */
      if(i==lgth&& lgth>0) iloc=gindr(&ll,&uu,mind,sorted,nx,kts[lgth-1]);
/* first knot */
      if(i==0   && lgth==0)iloc=gindx(&ll,&uu,nx,0);
      if(i==-1  && lgth==0)iloc=gindx(&ll,&uu,nx,1);
      if(i==-2  && lgth==0)iloc=gindx(&ll,&uu,nx,2);
/* in between knots */
      if(i>0    && i<lgth) iloc=gindm(&ll,&uu,mind,sorted,nx,kts[i-1],kts[i]);
/* possible location */
      if(iloc>=0){
         critnew=testbasis(new,newt,crit,data,i0,(*data).ncov,0,0,sorted[iloc]);
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

/* as long as the locations are different, do interval halving */
   do{
      if(sorted[uloc]>sorted[lloc]){
         iloc2=gindyr(uloc,bloc,sorted);
/* two search points, the upper one */
         if(iloc2>=0){
            crit2=testbasis(new,newt,crit,data,i0,(*data).ncov,0,
                                                               0,sorted[iloc2]);
         }
         else crit2=crit;

/* two search points, the lower one */
         iloc=gindyl(bloc,lloc,sorted);
         if(iloc>=0){
            critnew=testbasis(new,newt,crit2,data,i0,(*data).ncov,0,
                                                                0,sorted[iloc]);
         }
         else critnew=crit;
/* the middle one is the best, we call it quits */
         if(crit>=critnew && crit>=crit2){
            lloc=uloc;
         }
         else{
/* the lower search point is the best */
            if(critnew>crit2){
               uloc=bloc;
               bloc=iloc;
               crit=critnew;
            }
            else{
/* the upper search point is the best */
               lloc=bloc;
               bloc=iloc2;
               crit=crit2;
            }
         }
      }
   }while(sorted[uloc]>sorted[lloc]);

/* clean up and leave */
   return crit;
}
/******************************************************************************/

/* after another routine has decided to check the rao-criterion for a model
   with an added basis, this routine first adds the basis (addbasis), then
   it checks the criterion (critswap) - there are lots of possibilities to
   check. */

static double testbasis(new,newt,criterion,data,i0,j0,ki,kj,ti)
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
   int   position;

/* arg       - arguments for fct1 and fct2 in addbasis
   critswap  - computes rao and if there is improvement swaps the space
   addbasis  - adds a basis function
   tswapout  - removes a t-knot from a space
   tswapin   - adds a t-knot to a space
   position  - which position has the new t-knot */

/* most common occurences - preset for linear in covariates */
   arg[0]=-1.; 
   arg[1]=-1.; 
   arg[2]=-1.; 
   arg[3]=-1.;
/* for 1-d subspaces this is */

/* if we are adding a t-knot everything is different */
   if(j0==(*data).ncov && i0==(*data).ncov){ 
/* which is the knot */
      (*newt).knots[(*newt).nknots]=ti;
/* add the basisfunction */
      arg[0]=ti;
      arg[2]=(*newt).nknots;
      addbasis(i0,j0,arg,data,&((*newt).basis[(*newt).ndim]));
/* change the ordering of the space */
      position=tswapin(newt,(*data).ncov);
/* compute rao */
      criterion=critswap(newt,data,new,criterion,i0,j0,0);
/* change the ordering of the space back */
      tswapout(newt,(*data).ncov,position);
/* done */
      return criterion;
   } 

/* 1 covariate subdimension */
   if(j0==(*data).ncov && i0<(*data).ncov){
/* this is not the first (i.e. linear) space */
      if((*newt).sub[i0][j0].dim1>0){
/* what is the knot to be added */
         arg[0]=ti;  
         arg[2]=(*newt).sub[i0][(*data).ncov].dim1-1;
        (*newt).sub[i0][(*data).ncov].ktsc[(*newt).sub[i0][(*data).ncov].dim1-1]
                                                                            =ti;
      }
   }

/* a crossproduct subdimension */
   if(j0<(*data).ncov){
/* if it is the first one, it is linear*time-knot[ki] */
      if(ki>=0){
         arg[2]=ki;
         if(i0<(*data).ncov)
            arg[0]=(*newt).sub[j0][(*data).ncov].ktsc[(int)arg[2]];
         if(i0==(*data).ncov) arg[0]=(*newt).knots[(int)arg[2]];
      }
      if(kj>=0){
         arg[3]=kj;
         arg[1]=(*newt).sub[j0][(*data).ncov].ktsc[(int)arg[3]];
      }
      (*newt).sub[i0][j0].kts1[ki+1][kj+1]=1;
   }

   addbasis(i0,j0,arg,data,&((*newt).basis[(*newt).ndim]));
/* compute rao. possibly swap */
   criterion=critswap(newt,data,new,criterion,i0,j0,1);
   if(j0<(*data).ncov)(*newt).sub[i0][j0].kts1[ki+1][kj+1]=0;

/* done */
   return criterion;
}
/******************************************************************************/

/* After a new time-knot is introduced, we need to reorganize the indices. That
   is, all other time knots that are larger than the new knot are shifted one
   position. The output argument position is the position of the new basis
   function */

static int tswapin(spc,ncov)
struct space *spc;
int ncov;

/* spc  - the space that has to be reorganized 
   ncov - number of covariates */
{
   int position,i,j,k;

/* i,j      - counter
   position - position of the new knot 
   dsort    - sort a double array */

/* if there is only 1 knot, don't bother */
   if((*spc).nknots==0) return -1;

/* if the knot is in the right position, don't bother either */
   if((*spc).knots[(*spc).nknots]>(*spc).knots[(*spc).nknots-1]) return -1;

/* first find the correct position */
   position=0;
   for(i=0;i<(*spc).nknots;i++){
      if((*spc).knots[(*spc).nknots]>(*spc).knots[i]) position=i+1;
   }

/* the last basis function refers to the best position */
   (*spc).basis[(*spc).ndim].iknots=position;
   (*spc).basis[(*spc).ndim].t1=position;

/* check for all other basis functions whether the reference changes */
   for(i=0;i<(*spc).ndim;i++){
      if((*spc).basis[i].iknots>=position){
         ((*spc).basis[i].iknots)++;
         if((*spc).basis[i].b1==ncov) (*spc).basis[i].t1=(*spc).basis[i].iknots;
         else (*spc).basis[i].t2=(*spc).basis[i].iknots;
      }
   }

/* check up the subdimensions that are time * covariate */
   for(i=0;i<ncov;i++){
      for(j=0;j<=(*spc).sub[i][ncov].dim1;j++){
         for(k=(*spc).nknots;k>position;k--){
            (*spc).sub[ncov][i].kts1[k+1][j]=(*spc).sub[ncov][i].kts1[k][j];
         }
         (*spc).sub[ncov][i].kts1[position+1][j]=0;
      }
   }

/* finally, reorder the knots, just use sort */
   dsort((*spc).knots,(*spc).nknots+1);

   return position;
}
/******************************************************************************/

/* removes the time knot on position from a space */
static void tswapout(spc,ncov,position)

struct space *spc;
int ncov,position;

/* spc      - the space to remove a time-knot from
   ncov     - number of covariates
   position - position of the time knot to remove */
{
   int i,j,k;

/* if the position is -1 there is nothing to worry about */
   if(position== -1) return;

/* check out the basisfunctions */
   for(i=0;i<(*spc).ndim;i++){
      if((*spc).basis[i].iknots>position){
         ((*spc).basis[i].iknots)-=1;
         if((*spc).basis[i].b1==ncov) (*spc).basis[i].t1=(*spc).basis[i].iknots;
         else (*spc).basis[i].t2=(*spc).basis[i].iknots;
      }
   }

/* check up the subdimensions that are time * covariate */
   for(i=0;i<ncov;i++){
      for(j=0;j<=(*spc).sub[i][ncov].dim1;j++){
         for(k=position;k<(*spc).nknots;k++){
            (*spc).sub[ncov][i].kts1[k+1][j]=(*spc).sub[ncov][i].kts1[k+2][j];
         }
         (*spc).sub[ncov][i].kts1[(*spc).nknots+1][j]=0;
      }
   }

/* finally, reorder the knots */
   for(i=position;i<(*spc).nknots;i++) (*spc).knots[i]=(*spc).knots[i+1];
}
/******************************************************************************/

/* after another routine has decided to add a basis function, this routine
   actually adds the basis function */

static void addbasis(i0,j0,arg,data,basis)
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
   fct1,fct2 - functions to compute version */

/* update the scalar part of the basis function */
   (*basis).beta=0.;
   if(i0==(*data).ncov)(*basis).iknots=arg[2];
   else (*basis).iknots=-1;
   (*basis).b1=i0;
   (*basis).b2=j0;
   (*basis).t1=arg[2];
   (*basis).t2=arg[3];

/* values */
   if(i0!=(*data).ncov){

/* depends on 2 covariates */
      if(j0!=(*data).ncov){
         for(i=0;i<(*data).ndata;i++){
            if((*data).same[i]==1) (*basis).values[i]=(*basis).values[i-1];
            else (*basis).values[i]=
                           fct2((*data).cov[i0][i],(*data).cov[j0][i],arg);
            (*basis).values2[i]=(*basis).values[i];
         }
      }

/* depends on 1 covariate */
      else{
        for(i=0;i<(*data).ndata;i++){
            if((*data).same[i]==1) (*basis).values[i]=(*basis).values[i-1];
            else (*basis).values[i]=fct1((*data).cov[i0][i],arg,1);
            (*basis).values2[i]=(*basis).values[i];
         }
      }
   }
   else{

/* depends on 1 covariate and time */
      if(j0!=(*data).ncov){
         for(i=0;i<(*data).ndata;i++){
            if((*data).same[i]==1) (*basis).values[i]=(*basis).values[i-1];
            else (*basis).values[i]=fct1((*data).cov[j0][i],arg,2);
            (*basis).values2[i]=(*basis).values[i]*fct1((*data).times[i],arg,0);
         }
      }

/* depends on time only */
      else{
         for(i=0;i<(*data).ndata;i++){
            (*basis).values[i]=1.;
            (*basis).values2[i]=fct1((*data).times[i],arg,0);
         }
      }
   }
}
/******************************************************************************/

/* This function is used to compute the version component of a basis-function
   of one covariate (plus possibly time). If t[2(or3)]<0, the function is linear
   in x while if t[2(or 3)]>0 the function is proportional to (x-t[0(or 1)])+.*/

static double fct1(x,t,i)
double x,*t;
int i;

/* t - see description above
   x - value in dimension 1 
   i - elements 0 and 2 of t (i=1) or 1 and 3 (i=2) */

{
   if(i==0){
      if(t[2]<0)return (double) 1.;
      if(x>=t[0])return (double)0.;
      return t[0]-x;
   }
   if(i==1){
      if(t[2]<0)return x;
      if(x<=t[0])return (double)0.;
      return x-t[0];
   }
   else{
      if(t[3]<0)return x;
      if(x<=t[1])return (double)0.;
      return x-t[1];
   }
}
/******************************************************************************/

/* this function is used to compute the version component of a basis-function
   of two covariates. If t[2]<0/t[3]<0, the function is linear in x/y while
   if t[2]>0/t[3]>=0 the function is proportional to (x-t[1])+/(y-t[1])+. */

static double fct2(x,y,t)
double x,y,*t;

/* t - see description above
   x - value in dimension 1
   y - value in dimension 2  */

{
   if(t[2]> -0.5){
      x=x-t[0];
      if(x<0.)return (double)0.;
   }
   if((int)t[3]> -0.5){
      y=y-t[1];
      if(y<0.)return (double)0.;
   }
   return x*y;
}
/******************************************************************************/

/* after the space newt has been updated with extra new basisfunctions, it
   computes the criterion. If this is an improvement it copies the basis into
   new, then it restores newt. */

static double critswap(newt,data,new,criterion,i0,j0,ij)
double criterion;
struct datastruct *data;
int i0,j0,ij;
struct space *new,*newt;

/* criterion - best rao p-value up to now
   data      - the data
   i0,j0     - which subdimension is being altered
   ij        - which of the dimension attributes are to be altered
   new       - best model with additions
   newt      - model tested whether it is better */

{
   double crit,r1;
   int i,j;

/* crit      - present value of criterion
   grao      - computes rao-criterion
   swapspace - copies one space into another */

   for(i=0;i<(*newt).ndim;i++){
      (*newt).xtx[i][(*newt).ndim]=0.;
      for(j=0;j<(*data).ndata;j++){
         (*newt).xtx[i][(*newt).ndim]+=(*newt).basis[i].values2[j]*
                                       (*newt).basis[(*newt).ndim].values2[j];
      }
      (*newt).xtx[(*newt).ndim][i]=(*newt).xtx[i][(*newt).ndim];
   }
   (*newt).xtx[(*newt).ndim][(*newt).ndim]=0.;
   for(j=0;j<(*data).ndata;j++){
      (*newt).xtx[(*newt).ndim][(*newt).ndim]+=
                                  (*newt).basis[(*newt).ndim].values2[j]*
                                  (*newt).basis[(*newt).ndim].values2[j];
   }
   r1=condition((*newt).xtx,(*newt).ndim+1);
   if(r1<pow((double)10.0,(double)-13.))return criterion;

/* update the dimension parameters */

   ((*newt).ndim)+=1;
   if(ij==0)((*newt).nknots)+=1;
   else ((*newt).sub[i0][j0].dim1)+=1;

   crit=grao(newt,data);

/* if there is improvement, copy the space */
   if(crit>criterion){
      swapgspace(new,newt,(*data).ndata,(*data).ncov);
      criterion=crit;
   }

/* change back the dimensions */

   ((*newt).ndim)-=1;
   if(ij==0)((*newt).nknots)-=1;
   else ((*newt).sub[i0][j0].dim1)-=1;

   return criterion;  
}
/******************************************************************************/

/* this routine computes the extra elements of hessian and score then it
   computes rao - the routine is very much like the routine complog, which
   is part of Newton, except that it does not compute the log-likelihood 
   and it makes use of the fact that part of b0, b1 and b2 might be known and
   completely at the end, it computes the rao statistic.                      */

static double grao(spc,data)
struct space *spc;
struct datastruct *data;

/* spc   - the present model 
   data  - the data */
{
   double l0,raoc=0.,*scorecopy,r[3],rr[4],**hhh,*ss;
   int i,j,l,where,whereold=0,naaap,extra,again1,again2;

/* extra   - (*spc).ndim-1
   l0      - lower integration bound
   raoc    - rao-score statistic
   veint   - computes integral (see newton for this routine)
   upbasis - updates basis for one point  (see newton for this routine)
   i,j,l  - counter
   r       -integrals
   where   - in which interval is the datapoint 
   whereold- where was the previous case
   naaap   - is this case similar to the previous one
   hhh     - old, still useful, hessian integrals
   ss      - old, still useful, score integrals 
   again1  - how many follow up cases are exactly the same 
   again2  - how many of these have delta=1 */

/* initializations */
   extra=(*spc).ndim-1;
   (*spc).score[extra]=0.;
   for(j=0;j<(*spc).ndim;j++){
      (*spc).info[extra][j]=0.;
      (*spc).info[j][extra]=0.;
   }

/* allocation */
   hhh=raohhh;
   ss=raoss;

/* now circle the datapoints */
   for(i=0;i<(*data).ndata;i++){

/* get the again things */
      again1=1;
      again2=(*data).delta[i];
      for(j=i+1;j<(*data).ndata;j++){
         if((*data).same[j]==1 && (*data).times[i]==(*data).times[j]){
            again1++;
            if((*data).delta[j]==1)again2++;
         }
         else{
            j=(*data).ndata;
         }
      }
      if((*spc).basis[extra].iknots>-1){
         l0=(*spc).knots[(*spc).basis[extra].iknots];
         if((*data).times[i]>l0){
            for(j=i+again1;j<(*data).ndata;j++){
               if((*data).same[j]==1 && (*data).times[j]> l0)again1++;
               else j=(*data).ndata;
            }
         }
      } 

/* in which interval is the datapoint ? */
      where=(*spc).nknots;
      for(j=0;j<(*spc).nknots;j++){
         if((*spc).knots[j]>(*data).times[i]){
            where=j;
            j=(*spc).nknots;
         }
      }
      naaap=0;
      if((*data).same[i]==1 && whereold==where) naaap=1;
      whereold=where;

/* initialize basis */
      if(naaap==0){
         for(j=0;j<=(*spc).nknots;j++){
            (*spc).b0[j][0]=0.;
            (*spc).b1[j][0]=0.;
         }
      }
      (*spc).b2[0]=0.;

/* update spc.b2,spc.b0,spc.b1 per basisfunction */
      for(j=0;j<(*spc).ndim;j++){
         if(naaap==0){
            upbasis((*spc).knots,(*spc).nknots,(*spc).b0,(*spc).b1,(*spc).b2,i,
                     j+1,&((*spc).basis[j]),where,0);
         }
         else{
            (*spc).b2[j+1]=(*spc).basis[j].values2[i];
            (*spc).b2[0]+=(*spc).b2[j+1]*(*spc).basis[j].beta;
         }
      }

/* add the delta terms to spc.score */
      (*spc).score[extra]+=(*spc).b2[extra+1]*again2;

/* the numerical integrals, all quite straight forward */
      if(naaap==0){
         ss[extra]=0.;
         for(l=0;l<(*spc).ndim;l++) hhh[extra][l]=0.;
         for(j=0;j<where;j++){

/* the lower bound */
            if(j==0) l0=0.;
            else l0=(*spc).knots[j-1];

            veint(r,(*spc).b1[j][0],(*spc).b0[j][0],l0,(*spc).knots[j]);
            rr[0]=r[0]*(*spc).b0[j][extra+1];
            rr[1]=r[1]*(*spc).b0[j][extra+1];
            rr[2]=r[1]*(*spc).b1[j][extra+1]+rr[0];
            rr[3]=r[2]*(*spc).b1[j][extra+1]+rr[1];

/* score */
            ss[extra]+=rr[2];
            for(l=0;l<(*spc).ndim;l++){

/* hessian */
               hhh[extra][l]+=rr[3]*(*spc).b1[j][l+1]+rr[2]*(*spc).b0[j][l+1];
            }
         }
      }
      j=where;
      if(j==0) l0=0.;
      else l0=(*spc).knots[j-1];
      veint(r,(*spc).b1[j][0],(*spc).b0[j][0],l0,(*data).times[i]);
      rr[0]=r[0]*(*spc).b0[j][extra+1];
      rr[1]=r[1]*(*spc).b0[j][extra+1];
      rr[2]=r[1]*(*spc).b1[j][extra+1]+rr[0];
      rr[3]=r[2]*(*spc).b1[j][extra+1]+rr[1];
      if(again1==1){
/* score */
         (*spc).score[extra]-=ss[extra]+rr[2];
/* hessian */
         for(l=0;l<(*spc).ndim;l++){
            (*spc).info[extra][l]-=hhh[extra][l]+
                               rr[3]*(*spc).b1[j][l+1]+rr[2]*(*spc).b0[j][l+1];
         }
      }
      else{
/* score */
         (*spc).score[extra]-=again1*(ss[extra]+rr[2]);
/* hessian */
         for(l=0;l<(*spc).ndim;l++){
            (*spc).info[extra][l]-=again1*(hhh[extra][l]+
                               rr[3]*(*spc).b1[j][l+1]+rr[2]*(*spc).b0[j][l+1]);
         }
         i+=again1-1;
      }
   }

/* symmaterize the hessian */
   for(l=0;l<extra;l++) (*spc).info[l][extra]=(*spc).info[extra][l];

/* copy hessian and score, since we're going to destroy it */
   scorecopy=raoscorecopy;
   for(l=0;l<(*spc).ndim;l++) scorecopy[l]=(*spc).score[l];

/* compute rao */
   l=glusolve2((*spc).info,(*spc).ndim,scorecopy);
   if(l>0) for(l=0;l<(*spc).ndim;l++) raoc+=scorecopy[l]*(*spc).score[l];
   else raoc=0.;

   return -raoc;
}
static void getse(spc)
struct space *spc;
{
   int i,j;
   double **xx;
   xx=getsexx;
   for(i=0;i<(*spc).ndim;i++){
      for(j=0;j<(*spc).ndim;j++)xx[i][j]=(*spc).info[i][j];
   }
   gluinverse(xx,(*spc).ndim);
   for(i=0;i<(*spc).ndim;i++)(*spc).basis[i].se=sqrt(-xx[i][i]);
}
/******************************************************************************/
/* finds a new location in an interval (l,b) - that is the lower end might not
   have been tested yet */
static int gindyl(u,l,x)
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
/* finds a new location in an interval (b,u) - that is the upper end might not
   have been tested yet */
static int gindyr(u,l,x)
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
/* Finds a possible location for a knot on the interval (0,knot1)
   ll - lowest number we can search on in the future
   uu - highest number we can search on in the future
   mind minimum distance between knots
   x  - data
   nx - length of data
   knt- knot */

static int gindl(ll,uu,mind,x,nx,knt)
double *x,knt;
int nx,*ll,*uu,mind;
{

/* i  - utility
   zlocation - finds uu */

   int i;

   (*uu)=zlocation(0,x,nx,knt);
   if((*uu)<mind)return -1;
   i=((*uu)-1)/2;
   if((*uu)-i<mind+1)i=(*uu)-mind-1;
   *ll=0;
   *uu=(*uu)-mind-1;
   return i;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (knot-last,nx-1)
   ll - lowest number we can search on in the future
   uu - highest number we can search on in the future
   mind minimum distance between knots
   x  - data
   nx - length of data
   knt- knot */

static int gindr(ll,uu,mind,x,nx,knt)
double *x,knt;
int nx,*ll,*uu,mind;
{

/* i  - utility
   zlocation - finds ll */

   int i;

   (*ll)=zlocation(1,x,nx,knt);
   if(nx-1-(*ll)<mind)return -1;
   i=(nx+(*ll))/2;
   if(i-(*ll)<mind+1)i=(*ll)+mind+1;
   *uu=nx-1;
   *ll=(*ll)+mind+1;
   return i;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (0,nx-1)
   ll - lowest number we can search on in the future
   uu - highest number we can search on in the future
   nx - length of data */

static int gindx(ll,uu,nx,i)
int nx,*ll,*uu,i;
{
   if(i==0){
      *ll=0;
      *uu=nx/2;
      return nx/4;
   }
   if(i==1){
      *ll=nx/4;
      *uu=(3*nx)/4;
      return nx/2;
   }
      *ll=nx/2;
      *uu=nx-1;
      return (3*nx)/4;
}
/******************************************************************************/
/* Finds a possible location for a knot on the interval (k0,k1)
   ll - lowest number we can search on in the future
   uu - highest number we can search on in the future
   mind minimum distance between knots
   x  - data
   nx - length of data
   k0 - knot
   k1 - knot */

static int gindm(ll,uu,mind,x,nx,k0,k1)
double *x,k0,k1;
int nx,*ll,*uu,mind;
{
/* zlocation - finds ll */


   (*ll)=zlocation(1,x,nx,k0);
   (*uu)=zlocation(0,x,nx,k1);
   if((*uu)-(*ll)<2*mind+1)return -1;
   *uu=(*uu)-mind-1;
   *ll=(*ll)+mind+1;
   return ((*uu)+(*ll))/2;
}
/******************************************************************************/
/* finds the lowest (if what = 0) or the highest (if what = 1) index of x for
   which x==k  */

/* what - see above
   x    - data
   nx   - length data
   k    - see above */

static int zlocation(what,x,nx,k)
int nx,what;
double k,*x;
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

/* this file contains the Newton-Raphson iteration loop subroutine, and the
   functions which it needs */

/******************************************************************************/

/* this routine does the Newton-Raphson iteration loop */

static double newton(spc,data,precision,silent,oops)
int *oops,silent,precision;
struct datastruct *data; 
struct space *spc;

/* spc   - the present model
   data  - the data 
   precision - how precise (0-not, 1-very) 
   silent- should diagnostic output be printed */

{
   double zerror=0.0001,lnew,logl,*scp;
   int iter,ihalf,maxiter=100,j,*where;

/* zerror  - convergence criterion
   lnew   - new loglikelihood
   compall- routine to compute hessian, score and likelihood
   complog- routine to compute likelihood
   logl   - log-likelihood 
   scp    - copy of the score vector
   iter   - iteration counter
   ihalf  - how many times has there been step-halving
   maxiter- maximum number of iterations 
   j      - counter
   where  - where are the points 
   oops   - there are problems if oops is 1 on exit */

/* allocation/initialization */
   (*oops)=0; 
   scp=newtonscp;
   where=newtonwhere;
   zerror=.01;
/* if(precision==0) zerror=.01; */
   if(precision==2){
      ihalf=0;
      logl=0.;
      for(j=0;j<(*data).ndata;j++){
         ihalf+=(*data).delta[j];
         logl+=(*data).times[j];
      }
      (*spc).basis[0].beta=-log(logl/(double)ihalf);
      for(j=1;j<(*spc).ndim;j++)(*spc).basis[j].beta=0.;
   }

/* start iteration */
   for(iter=0;iter<maxiter;iter++){

/* compute score, info and old loglikelihood */
      logl=compall(spc,(*data).times,(*data).delta,(*data).ndata,
                                                      iter,where,(*data).same);
/*    (void)Rprintf("&& %f\n",logl);(void)fflush(stdout); */

/* copy score for possible later use */
      for(j=0;j<(*spc).ndim;j++) scp[j]=(*spc).score[j];

/* solve system */
      j=glusolve2((*spc).info,(*spc).ndim,(*spc).score);
      if(j==0){
         if(precision==1){
            (*oops)=1;
            return 0.;
         }
         nrerror("instable system during NR iterations");
      }
      ihalf=1; 
      do{

/* compute new loglikelihood */
         lnew=complog(spc,(*data).times,(*data).delta,
                                          (*data).ndata,0,where,(*data).same);
/*       (void)Rprintf("^^ %f\n",lnew);(void)fflush(stdout); */

/* step halving if required */
         if(lnew<logl-zerror){
     /*     (void)Rprintf("%20.14f %20.14f \n",lnew,logl); (void)fflush(stdout);  */

/* oops, too much step halving */
            if(ihalf>2048 || (ihalf >256 && precision==1)){
               if(precision==1){
                  (*oops)=1;
                  return 0.;
               }
               nrerror("too much step-halving");
            }

/* the actual halving */
            ihalf=ihalf*2;
            for(j=0;j<(*spc).ndim;j++) (*spc).score[j]=(*spc).score[j]/2;
         }
      }while(lnew<logl-zerror);

/* record the new solution */
      for(j=0;j<(*spc).ndim;j++){

/* add the step */
         (*spc).basis[j].beta-=(*spc).score[j];
      }

/* did we converge */
      if(precision==1 && humbertester(logl)+humbertester(lnew)!=6){
         (*oops)=1;
         return 0.;
      }
      if(lnew-logl<zerror) iter=maxiter+1000;
   }

/* did we finish because we converged? */
   if(iter<maxiter+500){
      nrerror("no convergence");
   }

/* get a fresh copy of score */
   if(silent!=1)Rprintf("|| logl= %.2f (nd=%d)\n",lnew,(*spc).ndim);
   (void)fflush(stdout);
   lnew=compall(spc,(*data).times,(*data).delta,(*data).ndata,iter,
                                                            where,(*data).same);

   return lnew;
}

/******************************************************************************/

/* this routine computes hessian score and log-likelihood */

static double compall(spc,data,delta,ndata,iter,where,same)
struct space *spc;
int ndata,*delta,iter,*where,*same;
double *data;

/* spc   - the present model
   data  - the data (only the times)
   delta - delta
   iter  - which iteration is this
   ndata - number of datapoints 
   where - in which interval is the datapoint */

{
   double r[3],l0,logl,*ss,**hhh,lala,rr[5];
   int i,j,k,l,naaap,whereold=0;

/* r      - intgrals
   l0     - lower integration bound
   logl   - loglikelihood
   ss     - integrals from the previous time for score.... still useful???
   hhh    - integrals from the previous time for hessian.... still useful???
   lala   - integrals from the previous time for logl.... still useful???
   upbasis- updates spc.b0, spc.b1 and spc.b2 for one point 
   veint  - computes integrals
   i,j,k,l- counters 
   whereold - where was the previous point
   naaap  - can we use some of the old stuff?   */

/* initializations  and allocations */
   ss=compallss;
   hhh=compallhhh;
   logl=0.;
   for(i=0;i<(*spc).ndim;i++){
      (*spc).score[i]=0.;
      for(j=0;j<(*spc).ndim;j++) (*spc).info[i][j]=0.;
   }

/* now circle the datapoints */
   for(i=0;i<ndata;i++){

/* in which interval is the datapoint ? */
      if(iter==0){
         where[i]=(*spc).nknots;
         for(j=0;j<(*spc).nknots;j++){
            if((*spc).knots[j]>data[i]){
               where[i]=j;
               j=(*spc).nknots;
            }
         }
      }

/* is this the same interval as the previous point (with the same covariates) */
      naaap=0;
      if(same[i]==1 && whereold==where[i]) naaap=1;
      whereold=where[i];

/* initialize basis */
      if(naaap==0){
         for(j=0;j<=(*spc).nknots;j++){
            (*spc).b0[j][0]=0.;
            (*spc).b1[j][0]=0.;
         }
      }
      (*spc).b2[0]=0.;
      
/* update spc.b2, spc.b0 and spc.b1  per basisfunction */
      for(j=0;j<(*spc).ndim;j++){
         if(naaap==0){
            upbasis((*spc).knots,(*spc).nknots,(*spc).b0,(*spc).b1,(*spc).b2,i,
                                    j+1,&((*spc).basis[j]),where[i],0);
         }
         else {
            (*spc).b2[j+1]=(*spc).basis[j].values2[i];
            (*spc).b2[0]+=(*spc).b2[j+1]*(*spc).basis[j].beta;
         }
      }

/* add the delta terms to spc.score and loglikelihood */
      if(delta[i]==1){
         logl+=(*spc).b2[0];
         for(j=0;j<(*spc).ndim;j++) (*spc).score[j]+=(*spc).b2[j+1];
      }

/* the numerical integrals, all quite straight forward */
      if(naaap==0){
         lala=0.;
         for(k=0;k<(*spc).ndim;k++){
            ss[k]=0.;
            for(l=0;l<=k;l++) hhh[k][l]=0.;
         }

/* for all the intervals between knots */
         for(j=0;j<where[i];j++){
/* the lower bound */
            if(j==0) l0=0.;
            else l0=(*spc).knots[j-1];

            veint(r,(*spc).b1[j][0],(*spc).b0[j][0],l0,(*spc).knots[j]);
            lala+=r[0];
            for(k=0;k<(*spc).ndim;k++){

/* score */
               rr[0]=r[0]*(*spc).b0[j][k+1];
               rr[1]=r[1]*(*spc).b0[j][k+1];
               rr[2]=r[1]*(*spc).b1[j][k+1]+rr[0];
               rr[3]=r[2]*(*spc).b1[j][k+1]+rr[1];
               ss[k]+=rr[2];
               for(l=0;l<=k;l++){

/* hessian */
                 hhh[k][l]+=rr[3]*(*spc).b1[j][l+1]+rr[2]*(*spc).b0[j][l+1];
               }
            }
         }
      }
/* for the interval between the last knot and the observation */
      j=where[i];
/* the lower bound */
      if(j==0) l0=0.;
      else l0=(*spc).knots[j-1];

/* the integrals */
      veint(r,(*spc).b1[j][0],(*spc).b0[j][0],l0,data[i]);

/* the loglikelihood */
      logl-=r[0]+lala;

/* score */
      for(k=0;k<(*spc).ndim;k++){
         rr[0]=r[0]*(*spc).b0[j][k+1];
         rr[1]=r[1]*(*spc).b0[j][k+1];
         rr[2]=r[1]*(*spc).b1[j][k+1]+rr[0];
         rr[3]=r[2]*(*spc).b1[j][k+1]+rr[1];
         (*spc).score[k]-=ss[k]+rr[2];
   
/* hessian */
         for(l=0;l<=k;l++){
            (*spc).info[k][l]-=hhh[k][l]+
                               rr[3]*(*spc).b1[j][l+1]+rr[2]*(*spc).b0[j][l+1];
         }
      }
   }

/* symmaterize the hessian */
   for(k=0;k<(*spc).ndim-1;k++) {
      for(l=k+1;l<(*spc).ndim;l++) (*spc).info[k][l]=(*spc).info[l][k];
   }
   return logl;
}

/******************************************************************************/

/* this routine updates basis0, basis1 and basis2 (spc.b0, spc.b1 and spc.b2)
   for one basisfunction for one datapoint, for the hessian/score/logl case */

static void upbasis(knots,nknots,basis0,basis1,basis2,idt,ifc,basf,where,il)
int nknots,idt,ifc,where,il;
double **basis0,**basis1,*basis2,*knots;
struct basisfunct *basf;

/* knots  - time-knots
   nknots - number of time-knots
   basis0 - first element constant term of lambda[between-2-knots,datapoint]
            others constant term basisfunction(j)[between-2-knots,datapoint]
   basis1 - basis0, but linear term
   basis2 - first element, lambda in a datapoint others, basisfunction(j) in a
            datapoint
   idt    - number of datapoint on which we are working
   ifc    - basisfunction on which we are working
   basf   - one basisfunction
   where  - in between which knots is the basisfunction 
   il     - should the lambda part be done (not always if called by rao) */

{
   int j;
   double x;

/* j      - counter
   x      - to save some typing (below) */

   x=(*basf).values[idt];

/* initializations */

   for(j=0;j<=nknots;j++){
      basis0[j][ifc]=0;
      basis1[j][ifc]=0;
   }
   basis2[ifc]=(*basf).values2[idt];

/* the two issues that determine the formulas are (i) is a t-knots involved
   in this basisfunction and (ii) where is the datapoint relative to this knot*/

/* no t-knot, basisfunction is constant */
   if((*basf).iknots== -1) for(j=0;j<=where;j++) basis0[j][ifc]=x;

/* 1 t-knot, basisfunction is 0 after knot, and has slope -1 before knot */
   else{
      for(j=0;j<=where && j<=(*basf).iknots;j++){
         basis0[j][ifc]=knots[(*basf).iknots]*x;
         basis1[j][ifc]= -x;
      }
   }

/* for the lambda part multiply by beta */
   if(il==0){
      for(j=0;j<=where;j++){
         basis0[j][0]+=basis0[j][ifc]*(*basf).beta;
         basis1[j][0]+=basis1[j][ifc]*(*basf).beta;
      }
      basis2[0]+=basis2[ifc]*(*basf).beta;
   }
}

/******************************************************************************/

/* this function computes the log-likelihood (log-density) of a vector (one)
   datapoint */

static double complog(spc,data,delta,ndata,iwhere,vwhere,same)
struct space *spc;
int ndata,*delta,iwhere,*vwhere,*same;
double *data;

/* spc   - the present model
   data  - the data (only time)
   delta - delta
   ndata - number of datapoints - if <=0: do only for one point 
   vwhere- where is the data 
   iwhere- is vwhere alradey computed
   same  - for which observations are the covariates the same as for previous */

{
   double basis2,*basis0,*basis1,l0,logl,rall;
   int i,j,where,ik1,ik2,naaap,whereold=0;

/* basis2 - lambda in a datapoint 
   basis0 - element constant term of lambda[between-2-knots,datapoint]
   basis1 - basis0, but linear term
   l0     - lower integration bound
   logl   - loglikelihood
   eint   - computes integral
   rall   - log-likelihood result of the previous observation
   upbasis2- updates basis for one point 
   i,j    - counter
   where  - in which interval is the datapoint 
   ik1,ik2- lower and upper bound on range of datapoints
   naaap  - can we use some of the stuff from the previous observation?
   whereold - where was the previous observation */

/* initializations */
   logl=0.;
   basis0=complogbasis0;
   basis1=complogbasis1;

/* for which datapoints? */
   ik1=0;
   ik2=ndata;

/* if ndata<0 - only for 1 point */
   if(ndata<=0){
      ik1=-ndata;
      ik2=ik1+1;
   }

/* circle the data */
   if(ndata>0){
      for(j=0;j<(*spc).ndim;j++) (*spc).basis[j].beta-=(*spc).score[j];
   }
   for(i=ik1;i<ik2;i++){

/* in which interval is the data? */
      if(iwhere==1){
         where=(*spc).nknots;
         for(j=0;j<(*spc).nknots;j++){
            if((*spc).knots[j]>data[i]){
               where=j;
               j=(*spc).nknots;
            }
         }
      }
      else{
         where=vwhere[i];
      }

/* can we use old stuff? */
      naaap=0;
      if(same[i]==1 && ndata>0 && whereold==where) naaap=1;
      whereold=where;

/* initialize */
      if(naaap==0){
         for(j=0;j<=(*spc).nknots;j++){
            basis0[j]=0.;
            basis1[j]=0.;
         }
      }
      basis2=0.;

/* update basis0, basis1 and basis2 */
      for(j=0;j<(*spc).ndim;j++){
         if(naaap==0){
            upbasis2((*spc).knots,basis0,basis1,&basis2,i,&((*spc).basis[j]),
                                                                where);
         }
         else{
            basis2+=(*spc).basis[j].values2[i]*(*spc).basis[j].beta;
         }
      }

/* the delta part */
      if(delta[i]==1) logl+=basis2;
 
/* the integrals */
      if(naaap==0){
         rall=0;

/* per interval between knots */
         for(j=0;j<where;j++){

/* lower bound */
            if(j==0) l0=0.;
            else l0=(*spc).knots[j-1];

/* integrals */
            rall+=eint(basis1[j],basis0[j],l0,(*spc).knots[j]);
         }
      }
/* from the last knot to the datapoint */
      if(where==0) l0=0.;
      else l0=(*spc).knots[where-1];
      logl-=rall+eint(basis1[where],basis0[where],l0,data[i]);
   }

/* clean up */
   if(ndata>0){
      for(j=0;j<(*spc).ndim;j++) (*spc).basis[j].beta+=(*spc).score[j];
   }
   return logl;
}

/******************************************************************************/

/* comprable to upbasis, but does less, since only the loglikelihood will
   be computed */

static void upbasis2(knots,basis0,basis1,basis2,idt,basf,where)
int idt,where;
double *basis0,*basis1,*basis2,*knots;
struct basisfunct *basf;

/* knots  - time-knots
   basis2 - lambda in a datapoint
   basis0 - element constant term of lambda[between-2-knots,datapoint]
   basis1 - basis0, but linear term
   idt    - number of datapoint on which we are working
   basf   - one basisfunction
   where  - in between which knots is the basisfunction */

{
   int j;
   double x;
/* j    - counter
   x    - save typing (see below) */

   x=(*basf).values[idt]*(*basf).beta;
   *basis2+=(*basf).values2[idt]*(*basf).beta;

/* no t-knots, basisfunction is constant */
   if((*basf).iknots== -1){
      for(j=0;j<=where;j++) basis0[j]+=x;
   }

/* 1 t-knot, basisfunction is 0 after knot, and has slope -1 before knot */
   else{
      for(j=0;j<=where&&j<=(*basf).iknots;j++){
         basis0[j]+=knots[(*basf).iknots]*x;
         basis1[j]-=x;
      }
   }
}

/******************************************************************************/

/*           u
            /
            | (b1*x+b2)
            |e          dx
            |
            /
           l            */

static double eint(b1,b2,l,u)
double b1,b2,l,u;

/* just work it out */
{
   double c1,c2;
   if(b1!=0.){
      c1=b1*u+b2;
      c2=b1*l+b2;

/* take the numerically most stable exponents */

      if(c1*c2<= 0.){
        return (exp(c1)-exp(c2))/b1;
      }
      if(fabs(c1)>fabs(c2)){
         return (exp(c1-c2)-1.)*exp(c2)/b1;
      }
      return (1.-exp(c2-c1))*exp(c1)/b1;
   }
   return (u-l)*exp(b2);
}
/******************************************************************************/

/*           u
            /
            |     2           (b1*x+b2)
            |(a1*x +a2*x+a3)e          dx
            |
            /
           l            */

static void veint(r,b1,b2,l,u)
double r[3],b1,b2,l,u;

/* just work it out */
{
   double d21,d22,d12,d11,c1,c2,xx,e1,e2;
   if(b1!=0.){
      xx=2./b1;
      d21=(u*(u-xx)+xx/b1);
      d22=(l*(l-xx)+xx/b1);
      d11=(u-1./b1);
      d12=(l-1./b1);
      c1=b1*u+b2;
      c2=b1*l+b2;

/* take the numerically most stable exponents */

      if(c1*c2<= 0.){
        e1=exp(c1)/b1;
        e2=exp(c2)/b1;
        r[2]=(d21*e1-d22*e2);
        r[1]=(d11*e1-d12*e2);
        r[0]=(e1-e2);
        return;
      }
      if(fabs(c1)>fabs(c2)){
         e1=exp(c1-c2);
         e2=exp(c2)/b1;
         r[2]=(d21*e1-d22)*e2;
         r[1]=(d11*e1-d12)*e2;
         r[0]=(e1-1.)*e2;
         return;
      }
      e1=exp(c2-c1);
      e2=exp(c1)/b1;
      r[2]=(d21-d22*e1)*e2;
      r[1]=(d11-d12*e1)*e2;
      r[0]=(1.-e1)*e2;
      return;
   }
   e1=exp(b2);
   r[0]=(u-l)*e1;
   r[1]=(u*u-l*l)*e1/2.;
   r[2]=(u*u*u-l*l*l)*e1/3.;
   return;
}
/******************************************************************************/
/* this routine searches all dimensions for something to remove */

static void gremdim(spc,ncov,ndata,silent)
struct space *spc;
int ncov,ndata,silent;

/* spc  - the model from which to remove something
   ncov - number of covariates
   ndata- number of observations 
   silent- should diagnostic output be printed?  */

{
   int i,j,k,n,bb1,bt1,nb1,nb2,nt1,nt2,bb2,is,bw;
   double criterion,wald,**x,*y,**xtx,*xty;

/* i,j,k  - counter
   n      - can this basisfunction be removed? 1=yes, 0=no
   bb1/bt1/bb1/bw - the b1/t1/b2 components and ranknumber of a basisfunction
            of the basisfunction which is the candidate to be removed.
   nb1/nb2/nt1/nt2 - the b1/b2/t1/t2 components of the present basis function
            under examination
   criterion - best wald criterion up to now
   wald   - criterion for the dimension now to be removed
   basisswap - cuts out a basisfunction;
   cleanupt/cleanup1 - removes basisfunctions and knots from spc 
   y      - fitted values before removal
   x      - inner products of basisfunctions versus data
   xtx    - x'x
   xty    - x'y and after solving, the new starting values
   is     - should we compute new starting values */

/* initializations, allocation of storage */
   criterion=pow((double)10.,(double)100.);
   y=remdimy;
   is=1;
   for(i=0;i<ndata;i++) y[i]=0.;

/* per basisfunction */
   for(j=0;j<(*spc).ndim;j++){
      if(fabs((*spc).basis[j].beta)>10000)is=0;
      if(is==1){
         for(i=0;i<ndata;i++){
            y[i]+=(*spc).basis[j].values2[i]*(*spc).basis[j].beta;
         } 
      } 
   }

/* invert the hessian */
   gluinverse((*spc).info,(*spc).ndim);
   for(i=0;i<(*spc).ndim;i++) (*spc).basis[i].se=sqrt(-(*spc).info[i][i]);

/* circle all the basisfunctions except for the first one */
   for(i=1;i<(*spc).ndim;i++){
      nb1=(*spc).basis[i].b1;
      nt1=(*spc).basis[i].t1;
      nb2=(*spc).basis[i].b2;
      nt2=(*spc).basis[i].t2;
      n=1;

/* check whether we can remove the basisfunction */
      for(j=1;j<(*spc).ndim;j++) if(j!=i){
/* if this is a 1-d basisfunction, the knot cannot be used anywhere else */
         if(nb2==ncov){
            if((*spc).basis[j].b1==nb1 && (*spc).basis[j].t1==nt1) n=0;
            if(nb1!=ncov&&(*spc).basis[j].b2==nb1&&(*spc).basis[j].t2==nt1) n=0;
            if(nt1==-1 && (*spc).basis[j].b1==nb1)n=0;
            if(nt1==-1 && (*spc).basis[j].b2==nb1)n=0;
         }
/* if this is a 2-d basisfunction and it is linear, all nots should be gone */
         else{
            if((*spc).basis[j].b1==nb1 && (*spc).basis[j].b2==nb2) {
               if(nt2==-1 && nt1==-1) n=0;
               if(nt1==-1 && nt2>= 0 && (*spc).basis[j].t2==nt2) n=0;
               if(nt1>= 0 && nt2==-1 && (*spc).basis[j].t1==nt1) n=0;
            }
         }
         if(n==0) j=(*spc).ndim;
      }

/* if we can still remove, compute wald */
      if(n>0){
         wald=-(*spc).basis[i].beta*(*spc).basis[i].beta/(*spc).info[i][i];

/* record if we improve */
         if(wald<criterion){
            bw=i;
            bb1=nb1;
            bb2=nb2;
            bt1=nt1;
            criterion=wald;
         }
      }
   }

/* now remove the worst dimension */
   (*spc).ndim-=1;
   if(silent!=1)uuu(spc,bb1,bb2,bt1,(*spc).basis[bw].t2,ncov,1);
   if(silent!=1)(void)Rprintf("(wald = %.2f) ",wald);
   basisswap(&((*spc).basis[bw]),&((*spc).basis[(*spc).ndim]),ndata);
/* it is a 1d basisfunction */
   if(bb2==ncov){
/* if it is a t-only dimension */
      if(bb1==ncov) cleanupt(bt1,spc,ncov);
/* if it is a covariate-only dimension */
      else cleanup1(bb1,bt1,spc,ncov);
   }
   else{
/* if it is a two variable dimension */
      (*spc).sub[bb1][bb2].dim1-=1;
   }

   if(is==1){
/* get xtx */
   xtx=remdimxtx;
   for(j=0;j<(*spc).ndim;j++){
      for(i=0;i<=j;i++){
         xtx[i][j]=0.;
         for(k=0;k<ndata;k++) xtx[i][j]+=
             ((*spc).basis[i].values2[k])*((*spc).basis[j].values2[k]);
      }
   }
   for(j=0;j<(*spc).ndim;j++){
      for(i=j+1;i<(*spc).ndim;i++) xtx[i][j]=xtx[j][i];
   }

/* get xty */
   xty=remdimxty;
   for(i=0;i<(*spc).ndim;i++){
      xty[i]=0.;
      for(k=0;k<ndata;k++) xty[i]+=((*spc).basis[i].values2[k])*y[k];
   }

/* solve */
   glusolve(xtx,(*spc).ndim,xty);

/* record */
   for(i=0;i<(*spc).ndim;i++) (*spc).basis[i].beta=xty[i];
   
   }
   return;
}
/******************************************************************************/
/* this routine removes a time knot from the space */

static void cleanupt(ki,spc,ncov)
struct space *spc;
int ki,ncov;

/* spc   - space to be reduced
   ki    - ranknumber of knot
   ncov  - number of covariates */
{
   int j;
/* j   - counter  */

/* change the number of t-knots */
   (*spc).nknots-=1;

/* change all the other pointers */
   if((*spc).nknots>0){

/* first in the basisfunctions */
      for(j=0;j<(*spc).ndim;j++){
         if((*spc).basis[j].b1==ncov && (*spc).basis[j].t1>ki){
            (*spc).basis[j].t1-=1;
            (*spc).basis[j].iknots-=1;
         }
      }

/* in the knots themselves */
      for(j=ki;j<(*spc).nknots;j++) (*spc).knots[j]=(*spc).knots[j+1];
   }
}
/******************************************************************************/
/* this routine removes a covariate knot from the space */

static void cleanup1(i0,ki,spc,ncov)
struct space *spc;
int i0,ki,ncov;

/* spc   - space to be reduced
   i0    - which covariate
   ki    - ranknumber of knot
   ncov  - number of covariates */
{
   int j;
/* j   - counters  */

/* change the number of knots */
   (*spc).sub[i0][ncov].dim1-=1;

/* change all the other pointers */
   if((*spc).sub[i0][ncov].dim1>0){

/* first in the basisfunctions */
      for(j=0;j<(*spc).ndim;j++){
         if((*spc).basis[j].b1==i0&&(*spc).basis[j].t1>ki)(*spc).basis[j].t1-=1;
         if((*spc).basis[j].b2==i0&&(*spc).basis[j].t2>ki)(*spc).basis[j].t2-=1;
      }

/* in the knots themselves */
      for(j=ki;j>(-1)&&j<(*spc).sub[i0][ncov].dim1;j++){
         (*spc).sub[i0][ncov].ktsc[j]=(*spc).sub[i0][ncov].ktsc[j+1];
      }
   }
}
/******************************************************************************/
/* copy one basisfunction into another */

static void basisswap(bs1,bs2,ndata)
struct basisfunct *bs1,*bs2;
int ndata;

/* bs1   - basisfunction out
   bs2   - basisfunction in
   ndata - number of datapoints */
{
   int i;
/* i     - counter */

/* routine should be obvious */
   (*bs1).b1=(*bs2).b1;
   (*bs1).b2=(*bs2).b2;
   (*bs1).t1=(*bs2).t1;
   (*bs1).t2=(*bs2).t2;
   (*bs1).beta=(*bs2).beta;
   (*bs1).iknots=(*bs2).iknots;
   for(i=0;i<ndata;i++){
      (*bs1).values[i]=(*bs2).values[i];
      (*bs1).values2[i]=(*bs2).values2[i];
   }
}

static void uuu(spc,b1,b2,t1,t2,ncov,ii)
struct space *spc;
int b1,b2,t1,t2,ncov,ii;
{
   if(ii==0)(void)Rprintf("added: ");
   else (void)Rprintf("removed: ");
   if(b1!=ncov)(void)Rprintf("(%d",b1+1);
   else (void)Rprintf("(T");
   if(b2==ncov)(void)Rprintf(")=(");
   else (void)Rprintf(",%d)=(",b2+1);
   if(b1==ncov)(void)Rprintf("%.2f",(*spc).knots[t1]);
   else {if(t1!=-1)(void)Rprintf("%.2f",(*spc).sub[b1][ncov].ktsc[t1]);
                   else (void)Rprintf("linear");}
   if(b2==ncov)(void)Rprintf(") ");
   else {if(t2!=-1)(void)Rprintf(",%.2f)",(*spc).sub[b2][ncov].ktsc[t2]);
                   else (void)Rprintf(",linear) ");}
}
/******************************************************************************/
/* this routine contains the S-version I-O program, and a function it needs  */


/******************************************************************************/

/* the S-I/O routine */

void sharex(ncov,ndata)
int *ncov,*ndata;
{
   *ndata=MAXSPACE;
   *ncov=MAXKNOTS;
   return;
}


void share(ncov,ndata,time,delta,xcov,penalty,mindis,ndmax,bbtt,cckk,vexclude,
          lins,silent,logs,fitter,ad,strt)

double *time;
double *penalty,*logs;
int *delta,*silent,*fitter,*ad,*lins,*strt;
int *mindis,*ncov,*ndmax,*ndata,*vexclude;
double *xcov,*bbtt,*cckk;

/* time   - the times
   penalty- alpha(bic)
   delta  - delta
   mindis - minimum distance between knots (in order statistics)
   ncov   - number of covariates
   ndmax  - maximum number of basisfunctions
   ndata  - number of datapoints
   vexclude - exclude in vector form
   xcov   - the covariates
   bbtt   - summarizes the model (basisfunctions)
   cckk   - summarizes the model (1d subspaces) 
   strt   - start with constant (0) ot linear (1) 
   fitter - was a starting model supplied? (1=yes, 0=no)
   ad     - is the best model of this size during adding or deletion
   silent - should diagnostic output be printed? (0=yes, 1=no) */
{
   int i,j,k,k2,k3,**exclude;
   struct space *best;
   struct datastruct *data;

/* i,k           - counters 
   exclude       - which terms to exclude
   best          - the final fitted model
   definegspace  - allocates storage for a space
   definedata    - allocates storage for data
   data          - structure containing all information about the data
   hare          - that is where it happens
   soutgspace    - arranges the output */

/* if ncov<0 we only want the parameters */
   if(*ncov<0){
      *ndata=MAXSPACE;
      *ncov=MAXKNOTS;
      return;
   }

/* arrange the data */
   data=definedata(*ndata,*ncov);
   (*data).ndata=*ndata;
   (*data).ncov=*ncov;
   for(i=0;i<(*ndata);i++)(*data).times[i]=time[i];
   for(i=0;i<(*ndata);i++)(*data).delta[i]=delta[i];
   exclude=igmatrix((*data).ncov+2,(*data).ncov+2);
   for(i=0;i<(*data).ncov;i++){
      for(k=0;k<(*data).ndata;k++) (*data).cov[i][k]=xcov[i*(*data).ndata+k];
   }
   hareallocer((*data).ndata);

/* check whether things are binary */
   for(k=0;k<(*data).ncov;k++){
      (*data).bincov[k]=1;
      for(i=0;i<(*data).ndata;i++){
         if((*data).cov[k][i]!=0. && (*data).cov[k][i]!=1){
            (*data).bincov[k]=0;
            i=(*data).ndata;
         }
      }
   }

/* initialize exclude */
   for(k=0;k<=(*data).ncov;k++){
      for(i=0;i<=(*data).ncov;i++)exclude[i][k]=0;
   }
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

/* compute the same vector - are the covariates the same */
   (*data).same[0]=0;
   for(i=1;i<(*data).ndata;i++){
      (*data).same[i]=1;
      for(k=0;k<(*data).ncov;k++){
         if((*data).cov[k][i]!=(*data).cov[k][i-1]){
            k=(*data).ncov;
            (*data).same[i]=0;
         }
      }
   }

/* make a space */
   best=definegspace((*data).ncov,(*data).ndata);

/* we should initialize the space if we started with a fit */
   if(*fitter>0){
      (*best).ndim=*fitter;
      (*best).nknots=cckk[0];
      for(i=0;i<(*best).nknots;i++) (*best).knots[i]=cckk[(i+1)*(*ncov+1)];
      for(i=0;i<(*ncov);i++){
         (*best).sub[i][*ncov].dim1=cckk[i+1]+1;
         for(k=0;k<(*best).sub[i][*ncov].dim1-1;k++){
            (*best).sub[i][*ncov].ktsc[k]=cckk[(k+1)*(*ncov+1)+1+i];
         }
      }
      for(i=0;i<(*best).ndim;i++){
         j=bbtt[i*6+0];
         k=bbtt[i*6+2];
         (*best).basis[i].t1=bbtt[i*6+1];
         (*best).basis[i].t2=bbtt[i*6+3];
         (*best).basis[i].beta=bbtt[i*6+4];
         if(j==0) (*best).basis[i].b1=(*ncov);
         else (*best).basis[i].b1=j-1;
         if(k<=0) (*best).basis[i].b2=(*ncov);
         else(*best).basis[i].b2=k-1;
         if((*best).basis[i].t1> -1)((*best).basis[i].t1)-=1;
         if((*best).basis[i].t2> -1)((*best).basis[i].t2)-=1;
         if((*best).basis[i].b1==(*ncov)){
            (*best).basis[i].iknots=(*best).basis[i].t1;
         }
         else{
            ((*best).basis[i].iknots)=-1;
         }
         if((*best).basis[i].b2!=(*ncov)){
             j=(*best).basis[i].b1;
             k=(*best).basis[i].b2;
             k2=(*best).basis[i].t1;
             k3=(*best).basis[i].t1;
             (*best).sub[j][k].kts1[k2+1][k3+1]=1;
         }
      }
      for(i=0;i<(*ncov);i++){
         if((*best).sub[i][*ncov].dim1==1){
            for(j=0;j<(*best).ndim;j++){
               if((*best).basis[j].b1==i) j=((*best).ndim)+100;
            }
            if(*strt==0){
               if(j<((*best).ndim)+50)(*best).sub[i][*ncov].dim1=0;
            }
            else{
               if(j<((*best).ndim)+50){
                  (*best).basis[(*best).ndim].b1=i;
                  (*best).basis[(*best).ndim].b2=(*ncov);
                  (*best).basis[(*best).ndim].t1=-1;
                  (*best).basis[(*best).ndim].t2=-1;
                  (*best).basis[(*best).ndim].iknots=-1;
                  (*best).basis[(*best).ndim].beta=0.;
                  ((*best).ndim)+=1;
               }
            }
         }
      }
      getvector(best,(*data).ndata,(*ncov),(*data).cov,(*data).times);
      *strt = -1;
   }

/* do the work */
   hare(best,data,*penalty,*ndmax,*mindis,exclude,*strt,*silent,logs,ad,lins);

/* admire the results */
   soutgspace(best,data,bbtt,cckk);
   *ndata=(*best).ndim; 
   return;
}

/******************************************************************************/

/* this is an output routine, it writes the matrices bbtt and cckk - which are
   given as output to S */

static void soutgspace(spc,data,bbtt,cckk)
struct space *spc;
struct datastruct *data;
double *cckk,*bbtt;

/* spc   - structure describing the model
   data  - data 
   bbtt  - matrix describing the basisfunctions
   cckk  - matrix describing the 1d subspaces */

{
   int i,j,k,l;

   l=MAXKNOTS+1;

/* cckk */

/* second line - 0, number of time knots, time knots */
   cckk[0]=(*spc).nknots;
   for(k=0;k<(*spc).nknots;k++) cckk[k+1]=(*spc).knots[k];
   for(k=(*spc).nknots+1;k<=MAXKNOTS;k++)cckk[k]=0.;

/* lines 3-ncov+2- covariate number, dimension of space (#kts+1), kts */
   for(j=0;j<(*data).ncov;j++){
      cckk[(j+1)*l]=(*spc).sub[j][(*data).ncov].dim1-1;
      for(i=0;i<(*spc).sub[j][(*data).ncov].dim1-1;i++){
         cckk[(j+1)*l+i+1]=(*spc).sub[j][(*data).ncov].ktsc[i];
      }
      for(k=(*spc).sub[j][(*data).ncov].dim1;k<=MAXKNOTS;k++)cckk[(j+1)*l+k]=0.;
   }

/* rest: basisfunctions  (bbtt) */
   for(j=0;j<(*spc).ndim;j++){
      if((*spc).basis[j].b1>=0){
         (*spc).basis[j].b1+=1;
         if((*spc).basis[j].b1>(*data).ncov) (*spc).basis[j].b1=0;
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
      bbtt[j*6+0]=(*spc).basis[j].b1;
      bbtt[j*6+1]=(*spc).basis[j].t1;
      bbtt[j*6+2]=(*spc).basis[j].b2;
      bbtt[j*6+3]=(*spc).basis[j].t2;
      bbtt[j*6+4]=(*spc).basis[j].beta;
      bbtt[j*6+5]=(*spc).basis[j].se;
   }
   return;
}

/******************************************************************************/
static void getvector(best,ndata,ncov,cov,data)
struct space *best;
int ndata,ncov;
double **cov,*data;

/* best  - the model
   ndata - number of datapoints
   ncov  - number of covariates
   cov   - covariates */
{
   int i,j,k,b1,b2,t1,t2;
   double xx,*vv,*ww;

/* i,j   - counter
   b1,b2,t1,t2 - b1,b2,t1,t2 for the present basisfunction
   vv,ww - save typing
   xx    - the second half */

/* circle the basisfunctions */
   for(i=0;i<(*best).ndim;i++){

/* allocate storage */
      vv=dgvector(ndata);
      ww=dgvector(ndata);

/* if it is time only it is easy */
      for(j=0;j<ndata;j++){
         b1=(*best).basis[i].b1;
         t1=(*best).basis[i].t1;
         if(b1==ncov){
            vv[j]=1.;
         }
         else{

/* then first take the first component of the basisfunction */
            vv[j]=cov[b1][j];

/* -1 means linear, otherwise it is piecewise linear */
            if(t1> -1){
               vv[j]-=(*best).sub[b1][ncov].ktsc[t1];
               if(vv[j]<0.) vv[j]=0.;
            }
         }
/* and then the second component of the basisfunction */
         b2=(*best).basis[i].b2;
         t2=(*best).basis[i].t2;
         if(b2!=ncov && b2!=-1){
            xx=cov[b2][j];

/* -1 means linear, otherwise it is piecewise linear */
            if(t2> -1){
               xx-=(*best).sub[b2][ncov].ktsc[t2];
               if(xx<0.) xx=0.;
            }
            vv[j]=vv[j]*xx;
         }
/* now values 2 */
         ww[j]=vv[j];
         t1=(*best).basis[i].t1;
         if(b1==ncov && t1>=0){
            xx= -data[j];
            xx+=(*best).knots[t1];
            if(xx<0.) xx=0.;
            ww[j]=ww[j]*xx;
         }
      }
      (*best).basis[i].values=vv;
      (*best).basis[i].values2=ww;
   }
   for(i=0;i<(*best).ndim;i++){
      for(j=i;j<(*best).ndim;j++){
         (*best).xtx[i][j]=0.;
         for(k=0;k<ndata;k++){
            (*best).xtx[i][j]+=(*best).basis[i].values2[k]*
                               (*best).basis[j].values2[k];
         }
      }
   }
   for(i=1;i<(*best).ndim;i++){
      for(j=0;j<i;j++) (*best).xtx[i][j]=(*best).xtx[j][i];
   }
   return;
}
/******************************************************************************/
/* this function allocates storage for a data structure */

static struct datastruct *definedata(ndata,ncov)
int ncov,ndata;
{
   struct datastruct *newdata;
   /* newdata=(struct datastruct *)S_alloc((long)1,sizeof(struct datastruct)); */
   newdata=(struct datastruct *)Salloc(1,struct datastruct);
   (*newdata).delta=igvector(ndata);
   (*newdata).bincov=igvector(ndata);
   (*newdata).same=igvector(ndata);
   (*newdata).times=dgvector(ndata);
   (*newdata).cov=dgmatrix(ncov+1,ndata);
   return newdata;
}
/******************************************************************************/

/* this function allocates storage for a space */

static struct space *definegspace(ncov,ndata)
int ncov,ndata;

/* ncov  - number of covariates */

{
   struct space *newspace;
   int i,j,k,l;

/* this routine is mainly a copy from the numerical recipes routines
   i,j,k       - counters
   newspace    - thing to be initialized
   definebasis - initializes an array of basisfunctions
   definedim   - initializes a matrix of subdimensions */

/* basic allocation */
   /* newspace=(struct space *)S_alloc((long)1,sizeof(struct space)); */
   newspace=(struct space *)Salloc(1,struct space);

/* the simple elements */
   (*newspace).knots=dgvector(MAXKNOTS);
   (*newspace).info=dgmatrix(MAXSPACE,MAXSPACE);
   (*newspace).score=dgvector(MAXSPACE);
   (*newspace).b0=dgmatrix(MAXKNOTS+1,MAXSPACE+1);
   (*newspace).b1=dgmatrix(MAXKNOTS+1,MAXSPACE+1);
   (*newspace).b2=dgvector(MAXSPACE+1);
   (*newspace).xtx=dgmatrix(MAXSPACE,MAXSPACE);
   (*newspace).ndim=0;
   (*newspace).nknots=0;
   (*newspace).aic=0.;

/* defines the basisfunctions, initializes */
   (*newspace).basis=definebasis();
   for(i=0;i<MAXSPACE;i++) {
      (*newspace).basis[i].values=dgvector(ndata);
      (*newspace).basis[i].values2=dgvector(ndata);
      (*newspace).basis[i].b1=-1;
      (*newspace).basis[i].b2=-1;
      (*newspace).basis[i].t1=-1;
      (*newspace).basis[i].t2=-1;
      (*newspace).basis[i].iknots=-1;
      (*newspace).basis[i].beta=0.;
      (*newspace).basis[i].se=0.;
      for(j=0;j<MAXSPACE;j++) (*newspace).xtx[i][j]=0.;
   }

/* defines the subdimensions */
   (*newspace).sub=definedim(ncov+1);

/* for the 2-covariate subdimensions */
   for(i=0;i<ncov;i++) for(j=i+1;j<ncov;j++){
      (*newspace).sub[i][j].kts1=iigmatrix(MAXKNOTS+1,MAXKNOTS+1);
      (*newspace).sub[i][j].dim1=0;
      for(k=0;k<MAXKNOTS+1;k++) for(l=0;l<MAXKNOTS+1;l++){
         (*newspace).sub[i][j].kts1[k][l]=0;
      }
   }

/* for the 1-covariate subdimensions */
   for(j=0;j<=ncov;j++){
      (*newspace).sub[j][ncov].ktsc=ddgvector(MAXKNOTS);
      (*newspace).sub[j][ncov].dim1=0;
   }

/* for the 1-covariate and time subdimensions */
   for(j=0;j<=ncov;j++){
      (*newspace).sub[ncov][j].kts1=iigmatrix(MAXKNOTS+1,MAXKNOTS+1);
      (*newspace).sub[ncov][j].dim1=0;
      for(k=0;k<MAXKNOTS+1;k++) for(l=0;l<MAXKNOTS+1;l++){
         (*newspace).sub[ncov][j].kts1[k][l]=0;
      }
   }

   return newspace;
}

/******************************************************************************/

/* this function allocates storage for an array of basisfunctions */
static struct basisfunct *definebasis()
{
   struct basisfunct *nb;

   nb=(struct basisfunct *)Salloc(MAXSPACE,struct basisfunct);
   return nb;
}

/******************************************************************************/

/* this function allocates storage for a matrix of subdimensions */
static struct subdim **definedim(ncov)
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
static int *igvector(l)
int l;
/* allocate an int vector with subscript range v[0...l] */
{
   int *v,i;
   v=(int *)Salloc(l+1,int);
   for(i=0;i<=l;i++)v[i]=0;
   return v;
}
/******************************************************************************/
static float *ddgvector(l)
int l;
/* allocate a float vector with subscript range v[0...l] */
{
   float *v;
   int i;
   v=(float *)Salloc(l+1,float);
   for(i=0;i<=l;i++)v[i]=0.;
   return v;
}
/******************************************************************************/
static double *dgvector(l)
int l;
/* allocate a double vector with subscript range v[0...l] */
{
   double *v;
   int i;
   v=(double *)Salloc(l+1,double);
   for(i=0;i<=l;i++)v[i]=0.;
   return v;
}
/******************************************************************************/
static short int **iigmatrix(r,c)
int r,c;
/* allocate an int matrix with subscript range m[0..r][0..c] */
{
   short int i,j,**m;
   m=(short int **) Salloc(r+1,short int*);
   for(i=0;i<=r;i++){
      m[i]=(short int *) Salloc(c+1,short int);
      for(j=0;j<=c;j++)m[i][j]=0;
   }
   return m;
}
/******************************************************************************/
static int **igmatrix(r,c)
int r,c;
/* allocate an int matrix with subscript range m[0..r][0..c] */
{
   int i,j,**m;
   m=(int **) Salloc(r+1,int*);
   for(i=0;i<=r;i++){
      m[i]=(int *) Salloc(c+1,int);
      for(j=0;j<=c;j++)m[i][j]=0;
   }
   return m;
}
/******************************************************************************/
static double **dgmatrix(r,c)
int r,c;
/* allocate a double matrix with subscript range m[0..r][0..c] */
{
   int i,j;
   double **m;
   m=(double **) Salloc(r+1,double*);
   for(i=0;i<=r;i++){
      m[i]=(double *) Salloc(c+1,double);
      for(j=0;j<=c;j++)m[i][j]=0.;
   }
   return m;
}
/******************************************************************************/
static void nrerror(error_text)
char error_text[];
{
   void exit();

   (void)error("%s\n this is serious!");
   exit(1);
}
/******************************************************************************/
static int humbertester(aa)
double aa;
/* if aa = -Inf: 0
      aa = +Inf: 1
      aa =  NaN: 2
      otherwise: 3 */
{
   int i1=0,i2=0,i3=0,i4=0;
   if(aa< 2.)i1=1;
   if(aa> 0.)i2=1;
   if(aa< pow(10.,500.))i3=1;
   if(aa> -pow(10.,500.))i4=1;
   if(i1+i2+i3+i4>=3)return 3;
   if(i2==1 && i4==1)return 1;
   if(i1==1 && i3==1)return 0;
   return 2;
}
/******************************************************************************/
static void glusolve(a,n,b)
int n;
double **a,*b;
{
   double aa[DIM5][DIM5],bb[DIM5];
   int kpvt[DIM5],info;
   int i,j;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
      bb[i]=b[i];
   }
   i=DIM5;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
   if(info!=0)nrerror("info in glusolve is not 0");
   F77_CALL(xdsisl)(aa,&i,&n,kpvt,bb);
   for(i=0;i<n;i++)b[i]=bb[i];
}
/******************************************************************************/
static int glusolve2(a,n,b)
int n;
double **a,*b;
{
   double aa[DIM5][DIM5],bb[DIM5];
   int kpvt[DIM5],info;
   int i,j;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
      bb[i]=b[i];
   }
   i=DIM5;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
   if(info!=0)return 0;
   F77_CALL(xdsisl)(aa,&i,&n,kpvt,bb);
   for(i=0;i<n;i++)b[i]=bb[i];
   return 1;
}
/******************************************************************************/
static void gluinverse(a,n)
int n;
double **a;
{
   double aa[DIM5][DIM5],bb[DIM5],det[2];
   int kpvt[DIM5],info,inert[3];
   int i,j;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
   }
   i=DIM5;
   j=1;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
   F77_CALL(xdsidi)(aa,&i,&n,kpvt,det,inert,bb,&j);
   for(i=0;i<n;i++){
      for(j=0;j<i;j++)a[i][j]=aa[i][j];
      for(j=i;j<n;j++)a[i][j]=aa[j][i];
   }
}
/******************************************************************************/
static void dsort(ra,n)
int n;
double *ra;
{
   double w[20000],w2[20000];
   int i;
   if(n>20000)nrerror("increase dim(w) in dsort");
   for(i=0;i<n;i++)w[i]=ra[i];
   i=1;
   F77_CALL(xssort)(w,w2,&n,&i);
   for(i=0;i<n;i++)ra[i]=w[i];
} 
/******************************************************************************/
static void sort(ra,rb,n)
int n;
double *ra,*rb;
{
   double w[20000],w2[20000];
   int i;
   if(n>20000)nrerror("increase dim(w) in dsort");
   for(i=0;i<n;i++)w[i]=rb[i];
   i=1;
   F77_CALL(xssort)(w,w2,&n,&i);
   for(i=0;i<n;i++)ra[i]=w[i];
} 
/******************************************************************************/
static double condition(a,n)
int n;
double **a;
{
   double aa[DIM5][DIM5],bb[DIM5],rcond;
   int kpvt[DIM5];
   int i,j;
   for(i=0;i<n;i++) for(j=0;j<n;j++)aa[i][j]=a[j][i];
   i=DIM5;
   j=1;
   F77_CALL(xdsico)(aa,&i,&n,kpvt,&rcond,bb);
   return rcond;
}
/******************************************************************************/
static void hareallocer(ndata)
int ndata;
{
   newtonwhere=igvector(ndata+1);
   searchsorted=dgvector(ndata+1);
   searchkts=dgvector(DIM5);
   searchsorted2=dgvector(ndata+1);
   remdimy=searchsorted;
   raoss=dgvector(DIM5);
   raoscorecopy=dgvector(DIM5);
   newtonscp=dgvector(DIM5);
   compallss=dgvector(DIM5);
   complogbasis0=dgvector(DIM5);
   complogbasis1=dgvector(DIM5);
   remdimxty=dgvector(DIM5);
   raohhh=dgmatrix(DIM5,DIM5);
   getsexx=dgmatrix(DIM5,DIM5);
   compallhhh=dgmatrix(DIM5,DIM5);
   remdimxtx=dgmatrix(DIM5,DIM5);
    
}
/******************************************************************************/
/* this function prints out the tables for hare.summary */

void ssumm(penalty,sample,logl,llogl,spcs,fcts,ndim,ncov)
int *sample,*llogl,*ncov,*ndim;
double *penalty,*logl,*spcs,*fcts;

/* penalty   - penalty term used
   sample    - sample size
   logl      - vector of loglikelihoods
               also where these ads or delete models
   llogl     - largest size
   spcs      - matrix containing information about the knots
   fcts      - matrix containing information about the basis functions
   ndim      - number of basis functions
   ncov      - number of covariates */
{
   int i,*k1,*k4,j,k;
   double *d1,*d2,*d3,*d4;

/* k1,k4,d1,d2,d3,d4 are the various columns of the first table 
   k1 - dimension
   k4 - add or del model
   d1 - loglikelihood
   d2 - bic
   d3 & d4 optimality bounds */

   k1=igvector(*llogl+3); d4=dgvector(*llogl+3); k4=igvector(*llogl+3);
   d1=dgvector(*llogl+3); d2=dgvector(*llogl+3); d3=dgvector(*llogl+3);
   j=0;

/* compute aic and loglikelihood */
   for(i=0;i<(*llogl);i++){
      if(logl[2*i]!=0){
         k1[j]=i+1; d1[j]=logl[2*i];
         k4[j]=logl[2*i+1];
         d2[j]=-2*d1[j]+(*penalty)*k1[j];
         d3[j]=-2.; d4[j]=-1.;
         j++;
      }
   }

/* compute optimality bounds */
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
   for(i=1;i<j-1;i++){
      if(d4[i]<d3[i]){
         d4[i]=-3.;
         d3[i]=-3.;
      }
   }
   (void)Rprintf("dim A/D   loglik       AIC        penalty \n");
   (void)Rprintf("                                min    max \n");
   k=0;
   for(i=0;i<j;i++){
      if(d2[i]<d2[k])k=i;
      if(k4[i]==0)(void)Rprintf("%3d Del %9.2f %9.2f",k1[i],d1[i],d2[i]);
      else (void)Rprintf("%3d Add %9.2f %9.2f",k1[i],d1[i],d2[i]);
      if(d3[i]>0)(void)Rprintf(" %7.2f",2*d3[i]);
      if(i==0)(void)Rprintf("     Inf");
      if(d3[i]<0 && i!=0 && i!=j-1)(void)Rprintf("      NA");
      if(i==j-1)(void)Rprintf("    0.00");
      if(d4[i]>0)(void)Rprintf(" %7.2f",2*d4[i]);
      if(d4[i]<0 && i!=j-1 && i!=0)(void)Rprintf("      NA");
      (void)Rprintf("\n");
   }
   (void)Rprintf("\nthe present optimal number of dimensions is %d.\n",k1[k]);
   if((int)exp(*penalty) == *sample){
      (void)Rprintf("penalty(AIC) was the default: BIC=log(samplesize): log(");
      (void)Rprintf("%d)=%.2f\n",*sample,*penalty);
   }
   else{
      (void)Rprintf("penalty(AIC) was %.2f",*penalty);
      (void)Rprintf(", the default (BIC), would have been %.2f.\n",
            log((double)*sample));
   }
   if(k1[0]>1){
      (void)Rprintf("models with fewer than %d dims ",k1[0]);
      (void)Rprintf("can be fitted, but they are not optimal for the\n");
      (void)Rprintf("present choice of penalty - choose penalty in ");
      (void)Rprintf("hare.fit larger to see these fits.\n");
   }
   (void)Rprintf("\n");
   (void)Rprintf("  dim1           dim2           beta        SE         Wald\n");
   for(i=0;i<(*ndim);i++){
      if(i==0) (void)Rprintf("Constant      ");
      else{
         if((int)fcts[i*6]==0)(void)Rprintf("Time");
         else (void)Rprintf("Co-%d",(int)fcts[i*6]);
         if((int)fcts[i*6+1]==0)(void)Rprintf("  linear  ");
         else (void)Rprintf(" %9.2g",spcs[(int)(fcts[i*6]+fcts[i*6+1]*(1+(*ncov)))]);
      }
      if((int)fcts[i*6+2]<=0)(void)Rprintf("               ");
      else{
         (void)Rprintf(" Co-%d",(int)fcts[i*6+2]);
         if((int)fcts[i*6+3]==0)(void)Rprintf("  linear  ");
        else (void)Rprintf(" %9.2g",spcs[(int)(fcts[i*6+2]+fcts[i*6+3]*(1+(*ncov)))]);
      }
      d1[0]=(fcts[i*6+4]/fcts[i*6+5]);
      (void)Rprintf(" %10.2g %10.2g %7.2f\n",fcts[i*6+4],fcts[i*6+5],d1[0]);
   }
}
/******************************************************************************/

/* this routine contains phare/dhare/qhare/share   = S routine                */

/******************************************************************************/

/* the C-I/O routine */

void sphare(ncov,ndim,ndata,xcov,ip,pdh,cckk,bbtt)
double *pdh,*xcov,*cckk,*bbtt;
int *ncov,*ndim,*ip,*ndata;

/* ncov          - number of covariates
   ndim          - dimensionality of the space
   ndata         - number of datapoints
   xcov          - covariates
   ip            - do we want quantiles (ip=1) or probabilities (ip=0)
                   or hazards (ip=2) or densities (ip=3)
   pdh           - probabilities/density/hazard/quantiles
   cckk          - the info about the 1d subspaces 
   bbtt          - the info about the basisfunctions */

{
   int i,j,k;
   struct space *best;
   double **cov,*qqq;

/* i,j,k         - counters     
   best          - the final fitted model
   hdefinegspace   - allocates storage for a space
   houtgspace     - gets the output (ip=0/2/3)
   poutgspace     - gets the output (ip=1)
   getvectors    - defines the vector part of the basisfunctions 
   cov           - covariates 
   qqq           - quantiles */
   
/* allocate space for the model */
   qqq=dgvector(*ndata);
   if(ncov==0) best=hdefinegspace(2,*ndata);
   else best=hdefinegspace(*ncov,*ndata);

/* essentially we need two copies */
   for(i=0;i<(*ndata);i++) qqq[i]=pdh[i];

/* read in the model - spaces */
   (*best).ndim=*ndim;
   (*best).nknots=cckk[0];
   for(i=0;i<(*best).nknots;i++) (*best).knots[i]=cckk[(i+1)*(*ncov+1)];
   for(i=0;i<(*ncov);i++){
      (*best).sub[i][*ncov].dim1=cckk[i+1]+1;
      for(k=0;k<(*best).sub[i][*ncov].dim1-1;k++){
         (*best).sub[i][*ncov].ktsc[k]=cckk[(k+1)*(*ncov+1)+1+i];
      }
   }

/* read in the model - basis functions */
   for(i=0;i<(*ndim);i++){
      j=bbtt[i*6+0];
      k=bbtt[i*6+2];
      (*best).basis[i].t1=bbtt[i*6+1];
      (*best).basis[i].t2=bbtt[i*6+3];
      (*best).basis[i].beta=bbtt[i*6+4];

/* adjust certain indices time becomes n (was 0) and covariates 0 to (n-1) 
   (was 1 to n) */
      if(j==0) (*best).basis[i].b1=(*ncov);
      else (*best).basis[i].b1=j-1;
      if(k<=0) (*best).basis[i].b2=(*ncov);
      else(*best).basis[i].b2=k-1;

/* adjust the knots by one too */
      ((*best).basis[i].t1)-=1;
      ((*best).basis[i].t2)-=1;
      if((*best).basis[i].b1==(*ncov)){
         (*best).basis[i].iknots=(*best).basis[i].t1;
      }
      else{
         ((*best).basis[i].iknots)=-1;
      }
   }

/* put the covariates in a matrix */
   cov=dgmatrix((*ncov)+1,*ndata);
   for(i=0;i<(*ncov);i++){
      for(j=0;j<(*ndata);j++) cov[i][j]=xcov[i*((*ndata))+j];
   }

/* get vectors  component of basisfunctions */
   getvectors(best,*ndata,*ncov,cov);

/* get probabilities/densities/hazard from quantiles case */
   if((*ip)!=1){
      houtgspace(best,pdh,qqq,*ndata,*ip);
   }

/* get quantiles from probabilities case */
   else{
       poutgspace(best,pdh,qqq,*ndata);
       for(i=0;i<(*ndata);i++)pdh[i]=qqq[i];
   }
}

/******************************************************************************/

/* this routine is used in the qhare case */

static void poutgspace(spc,ppp,qqq,ndata)
struct space *spc;
double *ppp,*qqq;
int ndata;

/* spc   - the model fitted
   ppp   - the probabilities (in)
   qqq   - the quantiles (out)
   ndata - number of datapoints */
{
   int i,j,k;
   double *lin,*con;

/* i,j,k,l - counters
   lin     - linear term[i] is between knots (i-1) and (i)
   con     - constant term[i] is between knots (i-1) and (i)
   getthosep - gets those q from those p */

/* allocation */
   lin=dgvector((*spc).nknots+1);
   con=dgvector((*spc).nknots+1);

/* gets the linear and constant term relevant to the next set of covariates */
   for(i=0;i<ndata;i++){

/* between each set of knots */
      for(j=0;j<(*spc).nknots+1;j++){
         lin[j]=0.;
         con[j]=0.;

/* for each basisfunction */
         for(k=0;k<(*spc).ndim;k++){

/* constant in time */
            if((*spc).basis[k].iknots==-1){
               con[j]+=(*spc).basis[k].values[i]*(*spc).basis[k].beta;
            }

/* piecewise linear in time */
            else{
               if(j<=(*spc).basis[k].t1){
                  lin[j]-=(*spc).basis[k].beta*(*spc).basis[k].values[i];
                  con[j]+=(*spc).basis[k].beta*(*spc).basis[k].values[i]*
                          (*spc).knots[(*spc).basis[k].t1];
               }
            }
         }
      }

/* for which cases is this relevant */
      k=i;
/* get those numbers */
      getthosep(lin,con,(*spc).nknots,(*spc).knots,ppp,qqq,i,k+1);
      i=k;
   }
   return;
}
/******************************************************************************/

/* this one geths the quantiles for some p's */

static void getthosep(lin,con,nk,kts,ppp,qqq,i1,i2)
double *lin,*con,*ppp,*qqq,*kts;
int i1,i2,nk;

/* lin - linear term
   con - constant term
   nk  - number of knots
   kts - knots
   ppp - probabilities
   qqq - quantiles
   i1  - first case in action
   i2  - first case after last case in action */

{
   int i,j,doing=i1;
   double rold=0.,lower,rnew;

/* i,j   - counter
   doing - which job are we working on
   rold  - probability at previous knot
   eint  - computes integral
   lower - lower integration boundary
   xeint - computes inverse integral 
   rnew  - pronanility at new knot */  

/* transform the probabilities */
   for(j=i1;j<i2;j++)qqq[j]=-log((double)(1.-ppp[j]));

/* compute the probability in the next knot */
   for(i=0;i<nk;i++){
      if(i==0)lower=0.;
      else lower=kts[i-1];
      rnew=rold+eint(lin[i],con[i],lower,kts[i]);

/* if we are too far, invert an integral */
      for(j=doing;j<i2;j++){
         if(qqq[j]<rnew){
            qqq[j]=xeint(lin[i],con[i],lower,qqq[j]-rold);
            doing++;
         }
         else{
            j=i2;
         }
      }
      rold=rnew;
   }

/* those after the last knot */
   for(j=doing;j<i2;j++) qqq[j]=xeint(lin[nk],con[nk],kts[nk-1],qqq[j]-rold);
}

/******************************************************************************/

/* this routine computes pdh from qqq */

static void houtgspace(spc,pdh,qqq,ndata,ip)
struct space *spc;
double *pdh,*qqq;
int ndata,ip;

/* spc   - structure describing the model
   pdh   - probabilities/densities/hazard
   qqq   - quantiles
   ndata - how much data
   ip    - ip=0: probabilities, ip=2: hazard, ip=3: density */

{
   int i,j;
   double s,*b1,*b0;

/* i,j - counters 
   s     - part of the hazard/density
   hcomplog - computes the log density 
   b0,b1 - for complog */
   
/* allocation */
   b0=dgvector((*spc).nknots+1);
   b1=dgvector((*spc).nknots+1);

/* for each datapoint */
   for(i=0;i<ndata;i++){

/* hazard or probability */
      if(ip!=3){
         pdh[i]=0.;

/* per basisfunction */
         for(j=0;j<(*spc).ndim;j++){
            s=(*spc).basis[j].values[i]*(*spc).basis[j].beta;
            if((*spc).basis[j].iknots<0){
               pdh[i]+=s;
            }
            else{
               if(qqq[i]<(*spc).knots[(*spc).basis[j].iknots]){
                  pdh[i]+=s*((*spc).knots[(*spc).basis[j].iknots]-qqq[i]);
               }
            }
         }
         pdh[i]=exp(pdh[i]);
      }

/* density or probability */
      if(ip!=2) s=hcomplog(spc,qqq[i],i,b0,b1);

/* density */
      if(ip==3) pdh[i]=s;

/* probability */
      if(ip==0) pdh[i]=1-s/pdh[i];
   }
   return;
}

/******************************************************************************/

/* this function computes the density in one datapoint */

static double hcomplog(spc,data,idt,basis0,basis1)
struct space *spc;
double data,*basis0,*basis1;
int idt;

/* spc   - the present model
   data  - the data (only time) 
   idt   - which datapoint 
   basis0 - element constant term of lambda[between-2-knots,datapoint]
   basis1 - basis0, but linear term */

{
   double basis2,l0,l1,logl;
   int j,where;

/* basis2 - lambda in a datapoint 
   l0     - lower integration bound
   l1     - upper integration bound
   logl   - loglikelihood
   eint   - computes integral
   upbasis3- updates basis for one point 
   j      - counter
   where  - in which interval is the datapoint */

/* initializations */
   logl=0.;

/* initialize */
   for(j=0;j<=(*spc).nknots;j++){
      basis0[j]=0.;
      basis1[j]=0.;
   }
   basis2=0.;

/* in which interval is the data? */
   where=(*spc).nknots;
   for(j=0;j<(*spc).nknots;j++){
      if((*spc).knots[j]>data){
         where=j;
         j=(*spc).nknots;
      }
   }

/* update basis0, basis1 and basis2 */
   for(j=0;j<(*spc).ndim;j++) upbasis3((*spc).knots,basis0,basis1,&basis2,idt,
                                                &((*spc).basis[j]),where,data);

/* the delta part */
   logl+=basis2;
 
/* the integrals */
   for(j=0;j<=where;j++){

/* lower bound */
      if(j==0) l0=0.;
      else l0=(*spc).knots[j-1];

/* upper bound */
      if(j==where) l1=data;
      else l1=(*spc).knots[j];

/* integrals */
      logl-=eint(basis1[j],basis0[j],l0,l1);
   }

/* clean up */
   return exp(logl);
}

/******************************************************************************/

/* comprable to upbasis, but does less, since only the loglikelihood will
   be computed */

static void upbasis3(knots,basis0,basis1,basis2,idt,basf,where,time)
int idt,where;
double *basis0,*basis1,*basis2,*knots,time;
struct basisfunct *basf;

/* knots  - time-knots
   basis0 - element constant term of lambda[between-2-knots,datapoint]
   basis1 - basis0, but linear term
   basis2 - lambda in a datapoint
   idt    - number of datapoint on which we are working
   basf   - one basisfunction
   where  - in between which knots is the basisfunction 
   time   - datapoint */

{
   int j;
   double x;

/* j    - counter
   x    - save typing (see below) */

   x=(*basf).values[idt]*(*basf).beta;

/* no t-knots, basisfunction is constant */
   if((*basf).iknots== -1){
      *basis2+=x;
      for(j=0;j<=where;j++) basis0[j]+=x;
   }

/* 1 t-knot, basisfunction is 0 after knot, and has slope -1 before knot */
   else{
      if(time <= knots[(*basf).iknots]){
         *basis2+=x*(knots[(*basf).iknots]-time);
      }
      for(j=0;j<=where&&j<=(*basf).iknots;j++){
         basis0[j]+=knots[(*basf).iknots]*x;
         basis1[j]-=x;
      }
   }
}

/******************************************************************************/

/* this function computes the vector part of a basisfunction */

static void getvectors(best,ndata,ncov,cov)
struct space *best;
int ndata,ncov;
double **cov;

/* best  - the model
   ndata - number of datapoints
   ncov  - number of covariates
   cov   - covariates */
{
   int i,j,b1,b2,t1,t2;
   double xx,*vv;

/* i,j   - counter
   b1,b2,t1,t2 - b1,b2,t1,t2 for the present basisfunction
   vv    - save typing
   xx    - the second half */

/* circle the basisfunctions */
   for(i=0;i<(*best).ndim;i++){

/* allocate storage */

/* if it is time only it is easy */
      vv=dgvector(ndata+5);
      for(j=0;j<ndata;j++){
         b1=(*best).basis[i].b1;
         t1=(*best).basis[i].t1;
         if(b1==ncov){
            vv[j]=1.;
         }
         else{

/* then first take the first component of the basisfunction */
            vv[j]=cov[b1][j];

/* -1 means linear, otherwise it is piecewise linear */
            if(t1> -1){
               vv[j]-=(*best).sub[b1][ncov].ktsc[t1];
               if(vv[j]<0.) vv[j]=0.;
            }
         }
/* and then the second component of the basisfunction */
         b2=(*best).basis[i].b2;
         t2=(*best).basis[i].t2;
         if(b2!=ncov && b2!=-1){
            xx=cov[b2][j];

/* -1 means linear, otherwise it is piecewise linear */
            if(t2> -1){
               xx-=(*best).sub[b2][ncov].ktsc[t2];
               if(xx<0.) xx=0.;
            }
            vv[j]=vv[j]*xx;
         }
      }
      (*best).basis[i].values=vv;
   }
   return;
}
/******************************************************************************/


/*                t
                 /
                 | (b1*x+b2)
solves for t  c= |e          dx
                 |
                 /
                l            */

static double xeint(b1,b2,l,c)
double b1,b2,l,c;

/* just work it out */
{
   if(b1!=0.){
      c=log(c*b1/exp(b2)+exp(b1*l))/b1;
      return c;
   }
   return l+c/exp(b2);
}
/******************************************************************************/


/* this function allocates storage for a space */

static struct space *hdefinegspace(ncov,ndata)
int ncov,ndata;

/* ncov  - number of covariates */

{
   struct space *newspace;
   int i,j;

/* this routine is mainly a copy from the numerical recipes routines
   i,j,k       - counters
   newspace    - thing to be initialized
   definebasis - initializes an array of basisfunctions
   definedim   - initializes a matrix of subdimensions */

/* basic allocation */
   newspace=(struct space *)Salloc(1,struct space); 
   if(!newspace) nrerror("allocation error in definegspace");

/* the simple elements */
   (*newspace).knots=dgvector(MAXKNOTS);
   (*newspace).info=dgmatrix(MAXSPACE,MAXSPACE);
   (*newspace).score=dgvector(MAXSPACE);
   (*newspace).b0=dgmatrix(MAXKNOTS+1,MAXSPACE+1);
   (*newspace).b1=dgmatrix(MAXKNOTS+1,MAXSPACE+1);
   (*newspace).b2=dgvector(MAXSPACE+1);

/* defines the basisfunctions, initializes */
   (*newspace).basis=hdefinebasis();
   for(i=0;i<MAXSPACE;i++) {
      (*newspace).basis[i].values=dgvector(ndata);
      (*newspace).basis[i].b1=-1;
      (*newspace).basis[i].b2=-1;
      (*newspace).basis[i].t1=-1;
      (*newspace).basis[i].t2=-1;
      (*newspace).basis[i].iknots=-1;
      (*newspace).basis[i].beta=0.;
      (*newspace).basis[i].se=0.;
   }

/* defines the subdimensions */
   (*newspace).sub=hdefinedim(ncov+1);

/* for the 2-covariate subdimensions */
   for(i=0;i<ncov;i++) for(j=i+1;j<ncov;j++){
      (*newspace).sub[i][j].dim1=0;
   }

/* for the 1-covariate subdimensions */
   for(j=0;j<=ncov;j++){
      (*newspace).sub[j][ncov].ktsc=ddgvector(MAXKNOTS);
      (*newspace).sub[j][ncov].dim1=0;
   }

/* for the 1-covariate and time subdimensions */
   for(j=0;j<=ncov;j++){
      (*newspace).sub[i][j].dim1=0;
   }

   return newspace;
}

/******************************************************************************/

/* this function allocates storage for an array of basisfunctions */
static struct basisfunct *hdefinebasis()
{
   struct basisfunct *nb;

   nb=(struct basisfunct *)Salloc(MAXSPACE,struct basisfunct); 
   return nb;
}

/******************************************************************************/

/* this function allocates storage for a matrix of subdimensions */
static struct subdim **hdefinedim(ncov)
int ncov;
{
   struct subdim **newdim;
   int i;
   newdim=(struct subdim **)Salloc(ncov,struct subdim *);
   for(i=0;i<=ncov;i++)
      newdim[i]=(struct subdim *)Salloc(ncov,struct subdim);
   return newdim;
}
/******************************************************************************/
