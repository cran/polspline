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

#define MAXKNOTS 35
#define HLENGTH MAXKNOTS+5

void F77_NAME(xdsifa)(double[][HLENGTH], int *, int *, int *, int *);
void F77_NAME(xdsisl)(double[][HLENGTH], int *, int *, int *, double *);
void F77_NAME(xdsidi)(double[][HLENGTH], int *, int *, int *, double *, int *, double *, int *);
void F77_NAME(xdgefa)(double[][HLENGTH], int *, int *, int *, int *);
void F77_NAME(xdgesl)(double[][HLENGTH], int *, int *, int *, double *, int *);


static double hmylog();


/* MAXKNOTS is the maximum number of knots in a model
   HLENGTH   is the generic vector length */

static void hlusolve(),hluinverse(),hlusolve2();
static int *ihvector(),**ihmatrix();
static double *dhvector(),**dhmatrix(),***dstriparray();
static double *wkddd,*wkvec1,*wkvec2,**wkmat1,*wkphi,**wkmat,*wkphi2,*wkxx,*wkcand;
static double *wkmasterpt,*wkphi3,*wkse3,*wkphi4,**wkhh,**wkpowdat,**wkpowvec;
static double **wkinfo2,*wkscore2,*wkscore3,*wknewbas,*wknewdata,*wkphi7,**wkmat33;
static double *wksorted;
static void heft(),hstart2();
static void intprep(),getcoef(),start(),thetaswap(),tossit(),nstart(),getcoefx();
static void hiter(),hknotplace(),dubmodel(),hetse(),allocer();
static int add(),hlocation();
static void midblob(),basis(),lgrange();
static int hopplus();
static double summer(),ilambda(),xlambda();
static int step();
static double summer2(),lambda();
static int step2();
static void hremoveknot(),getse2();
static void thetaform(),newnew();
static int hindyl(),hindyr(),hindl(),hindr(),hindx(),hindm();
static double hrao();

struct model {
   int nk,*iknots,**icoef,nk1,*ad;
   double *knots,*theta,**coef2,***coef3,aic,*score,**hessian,*logl,tailse[2];
   double *basvec,**basmat,*mult,*tails,*yknots,ll;
};

/* nk     - number of knots in a model
   iknots - which of the potential knots in yknots are a member of knots
            (1=yes, 0=no)
   icoef  - does basisfunction i exist in interval j (related to coef2/coef3)
   nk1    - largest number of knots fitted
   ad     - vector of 0/1/2 which indicate whether the best model of the
            corresponding dimension was not fit (2), fit during addition (0)
            or deletion (1)
   knots  - the knots
   theta  - theta-hat
   coef2  - representation of the basis into truncated polynomial representation
   coef3  - alternative representation of the basis into truncated polynomials
   aic    - aikaike criterion of this model
   score  - score function
   hessian- hessian
   logl   - vector: for those elements of ad that are 0 or 1 the loglikelihood
            of the corrsponding model
   tailse - standard errors of the log terms
   basvec - used for the numerical integration
   basmat - used for the numerical integration
   mult   - used for the numerical integration
   tails  - vector indicating the status of the tail basisfunctions
   yknots - all knots ever used during the analysis, knots is a subset of yknots
   ll                                                                         */

struct datas {
   int nd,*delta;
   double *data,cc,**basdata1,**basdata2;
};

/* nd     - length of the data
   delta  - vector with delta values
   data   - sorted vector of observations
   cc     - c=number
   basdata1 - used for numerical integrations
   basdata2 - used for numerical integrations  */
static struct model *makemodel();
static struct datas *makedata();
/******************************************************************************/
/* this routine looks ugly - and is ugly. It only relates the S-variables to
   the C-variables - almost all should be self explanatory                    */
void sheftx(nx)
int *nx;
{
   *nx=HLENGTH;
   return;
}

void sheft(nx,data,delta,nkstart,knots,alpha,tails,iauto,logl,theta,iknots,zerror,cc,
      nkmax,ad,mindist)
int *nx,*delta,*nkmax,*nkstart,*zerror,*iknots,*iauto,*ad,*mindist;
double *data,*alpha,*theta,*cc,*logl,*knots,*tails;

{
   int i;
   struct model *mod1;
   struct datas *dat;

/* if nx<1 we are only interested in the setting of HLENGTH */
   if(*nx<1){
      *nx=HLENGTH;
      return;
   }

/* in */
   dat=makedata(*nx);
   (*dat).data=data;
   (*dat).delta=delta;
   (*dat).cc=*cc;
   mod1=makemodel();
   for(i=0;i<HLENGTH;i++) (*mod1).knots[i]=knots[i];
   (*mod1).tails[0]=tails[0];
   (*mod1).tails[1]=tails[1];
   (*mod1).tails[2]=tails[2];
   (*mod1).tails[3]=tails[3];
   (*mod1).tails[4]=tails[4];

/* do it */
   heft(dat,*nkstart,*alpha,mod1,*iauto,zerror,*nkmax,*mindist);
   if((*nkstart)< -900)return;

/* out */
   *nkmax=(*mod1).nk1;
   *nkstart=(*mod1).nk1;
   for(i=0;i<HLENGTH;i++){
      iknots[i]=(*mod1).iknots[i];
      knots[i]=(*mod1).yknots[i];
      theta[i]=(*mod1).theta[i];
      logl[i]=(*mod1).logl[i];
      ad[i]=(*mod1).ad[i];
   }
   tails[0]=(*mod1).tails[0];
   tails[1]=(*mod1).tailse[0];
   tails[2]=(*mod1).tails[2];
   tails[3]=(*mod1).tailse[1];
   tails[4]=(*mod1).tails[4];
   return;
}
static struct model *makemodel()
/* allocates storage for a model */
{
   int i;
   struct model *m1;
   m1=(struct model *)Salloc(1,struct model);
   (*m1).aic=pow(10.,100.);
   (*m1).nk=0;
   (*m1).nk1=0;
   (*m1).ll=0.;
   (*m1).tailse[0]=0.;
   (*m1).tailse[1]=0.;
   (*m1).iknots=ihvector(HLENGTH);
   for(i=0;i<HLENGTH;i++)(*m1).iknots[i]=1;
   (*m1).tails=dhvector(5);
   for(i=0;i<5;i++)(*m1).tails[i]=0.;
   (*m1).icoef=ihmatrix(HLENGTH,HLENGTH);
   (*m1).knots=dhvector(HLENGTH);
   (*m1).yknots=dhvector(HLENGTH);
   (*m1).logl=dhvector(HLENGTH);
   (*m1).ad=ihvector(HLENGTH);
   for(i=0;i<HLENGTH;i++)(*m1).ad[i]=2;
   (*m1).theta=dhvector(HLENGTH);
   (*m1).coef2=dhmatrix(HLENGTH,HLENGTH);
   (*m1).coef3=dstriparray(HLENGTH,4,HLENGTH);
   (*m1).score=dhvector(HLENGTH);
   (*m1).hessian=dhmatrix(HLENGTH,HLENGTH);
   return m1;
}
static struct datas *makedata(i)
/* allocates storage for a data-datastructure with i observations */
int i;
{
   struct datas *d1;
   d1=(struct datas *)Salloc(1,struct datas);
   /* if(!d1) hrerror("allocation error in makedata()"); */
   (*d1).nd=i;
   (*d1).cc=0.;
   (*d1).delta=ihvector(i);
   (*d1).data=dhvector(i);
   (*d1).basdata1=dhmatrix(i,HLENGTH);
   (*d1).basdata2=dhmatrix(i,HLENGTH);
   for(i=i-1;i>=0;i--)(*d1).delta[i]=1;
   return d1;
}
/******************************************************************************/
static double ***dstriparray(r,c,s)
int r,c,s;
{
   int i,j,k;
   double ***m;

   m=(double ***) Salloc(r+1,double**);
   for(i=0;i<=r;i++) {
      m[i]=(double **)Salloc(c+1,double*);
      for(j=0;j<=c;j++){
         m[i][j]=(double *)Salloc(s+1,double);
         for(k=0;k<=s;k++)m[i][j][k]=0.;
      }
   }
   return m;
}
/******************************************************************************/
static int *ihvector(l)
int l;
/* allocate an int vector with subscript range v[0...l] */
{
   int i,*v;
   v=(int *)Salloc(l+1,int);
   for(i=0;i<=l;i++)v[i]=0;
   return v;
}
/******************************************************************************/
static double *dhvector(l)
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
static int **ihmatrix(r,c)
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
static double **dhmatrix(r,c)
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

/* this is the main loop */

static void heft(dat,nkstart,alpha,mod1,iauto,zerror,nkmax,mind)

struct datas *dat;
struct model *mod1;
int nkstart,iauto,nkmax,*zerror,mind;
double alpha;

/* dat    - the data
   nkstart- starting number of knots
   alpha  - penalty parameter (bic)
   mod1   - the working model
   iauto  - 0 - fully automatic knots
            2 - user chooses knots
   zerror  - error conditions 
   nkmax  - maximum number of knots 
   mind   - minimum distance between knots */

{
   int i00;

   struct model *modmin,*modold,*makemodel();

   int nint1=20,nint2=50,nint,nintx1=2,nintx2=6,nintx=0,i,j,addi=0,nkmax2,ndd;
   double r,newk,lold[HLENGTH],*ddd;

/* modmin - model with minimum aic
   modold - old model
   nkmax2 - copy of nkmax on entrance
   nint   - number of integration points active
   nint1  - number of integration points low precision
   nint2  - number of integration points low precision
   nintx  - number of times more points before first knot and in last interval
   nintx1,nintx2 - to nintx as nint2 and nint1 are to nint
   i,j    - counter
   addi   - are we adding (0=no, 1=yes, 2=just gave up)
   r      - utility
   newk   - knot to be added
   lold   - old loglikelihoods - in case we don't want to add anymore */

/* initialize */
   for(i=0;i<HLENGTH;i++) lold[i]=-pow((double)10.,(double)99.);

/* allocate memory */
   nkmax2=nkmax;
   modmin=makemodel();
   (*modmin).aic=pow((double)10.,(double)99.);
   modold=makemodel();
   i00=(HLENGTH+15+2*nintx2)*nint2;
   (*mod1).basvec=dhvector(i00);
   (*mod1).mult=dhvector(i00);
   (*mod1).basmat=dhmatrix(i00,HLENGTH);
   allocer((*dat).nd,i00);
   ddd=wkddd;

   ndd=0;
   for(i=0;i<(*dat).nd;i++){
      if((*dat).delta[i]==1){
         ddd[ndd]=(*dat).data[i];
         ndd++;
      }
   }
/* cc is the upper quartile */
   if((*dat).cc< -90000.){
      r=(.75)*ndd-0.75;
      i=floor(r);
      r=r-i;
      (*dat).cc=(1-r)*ddd[i]+r*ddd[i+1];
   }

/* only deletion */
   if(iauto > 0){

/* place the knots */
      hknotplace(&nkstart,mod1);
      if(nkstart< -900)return;
      if(nkmax==nkstart)addi=0;
/* positioned */
      else addi=1;
   }

/* addition and deletion compute the maximum number of knots */
   if(nkmax==0){
      r = 5.*pow((double)((*dat).nd),0.2);
      if(r>29.9)r=29.9;
      if((*dat).nd<=60)r=(*dat).nd/5.;
      if(r<1.5)r=1.5;
      nkmax=ceil(r);
   }

/* place knots */
   if(iauto==0){
      r=(.25)*ndd-0.25;
      i=floor(r);
      r=r-i;
      (*mod1).knots[0]=(1-r)*ddd[i]+r*ddd[i+1];
/* if we are constant in the left tail, we start with 2 knots */
      if((*mod1).tails[4]<0.5){
         (*mod1).nk = 2;
         nkstart=2;
         r=(.75)*ndd-0.75;
         i=floor(r);
         r=r-i;
         (*mod1).knots[1]=(1-r)*ddd[i]+r*ddd[i+1];
         if((*mod1).knots[1]==(*mod1).knots[0]){
            i=hlocation(1,ddd,ndd,(*mod1).knots[0]);
            if(i==ndd-1){
               (void)Rprintf("too few distinct data: 1st quart=max\n");
               nkstart=-998;
               return;
            }
            (*mod1).knots[1]=ddd[i+1];
         }
      }
      else{
/* if we are linear in the left tail, we start with 3 knots */
         r=(.5)*ndd-0.5;
         i=floor(r);
         r=r-i;
         (*mod1).knots[1]=(1-r)*ddd[i]+r*ddd[i+1];
         r=(.75)*ndd-0.75;
         i=floor(r);
         r=r-i;
         (*mod1).knots[2]=(1-r)*ddd[i]+r*ddd[i+1];
         (*mod1).nk = 3;
         nkstart=3;
         if((*mod1).knots[1]==(*mod1).knots[0]){
            i=hlocation(0,ddd,ndd,(*mod1).knots[0]);
            if(i==0){
               (void)Rprintf("too few distinct data: median=min data\n");
               nkstart=-998;
               return;
            }
            (*mod1).knots[0]=ddd[i-1];
         }
         if((*mod1).knots[1]==(*mod1).knots[2]){
            i=hlocation(1,ddd,ndd,(*mod1).knots[1]);
            if(i==ndd-1){
               (void)Rprintf("too few distinct data: median=max data\n");
               nkstart=-998;
               return;
            }
            (*mod1).knots[2]=ddd[i+1];
         }
      }
      addi=1;
   }
   if(zerror[6]==37){
      (void)Rprintf("starting knots at ");
      for(i=0;i<nkstart;i++)(void)Rprintf("%.2f ",(*mod1).knots[i]);
      (void)Rprintf("\n");
   }

/* the knot addition loop starts here */
   if(addi==1){
      do{
         getcoef(mod1);

/* compute the integration multipliers  */
         nint=nint1;
         nintx=nintx1;
         if((*mod1).nk< -5){
            nint=nint2;
            nintx=nintx2;
         }
         intprep(&nint,nintx,mod1,dat,1); /* should get basis in it */

/* compute the starting values */
         nstart(mod1,dat,nkstart);

/* fit the model */
         hiter(mod1,dat,zerror,nint,1);

/* if zerror[1]>0 there were problems */
         if(zerror[1]>0 && (*mod1).nk<6){
            (void)Rprintf("sorry - can't recover with so few knots (%d)\n",(*mod1).nk);
         /* if(zerror[0]==0) exit(1);
            else  return; */
            return;
         }
/* if zerror[1]=2 we might be helped by starting to remove */
         if(zerror[1]==2){
            nkstart=(*mod1).nk-1;
            dubmodel(mod1,modold);
            (void)Rprintf("trying to start removing knots.....\n"); 
            addi=2;
         }

/* record the old fit, justin case */
         dubmodel(modold,mod1);

/* no reason to keep on adding */
         lold[(*mod1).nk]=(*mod1).ll;
         if(nkmax!=nkstart && nkmax2==0){
            for(i=2;i<(*mod1).nk-2;i++){
               if((*mod1).ll-lold[i]<((*mod1).nk-i-2.)/2.+0.5){
                  nkmax=(*mod1).nk;
                  addi=2;
               }
            }   
         }

/* have we added enough? */
         if(addi==1 && nkmax==(*mod1).nk) addi=2;
         hetse(mod1,modmin,alpha);

/* post-processing addition */
         if(addi==1){
            newk=add(mod1,dat,nint,zerror,modmin,mind);
/* oops, cannot add anymore */
            if(newk<0) {
               addi=2;
               nkmax=(*mod1).nk;
            }
         }
      }while(addi==1);
   }

/* record where we start from */
   for(i=0;i<HLENGTH;i++) (*mod1).yknots[i]=(*mod1).knots[i];
   for(i=0;i<HLENGTH;i++) (*modmin).yknots[i]=(*mod1).knots[i];
   (*modmin).nk1=nkmax;
   (*mod1).nk1=nkmax;

/* the knot removal loop starts here */
   do{
      getcoef(mod1);

/* compute the coefficients of the basisfunctions */
      if(addi==2){
         nint=nint2;
         nintx=nintx2;
      }
      else{
         for(i=0;i<(*dat).nd;i++){
            for(j=0;j<=(*mod1).nk+1;j++){
               (*dat).basdata2[i][j]=(*dat).basdata1[i][j];
            }
         }
      }
      intprep(&nint,nintx,mod1,dat,addi); /* basis , intprep only if addi==2 */

/* compute the starting values */
      if(addi==0) start(mod1,dat);
      addi=0;

/* fit the model */
      hiter(mod1,dat,zerror,nint,0);
      if(zerror[1]>0){
         zerror[1]=0;
         nstart(mod1,dat,(*mod1).nk);
         hiter(mod1,dat,zerror,nint,1);
         if(zerror[1]>0){
            (void)Rprintf("sorry - cannot recover during removal fase..\n");
         /* if(zerror[0]==0) exit(1);
            else  return; */
         }
      }

/* post-processing - removal */
      tossit(mod1,modmin,alpha,zerror);
   }while(((*mod1).nk>1 && (*mod1).tails[4]<0.5 )|| (*mod1).nk>2);

/* send the correct model back */
   dubmodel(mod1,modmin);
   for(i=0;i<HLENGTH;i++)(*mod1).ad[i]=(*modmin).ad[i];

/* get theta in the power basis format */
   thetaswap(mod1);
   return;
}
/******************************************************************************/
/* checks the knots */


static void hknotplace(nkstart,mod1)

int *nkstart;
struct model *mod1;

/* nkstart - starting number of knots
   mod1    - model */

{
   int i,i1,jj;
   double r;

/* i    - counter
   i1   - utility
   r    - utility */

/* check for negative knots */
   i1=0;
   jj=1;
   if((*mod1).knots[0]<=0.){
      (void)Rprintf("*** first knot <= 0 ***\n");
      jj=-999;
   }

/* check for knots out of sequence and double knots */
   if(jj>0)for(i=1;i<*nkstart;i++){
      if((*mod1).knots[i]>(*mod1).knots[i1]){
         i1++;
         (*mod1).knots[i1]=(*mod1).knots[i];
      }
      else{
         if((*mod1).knots[i]<(*mod1).knots[i1]){
            (void)Rprintf("** knots not in sequence **\n");
            jj = -999;
         }
         if((*mod1).knots[i]==(*mod1).knots[i1]){
            (void)Rprintf("*** warning, knot %d is double: removed ***\n",i);
         }
      }
   }

/* how many knots are left */
   if(jj>0){
      *nkstart = i1+1;

/* copy the knots in yknots */
      for(i=0;i<*nkstart;i++)(*mod1).yknots[i]=(*mod1).knots[i];
      r=1.;
      for(i=1;i<*nkstart;i++){
         if((*mod1).knots[i]/(*mod1).knots[i-1]>r){
            r=(*mod1).knots[i]/(*mod1).knots[i-1];
         }
      }
      if(r>4000.){
         (void)Rprintf( "*** warning: max knot-ratio is %e - answers inaccurate ***\n",r);
      }
      (*mod1).nk=*nkstart;
   }
   else{
      *nkstart= -999;
   }
}
/******************************************************************************/
/* This function computes the coefficients of the basis functions from the
knots the basis funcftions are G2-G(p-1), where (p=K+1). G2=B1.
Basis functions B(2)-B(nk-3) are multiples of B-splines.  Further the
coefficients are choosen such that the quadratic and  cubic terms in both
tails are 0;  this leads to differnt basis functions for B(1), B(nk-2)
and B(nk-1). B(1) is linear left of the first knot. B(nk-2) is constant to
the right of the last knot and B(nk-1) is constant 1 everywhere */

static void getcoef(m1)
struct model *m1;
{
   getcoefx((*m1).coef2,(*m1).coef3,(*m1).knots,(*m1).icoef,(*m1).nk);
}

static void getcoefx(coef2,coef3,knots,icoef,nk)

double **coef2,***coef3,*knots;
int nk,**icoef;

/* coef2 - first index: basis function number-1,
           second index: 0:1, 1:x, 2:(x-t1)+^3,  3:(x-t2)+^3, 4:(x-t3)+^3,.....
   coef3 - between knot(i) and knot(i+1) the coef of x^power of basisfct(j)
           first index: basis function number-1 (j-1)
           second index: power of x
           third index: interval (i)  
   icoef - does basisfunction i exist in interval j?
   knots - knots
   nk    - number of knots */

{
   int i,j,k;
   double z0,z1;

/* i j k       - counter
   z0,z1       - value of constants of two succesive basisfunctions  */

/* Initializations  */
   for(i=0; i<nk-1; i++){
      for(j=0; j<nk+2; j++){
         coef2[i][j]=0.;
         icoef[i][j]=0;
         for(k=0; k<4; k++) coef3[i][k][j]=0.;
      }
   }

/* The coefficients for the two tail basis functions are easy to compute  */
   if(nk > 2){
      coef2[0][2] = 1.;
      coef2[0][3] = (knots[0]-knots[2]) / (knots[2]-knots[1]);
      coef2[0][4] = (knots[1]-knots[0]) / (knots[2]-knots[1]);
      coef2[0][1] = -3. * (pow(knots[0],2.) + coef2[0][3] * pow(knots[1],2.)
                                            + coef2[0][4] * pow(knots[2],2.));
      coef2[0][0] = - knots[nk-1] * coef2[0][1]
                    - coef2[0][2] * pow((knots[nk-1]-knots[0]),3.)
                    - coef2[0][3] * pow((knots[nk-1]-knots[1]),3.)
                    - coef2[0][4] * pow((knots[nk-1]-knots[2]),3.);
      coef2[0][5] = 0.;
   }
   coef2[nk-2][0] = 1.;

/* we first create basis functions that are 0 before knot[i] and constant
   after knot [i+3]  */

   if(nk > 3){
      for(i=1;i<nk-2;i++){
         coef2[i][i+1] = 1.;
         coef2[i][i+4] = (knots[i+1]-knots[i-1])*(knots[i-1]-knots[i])
                         /((knots[i+1]-knots[i+2])*(knots[i]-knots[i+2]));
         coef2[i][i+3] = (coef2[i][i+4]*(knots[i]-knots[i+2])
                         +knots[i]-knots[i-1])/(knots[i+1]-knots[i]);
         coef2[i][i+2] = -1.-coef2[i][i+3]-coef2[i][i+4];
      }
   }

/* In the following part we subtract a number of times one basis
   function from another - so that basis function i becomes 0 after knot[i+4] */
   if(nk > 4){
      for(i=1;i<nk-3;i++){
         z0 = 0.;
         z1 = 0.;
         for(j=2;j<nk+1;j++){
            z0 += coef2[i][j]   * pow((knots[nk-1]-knots[j-2]),3.);
            z1 += coef2[i+1][j] * pow((knots[nk-1]-knots[j-2]),3.);
         }
         for(j=2; j<nk+2; j++) coef2[i][j] += - (z0 / z1) * coef2[i+1][j];
      }
   }

/* Now the coef3 matrix. First basis function 1. */
   if(nk>2){
      for(k=0; k<3; k++){
         coef3[0][1][k]=coef2[0][1];
         coef3[0][0][k]=coef2[0][0];
         icoef[0][k]=1;
      }

/* The rest is a bit tricking with the correct indices */
      for(i=0;i<nk-2;i++){
         for(j=i;j<i+4;j++){
            for(k=i+1;k<j+2;k++){
               if(j>0 && j<nk+1 && (i!=0 || j!=3)){
                  if(k != 1){
                     coef3[i][0][j] += -coef2[i][k]*pow(knots[k-2],3.);
                     coef3[i][1][j] += 3.*coef2[i][k]*pow(knots[k-2],2.);
                     coef3[i][2][j] += -3.*coef2[i][k]*knots[k-2];
                     coef3[i][3][j] += coef2[i][k];
                     icoef[i][j]=1;
                  }
               }
            }
         }
      }
   }

/* initialize the constant basis */
   for(j=0;j<nk+1;j++){
      coef3[nk-2][0][j]=1.;
      icoef[nk-2][j]=1;
   }
}

/******************************************************************************/

/* computes the starting values removal stage - L2 projection on a smaller
   space  */

static void start(mod1,dat)
struct model *mod1;
struct datas *dat;
{
   hstart2((*mod1).theta,(*mod1).nk,(*dat).basdata1,(*dat).basdata2,
          (*dat).nd,(*mod1).tails);
}

static void hstart2(theta,nk,basdata1,basdata2,nx,tails)

int nk,nx;
double *theta,**basdata1,**basdata2,*tails;

/* theta     - theta 
   nk        - present number of knots
   basdata1  - present basis functions in datapoints
   basdata2  - previous basis functions in datapoints
   nx        - number of data points 
   tails     - which tail basis functions are included ? */

{
   int i,j,k;
   double **mat1,*vec2,*vec1;

/* i,j,k   - counter
   mat1    - X matrix
   vec2    - Y 
   vec1    - XtY */

   if(tails[0]>0) theta[0]=tails[1];
   if(tails[2]>0) theta[nk]=tails[3];

/* first allocate some storage space */
   mat1=wkmat1;
   vec1=wkvec1;
   vec2=wkvec2;

/* compute the fitted values = Y */
   for(i=0;i<nx;i++){
      vec2[i]=0.;
      for(j=1;j<=nk;j++) vec2[i]+=theta[j]*basdata2[i][j];
   }

/* compute XtX */
   for(i=1;i<nk;i++){
      for(j=i;j<nk;j++){
         mat1[i-1][j-1]=0.;
         for(k=0;k<nx;k++) mat1[i-1][j-1] += basdata1[k][i]*basdata1[k][j];
      }
   }

/* make XtX symmetric */
   for(i=2;i<nk;i++){
      for(j=1;j<i;j++) mat1[i-1][j-1]=mat1[j-1][i-1];
   } 

/* Compute XtY */
   for(i=1;i<nk;i++){
      vec1[i-1]=0.;
      for(k=0;k<nx;k++) vec1[i-1] += basdata1[k][i]*vec2[k];
   }

/* if there is no linear term */
   if(tails[4]>0.5){
      for(i=1;i<nk;i++){
         mat1[0][i]=0.;
         mat1[i][0]=0.;
      }
      mat1[0][0]=1.;
      vec1[0]=0.;
   }

/* solve the system */
   i=0;
   hlusolve2(mat1,nk-1,vec1,&i);
   for(i=1;i<nk;i++) theta[i]=vec1[i-1];
   theta[nk]=theta[nk+1];
   for(i=0;i<nx;i++){
      vec2[i]=0.;
      for(j=1;j<nk;j++) vec2[i]+=theta[j]*basdata1[i][j];
   }
}

/******************************************************************************/

/* starting values after a knot was added - L2 projection on a larger space */

static void nstart(mod1,dat,nkstart)
struct model *mod1;
struct datas *dat;
int nkstart;

/* mod1  - the model
   dat   - the data
   nkstart starting number of knots */

{
   double **mat,*phi,r;
   int i,j,nk=(*mod1).nk;

/* mat    - relates coefficiets in one basis with the other basis
   phi    - first old later new values
   i,j    - counter
   nk     - see above */

/* if this is the first time we compute a one parameter estimate */
   if(nkstart==nk){
      for(i=0;i<=nk;i++)(*mod1).theta[i]=0;
      if((*mod1).tails[0]>0) (*mod1).theta[0]=(*mod1).tails[1];
      if((*mod1).tails[2]>0) (*mod1).theta[nk]=(*mod1).tails[3];
      r=0;
      j=0;
      for(i=0;i<(*dat).nd;i++){
         r+=(*dat).data[i];
         j+=(*dat).delta[i];
      }
      r = r/(double)j;
      (*mod1).theta[nk-1]= -hmylog(r);
      return;
   }

/* things are on a power basis and should get to the real basis */
   phi=wkphi;
   mat=wkmat;
   for(j=0; j<nk+2; j++){
      for(i=0;i<nk-1;i++) mat[j][i]=(*mod1).coef2[i][j];
   }
   for(j=0;j<nk+2;j++) phi[j]=(*mod1).theta[j];
   j=0;
   hlusolve2(mat,nk-1,phi,&j);
   (*mod1).theta[0]=(*mod1).theta[nk+2];
   (*mod1).theta[nk]=(*mod1).theta[nk+3];
   for(j=0;j<nk-1;j++) (*mod1).theta[j+1]=phi[j];
}
/******************************************************************************/

/* This function changes theta from the basisfunction representation into the
   truncated power basis representation */

static void thetaswap(mod1)
struct model *mod1;

{
   double *phi;
   int i,j;

/* phi - theta for (part of the) power basis */

   phi=wkphi2;
   for(i=0;i<HLENGTH;i++) phi[i]=0;
   for(j=0;j<(*mod1).nk-1;j++){
      for(i=0;i<(*mod1).nk+2;i++)phi[i]+=(*mod1).theta[j+1]*(*mod1).coef2[j][i];
   }
   (*mod1).theta[(*mod1).nk1+3]=(*mod1).theta[(*mod1).nk];
   (*mod1).theta[(*mod1).nk1+2]=(*mod1).theta[0];
   (*mod1).theta[(*mod1).nk1+1]=phi[1];
   (*mod1).theta[(*mod1).nk1]=phi[0];
   j=1;
   for(i=0;i<(*mod1).nk1;i++){
      (*mod1).theta[i]=0.;
      if((*mod1).iknots[i]==1){
         j++;
         (*mod1).theta[i]=phi[j];
      }
   }
}
/******************************************************************************/
/* these routines compute the multipliers for the integrals */


static void intprep(nint,nintx,mod1,dat,what)

int *nint,nintx,what;
struct model *mod1;
struct datas *dat;

/* nint    - number of integration points per knot-interval
   nintx   - extra integration points before first knot and in last interval
   mod1    - the model
   dat     - the data
   what    - if 0 prepare for high-precision
             if 1 prepare for low precision
             if 2 do only basis */
{
   int i,j=0,k=0,where=0,npart=*nint,nkplus=0,nx=(*dat).nd,hopla=1;
   double *masterpt,*xx,cr;

/* xx - data[0],data[nx-1] andd some other quantiles
   nkplus  - nr of distinct elements in (knots,xx)
   masterpt- distinct elements of (knots,data[0],xx)
   cr      - we don't need extra points closer together than this.
   npart   - nint per part
   i,j     - counter
   where   - coefficient of the next datapoint to be covered 
   hopla   - keep track of length of xx
   nx      - (*dat).nd
   k       - was the last point of masterpt a knot or a xx? */

   if(what!=0){
      masterpt=wkmasterpt;
      xx=wkxx;

/* select datapoints that are potential masterpoints */
      xx[0]=(*dat).data[nx/100];
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/25));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/15));
      hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/10));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/8));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/6));
      hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/4));
      hopla=hopplus(xx,(*dat).data,hopla,(int)(nx/2));
      hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-nx/4));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-nx/6));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-nx/8));
      hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-nx/10));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-nx/15));
      if(what==2) hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-nx/25));
      hopla=hopplus(xx,(*dat).data,hopla,(int)(nx-1));
   
/* points shouldn't be too close together, cr is the minimum distance */
      cr=(xx[1]-xx[0])/1.5;
      for(i=1;i<hopla-1;i++) if((xx[i+1]-xx[i])/1.5<cr) cr=(xx[i+1]-xx[i])/1.5;
      i=0;

/* merge the knots and the xx into masterpt */
      do{
         if(xx[i]<(*mod1).knots[j] || j==(*mod1).nk){
            masterpt[nkplus]=xx[i];
            nkplus++;
            if(k==1 && masterpt[nkplus-1]-masterpt[nkplus-2]<cr){
               if(i<hopla-1 || (*mod1).knots[(*mod1).nk-1]>=(*dat).data[nx-3]){
                  nkplus--;
               }
               else k=2;
            }
            else k=2;
            i++;
         }
         else{
            masterpt[nkplus]=(*mod1).knots[j];
            nkplus++;
            if(k==2 && masterpt[nkplus-1]-masterpt[nkplus-2]<cr){
               if(j>0 || i>0 || (*mod1).knots[0]<(*dat).data[2] || i<1){
                  nkplus--;
                  masterpt[nkplus-1]=masterpt[nkplus];
               }
            }
            j++;
            k=1;
         }
      }while(i<hopla || j<(*mod1).nk);
      masterpt[nkplus-1]=(*dat).data[nx-1];
      if(masterpt[0]==(double)0) masterpt[0]=masterpt[1]/2.;
      *nint=(nkplus+2*nintx)*(*nint)+1; 

/* first compute basvec */
      (*mod1).basvec[0]=1./(21.365*(double)((nintx+1)*npart))*masterpt[0];
      (*mod1).basvec[1]=1./(2.0904*(double)((nintx+1)*npart))*masterpt[0];
      for(i=2;i<=npart*(nintx+1);i++){
         (*mod1).basvec[i]=
            ((double)(i-1)/(double)((nintx+1)*npart-1))*masterpt[0];
      }
      if(nkplus>2){
         for(j=1;j<nkplus-1;j++){
            for(i=1;i<=npart;i++){
               (*mod1).basvec[(j+nintx)*npart+i]=((double)i/(double)npart)*
                             (masterpt[j]-masterpt[j-1])+masterpt[j-1];
            }
         }
      }
      for(i=1;i<=npart*(nintx+1);i++){
         (*mod1).basvec[i+(nkplus+nintx-1)*npart]=
            ((double)i/(double)((nintx+1)*npart))
                    *(masterpt[nkplus-1]-masterpt[nkplus-2])+masterpt[nkplus-2];
      }

/* initialize the multipliers */
      for(i=0;i<(npart*(nkplus+2*nintx)+1);i++) (*mod1).mult[i]=0.;

/* integrate per interval in between two integration points */
      for(j=0;j<(nkplus+2*nintx);j++){
         for(i=0;i<npart;i++) midblob(&where,j,i,npart,
            (*mod1).basvec,(*mod1).mult,nx,(*dat).data);
      }
   }
   basis((*dat).data,(*dat).nd,(*mod1).knots,(*mod1).nk,(*dat).basdata1,
            (*dat).cc,(*mod1).icoef,(*mod1).coef3);
   basis((*mod1).basvec,*nint,(*mod1).knots,(*mod1).nk,(*mod1).basmat,
            (*dat).cc,(*mod1).icoef,(*mod1).coef3);
}
/******************************************************************************/

/* this one integrates between two integration points */

static void midblob(where,j,i,nint,basvec,mult,nx,data)

int nint,j,i,*where,nx;
double *basvec,*mult,*data;

/* see all above */

{
   double x1,x2,x3,x4,zmin,zmax;
   int i1,i2,i3,i4;

/* x1,x2,x3,x4 - the four points that are used to interpolate
   i1,i2,i3,i4 - their indices
   zmin,zmax   - boundaries of the part of the integral that still has to be
                 done  */

/* set zmin and zmax */
   zmin=basvec[j*nint+i];
   zmax=basvec[j*nint+i+1];
   if(j*nint+i==0)zmin=0.;

/* select the interpolation points */
   if(i==0)i=1;
   if(i==nint-1)i=nint-2;
   i1=j*nint+i-1;
   i2=j*nint+i;
   i3=j*nint+i+1;
   i4=j*nint+i+2;
   x1=basvec[i1];
   x2=basvec[i2];
   x3=basvec[i3];
   x4=basvec[i4];

/* figure out what the next upper limit is, and integrate until there */
   do{
      if(data[*where]>=zmax || *where >= nx){
         lgrange(x1,x2,x3,x4,(double)(nx-*where),zmin,zmax,mult,i4);
         lgrange(x1,x2,x4,x3,(double)(nx-*where),zmin,zmax,mult,i3);
         lgrange(x1,x3,x4,x2,(double)(nx-*where),zmin,zmax,mult,i2);
         lgrange(x2,x3,x4,x1,(double)(nx-*where),zmin,zmax,mult,i1);
         zmin=zmax+1.; 
         if(data[*where]==zmax) (*where)++;
      }
      else{
         lgrange(x1,x2,x3,x4,(double)(nx-*where),zmin,data[*where],mult,i4);
         lgrange(x1,x2,x4,x3,(double)(nx-*where),zmin,data[*where],mult,i3);
         lgrange(x1,x3,x4,x2,(double)(nx-*where),zmin,data[*where],mult,i2);
         lgrange(x2,x3,x4,x1,(double)(nx-*where),zmin,data[*where],mult,i1);
         zmin=data[*where];
         (*where)++;
      }
   }while(zmin<zmax);
}

/******************************************************************************/
/* this computes the integral shown below, for q=0,1 and 2. It adds the
   results to row j of mult */

static void lgrange(a,b,c,d,n,u,v,mult,j)

double a,b,c,d,n,u,v,*mult;
int j;

/*  a,b,c,d,n,u,v - see below
    j             - row of mult to which the integral should be added
    mult          - multipliers

       v
      /
      | n*(x-a)*(x-b)*(x-c)        q
      | -------------------(log(x))  dx
      |  (d-a)*(d-b)*(d-c)
      /
     u                                     */
               
{
   double cc[5],uu[6],vv[6];
   int i;

/* prepare the coef */
   n=n/((d-a)*(d-b)*(d-c));
   cc[4]= n;
   cc[3]= -n*(a+b+c);
   cc[2]= n*(a*b+a*c+b*c);
   cc[1]= -a*b*c*n;

/* prepare the lower bound */
   uu[1]=u;
   for(i=2;i<5;i++) uu[i]=u*uu[i-1]*(double)(i-1)/(double)i;

/* prepare the upper bound */
   vv[1]=v;
   for(i=2;i<5;i++) vv[i]=v*vv[i-1]*(double)(i-1)/(double)i;

/* compute the integrals, note that the second index is q */
   for(i=1;i<5;i++) mult[j]-=(vv[i]-uu[i])*cc[i];
}     

/******************************************************************************/
/* silly, but needed in intprep */
static int hopplus(xx,data,i,j)
double *xx,*data;
int i,j;
{
   xx[i]=data[j];
   if(xx[i]>xx[i-1])return i+1;
   return i;
}
/******************************************************************************/
static void basis(x,nx,knots,nk,basmat,cc,icoef,coef3)

double *x,*knots,**basmat,cc,***coef3;
int nx,nk,**icoef;

/* x      - sorted vector of data points, in which the basisfunctions are to be
            computed
   nx     - length of x
   knots  - vector of knots
   nk     - length of knots
   basmat - to be the basisfunctions in each point
   cc     - the cc number for the log-terms.  
   icoef  - does basis function [i] exist in interval t(j-1)-t(j), it is in
            icoef[i][j];
   coef3  - coefficient of x^j for basis function [i] in interval t(k-1)-t(k).*/

{
   int i,j,where=0;
/* where indicates in between which two knots a point is */

/* inialize */
   for(i=0;i<nx;i++){
      for(j=1;j<nk;j++) basmat[i][j]=0.;
   }

   for(i=0;i<nx;i++){

/* the two log basis functions */
      if(x[i]>0) basmat[i][0]=hmylog(x[i]/(x[i]+cc));
      basmat[i][nk]=hmylog(x[i]+cc);

/* find where the knot is */
      if(x[i]>knots[where] && where<nk){
         do{
            where++;
         } while(x[i]>knots[where] && where<nk);
      }
      basmat[i][nk+1]=0.;
      basmat[i][nk+2]=0.;
      for(j=1;j<nk-1;j++){
         if(basmat[i][nk+1]<0.5 && icoef[j-1][where]!=0){
            basmat[i][nk+1]=j;
            j=nk+10;
         }
      }
      for(j=nk-2;j>0;j--){
         if(basmat[i][nk+2]<0.5 && icoef[j-1][where]!=0){
            basmat[i][nk+2]=j;
            j=0;
         }
      }

/* update the other basis functions */
      for(j=1;j<nk;j++){
         if(icoef[j-1][where]!=0){
            basmat[i][j]=((coef3[j-1][3][where]*x[i]+coef3[j-1][2][where])
                                               *x[i]+coef3[j-1][1][where])
                                               *x[i]+coef3[j-1][0][where];
         }
      }
   }
}
/******************************************************************************/
/* the main Newton Raphson loop */ 


static void hiter(mod1,dat,zerror,nint,what)

struct model *mod1;
struct datas *dat;
int *zerror,nint,what;

/* mod1   - model
   dat    - data
   nint   - number of integration points
   zerror  - zerror conditions 
   what   - 0 if we are deleting and could have another shot at the starting
            values, 1 else */

{
   double ldif=0.;
   int i,j,ctr,itails[3],status;


/* ldif   - lnew - lold
   i,j,k  - counter
   ctr    - iteration counter */

/* status
   0
   1
   2   - left tail allert (itails[0]=2)
   3
   4   - right tail allert (itails[2]=2) */
   
/* itails[0]:
  -1 - converged with itails[0]==2, let's now try....
   0 - left log term included and operational
   1 - left log term not included (or user fixed)
   2 - left log term protection against -1 */
   if((*mod1).tails[0]>0.5) itails[0]=1;
   else itails[0]=0;

/* itails[1]:
   0 - linear left term included
   1 - linear left term not included */
   if((*mod1).tails[4]>0.5) itails[1]=1;
   else itails[1]=0;

/* itails[2]:
   0 - right term included and operational
   1 - right term not included (or user fixed)
   2 - right term protection against -1 */
   if((*mod1).tails[2]>0.5) itails[2]=1;
   else{
      itails[2]=0;
      if((*mod1).theta[(*mod1).nk]< -0.999) itails[2]=2;
   }

/* iterations start */
   for(ctr=1;ctr<500;ctr++){

/* one step */
      status=step(dat,mod1,itails,&ldif,nint,zerror,what);

/* problems in the right tail */
      if((*mod1).theta[(*mod1).nk]<-1){
         if(zerror[0]==0){
            warning("*** warning: right tail adjustment ***\n");
         }
         (*mod1).theta[(*mod1).nk]=-1;
         itails[2]=2;
         status=4;
      }

/* serious problems */
      if(status==1 || status==3 || (status==2 && itails[0]==-1)){
         zerror[1]=1;
         return;
      }

/* problems in the left tail */
      if(status==2 && itails[0]==0){
         (*mod1).theta[0]=-0.8;
         itails[0]=2;
      }

/* is there convergence (or pseudo convergence) */
      if(status==0 && ldif<0.0000001){

/* we are done */
         if(itails[0]<2 && itails[2]<2) ctr+=10000;

/* we might be done, we have converged in a subspace */
         if(itails[2]==2){
            ldif=summer(mod1,2,nint,dat);
            if((*mod1).score[(*mod1).nk]>0.) itails[2]=0;
            else ctr+=10000; 
         }

/* we have converged in a subspace and take it from there */
         if(itails[0]==2) itails[0]=-1;
      }
   }

/* if ctr<1000 there was no convergence */
   if(ctr<1000){
      zerror[1]=2;
      (void)Rprintf("*** zerror: no convergence ***\n");
      return;
   }

/* final bookkeeping */
   if((*mod1).tails[2]<2.5) ldif=summer(mod1,2,nint,dat);

/* adjust for fixed tail thetas */
   for(i=0;i<3;i++){
      if(itails[i]!=0){
         if(i==2) i=(*mod1).nk;
         for(j=0;j<=(*mod1).nk;j++){
            (*mod1).hessian[j][i]=0;
            (*mod1).hessian[i][j]=0;
         }
         (*mod1).score[i]=0.;
         (*mod1).hessian[i][i]=-1.;
      }
   }

   if(zerror[6]==37 || zerror[0]==0){
      (void)Rprintf("logl= %.2f ",ldif);
      (void)Rprintf("(nk = %d)\n",(*mod1).nk);
   }
     
   (*mod1).ll=ldif;
   return;
}
/******************************************************************************/
/* computes l(), S() and H() */

static double summer(mod1,what,nint,dat)
struct model *mod1;
int nint,what;
struct datas *dat;
{
   return summer2((*mod1).score,(*mod1).hessian,what,(*mod1).nk,(*dat).nd,nint,
        (*mod1).theta,(*dat).basdata1,(*mod1).basmat,(*dat).delta,(*mod1).mult);
}

static double summer2(score,hessian,what,nk,ndata,nint,theta,basdata,basint,delta,mult)

double *score,**hessian,*theta,**basdata,**basint,*mult;
int what,nk,ndata,nint,*delta;

/* score   - score function
   hessian - hessian
   what    - 0: just logl, 1: also score, 2: also hessian;
   nk      - number of knots
   ndata   - number of datapoints
   nint    - number of integration points
   theta   - theta (see above)
   basdata - basisfunctions in datapoints
   basint  - basisfunctions in integration points
   delta   - delta for data points
   mult    - multipliers for p=0/1/2 in integration points */

{
   double logl=0.,lam,lm0,lm1;

/* logl - loglikelihood
   lam  - lambda or exp(lambda) 
   i,j,k- counters */

/* initializations */
   int i,j,k;
   if(what>0){
      for(i=0;i<=nk;i++){
         score[i]=0.;
         if(what>1) for(j=0;j<=nk;j++) hessian[i][j]=0.;
      }
   }

/* the integral part, anything related to basisfunction 1 goes different */
   for(i=0;i<nint;i++){
      lam = exp(lambda(nk,basint,theta,i));
      lm0 = lam*mult[i];
      logl += lm0;
      if(what >0){
         score[0] += lm0*basint[i][0];
         score[nk-1] += lm0*basint[i][nk-1];
         score[nk] += lm0*basint[i][nk];
         for(j=(int)basint[i][nk+1];j<=(int)basint[i][nk+2] && j>0;j++){ 
            score[j] += lm0*basint[i][j];
         }
         if(what >1){
            lm1=lm0*basint[i][nk];
            for(k=0;k<=nk;k++) hessian[k][nk] += lm1*basint[i][k];
            lm1=lm0*basint[i][nk-1];
            for(k=0;k<=nk-1;k++) hessian[k][nk-1] += lm1*basint[i][k];
            lm1=lm0*basint[i][0];
            hessian[0][0] += lm1*basint[i][0];
            for(j=(int)basint[i][nk+1];j<=(int)basint[i][nk+2] && j>0;j++){ 
               lm1=lm0*basint[i][j];
               for(k=0;k<=j;k++) hessian[k][j] += lm1*basint[i][k];
            }
         }
      }
   }

/* symmatrize the hessian */
   for(j=0;j<nk;j++) for(k=j+1;k<=nk;k++) hessian[k][j] = hessian[j][k];

/* the delta - data part */
   for(i=0;i<ndata;i++){
      if(delta[i]==1){
         lam = lambda(nk,basdata,theta,i);
         logl += lam;
         if(what >0) for(j=0;j<=nk;j++) score[j] += basdata[i][j];
      }
  }
  return logl;
}
/******************************************************************************/

/* this routine computes lambda(theta) */

static double lambda(nk,basis,theta,which)

double **basis,*theta;
int nk,which;

/* nk        - number of knots
   theta[k]  - theta of B(k), (for k=1....k-1)
   theta[0]  - theta of G(1)
   theta[nk] - theta of G(p) 
   which     - see next line
   basis     - matrix with in position [which][i] basisfunction i in which */

{
   int k;
   double r=0;
   r=theta[0]*basis[which][0]+theta[nk]*basis[which][nk]
        +theta[nk-1]*basis[which][nk-1];
   for(k=(int)basis[which][nk+1];k<=(int)basis[which][nk+2] && k>0;k++){ 
       r += theta[k]*basis[which][k];
   }
   return r;
}
/******************************************************************************/
/* this routine does one Newton Raphson step */

static int step(dat,mod1,itails,ldif,nint,zerror,what)
struct model *mod1;
struct datas *dat;
int *itails,nint,*zerror,what;
double *ldif;
{
   return step2((*dat).nd,(*dat).delta,nint,(*mod1).basmat,(*mod1).mult,
          (*mod1).theta,(*mod1).nk,(*dat).basdata1,(*mod1).hessian,zerror,
          (*mod1).score,itails,ldif,what);
}

static int step2(nx,delta,nint,basmat,mult,theta,nk,basdata,hessian,zerror,score,
         itails,ldif,what)
int nx,*delta,nint,nk,*zerror,*itails,what;
double **basmat,*mult,*theta,**basdata,**hessian,*score,*ldif;

/* nx     - sample size
   delta  - censoring (0=yes, 1=no)
   itails - status of the three tail thetas
   ldif   - returns the difference between the likelihoods
   nint   - number of integration points during first part of iteration
   basmat - basis functions in integration points
   mult   - integration multipliers in integration points
   theta  - theta
   score  - score function
   nk     - present number of knots
   basdata- basis functions in datapoints
   hessian- hessian of present solution
   zerror  - zerror conditions */

{
   double *cand,lnew=0.,r,lold;
   int i,j;

/* i,k    - counter
   lold   - old log-likelihood
   lnew   - new log-likelihood
   cand   - candidate for theta
   r      - utility */

/* allocate memmory */
   cand=wkcand;
/* compute likelihood, score and hessian */
   lold=summer2(score,hessian,2,nk,nx,nint,theta,basdata,basmat,delta,mult);

/* fix some things if thetas are fixed */
   for(i=0;i<3;i++){
      if(itails[i]>0){
         if(i==2) i=nk;
         for(j=0;j<=nk;j++){
            hessian[j][i]=0;
            hessian[i][j]=0;
         }
         score[i]=0.;
         hessian[i][i]=1.;
      }
   }

/* solve the system */
   i=1;
   hlusolve(hessian,nk+1,score,&i);
   if(i==-1){
      if(what==1)(void)Rprintf("*** oops, an unstable system ***\n");
      return 1;
   }

/* if the left theta is free, it shouldn't become smaller than -1 */
   if(itails[0]<=0){
      r= -theta[0]-1.;
      if(r>-score[0]){
         r=1./pow(1.5,ceil(hmylog(-score[0]/r)/hmylog(1.5)));
         if(zerror[0]==0){
            warning("*** warning: step (-1) halving(%e) ***\n",r);
         }
         if(r<0.0001 && itails[0]>=0){
            if(zerror[0]==0){
               if(what==1)(void)Rprintf("*** warning: too much step halving ***\n");
            }
            return 2.;
         }
         for(i=0;i<=nk;i++)score[i]=score[i]*r;
      }
   }

/* step halving to increase loglikelihood */
   r=2.;
   i= -1;
   do{
      r=r/2.;
      i++;
      for(j=0;j<=nk;j++) cand[j]=theta[j]-r*score[j];
      if(zerror[0]==0 && i>0){
         warning("*** warning: step (ll) halving(%e,%e)***\n",lold,lnew);
      }
      if(r<0.000000001){
         if(what==1)(void)Rprintf("*** warning: too much step halving ***\n");
         return 3;
      }
      lnew=summer2(score,hessian,0,nk,nx,nint,cand,basdata,basmat,delta,mult);
      *ldif=lnew-lold;
   }while(*ldif<-0.00000001 && r>0);

/* record the solution */
   if(r>0){
      if(cand[nk]<-1.02){
         r=(-1.02-theta[nk])/(cand[nk]-theta[nk]);
         for(j=0;j<=nk;j++) theta[j]=r*cand[j]+(1.-r)*theta[j];
      }
      else{
         for(j=0;j<=nk;j++) theta[j]=cand[j];
      }
   }
   return 0;
}
/******************************************************************************/
/* this routine does the post-processing in the case of knot removal */

static void tossit(mod1,modmin,alpha,zerror)

struct model *mod1,*modmin;
double alpha;
int *zerror;
       
/* mod1   - present model
   modmin - minimum aic model
   alpha  - alpha (AIC)
   zerror  - zerror status */

{

/* record things like loglikelihood, check whether we improved */
   (*mod1).aic= -2.*(*mod1).ll + alpha * ((*mod1).nk+1);
   if((*mod1).tails[4]>0.5) (*mod1).aic-=alpha;
   if((*mod1).ll>(*mod1).logl[(*mod1).nk] || (*modmin).ad[(*mod1).nk]==2){
      (*mod1).logl[(*mod1).nk]=(*mod1).ll;
      (*modmin).logl[(*mod1).nk]=(*mod1).ll;
      (*modmin).ad[(*mod1).nk]=1;
   }

/* did we improve */
   if((*mod1).aic <= (*modmin).aic) dubmodel(modmin,mod1);
   else if( -2.*(*mod1).ll + alpha > (*modmin).aic) (*mod1).nk=0;

/* figure out which knot to remove (and update knots and iknots) */
   hremoveknot(mod1,zerror);
   if((*modmin).nk == (*mod1).nk+1 &&(*modmin).ad[(*mod1).nk]==1){
      (*modmin).tailse[0]=(*mod1).tailse[0];
      (*modmin).tailse[1]=(*mod1).tailse[1];
   } 
}

/******************************************************************************/
/* selects which knot to remove */


static void hremoveknot(mod1,zerror)

struct model *mod1;
int *zerror;

/* mod1 - model
   zerror- zerror status */
{
   double ratmax=0.,*se,*phi;
   int i,j,k,irmax=1,nk;

/* i j k      - counters
   phi        - linear combination of thetas
   se         - standard zerrors of phi
   ratmax     - maximum ratio se/phi
   irmax      - index of maximum ratio   
   nk         - (*mod1).nk */

/* allocate storage */
   se=wkse3;
   phi=wkphi3;

   ((*mod1).nk) += -1; 
   nk=(*mod1).nk;

/* Take linear combinations of theta such that phi is theta(phi) for
   the truncated power basis. (Which is not a basis.) */
   for(i=0;i<=nk;i++){
      phi[i] = 0.;
      for(j=0;j<nk;j++) phi[i] +=(*mod1).theta[j+1]*(*mod1).coef2[j][i+2];
      phi[i]=fabs(phi[i]);
   }

/* in case there is no left log term */
   if((*mod1).tails[0]>0.5){
      (*mod1).hessian[0][0]=-1.;
      for(j=1;j<nk+2;j++){
         (*mod1).hessian[0][j]=0.;
         (*mod1).hessian[j][0]=0.;
      }
   }
/* in case there is no right log term */
   if((*mod1).tails[2]>0.5 || (*mod1).theta[nk+1]<= -0.999999){
      for(j=0;j<nk+2;j++){
         (*mod1).hessian[nk+1][j]=0.;
         (*mod1).hessian[j][nk+1]=0.;
      }
      (*mod1).hessian[nk+1][nk+1]=-1.;
   }
/* in case there is no left-linear term */
   if((*mod1).tails[4]>0.5){
      for(j=0;j<nk+2;j++){
         (*mod1).hessian[1][j]=0.;
         (*mod1).hessian[j][1]=0.;
      }
      (*mod1).hessian[1][1]=-1.;
   }


/* Invert the information matrix, giving the covariance matrix for theta */
   hluinverse((*mod1).hessian,(int)(nk+2));

/* the Standard errors for the tail things */
   if((*mod1).tails[0]>0.5) (*mod1).tailse[0]=0.;
   else  (*mod1).tailse[0]=sqrt(-(*mod1).hessian[0][0]);
   if((*mod1).tails[2]>0.5 || (*mod1).theta[nk+1]<= -1.) (*mod1).tailse[1]=0.;
   else  (*mod1).tailse[1]=sqrt(-(*mod1).hessian[nk+1][nk+1]);

/* we are done */
   if(nk==1 || (nk==2 && (*mod1).tails[4]>0.5)) return;

/* in case there is no left-linear term */
   if((*mod1).tails[4]>0.5){
      for(j=0;j<nk+2;j++){
         (*mod1).hessian[1][j]=0.;
         (*mod1).hessian[j][1]=0.;
      }
      (*mod1).hessian[1][1]=0.;
   }

/* Take linear combinations, to get the standard errors of phi      */
   if(nk>3 || (nk==2 && (*mod1).tails[4]<0.5)){
      for(i=0;i<nk+1;i++){
         se[i] = 0.;
         for(j=0;j<nk;j++){
            for(k=0;k<nk;k++){
               se[i]-=(*mod1).coef2[j][i+2]*(*mod1).coef2[k][i+2] 
                                           *(*mod1).hessian[j+1][k+1];
            }
         }
/* not really correct, but it saves numerical trouble */
         se[i] = sqrt(fabs(se[i]));
   
/* Select for which knot se/phi takes it maximal value */
         if(se[i] > phi[i] * ratmax){
            ratmax = se[i] / phi[i];
            irmax = i;
         }
      }
   }
   else irmax=1;

/* update iknots */
   j=0;
   for(i=0;i<HLENGTH;i++){
      if((*mod1).iknots[i]==1){
         if(j==irmax){
            (*mod1).iknots[i]=0;
            i=HLENGTH;
         }
         j++;
      }
   }
   if(zerror[6]==37 && ratmax!=0){
      (void)Rprintf("knot at %.2f removed ", (*mod1).knots[irmax]);
      if(ratmax!=0) (void)Rprintf("(wald = %.2f) || ",1./(ratmax*ratmax));
   }

/* update knots */
   if(irmax<nk){
      for(i=irmax;i<nk;i++)(*mod1).knots[i]=(*mod1).knots[i+1];
   }
   return;
}
/******************************************************************************/
/* this routine checks whether the model is better , gets the SE's for the log
   terms */

static void hetse(mod1,modmin,alpha)

struct model *mod1,*modmin;
double alpha;
       
/* mod1   - present model
   modmin - minimum aic model
   alpha  - alpha (AIC) */

{

   (*mod1).aic= -2.*(*mod1).ll + alpha * ((*mod1).nk+1);
   if((*mod1).tails[4]>0.5) (*mod1).aic-=alpha;
   (*mod1).logl[(*mod1).nk]=(*mod1).ll;
   (*modmin).ad[(*mod1).nk]=0; 
   (*modmin).logl[(*mod1).nk]=(*mod1).ll;

/* did we improve */
   if((*mod1).aic <= (*modmin).aic){
      dubmodel(modmin,mod1); 
/* get the se */
      getse2(mod1);
      (*modmin).tailse[0]=(*mod1).tailse[0];
      (*modmin).tailse[1]=(*mod1).tailse[1];
   } 
}

/******************************************************************************/
/* finds the SEs */


static void getse2(mod1)

struct model *mod1;

/* mod1 - model
   zerror- zerror status */
{
   double *phi,**hh;
   int i,j,nk;

/* i j k      - counters
   phi        - linear combination of thetas
   se         - standard zerrors of phi
   ratmax     - maximum ratio se/phi
   irmax      - index of maximum ratio   
   nk         - (*mod1).nk */

/* allocate storage */
   phi=wkphi4;
   hh=wkhh;

   nk=(*mod1).nk-1;
   for(j=0;j<HLENGTH;j++){
      for(i=0;i<HLENGTH;i++) hh[i][j]=(*mod1).hessian[i][j];
   }

/* Take linear combinations of theta such that phi is theta(phi) for
   the truncated power basis. (Which is not a basis.) */
   for(i=0;i<=nk;i++){
      phi[i] = 0.;
      for(j=0;j<nk;j++) phi[i] +=(*mod1).theta[j+1]*(*mod1).coef2[j][i+2];
      phi[i]=fabs(phi[i]);
   }

/* in case there is no left log term */
   if((*mod1).tails[0]>0.5){
      hh[0][0]=-1.;
      for(j=1;j<nk+2;j++){
         hh[0][j]=0.;
         hh[j][0]=0.;
      }
   }
/* in case there is no right log term */
   if((*mod1).tails[2]>0.5 || (*mod1).theta[nk+1]<= -0.999999){
      for(j=0;j<nk+2;j++){
         hh[nk+1][j]=0.;
         hh[j][nk+1]=0.;
      }
      hh[nk+1][nk+1]=-1.;
   }
/* in case there is no left-linear term */
   if((*mod1).tails[4]>0.5){
      for(j=0;j<nk+2;j++){
         hh[1][j]=0.;
         hh[j][1]=0.;
      }
      hh[1][1]=-1.;
   }


/* Invert the information matrix, giving the covariance matrix for theta */
   hluinverse(hh,(int)(nk+2));

/* the Standard errors for the tail things */
   if((*mod1).tails[0]>0.5) (*mod1).tailse[0]=0.;
   else  (*mod1).tailse[0]=sqrt(-hh[0][0]);
   if((*mod1).tails[2]>0.5 || (*mod1).theta[nk+1]<= -1.) (*mod1).tailse[1]=0.;
   else  (*mod1).tailse[1]=sqrt(-hh[nk+1][nk+1]);

   return;
}
/******************************************************************************/
/* this routine figures out where to add a knot  using the Rao criterion */

static int add(mod1,dat,nint,zerror,modmin,mind)
struct model *mod1,*modmin;
struct datas *dat;
int nint,*zerror,mind;

/* mod1   - present model
   dat    - data                     
   nint   - number of integration points
   mind   - minimum distance between knots
   zerror  - zerror status 
   modmin - best model up to now
   mind   - minimum distance between knots */

{
   int i,j,ipowdat[3],ipowvec[3],besti=-1,ll=0,uu=0,nowloc2,bestloc=-1,nx;
   int loloc=0,uploc=0,nowloc1=0;
   double bestrao=-1.,nowrao1,nowrao2,**powdat,**powvec,*sorted;

/* bestrao - bestrao statistic up to now
   nowrao1 - new rao statistic 
   besti  - in between which two knots is the new one                         
   nowrao2 - another new rao statistic 
   bestloc - location of best rao statistic
   newloc1 - location of new rao statistic
   newloc2 - another location of new rao statistic
   i       - counter
   loloc   - smallest possible location
   uploc   - largest possible location
   ll      - potential loloc
   uu      - potential uploc
   find..  - find various locations
   j       - stopper 
   powvec  - piecewise polynomial products, used by hrao
   powdat  - piecewise polynomial products, used by hrao
   rao      - computes a rao-statistic */

/* powvec and powdat are the piecewise polynomial products */
   sorted=wksorted;
   nx=0;
   for(i=0;i<(*dat).nd;i++){
      if((*dat).delta[i]==1){
         sorted[nx]=(*dat).data[i];
         nx++;
      }
   }
   powvec=wkpowvec;
   powdat=wkpowdat;
   if((*mod1).nk!=2){
      for(j=0;j<3;j++){
         ipowvec[j]=-1;
         for(i=0;i<nint;i++){
            if((*mod1).basvec[i]>(*mod1).knots[(*mod1).nk-3+j]){
               powvec[i][j]=pow((*mod1).basvec[i]-(*mod1).knots[(*mod1).nk-3+j],
                                                                    (double)3.);
            }
            else{
               powvec[i][j]=0.;
               ipowvec[j]=i;
            }
         }
      }
      for(j=0;j<3;j++){
         ipowdat[j]=-1;
         for(i=0;i<(*dat).nd;i++){
            if((*dat).data[i]>(*mod1).knots[(*mod1).nk-3+j]){
               powdat[i][j]=pow((*dat).data[i]-(*mod1).knots[(*mod1).nk-3+j],
                                                                    (double)3.);
            }
            else{
               powdat[i][j]=0.;
               ipowdat[j]=i;
            }
         }
      }
   }

/* find the interval */
   for(i=0;i<=(*mod1).nk;i++){

/* before first knot */
      if(i==0 && (*mod1).nk>0)
         nowloc1=hindl(&ll,&uu,mind,sorted,nx,(*mod1).knots[0]);

/* after last knot */
      if(i==(*mod1).nk && (*mod1).nk>0) nowloc1=hindr(&ll,&uu,mind,
                                         sorted,nx,(*mod1).knots[(*mod1).nk-1]);

/* first knot */
      if(i==0 && (*mod1).nk==0) nowloc1=hindx(&ll,&uu,nx);

/* in between knots */
      if(i>0 && i<(*mod1).nk) nowloc1=hindm(&ll,&uu,mind,
                                 sorted,nx,(*mod1).knots[i-1],(*mod1).knots[i]);

/* possible location */
      if(nowloc1>=0){
         nowrao1=hrao(mod1,dat,sorted[nowloc1],nint,powvec,powdat,ipowvec,
                                                                       ipowdat);
         if(nowrao1>bestrao){
            loloc=ll;
            uploc=uu;
            bestloc=nowloc1;
            bestrao=nowrao1;
            besti=i;
         }
      }
   }
   if(bestloc<0)return -1;

/* as long as the locations are different, do interval halving */
   do{
      if(sorted[uploc]>sorted[loloc]){
         nowloc2=hindyr(uploc,bestloc,sorted);

/* two search points, the upper one */
         if(nowloc2>=0){
            nowrao2=hrao(mod1,dat,sorted[nowloc2],nint,powvec,powdat,ipowvec,
                                                                       ipowdat);
         }
         else nowrao2=bestrao;

/* two search points, the lower one */
         nowloc1=hindyl(bestloc,loloc,sorted);
         if(nowloc1>=0){
            nowrao1=hrao(mod1,dat,sorted[nowloc1],nint,powvec,powdat,ipowvec,
                                                                       ipowdat);
         }
         else nowrao1=bestrao;

/* the middle one is the best, we call it quits */
         if(bestrao>=nowrao2 && bestrao>=nowrao1){
            loloc=uploc;
         }
         else{

/* the lower search point is the best */
            if(nowrao1>bestrao){
               uploc=bestloc;
               bestloc=nowloc1;
               bestrao=nowrao1;
            }

/* the upper search point is the best */
            else{
               loloc=bestloc;
               bestloc=nowloc2;
               bestrao=nowrao2;
            }
         }
      }
   }while(sorted[uploc]>sorted[loloc]);

/* failure */
   if(bestloc<0)return bestloc;

/* done record the new knot in its correct position */
   if(besti==(*mod1).nk) 
          (*mod1).knots[(*mod1).nk]=sorted[bestloc];
   else{
      for(i=(*mod1).nk;i>besti;i=i-1){
         (*mod1).knots[i]=(*mod1).knots[i-1];
         (*modmin).iknots[i]=(*modmin).iknots[i-1];
      }
      (*mod1).knots[besti]=sorted[bestloc];
      (*modmin).iknots[besti]=0;
   }
   ((*mod1).nk)++;
   thetaform(mod1,besti);
   if(zerror[6]==37){
      (void)Rprintf("knot added at %.2f ",sorted[bestloc]);
      (void)Rprintf("(rao = %.2f) || ",bestrao);
   }
   return bestloc;
}

/******************************************************************************/
/* these routines compute the rao-score statistic in cand */

static double hrao(mod1,dat,cand,nint,powvec,powdat,ipowvec,ipowdat)
double **powvec,**powdat,cand;
struct model *mod1;
struct datas *dat;
int nint,ipowvec[3],ipowdat[3];

{
   double **info2,*score2,*score3,*newbas,*newdata,lm0,lam,r;
   int i,j,nk=(*mod1).nk;

/* info2    - larger copy of hessian
   score2   - larger copy of score
   score3   - another larger copy of score
   newbas   - new basis function in integration poins
   newdata  - new basis function in data points
   lam      - lambda or exp(lambda)
   lm0      - multiplier times lam
   r        - value of rao
   i,j      - counter 
   nk       - (*mod1).nk 
   newnew   - computes newdata and newbas */

/* allocate memmory */
   info2=wkinfo2;
   score2=wkscore2;
   score3=wkscore3;
   newbas=wknewbas;
   newdata=wknewdata;

/* copy score and info in score2 and info2 */
   score2[nk+1]=0.;
   info2[nk+1][nk+1]=0.;
   for(i=0;i<=nk;i++){
      score2[i]=(*mod1).score[i];
      info2[i][nk+1]=0.;
      info2[nk+1][i]=0.;
      for(j=0;j<=nk;j++){
         info2[i][j]=(*mod1).hessian[i][j];
      }
   }

/* compute newdata and newbas */
   newnew((*mod1).knots,nk,cand,newbas,newdata,nint,dat,
         (*mod1).basvec,powdat,powvec,ipowdat,ipowvec);

/* compute the extra row of info and extra element of score - compare mint */
   for(i=0;i<nint;i++){
      lam=exp(lambda(nk,(*mod1).basmat,(*mod1).theta,i));
      lm0=lam*(*mod1).mult[i]*newbas[i];
      score2[nk+1]+=lm0;
      info2[nk+1][nk+1]+=lm0*newbas[i];
      info2[0][nk+1]+=lm0*(*mod1).basmat[i][0];
      info2[nk-1][nk+1]+=lm0*(*mod1).basmat[i][nk-1];
      info2[nk][nk+1]+=lm0*(*mod1).basmat[i][nk];
      for(j=(int)((*mod1).basmat[i][nk+1]);
         j<=(int)((*mod1).basmat[i][nk+2]) && j>0;j++){
         info2[j][nk+1]+=lm0*(*mod1).basmat[i][j];
      }
   }

/* add the delta part to the score function */
   for(i=0;i<(*dat).nd;i++){
      if((*dat).delta[i]==1) score2[nk+1]+=newdata[i];
   }

/* left tail peculiarities */
   if((*mod1).tails[0]>0.5 || (*mod1).theta[0]<-0.999){
      score2[0]=0.;
      info2[0][0]=-1.;
      for(i=1;i<=nk+1;i++){
         info2[0][i]=0.;
         info2[i][0]=0.;
      }
   }

/* more left tail peculiarities */
   if((*mod1).tails[4]>0.5){
      score2[1]=0.;
      for(i=0;i<=nk+1;i++){
         info2[1][i]=0.;
         info2[i][1]=0.;
      }
      info2[1][1]=-1.;
   }

/* right tail peculiarities */
   if((*mod1).tails[2]>0.5 || (*mod1).theta[nk]<-0.999){
      score2[nk]=0.;
      for(i=0;i<=nk+1;i++){
         info2[nk][i]=0.;
         info2[i][nk]=0.;
      }
      info2[nk][nk]=-1.;
   }

/* symmaterize  and copy in score3 */
   for(j=0;j<=nk+1;j++){
      info2[nk+1][j]=info2[j][nk+1];
      score3[j]=score2[j];
   }

/* compute rao in 2-steps, solving a system and an inner product */
   i=0;
   hlusolve(info2,nk+2,score2,&i);
   r=0.;
   for(i=0;i<nk+2;i++) r+=score2[i]*score3[i];

   return -r;
}
/******************************************************************************/
/* computes the new basisfunction */

static void newnew(knots,nk,cand,newbas,newdata,nint,dat,basvec,
         powdat,powvec,ipowdat,ipowvec)
double *knots,cand,*newbas,*newdata,*basvec,**powdat,**powvec;
int nk,nint,ipowdat[3],ipowvec[3];
struct datas *dat;

/* see all above */

{
   double **mat,vec[3];
   int i,j;

/* i,j - counter
   mat - coefficients of new basisfunction
   vec - coefficients of new basisfunction */

   mat=wkmat33;

/* the basis function is slightly different if there were only two
   knots before */
   if(nk>2){

/* compute the coefficients */
      vec[0]=-1.;
      vec[1]=-cand;
      vec[2]=-cand*cand;
      mat[0][0]=1.;
      mat[0][1]=1.;
      mat[0][2]=1.;
      mat[1][0]=knots[nk-3];
      mat[1][1]=knots[nk-2];
      mat[1][2]=knots[nk-1];
      mat[2][0]=knots[nk-3]*knots[nk-3];
      mat[2][1]=knots[nk-2]*knots[nk-2];
      mat[2][2]=knots[nk-1]*knots[nk-1];
      i=0;
      hlusolve2(mat,(int)3,vec,&i);

/* compute in integration points */
      for(i=0;i<nint;i++){
         newbas[i]=0.;
         if(basvec[i]>cand) newbas[i]=pow((basvec[i]-cand),(double)3.);
      }
      for(j=0;j<3;j++){
         if(ipowvec[j]<nint-1){
            for(i=ipowvec[j]+1;i<nint;i++) newbas[i]+=vec[j]*powvec[i][j];
         }
      }

/* compute in basis points */
      for(i=0;i<(*dat).nd;i++){
         newdata[i]=0.;
         if((*dat).delta[i]==1){
            if((*dat).data[i]>cand){
               newdata[i]=pow(((*dat).data[i]-cand),(double)3.);
            }
         }
      }
      for(j=0;j<3;j++){
         if(ipowdat[j]<(*dat).nd-1){
            for(i=ipowdat[j]+1;i<(*dat).nd;i++){
               if((*dat).delta[i]==1) newdata[i]+=vec[j]*powdat[i][j];
            }
         }
      }
   }

/* if there were only two knots */
   if(nk==2){

/* compute coefficients */
      vec[0]=(knots[1]-cand)/(knots[0]-knots[1]);
      vec[1]=(cand-knots[0])/(knots[0]-knots[1]);

/* compute in data points */
      for(i=0;i<(*dat).nd;i++){
         newdata[i]=0.;
         if((*dat).delta[i]==1){
            if((*dat).data[i]<cand){
               newdata[i]=pow(((*dat).data[i]-cand),(double)3.);
            }
            if((*dat).data[i]<knots[1]){
               newdata[i]+=vec[1]*pow(((*dat).data[i]-knots[1]),(double)3.);
               if((*dat).data[i]<knots[0]){
                  newdata[i]+=vec[0]*pow(((*dat).data[i]-knots[0]),(double)3.);
               }
            }
         }
      }

/* compute in integration points */
      for(i=0;i<nint;i++){
         newbas[i]=0.;
         if(basvec[i]<cand) newbas[i]=pow((basvec[i]-cand),(double)3.);
         if(basvec[i]<knots[1]){
            newbas[i]+=vec[1]*pow((basvec[i]-knots[1]),(double)3.);
            if(basvec[i]<knots[0]){
               newbas[i]+=vec[0]*pow((basvec[i]-knots[0]),(double)3.);
            }
         }
      }
   }
}
/******************************************************************************/
static void thetaform(mod1,besti)

struct model *mod1;
int besti;

/* mod1  - present model
   besti - index of new basisfunction */

{
   double *phi;
   int i,j;

/* phi - theta for (part of the) power basis */

   phi=wkphi7;
   for(i=0;i<HLENGTH;i++) phi[i]=0;
   for(j=0;j<(*mod1).nk-2;j++){
      for(i=0;i<(*mod1).nk+1;i++)phi[i]+=(*mod1).theta[j+1]*(*mod1).coef2[j][i];
   }
   (*mod1).theta[(*mod1).nk+3]=(*mod1).theta[(*mod1).nk-1];
   (*mod1).theta[(*mod1).nk+2]=(*mod1).theta[0];
   (*mod1).theta[1]=phi[1];
   (*mod1).theta[0]=phi[0];
   (*mod1).theta[besti+2]=0.;
   for(i=0;i<=(*mod1).nk-1;i++){
      if(i<besti) (*mod1).theta[i+2]=phi[i+2];
      if(i>besti) (*mod1).theta[i+2]=phi[i+1];
   }
}
/******************************************************************************/
/* finds a new location in an interval (l,b) - that is the lower end might not
   have been tested yet */
static int hindyl(u,l,x)
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
static int hindyr(u,l,x)
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

static int hindl(ll,uu,mind,x,nx,knt)
double *x,knt;
int nx,*ll,*uu,mind;
{

/* i  - utility
   hlocation - finds uu */

   int i;

   (*uu)=hlocation(0,x,nx,knt);
   if((*uu)<mind-1)return -1;
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

static int hindr(ll,uu,mind,x,nx,knt)
double *x,knt;
int nx,*ll,*uu,mind;
{

/* i  - utility
   hlocation - finds ll */

   int i;

   (*ll)=hlocation(1,x,nx,knt);
   if(nx-1-(*ll)<=mind)return -1;
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

static int hindx(ll,uu,nx)
int nx,*ll,*uu;
{
   *ll=0;
   *uu=nx-1;
   return nx/2;
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

static int hindm(ll,uu,mind,x,nx,k0,k1)
double *x,k0,k1;
int nx,*ll,*uu,mind;
{
/* hlocation - finds ll */


   (*ll)=hlocation(1,x,nx,k0);
   (*uu)=hlocation(0,x,nx,k1);
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

static int hlocation(what,x,nx,k)
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
static void dubmodel(m2,m1)
/* copies model m1 into model m2 */
struct model *m1,*m2;
{
   int i1,i2;
   (*m2).tailse[0]=(*m1).tailse[0];
   (*m2).tailse[1]=(*m1).tailse[1];
   (*m2).ll=(*m1).ll;
   (*m2).nk=(*m1).nk;
   (*m2).nk1=(*m1).nk1;
   (*m2).aic=(*m1).aic;
   for(i1=0;i1<HLENGTH;i1++){
      (*m2).iknots[i1]=(*m1).iknots[i1];
      (*m2).logl[i1]=(*m1).logl[i1];
      (*m2).knots[i1]=(*m1).knots[i1];
      (*m2).yknots[i1]=(*m1).yknots[i1];
      (*m2).theta[i1]=(*m1).theta[i1];
      for(i2=0;i2<HLENGTH;i2++){
         (*m2).icoef[i1][i2]=(*m1).icoef[i1][i2];
         (*m2).coef2[i1][i2]=(*m1).coef2[i1][i2];
      }
   }
   for(i1=0;i1<5;i1++) (*m2).tails[i1]=(*m1).tails[i1];
}
/******************************************************************************/
static void hlusolve(a,n,b,k)
int n,*k;
double **a,*b;
{
   double aa[HLENGTH][HLENGTH],bb[HLENGTH];
   int kpvt[HLENGTH],info;
   int i,j;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
      bb[i]=b[i];
   }
   i=HLENGTH;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
   (*k)=0;
   if(info!=0)(*k)= -1;
   F77_CALL(xdsisl)(aa,&i,&n,kpvt,bb);
   for(i=0;i<n;i++)b[i]=bb[i];
}
/******************************************************************************/
static void hlusolve2(a,n,b,k)
int n,*k;
double **a,*b;
{
   double aa[HLENGTH][HLENGTH],bb[HLENGTH];
   int kpvt[HLENGTH],info;
   int i,j;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
      bb[i]=b[i];
   }
   i=HLENGTH;
   F77_CALL(xdgefa)(aa,&i,&n,kpvt,&info);
   (*k)=0;
   if(info!=0)(*k)= -1;
   j=0;
   F77_CALL(xdgesl)(aa,&i,&n,kpvt,bb,&j);
   for(i=0;i<n;i++)b[i]=bb[i];
}
/******************************************************************************/
static void hluinverse(a,n)
int n;
double **a;
{
   double aa[HLENGTH][HLENGTH],bb[HLENGTH],det[2]; int inert[3];
   int kpvt[HLENGTH],info;
   int i,j;
   for(i=0;i<n;i++){
      for(j=0;j<n;j++)aa[i][j]=a[j][i];
   }
   i=HLENGTH;
   j=1;
   F77_CALL(xdsifa)(aa,&i,&n,kpvt,&info);
   F77_CALL(xdsidi)(aa,&i,&n,kpvt,det,inert,bb,&j);
   for(i=0;i<n;i++){
      for(j=0;j<i;j++)a[i][j]=aa[i][j];
      for(j=i;j<n;j++)a[i][j]=aa[j][i];
   }
}
/******************************************************************************/
/* computes heft probabilities or quantiles */

void heftpq(knots,cc,thetak,thetal,thetap,what,pp,qq,nk,np)
int *what,*np,*nk;
double *knots,*cc,*thetak,*thetal,*thetap,*pp,*qq;

/* knots   - knots
   cc      - median of data
   thetak  - theta related to (x-t)+^3
   thetal  - theta related to log(x/(x+cc)) and log(x+cc)
   thetap  - theta related to 1 and x
   what    - if 0: get quantiles from probabilities
             if 1: get probabilities from quantiles
   pp      - probabilities
   qq      - quantiles
   nk      - number of knots
   np      - number of points of interest */

{
   double x=0.,z=0.,zk,xl=0.,zl=0.,dpp=30.;
   int i=0,j,l=0;

/* x   - last point until where the complete integration has been carried out
   z   - value of  -log(1-p) at x
   xl  - next integration point on one tenth of knot distant
   zl  - value of  -log(1-p) at xl
   zk  - value of  -log(1-p) at next knot
   i   - counter (for intervals)
   j   - counter (for points of interest)
   l   - keeps track how far we are inbetween knots */

/* get probabilities from quantiles */
   if(*what==1){
      for(j=0;j<*np;j++){

/* extreme cases */
         if(qq[j]<0.) pp[j]=0.;

/* first integrate until the closest knot before qq[j], start at last point */
         else{
            if(qq[j]>knots[i] && i< *nk){
               do{
                  z+=ilambda(knots,*cc,thetak,thetal,thetap,*nk,x,knots[i],i);
                  x=knots[i];
                  i++;
               } while(qq[j]>knots[i] && i< *nk);
            }

/* then integrate to qq[j] */
            z+=ilambda(knots,*cc,thetak,thetal,thetap,*nk,x,qq[j],i);
            pp[j]=1-exp(-z);
            x=qq[j];
         }
      }
   }

/* get quantiles from probabilities */
   else{

/* first compute -log(1-p) in first knot */
      zk=ilambda(knots,*cc,thetak,thetal,thetap,*nk,(double)0.,knots[0],0);

      for(j=0;j<*np;j++){
         if(pp[j]>0 && pp[j]<1){
            pp[j]= -hmylog(1-pp[j]);

/* first integrate until the closest knot before pp[j] */
            if(pp[j]>zk && i < *nk){
               do{
                  z=zk;
                  x=knots[i];
                  i++;
                  zk+=ilambda(knots,*cc,thetak,thetal,thetap,*nk,x,knots[i],i);
                  xl=x;
                  zl=0;
                  l=0;
               } while(pp[j]>zk && i<*nk);
            }

/* then takes steps of one tenth of the interval */
            if(pp[j]>z+zl){
               do{
                  l++;
                  if(i<*nk && i>0){
                     x=xl;
                     z+=zl;
                     xl=(double)l/dpp*knots[i]+(dpp-l)/dpp*knots[i-1];
                  }
                  if(i==0){
                     x=xl;
                     z+=zl;
                     xl=(double)l/dpp*knots[i];
                  }

/* outside the most extreme knot, we double the distance */
                  if(i==*nk){
                     z+=zl;
                     x=xl;
                     xl=2.*(x-knots[*nk-2])+knots[*nk-2];
                  }
                  zl=ilambda(knots,*cc,thetak,thetal,thetap,*nk,x,xl,i);
               } while(pp[j]>z+zl);
            } 

/* linear interpolate further */
            qq[j]=x+(pp[j]-z)/zl*(xl-x);
         } 
      }
   }
}

/******************************************************************************/

static double xlambda(knots,cc,thetak,thetal,thetap,nk,x)
double *knots,cc,*thetak,*thetal,*thetap,x;
int nk;

/* computes exp(lambda(x)), all quantities, see above */

{
   double y;
   int i;
   if(x>0){
      y=thetap[0]+x*thetap[1]+hmylog(x+cc)*thetal[1]+hmylog(x/(x+cc))*thetal[0];
      for(i=0;i<nk && x>knots[i];i++) 
         y+=(x-knots[i])*(x-knots[i])*(x-knots[i])*thetak[i];
      return exp(y);
   }

/* if x is 0 forget about log(x/(x+c)) */
   else{
      y=thetap[0]+x*thetap[1]+hmylog(x+cc)*thetal[1];
      for(i=0;i<nk && x>knots[i];i++)
         y+=(x-knots[i])*(x-knots[i])*(x-knots[i])*thetak[i];
      return exp(y);
   }
}

/******************************************************************************/

static double ilambda(knots,cc,thetak,thetal,thetap,nk,z1,z2,i)
double *knots,cc,*thetak,*thetal,*thetap,z1,z2;
int nk,i;

/* integrates exp(lambda(x)) from z1 to z2, which is between knot[i-1] and
   knot[i] (knot[-1]=0, knot[nk]=infty) */
{
   double r1,r2,y[60],w[60],f;
   int k;
   r1 = (z2-z1)/2.;
   r2 = (z2+z1)/2.;

/* Gaussian quadrature - see Abromowitz and Stegun */

   y[0] = 0.125233408511469 * r1; w[0] = 0.249147045813403 * r1;
   y[1] = 0.367831498998180 * r1; w[1] = 0.233492536538355 * r1;
   y[2] = 0.587317954286617 * r1; w[2] = 0.203167426723066 * r1;
   y[3] = 0.769902674194305 * r1; w[3] = 0.160078328543346 * r1;
   y[4] = 0.904117256370475 * r1; w[4] = 0.106939325995318 * r1;
   y[5] = 0.981560634246719 * r1; w[5] = 0.047175336386512 * r1;
   k=6;  /*
   w[ 0]=  0.00178328072169643 * r1; y[0 ]=  0.99930504173577217 * r1;
   w[ 1]=  0.00414703326056247 * r1; y[1 ]=  0.99634011677195533 * r1;
   w[ 2]=  0.00650445796897836 * r1; y[2 ]=  0.99101337147674429 * r1;
   w[ 3]=  0.00884675982636395 * r1; y[3 ]=  0.98333625388462598 * r1;
   w[ 4]=  0.01116813946013113 * r1; y[4 ]=  0.97332682778991098 * r1;
   w[ 5]=  0.01346304789671864 * r1; y[5 ]=  0.96100879965205377 * r1;
   w[ 6]=  0.01572603047602472 * r1; y[6 ]=  0.94641137485840277 * r1;
   w[ 7]=  0.01795171577569734 * r1; y[7 ]=  0.92956917213193957 * r1;
   w[ 8]=  0.02013482315353021 * r1; y[8 ]=  0.91052213707850282 * r1;
   w[ 9]=  0.02227017380838325 * r1; y[9 ]=  0.88931544599511414 * r1;
   w[10]=  0.02435270256871087 * r1; y[10]=  0.86599939815409277 * r1;
   w[11]=  0.02637746971505466 * r1; y[11]=  0.84062929625258032 * r1;
   w[12]=  0.02833967261425948 * r1; y[12]=  0.81326531512279754 * r1;
   w[13]=  0.03023465707240248 * r1; y[13]=  0.78397235894334139 * r1;
   w[14]=  0.03205792835485155 * r1; y[14]=  0.75281990726053194 * r1;
   w[15]=  0.03380516183714161 * r1; y[15]=  0.71988185017161088 * r1;
   w[16]=  0.03547221325688239 * r1; y[16]=  0.68523631305423327 * r1;
   w[17]=  0.03705512854024005 * r1; y[17]=  0.64896547125465731 * r1;
   w[18]=  0.03855015317861563 * r1; y[18]=  0.61115535517239328 * r1;
   w[19]=  0.03995374113272034 * r1; y[19]=  0.57189564620263400 * r1;
   w[20]=  0.04126256324262353 * r1; y[20]=  0.53127946401989457 * r1;
   w[21]=  0.04247351512365359 * r1; y[21]=  0.48940314570705296 * r1;
   w[22]=  0.04358372452932345 * r1; y[22]=  0.44636601725346409 * r1;
   w[23]=  0.04459055816375657 * r1; y[23]=  0.40227015796399163 * r1;
   w[24]=  0.04549162792741814 * r1; y[24]=  0.35722015833766813 * r1;
   w[25]=  0.04628479658131442 * r1; y[25]=  0.31132287199021097 * r1;
   w[26]=  0.04696818281621002 * r1; y[26]=  0.26468716220876742 * r1;
   w[27]=  0.04754016571483031 * r1; y[27]=  0.21742364374000708 * r1;
   w[28]=  0.04799938859645831 * r1; y[28]=  0.16964442042399283 * r1;
   w[29]=  0.04834476223480295 * r1; y[29]=  0.12146281929612056 * r1;
   w[30]=  0.04857546744150343 * r1; y[30]=  0.07299312178779904 * r1;
   w[31]=  0.04869095700913972 * r1; y[31]=  0.02435029266342443 * r1;
   k=32; */
   f=0.;
   for(i=0;i<k;i++){
      f+= w[i]*(xlambda(knots,cc,thetak,thetal,thetap,nk,r2-y[i])
               +xlambda(knots,cc,thetak,thetal,thetap,nk,r2+y[i]));
   }
   return f;
}
static void allocer(nd,i00)
int nd,i00;
{
   wkddd=dhvector(nd);
   wkvec2=wkddd;
   wkmat1=dhmatrix(HLENGTH,HLENGTH);
   wkvec1=dhvector(HLENGTH);
   wkphi=dhvector(HLENGTH+2);
   wkmat=dhmatrix(HLENGTH+2,HLENGTH+2);
   wkphi2=dhvector(HLENGTH);
   wkmasterpt = dhvector(HLENGTH+100);
   wkxx = dhvector(HLENGTH+100);
   wkcand = dhvector(HLENGTH);
   wkse3=dhvector(HLENGTH);
   wkphi3=dhvector(HLENGTH);
   wkphi7=dhvector(HLENGTH);
   wkphi4=dhvector(HLENGTH);
   wkhh=dhmatrix(HLENGTH,HLENGTH);
   wkpowdat=dhmatrix(nd,2);
   wksorted=wkddd;
   wkpowvec=dhmatrix(i00,2);
   wkinfo2=dhmatrix(HLENGTH,HLENGTH);
   wkscore2=dhvector(HLENGTH);
   wkscore3=dhvector(HLENGTH);
   wknewbas=dhvector(i00);
   wknewdata=dhvector(nd);
   wkmat33=dhmatrix(3,3);
}
static double hmylog(x)
double x;
{
if(x < 10.e-250)return (double)(-575.64627);
else return log(x);
}

