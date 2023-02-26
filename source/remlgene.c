/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//REML with one gene - have not added Bayesian version

///////////////////////////

double gene_reml(int ns, int num_fixed, double *Yadj, double detZTZ, double *X, int Xtotal, double *U, double *E2, double Xsum, double *Xnss, double *Xrhos, double *Xsqs, double *stats, double tol, int maxiter)
{
int i, j, k, count, best, one=1, flag;
double nfree, scale, sum, value, alpha, beta;
double lambda, lambda2, lambdadiff, her, likenull, like, like2, likeold, diff;

double S1, S2, S3, T1, T2, T3, gam, deriv, dderiv, dderiv2;
double YTCY, *XTCY, *D, *E;

double begin[19]={0,1e-6,5e-6,.00001,.00005,.0001,.0002,.0005,.001,.002,.005,.01,.02,.05,.1,.2,.5,.99,.995};


//set nfree and scale (how much to scale up YTCY, XTCX, detZTZ and E (but not U) when using summary statistics)
if(Xnss==NULL)	//using real data
{
nfree=ns-num_fixed;
scale=1;
}
else	//using summaries
{
sum=0;for(j=0;j<Xtotal;j++){sum+=Xnss[j];}
nfree=sum/Xtotal;
scale=nfree/ns;
}

//allocate variables - will always have Xtotal>0

XTCY=malloc(sizeof(double)*Xtotal);
D=malloc(sizeof(double)*Xtotal);
E=malloc(sizeof(double)*Xtotal);

//get YTCY and likenull
YTCY=0;for(i=0;i<ns;i++){YTCY+=scale*pow(Yadj[i],2);}
likenull=-.5*(nfree+nfree*log(2*M_PI*YTCY/nfree)+detZTZ+log(scale));

//get XTCY
if(Xnss==NULL)	//using real data
{
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &Xtotal, &alpha, X, &ns, Yadj, &one, &beta, XTCY, &one);
}
else	//using summaries
{
for(j=0;j<Xtotal;j++){XTCY[j]=Xrhos[j]*pow(scale*Xsqs[j]*YTCY,.5);}
}

//get D = UTXTCY
alpha=1.0;beta=0.0;
dgemv_("T", &Xtotal, &Xtotal, &alpha, U, &Xtotal, XTCY, &one, &beta, D, &one);

//load up E
for(j=0;j<Xtotal;j++){E[j]=E2[j]*scale;}

//////////////////////////

//likelihood -.5 [nfree (1 + log(2pi gamma/nfree)) + log |E+lambda| - Xtotal log(lambda) + detZTZ]
//where gamma is YTCY - DT (E + lambda)^-1 D

flag=0;	//indicates whether there is a problem with the gene (negative gamma)

if(Xnss!=NULL)	//check least-squares gamma (ie lambda=0) is sensible
{
S1=0;
for(j=0;j<Xtotal;j++)
{
if(E[j]!=0){S1+=pow(D[j],2)/E[j];};
}
gam=YTCY-S1;

if(S1>=YTCY||S1<=0){flag=1;}
}

//starting heritability is zero
best=0;
like=likenull;

if(flag==0)	//see if we can find a better starting value
{
for(k=1;k<19;k++)
{
her=begin[k];
lambda2=Xsum*(1-her)/her;

S1=0;T1=0;
for(j=0;j<Xtotal;j++)
{
S1+=pow(D[j],2)/(E[j]+lambda2);
if(E[j]+lambda2>0){T1+=log(E[j]+lambda2);}
}
gam=YTCY-S1;

if(gam<=0)	//should now never happen
{printf("Error, negative gamma for starting heritabilities, please tell Doug %f %f | %f and her %f perm %d %f\n\n", S1, T1, gam, her, (stats==NULL), begin[best]);exit(1);}

like2=-.5*(nfree+nfree*log(2*M_PI*gam/nfree)+T1-Xtotal*log(lambda2)+detZTZ+log(scale));

if(like2>like)	//have found a better value
{best=k;like=like2;}
else	//will only accept the first mode, so can break
{break;}
}
}

////////

if(best==0||best==18)	//starting her is either 0 or (approx) 1, so will not iterate
{
if(best==0)	//best her is 0 - final estimate will be 1e-6
{her=1e-6;}
if(best==18)	//best her is 0.995 - final estimate will be 0.99
{her=0.99;}
like=likenull;

if(stats!=NULL)	//fill stats
{
stats[0]=her;
stats[1]=-9999;
stats[2]=likenull;
stats[3]=like;
stats[4]=2*(like-likenull);
if(stats[4]>0){stats[5]=.5*erfc(pow(stats[4],.5)*M_SQRT1_2);}
else{stats[5]=.75;}
}
}
else	//best her was not one of the extremes, so will iterate
{
her=begin[best]; 
lambda=Xsum*(1-her)/her;

//reset first and last values of begin so that her does not go outside [1e-6,0.99]
begin[0]=1e-6;
begin[18]=0.99;

count=0;
while(1)
{
//get likelihood and derivatives
S1=0;S2=0;S3=0;T1=0;T2=0;T3=0;
for(j=0;j<Xtotal;j++)
{
S1+=pow(D[j],2)/(E[j]+lambda);
S2+=pow(D[j],2)*pow(E[j]+lambda,-2);
S3+=pow(D[j],2)*pow(E[j]+lambda,-3);
if(E[j]+lambda>0){T1+=log(E[j]+lambda);}
T2+=1/(E[j]+lambda);
T3+=pow(E[j]+lambda,-2);
}
gam=YTCY-S1;

deriv=-.5*nfree/gam*S2-.5*T2+.5*Xtotal/lambda;
dderiv=-.5*nfree*pow(gam,-2)*pow(S2,2)+nfree/gam*S3+.5*T3-.5*Xtotal*pow(lambda,-2);

//see if breaking
if(count>0)	//test for convergence
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;
if(count==maxiter){printf("Warning, REML failed to converge within %d iterations - this is only a problem if it happens very often\n", maxiter);break;}

//get proposed move, ensuring heritability remains within neighbouring start points
lambdadiff=-deriv/dderiv;

her=Xsum/(Xsum+lambda+lambdadiff);
if(her<begin[best-1]){lambdadiff=Xsum*(1-begin[best-1])/begin[best-1]-lambda;}
if(her>begin[best+1]){lambdadiff=Xsum*(1-begin[best+1])/begin[best+1]-lambda;}

value=1;
while(value>0.0001)
{
//get likelihood based on moving value*lambdadiff
lambda2=lambda+value*lambdadiff;
her=Xsum/(Xsum+lambda2);
S1=0;T1=0;
for(j=0;j<Xtotal;j++)
{
S1+=pow(D[j],2)/(E[j]+lambda2);
if(E[j]+lambda2>0){T1+=log(E[j]+lambda2);}
}
gam=YTCY-S1;

if(gam<=0)	//should now never happen
{printf("Error, negative gamma, please tell Doug %f %f | %f and her %f perm %d %f\n\n", S1, T1, gam, her, (stats==NULL), begin[best]);exit(1);}

like2=-.5*(nfree+nfree*log(2*M_PI*gam/nfree)+T1-Xtotal*log(lambda2)+detZTZ+log(scale));

//see whether to accept move, or switch to a smaller move
if(like2>like-tol){lambda=lambda2;like=like2;break;}
else{value*=.5;}
}

count++;
}	//end of while loop

if(stats!=NULL)	//save statistics
{
her=Xsum/(Xsum+lambda);
dderiv2=pow(scale*Xsum,2)*pow(her,-4)*dderiv+2*scale*Xsum*pow(her,-3)*deriv;

stats[0]=her;
if(dderiv2<0){stats[1]=pow(-dderiv2,-.5);}
else{stats[1]=-9999;}
stats[2]=likenull;
stats[3]=like;
stats[4]=2*(like-likenull);
if(stats[4]>0){stats[5]=.5*erfc(pow(stats[4],.5)*M_SQRT1_2);}
else{stats[5]=.75;}
}	//end of saving

}	//end of iterating

free(XTCY);free(D);free(E);

return(2*(like-likenull));
}	//end of reml_lite

///////////////////////////

