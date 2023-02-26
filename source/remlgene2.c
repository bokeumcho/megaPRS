/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//REML with one gene - have not added Bayesian version

///////////////////////////

double gene_reml(int ns, double *Y2, int num_fixed, double detZTZ, double *X, int Xtotal, double *sqdevs, double Xsum, double *Xnss, double *Xrhos, double *stats, double tol, int maxiter)
{
int i, j, j2, k, count, best, one=1, lwork, info, flag;
double nfree, scale, sum, value, alpha, beta, wkopt, *work;
double lambda, lambda2, lambdadiff, her, likenull, like, like2, likeold, diff;

double S1, S2, S3, T1, T2, T3, gam, deriv, dderiv, dderiv2;
double YTCY, *XTCY, *D;

double begin[19]={0,1e-6,5e-6,.00001,.00005,.0001,.0002,.0005,.001,.002,.005,.01,.02,.05,.1,.2,.5,.99,.995};


//set nfree and scale (how much to scale up)
if(Xnss==NULL)	//using real data
{
nfree=ns-num_fixed;
scale=1;
}
else	//using summary statistics
{
sum=0;for(j=0;j<Xtotal;j++){sum+=Xnss[j];}
nfree=sum/Xtotal;
scale=nfree/ns;
}

//allocate variables - will always have Xtotal>0

XTCY=malloc(sizeof(double)*Xtotal);
D=malloc(sizeof(double)*Xtotal);

//get YTCY
if(Xnss==NULL)	//using real data
{
YTCY=0;for(i=0;i<ns;i++){YTCY+=pow(Y2[i],2);}
}
else	//using summary statistics - assume variance of phenotypes is one
{
YTCY=nfree;
}

//get likelihood
likenull=-.5*(nfree+nfree*log(2*M_PI*YTCY/nfree)+detZTZ);

//get XTCX (when using summary statistics, this equals covar(X)*nfree)
alpha=scale;beta=0.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &ns, &alpha, X, &ns, X, &ns, &beta, XTCX, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &num_fixed, &alpha, ZTX, &num_fixed, ZTZZTX, &num_fixed, &beta, XTCX, &Xtotal);

//get XTCY - only use summary statistics for real analysis (not when permuting)
if(Xnss==NULL)	//use real data
{
alpha=scale;beta=0.0;
dgemv_("T", &ns, &Xtotal, &alpha, X, &ns, Y, &one, &beta, XTCY, &one);
alpha=-1.0;beta=1.0;
dgemv_("T", &num_fixed, &Xtotal, &alpha, ZTX, &num_fixed, ZTZZTY, &one, &beta, XTCY, &one);
}
else	//from summaries - the diagonal of XTCX is the sum of squares
{
for(j=0;j<Xtotal;j++){XTCY[j]=Xrhos[j]*pow(XTCX[j+j*Xtotal]*YTCY,.5);}
}

if(shrink!=1)	//alter XTCX
{
if(shrink<1)	//deflate off diagonal terms
{
for(j=0;j<Xtotal;j++)
{
for(j2=j+1;j2<Xtotal;j2++){XTCX[j+j2*Xtotal]*=shrink;XTCX[j2+j*Xtotal]*=shrink;}
}
}
if(shrink>1)	//inflate diagonal terms
{
for(j=0;j<Xtotal;j++){XTCX[j+j*Xtotal]*=shrink;}
}
}

//decomp XTCX
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){U[j+j2*Xtotal]=XTCX[j+j2*Xtotal];}
}
lwork=-1;
dsyev_("V", "U", &Xtotal, U, &Xtotal, E, &wkopt, &lwork, &info);
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &Xtotal, U, &Xtotal, E, work, &lwork, &info);
free(work);

/*
if(strip>0)	//set end proportion of E to zero
{
sum=0;for(j=0;j<Xtotal;j++){sum+=E[j];}
sum2=0;
for(j=Xtotal-1;j>=0;j--)
{
sum2+=E[j];
if(sum2/sum>1-strip){E[j]=0;}
}
}
*/

//get D = UTXTCY
alpha=1.0;beta=0.0;
dgemv_("T", &Xtotal, &Xtotal, &alpha, U, &Xtotal, XTCY, &one, &beta, D, &one);

//////////////////////////

//likelihood -.5 [nfree (1 + log(2pi gamma/nfree)) + log |E+lambda| - Xtotal log(lambda) + detZTZ]
//where gamma is YTCY - DT (E + lambda)^-1 D

flag=0;	//indicates whether there is a problem with the gene (negative gamma)

if(Xnss!=NULL)	//compute gamma corresponding to least squares estimate (ie lambda=0)
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

like2=-.5*(nfree+nfree*log(2*M_PI*gam/nfree)+T1-Xtotal*log(lambda2)+detZTZ);

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

like2=-.5*(nfree+nfree*log(2*M_PI*gam/nfree)+T1-Xtotal*log(lambda2)+detZTZ);

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

free(ZTY);free(ZTX);free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTZZTY);free(ZTZZTX);
free(XTCY);free(XTCX);free(U);free(E);free(D);

return(2*(like-likenull));
}	//end of reml_lite

///////////////////////////

