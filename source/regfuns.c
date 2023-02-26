/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Regression functions

///////////////////////////

double get_factor(double *Y, int ns, double prev, double ascer, char *outfile)
{
int i;
double sum, mean, value, factor;

char filename[500];
FILE *output;


if(ascer==-9999)	//get mean from phenotypes
{
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
mean=sum/ns;
}
else{mean=ascer;}

value=pow(2*M_PI,-.5)*exp(-.5*pow(normal_inv(1-prev),2));
factor=pow(prev*(1-prev),2)/mean/(1-mean)*pow(value,-2);
printf("Prevalence: %.6f, Ascertainment: %.2f, Scaling factor (obs->liab): %.2f\n\n", prev, mean, factor);

if(outfile!=NULL)
{
sprintf(filename,"%s.factor", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Prevalence %.6f\nAscertainment %.6f\nScaling_Factor_Observed->Liability %.6f\nScaling_Factor_Liability->Observed %.6f\n", prev, mean, factor, 1.0/factor);
fclose(output);
}

return(factor);
}	//end of get_factor

///////////////////////////

void reg_covar_lin(double *Y, double *Z, int ns, int num_covars, int num_tops, double *thetas, double *thetasds, double *thetapvas, double *Yadj, int type, double *covher, double *topher)
//regress Y on Z linearly, then maybe fill Yadj with residuals and get fixed effect heritabilities
//type=0 - dont standardize residuals, type=1 - do
{
int i, j, num_fixed, one=1;
double sum, sumsq, value, origss, fixedss, covarss, alpha, beta;
double *ZTZ, *ZTZ2, *ZTZ3, *ZTY, *thetas2;


num_fixed=num_covars+num_tops;

ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTY=malloc(sizeof(double)*num_fixed);

//get sumsq and original ss
sum=0;sumsq=0;
for(i=0;i<ns;i++){sum+=Y[i];sumsq+=pow(Y[i],2);}
origss=sumsq-sum/ns*sum;

//solve with all fixed, and get fixed ss
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, thetas, &one);

fixedss=sumsq;for(j=0;j<num_fixed;j++){fixedss-=ZTY[j]*thetas[j];}

if(thetasds!=NULL)	//get sds and pvalues
{
for(j=0;j<num_fixed;j++)
{
thetasds[j]=pow(ZTZ[j+j*num_fixed]*fixedss/(ns-num_fixed),.5);
thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);
}
}

if(Yadj!=NULL)	//put residuals into Yadj
{
for(i=0;i<ns;i++){Yadj[i]=Y[i];}
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Yadj, &one);

if(type==1)	//standardize
{
value=pow(fixedss/ns,-.5);
for(i=0;i<ns;i++){Yadj[i]=Yadj[i]*value;}
}
}

if(covher!=NULL)	//solve for just covariates, and get covher and topher
{
if(num_tops==0)	//already done
{covarss=fixedss;}
else
{
thetas2=malloc(sizeof(double)*num_covars);
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_covars);
dgemv_("T", &ns, &num_covars, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
(void)eigen_invert(ZTZ, num_covars, ZTZ2, -1, ZTZ3, 1);
dgemv_("N", &num_covars, &num_covars, &alpha, ZTZ, &num_covars, ZTY, &one, &beta, thetas2, &one);

covarss=sumsq;for(j=0;j<num_covars;j++){covarss-=ZTY[j]*thetas2[j];}
free(thetas2);
}

*covher=(origss-covarss)/origss;
*topher=(covarss-fixedss)/covarss;
}

free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTY);
}	//end of reg_covar_lin	

////////

void reg_covar_matrix(double *X, double *Z, int ns, int length, int num_fixed)
{
//replace X with residual from regressing on Z
double alpha, beta;
double *ZTZ, *ZTZ2, *ZTX;


ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTX=malloc(sizeof(double)*num_fixed*length);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
dgemm_("T", "N", &num_fixed, &length, &ns, &alpha, Z, &ns, X, &ns, &beta, ZTX, &num_fixed);

(void)eigen_invert(ZTZ, num_fixed, ZTZ2, length, ZTX, 1);

alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &length, &num_fixed, &alpha, Z, &ns, ZTX, &num_fixed, &beta, X, &ns);

free(ZTZ);free(ZTZ2);free(ZTX);
}

////////

void reg_single_lin(double *Y, double *Z, int ns, int num_fixed, double *stats)
{
int i, j, one=1;
double sumsq, alpha, beta;
double *ZTZ, *ZTZ2, *ZTZ3, *ZTY, *thetas;


ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTY=malloc(sizeof(double)*num_fixed);
thetas=malloc(sizeof(double)*num_fixed);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, thetas, &one);

sumsq=0;
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
for(j=0;j<num_fixed;j++){sumsq-=ZTY[j]*thetas[j];}

stats[0]=thetas[num_fixed-1];
stats[1]=pow(ZTZ[num_fixed*num_fixed-1]*sumsq/(ns-num_fixed),.5);
stats[2]=stats[0]/stats[1];
stats[3]=erfc(fabs(stats[0]/stats[1])*M_SQRT1_2);

free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTY);free(thetas);
}	//end of reg_single_lin	

///////////////////////////

void reg_covar_log(double *Y, double *Z, int ns, int num_covars, int num_tops, double *thetas, double *thetasds, double *thetapvas, double *Yadj, double *covher, double *topher, double prev, double tol, int maxiter)
//regress Y on Z logistically, then maybe fill Yadj with function of liablities and get fixed effect heritabilities
{
int i, j, j2, count, count2, num_fixed, one=1;
double sum, sum2, sumsq, sumsq2, var, var2, value, value2;
double ascer, like, likeold, diff, diff2, fixedvar, covarvar, alpha, beta;
double *Zthetas, *probs, *AI, *AI2, *AI3, *BI, *Ks, *thresh, *regs, *thetas2;


num_fixed=num_covars+num_tops;

Zthetas=malloc(sizeof(double)*ns);
probs=malloc(sizeof(double)*ns);
AI=malloc(sizeof(double)*num_fixed*num_fixed);
AI2=malloc(sizeof(double)*num_fixed);
AI3=malloc(sizeof(double)*num_fixed*num_fixed);
BI=malloc(sizeof(double)*num_fixed);

if(prev!=-9999)
{
Ks=malloc(sizeof(double)*ns);
thresh=malloc(sizeof(double)*ns);
regs=malloc(sizeof(double)*ns);
}

//get ascer
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
ascer=sum/ns;

//solve for all fixed - first covariate starts at ascer, rest are zero
thetas[0]=log(ascer/(1-ascer));
for(j=1;j<num_fixed;j++){thetas[j]=0;}

count=0;
while(1)
{
//get Zthetas then probs
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Zthetas, &one);
for(i=0;i<ns;i++){probs[i]=pow(1+exp(-Zthetas[i]),-1);}

//now likelihood and diff
like=0;for(i=0;i<ns;i++){like+=Y[i]*log(probs[i])+(1-Y[i])*log(1-probs[i]);}
if(count>1){diff2=diff;}
if(count>0){diff=like-likeold;}
likeold=like;

//AI = -2nd deriv = sum(pi(1-pi)ZiZi) and Bi = 1st deriv = sum(Yi-pi)Zi
for(j=0;j<num_fixed;j++)
{
BI[j]=0;for(i=0;i<ns;i++){BI[j]+=(Y[i]-probs[i])*Z[i+j*ns];}
AI[j+j*num_fixed]=0;for(i=0;i<ns;i++){AI[j+j*num_fixed]+=probs[i]*(1-probs[i])*pow(Z[i+j*ns],2);}
for(j2=j+1;j2<num_fixed;j2++)
{
AI[j+j2*num_fixed]=0;for(i=0;i<ns;i++){AI[j+j2*num_fixed]+=probs[i]*(1-probs[i])*Z[i+j*ns]*Z[i+j2*ns];}
AI[j2+j*num_fixed]=AI[j+j2*num_fixed];
}
}

//invert AI
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);

if(count>1)	//will have set diff and diff2 so can check if converged
{
if(fabs(diff)<tol&&fabs(diff2)<tol){break;}
}
if(count==maxiter){printf("Error, logistic regression did not converge after %d iterations; consider reducing the number of fixed effects, or using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %.6f)\n\n", maxiter, tol);exit(1);}

//update thetas using newton raphson
alpha=1.0;beta=1.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, AI, &num_fixed, BI, &one, &beta, thetas, &one);

count++;
}

//get sds and pvalues
for(j=0;j<num_fixed;j++)
{
if(AI[j+j*num_fixed]>0)
{thetasds[j]=pow(AI[j+j*num_fixed],.5);
thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);}
else{thetasds[j]=-9999;thetapvas[j]=-9999;}
}

if(prev!=-9999)	//get thresholds, then fill Yadj and/or get covarvar and topher 
{
//get Ks (population probs of being case given covariates), thresholds and regressions - probs is set from above
value=prev*(1-ascer)/ascer/(1-prev);
for(i=0;i<ns;i++){Ks[i]=value*probs[i]/(1+value*probs[i]-probs[i]);thresh[i]=normal_inv(1-Ks[i]);}

for(i=0;i<ns;i++)
{
value=pow(2*M_PI,-.5)*exp(-.5*pow(thresh[i],2))*(1-probs[i]*(ascer-prev)/ascer/(1-prev));
value2=pow(probs[i]*(1-probs[i]),.5)*(Ks[i]+(1-Ks[i])*prev*(1-ascer)/ascer/(1-prev));
regs[i]=value/value2;
}

if(Yadj!=NULL)	//set Yadj to (Y-probs)/root(probs*(1-probs)), divided by regressions
{
for(i=0;i<ns;i++){Yadj[i]=(Y[i]-probs[i])*pow(probs[i]*(1-probs[i]),-.5)/regs[i];}
}

if(covher!=NULL)	//get covher and topher
{
//get variance due to all fixed effects
sum=0;sum2=0;sumsq=0;sumsq2=0;count=0;count2=0;
for(i=0;i<ns;i++)
{
if(Y[i]==0){sum+=thresh[i];sumsq+=pow(thresh[i],2);count++;}
else{sum2+=thresh[i];sumsq2+=pow(thresh[i],2);count2++;}
}
var=sumsq/count-pow(sum/count,2);
var2=sumsq2/count2-pow(sum2/count2,2);
fixedvar=prev*(1-prev)*pow(sum2/count2-sum/count,2)+prev*var2+(1-prev)*var;

//get variance for just covariates
if(num_tops==0)	//already done
{covarvar=fixedvar;}
else	//must estimate coefficients for just covariates
{
thetas2=malloc(sizeof(double)*num_covars);
thetas2[0]=log(ascer/(1-ascer));
for(j=1;j<num_covars;j++){thetas2[j]=0;}

count=0;
while(1)
{
//get Zthetas then probs
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_covars, &alpha, Z, &ns, thetas2, &one, &beta, Zthetas, &one);
for(i=0;i<ns;i++){probs[i]=pow(1+exp(-Zthetas[i]),-1);}

//now likelihood and diff
like=0;for(i=0;i<ns;i++){like+=Y[i]*log(probs[i])+(1-Y[i])*log(1-probs[i]);}
if(count>1){diff2=diff;}
if(count>0){diff=like-likeold;}
likeold=like;

//AI = -2nd deriv = sum(pi(1-pi)ZiZi) and Bi = 1st deriv = sum(Yi-pi)Zi
for(j=0;j<num_covars;j++)
{
BI[j]=0;for(i=0;i<ns;i++){BI[j]+=(Y[i]-probs[i])*Z[i+j*ns];}
AI[j+j*num_covars]=0;for(i=0;i<ns;i++){AI[j+j*num_covars]+=probs[i]*(1-probs[i])*pow(Z[i+j*ns],2);}
for(j2=j+1;j2<num_covars;j2++)
{
AI[j+j2*num_covars]=0;for(i=0;i<ns;i++){AI[j+j2*num_covars]+=probs[i]*(1-probs[i])*Z[i+j*ns]*Z[i+j2*ns];}
AI[j2+j*num_covars]=AI[j+j2*num_covars];
}
}

//invert AI
(void)eigen_invert(AI, num_covars, AI2, -1, AI3, 1);

if(count>1)	//will have set diff and diff2 so can check if converged
{
if(fabs(diff)<tol&&fabs(diff2)<tol){break;}
}
if(count==maxiter){printf("Error, logistic regression did not converge after %d iterations; consider reducing the number of fixed effects, or using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %.6f)\n\n", maxiter, tol);exit(1);}

//update thetas using newton raphson
alpha=1.0;beta=1.0;
dgemv_("N", &num_covars, &num_covars, &alpha, AI, &num_covars, BI, &one, &beta, thetas2, &one);

count++;
}

//get Ks (population probs of being case) and corresponding thresholds
value=prev*(1-ascer)/ascer/(1-prev);
for(i=0;i<ns;i++){Ks[i]=value*probs[i]/(1+value*probs[i]-probs[i]);thresh[i]=normal_inv(1-Ks[i]);}

//get variance
sum=0;sum2=0;sumsq=0;sumsq2=0;count=0;count2=0;
for(i=0;i<ns;i++)
{
if(Y[i]==0){sum+=thresh[i];sumsq+=pow(thresh[i],2);count++;}
else{sum2+=thresh[i];sumsq2+=pow(thresh[i],2);count2++;}
}
var=sumsq/count-pow(sum/count,2);
var2=sumsq2/count2-pow(sum2/count2,2);
covarvar=prev*(1-prev)*pow(sum2/count2-sum/count,2)+prev*var2+(1-prev)*var;

free(thetas2);
}	//end of solving for just covariates

*covher=covarvar/(1+fixedvar);
*topher=(fixedvar-covarvar)/(1+fixedvar-covarvar);
}
}	//end of prev!=-9999

free(Zthetas);free(probs);free(AI);free(AI2);free(AI3);free(BI);
if(prev!=-9999){free(Ks);free(thresh);free(regs);}
}	//end of reg_covar_log

////////

void reg_single_log(double *Y, double *Z, int ns, int num_fixed, double *stats, double tol, int maxiter)
//regress Y on Z logistically
{
int i, j, j2, count, one=1;
double sum, ascer, like, likeold, diff, diff2, alpha, beta;
double *thetas, *Zthetas, *probs, *AI, *AI2, *AI3, *BI;


thetas=malloc(sizeof(double)*num_fixed);
Zthetas=malloc(sizeof(double)*ns);
probs=malloc(sizeof(double)*ns);
AI=malloc(sizeof(double)*num_fixed*num_fixed);
AI2=malloc(sizeof(double)*num_fixed);
AI3=malloc(sizeof(double)*num_fixed*num_fixed);
BI=malloc(sizeof(double)*num_fixed);

//get ascer
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
ascer=sum/ns;

//solve for all fixed - first covariate starts at ascer, rest are zero
thetas[0]=log(ascer/(1-ascer));
for(j=1;j<num_fixed;j++){thetas[j]=0;}

count=0;
while(1)
{
//get Zthetas then probs
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Zthetas, &one);
for(i=0;i<ns;i++){probs[i]=pow(1+exp(-Zthetas[i]),-1);}

//now likelihood and diff
like=0;for(i=0;i<ns;i++){like+=Y[i]*log(probs[i])+(1-Y[i])*log(1-probs[i]);}
if(count>1){diff2=diff;}
if(count>0){diff=like-likeold;}
likeold=like;

//AI = -2nd deriv = sum(pi(1-pi)ZiZi) and Bi = 1st deriv = sum(Yi-pi)Zi
for(j=0;j<num_fixed;j++)
{
BI[j]=0;for(i=0;i<ns;i++){BI[j]+=(Y[i]-probs[i])*Z[i+j*ns];}
AI[j+j*num_fixed]=0;for(i=0;i<ns;i++){AI[j+j*num_fixed]+=probs[i]*(1-probs[i])*pow(Z[i+j*ns],2);}
for(j2=j+1;j2<num_fixed;j2++)
{
AI[j+j2*num_fixed]=0;for(i=0;i<ns;i++){AI[j+j2*num_fixed]+=probs[i]*(1-probs[i])*Z[i+j*ns]*Z[i+j2*ns];}
AI[j2+j*num_fixed]=AI[j+j2*num_fixed];
}
}

//invert AI
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);

if(count>1)	//will have set diff and diff2 so can check if converged
{
if(fabs(diff)<tol&&fabs(diff2)<tol){break;}
}
if(count==maxiter){printf("Warning, logistic regression did not converge after %d iterations; consider reducing the number of fixed effects, or using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %.6f)\n\n", maxiter, tol);}

//update thetas using newton raphson
alpha=1.0;beta=1.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, AI, &num_fixed, BI, &one, &beta, thetas, &one);

count++;
}

stats[0]=thetas[num_fixed-1];
if(AI[num_fixed*num_fixed-1]>0)
{
stats[1]=pow(AI[num_fixed*num_fixed-1],.5);
stats[2]=stats[0]/stats[1];
stats[3]=erfc(fabs(stats[0]/stats[1])*M_SQRT1_2);
}
else{stats[1]=-9999;stats[2]=-9999;stats[3]=-9999;}

free(thetas);free(Zthetas);free(probs);free(AI);free(AI2);free(AI3);free(BI);
}	//end of reg_single_log

///////////////////////////

double get_YTdata(int j, double *residuals, int dtype, unsigned char **data_char, int nt, float *speedstarts, float *speedscales, double *centres)	//do not currently use
//could speed this up a bit if sure residuals have sum zero
{
int i, readint;
double sum;


sum=0;
if(dtype==1)	//bed format - actual values stored (3 means missing)
{
for(i=0;i<nt;i++)
{
readint=(int)((data_char[j][i/4] >> (2*(i%4))) & 3);
if(readint!=3){sum+=(readint-centres[j])*residuals[i];}
}
}
else	//short speed format
{
for(i=0;i<nt;i++)
{
readint=(int)data_char[j][i];
if(readint!=255){sum+=(speedstarts[j]+speedscales[j]*readint-centres[j])*residuals[i];}
}
}

return(sum);
}

////////

void update_residuals(int j, double change, double *residuals, int dtype, unsigned char **data_char, int nt, float *speedstarts, float *speedscales, double *centres)	//do not currently use
{
int i, readint;

if(change!=0)
{
if(dtype==1)	//bed format
{
for(i=0;i<nt;i++)
{
readint=(int)((data_char[j][i/4] >> (2*(i%4))) & 3);
if(readint!=3){residuals[i]-=(readint-centres[j])*change;}
}
}
else	//short speed format
{
for(i=0;i<nt;i++)
{
readint=(int)data_char[j][i];
if(readint!=255){residuals[i]-=(speedstarts[j]+speedscales[j]*readint-centres[j])*change;}
}
}
}
}

///////////////////////////

void set_lambdas(int p, double *lambdas, double *lambdas2, double *lambdas3, double *lambdas4, int length, double *exps, double varphen, double her, double *trylams, double *tryps, double *tryp2s, double *tryp3s, double *tryp4s, double *tryf2s, double enalpha, int type)
{
//type=1 - lasso-sparse, type=2 - lasso-shrink, type=3 - ridge, type=4 - bolt, type=5 bayesr-sparse, type=6 - bayesr-shrink
//type=7 - elastic net
//have already set all lambdas to zero
int j;
double value, value2;


if(type==1)	//lasso-sparse - lambda is set to match the lassosum paper
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[j+p*length]=trylams[p]*pow(length*exps[j],-.5);}
}
}

if(type==2)	//lasso-shrink - exp(betaj^2) = exps*varphen*her = 2/lambdas^2
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[j+p*length]=pow(2.0/exps[j]/her/varphen,.5);}
}
}

if(type==3)	//ridge - exp(betaj^2) = exps*varphen*her = lambdas
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[j+p*length]=exps[j]*her*varphen;}
}
}

if(type==4)	//bolt - exp(betaj^2) = exps*varphen*her = plambdas + (1-p)lambdas2 - might have f2=0, in which case lambdas2=0
{
value=(1-tryf2s[p])/tryps[p];
value2=tryf2s[p]/tryp2s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[j+p*length]=exps[j]*her*varphen*value;lambdas2[j+p*length]=exps[j]*her*varphen*value2;}
}
}

if(type==5)	//bayesr-sparse - exp(betaj^2) = exps*varphen*her = p2lambdas2 + p3lambdas3 + p4lambdas4, where lambdas2=lambdas4/100, lambdas3=lambdas4/10
{
value=tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0)
{
value2=exps[j]*her*varphen/value;
lambdas2[j+p*length]=value2/100;
lambdas3[j+p*length]=value2/10;
lambdas4[j+p*length]=value2;
}
}
}

if(type==6)	//bayesr-shrink - exp(betaj^2) = exps*varphen*her = plambdas + p2lambdas2 + p3lambdas3 + p4lambdas4
{
value=tryps[p]/1000+tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0)
{
value2=exps[j]*her*varphen/value;
lambdas[j+p*length]=value2/1000;
lambdas2[j+p*length]=value2/100;
lambdas3[j+p*length]=value2/10;
lambdas4[j+p*length]=value2;
}
}
}
}	//end of set_lambdas

///////////////////////////

double get_postmean(double sum, double lam, double lam2, double lam3, double lam4, double dsq, double resvar, double pp, double pp2, double pp3, double pp4, double *pen, int type, double *mcmc)
{
//type=1 - lasso-sparse, type=2 - lasso-shrink, type=3 - ridge, type=4 - bolt, type=5 bayesr-sparse, type=6 - bayesr-shrink
//when using gaussians, means are sum/A, vars are resvar/A, where A is (dsq+resvar/lam), mults are prop*sqrt(var/lam)*exp(mean^2/2var)

double mean, mean2, mean3, mean4, var, var2, var3, var4, frac, frac2, frac3, frac4, area, area2, trun, trun2, brun, brun2;
double max, value, value2, value3, value4, postmean, postvar, unifrand;

if(type==1)	//lasso with posterior mode
{
if(fabs(sum)>resvar*lam)
{
if(sum>0){postmean=(sum-resvar*lam)/dsq;}
else{postmean=(sum+resvar*lam)/dsq;}
}
else{postmean=0;}

//will not be using penalty, and can not perform mcmc
}

if(type==2)	//lasso with posterior mean - divide prior as positive exp then negative exp
{
//these are the means and variances of the non-truncated distributions
mean=(sum-resvar*lam)/dsq;
mean2=(sum+resvar*lam)/dsq;
var=resvar/dsq;
var2=resvar/dsq;

//for positive normal, area (denominator) is 1-Phi(-mean/var^.5) = Phi(mean/var^.5) = .5*erfc(-mean/var^.5 root(1/2))
//for negative normal, area (denominator) is Phi(-mean2/var2^.5) = .5*erfc(mean2/var2^.5 root(1/2))
area=.5*erfc(-mean*pow(var,-.5)*M_SQRT1_2);
area2=.5*erfc(mean2*pow(var2,-.5)*M_SQRT1_2);

//for means and variances, require the following
value=exp(-.5*pow(-mean,2)/var)*pow(2*M_PI,-.5);	//on wiki page, this is phi(alpha)
value2=exp(-.5*pow(-mean2,2)/var2)*pow(2*M_PI,-.5);	//on wiki page, this is phi(beta)

//get means of truncated distributions (if area is zero, values are irrelevant)
trun=mean+value/area*pow(var,.5);
trun2=mean2-value2/area2*pow(var2,.5);

//get variances of truncated distributions (again, values are irrelevant if area zero)
brun=var*(1-mean*pow(var,-.5)*value/area-pow(value/area,2));
brun2=var2*(1+mean2*pow(var2,-.5)*value2/area2-pow(-value2/area2,2));

max=pow(mean,2)/var;
if(pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
//here ignore pow(var,.5) in both lines
value=area*exp(.5*pow(mean,2)/var-.5*max);
value2=area2*exp(.5*pow(mean2,2)/var2-.5*max);

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*trun+frac2*trun2;
postvar=frac*(brun+pow(trun,2))+frac2*(brun2+pow(trun2,2))-pow(postmean,2);

if(postmean!=postmean||isinf(postmean))	//seems caused because lambda very large, so just set to zero
{postmean=0;}

//will not be using penalty

//have not worked out mcmc yet
}

if(type==3)	//ridge - just one gaussian
{
postmean=sum/(dsq+resvar/lam);
postvar=resvar/(dsq+resvar/lam);

if(postmean!=postmean){printf("sum %f dsp %f resvar %f lam %f\n", sum, dsq, resvar, lam);exit(1);}

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
*pen-=.5*(1+log(postvar/lam)-(postvar+pow(postmean,2))/lam);
}

if(mcmc!=NULL){*mcmc=postmean+rnorm_safe()*pow(postvar,.5);}
}

if(type==4)	//bolt - note that lam2 can be zero (but pp always positive)
{
if(lam2>0)	//have two gaussians
{
mean=sum/(dsq+resvar/lam);
mean2=sum/(dsq+resvar/lam2);
var=resvar/(dsq+resvar/lam);
var2=resvar/(dsq+resvar/lam2);

max=pow(mean,2)/var;
if(pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
value=pp*pow(var/lam,.5)*exp(.5*pow(mean,2)/var-.5*max);
value2=pp2*pow(var2/lam2,.5)*exp(.5*pow(mean2,2)/var2-.5*max);

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*mean+frac2*mean2;
postvar=frac*(var+pow(mean,2))+frac2*(var2+pow(mean2,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
*pen-=.5*frac2*(1+log(var2/lam2)-(var2+pow(mean2,2))/lam2);
}

if(mcmc!=NULL)
{
unifrand=(double)rand()/RAND_MAX;
if(unifrand<frac){*mcmc=mean+rnorm_safe()*pow(var,.5);}
else{*mcmc=mean2+rnorm_safe()*pow(var2,.5);}
}
}
else	//have gaussian and point mass
{
mean=sum/(dsq+resvar/lam);
var=resvar/(dsq+resvar/lam);

max=pow(mean,2)/var;
if(0>max){max=0;}
value=pp*pow(var/lam,.5)*exp(.5*pow(mean,2)/var-.5*max);
value2=pp2*exp(-.5*max);

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*mean;
postvar=frac*(var+pow(mean,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
}

if(mcmc!=NULL)
{
unifrand=(double)rand()/RAND_MAX;
if(unifrand<frac){*mcmc=mean+rnorm_safe()*pow(var,.5);}
else{*mcmc=0;}
}
}
}

if(type==5)	//bayesr with three gaussians and a point mass (pp, pp2, pp3 and pp4 can each be zero)
{
mean2=sum/(dsq+resvar/lam2);
mean3=sum/(dsq+resvar/lam3);
mean4=sum/(dsq+resvar/lam4);
var2=resvar/(dsq+resvar/lam2);
var3=resvar/(dsq+resvar/lam3);
var4=resvar/(dsq+resvar/lam4);

max=0;
if(pp2>0&&pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
if(pp3>0&&pow(mean3,2)/var3>max){max=pow(mean3,2)/var3;}
if(pp4>0&&pow(mean4,2)/var4>max){max=pow(mean4,2)/var4;}

value=0;value2=0;value3=0;value4=0;
if(pp>0){value=pp*exp(-.5*max);}
if(pp2>0){value2=pp2*pow(var2/lam2,.5)*exp(.5*pow(mean2,2)/var2-.5*max);}
if(pp3>0){value3=pp3*pow(var3/lam3,.5)*exp(.5*pow(mean3,2)/var3-.5*max);}
if(pp4>0){value4=pp4*pow(var4/lam4,.5)*exp(.5*pow(mean4,2)/var4-.5*max);}

frac=value/(value+value2+value3+value4);
frac2=value2/(value+value2+value3+value4);
frac3=value3/(value+value2+value3+value4);
frac4=value4/(value+value2+value3+value4);
postmean=frac2*mean2+frac3*mean3+frac4*mean4;
postvar=frac2*(var2+pow(mean2,2))+frac3*(var3+pow(mean3,2))+frac4*(var4+pow(mean4,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
if(frac3>0){*pen+=frac3*log(frac3/pp3);}
if(frac4>0){*pen+=frac4*log(frac4/pp4);}
*pen-=.5*frac2*(1+log(var2/lam2)-(var2+pow(mean2,2))/lam2);
*pen-=.5*frac3*(1+log(var3/lam3)-(var3+pow(mean3,2))/lam3);
*pen-=.5*frac4*(1+log(var4/lam4)-(var4+pow(mean4,2))/lam4);
}

if(mcmc!=NULL)
{
unifrand=(double)rand()/RAND_MAX;
if(unifrand<frac){*mcmc=0;}
if(unifrand>=frac&&unifrand<frac+frac2){*mcmc=mean2+rnorm_safe()*pow(var2,.5);}
if(unifrand>=frac+frac2&&unifrand<frac+frac2+frac3){*mcmc=mean3+rnorm_safe()*pow(var3,.5);}
if(unifrand>=frac+frac2+frac3){*mcmc=mean4+rnorm_safe()*pow(var4,.5);}
}
}

if(type==6)	//bayesr with four gaussians
{
mean=sum/(dsq+resvar/lam);
mean2=sum/(dsq+resvar/lam2);
mean3=sum/(dsq+resvar/lam3);
mean4=sum/(dsq+resvar/lam4);
var=resvar/(dsq+resvar/lam);
var2=resvar/(dsq+resvar/lam2);
var3=resvar/(dsq+resvar/lam3);
var4=resvar/(dsq+resvar/lam4);

//pp will always be >0
max=pow(mean,2)/var;
if(pp2>0&&pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
if(pp3>0&&pow(mean3,2)/var3>max){max=pow(mean3,2)/var3;}
if(pp4>0&&pow(mean4,2)/var4>max){max=pow(mean4,2)/var4;}

value=0;value2=0;value3=0;value4=0;
if(pp>0){value=pp*exp(-.5*max);}
if(pp2>0){value2=pp2*pow(var2/lam2,.5)*exp(.5*pow(mean2,2)/var2-.5*max);}
if(pp3>0){value3=pp3*pow(var3/lam3,.5)*exp(.5*pow(mean3,2)/var3-.5*max);}
if(pp4>0){value4=pp4*pow(var4/lam4,.5)*exp(.5*pow(mean4,2)/var4-.5*max);}

frac=value/(value+value2+value3+value4);
frac2=value2/(value+value2+value3+value4);
frac3=value3/(value+value2+value3+value4);
frac4=value4/(value+value2+value3+value4);
postmean=frac*mean+frac2*mean2+frac3*mean3+frac4*mean4;
postvar=frac*(var+pow(mean,2))+frac2*(var2+pow(mean2,2))+frac3*(var3+pow(mean3,2))+frac4*(var4+pow(mean4,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
if(frac3>0){*pen+=frac3*log(frac3/pp3);}
if(frac4>0){*pen+=frac4*log(frac4/pp4);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
*pen-=.5*frac2*(1+log(var2/lam2)-(var2+pow(mean2,2))/lam2);
*pen-=.5*frac3*(1+log(var3/lam3)-(var3+pow(mean3,2))/lam3);
*pen-=.5*frac4*(1+log(var4/lam4)-(var4+pow(mean4,2))/lam4);
}

if(mcmc!=NULL)
{
unifrand=(double)rand()/RAND_MAX;
if(unifrand<frac){*mcmc=mean+rnorm_safe()*pow(var,.5);}
if(unifrand>=frac&&unifrand<frac+frac2){*mcmc=mean2+rnorm_safe()*pow(var2,.5);}
if(unifrand>=frac+frac2&&unifrand<frac+frac2+frac3){*mcmc=mean3+rnorm_safe()*pow(var3,.5);}
if(unifrand>=frac+frac2+frac3){*mcmc=mean4+rnorm_safe()*pow(var4,.5);}
}
}

return(postmean);
}	//end of get_postmean

///////////////////////////

