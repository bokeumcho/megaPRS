/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Construct lasso, ridge, bolt and bayesr models from summary statistics - assumes phenotype has variance one

///////////////////////////

if(strcmp(sumsfile,"blank")!=0)	//neff is the average sample size for the main summaries
{
sum=0;for(j=0;j<data_length;j++){sum+=nss[j];}
neff=sum/data_length;
}

if(strcmp(pseudostem,"blank")!=0)	//neff2 is the average sample size for the training summaries
{
sum=0;for(j=0;j<data_length;j++){sum+=nss2[j];}
neff2=sum/data_length;
}

if(num_blocks>0)	//neff3 is the average sample size for jackknives - will already have neff
{
neff3=neff*(num_blocks-1)/num_blocks;
}

if(strcmp(pseudostem,"blank")!=0)	//set ldpreds, which indicates which predictors to use for testing
{
ldpreds=malloc(sizeof(int)*data_length);
for(j=0;j<data_length;j++){ldpreds[j]=1;}

if(strcmp(ldfile,"blank")!=0)	//exclude predictors in file
{
count=countrows(ldfile);
printf("Reading list of %d high-LD predictors from %s\n", count, ldfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(ldfile, wantpreds, count, NULL, 1, 0);

indexer=malloc(sizeof(int)*data_length);
count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Warning, none of the predictors are in the data\n\n");}
if(count2<count){printf("Warning, only %d of these are in the data\n\n", count2);}
if(count2==count){printf("All of these are in the data\n\n");}

for(j=0;j<count2;j++){ldpreds[indexer[j]]=0;}
}
}

//set num_try
if(strcmp(bestfile,"blank")!=0)	//considering a single model
{
num_try=1;
}
else	//considering multiple models
{
if(strcmp(fracfile,"blank")==0)	//defaults vary according to model
{
if(ptype==0)	//lasso-shrink, ridge, bolt and bayesr (must have multiher=1 and basicbolt=0)
{num_try=11+11+131+83;}
if(ptype==1)	//lasso-sparse
{num_try=80;}
if(ptype==2)	//lasso-shrink
{num_try=1+10*multiher;}
if(ptype==3)	//ridge
{num_try=1+10*multiher;}
if(ptype==4)
{num_try=132-114*basicbolt-110*ldpred;}
if(ptype==5)	//bayesr-sparse
{num_try=84;}
if(ptype==6)	//bayesr-shrink
{num_try=84;}
}
else	//will set models based on values in fracfile
{
num_try=countrows(fracfile);
}
}

//allocate parameters
trytypes=malloc(sizeof(int)*num_try);
tryhers=malloc(sizeof(double)*num_try);
trylams=malloc(sizeof(double)*num_try);
tryscales=malloc(sizeof(double)*num_try);
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

//set tryscales to 1, and hers, lams, ps & f2s to NA, then will change when required
for(p=0;p<num_try;p++)
{tryscales[p]=1.0;tryhers[p]=-9999;trylams[p]=-9999;tryps[p]=-9999;tryp2s[p]=-9999;tryp3s[p]=-9999;tryp4s[p]=-9999;tryf2s[p]=-9999;}

if(strcmp(bestfile,"blank")!=0)	//have already checked size of bestfile
{
//open bestfile, skip the header row and get model
if((input=fopen(bestfile,"r"))==NULL)
{printf("Error opening %s\n\n",bestfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
if(fscanf(input, "%d %s ", &readint, readstring)!=2)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
if(strcmp(readstring,"lasso-sparse")!=0&&strcmp(readstring,"lasso")!=0&&strcmp(readstring,"ridge")!=0&&strcmp(readstring,"bolt")!=0 &&strcmp(readstring,"bayesr")!=0&&strcmp(readstring,"bayesr-shrink")!=0)
{printf("Error, Element 2 of Row 2 of %s should be \"lasso\", \"ridge\", \"bolt\", \"bayesr\", \"bayesr-shrink\" or \"lasso-sparse\" (not %s), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, readstring);exit(1);}

if(strcmp(readstring,"lasso-sparse")==0)	//lasso-sparse - read lambda and s
{
trytypes[0]=1;
if(fscanf(input, "NA %lf %lf ", trylams, tryscales)!=2){printf("Error reading lambda and s from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"lasso")==0)	//lasso-shrink - read her scalings
{
trytypes[0]=2;
if(fscanf(input, "%lf ", tryhers)!=1){printf("Error reading heritability scaling from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"ridge")==0)	//ridge - read her scalings
{
trytypes[0]=3;
if(fscanf(input, "%lf ", tryhers)!=1){printf("Error reading heritability scaling from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"bolt")==0)	//bolt - read p and f2 (and set p2)
{
trytypes[0]=4;
if(fscanf(input, "%lf NA NA %lf %lf ", tryhers, tryps, tryf2s)!=3){printf("Error reading heritability scaling, p and f2 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
tryp2s[0]=1-tryps[0];
}

if(strcmp(readstring,"bayesr")==0)	//bayesr-sparse - read p, p2, p3 and p4
{
trytypes[0]=5;
if(fscanf(input, "%lf NA NA NA NA %lf %lf %lf %lf ", tryhers, tryps, tryp2s, tryp3s, tryp4s)!=5){printf("Error reading heritability scaling, p, p2, p3 and p4 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"bayesr-shrink")==0)	//bayesr-shrink - read p, p2, p3 and p4
{
trytypes[0]=6;
if(fscanf(input, "%lf NA NA NA NA %lf %lf %lf %lf ", tryhers, tryps, tryp2s, tryp3s, tryp4s)!=5){printf("Error reading heritability scaling, p, p2, p3 and p4 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

fclose(input);
}
else
{
if(strcmp(fracfile,"blank")==0)	//using defaults
{
count=0;

if(ptype==0)	//lasso-shrink, ridge, bolt and bayesr-sparse (must have multiher=1, basicbolt=0 and ldpred=0)
{
//lasso-shrink and ridge
for(p=5;p<16;p++){trytypes[count]=2;tryhers[count]=p*0.1;count++;}
for(p=5;p<16;p++){trytypes[count]=3;tryhers[count]=p*0.1;count++;}

//bolt - must skip ridge
loads=malloc(sizeof(double)*12);
loads2=malloc(sizeof(double)*11);
loads[0]=.5;loads[1]=.45;loads[2]=.4;loads[3]=.35;loads[4]=.3;loads[5]=.25;
loads[6]=.2;loads[7]=.15;loads[8]=.1;loads[9]=.05;loads[10]=0.02;loads[11]=0.01;
loads2[0]=.5;loads2[1]=.45;loads2[2]=.4;loads2[3]=.35;loads2[4]=.3;loads2[5]=.25;
loads2[6]=.2;loads2[7]=.15;loads2[8]=.1;loads2[9]=.05;loads2[10]=0;
for(j=0;j<12;j++){
for(j2=0;j2<11;j2++){
if(j+j2>0)
{trytypes[count]=4;tryhers[count]=1.0;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;}
}}
free(loads);free(loads2);

//bayesr-sparse - p4+p3+p2 must be positive
loads=malloc(sizeof(double)*7);
loads[0]=0;loads[1]=.005;loads[2]=.01;loads[3]=.02;loads[4]=.05;loads[5]=.1;loads[6]=.2;
for(j=0;j<7;j++){
for(j2=j;j2<7;j2++){
for(j3=j2;j3<7;j3++){
if(j+j2+j3>0)
{trytypes[count]=5;tryhers[count]=1.0;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;}
}}}
free(loads);
}

if(ptype==1)	//lasso-sparse
{
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.9;trylams[count]=0.001*pow(100,(float)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.5;trylams[count]=0.001*pow(100,(float)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.2;trylams[count]=0.001*pow(100,(float)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.1;trylams[count]=0.001*pow(100,(float)p/19);count++;}
}

if(ptype==2)	//lasso-shrink
{
if(multiher==1)	//ok if scaled her exceeds one
{
for(p=5;p<16;p++){trytypes[count]=2;tryhers[count]=p*0.1;count++;}
}
else{trytypes[count]=3;tryhers[count]=1.0;count++;}
}

if(ptype==3)	//ridge
{
if(multiher==1)	//ok if scaled her exceeds one
{
for(p=5;p<16;p++){trytypes[count]=3;tryhers[count]=p*0.1;count++;}
}
else{trytypes[count]=3;tryhers[count]=1.0;count++;}
}

if(ptype==4)	//bolt
{
if(basicbolt==0&&ldpred==0)	//132 pairs
{
loads=malloc(sizeof(double)*12);
loads2=malloc(sizeof(double)*11);
loads[0]=.5;loads[1]=.45;loads[2]=.4;loads[3]=.35;loads[4]=.3;loads[5]=.25;
loads[6]=.2;loads[7]=.15;loads[8]=.1;loads[9]=.05;loads[10]=0.02;loads[11]=0.01;
loads2[0]=.5;loads2[1]=.45;loads2[2]=.4;loads2[3]=.35;loads2[4]=.3;loads2[5]=.25;
loads2[6]=.2;loads2[7]=.15;loads2[8]=.1;loads2[9]=.05;loads2[10]=0;
for(j=0;j<12;j++){
for(j2=0;j2<11;j2++){
trytypes[count]=4;tryhers[count]=1.0;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}
else
{
if(basicbolt==1)	//just 18 pairs
{
loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
for(j=0;j<6;j++){
for(j2=0;j2<3;j2++){
trytypes[count]=4;tryhers[count]=1.0;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}
if(ldpred==1)	//just 22 pairs
{
loads=malloc(sizeof(double)*22);
loads[0]=1;loads[1]=.95;loads[2]=.9;loads[3]=.85;loads[4]=.8;loads[5]=.75;loads[6]=.7;loads[7]=.65;
loads[8]=.6;loads[9]=.55;loads[10]=.5;loads[11]=.45;loads[12]=.4;loads[13]=.35;loads[14]=.3;
loads[15]=.25;loads[16]=.2;loads[17]=.15;loads[18]=.1;loads[19]=.05;loads[20]=.02;loads[21]=.01;
for(j=0;j<22;j++)
{trytypes[count]=4;tryhers[count]=1.0;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=0;count++;}
free(loads);
}
}
}

if(ptype==5)	//bayesr-sparse (can not have all ps zero) - add in ridge model
{
loads=malloc(sizeof(double)*7);
loads[0]=0;loads[1]=.005;loads[2]=.01;loads[3]=.02;loads[4]=.05;loads[5]=.1;loads[6]=.2;
for(j=0;j<7;j++){
for(j2=j;j2<7;j2++){
for(j3=j2;j3<7;j3++){
if(j+j2+j3>0)
{trytypes[count]=5;tryhers[count]=1.0;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;}
}}}
free(loads);
trytypes[count]=5;tryhers[count]=1.0;tryps[count]=0;tryp2s[count]=0;tryp3s[count]=0;tryp4s[count]=1;count++;
}

if(ptype==6)	//bayesr-shrink - already includes ridge model
{
loads=malloc(sizeof(double)*7);
loads[0]=0;loads[1]=.005;loads[2]=.01;loads[3]=.02;loads[4]=.05;loads[5]=.1;loads[6]=.2;
for(j=0;j<7;j++){
for(j2=0;j2<7;j2++){
for(j3=0;j3<7;j3++){
if(j2>=j&&j3>=j2)
{
trytypes[count]=6;tryhers[count]=1.0;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}

if(count!=num_try){printf("Doug error %d and %d\n", num_try, count);exit(1);}
}	//end of setting using defaults
else	//read fracfile
{
for(p=0;p<num_try;p++){trytypes[p]=ptype;}

if(ptype==2)	//lasso-shrink - read her scalings
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0){printf("Error, heritability scaling in Row %d of %s is negative (%.4f)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}
}

if(ptype==3)	//ridge - read her
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0){printf("Error, heritability scaling in Row %d of %s is negative (%.4f)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}
}

if(ptype==4)	//bolt - read p and f2
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<=0){printf("Error, probability in Row %d of %s is non-positive (%.4f)\n\n", p+1, fracfile, tryps[p]);exit(1);}
if(tryps[p]>1){printf("Error, probability in Row %d of %s is greater than one (%.4f)\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
read_values(fracfile, tryf2s, num_try, NULL, 2, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryf2s[p]<0){printf("Error, fraction in Row %d of %s is negative (%.4f)\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
if(tryf2s[p]>=1){printf("Error, fraction in Row %d of %s not less than one (%.4f)\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryhers[p]=1.0;tryp2s[p]=1-tryps[p];}
}

if(ptype==5||ptype==6)	//bayes-sparse or bayesr-shrink - read p, p2, p3 and p4
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
read_values(fracfile, tryp2s, num_try, NULL, 2, 0, 0);
read_values(fracfile, tryp3s, num_try, NULL, 3, 0, 0);
read_values(fracfile, tryp4s, num_try, NULL, 4, 0, 0);

for(p=0;p<num_try;p++)
{
sum=tryps[p]+tryp2s[p]+tryp3s[p]+tryp4s[p];
if(sum<0.99||sum>1.01)
{printf("Error, values in Row %d of %s sum to %.4f (not one)\n\n", p+1, fracfile, sum);exit(1);}
tryps[p]=tryps[p]/sum;
tryp2s[p]=tryp2s[p]/sum;
tryp3s[p]=tryp3s[p]/sum;
tryp4s[p]=tryp4s[p]/sum;
}
for(p=0;p<num_try;p++){tryhers[p]=1.0;}
}
}	//end of reading fracfile
}

//save parameters
sprintf(filename,"%s.parameters",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Model\tType\tHer_Scaling\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\n");
for(p=0;p<num_try;p++)
{
if(trytypes[p]==1){fprintf(output, "%d\tlasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, trylams[p], tryscales[p]);}
if(trytypes[p]==2){fprintf(output, "%d\tlasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryhers[p]);}
if(trytypes[p]==3){fprintf(output, "%d\tridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryhers[p]);}
if(trytypes[p]==4){fprintf(output, "%d\tbolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\n", p+1, tryhers[p], tryps[p], tryf2s[p]);}
if(trytypes[p]==5){fprintf(output, "%d\tbayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryhers[p], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==6){fprintf(output, "%d\tbayesr-shrink\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryhers[p], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
}
fclose(output);

////////////////

//read headers from cors
usedpreds=malloc(sizeof(int)*num_preds);
maxnums=malloc(sizeof(int)*num_preds);
scumsums=malloc(sizeof(size_t)*num_preds);

sprintf(filename,"%s.cors.bin", corsfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read indicators and check compatible with keeppreds_use
fseeko(input, 0, SEEK_SET);
if(fread(usedpreds, sizeof(int), num_preds, input)!=num_preds)
{printf("Error reading predictor indicators from %s\n\n", filename);exit(1);}
for(j=0;j<data_length;j++)
{
if(usedpreds[keeppreds_use[j]]!=1){printf("Error, Predictor %s was not included when calculating predictor-predictor correlations; either remake the correlations or exclude predictors not included\n\n", preds[j]);exit(1);}
}

//read maxnums, and use to set cumsums
fseeko(input, (off_t)sizeof(int)*num_preds, SEEK_SET);
if(fread(maxnums, sizeof(int), num_preds, input)!=num_preds)
{printf("Error reading counts from %s\n\n", filename);exit(1);}
scumsums[0]=0;for(j=1;j<num_preds;j++){scumsums[j]=scumsums[j-1]+maxnums[j-1];}

fclose(input);

////////

//allocate variables
scount=0;for(j=0;j<data_length;j++){scount+=maxnums[keeppreds[j]];}
value=(double)scount/1024/1024/1024*8;
if(value>1){printf("Warning, to read the correlations requires %.1f Gb\n\n", value);}

//total is max of num_try and num_blocks
total=num_try;
if(num_blocks>num_try){total=num_blocks;}
anal_warn(data_length, num_blocks+4*num_try+3*total);

bigs=malloc(sizeof(int*)*data_length);
rjks=malloc(sizeof(float*)*data_length);
for(j=0;j<data_length;j++)
{
bigs[j]=malloc(sizeof(int)*maxnums[keeppreds_use[j]]);
rjks[j]=malloc(sizeof(float)*maxnums[keeppreds_use[j]]);
}

actnums=malloc(sizeof(int)*data_length);
rjksums=malloc(sizeof(float)*data_length);

datasqs=malloc(sizeof(double)*data_length);
YTdata=malloc(sizeof(double)*data_length);

if(num_blocks>0){YTdata2=malloc(sizeof(double)*data_length*num_blocks);}

exps=malloc(sizeof(double)*data_length);

lambdas=malloc(sizeof(double)*data_length*num_try);
lambdas2=malloc(sizeof(double)*data_length*num_try);
lambdas3=malloc(sizeof(double)*data_length*num_try);
lambdas4=malloc(sizeof(double)*data_length*num_try);

effs=malloc(sizeof(double)*data_length*total);
effs2=malloc(sizeof(double)*data_length*total);
convs=malloc(sizeof(int)*data_length*total);

if(strcmp(pseudostem,"blank")!=0){predcors=malloc(sizeof(double)*num_try);}

////////

//reopen cors
sprintf(filename,"%s.cors.bin", corsfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}

//read means, allowing for filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(centres+j, sizeof(double), 1, input)!=1)
{printf("Error reading mean for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}

//now scalings, allowing for filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(mults+j, sizeof(double), 1, input)!=1)
{printf("Error reading scaling for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}

//and variances, allowing for filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*2, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*2+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(sqdevs+j, sizeof(double), 1, input)!=1)
{printf("Error reading variance for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}

//finally read correlations - first with crude filtering
fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*3, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(maxnums[keeppreds_use[j]]>0)
{
if(keeppreds_use[j]!=current)
{fseeko(input, (off_t)sizeof(int)*num_preds*2+sizeof(double)*num_preds*3+(sizeof(int)+sizeof(float))*scumsums[keeppreds_use[j]], SEEK_SET);}
if(fread(bigs[j], sizeof(int), maxnums[keeppreds_use[j]], input)!=maxnums[keeppreds_use[j]])
{printf("Error reading indexes for Predictor %d from %s\n\n", j+1, filename);exit(1);}
if(fread(rjks[j], sizeof(float), maxnums[keeppreds_use[j]], input)!=maxnums[keeppreds_use[j]])
{printf("Error reading correlations for Predictor %d from %s\n\n", j+1, filename);exit(1);}
current=keeppreds_use[j]+1;
}
}
fclose(input);

//set actnums and rjksums, and replace rjks with shrink * covjk (this also finishes any filtering)
for(j=0;j<num_preds;j++){usedpreds[j]=-1;}
for(j=0;j<data_length;j++){usedpreds[keeppreds_use[j]]=j;}

for(j=0;j<data_length;j++)
{
count=0;
sum=1;
for(j2=0;j2<maxnums[keeppreds_use[j]];j2++)
{
j3=bigs[j][j2];
if(usedpreds[j3]!=-1)	//predictor remains (and has been mapped to usedpreds[j3])
{
sum+=pow(rjks[j][j2],2);
bigs[j][count]=usedpreds[j3];
rjks[j][count]=shrink*rjks[j][j2]*pow(sqdevs[j],.5)*pow(sqdevs[usedpreds[j3]],.5);
count++;
}
}
actnums[j]=count;
rjksums[j]=sum;
}

if(num_blocks>0)	//put jackknife XTY into YTdata2
{
//get indexes of predictors we want
sprintf(filename,"%s.jacks.predictors", jackstem);
count=countrows(filename);
wantpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*data_length);
datarands=malloc(sizeof(double)*count);

read_strings(filename, wantpreds, count, NULL, 1, 0);
count2=find_strings(preds, data_length, wantpreds, count, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, %s contains none of the %d predictors\n\n", filename, data_length);exit(1);}
if(count2<data_length){printf("Error, %s contains only %d of the %d predictors\n\n", filename, count2, data_length);exit(1);}

//check ns are reasonable
read_values(filename, datarands, count, NULL, 2, 0, 0);
for(j=0;j<data_length;j++)
{
if(fabs(datarands[indexer[j]]-nss[j])>.1)
{printf("Error, sample size for Predictor %s in %s (%.1f) does not match that in %s (%.1f)\n\n", preds[j], filename, datarands[indexer[j]], sumsfile, nss[j]);exit(1);}
}

//open bin
sprintf(filename2,"%s.jacks.bin", jackstem);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

//read rhos and convert to XTY (will already have sqdevs and neff3)
fseeko(input2, 0, SEEK_SET);
for(p=0;p<num_blocks;p++)
{
if(fread(datarands, sizeof(double), count, input2)!=count)
{printf("Error reading values from %s\n\n", filename2);exit(1);}
for(j=0;j<data_length;j++){YTdata2[(size_t)p*data_length+j]=datarands[indexer[j]]*pow(sqdevs[j],.5)*neff3;}
}
fclose(input2);

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);free(datarands);
}

////////

//exps should contain expected squared effect size of each predictor divided by hers
//betaj^2/hers = ej var(Y) / var(Xj), where ej is prop of her contributed by predictor j

//start by getting expected heritabilities (set to zero for trivial predictors)
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

//scale exps so they sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}

//now divide by (expected) predictor variance
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(hwestand==1){exps[j]*=pow(centres[j]*(1-centres[j]/2),-1);}
else{exps[j]*=pow(sqdevs[j],-1);}
}
}

//set lambdas
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++){lambdas[(size_t)p*data_length+j]=0;lambdas2[(size_t)p*data_length+j]=0;lambdas3[(size_t)p*data_length+j]=0;lambdas4[(size_t)p*data_length+j]=0;}

set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, 1.0, tryhers[p]*her, trylams, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, trytypes[p]);
}

//blank progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

///////////////////////

if(strcmp(pseudostem,"blank")!=0)	//solve using training summaries - var(Y)=1, so YTY=neff2 - then test using test summaries
{
//screen and file print
if(num_try==1){printf("Estimating effect sizes for one model using summary statistics in %s.train.summaries\n\n", pseudostem);}
else{printf("Estimating effect sizes for %d models using summary statistics in %s.train.summaries (if using multiple cores, models will finish in a random order)\n\n", num_try, pseudostem);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using summary statistics in %s.train.summaries\n", num_try, pseudostem);
fprintf(output, "Model\tType\tHeritability\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\tNum_Predictors_Failed\n");
fclose(output);

//get XTX and XTY for each predictor
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*neff2;YTdata[j]=rhos2[j]*pow(sqdevs[j],.5)*neff2;}

#pragma omp parallel for private(p,j,j2,count,sum,window_kb2,bitstart,bitstart2,bitend,sumsq,sumsq2) schedule(dynamic, 1)
for(p=0;p<num_try;p++)
{
//set starting effects
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
else	//effect sizes begin at marginal postmeans under own model, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL)/rjksums[j];}
}
else{effs[(size_t)p*data_length+j]=0;}
}

//set convs to zero
for(j=0;j<data_length;j++){convs[(size_t)p*data_length+j]=0;}

window_kb2=window_kb/segments;
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb2||chr[bitend]!=chr[bitstart]){break;}
}

//calculate ess for bit
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff2*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

//save effect sizes for bit, in case fail to converge
for(j=bitstart;j<bitend;j++){effs2[(size_t)p*data_length+j]=effs[(size_t)p*data_length+j];}

count=0;
while(1)
{
count++;

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j];
for(j2=0;j2<actnums[j];j2++){sum-=tryscales[p]*rjks[j][j2]*neff2*effs[(size_t)p*data_length+bigs[j][j2]];}

//update effect size
if(trytypes[p]==1){effs[(size_t)p*data_length+j]=get_postmean(sum, lambdas[(size_t)p*data_length+j]*neff2, -9999, -9999, -9999, datasqs[j], 1.0, -9999, -9999, -9999, -9999, NULL, 1, NULL);}	
else{effs[(size_t)p*data_length+j]=get_postmean(sum, lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL);}
}
}	//end of j loop

//save old then calculate new ess for bit
sumsq2=sumsq;
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff2*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

if(fabs(sumsq-sumsq2)/neff2<tol)	//converged
{
for(j=bitstart;j<bitend;j++){convs[(size_t)p*data_length+j]=1;}

if(window_kb2==window_kb)	//did a full window, so move start forwards window_kb/segments
{
for(bitstart2=bitstart+1;bitstart2<data_length;bitstart2++)
{
if(cmbp[bitstart2]-cmbp[bitstart]>1000*window_kb/segments||chr[bitstart2]!=chr[bitstart]){break;}
}
if(bitstart2==bitend)	//must have been at end of chromosome
{window_kb2=window_kb/segments;}
bitstart=bitstart2;
}
else	//did a small window, so start stays same, but increase window wise
{window_kb2*=2;}

break;
}

if(count==maxiter)	//failed (almost certainly due to end of bit) - undo iterations, skip rest of bit, and reduce window size 
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[(size_t)p*data_length+j];}
bitstart=bitend;
window_kb2=window_kb/segments;
break;
}
}	//end of inside while loop

}	//end of outside while loop

//count how many failed to converged
count=0;for(j=0;j<data_length;j++){count+=(exps[j]>0&&convs[(size_t)p*data_length+j]==0);}

//print update
#pragma omp critical
{
printf("Constructed Model %d: ", p+1);
if(trytypes[p]==1){printf("lasso-sparse, lambda %.4f, scale %.2f ", trylams[p], tryscales[p]);}
if(trytypes[p]==2){printf("lasso, heritability %.4f ", tryhers[p]*her);}
if(trytypes[p]==3){printf("ridge, heritability %.4f ", tryhers[p]*her);}
if(trytypes[p]==4){printf("bolt, heritability %.4f, p %.4f, f2 %.4f ", tryhers[p]*her, tryps[p], tryf2s[p]);}
if(trytypes[p]==5){printf("bayesr, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==6){printf("bayesr-shrink, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(count==0){printf(" - effect sizes converged for all predictors\n");}
else{printf(" - effect sizes failed to converge for %d predictors\n", count);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "%d\t", p+1);
if(trytypes[p]==1){fprintf(output, "lasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", trylams[p], tryscales[p], count);}
if(trytypes[p]==2){fprintf(output, "lasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, count);}
if(trytypes[p]==3){fprintf(output, "ridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, count);}
if(trytypes[p]==4){fprintf(output, "bolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, tryps[p], tryf2s[p], count);}
if(trytypes[p]==5){fprintf(output, "bayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], count);}
if(trytypes[p]==6){fprintf(output, "bayesr-shrink\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], count);}
fclose(output);
}
}	//end of p loop
printf("\n");

//test the models, and record the best (note that if exps[j]=0, the predictor will have effect zero)
printf("Testing the models using summary statistics in %s.test.summaries\n", pseudostem);

#pragma omp parallel for private(p,j,value,value2) schedule(dynamic, 1)
for(p=0;p<num_try;p++)
{
//get cov(score,Y)/SD(Y) = weighted sum of cov(Xj,Y)/SD(Y) - cov(Xj,Y)/SD(Y) = rhos[j]*SD(Xj)
value=0;
for(j=0;j<data_length;j++)
{
if(ldpreds[j]==1){value+=effs[(size_t)p*data_length+j]*rhos3[j]*pow(sqdevs[j],.5);}
}

//now get var(score) = sum betaj betak cov(Xj,Xk)
value2=0;
for(j=0;j<data_length;j++)
{
if(ldpreds[j]==1)
{
value2+=pow(effs[(size_t)p*data_length+j],2)*sqdevs[j];
for(j2=0;j2<actnums[j];j2++)
{
if(ldpreds[bigs[j][j2]]==1){value2+=effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]]*rjks[j][j2];}
}
}}

//compute correlation (if possible)
if(value2>0){predcors[p]=value*pow(value2,-.5);}
else{predcors[p]=-9999;}
}

//find best
best=-1;
for(p=0;p<num_try;p++)
{
if(predcors[p]!=-9999)
{
if(best==-1){best=p;value=predcors[p];}
if(predcors[p]>value){best=p;value=predcors[p];}
}
}

if(best==-1)	//not possible to compute a correlation for any models
{printf("Error, it was not possible to compute a correlation for any of the models (suggesting they all have zero effect sizes)\n\n");exit(1);}

//save correlations
sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Model\tCorrelation\n");

for(p=0;p<num_try;p++)
{
if(predcors[p]!=-9999){fprintf(output2,"%d\t%.4f\n", p+1, predcors[p]);}
else{fprintf(output2,"%d\tNA\n", p+1);}
}
fclose(output2);

//save the best model
sprintf(filename3,"%s.best",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Model\tType\tHer_Scaling\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\n");

if(trytypes[best]==1){fprintf(output3, "%d\tlasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, trylams[best], tryscales[best]);}
if(trytypes[best]==2){fprintf(output3, "%d\tlasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, tryhers[best]);}
if(trytypes[best]==3){fprintf(output3, "%d\tridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, tryhers[best]);}
if(trytypes[best]==4){fprintf(output3, "%d\tbolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\n", best+1, tryhers[best], tryps[best], tryf2s[best]);}
if(trytypes[best]==5){fprintf(output3, "%d\tbayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", best+1, tryhers[best], tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(trytypes[best]==6){fprintf(output3, "%d\tbayesr-shrink\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", best+1, tryhers[best], tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
fclose(output3);

printf("The estimated correlations of models are saved in %s, while the parameters corresponding to the best model are saved in %s\n\n", filename2, filename3);
}	//end of using training and test summary statistics

////////

if(strcmp(sumsfile,"blank")!=0)	//solve for main statistics - var(Y)=1, so YTY=neff
{
//set start and end
if(strcmp(pseudostem,"blank")!=0)	//test only best model
{start=best;end=best+1;}
else	//test all models
{start=0;end=num_try;}

//screen and file print
if(strcmp(pseudostem,"blank")!=0){printf("Estimating effect sizes using the best parameters and summary statistics in %s\n\n", sumsfile);}
else{printf("Estimating effect sizes for %d models using summary statistics in %s (if using multiple cores, models will finish in a random order)\n\n", end-start, sumsfile);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using summary statistics in %s\n", end-start, sumsfile);
fprintf(output, "Model\tType\tHeritability\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\tNum_Predictors_Failed\n");
fclose(output);

//get XTX and XTY for each predictor
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*neff;YTdata[j]=rhos[j]*pow(sqdevs[j],.5)*neff;}

#pragma omp parallel for private(p,j,j2,count,sum,window_kb2,bitstart,bitstart2,bitend,sumsq,sumsq2) schedule(dynamic, 1)
for(p=start;p<end;p++)
{
//set starting effects
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(trytypes[p]==1)	//lasso-sparse - set to zero(size_t)p*data_length+j
{effs[(size_t)p*data_length+j]=0;}
else	//effect sizes begin at marginal postmeans under own model, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL)/rjksums[j];}
}
else{effs[(size_t)p*data_length+j]=0;}
}

//set convs to zero
for(j=0;j<data_length;j++){convs[(size_t)p*data_length+j]=0;}

window_kb2=window_kb/segments;
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb2||chr[bitend]!=chr[bitstart]){break;}
}

//calculate ess for bit
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

//save effect sizes for bit, in case fail to converge
for(j=bitstart;j<bitend;j++){effs2[(size_t)p*data_length+j]=effs[(size_t)p*data_length+j];}

count=0;
while(1)
{
count++;

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j];
for(j2=0;j2<actnums[j];j2++){sum-=tryscales[p]*rjks[j][j2]*neff*effs[(size_t)p*data_length+bigs[j][j2]];}

//update effect size
if(trytypes[p]==1){effs[(size_t)p*data_length+j]=get_postmean(sum, lambdas[(size_t)p*data_length+j]*neff, -9999, -9999, -9999, datasqs[j], 1.0, -9999, -9999, -9999, -9999, NULL, 1, NULL);}
else
{effs[(size_t)p*data_length+j]=get_postmean(sum, lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL);
}
}
}	//end of j loop

//save old then calculate new ess for bit
sumsq2=sumsq;
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

if(fabs(sumsq-sumsq2)/neff<tol)	//converged
{
for(j=bitstart;j<bitend;j++){convs[(size_t)p*data_length+j]=1;}

if(window_kb2==window_kb)	//did a full window, so move start forwards window_kb/segments
{
for(bitstart2=bitstart+1;bitstart2<data_length;bitstart2++)
{
if(cmbp[bitstart2]-cmbp[bitstart]>1000*window_kb/segments||chr[bitstart2]!=chr[bitstart]){break;}
}
if(bitstart2==bitend)	//must have been at end of chromosome
{window_kb2=window_kb/segments;}
bitstart=bitstart2;
}
else	//did a small window, so start stays same, but increase window size
{window_kb2*=2;}

break;
}

if(count==maxiter)	//failed (almost certainly due to end of bit) - undo iterations, skip rest of bit, and reduce window size 
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[(size_t)p*data_length+j];}
bitstart=bitend;
window_kb2=window_kb/segments;
break;
}
}	//end of inside while loop

}	//end of outside while loop

//count how many failed to converged
count=0;for(j=0;j<data_length;j++){count+=(exps[j]>0&&convs[(size_t)p*data_length+j]==0);}

//print update
#pragma omp critical
{
printf("Constructed Model %d: ", p-start+1);
if(trytypes[p]==1){printf("lasso-sparse, lambda %.4f, scale %.4f ", trylams[p], tryscales[p]);}
if(trytypes[p]==2){printf("lasso, heritability %.4f ", tryhers[p]*her);}
if(trytypes[p]==3){printf("ridge, heritability %.4f ", tryhers[p]*her);}
if(trytypes[p]==4){printf("bolt, heritability %.4f, p %.4f, f2 %.4f ", tryhers[p]*her, tryps[p], tryf2s[p]);}
if(trytypes[p]==5){printf("bayesr, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==6){printf("bayesr-shrink, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(count==0){printf("- effect sizes converged for all predictors\n");}
else{printf("- effect sizes failed to converge for %d predictors\n", count);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "%d\t", p+1);
if(trytypes[p]==1){fprintf(output, "lasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\t%d\n", trylams[p], tryscales[p], count);}
if(trytypes[p]==2){fprintf(output, "lasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, count);}
if(trytypes[p]==3){fprintf(output, "ridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, count);}
if(trytypes[p]==4){fprintf(output, "bolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\t%d\n", tryhers[p]*her, tryps[p], tryf2s[p], count);}
if(trytypes[p]==5){fprintf(output, "bayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], count);}
if(trytypes[p]==6){fprintf(output, "bayesr-shrink\t%.4f\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], count);}
fclose(output);
}
}	//end of p loop
printf("\n");

//save models
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre");
for(p=start;p<end;p++){fprintf(output2," Model%d", p+1);}
fprintf(output2,"\n");
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=start;p<end;p++){fprintf(output2," %.4e", effs[(size_t)p*data_length+j]);}
fprintf(output2,"\n");
}
}
fclose(output2);

if(end-start==1){printf("Model saved in %s\n\n", filename2);}
else{printf("Models saved in %s\n\n", filename2);}
}	//end of using sumsfile

////////

if(num_blocks>0)	//now get jackknife models
{
//ensure best model is set - if using pseudos, will already be set, else will have num_try=1
if(strcmp(pseudostem,"blank")==0){best=0;}

//screen and file print
printf("Estimating effect sizes for %d jackknife models (if using multiple cores, models will finish in a random order)\n\n", num_blocks);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d jackknife models\n", num_blocks);
fprintf(output, "Model\tType\tHeritability\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\tNum_Predictors_Failed\n");
fclose(output);

//get XTX for each predictor (already have XTY)
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*neff3;}

#pragma omp parallel for private(p,j,j2,count,sum,window_kb2,bitstart,bitstart2,bitend,sumsq,sumsq2) schedule(dynamic, 1)
for(p=0;p<num_blocks;p++)
{
//set starting effects
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
else	//effect sizes begin at marginal postmeans under own model, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[(size_t)best*data_length+j], lambdas2[(size_t)best*data_length+j], lambdas3[(size_t)best*data_length+j], lambdas4[(size_t)best*data_length+j], datasqs[j], 1.0, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], NULL, trytypes[best], NULL)/rjksums[j];}
}
else{effs[(size_t)p*data_length+j]=0;}
}

//set convs to zero
for(j=0;j<data_length;j++){convs[(size_t)p*data_length+j]=0;}

window_kb2=window_kb/segments;
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb2||chr[bitend]!=chr[bitstart]){break;}
}

//calculate ess for bit
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata2[(size_t)p*data_length+j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff3*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

//save effect sizes for bit, in case fail to converge
for(j=bitstart;j<bitend;j++){effs2[(size_t)p*data_length+j]=effs[(size_t)p*data_length+j];}

count=0;
while(1)
{
count++;

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata2[(size_t)p*data_length+j];
for(j2=0;j2<actnums[j];j2++){sum-=tryscales[best]*rjks[j][j2]*neff3*effs[(size_t)p*data_length+bigs[j][j2]];}

//update effect size
if(trytypes[best]==1){effs[(size_t)p*data_length+j]=get_postmean(sum, lambdas[(size_t)best*data_length+j]*neff3, -9999, -9999, -9999, datasqs[j], 1.0, -9999, -9999, -9999, -9999, NULL, 1, NULL);}
else
{effs[(size_t)p*data_length+j]=get_postmean(sum, lambdas[(size_t)best*data_length+j], lambdas2[(size_t)best*data_length+j], lambdas3[(size_t)best*data_length+j], lambdas4[(size_t)best*data_length+j], datasqs[j], 1.0, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], NULL, trytypes[best], NULL);}
}
}	//end of j loop

//save old then calculate new ess for bit
sumsq2=sumsq;
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata2[(size_t)p*data_length+j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff3*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

if(fabs(sumsq-sumsq2)/neff3<tol)	//converged
{
for(j=bitstart;j<bitend;j++){convs[(size_t)p*data_length+j]=1;}

if(window_kb2==window_kb)	//did a full window, so move start forwards window_kb/segments
{
for(bitstart2=bitstart+1;bitstart2<data_length;bitstart2++)
{
if(cmbp[bitstart2]-cmbp[bitstart]>1000*window_kb/segments||chr[bitstart2]!=chr[bitstart]){break;}
}
if(bitstart2==bitend)	//must have been at end of chromosome
{window_kb2=window_kb/segments;}
bitstart=bitstart2;
}
else	//did a small window, so start stays same, but increase window size
{window_kb2*=2;}

break;
}

if(count==maxiter)	//failed (almost certainly due to end of bit) - undo iterations, skip rest of bit, and reduce window size 
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[(size_t)p*data_length+j];}
bitstart=bitend;
window_kb2=window_kb/segments;
break;
}
}	//end of inside while loop
}	//end of outside while loop

//count how many failed to converged
count=0;for(j=0;j<data_length;j++){count+=(exps[j]>0&&convs[(size_t)p*data_length+j]==0);}

//print update
#pragma omp critical
{
printf("Constructed Model %d: ", p+1);
if(trytypes[best]==1){printf("lasso-sparse, lambda %.4f, scale %.4f ", trylams[best], tryscales[best]);}
if(trytypes[best]==2){printf("lasso, heritability %.4f ", tryhers[best]*her);}
if(trytypes[best]==3){printf("ridge, heritability %.4f ", tryhers[best]*her);}
if(trytypes[best]==4){printf("bolt, heritability %.4f, p %.4f, f2 %.4f ", tryhers[best]*her, tryps[best], tryf2s[best]);}
if(trytypes[best]==5){printf("bayesr, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[best]*her, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(trytypes[best]==6){printf("bayesr-shrink, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[best]*her, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(count==0){printf("- effect sizes converged for all predictors\n");}
else{printf("- effect sizes failed to converge for %d predictors\n", count);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "%d\t", p+1);
if(trytypes[best]==1){fprintf(output, "lasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\t%d\n", trylams[best], tryscales[best], count);}
if(trytypes[best]==2){fprintf(output, "lasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[best]*her, count);}
if(trytypes[best]==3){fprintf(output, "ridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[best]*her, count);}
if(trytypes[best]==4){fprintf(output, "bolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\t%d\n", tryhers[best]*her, tryps[best], tryf2s[best], count);}
if(trytypes[best]==5){fprintf(output, "bayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[best]*her, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], count);}
if(trytypes[best]==6){fprintf(output, "bayesr-shrink\t%.4f\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[best]*her, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], count);}
fclose(output);
}
}	//end of p loop
printf("\n");

//save models
sprintf(filename2,"%s.jacks",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre");
for(p=0;p<num_blocks;p++){fprintf(output2," Model%d", p+1);}
fprintf(output2,"\n");
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=0;p<num_blocks;p++){fprintf(output2," %.4e", effs[(size_t)p*data_length+j]);}
fprintf(output2,"\n");
}
}
fclose(output2);

printf("Models saved in %s\n\n", filename2);
}	//end of num_blocks>0

////////

if(strcmp(pseudostem,"blank")!=0){free(ldpreds);}
free(trytypes);free(tryhers);free(trylams);free(tryscales);free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);free(tryf2s);
free(usedpreds);free(maxnums);free(scumsums);
for(j=0;j<data_length;j++){free(bigs[j]);free(rjks[j]);}free(bigs);free(rjks);
free(actnums);free(rjksums);
free(datasqs);free(YTdata);
if(num_blocks>0){free(YTdata2);}
free(exps);
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(effs);free(effs2);free(convs);
if(strcmp(pseudostem,"blank")!=0){free(predcors);}

///////////////////////////

//made a basic mcmc solver
/*
if(mcmc>0)	//try mcmc - must have one model, and for moment will only have mcmc=1
{
//ensure best model is set - if using pseudos, will already be set, else will have num_try=1
if(strcmp(pseudostem,"blank")==0){best=0;}

//set start and end
start=best;end=best+1;

//screen and file print
if(end-start==1){printf("Estimating effect sizes via MCMC for one model using summary statistics in %s\n\n", sumsfile);}
else{printf("Estimating effect sizes via MCMC for %d models using summary statistics in %s (if using multiple cores, models will finish in a random order)\n\n", end-start, sumsfile);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models via MCMC using summary statistics in %s\n", end-start, sumsfile);
fprintf(output, "Model\tType\tHeritability\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\tNum_Predictors_Failed\n");
fclose(output);

//get XTX and XTY for each predictor
for(j=0;j<data_length;j++){datasqs[j]=sqdevs[j]*neff;YTdata[j]=rhos[j]*pow(sqdevs[j],.5)*neff;}

#pragma omp parallel for private(p,j,j2,count,sum,window_kb2,bitstart,bitstart2,bitend,sumsq,sumsq2) schedule(dynamic, 1)
for(p=start;p<end;p++)
{
//set starting effects - can not have lasso
for(j=0;j<data_length;j++)
{
if(exps[j]>0){effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL)/rjksums[j];}
else{effs[(size_t)p*data_length+j]=0;}
}

//set convs to zero
for(j=0;j<data_length;j++){convs[(size_t)p*data_length+j]=0;}

//will store final estimates in effs2
for(j=0;j<data_length;j++){effs2[(size_t)p*data_length+j]=0;}

//will always move window_kb
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb||chr[bitend]!=chr[bitstart]){break;}
}

for(count=0;count<100+500;count++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j];
for(j2=0;j2<actnums[j];j2++){sum-=tryscales[p]*rjks[j][j2]*neff*effs[(size_t)p*data_length+bigs[j][j2]];}

//update effect size - can not have lasso
value=get_postmean(sum, lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], effs+(size_t)p*data_length+j);

if(count>100)	//saving effect
{effs2[(size_t)p*data_length+j]+=value/500;}
}
}	//end of j loop

//calculate ess for bit
sumsq=0;
for(j=bitstart;j<bitend;j++)
{
sumsq+=2*YTdata[j]*effs[(size_t)p*data_length+j]-datasqs[j]*pow(effs[(size_t)p*data_length+j],2);
for(j2=0;j2<actnums[j];j2++)
{
if(bigs[j][j2]>=bitstart&&bigs[j][j2]<bitend)
{sumsq-=rjks[j][j2]*neff*effs[(size_t)p*data_length+j]*effs[(size_t)p*data_length+bigs[j][j2]];}
}
}

if(sumsq<0||sumsq>neff)	//bit has failed - effects return to starting values
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0){effs2[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[(size_t)p*data_length+j], lambdas2[(size_t)p*data_length+j], lambdas3[(size_t)p*data_length+j], lambdas4[(size_t)p*data_length+j], datasqs[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL)/rjksums[j];}
else{effs2[(size_t)p*data_length+j]=0;}

convs[(size_t)p*data_length+j]=1;
}
break;
}

}	//end of count loop

bitstart=bitend;
}	//end of outside while loop

//count how many failed to converged
count=0;for(j=0;j<data_length;j++){count+=(exps[j]>0&&convs[(size_t)p*data_length+j]==0);}

//print update
#pragma omp critical
{
printf("Constructed Model %d: ", p+1);
if(trytypes[p]==1){printf("lasso-sparse, lambda %.4f, scale %.4f ", trylams[p], tryscales[p]);}
if(trytypes[p]==2){printf("lasso, heritability %.4f ", tryhers[p]*her);}
if(trytypes[p]==3){printf("ridge, heritability %.4f ", tryhers[p]*her);}
if(trytypes[p]==4){printf("bolt, heritability %.4f, p %.4f, f2 %.4f ", tryhers[p]*her, tryps[p], tryf2s[p]);}
if(trytypes[p]==5){printf("bayesr, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==6){printf("bayesr-shrink, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f ", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(count==0){printf("- effect sizes converged for all predictors\n");}
else{printf("- effect sizes failed to converge for %d predictors\n", count);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "%d\t", p+1);
if(trytypes[p]==1){fprintf(output, "lasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\t%d\n", trylams[p], tryscales[p], count);}
if(trytypes[p]==2){fprintf(output, "lasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, count);}
if(trytypes[p]==3){fprintf(output, "ridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d\n", tryhers[p]*her, count);}
if(trytypes[p]==4){fprintf(output, "bolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\t%d\n", tryhers[p]*her, tryps[p], tryf2s[p], count);}
if(trytypes[p]==5){fprintf(output, "bayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], count);}
if(trytypes[p]==6){fprintf(output, "bayesr-shrink\t%.4f\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", tryhers[p]*her, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], count);}
fclose(output);
}
}	//end of p loop
printf("\n");

//save models
sprintf(filename2,"%s.mcmc",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre");
for(p=start;p<end;p++){fprintf(output2," Model%d", p+1);}
fprintf(output2,"\n");
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=start;p<end;p++){fprintf(output2," %.4e", effs2[(size_t)p*data_length+j]);}
fprintf(output2,"\n");
}
}
fclose(output2);

if(end-start==1){printf("Model saved in %s\n\n", filename2);}
else{printf("Models saved in %s\n\n", filename2);}
}	//end of mcmc
*/

