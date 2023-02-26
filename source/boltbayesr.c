/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Perform ridge, bolt and bayesr

///////////////////////////

//set num_try
if(strcmp(fracfile,"blank")==0)	//defaults vary according to model
{
if(mode==151)	//ridge
{num_try=1+4*multiher;}
if(mode==152)	//bolt
{
if(ldpred==0){num_try=18;}
else{num_try=22;}
}
if(mode==153)	//bayesr
{
if(fullspace==0){num_try=35;}
else{num_try=125;}
}
}
else	//will set models based on values in fracfile
{
num_try=countrows(fracfile);
}

//allocate parameters
tryhers=malloc(sizeof(double)*num_try);
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

//set tryhers to 1, and other values to NA, then will change when required
for(p=0;p<num_try;p++){tryhers[p]=1.0;tryps[p]=-9999;tryp2s[p]=-9999;tryp3s[p]=-9999;tryp4s[p]=-9999;tryf2s[p]=-9999;}

if(strcmp(fracfile,"blank")==0)	//using defaults
{
if(mode==151)	//ridge
{
if(multiher==1)	//make sure scaled her does not exceed one
{
count=0;
for(p=8;p<12;p++)
{
tryhers[count]=p*0.1;
if(tryhers[count]*her>=1){tryhers[count]=0.99/her;}
count++;
}
}
else{tryhers[0]=1.0;}
}

if(mode==152)	//bolt
{
if(ldpred==0)
{
loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
count=0;
for(j=0;j<6;j++)
{
for(j2=0;j2<3;j2++){tryps[count]=loads[j];tryp2s[count]=1-loads[count];tryf2s[count]=loads2[j2];count++;}
}
free(loads);free(loads2);
}
else
{
loads=malloc(sizeof(double)*21);
loads[0]=.95;loads[1]=.9;loads[2]=.85;loads[3]=.8;loads[4]=.75;loads[5]=.7;loads[6]=.65;
loads[7]=.6;loads[8]=.55;loads[9]=.5;loads[10]=.45;loads[11]=.4;loads[12]=.35;loads[13]=.3;
loads[14]=.25;loads[15]=.2;loads[16]=.15;loads[17]=.1;loads[18]=.05;loads[19]=.02;loads[20]=.01;
tryps[0]=.5;tryp2s[0]=.5;tryf2s[0]=.5;
count=1;
for(j=0;j<21;j++){tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=0;count++;}
free(loads);
}
}

if(mode==153)	//bayesr - when pointmass=1 (sparse version), not allowed to have ps4, ps3 and ps2 all zero
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;
count=0;
for(j=0;j<5;j++){
for(j2=0;j2<5;j2++){
for(j3=0;j3<5;j3++){
if((j+j2+j3>0||pointmass==0)&&((j2>=j&&j3>=j2)||fullspace==1))
{
tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
tryps[count]=0;tryp2s[count]=0;tryp3s[count]=0;tryp4s[count]=1;count++;
free(loads);
}
}
else	//read fracfile
{
if(mode==151)	//ridge - read her scaling
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0){printf("Error, heritability scaling in Row %d of %s is negative (%.4f)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
if(tryhers[p]*her>=1){tryhers[p]=0.99/her;}
}
}

if(mode==152)	//bolt - read p and f2
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
for(p=0;p<num_try;p++){tryp2s[p]=1-tryps[p];}
}

if(mode==153)	//bayes-sparse or bayesr-shrink - read p, p2, p3 and p4
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
}
}	//end of reading fracfile

//save parameters
sprintf(filename2,"%s.parameters",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(mode==151)
{
fprintf(output2, "Model\tHer_Scaling\n");
for(p=0;p<num_try;p++){fprintf(output2, "%d\t%.4f\n", p+1, tryhers[p]);}
}
if(mode==152)
{
fprintf(output2, "Model\tp\tf2\n");
for(p=0;p<num_try;p++){fprintf(output2, "%d\t%.4f\t%.4f\n", p+1, tryps[p], tryf2s[p]);}
}
if(mode==153)
{
fprintf(output2, "Model\tp1\tp2\tp3\n");
for(p=0;p<num_try;p++){fprintf(output2, "%d\t%.4f\t%.4f\t%.4f\n", p+1, tryp4s[p], tryp3s[p], tryp2s[p]);}
}
fclose(output2);

////////

//restage indicates which stage we start on (-1 - cv, -2 - many finals, >=0 - best final)
//recount is number of iterations at start

if(restart==0)	//normal start - make a blank save file, and set restage and recount
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=-1;fwrite(&readint, sizeof(int), 1, output2);
readint=0;fwrite(&readint, sizeof(int), 1, output2);
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);

if(skipcv==0){restage=-1;}
else{restage=-2;}
recount=0;
}
else	//restarting from previous run
{
//set readstring, for now and later
if(mode==151){sprintf(readstring,"\"--ridge\"");}
if(mode==152){sprintf(readstring,"\"--bolt\"");}
if(mode==153){sprintf(readstring,"\"--bayesr\"");}

printf("Reading details from previous run\n");

//check .save file exists and read restage and recount
sprintf(filename2,"%s.save.bin", outfile);
if(just_check(filename2)!=0)
{printf("Error reading %s; make sure the output filename is the same as when you first ran %s\n\n", filename2, readstring);exit(1);}

if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&restage, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&recount, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(restage==-1&&skipcv==1){printf("Error, the previous run performed cross-validation; make sure the cross-validation options are the same as when you first ran %s\n\n", readstring);exit(1);}
if(restage==-2&&skipcv==0){printf("Error, the previous run did not perform cross-validation; make sure the cross-validation options are the same as when you first ran %s\n\n", readstring);exit(1);}

if(recount==0){printf("Error, the previous run had not completed at least ten iterations (if running on a cluster, make sure to allocate sufficient memory and time)\n\n");exit(1);}

//check size
if(restage==-1)	//solving for num_try models
{
fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use)*num_try)
{printf("Error, %s should have size %jd (not %jd); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, (off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use)*num_try, ftello(input2), readstring);exit(1);}
}
if(restage>0)	//solving for 1+num_blocks models
{
fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use)*(1+num_blocks))
{printf("Error, %s should have size %jd (not %jd); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, (off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use), ftello(input2), readstring);exit(1);}
}
if(restage==-2)	//solving for num_try models
{
fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use)*num_try)
{printf("Error, %s should have size %jd (not %jd); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, (off_t)sizeof(int)*2+sizeof(double)*(1+data_length+num_samples_use), ftello(input2), readstring);exit(1);}
}
fclose(input2);

if(restage>=0)	//set best
{best=restage;}

printf("Will be restarting from Iteration %d\n\n", recount+1);
}

////////

//set num_train, num_test and keepboth (indexes training then test samples, relative to ids)
keepboth=malloc(sizeof(int)*num_samples_use);
keepboth2=malloc(sizeof(int)*num_samples_use);

if(skipcv==0)	//sort cv samples
{
if(cvprop!=-9999)
{
num_test=cvprop*num_samples_use;
num_train=num_samples_use-num_test;

if(restart==0)	//pick samples at random then save
{
for(i=0;i<num_samples_use;i++){keepboth[i]=i;}
permute_int(keepboth,num_samples_use);

sprintf(filename2,"%s.cv.samples",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(i=num_train;i<num_samples_use;i++){fprintf(output2,"%s %s\n", ids1[keepboth[i]], ids2[keepboth[i]]);}
fclose(output2);

sprintf(filename2,"%s.cv.index",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(i=0;i<num_samples_use;i++){fprintf(output2,"%d\n", keepboth[i]);}
fclose(output2);
}
else	//read samples used in first run - will have already set readstring
{
sprintf(filename2,"%s.cv.index", outfile);
if(just_check(filename2)!=0)
{printf("Error reading %s; make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, readstring);exit(1);}
if(countrows(filename2)!=num_samples_use)
{printf("Error, %s should have %d rows (not %d); make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename2, num_test, countrows(filename2), readstring);exit(1);}

if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{
if(fscanf(input2,"%d ", keepboth+i)!=1)
{printf("Error reading row %d of %s\n\n", i+1, filename2);exit(1);}
}
fclose(input2);
}
}	//end of cvprop!=-9999

if(strcmp(bvsfile,"blank")!=0)
{
count=countrows(bvsfile);
printf("Reading list of %d cross-validation samples from %s\n", count, bvsfile);
wantids=malloc(sizeof(char*)*count);
read_ids(bvsfile, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
num_test=find_strings(ids3, num_samples_use, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(num_test==0){printf("Error, can not find any of these samples in the data\n\n");exit(1);}
if(num_test<count){printf("Warning, only %d of these are in the data\n", num_test);}
num_train=num_samples_use-num_test;

usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(i=0;i<num_test;i++){usedids[indexer[i]]=1;}
count2=0;
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==0){keepboth[count2]=i;count2++;}
}
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==1){keepboth[count2]=i;count2++;}
}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);free(usedids);
}

printf("Will be using %d samples to train and %d to test\n\n", num_train, num_test);
if(num_train<3){printf("Error, unable to continue with fewer than three training samples\n\n");exit(1);}
if(num_test<3){printf("Error, unable to continue with fewer than three test samples\n\n");exit(1);}
}
else	//will not use num_train or num_test
{
num_train=num_samples_use;
num_test=0;
for(i=0;i<num_samples_use;i++){keepboth[i]=i;}
}

//set keepboth2 (indexes all samples)
for(i=0;i<num_samples_use;i++){keepboth2[i]=keepsamps[keepboth[i]];}

////////

//allocate other variables
if(dtype==1){count=(num_samples_use-1)/4+1;}
else{count=num_samples_use;}
value=(double)data_length/1024/1024/1024*count;
if(value>1){printf("Warning, to read the data requires %.1f Gb (this can not be reduced)\n\n", value);}

//total is max of num_try and 1+num_blocks
total=num_try;
if(1+num_blocks>num_try){total=1+num_blocks;}
anal_warn(data_length, 4*num_try+2*total);

data_char=malloc(sizeof(unsigned char*)*data_length);
for(j=0;j<data_length;j++){data_char[j]=malloc(sizeof(unsigned char)*count);}
if(dtype==4){speedstarts=malloc(sizeof(float)*data_length);speedscales=malloc(sizeof(float)*data_length);}

Y=malloc(sizeof(double)*num_samples_use);
Yadj=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
datasqs=malloc(sizeof(double)*data_length);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

exps=malloc(sizeof(double)*data_length);

lambdas=malloc(sizeof(double)*data_length*num_try);
lambdas2=malloc(sizeof(double)*data_length*num_try);
lambdas3=malloc(sizeof(double)*data_length*num_try);
lambdas4=malloc(sizeof(double)*data_length*num_try);

effs=malloc(sizeof(double)*data_length*total);
residuals=malloc(sizeof(double)*num_samples_use*total);

pens=malloc(sizeof(double)*total);
likes=malloc(sizeof(double)*total);
likesold=malloc(sizeof(double)*total);

if(num_blocks>0){order=malloc(sizeof(int)*num_samples_use);}

//variables for speeding up code
anal_warn(bitsize,(num_samples_use+data_length));
data=malloc(sizeof(double)*(num_samples_use*bitsize+4));
cors=malloc(sizeof(double)*bitsize*data_length);
YTdata=malloc(sizeof(double)*bitsize*total);
changes=malloc(sizeof(double)*bitsize*total);
if(dtype==1)
{
value=(double)data_length/1024/1024/1024*8*256*4;
if(value>1){printf("Warning, to create lookup tables requires %.1f Gb (this can not be reduced)\n\n", value);}

bytelookup=malloc(sizeof(double*)*data_length);
for(j=0;j<data_length;j++){bytelookup[j]=malloc(sizeof(double)*256*4);}
}

//will be reading data in bits
bittotal=(data_length-1)/bitsize+1;

////////

//read in (all) the data - will be using either bed format or short SPEED format (dtype 1 or 4)
if(dtype==1)
{(void)read_bed_full(datafile, data_char, num_samples_use, keepboth2, data_length, keeppreds_use, num_samples, num_preds, maxthreads);}
else
{(void)read_speed_full(datafile, speedstarts, speedscales, data_char, NULL, num_samples_use, keepboth2, data_length, keeppreds_use, num_samples, num_preds, nonsnp, maxthreads);}

//fill Y and Z, get thetas and residuals
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[keepboth[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[keepboth[i]+j*num_samples_use];}
}
reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);

//save
sprintf(filename2,"%s.coeff", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Component Effect SD P\n");
fprintf(output2, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output2, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output2);

//get variance of residuals (Yadj will have mean zero)
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(Yadj[i],2);}
varphen=sumsq/num_samples_use;

////////

//get means and variances - mults indicates which predictors non-trivial (will subsequently use exps)
printf("Computing predictor means and variances\n\n");
wcount=0;
#pragma omp parallel for private(j,sum,sumsq,indcount,i,readint,mean,var) schedule (static)
for(j=0;j<data_length;j++)
{
if(dtype==1)	//bed format - actual values stored (3 means missing)
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
readint=(int)((data_char[j][i/4] >> (2*(i%4))) & 3);
if(readint!=3){sum+=readint;sumsq+=pow(readint,2);indcount++;}
}
}
else	//short speed format - values are speedstart + value * speedscale
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
readint=(int)data_char[j][i];
if(readint!=255){sum+=readint;sumsq+=pow(readint,2);indcount++;}
}
sumsq=indcount*pow(speedstarts[j],2)+2*speedstarts[j]*speedscales[j]*sum+pow(speedscales[j],2)*sumsq;
sum=indcount*speedstarts[j]+speedscales[j]*sum;
}

if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=-9999;var=0;}

centres[j]=mean;
if(var>0){mults[j]=1;}
else	//trivial
{
mults[j]=-9999;
#pragma omp critical
{
if(wcount<5){printf("Warning, Predictor %s is trivial (takes at most one non-missing value) and will be ignored\n", preds[j]);}
wcount++;
}
}
sqdevs[j]=var*indcount/num_samples_use;
}	//end of j loop
if(wcount==data_length){printf("Error, all predictors are trivial\n\n");exit(1);}
if(wcount>5){printf("In total, %d predictors are trivial\n", wcount);}
if(wcount>0){printf("\n");}

if(dtype==1)	//get look up for each predictor
{
printf("Constructing lookup table for each predictor\n\n");
#pragma omp parallel for private(j,k,readint) schedule (static)
for(j=0;j<data_length;j++)
{
for(k=0;k<256;k++)
{
readint=(k & 3);
if(readint!=3){bytelookup[j][k*4]=readint-centres[j];}
else{bytelookup[j][k*4]=0;}
readint=((k >> 2) & 3);
if(readint!=3){bytelookup[j][k*4+1]=readint-centres[j];}
else{bytelookup[j][k*4+1]=0;}
readint=((k >> 4) & 3);
if(readint!=3){bytelookup[j][k*4+2]=readint-centres[j];}
else{bytelookup[j][k*4+2]=0;}
readint=((k >> 6) & 3);
if(readint!=3){bytelookup[j][k*4+3]=readint-centres[j];}
else{bytelookup[j][k*4+3]=0;}
}
}
}

////////

//exps should contain expected squared effect size of each predictor divided by (varphen * hers)
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
for(j=0;j<data_length;j++){lambdas[j+p*data_length]=0;lambdas2[j+p*data_length]=0;lambdas3[j+p*data_length]=0;lambdas4[j+p*data_length]=0;}

if(mode==151)	//ridge
{set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, varphen, tryhers[p]*her, NULL, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, 3);}
if(mode==152)	//bolt
{set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, varphen, tryhers[p]*her, NULL, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, 4);}
if(mode==153)	//bayesr
{set_lambdas(p, lambdas, lambdas2, lambdas3, lambdas4, data_length, exps, varphen, tryhers[p]*her, NULL, tryps, tryp2s, tryp3s, tryp4s, tryf2s, -9999, 5+(pointmass==0));}
}

//blank or check progress file
sprintf(filename,"%s.progress",outfile);
if(restart==0)
{
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);
}
else
{
if(just_check(filename)!=0)
{printf("Error reading %s; make sure the options you use now are the same as when you first ran %s (except for adding \"--restart YES\")\n\n", filename, readstring);exit(1);}
}

////////

if(restage==-1)	//make and test training models
{
if(num_try==1)	//just need to set effect sizes to zero and best
{
printf("Can skip training phase because there is only set of parameters\n\n");
for(j=0;j<data_length;j++){effs[j]=0;}
for(i=0;i<num_samples_use;i++){residuals[i]=Yadj[i];}
best=0;
}
else
{
//get predictor-predictor covariances, from which can extract XTX
printf("Computing predictor-predictor covariances across training samples\n\n");
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_train, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+(j-bitstart)];}
}

if(restart==0)	//start with effect sizes zero
{
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++){effs[j+p*data_length]=0;}
for(i=0;i<num_samples_use;i++){residuals[i+p*num_samples_use]=Yadj[i];}
}
}
else	//read likelihoods, effect sizes and residuals from file
{
sprintf(filename2,"%s.save.bin", outfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(fread(likes, sizeof(double), num_try, input2)!=num_try)
{printf("Error reading model likelihoods from %s\n\n", filename2);exit(1);}
for(p=0;p<num_try;p++)
{
if(fread(effs+p*data_length, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading effects sizes for Model %d from %s\n\n", p+1, filename2);exit(1);}
if(fread(residuals+p*num_samples_use, sizeof(double), num_samples_use, input2)!=num_samples_use)
{printf("Error reading residuals for Model %d from %s\n\n", p+1, filename2);exit(1);}
}
fclose(input2);
}

//screen and file print (no need to reprint headers if restart=1)
printf("Estimating effect sizes for %d models using training samples\nIter\tNum_Con\tMax_Diff\tTarget\n", num_try);

if(restart==0)
{
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using training samples\n", num_try);
fprintf(output,"Iter\tNum_Converged\tMax_Difference\tTarget\n");
fclose(output);
}

//iterate effect sizes
count=recount;
while(1)
{
count++;

for(p=0;p<num_try;p++){pens[p]=0;}	//use this to keep track of likelihood terms
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//get t(X) residuals
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_try, &num_train, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitsize);

#pragma omp parallel for private(p,j,j2,sum,postmean) schedule(dynamic, 1)
for(p=0;p<num_try;p++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XTresiduals
sum=effs[j+p*data_length]*datasqs[j]+YTdata[(j-bitstart)+p*bitsize];

if(mode==151)	//ridge
{postmean=get_postmean(sum, lambdas[j+p*data_length], -9999, -9999, -9999, datasqs[j], (1-tryhers[p]*her)*varphen, -9999, -9999, -9999, -9999, pens+p, 3, NULL);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], -9999, -9999, datasqs[j], (1-tryhers[p]*her)*varphen, tryps[p], tryp2s[p], -9999, -9999, pens+p, 4, NULL);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], (1-tryhers[p]*her)*varphen, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], pens+p,  5+(pointmass==0), NULL);}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitsize]=postmean-effs[j+p*data_length];
effs[j+p*data_length]=postmean;
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitsize]-=changes[j-bitstart+p*bitsize]*cors[(size_t)j*bitsize+j2-bitstart];}
}
}	//end of j loop
}	//end of p loop

//update residuals (for all samples)
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_try, &bitlength, &alpha, data, &num_samples_use, changes, &bitsize, &beta, residuals, &num_samples_use);
}	//end of bit loop

for(p=0;p<num_try;p++)	//approx likelihoods
{
if(count>1){likesold[p]=likes[p];}
sumsq=0;for(i=0;i<num_train;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
likes[p]=-.5*num_train*log(2*M_PI*(1-tryhers[p]*her)*varphen)-.5*sumsq/(1-tryhers[p]*her)/varphen-pens[p];
}

if(count%10==0)	//save in case need to restart
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=-1;fwrite(&readint, sizeof(int), 1, output2);
readint=count;fwrite(&readint, sizeof(int), 1, output2);
fwrite(likes, sizeof(double), num_try, output2);
for(p=0;p<num_try;p++)
{
fwrite(effs+p*data_length, sizeof(double), data_length, output2);
fwrite(residuals+p*num_samples_use, sizeof(double), num_samples_use, output2);
}
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);
}

if(count==1)	//just print update
{
printf("%d\t0\tNA\t%.4f\n", count, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t0\tNA\t%.4f\n", count, tol);
fclose(output);
}
else	//get number converged, largest difference, print update, then see if breaking
{
cflag=0;
for(p=0;p<num_try;p++)
{
cflag+=(fabs(likes[p]-likesold[p])<tol);
if(p==0){diff=fabs(likes[p]-likesold[p]);}
if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

printf("%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);
fclose(output);

if(diff<tol){break;}
}

if(count>=maxiter){printf("\nWarning, Variational Bayes did not converge after %d iterations (%d out of %d converged, largest difference in likelihood %f, tolerance %f); consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance\n", count, cflag, num_try, diff, tol);break;}
}	//end of while loop
printf("\n");

//measure and save mse, recording best
printf("Measuring accuracy of each model\n");
sprintf(filename2,"%s.mse",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(mode==151){fprintf(output2,"Model\tHer_Scaling\tMean_Squared_Error\tApprox_Likelihood\n");}
if(mode==152){fprintf(output2,"Model\tp\tf2\tMean_Squared_Error\tApprox_Likelihood\n");}
if(mode==153){fprintf(output2,"Model\tp1\tp2\tp3\tMean_Squared_Error\tApprox_Likelihood\n");}

for(p=0;p<num_try;p++)
{
sumsq=0;for(i=num_train;i<num_samples_use;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
value=sumsq/num_test;

if(mode==151)
{printf("Model %d: heritability %.4f, mean squared error %.4f\n", p+1, tryhers[p]*her, value);
fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\n", p+1, tryhers[p], value, likes[p]);}
if(mode==152)
{printf("Model %d: p %.4f, f2 %.4f, mean squared error %.4f\n", p+1, tryps[p], tryf2s[p], value);
fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryps[p], tryf2s[p], value, likes[p]);}
if(mode==153)
{printf("Model %d: p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f, mean squared error %.4f\n", p+1, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], value);
fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], value, likes[p]);}

if(p==0){value2=value;best=0;}
if(value<value2){value2=value;best=p;}
}
printf("\n");
fclose(output2);

//measure and save correlation (currently only for interest, to see how it compares to mse)
sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(mode==151){fprintf(output2,"Model\tHer_Scaling\tCorrelation\tProp_Variance_Exp\n");}
if(mode==152){fprintf(output2,"Model\tp\tf2\tCorrelation\tProp_Variance_Exp\n");}
if(mode==153){fprintf(output2,"Model\tp1\tp2\tp3\tCorrelation\tProp_Variance_Exp\n");}

for(p=0;p<num_try;p++)
{
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;value2=0;
for(i=num_train;i<num_samples_use;i++)
{
sum+=Yadj[i]-residuals[i+p*num_samples_use];
sum2+=Yadj[i];
sumsq+=pow(Yadj[i]-residuals[i+p*num_samples_use],2);
sumsq2+=pow(Yadj[i],2);
sumsq3+=(Yadj[i]-residuals[i+p*num_samples_use])*Yadj[i];
value2+=pow(residuals[i+p*num_samples_use],2);
}
mean=sum/num_train;mean2=sum2/num_train;
value=(sumsq3-num_train*mean*mean2)/pow(sumsq-num_train*mean*mean,.5)/pow(sumsq2-num_train*mean2*mean2,.5);

if(mode==151){fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\n", p+1, tryhers[p], value, 1-value2/sumsq2);}
if(mode==152){fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryps[p], tryf2s[p], value, 1-value2/sumsq2);}
if(mode==153){fprintf(output2,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], value, 1-value2/sumsq2);}
}
printf("\n");
fclose(output2);
}
}	//end of making and testing training models

////////

if(skipcv==0)	//make one model (and maybe jackknife models) using all samples - will already have best
{
total2=1+num_blocks;
if(num_blocks>0)	//get num_train2 and initialize order
{
num_train2=(int)(subprop*num_samples_use);
for(i=0;i<num_samples_use;i++){order[i]=i;}
}

//get predictor-predictor correlations, from which can extract XTX
printf("Computing predictor-predictor covariances across all samples\n\n");
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+(j-bitstart)];}
}

if(restage==-1)
{
for(p=0;p<total2;p++)	//set effects and residuals for all models to those from best model
{
for(j=0;j<data_length;j++){effs[j+p*data_length]=effs[j+best*data_length];}
for(i=0;i<num_samples_use;i++){residuals[i+p*num_samples_use]=residuals[i+best*num_samples_use];}
}

for(p=1;p<total2;p++)	//adjust residuals of jackknife models (zero some phenotypes, scale remainder)
{
permute_int(order,num_samples_use);
for(i=0;i<num_train2;i++){residuals[order[i]+p*num_samples_use]+=Yadj[order[i]]*(1-subprop)/subprop;}
for(i=num_train2;i<num_samples_use;i++){residuals[order[i]+p*num_samples_use]-=Yadj[order[i]];}
}
}
else	//read effect sizes, residuals and likelihood from file
{
sprintf(filename2,"%s.save.bin", outfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(total2==1)
{
if(fread(likes, sizeof(double), 1, input2)!=1)
{printf("Error reading likelihood from %s\n\n", filename2);exit(1);}
if(fread(effs, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading effects sizes from %s\n\n", filename2);exit(1);}
if(fread(residuals, sizeof(double), num_samples_use, input2)!=num_samples_use)
{printf("Error reading residuals from %s\n\n", filename2);exit(1);}
}
else
{
if(fread(likes, sizeof(double), total2, input2)!=total2)
{printf("Error reading likelihoods from %s\n\n", filename2);exit(1);}
for(p=0;p<total2;p++)
{
if(fread(effs+p*data_length, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading effects sizes for Model %d from %s\n\n", p+1, filename2);exit(1);}
if(fread(residuals+p*num_samples_use, sizeof(double), num_samples_use, input2)!=num_samples_use)
{printf("Error reading residuals for Model %d from %s\n\n", p+1, filename2);exit(1);}
}
}
fclose(input2);
}

//screen and file print (no need to reprint headers if restage!=1)
if(mode==151){printf("Estimating effect sizes for best-fitting model (heritability %.4f)", tryhers[best]*her);}
if(mode==152){printf("Estimating effect sizes for best-fitting model (p %.4f, f2 %.4f)", tryps[best], tryf2s[best]);}
if(mode==153){printf("Estimating effect sizes for best-fitting model (p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f)", tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(total2==1){printf("\nIter\tDiff\tTarget\n");}
else{printf(", and for %d jackknife models\nIter\tNum_Con\tDiff\tTarget\n", num_blocks);}

if(restage==-1)
{
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}

if(total2==1){fprintf(output,"Constructing best-fitting model using all samples\nIter\tConverged\tDifference\tTarget\n");}
else{fprintf(output,"Constructing best-fitting model and %d jackknife models using all samples\nIter\tNum_Converged\tMax_Difference\tTarget\n", num_blocks);}
fclose(output);
}

//iterate effect sizes
if(restage==-1){count=0;}
else{count=recount;}
while(1)
{
count++;

for(p=0;p<total2;p++){pens[p]=0;}
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//get t(X) residuals
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total2, &num_samples_use, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitsize);

#pragma omp parallel for private(p,j,j2,sum,postmean) schedule(dynamic, 1)
for(p=0;p<total2;p++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XTresiduals
sum=effs[j+p*data_length]*datasqs[j]+YTdata[(j-bitstart)+p*bitsize];

if(p==0)	//normal model
{
if(mode==151)	//ridge
{postmean=get_postmean(sum, lambdas[j+best*data_length], -9999, -9999, -9999, datasqs[j], (1-tryhers[best]*her)*varphen, -9999, -9999, -9999, -9999, pens+p, 3, NULL);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[j+best*data_length], lambdas2[j+best*data_length], -9999, -9999, datasqs[j], (1-tryhers[best]*her)*varphen, tryps[best], tryp2s[best], -9999, -9999, pens+p, 4, NULL);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[j+best*data_length], lambdas2[j+best*data_length], lambdas3[j+best*data_length], lambdas4[j+best*data_length], datasqs[j], (1-tryhers[best]*her)*varphen, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], pens+p,  5+(pointmass==0), NULL);}
}
else	//values are too high by 1/subprop
{
if(mode==151)	//ridge
{postmean=get_postmean(sum*subprop, lambdas[j+best*data_length], -9999, -9999, -9999, datasqs[j]*subprop, (1-tryhers[best]*her)*varphen, -9999, -9999, -9999, -9999, pens+p, 3, NULL);}
if(mode==152)	//bolt
{postmean=get_postmean(sum*subprop, lambdas[j+best*data_length], lambdas2[j+best*data_length], -9999, -9999, datasqs[j]*subprop, (1-tryhers[best]*her)*varphen, tryps[best], tryp2s[best], -9999, -9999, pens+p, 4, NULL);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum*subprop, lambdas[j+best*data_length], lambdas2[j+best*data_length], lambdas3[j+best*data_length], lambdas4[j+best*data_length], datasqs[j]*subprop, (1-tryhers[best]*her)*varphen, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], pens+p,  5+(pointmass==0), NULL);}
}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitsize]=postmean-effs[j+p*data_length];
effs[j+p*data_length]=postmean;
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitsize]-=changes[j-bitstart+p*bitsize]*cors[(size_t)j*bitsize+j2-bitstart];}
}
}	//end of j loop
}	//end of p loop

//update residuals
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total2, &bitlength, &alpha, data, &num_samples_use, changes, &bitsize, &beta, residuals, &num_samples_use);
}	//end of bit loop

for(p=0;p<total2;p++)	//approx likelihoods
{
if(count>1){likesold[p]=likes[p];}
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
likes[p]=-.5*num_samples_use*log(2*M_PI*(1-tryhers[best]*her)*varphen)-.5*sumsq/(1-tryhers[best]*her)/varphen-pens[p];
}

if(count%10==0)	//save in case need to restart
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=best;fwrite(&readint, sizeof(int), 1, output2);
readint=count;fwrite(&readint, sizeof(int), 1, output2);
fwrite(likes, sizeof(double), total2, output2);
for(p=0;p<total2;p++)
{
fwrite(effs+p*data_length, sizeof(double), data_length, output2);
fwrite(residuals+p*num_samples_use, sizeof(double), num_samples_use, output2);
}
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);
}

if(count==1)	//just print update
{
if(total2==1){printf("%d\tNA\t%.4f\n", count, tol);}
else{printf("%d\t0\tNA\t%.4f\n", count, tol);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(total2==1){fprintf(output,"%d\tNA\t%.4f\n", count, tol);}
else{fprintf(output, "%d\t0\tNA\t%.4f\n", count, tol);}
fclose(output);
}
else	//get number converged, difference for first model, print update, then see if breaking
{
cflag=0;
for(p=0;p<total2;p++)
{
cflag+=(fabs(likes[p]-likesold[p])<tol);
if(p==0){diff=fabs(likes[p]-likesold[p]);}
//if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

if(total2==1){printf("%d\t%.4f\t%.4f\n", count, diff, tol);}
else{printf("%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(total2==1){fprintf(output,"%d\t%.4f\t%.4f\n", count, diff, tol);}
else{fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);}
fclose(output);

if(diff<tol){break;}
}

if(count>=maxiter){printf("\nWarning, Variational Bayes did not converge after %d iterations (%d out of %d converged, largest difference in likelihood %f, tolerance %f); consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance\n", count, cflag, total2, diff, tol);break;}
}	//end of while loop
printf("\n");

//save effects
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect");
for(p=0;p<num_blocks;p++){fprintf(output2," Jack%d", p+1);}
fprintf(output2,"\n");

for(j=0;j<num_tops;j++)
{
fprintf(output2, "%s %c %c %.6f", tpreds[j], tal1[j], tal2[j], tcentres[j]);
for(p=0;p<total2;p++){fprintf(output2," %.4e", thetas[num_covars+j]);}
fprintf(output2,"\n");
}
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=0;p<total2;p++){fprintf(output2," %.4e", effs[j+p*data_length]);}
fprintf(output2,"\n");
}
}
fclose(output2);

printf("Best-fitting model saved in %s\n\n", filename2);
}	//end of skipcv=0

////////

if(skipcv==1)	//make all models using all samples
{
//get predictor-predictor correlations, from which can extract XTX
printf("Computing predictor-predictor covariances across all samples\n\n");
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+(j-bitstart)];}
}

if(restart==0)	//start with effect sizes zero
{
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++){effs[j+p*data_length]=0;}
for(i=0;i<num_samples_use;i++){residuals[i+p*num_samples_use]=Yadj[i];}
}
}
else	//read likelihoods, effect sizes and residuals from file
{
sprintf(filename2,"%s.save.bin", outfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading first value of %s\n\n", filename2);exit(1);}
if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading second value of %s\n\n", filename2);exit(1);}

if(fread(likes, sizeof(double), num_try, input2)!=num_try)
{printf("Error reading model likelihoods from %s\n\n", filename2);exit(1);}
for(p=0;p<num_try;p++)
{
if(fread(effs+p*data_length, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading effects sizes for Model %d from %s\n\n", p+1, filename2);exit(1);}
if(fread(residuals+p*num_samples_use, sizeof(double), num_samples_use, input2)!=num_samples_use)
{printf("Error reading residuals for Model %d from %s\n\n", p+1, filename2);exit(1);}
}
fclose(input2);
}

//screen and file print (no need to reprint headers if restart=1)
printf("Estimating effect sizes for %d models using all samples\nIter\tNum_Con\tMax_Diff\tTarget\n", num_try);

if(restart==0)
{
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using training samples\n", num_try);
fprintf(output,"Iter\tNum_Converged\tMax_Difference\tTarget\n");
fclose(output);
}

//iterate effect sizes
if(restart==0){count=0;}
else{count=recount;}
while(1)
{
count++;

for(p=0;p<num_try;p++){pens[p]=0;}
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//get t(X) residuals
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_try, &num_samples_use, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitsize);

#pragma omp parallel for private(p,j,j2,sum,postmean) schedule(dynamic, 1)
for(p=0;p<num_try;p++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XTresiduals
sum=effs[j+p*data_length]*datasqs[j]+YTdata[(j-bitstart)+p*bitsize];

if(mode==151)	//ridge
{postmean=get_postmean(sum, lambdas[j+p*data_length], -9999, -9999, -9999, datasqs[j], (1-tryhers[p]*her)*varphen, -9999, -9999, -9999, -9999, pens+p, 3, NULL);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], -9999, -9999, datasqs[j], (1-tryhers[p]*her)*varphen, tryps[p], tryp2s[p], -9999, -9999, pens+p, 4, NULL);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[j+p*data_length], lambdas2[j+p*data_length], lambdas3[j+p*data_length], lambdas4[j+p*data_length], datasqs[j], (1-tryhers[p]*her)*varphen, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], pens+p,  5+(pointmass==0), NULL);}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitsize]=postmean-effs[j+p*data_length];
effs[j+p*data_length]=postmean;
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitsize]-=changes[j-bitstart+p*bitsize]*cors[(size_t)j*bitsize+j2-bitstart];}
}
}	//end of j loop
}	//end of p loop

//update residuals
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_try, &bitlength, &alpha, data, &num_samples_use, changes, &bitsize, &beta, residuals, &num_samples_use);
}	//end of bit loop

for(p=0;p<num_try;p++)	//approx likelihoods
{
if(count>1){likesold[p]=likes[p];}
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(residuals[i+p*num_samples_use],2);}
likes[p]=-.5*num_samples_use*log(2*M_PI*(1-tryhers[p]*her)*varphen)-.5*sumsq/(1-tryhers[p]*her)/varphen-pens[p];
}

if(count%10==0)	//save in case need to restart
{
sprintf(filename2,"%s.save.part", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readint=-2;fwrite(&readint, sizeof(int), 1, output2);
readint=count;fwrite(&readint, sizeof(int), 1, output2);
fwrite(likes, sizeof(double), num_try, output2);
for(p=0;p<num_try;p++)
{
fwrite(effs+p*data_length, sizeof(double), data_length, output2);
fwrite(residuals+p*num_samples_use, sizeof(double), num_samples_use, output2);
}
fclose(output2);
sprintf(cmd, "mv %s.save.part %s.save.bin", outfile, outfile);
system(cmd);
}

if(count==1)	//just print update
{
printf("%d\t0\tNA\t%.4f\n", count, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t0\tNA\t%.4f\n", count, tol);
fclose(output);
}
else	//get number converged, largest difference, print update, then see if breaking
{
cflag=0;
for(p=0;p<num_try;p++)
{
cflag+=(fabs(likes[p]-likesold[p])<tol);
if(p==0){diff=fabs(likes[p]-likesold[p]);}
if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

printf("%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol);
fclose(output);

if(diff<tol){break;}
}

if(count>=maxiter){printf("\nWarning, Variational Bayes did not converge after %d iterations (%d out of %d converged, largest difference in likelihood %f, tolerance %f); consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance\n", count, cflag, num_try, diff, tol);break;}
}	//end of while loop
printf("\n");

//save effects
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre");
for(p=0;p<num_try;p++){fprintf(output2," Effect%d", p+1);}
fprintf(output2,"\n");

for(j=0;j<num_tops;j++)
{
fprintf(output2, "%s %c %c %.6f", tpreds[j], tal1[j], tal2[j], tcentres[j]);
for(p=0;p<num_try;p++){fprintf(output2," %.4e", thetas[num_covars+j]);}
fprintf(output2,"\n");
}
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=0;p<num_try;p++){fprintf(output2," %.4e", effs[j+p*data_length]);}
fprintf(output2,"\n");
}
}
fclose(output2);

printf("Models saved in %s\n\n", filename2);
}	//end of skipcv=0

////////

free(tryhers);free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);free(tryf2s);
free(keepboth);free(keepboth2);
for(j=0;j<data_length;j++){free(data_char[j]);}free(data_char);
if(dtype==4){free(speedstarts);free(speedscales);}
free(Y);free(Yadj);free(Z);free(datasqs);
free(thetas);free(thetasds);free(thetapvas);
free(exps);
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(effs);free(residuals);
free(pens);free(likes);free(likesold);
free(data);free(cors);free(YTdata);free(changes);
if(dtype==1){for(j=0;j<data_length;j++){free(bytelookup[j]);}free(bytelookup);}

///////////////////////////

