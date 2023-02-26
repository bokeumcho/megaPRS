/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for reading stuff
//general syntax - thing to read, objects to be filled, length, keeps, other stuff req to read, misc

///////////////////////////

int read_bimfile(char *bimfile, int *chr, char **preds, double *cm, double *bp, char *al1, char *al2, int length, int *keeppreds_use, int type, double window_cm, int flag)
//type=0 - keep quiet, type=1 - warn about everything
//flag=0 - positions must be in order, flag=1 - do not
{
int j, count, found, want, bcount, lcount, mcount, ncount, scount;

char *rc, *rs, *rm, *rbp, *ra1, *ra2;

FILE *input;

rc=malloc(sizeof(char)*10000000);rs=malloc(sizeof(char)*10000000);
rm=malloc(sizeof(char)*10000000);rbp=malloc(sizeof(char)*10000000);
ra1=malloc(sizeof(char)*10000000);ra2=malloc(sizeof(char)*10000000);


if(countcols(bimfile)!=6)
{printf("Error, %s should have six columns (not %d); columns should be chr, name, genetic distance, bp, A1, A2\n\n", bimfile, countcols(bimfile));exit(1);}
count=countrows(bimfile);

if((input=fopen(bimfile,"r"))==NULL)
{printf("Error opening %s\n\n",bimfile);exit(1);}

found=0;bcount=0;mcount=0;ncount=0;scount=0;lcount=0;
for(j=0;j<count;j++)
{
want=j;
if(keeppreds_use!=NULL){want=keeppreds_use[found];}
if(fscanf(input, "%s %s %s %s %s %s ", rc, rs, rm, rbp, ra1, ra2)!=6)
{printf("Error reading Row %d of %s\n\n", j+1, bimfile);exit(1);}

if(j==want)	//will be using row
{
//save chr
chr[found]=-1;
if(strcmp(rc,"X")==0){chr[found]=23;}
if(strcmp(rc,"Y")==0){chr[found]=24;}
if(strcmp(rc,"XY")==0){chr[found]=25;}
if(strcmp(rc,"MT")==0){chr[found]=26;}
if(strcmp(rc,"0")==0){chr[found]=0;}

if(chr[found]==-1)	//not found yet
{
chr[found]=atoi(rc);

if(chr[found]<=0)	//so not valid
{printf("Error, Predictor %s has chromosome %s; all values must be a positive integer, X (23), Y (24), XY (25), MT (26) or 0\n\n", rs, rc);exit(1);}
if(chr[found]>26)	//unexpectedly large
{
if(type==1)
{
if(bcount<5){printf("Warning, Predictor %s has chromosome %s (the largest expected for humans is 26)\n", rs, rc);}
bcount++;
}}
}

if(found>0&&flag==0)	//check not smaller than previous one
{
if(chr[found]<chr[found-1])
{printf("Error, chromosome for %s (%d) is lower than that for %s (%d)\n\n", rs, chr[found], preds[found-1], chr[found-1]);exit(1);}
}

//save name
copy_string(preds,found,rs);

//save cm and check whether negative
cm[found]=atof(rm);
if(cm[found]<0)
{
if(type==1)
{
if(mcount<5){printf("Warning, Predictor %s has a negative genetic distance (%s)\n", rs, rm);}
mcount++;
}}

if(window_cm!=-9999&&found>0&&flag==0)	//check cm consistent with previous one
{
if(chr[found]==chr[found-1]&&cm[found]<cm[found-1])
{printf("Error, genetic distance for %s (%.2f) is lower than that for %s (%.2f)\n", rs, cm[found],  preds[found-1], cm[found-1]);exit(1);}
}

//save bp and check
bp[found]=atof(rbp);
if(bp[found]<=0)
{
if(type==1)
{
if(ncount<5){printf("Warning, Predictor %s has a non-positive basepair (%s) and will be ignored\n", rs, rbp);}
ncount++;
}}

if(found>0&&bp[found]>0&&flag==0)	//check bp consistent with previous one and for duplicates
{
if(chr[found]==chr[found-1]&&bp[found]<bp[found-1])
{printf("Error, basepair for %s (%.2f) is lower than that for %s (%.2f)\n", rs, bp[found],  preds[found-1], bp[found-1]);exit(1);}
if(chr[found]==chr[found-1]&&bp[found]==bp[found-1])
{
if(type==1)
{
if(scount<5){printf("Warning, %s and %s have the same basepair (%.2f)\n", preds[found-1], preds[found], bp[found]);}
scount++;
}}
}

//check then load alleles
if(strlen(ra1)+strlen(ra2)>2)
{
if(type==1)
{
if(lcount<5){printf("Warning, Predictor %s has multi-character alleles (%s and %s) and will be ignored\n", preds[found], ra1, ra2);}
lcount++;
}
bp[found]=-1;
}
al1[found]=ra1[0];al2[found]=ra2[0];

found++;
if(found==length){break;}
}	//end of found
}	//end of j loop
fclose(input);

if(window_cm!=-9999)	//will be using maps, so check there are some valid values
{
count=0;for(j=0;j<found;j++){count+=(cm[j]!=0);}
if(count==0)
{printf("\nError, Column 3 of %s should contain genetic distances\n\n", bimfile);exit(1);}
}

if(bcount>5){printf("In total %d predictors have chromosome values greater than 26\n", bcount);}
if(mcount>5){printf("In total %d predictors have negative genetic distances\n", mcount);}
if(ncount>5){printf("In total %d predictors have non-positive basepairs\n", ncount);}
if(scount>5){printf("In total %d pairs of predictors have the same basepair\n", scount);}
if(lcount>5){printf("In total %d predictors have multi-character alleles\n", lcount);}
//if(bcount+mcount+ncount+scount+lcount>0){printf("\n");}

free(rc);free(rs);free(rm);free(rbp);free(ra1);free(ra2);
return(0);
}	//end of read_bimfile

///////////////////////////

int read_genfile(char *genfile, char **preds, double *bp, char *al1, char *al2, int num_samples, int genskip, int genheaders, int genprobs, int maxpreds, int type)
//type=0 - keep quiet, type=1 warn about everything
{
int found, ncount, scount, lcount, size;

char *rs, *rs2, *rbp, *ra1, *ra2;
char *gzbuffer;

gzFile inputgz;

rs=malloc(sizeof(char)*10000000);rs2=malloc(sizeof(char)*10000000);
rbp=malloc(sizeof(char)*10000000);ra1=malloc(sizeof(char)*10000000);ra2=malloc(sizeof(char)*10000000);


//open, which checks size
open_datagz(&inputgz, genfile, num_samples, genskip, genheaders, genprobs);

//now can read
size=30000000+num_samples*(genprobs*20+(genprobs==0)*4);
gzbuffer=malloc(sizeof(char)*size);

found=0;ncount=0;scount=0;lcount=0;
while(gzgets(inputgz,gzbuffer,size)!=NULL)
{
if(strlen(gzbuffer)==size-1)
{printf("Error reading %s; Row %d is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", genfile, genskip+1, (int)strlen(gzbuffer), size);exit(1);}
if(genskip+found+1%20000==0){printf("Reading Row %d\n", genskip+found+1);}

if(genheaders==3)
{
if(sscanf(gzbuffer, "%s %s %s", rs, ra1, ra2)!=3)
{printf("Error reading details on Row %d\n\n", genskip+found+1);exit(1);}
}
if(genheaders==4)
{
if(sscanf(gzbuffer, "%s %s %s %s ", rs, rbp, ra1, ra2)!=4)
{printf("Error reading details on Row %d\n\n", genskip+found+1);exit(1);}
}
if(genheaders==5)
{
if(sscanf(gzbuffer, "%s %s %s %s %s ", rs2, rs, rbp, ra1, ra2)!=5)
{printf("Error reading details on Row %d\n\n", genskip+found+1);exit(1);}
}

if(found<maxpreds)	//saving
{
//save name
copy_string(preds,found,rs);

if(genheaders==4||genheaders==5)	//save bp and check
{
bp[found]=atof(rbp);
if(bp[found]<=0)
{
if(type==1)
{
if(ncount<5){printf("Warning, Predictor %s has a non-positive basepair (%s) and will ignored\n", rs, rbp);}
ncount++;
}}

if(found>0&&bp[found]>0)	//check bp consistent with previous one and for duplicates
{
if(bp[found]<bp[found-1])
{printf("Error, basepair for %s (%.2f) is lower than that for %s (%.2f)\n", rs, bp[found],  preds[found-1], bp[found-1]);exit(1);}
if(bp[found]==bp[found-1])
{
if(type==1)
{
if(scount<5){printf("Warning, %s and %s have the same basepair (%.2f)\n", preds[found-1], preds[found], bp[found]);}
scount++;
}}
}
}	//end of have bp
else
{bp[found]=found;}

//check then load alleles
if(strlen(ra1)+strlen(ra2)>2)
{
if(type==1)
{
if(lcount<5){printf("Warning, Predictor %s has multi-character alleles (%s and %s) and will be ignored\n", preds[found], ra1, ra2);}
lcount++;
}
bp[found]=-1;
}

al1[found]=ra1[0];al2[found]=ra2[0];
}
found++;
}	//end of while loop
gzclose(inputgz);

if(found<=maxpreds)
{
if(ncount>5){printf("In total %d predictors have non-positive basepairs\n", ncount);}
if(scount>5){printf("In total %d pairs of predictors have the same basepair\n", scount);}
if(lcount>5){printf("In total %d predictors have multi-character alleles\n", lcount);}
}

free(rs);free(rs2);free(rbp);free(ra1);free(ra2);free(gzbuffer);
return(found);
}	//end of read_genfile

///////////////////////////

int read_weightsfile(char *weightsfile, double *weights, int length, char **preds, int *keeppreds_use, int num_preds, char **allpreds, char *bimfile)
{
int j, count, count2, count3, head;
int *indexer, *indexer2;

char **gotpreds;
double *weightsb;


count=countcols(weightsfile);

if(count==2)	//read all in, then match up
{
head=check_head(weightsfile,"Predictor","SNP",0);
count2=countrows(weightsfile)-head;

gotpreds=malloc(sizeof(char*)*count2);
weightsb=malloc(sizeof(double)*count2);
indexer=malloc(sizeof(int)*length);
indexer2=malloc(sizeof(int)*length);

read_strings(weightsfile, gotpreds, count2, NULL, 1, head);
read_values(weightsfile, weightsb, count2, NULL, 2, head, 0);
count3=find_strings(preds, length, gotpreds, count2, indexer, indexer2, NULL, weightsfile, NULL, NULL, 3);
if(count3==0)
{printf("Error, %s does not contain weights for any of the %d predictors\n\n", weightsfile, length);exit(1);}
if(count3<length)
{printf("Warning, %s contains weights for only %d of the %d predictors; the remaining predictors will be given weight zero\n\n", weightsfile, count3, length);}
for(j=0;j<length;j++){weights[j]=0;}
for(j=0;j<count3;j++){weights[indexer[j]]=weightsb[indexer2[j]];}
for(j=0;j<count2;j++){free(gotpreds[j]);};free(gotpreds);
free(weightsb);free(indexer);free(indexer2);
}

if(count==6)	//weightsfile should align with bimfile
{
count2=countrows(weightsfile)-1;
if(count2!=num_preds)
{printf("Error, the number of predictors in %s (%d) does not match the number of rows of %s (%d), suggesting these weights correspond to a different dataset\nTo use these weights regardless, create a two-column file containing just the predictor names and weights (predictors with zero weights can be excluded), for example, using the command:\nawk < %s \'($2>0){print $1, $2}\' > weights.short\n\n", weightsfile, count2, bimfile, num_preds, weightsfile);exit(1);}

gotpreds=malloc(sizeof(char*)*count2);
read_strings(weightsfile, gotpreds, count2, NULL, 1, 1);
for(j=0;j<count2;j++)
{
if(strcmp(allpreds[j], gotpreds[j])!=0)
{printf("Error, %s does not align with %s, they have different predictors (%s and %s) on Row %d, indicating that these weights correspond to a different dataset\nTo use these weights regardless, create a two-column file containing just the predictor names and weights (predictors with zero weights can be excluded), for example, using the command:\nawk < %s \'($2>0){print $1, $2}\' > weights.short\n\n", weightsfile, bimfile, gotpreds[j], allpreds[j], j+1, weightsfile);exit(1);}
}
for(j=0;j<count2;j++){free(gotpreds[j]);}free(gotpreds);
read_values(weightsfile, weights, length, keeppreds_use, 2, 1, 0);
}

//check for -1 and negative
for(j=0;j<length;j++)
{
if(weights[j]==-1&&count==6)
{printf("Error reading %s; Predictor %s has weight -1, indicating this predictor was not included when calculating weights, so can not be used in subsequent applications\nTo use these weights regardless, create a two-column file containing just the predictor names and weights (predictors with zero weights can be excluded), for example, using the command:\nawk < %s \'($2>0){print $1, $2}\' > weights.short\n\n", weightsfile, preds[j], weightsfile);exit(1);}
if(weights[j]<0)
{printf("Error reading %s; Predictor %s has weight %.6f (all weights should be non-negative)\n\n", weightsfile, preds[j], weights[j]);exit(1);}
}

return(0);
}	//end of read_weightsfile

////////

void read_infosfile(char *infosfile, double *infos, int length, char **preds)
{
int j, count, count2, head;
int *indexer, *indexer2;

char **gotpreds;
double *infosb;


head=check_head(infosfile,"Predictor","SNP",0);
count=countrows(infosfile)-head;

gotpreds=malloc(sizeof(char*)*count);
infosb=malloc(sizeof(double)*count);
indexer=malloc(sizeof(int)*length);
indexer2=malloc(sizeof(int)*length);

read_strings(infosfile, gotpreds, count, NULL, 1, head);
read_values(infosfile, infosb, count, NULL, 2, head, 0);
count2=find_strings(preds, length, gotpreds, count, indexer, indexer2, NULL, infosfile, NULL, NULL, 3);

if(count2==0)
{printf("Error, %s does not contain infos for any of the %d predictors\n\n", infosfile, length);exit(1);}
if(count2<length)
{printf("Warning, %s contains infos for only %d of the %d predictors; the remaining %d predictors will be given info zero\n\n", infosfile, count2, length, length-count2);}
for(j=0;j<length;j++){infos[j]=0;}
for(j=0;j<count2;j++){infos[indexer[j]]=infosb[indexer2[j]];}

//check for invalid
for(j=0;j<length;j++)
{
if(infos[j]<0||infos[j]>1.01)
{printf("Error reading %s; Predictor %s has info %.6f (infos should be within [0,1])\n\n", infosfile, preds[j], infos[j]);exit(1);}
}

for(j=0;j<count;j++){free(gotpreds[j]);};free(gotpreds);
free(infosb);free(indexer);free(indexer2);
}	//end of read_infosfile

///////////////////////////

void read_sumsfile(char *sumsfile, double *nss, double *chis, double *rhos, int length, char **preds, char *al1, char *al2, char *filename, int amb, int fixn, double scaling, int type, int checksums)
{
//type=0 - normal, type=1 - sum-hers or sum-cors with plet=1, type=2 - sum-cors with plet=0
//for types 0 and 2, need signed statistics, for 1 unsigned
//for type 0, must have all predictions; for types 1 and 2, not required when checksums=0
int j, j2, j3, k, count, count2, count3, flag, ccount, acount, pcount;

int cols[6], zsp, *indexer, *indexer2;
int *signsb;
double *nssb;
double *chisb;
char **gotpreds, *gotal1, *gotal2;

char readchar, **readlist, *rs;

FILE *input;

rs=malloc(sizeof(char)*10000000);


//read all summaries, doing some basic checks - have already checked have required columns

count2=countcols(sumsfile);
if(find_head("Predictor", sumsfile, count2)!=-1){cols[0]=find_head("Predictor", sumsfile, count2);}
else{cols[0]=find_head("SNP", sumsfile, count2);}
cols[1]=find_head("A1", sumsfile, count2);
cols[2]=find_head("A2", sumsfile, count2);

if(type==0||type==2)	//need signed test statistics
{
if(find_head("Z", sumsfile, count2)!=-1){cols[3]=find_head("Z", sumsfile, count2);cols[5]=cols[3];zsp=1;}
else
{
if(find_head("Stat", sumsfile, count2)!=-1){cols[3]=find_head("Stat", sumsfile, count2);cols[5]=find_head("Direction", sumsfile, count2);zsp=2;}
else{cols[3]=find_head("P", sumsfile, count2);cols[5]=find_head("Direction", sumsfile, count2);zsp=3;}
}
}
else	//need unsigned test statistics
{
if(find_head("Z", sumsfile, count2)!=-1){cols[3]=find_head("Z", sumsfile, count2);zsp=1;}
else
{
if(find_head("Stat", sumsfile, count2)!=-1){cols[3]=find_head("Stat", sumsfile, count2);zsp=2;}
else{cols[3]=find_head("P", sumsfile, count2);zsp=3;}
}
}

if(fixn==-9999)
{
if(find_head("n", sumsfile, count2)!=-1){cols[4]=find_head("n", sumsfile, count2);}
else{cols[4]=find_head("N", sumsfile, count2);}
}

////////

count=countrows(sumsfile)-1;
gotpreds=malloc(sizeof(char*)*count);
gotal1=malloc(sizeof(char)*count);
gotal2=malloc(sizeof(char)*count);
nssb=malloc(sizeof(double)*count);
signsb=malloc(sizeof(int)*count);
chisb=malloc(sizeof(double)*count);

readlist=malloc(sizeof(char*)*count2);
for(j=0;j<count2;j++){readlist[j]=malloc(sizeof(char)*500);}

if((input=fopen(sumsfile,"r"))==NULL)
{printf("Error opening %s\n\n", sumsfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

for(j=0;j<count;j++)
{
for(k=0;k<count2;k++)
{
if(fscanf(input, "%s ", rs)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", k+1, j+2, sumsfile);exit(1);}
if(strlen(rs)>500){free(readlist[k]);copy_string(readlist,k,rs);}
else{strcpy(readlist[k],rs);}
}

if(strlen(readlist[cols[1]])+strlen(readlist[cols[2]])>2)
{printf("Error, Predictor %s has alleles %s and %s (both alleles must be single characters)\n\n", readlist[cols[0]], readlist[cols[1]], readlist[cols[2]]);exit(1);}

if(strcmp(readlist[cols[1]],readlist[cols[2]])==0)
{printf("Error, the two alleles for Predictor %s are the same (%s)\n\n", readlist[cols[0]], readlist[cols[1]]);exit(1);}

//copy in names and alleles
copy_string(gotpreds,j,readlist[cols[0]]);
gotal1[j]=readlist[cols[1]][0];
gotal2[j]=readlist[cols[2]][0];

if(zsp==1)	//using Z
{
if(sscanf(readlist[cols[3]], "%[0-9eE.+-] ", rs)==1){chisb[j]=pow(atof(rs),2);}
else
{
if(strcmp(readlist[cols[3]],"NA")==0)
{printf("Error, the Z statistic for Predictor %s is NA\n\n", gotpreds[j]);exit(1);}
else{printf("Error, the Z statistic for Predictor %s is unrecognisable (%s)\n\n", gotpreds[j], readlist[cols[3]]);exit(1);}
}
}
if(zsp==2)	//using Stat
{
if(sscanf(readlist[cols[3]], "%[0-9eE.+-] ", rs)==1)
{
if(atof(rs)<0)
{printf("Error, the chi-squared statistic for Predictor %s should be non-negative (not %s); make sure %s provides chi-squared test statistics (not Z statistics)\n\n", gotpreds[j], rs, sumsfile);exit(1);}
chisb[j]=atof(rs);
}
else
{
if(strcmp(readlist[cols[3]],"NA")==0)
{printf("Error, the Chi-Squared statistic for Predictor %s is NA\n\n", gotpreds[j]);exit(1);}
else{printf("Error, the Chi-Squared statistic for Predictor %s is unrecognisable (%s)\n\n", gotpreds[j], readlist[cols[3]]);exit(1);}
}
}
if(zsp==3)	//using P
{
if(sscanf(readlist[cols[3]], "%[0-9eE.+-] ", rs)==1)
{
if(atof(rs)<0&&atof(rs)>1)
{printf("Error, the p-value for Predictor %s should be within [0,1] (not %s)\n\n", gotpreds[j], rs);exit(1);}
if(atof(rs)<1e-300){chisb[j]=-1;}	//will warn below if predictor used
else{chisb[j]=pow(normal_inv(atof(rs)/2),2);}
}
else
{
if(strcmp(readlist[cols[3]],"NA")==0)
{printf("Error, the p-value for Predictor %s is NA\n\n", gotpreds[j]);exit(1);}
else{printf("Error, the p-value for Predictor %s is unrecognisable (%s)\n\n", gotpreds[j], readlist[cols[3]]);exit(1);}
}
}

if(fixn==-9999)
{
if(sscanf(readlist[cols[4]], "%[0-9eE.+-] ", rs)==1)
{
if(atof(rs)<0.5){printf("Error, Predictor %s has non-positive sample size (%s)\n\n", gotpreds[j], rs);exit(1);}
nssb[j]=atof(rs);
}
else
{
if(strcmp(readlist[cols[4]],"NA")==0)
{printf("Error, the sample size for Predictor %s is NA\n\n", gotpreds[j]);exit(1);}
else
{printf("Error, the sample size for Predictor %s is unrecognisable (%s)\n\n", gotpreds[j], readlist[cols[3]]);exit(1);}
}
}
else{nssb[j]=fixn;}

if(type==0||type==2)	//need signs
{
if(sscanf(readlist[cols[5]], "%[0-9eE.+-] ", rs)==1)
{
if(atof(rs)>=0){signsb[j]=1;}
else{signsb[j]=-1;}
}
else
{
if(strcmp(readlist[cols[5]],"NA")==0)
{printf("Error, the direction for Predictor %s is NA\n\n", gotpreds[j]);exit(1);}
else{printf("Error, the direction for Predictor %s is unrecognisable (%s)\n\n", gotpreds[j], readlist[cols[5]]);exit(1);}
}
}
else{signsb[j]=1;}

}
fclose(input);

for(j=0;j<count2;j++){free(readlist[j]);}free(readlist);

////////

//line up
indexer=malloc(sizeof(int)*length);
indexer2=malloc(sizeof(int)*length);
count2=find_strings(preds, length, gotpreds, count, indexer, indexer2, NULL, sumsfile, NULL, NULL, 3);

//load up, keeping track of consistent, ambiguous, and (if used) zero pvalues
for(j=0;j<length;j++){nss[j]=-9999;chis[j]=-9999;rhos[j]=-9999;}

count3=0;
ccount=0;acount=0;pcount=0;
for(j=0;j<count2;j++)
{
j2=indexer[j];
j3=indexer2[j];
flag=0;

if((al1[j2]!=gotal1[j3]&&al1[j2]!=gotal2[j3])||(al2[j2]!=gotal1[j3]&&al2[j2]!=gotal2[j3]))	//inconsistent
{
flag=1;
if(ccount<5){printf("Warning, Predictor %s will be ignored as its alleles (%c and %c) are not consistent with those in %s (%c and %c)\n", gotpreds[j3], gotal1[j3], gotal2[j3], filename, al1[j2], al2[j2]);}
ccount++;
}
else	//alleles consistent
{
if((int)gotal1[j3]+(int)gotal2[j3]==138||(int)gotal1[j3]+(int)gotal2[j3]==149)	//ambiguous
{
if(amb==0)	//will ignore
{
flag=1;
if(acount<5){printf("Warning, Predictor %s will be ignored as it has ambiguous alleles (%c and %c)\n", gotpreds[j3], al1[j2], al2[j2]);}
}
acount++;
}
}

if(chisb[j3]==-1&&flag==0)
{
if(pcount<5){printf("Warning, Predictor %s has p-value below 1e-300; will be replaced with 1e-300\n", gotpreds[j3]);}
chisb[j3]=pow(normal_inv(1e-300/2),2);
pcount++;
}

if(flag==0)	//can use
{
nss[j2]=nssb[j3];
chis[j2]=chisb[j3]*scaling;
rhos[j2]=signsb[j3]*pow(chisb[j3]/(chisb[j3]+nssb[j3]),.5);	//should really adjust demon for num covar in gwas
if(al1[j2]!=gotal1[j3]){rhos[j2]*=-1;}
count3++;
}
}

if(ccount>5||(acount>5&&amb==0)||pcount>5){printf("\n");}
if(ccount>5){printf("In total, %d predictors have inconsistent alleles\n", ccount);}
if(amb==0)
{
if(acount>5){printf("In total, %d predictors have ambiguous alleles\n", acount);}
if(acount>0){printf("If you are confident that ambiguous predictors are correctly oriented, you can retain them by adding \"--allow-ambiguous YES\"\n");}
}
if(amb==1&&acount>0){printf("In total, %d predictors have ambiguous alleles\n", acount);}
if(pcount>5){printf("In total %d predictors have p-value zero\n", pcount);}

if(count3==0){printf("\nError, %s contains summary statistics for none of the %d predictors\n\n", sumsfile, length);exit(1);}
if(count3<length)
{
if(type==0){printf("\nError, %s contains summary statistics for only %d of the %d predictors; either correct the file or use \"--extract\" or \"--exclude\" to exclude predictors without summary statistics\n\n", sumsfile, count3, length);exit(1);}

if(checksums==1)
{printf("\nError, %s contains summary statistics for only %d of the %d predictors; to continue despite this error, use \"--check-sums NO\"\n\n", sumsfile, count3, length);exit(1);}
printf("Warning, %s contains summary statistics for only %d of the %d predictors\n", sumsfile, count3, length);
}
else{printf("\nHave found summary statistics for all %d predictors\n", length);}

for(j=0;j<count;j++)
{free(gotpreds[j]);};
free(gotpreds);
free(gotal1);free(gotal2);free(nssb);free(signsb);free(chisb);
free(indexer);free(indexer2);

free(rs);
}	//end of read_sumsfile

///////////////////////////

int read_tagfile(char *tagfile, char **catlabels, char **spreds, char *sal1, char *sal2, double *stags, double **svars, double **ssums, int num_parts, int length, int parttype, char *altfile, char *bpredfile, char *cpredfile, int tagone)
{
int j, q, q2, count, count2, head;
int *predorder, *usedpreds, *indexer, *indexer2;
double *avalues;

char readstring[500], **wantpreds, *rs;

FILE *input;

rs=malloc(sizeof(char)*10000000);


printf("Reading details for %d predictors from %s\n\n", length, tagfile);

//first read the whole file

if((input=fopen(tagfile,"r"))==NULL)
{printf("Error opening %s\n\n",tagfile);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", rs, readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=9)
{printf("Error reading Row 1 of %s\n\n", tagfile);exit(1);}
if(strcmp(rs,"Predictor")!=0)
{printf("Error, %s should begin \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, readstring);exit(1);}

//save labels
for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%s ", catlabels[q])!=1)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile);exit(1);}
}

for(j=0;j<length;j++)
{
if(fscanf(input, "%s %c %c %s %lf %s %s %s %s ", rs, sal1+j, sal2+j, readstring, stags+j, readstring, readstring, readstring, readstring)!=9)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", j+2, tagfile);exit(1);}
copy_string(spreds, j, rs);

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%lf ", svars[q]+j)!=1)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", j+2, tagfile);exit(1);}
}
}

//now penultimate num_part rows
for(q2=0;q2<num_parts;q2++)
{
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", length+2+q2, tagfile);exit(1);}

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%lf ", ssums[q]+q2)!=1)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", length+2+q2, tagfile);exit(1);}
}
}

//and the last row
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", length+2+num_parts, tagfile);exit(1);}
if(strcmp(readstring,"The")==0)
{printf("Error, the tagging file %s was made using an older version of LDAK; please remake using this version\n\n", tagfile);exit(1);}
if(strcmp(readstring,"There")!=0)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", length+2+num_parts, tagfile);exit(1);}

if(fscanf(input, "%s %d %s %s %s %s %s %s ", readstring, &count, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", length+2+num_parts, tagfile);exit(1);}

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%lf ", ssums[q]+num_parts)!=1)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", length+2+num_parts, tagfile);exit(1);}
}

fclose(input);

//for printing need proportion of predictors in each category
for(q=0;q<num_parts;q++){ssums[q][num_parts+1]=ssums[q][num_parts]/count;}

//old code was exact for annotations, and approximate for partitions
/*
if(parttype==0)	//using annotations - ssums[num_parts-1][num_parts] is total number of predictors
{
for(q=0;q<num_parts;q++){ssums[q][num_parts+1]=ssums[q][num_parts]/ssums[num_parts-1][num_parts];}
}
else	//using partitions - sum of ssums[][num_parts] is (approx) total number of predictors
{
sum=0;for(q=0;q<num_parts;q++){sum+=ssums[q][num_parts];}
for(q=0;q<num_parts;q++){ssums[q][num_parts+1]=ssums[q][num_parts]/sum;}
}
*/

////////

if(strcmp(altfile,"blank")!=0)	//reduce to regresion predictors and update tags
{
head=check_head(altfile,"Predictor","SNP",0);
count=countrows(altfile)-head;
printf("Reading %d regression predictors and taggings from %s\n", count, altfile);
wantpreds=malloc(sizeof(char*)*count);
avalues=malloc(sizeof(double)*count);
read_strings(altfile, wantpreds, count, NULL, 1, head);
read_values(altfile, avalues, count, NULL, 2, head, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(spreds, length, wantpreds, count, indexer, indexer2, tagfile, altfile, NULL, NULL, 3);
if(count2==0){printf("Error, none of these are in %s\n\n", tagfile);exit(1);}
if(count2<count){printf("Warning, only %d of these are in %s\n", count2, tagfile);}
for(j=0;j<count2;j++){stags[indexer[j]]=avalues[indexer2[j]];}

//squeeze down
for(j=0;j<count2;j++)
{
if(indexer[j]!=j)
{
free(spreds[j]);copy_string(spreds,j,spreds[indexer[j]]);
sal1[j]=sal1[indexer[j]];
sal2[j]=sal2[indexer[j]];
stags[j]=stags[indexer[j]];
for(q=0;q<num_parts;q++){svars[q][j]=svars[q][indexer[j]];}
}}
for(j=count2;j<length;j++){free(spreds[j]);}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);free(indexer2);free(avalues);
length=count2;
}

if(tagone==1)	//set tags to one - hopefully using pruned predictors
{
for(j=0;j<length;j++){stags[j]=1;}
}

if(strcmp(bpredfile,"blank")!=0||strcmp(cpredfile,"blank")!=0)	//extract
{
predorder=malloc(sizeof(int)*length);
check_dups(spreds,length,tagfile,predorder,1);

usedpreds=malloc(sizeof(int)*length);
for(j=0;j<length;j++){usedpreds[j]=1;}

count=extraction(usedpreds, length, spreds, NULL,predorder, bpredfile, cpredfile, -9999, "blank", tagfile);
if(count<length){printf("Will use %d of the predictors in %s\n\n", count, tagfile);}
else{printf("Will use all %d predictors in %s\n\n", length, tagfile);}

//squeeze down
count=0;
for(j=0;j<length;j++)
{
if(usedpreds[j]==1)
{
if(count!=j)
{
free(spreds[count]);copy_string(spreds,count,spreds[j]);
sal1[count]=sal1[j];
sal2[count]=sal2[j];
stags[count]=stags[j];
for(q=0;q<num_parts;q++){svars[q][count]=svars[q][j];}
}
count++;
}}
for(j=count;j<length;j++){free(spreds[j]);}
length=count;

free(predorder);free(usedpreds);
}

free(rs);
return(length);
}	//read_tagfile

///////////////////////////

int read_genefile(char *genefile, char **gnames, int *gchr, double *gbp1, double *gbp2, int *gstrand)
{
int j, count, count2, bcount, zcount, head;

char readchar, *rs, *rc, *rbp1, *rbp2;

FILE *input;

rs=malloc(sizeof(char)*10000000);rc=malloc(sizeof(char)*10000000);
rbp1=malloc(sizeof(char)*10000000);rbp2=malloc(sizeof(char)*10000000);

head=check_head(genefile,"Gene","Name",1);
count=countrows(genefile)-head;
count2=countcols(genefile);
printf("Reading details for %d genes from %s; LDAK assumes this file is in BED (Browser Extensible Data) format, which uses 0-start, half-open coordinates (for example, if Gene ABC is on Chromosome 7 and contains basepairs 1-10, the corresponding row would be \"ABC 7 0 10\", \"ABC 7 0 10 +\" or \"ABC 7 0 10 -\")\n\n", count, genefile);

//open and skip the header row if present
if((input=fopen(genefile,"r"))==NULL)
{printf("Error opening %s\n\n",genefile);exit(1);}
if(head==1)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

bcount=0;zcount=0;
for(j=0;j<count;j++)
{
gstrand[j]=1;
if(count2==4)
{
if(fscanf(input, "%s %s %s %s ", rs, rc, rbp1, rbp2)!=4)
{printf("Error reading Row %d of %s\n\n", j+1, genefile);exit(1);}
}
else
{
if(fscanf(input, "%s %s %s %s %c ", rs, rc, rbp1, rbp2, &readchar)!=5)
{printf("Error reading Row %d of %s\n\n", j+1, genefile);exit(1);}
if(readchar!='+'&&readchar!='-')
{printf("Error reading Row %d of %s; the final column should contain \"+\" or \"-\" (not %c)\n\n", j+1, genefile, readchar);exit(1);}
if(readchar=='-'){gstrand[j]=-1;}
}

copy_string(gnames, j, rs);

gchr[j]=-1;
if(strcmp(rc,"X")==0){gchr[j]=23;}
if(strcmp(rc,"Y")==0){gchr[j]=24;}
if(strcmp(rc,"XY")==0){gchr[j]=25;}
if(strcmp(rc,"MT")==0){gchr[j]=26;}
if(strcmp(rc,"0")==0){gchr[j]=0;}

if(gchr[j]==-1)	//not found so far
{
gchr[j]=atoi(rc);

if(gchr[j]<=0)	//so not valid
{printf("Error reading %s; Gene %s has chromosome %s\nValue must be positive integer, X (23), Y (24), XY (25), MT (26) or 0\n\n", genefile, rs, rc);exit(1);}
if(gchr[j]>26)	//so unexpectedly large
{
if(bcount<5)
{printf("Warning, Gene %s has chromosome %s (the largest expected for humans is 26)\n", rs, rc);}
bcount++;
}
}

if(j>0)	//check not smaller than previous one
{
if(gchr[j]<gchr[j-1])
{printf("Error, chromosome for Gene %s (%d) is lower than that for Gene %s (%d); genes should be in genomic order\n\n", gnames[j], gchr[j], gnames[j-1], gchr[j-1]);exit(1);}
}

gbp1[j]=atof(rbp1)+1;	//this allows for fact that bed format starts at zero
if(gbp1[j]<=0)
{printf("Error reading %s; Gene %s has a non-positive start basepair (%s)\n\n", genefile, rs, rbp1);exit(1);}
gbp2[j]=atof(rbp2);
if(gbp2[j]<=0)
{printf("Error reading %s; Gene %s has a non-positive end basepair (%s)\n\n", genefile, rs, rbp2);exit(1);}

if(strcmp(rbp1,rbp2)==0)
{printf("Error reading %s; Gene %s has equal basepairs (%s)\nNote that the file should be 0-based; for example, if a gene contains only Basepair 1000, its entry should have start 999 and end 1000\n\n", genefile, rs, rbp1);exit(1);}
if(gbp1[j]-1>gbp2[j])
{printf("Error reading %s; Gene %s has a higher start basepair (%s) than end basepair (%s)\n\n", genefile, rs, rbp1, rbp2);exit(1);}

if(gbp1[j]==gbp2[j])
{
if(zcount<5)
{printf("Warning, Gene %s has length one\n", rs);}
zcount++;
}

if(j>0)	//check start bp consistent with previous one (end bps do not need to be consistent)
{
if(gchr[j]==gchr[j-1]&&gbp1[j]<gbp1[j-1])
{printf("Error, Gene %s has a lower start basepair (%.2f) than Gene %s (%.2f); genes should be in genomic order\n\n", gnames[j], gbp1[j], gnames[j-1], gbp1[j-1]);exit(1);}
}
}	//end of j loop
fclose(input);

if(bcount>5){printf("In total %d genes have chromosome values greater than 26\n", bcount);}
if(zcount>5){printf("In total %d genes have length one\n", zcount);}
if(bcount+zcount>0){printf("\n");}

//check whether gene names are unique

free(rs);free(rc);free(rbp1);free(rbp2);
return(0);
}	//end of read_genefile

///////////////////////////

int read_centresfile(char *centresfile, double *centres, int length, char **preds, char *al1, char *al2, char *bimfile)
{
int j, j2, count, count2, found, head;
int *indexer, *indexer2;

char **gotpreds;

char readchar, *rs, *ra1, *ra2;

FILE *input;

rs=malloc(sizeof(char)*10000000);
ra1=malloc(sizeof(char)*10000000);ra2=malloc(sizeof(char)*10000000);


//start by getting indexes
head=check_head(centresfile,"Predictor","SNP",0);
count=countrows(centresfile)-head;

gotpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*length);
indexer2=malloc(sizeof(int)*length);

read_strings(centresfile, gotpreds, count, NULL, 1, head);
count2=find_strings(gotpreds, count, preds, length, indexer, indexer2, centresfile, NULL, NULL, NULL, 3);
if(count2==0)
{printf("Error, %s does not contain centres for any of the %d predictors\n\n", centresfile, length);exit(1);}
if(count2<length)
{printf("Error, %s contains centres for only %d of the %d predictors\n", centresfile, count2, length);exit(1);}

//now read - open centresfile and skip the header row if present
if((input=fopen(centresfile,"r"))==NULL)
{printf("Error opening %s\n\n", centresfile);exit(1);}
if(head==1)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

found=0;
for(j=0;j<count;j++)
{
if(j==indexer[found])	//will be using
{
j2=indexer2[found];

//sort name and alleles
if(fscanf(input, "%s %s %s ", rs, ra1, ra2)!=3)
{printf("Error reading Row %d of %s\n\n", j+1, centresfile);exit(1);}

if(strlen(ra1)+strlen(ra2)>2)
{printf("Error, Predictor %s has multi-character alleles (%s and %s)\n\n", rs, ra1, ra2);exit(1);}
if(ra1[0]==ra2[0])
{printf("Error, both alleles of Predictor %s are %c\n\n", rs, ra1[0]);exit(1);} 

if((ra1[0]!=al1[j2]&&ra1[0]!=al2[j2])||(ra2[0]!=al1[j2]&&ra2[0]!=al2[j2]))
{printf("Error reading %s; alleles for predictor %s (%c and %c) are not consistent with those in %s (%c and %c)\n\n", centresfile, rs, ra1[0], ra2[0], bimfile, al1[j2], al2[j2]);exit(1);}

//read centre then make sure at end of row
if(fscanf(input, "%[0-9eE.+-]%c", rs, &readchar)!=2)
{
if(fscanf(input, "%s%c", rs, &readchar)!=2)
{printf("Error reading centre of Row %d of %s\n\n", head+j+1, centresfile);exit(1);}
if(strcmp(rs,"NA")==0)
{printf("Error reading %s; centre of Row %d is NA\n\n", centresfile,  head+j+1);exit(1);}
else
{printf("Error reading %s; centre of Row %d is unrecognisable (%s)\n\n", centresfile,  head+j+1, rs);exit(1);}
}
while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//load up
if(al1[j2]==ra1[0]){centres[j2]=atof(rs);}
else{centres[j2]=2-atof(rs);}

found++;
if(found==length){break;}
}
else	//skip row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}	//end of j loop

return(0);
}	//end of read_centresfile

////////

int read_pvafile(char *pvafile, double *pvalues, int length, char **preds)
{
int j, count, count2, head;
int *indexer, *indexer2;

char **gotpreds;
double *pvaluesb;


head=check_head(pvafile,"Predictor","SNP",0);
count=countrows(pvafile)-head;

gotpreds=malloc(sizeof(char*)*count);
pvaluesb=malloc(sizeof(double)*count);
indexer=malloc(sizeof(int)*length);
indexer2=malloc(sizeof(int)*length);

read_strings(pvafile, gotpreds, count, NULL, 1, head);
read_values(pvafile, pvaluesb, count, NULL, 2, head, 0);
count2=find_strings(preds, length, gotpreds, count, indexer, indexer2, NULL, pvafile, NULL, NULL, 3);
if(count2==0)
{printf("Error, %s does not contain p-values for any of the %d predictors\n\n", pvafile, length);exit(1);}
if(count2<length)
{printf("Warning, %s contains p-values for only %d of the %d predictors\n\n", pvafile, count2, length);}
for(j=0;j<length;j++){pvalues[j]=2;}
for(j=0;j<count2;j++){pvalues[indexer[j]]=pvaluesb[indexer2[j]];}

for(j=0;j<count;j++){free(gotpreds[j]);};free(gotpreds);
free(pvaluesb);free(indexer);free(indexer2);

//check for invalid
for(j=0;j<length;j++)
{
if(pvalues[j]!=2&&(pvalues[j]<0||pvalues[j]>1))
{printf("Error, Predictor %s has p-value %.6f (p-values should be within [0,1])\n\n", preds[j], pvalues[j]);exit(1);}
}

return(0);
}	//end of read_pvafile

////////

int read_indhers(char *indhers, double *weights, int length, char **preds)
{
int j, count, count2, head;
int *indexer, *indexer2;

char **gotpreds;
double *weightsb;


head=check_head(indhers,"Predictor","SNP",0);
count=countrows(indhers)-head;

gotpreds=malloc(sizeof(char*)*count);
weightsb=malloc(sizeof(double)*count);
indexer=malloc(sizeof(int)*length);
indexer2=malloc(sizeof(int)*length);

read_strings(indhers, gotpreds, count, NULL, 1, head);
read_values(indhers, weightsb, count, NULL, 2, head, 0);
count2=find_strings(preds, length, gotpreds, count, indexer, indexer2, NULL, indhers, NULL, NULL, 3);
if(count2==0)
{printf("Error, %s does not contain heritabilities for any of the %d predictors\n\n", indhers, length);exit(1);}
if(count2<length)
{printf("Warning, %s contains heritabilities for only %d of the %d predictors; the remaining predictors will be given heritability zero\n\n", indhers, count2, length);}
for(j=0;j<length;j++){weights[j]=0;}
for(j=0;j<count2;j++){weights[indexer[j]]=weightsb[indexer2[j]];}
for(j=0;j<count;j++){free(gotpreds[j]);};free(gotpreds);
free(weightsb);free(indexer);free(indexer2);

return(0);
}	//end of read_indhers

///////////////////////////

int read_regions(int **regindex, int *rkeeppreds, int num_regs, char *regpref, int num_preds, char **allpreds, int *predorder, int num_preds_use, int *keeppreds)
{
int j, r, count, count2, count3, wcount;
int *usedpreds;

char **wantpreds;

char filename[500];


usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds[j]]=1;}

for(r=0;r<num_regs;r++)
{
sprintf(filename,"%s%d",regpref, r+1);
count=countrows(filename);
printf("Reading %d region predictors from %s\n", count, filename);
regindex[r]=malloc(sizeof(int)*(1+count));
wantpreds=malloc(sizeof(char*)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);

count2=find_strings(allpreds, num_preds, wantpreds, count, regindex[r]+1, NULL, NULL, filename, predorder, NULL, 3);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the data\n", count2);}

//allow for predictor filterings
count3=0;
for(j=0;j<count2;j++)
{
if(usedpreds[regindex[r][1+j]]==1){regindex[r][1+count3]=regindex[r][1+j];count3++;}
}
if(count3==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
if(count3<count2){printf("Warning, will be using only %d of these\n", count3);}
regindex[r][0]=count3;

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
}	//end of r loop
printf("\n");

//see which predictors being used
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(r=0;r<num_regs;r++)
{
for(j=0;j<regindex[r][0];j++){usedpreds[regindex[r][1+j]]++;}
}

//fill up rkeeppreds then squeeze down regindex
count=0;wcount=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]>1)
{
if(wcount<5){printf("Warning, Predictor %s is in more than one region\n", allpreds[j]);}
wcount++;
}
if(usedpreds[j]>0){rkeeppreds[count]=j;usedpreds[j]=count;count++;}
}
if(wcount>5){printf("In total, %d predictors are in more than one region\n", wcount);}
if(wcount>0){printf("\n");}

count2=0;
for(r=0;r<num_regs;r++)
{
for(j=0;j<regindex[r][0];j++){regindex[r][1+j]=usedpreds[regindex[r][1+j]];}
count2+=regindex[r][0];
}

if(num_regs>1)
{
if(count==count2){printf("The %d regions contain %d predictors\n\n", num_regs, count);}
else{printf("The %d regions contain %d predictors (%d unique)\n\n", num_regs, count2, count);}
}

free(usedpreds);

return(count);
}	//end of read_regions

///////////////////////////

void check_respfile(char *respfile, int *usedids, int num_samples, char **allids3, int num_resps_use, int *keepresps, int num_resps, double missingvalue)
{
int i, m, count, count2, count3, head;

int *founds, *indexer, *indexer2;
char **wantids;

double *resptemp;

char readchar, readstring[500];

FILE *input;


//want to tally how many responses present for each sample
founds=malloc(sizeof(int)*num_samples);
for(i=0;i<num_samples;i++){founds[i]=0;}

//see which individuals are available
head=check_head_ids(respfile,0);
count=countrows(respfile)-head;
printf("Checking responses for %d samples from %s\n", count, respfile);
wantids=malloc(sizeof(char*)*count);
read_ids(respfile, NULL, NULL, wantids, count, NULL, head);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, allids3, num_samples, indexer, indexer2, respfile, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find any of these samples\n\n");exit(1);}

//open respfile and skip the header row if present
if((input=fopen(respfile,"r"))==NULL)
{printf("Error opening %s\n\n", respfile);exit(1);}
if(head==1)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

resptemp=malloc(sizeof(double)*num_resps);
count3=0;
for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", readstring, readstring)!=2)
{printf("Error reading IDs on row %d of %s\n\n", i+1, respfile);exit(1);}

if(i==indexer[count3])	//using this row
{
for(m=0;m<num_resps;m++)
{
if(fscanf(input, "%[0-9eE.+-] ", readstring)==1){resptemp[m]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 3+m, head+i+1, respfile);exit(1);}
if(strcmp(readstring,"NA")==0){resptemp[m]=missingvalue;}
else{printf("Error reading %s; Element %d of Row %d is unrecognisable (%s)\n\n", respfile, 3+m, head+i+1, readstring);exit(1);}
}
}
for(m=0;m<num_resps_use;m++){founds[indexer2[count3]]+=(resptemp[keepresps[m]]!=missingvalue);}
count3++;
if(count3==count2){break;}
}
else	//not using this row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(resptemp);

//see if valid and whether losing any
for(i=0;i<num_samples;i++)
{
if(usedids[i]==1&&founds[i]!=0&&founds[i]!=num_resps_use)	//so must have num_resps_use>1
{printf("Error, to analyse multiple phenotypes, each sample must have either all phenotypes present or all missing; you should instead either analyse each phenotype separately, or use \"--dentist YES\" to pad missing values\n\n");exit(1);}
}

count=0;for(i=0;i<num_samples;i++){count+=usedids[i];}
count2=0;for(i=0;i<num_samples;i++){usedids[i]=(usedids[i]==1&&founds[i]>0);count2+=usedids[i];}
if(count2==0){printf("Error, all phenotypes are missing\n\n");exit(1);}
if(count2<count)
{printf("Due to missing phenotypic values, the number of samples is reduced from %d to %d\n", count, count2);}

free(founds);
}	//end of check_respfile

////////

int read_respfile(char *respfile, double *resp, int ns, char **ids3, int num_resps_use, int *keepresps, int num_resps, int *respcounts, double missingvalue, double binary, int pad)
{
int i, m, count, count2, count3, indcount, head;
int n0, n1, n2, n9, nm, nr;
double sum, sumsq, mean, var;

int *indexer, *indexer2;
char **wantids;

double *resptemp;

char readchar, readstring[500];

FILE *input;


//see which individuals are available
head=check_head_ids(respfile,0);
count=countrows(respfile)-head;
printf("Reading responses for %d samples from %s\n", count, respfile);
wantids=malloc(sizeof(char*)*count);
read_ids(respfile, NULL, NULL, wantids, count, NULL, head);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, respfile, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find any of these samples\n\n");exit(1);}

//set responses to missing
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<ns;i++){resp[i+m*ns]=missingvalue;}
}

//open respfile and skip the header row if present
if((input=fopen(respfile,"r"))==NULL)
{printf("Error opening %s\n\n", respfile);exit(1);}
if(head==1)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

resptemp=malloc(sizeof(double)*num_resps);
count3=0;
for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", readstring, readstring)!=2)
{printf("Error reading IDs on row %d of %s\n\n", i+1, respfile);exit(1);}

if(i==indexer[count3])	//using this row
{
for(m=0;m<num_resps;m++)
{
if(fscanf(input, "%[0-9eE.+-] ", readstring)==1){resptemp[m]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 3+m, head+i+1, respfile);exit(1);}
if(strcmp(readstring,"NA")==0){resptemp[m]=missingvalue;}
else
{printf("Error reading %s; Element %d of Row %d is unrecognisable (%s)\n\n", respfile, 3+m, head+i+1, readstring);exit(1);}
}
}
for(m=0;m<num_resps_use;m++){resp[indexer2[count3]+m*ns]=resptemp[keepresps[m]];}
count3++;
if(count3==count2){break;}
}
else	//not using this row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(resptemp);

//perform some checks
for(m=0;m<num_resps_use;m++)
{
n0=0;n1=0;n2=0;n9=0;nm=0;
for(i=0;i<ns;i++)
{
if(resp[i+m*ns]==0){n0++;}
if(resp[i+m*ns]==1){n1++;}
if(resp[i+m*ns]==2){n2++;}
if(resp[i+m*ns]==-9){n9++;}
if(resp[i+m*ns]==missingvalue){nm++;}
}
nr=ns-n0-n1-n2-nm;

if(n9>0)
{printf("Warning, some samples have value -9 for Phenotype %d; LDAK does NOT treat these as missing\n\n", keepresps[m]+1);}

if(binary==1&&nr>0)
{printf("Error, Phenotype %d is not binary\n\n", keepresps[m]+1);exit(1);}
if(binary==1&&n0>0&&n1>0&&n2>0)
{printf("Error, Phenotype %d takes three different values (note that LDAK does NOT treat 0 as missing)\n\n", keepresps[m]+1);exit(1);}

if(nr==0&&n2==0)	//cases=1, controls=0
{
if(m<10)
{
printf("Response %d in %s is binary, with %d cases and %d controls", keepresps[m]+1, respfile, n1, n0);
if(nm>0){printf(", %d missing", nm);}
printf("\n");
}
}
if(nr==0&&n2>0)	//cases=2, controls=1 - subtract one to get cases=1, controls=0
{
if(m<10)
{
printf("Response %d in %s is binary, with %d cases and %d controls", keepresps[m]+1, respfile, n2, n1);
if(nm>0){printf(" (and %d missing)", nm);}
printf("\n");
}
for(i=0;i<ns;i++)
{
if(resp[i+m*ns]!=missingvalue){resp[i+m*ns]--;}
}
}
}	//end of m loop

//check counts and variances, and possibly pad
count=0;
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;indcount=0;mean=0;var=0;
for(i=0;i<ns;i++)
{
if(resp[i+m*ns]!=missingvalue)
{sum+=resp[i+m*ns];sumsq+=pow(resp[i+m*ns],2);indcount++;}
}

if(indcount<3)
{printf("Error reading %s; Phenotype %d has only %d non-missing values\n\n", respfile, keepresps[m]+1, indcount);exit(1);}

mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0)
{printf("Error reading %s; Phenotype %d has variance 0\n\n", respfile, keepresps[m]+1);exit(1);}
respcounts[m]=indcount;

if(pad==1&&indcount<ns)	//pad remainder
{
for(i=0;i<ns;i++)
{
if(resp[i+m*ns]==missingvalue){resp[i+m*ns]=mean;}
}
count++;
}

if(pad==1&&m<10)
{
if(indcount<ns)
{printf("%d out of %d samples have values for Phenotype %d (remainder will be set to the mean)\n", indcount, ns, m+1);}
else
{printf("All %d samples have values for Phenotype %d (no need to pad)\n", ns,  m+1);}
}
}
printf("\n");

if(num_resps_use>1&&pad==1)
{printf("In total, %d of the %d phenotypes have missing values padded\n\n", count, num_resps_use);}

return(0);
}	//end of read_respfile

///////////////////////////

int read_covarfile(char *covarfile, double *covar, int ns, char **ids3, int num_covars, double missingvalue)
{
int i, j, j2, count, count2, count3, indcount, head;
double sum, sumsq, mean, var, alpha, beta;

int *indexer, *indexer2;
char **wantids;
double *stands, *cors;

char readchar, readstring[500];

FILE *input;


//see which individuals are available (now require all to be present)
head=check_head_ids(covarfile,0);
count=countrows(covarfile)-head;
printf("Reading %d covariates for %d samples from %s\n", num_covars, count, covarfile);
wantids=malloc(sizeof(char*)*count);
read_ids(covarfile, NULL, NULL, wantids, count, NULL, head);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, covarfile, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find any of these samples\n\n");exit(1);}
if(count2<ns){printf("Error, %s only contains covariates for %d of the %d samples\n\n", covarfile, count2, ns);exit(1);}

//set covariates to missing
for(j=0;j<num_covars;j++)
{
for(i=0;i<ns;i++){covar[i+j*ns]=missingvalue;}
}

//open covarfile and skip the header row if present
if((input=fopen(covarfile,"r"))==NULL)
{printf("Error opening %s\n\n", covarfile);exit(1);}
if(head==1)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

count3=0;
for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", readstring, readstring)!=2)
{printf("Error reading IDs on row %d of %s\n\n", i+1, covarfile);exit(1);}

if(i==indexer[count3])	//using this row
{
for(j=0;j<num_covars;j++)
{
if(fscanf(input, "%[0-9eE.+-] ", readstring)==1)
{covar[indexer2[count3]+j*ns]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 3+j, head+i+1, covarfile);exit(1);}
if(strcmp(readstring,"NA")!=0){printf("Error reading %s; Element %d of Row %d is unrecognisable (%s)\n\n", covarfile, 3+j, head+i+1, readstring);exit(1);}
}
}
count3++;
if(count3==count2){break;}
}
else	//not using this row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}
fclose(input);

////////

//get counts, means and variances, then set missing to mean and fill stands
stands=malloc(sizeof(double)*ns*num_covars);
count2=0;
for(j=0;j<num_covars;j++)
{
sum=0;sumsq=0;indcount=0;mean=0;var=0;
for(i=0;i<ns;i++)
{
if(covar[i+j*ns]!=missingvalue)
{sum+=covar[i+j*ns];sumsq+=pow(covar[i+j*ns],2);indcount++;}
}

if(indcount==0){printf("Error reading %s; Covariate %d has no non-missing values\n\n", covarfile, j+1);exit(1);}
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0){printf("Error reading %s; Covariate %d has variance 0\n\n", covarfile, j+1);exit(1);}

for(i=0;i<ns;i++)
{
if(covar[i+j*ns]==missingvalue){covar[i+j*ns]=mean;count2++;}
stands[i+j*ns]=(covar[i+j*ns]-mean)*pow(var*indcount,-.5);
}}

if(count2>0){printf("There are %d missing values; these are replaced by the mean value of the corresponding covariate\n", count2);}
printf("\n");

//check correlations
cors=malloc(sizeof(double)*num_covars*num_covars);
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, stands, &ns, stands, &ns, &beta,  cors, &num_covars);

count2=0;
for(j=0;j<num_covars;j++)
{
for(j2=0;j2<j;j2++)
{
if(pow(cors[j+j2*num_covars],2)>.99)
{printf("Error, Columns %d and %d in %s are (almost) linearly dependent (correlation %.3f)\n", j2+3, j+3, covarfile, cors[j+j2*num_covars]);count2++;}
}
}
if(count2>0){printf("\n");exit(1);}


for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(stands);free(cors);

return(0);
}	//end of read_covarfile

///////////////////////////

int read_envfile(char *envfile, double *covar, int ns, char **ids3, int num_envs, int discenv, char **ids1, char **ids2)
{
int i, j, count, count2, count3, head, lwork, *iwork, info, one=1;
double sum, sumsq, mean, var, *mat, *S, wkopt, *work;

int *indexer, *indexer2, *usedids;
char **wantids;

char readchar, readstring[500];

FILE *input;


//see which individuals are available (now require all to be present)
head=check_head_ids(envfile,0);
count=countrows(envfile)-head;
printf("Reading %d environmental variables for %d samples from %s\n", num_envs, count, envfile);
wantids=malloc(sizeof(char*)*count);
read_ids(envfile, NULL, NULL, wantids, count, NULL, head);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, envfile, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find any of these samples\n\n");exit(1);} 
if(count2<ns){printf("Error, %s only contains environmental variables for %d of the %d samples\n\n", envfile, count2, ns);exit(1);}

//open envfile and skip the header row if present
if((input=fopen(envfile,"r"))==NULL)
{printf("Error opening %s\n\n", envfile);exit(1);}
if(head==1)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

count2=0;
for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", readstring, readstring)!=2)
{printf("Error reading IDs on row %d of %s\n\n", i+1, envfile);exit(1);}

if(i==indexer[count2])	//using this row
{
for(j=0;j<num_envs;j++)
{
if(fscanf(input, "%[0-9eE.+-] ", readstring)==1)
{covar[indexer2[count2]+j*ns]=atof(readstring);}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 3+j, head+i+1, envfile);exit(1);}
if(strcmp(readstring,"NA")==0){printf("Error reading %s; Element %d of Row %d is NA\n\n", envfile, 3+j, head+i+1);exit(1);}
else{printf("Error reading %s; Element %d of Row %d is unrecognisable (%s)\n\n", envfile, 3+j, head+i+1, readstring);exit(1);}
}
}
count2++;
if(count2==ns){break;}
}
else	//not using this row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}
fclose(input);

////////

//check variances
for(j=0;j<num_envs;j++)
{
sum=0;sumsq=0;mean=0;var=0;
for(i=0;i<ns;i++){sum+=covar[i+j*ns];sumsq+=pow(covar[i+j*ns],2);}
mean=sum/ns;
var=sumsq/ns-pow(mean,2);
if(var==0){printf("Error reading %s; Variable %d has variance 0\n\n", envfile, j+1);exit(1);}
}

//compute rank
S=malloc(sizeof(double)*num_envs);
mat=malloc(sizeof(double)*ns*num_envs);
for(j=0;j<num_envs;j++)
{
for(i=0;i<ns;i++){mat[i+j*ns]=covar[i+j*ns];}
}

lwork=-1;
iwork=malloc(sizeof(int)*8*num_envs);
dgesdd_("N", &ns, &num_envs, mat, &ns, S, NULL, &one, NULL, &one, &wkopt, &lwork, iwork, &info);
if(info!=0){printf("Rank error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

work=malloc(sizeof(double)*lwork);
dgesdd_("N", &ns, &num_envs, mat, &ns, S, NULL, &one, NULL, &one, work, &lwork, iwork, &info);
if(info!=0){printf("Rank error 2; please tell Doug info %d\n\n", info);exit(1);}
free(iwork);free(work);

count2=0;for(j=0;j<num_envs;j++){count2+=(fabs(S[j])>1e-6);}
if(count2<num_envs){printf("Error, the variables are linearly dependent (have rank %d)\n\n", count2);exit(1);} 

if(discenv==1)	//check binary, and all individuals in one group
{
printf("Checking the environmental variables partition the samples into subgroups\n");
usedids=malloc(sizeof(int)*ns);
for(i=0;i<ns;i++){usedids[i]=0;}

for(j=0;j<num_envs;j++)
{
for(i=0;i<ns;i++)
{
if(covar[i+j*ns]!=0&&covar[i+j*ns]!=1)
{printf("Error, Variable %d takes value %.6f (all variables should be 0 or 1)\n\n", j+1, covar[i+j*ns]);exit(1);}
usedids[i]+=covar[i+j*ns];
}
}

count2=0;count3=0;
for(i=0;i<ns;i++){count2+=(usedids[i]==0);count3+=(usedids[i]>1);}
if(count2>0){printf("Error, %d samples are not in a subgroup; consider adding an extra subgroup containing these\n", count2);}
if(count3>0){printf("Error, %d samples are in more than one subgroup\n", count3);}
if(count2+count3>0){printf("\n");exit(1);}
free(usedids);
}

printf("\n");

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(S);free(mat);

return(0);
}	//end of read_envfile

////////

int read_sampwfile(char *sampwfile, double *sweights, int ns, char **ids3, char **ids1, char **ids2)
{
int i, count, count2, head;
double *sweightsb;

int *indexer;
char **wantids;


//see which individuals are available (now require all to be present)
head=check_head_ids(sampwfile,0);
count=countrows(sampwfile)-head;
printf("Reading weightings for %d samples from %s\n", count, sampwfile);
wantids=malloc(sizeof(char*)*count);
sweightsb=malloc(sizeof(double)*count);
read_ids(sampwfile, NULL, NULL, wantids, count, NULL, head);
read_values(sampwfile, sweightsb, count, NULL, 3, head, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, NULL, indexer, sampwfile, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find any of these samples\n\n");exit(1);}
if(count2<ns){printf("Error, %s only contains weightings for %d of the %d samples (if this is intentional, use \"--keep\" or \"--remove\" to exclude the missing samples from the analysis)\n\n", sampwfile, count2, ns);exit(1);}

for(i=0;i<ns;i++){sweights[i]=sweightsb[indexer[i]];}

for(i=0;i<ns;i++)
{
if(sweights[i]<=0){printf("Error, Sample %s %s has weighting %.4f (all weightings must be positive)\n\n", ids1[i], ids2[i], sweights[i]);exit(1);}
}
printf("\n");

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(sweightsb);

return(0);
}	//end of read_sampwfile

///////////////////////////

int read_tops(char *datafile, double *covar, int ns, char **ids3, int num_tops, int *tkeeppreds, double *tcentres, double *tvars, char *famfile, int famhead, int dtype, int binary, int genskip, int genheaders, int genprobs, int num_preds, char **allpreds, double missingvalue, int nonsnp, char *sumsfile)
{
int i, j, j2, count, count2, count3, indcount;
double sum, sumsq, mean, var, value, alpha, beta;
 
int *indexer, *indexer2;
char **wantids;
double *datatemp, *stands, *cors;

gzFile datainputgz;


//have already checked all top predictors present
printf("Reading values for %d top predictors\n", num_tops);
count=countrows(famfile)-famhead;
wantids=malloc(sizeof(char*)*count);
read_ids(famfile, NULL, NULL, wantids, count, NULL, famhead);
indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);

//read in
datatemp=malloc(sizeof(double)*count2*num_tops);
if(binary==0){open_datagz(&datainputgz, datafile, count, genskip, genheaders, genprobs);}
(void)read_data_fly(datafile, dtype, datatemp, NULL, count2, indexer, 0, num_tops, tkeeppreds, datainputgz, 0, count, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp,1);
if(binary==0){gzclose(datainputgz);}

//fill up
for(j=0;j<num_tops;j++)
{
for(i=0;i<ns;i++){covar[i+j*ns]=missingvalue;}
for(i=0;i<count2;i++){covar[indexer2[i]+j*ns]=datatemp[i+j*count2];}
}

//get counts, means and variances, then centre and fill stands
stands=malloc(sizeof(double)*ns*num_tops);
count3=0;
for(j=0;j<num_tops;j++)
{
sum=0;sumsq=0;indcount=0;mean=0;var=0;
for(i=0;i<ns;i++)
{
if(covar[i+j*ns]!=missingvalue)
{sum+=covar[i+j*ns];sumsq+=pow(covar[i+j*ns],2);indcount++;}
}

if(indcount==0)
{printf("Error, Predictor %s has no non-missing values\n\n", allpreds[tkeeppreds[j]]);exit(1);}
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0){printf("Error, Predictor %s has variance 0\n\n", allpreds[tkeeppreds[j]]);exit(1);}

for(i=0;i<ns;i++)
{
if(covar[i+j*ns]!=missingvalue){covar[i+j*ns]-=mean;}
else{covar[i+j*ns]=0;count3++;}
stands[i+j*ns]=covar[i+j*ns]*pow(var*indcount,-.5);
}

tcentres[j]=mean;
tvars[j]=var*indcount/ns;
}	//end of j loop

if(count3>0){printf("There are %d missing values; these are replaced by the mean value of the corresponding predictor\n", count3);}

//check correlations
cors=malloc(sizeof(double)*num_tops*num_tops);
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_tops, &num_tops, &ns, &alpha, stands, &ns, stands, &ns, &beta,  cors, &num_tops);

count2=0;count3=0;value=0;
for(j=0;j<num_tops;j++)
{
for(j2=0;j2<j;j2++)
{
if(pow(cors[j+j2*num_tops],2)>.99)
{
if(count2<5)
{printf("Error, Predictors %s and %s are (almost) linearly dependent (correlation %.3f)\n", allpreds[tkeeppreds[j2]], allpreds[tkeeppreds[j]], cors[j+j2*num_tops]);}
count2++;
}
if(pow(cors[j+j2*num_tops],2)>value){value=pow(cors[j+j2*num_tops],2);}
count3+=(pow(cors[j+j2*num_tops],2)>.95);
}}
if(count2>5){printf("In total, there are %d pairs of highly correlated predictors\n", count2);}
if(count2>0)
{
if(strcmp(sumsfile,"blank")!=0)
{printf("You should first prune the top predictors, for example, using the LDAK commands \"--thin\" with option \"--window-prune 0.1\"\n\n");exit(1);}
else
{printf("You should first prune the top predictors, for example, using the LDAK commands \"--thin\" with option \"--window-prune 0.5\"\n\n");exit(1);}
}

if(num_tops>1)
{
printf("The highest correlation squared between pairs of top predictors is %.6f (%d are above 0.95)\n", value, count3);
if(value>0.1&&strcmp(sumsfile,"blank")!=0)
{printf("Warning, when using summary statistics, we recommend no pair of top predictors has correlation squared higher than 0.1\n");}
}
printf("\n");

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(datatemp);free(stands);free(cors);

return(0);
}

///////////////////////////

int get_effects(double **effects, double bivar, int *keeppreds_use, int num_phenos, int num_causals, char *causalsfile, char *effectsfile, int num_preds_use, int *keeppreds, int num_preds, char **allpreds, int *predorder)
{
int j, m, count;
int **causals, *order, *usedpreds, *revs;
double unifrand, value, value2, mat[4], mat2[2];

char **wantpreds, *rs;

FILE *input;


//get causals
causals=malloc(sizeof(int*)*num_phenos);
for(m=0;m<num_phenos;m++){causals[m]=malloc(sizeof(int)*num_causals);}

if(strcmp(causalsfile,"blank")!=0)	//read in and match against allpreds (not using bivar)
{
printf("Reading causal predictors from %s\n", causalsfile);
rs=malloc(sizeof(char)*10000000);
wantpreds=malloc(sizeof(char*)*num_causals);

if((input=fopen(causalsfile,"r"))==NULL)
{printf("Error opening %s\n\n", causalsfile);exit(1);}
for(m=0;m<num_phenos;m++)
{
for(j=0;j<num_causals;j++)
{
if(fscanf(input, "%s ", rs)!=1)
{printf("Error reading Predictor %d on Row %d of %s\n\n", j+1, m+1, causalsfile);exit(1);}
copy_string(wantpreds,j,rs);
}
count=find_strings(wantpreds, num_causals, allpreds, num_preds, NULL, causals[m], causalsfile, NULL, NULL, predorder, 3);
if(count<num_causals)
{printf("Error reading %s; can only find %d of the %d predictors in Row %d\n\n", causalsfile, count, num_causals, m+1);exit(1);}
for(j=0;j<count;j++){free(wantpreds[j]);}
}
fclose(input);

free(rs);free(wantpreds);
}
else	//pick at random
{
if(num_causals==1)
{
for(m=0;m<num_phenos;m++)
{unifrand=(double)rand()/RAND_MAX;causals[m][0]=(int)(unifrand*num_preds_use);}
}
else
{
order=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){order[j]=j;}
for(m=0;m<num_phenos;m++)
{
permute_int(order,num_preds_use);
//tweaked here for correlated effects
//qsort(order, num_causals, sizeof(int), compare_int);
for(j=0;j<num_causals;j++){causals[m][j]=keeppreds[order[j]];}
}
free(order);
}

if(bivar!=-9999)	//must ensure pairs of causals are the same
{
for(m=0;m<num_phenos/2;m++)
{
for(j=0;j<num_causals;j++){causals[m*2+1][j]=causals[m*2][j];}
}
}
}

//work out which used
usedpreds=malloc(sizeof(int)*num_preds);
revs=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(m=0;m<num_phenos;m++)
{
for(j=0;j<num_causals;j++){usedpreds[causals[m][j]]++;}
}

count=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]>0){keeppreds_use[count]=j;revs[j]=count;count++;}
}

////////

//now get effects
for(m=0;m<num_phenos;m++)
{
effects[m]=malloc(sizeof(double)*count);
for(j=0;j<count;j++){effects[m][j]=0;}
}

if(strcmp(effectsfile,"blank")!=0)	//read in 
{
if((input=fopen(effectsfile,"r"))==NULL)
{printf("Error opening %s\n\n", effectsfile);exit(1);}
for(m=0;m<num_phenos;m++)
{
for(j=0;j<num_causals;j++)
{
if(fscanf(input, "%lf ", effects[m]+revs[causals[m][j]])!=1)
{printf("Error reading Effect %d on Row %d of %s\n\n", j+1, m+1, effectsfile);exit(1);}
}
}
fclose(input);
}
else	//pick at random
{
if(bivar==-9999)	//effects independent
{
for(m=0;m<num_phenos;m++)
{
for(j=0;j<num_causals;j++)
{
effects[m][revs[causals[m][j]]]=rnorm_safe();
//tweaked here for correlated effects
//if(j%2==1){effects[m][revs[causals[m][j]]]=effects[m][revs[causals[m][j-1]]];}
}
}
}
else	//jointly sample effect sizes for pairs of traits
{
//get "cholesky"
if(bivar==-1){mat[0]=1;mat[1]=-1;mat[2]=0;mat[3]=0;}
if(bivar==1){mat[0]=1;mat[1]=1;mat[2]=0;mat[3]=0;}
if(bivar!=-1&&bivar!=1)	
{
mat[0]=1;mat[1]=bivar;mat[2]=bivar;mat[3]=1;
eigen_invert(mat, 2, mat2, 0, NULL, 0);
}

for(m=0;m<num_phenos/2;m++)
{
for(j=0;j<num_causals;j++)
{
value=rnorm_safe();
value2=rnorm_safe();
effects[m*2][revs[causals[m*2][j]]]=mat[0]*value+mat[2]*value2;
effects[m*2+1][revs[causals[m*2+1][j]]]=mat[1]*value+mat[3]*value2;
}
}
}	//end of bivar!=-9999
}	//end of picking at random

for(m=0;m<num_phenos;m++){free(causals[m]);}free(causals);
free(usedpreds);free(revs);

return(count);
}	//end of get_effects

///////////////////////////

int get_her_model(int num_parts, char *partpref, int **pindexes, double **pweights, int length, int *keeppreds_use, int num_preds, char **allpreds, int *predorder, int parttype, int backpart, int allone)
{
int j, j2, q, count, count2, count3;
int addpart, *usedpreds, *indexer, *indexer2;

double *avalues;
char **wantpreds;

char filename[500];


for(q=0;q<num_parts+1;q++)
{
pindexes[q]=malloc(sizeof(int)*length);
pweights[q]=malloc(sizeof(double)*length);
for(j=0;j<length;j++){pindexes[q][j]=0;pweights[q][j]=0;}
}

if(strcmp(partpref,"blank")==0)	//simple case of just the base
{
addpart=1;
for(j=0;j<length;j++){pindexes[0][j]=1;pweights[0][j]=1;}
}
else	//fill predictor lists and sort base/background
{
usedpreds=malloc(sizeof(int)*num_preds);

for(q=0;q<num_parts;q++)
{
sprintf(filename, "%s%d", partpref, q+1);
count=countrows(filename);
count3=countcols(filename);
if(count3==2){printf("Reading %d predictors and weights from %s\n", count, filename);}
else{printf("Reading %d predictors from %s\n", count, filename);}
wantpreds=malloc(sizeof(char*)*count);
avalues=malloc(sizeof(double)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);
if(count3==2){read_values(filename, avalues, count, NULL, 2, 0, 0);}
else
{
for(j=0;j<count;j++){avalues[j]=1;}
}

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(allpreds, num_preds, wantpreds, count, indexer, indexer2, NULL, filename, predorder, NULL, 3);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the data\n", count2);}

for(j=0;j<num_preds;j++){usedpreds[j]=-1;}
for(j=0;j<count2;j++){usedpreds[indexer[j]]=indexer2[j];}
count2=0;
for(j=0;j<length;j++)
{
j2=keeppreds_use[j];
if(usedpreds[j2]!=-1)	//then the j2th predictor is in the partition
{pindexes[q][j]=1;pweights[q][j]=avalues[usedpreds[j2]];count2++;}
}
if(count2==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);free(indexer2);free(avalues);
}
printf("\n");

//see which predictors are being used
for(j=0;j<length;j++)
{
usedpreds[j]=0;for(q=0;q<num_parts;q++){usedpreds[j]+=pindexes[q][j];}
}
count=0;count2=0;
for(j=0;j<length;j++){count+=(usedpreds[j]>0);count2+=(usedpreds[j]==1);}

if(parttype==0)	//checks annotations are sensible and add base
{
if(count2==length){printf("Warning, all %d predictors are in exactly one annotation; consider using \"--partition-number\" and \"--partition-prefix\" instead of \"--annotation-number\" and \"--annotation-prefix\"\n\n", length);}

addpart=1;
for(j=0;j<length;j++){pindexes[num_parts][j]=1;pweights[num_parts][j]=1;}
}
else	//checks for partitions (adjust weights if necessary and decide background)
{
if(count>count2)	//some predictors are in more than one partition
{
if(allone==0)	//share across overlapping predictors
{
printf("There are %d predictors in more than one partition; a predictor in X partitions will be allocated 1/X to each (to instead allocate predictors fully to each partition to which they belong, use \"--all-one\" YES)\n\n", count-count2);
for(j=0;j<length;j++)
{
if(usedpreds[j]>0)	//must be non-background
{
for(q=0;q<num_parts;q++){pweights[q][j]=pweights[q][j]/usedpreds[j];}
}}
}
else
{printf("There are %d predictors in more than one partition; each predictor will be allocated fully to all partition it belongs (the alternative is to use \"--all-one\" NO, in which case a predictor in X categories will be allocated 1/X to each)\n\n", count-count2);}
}
else{printf("No predictor is in more than one partition\n");}

if(count==length)	//no extra predictors
{
printf("All %d predictors are in at least one partition\n\n", length);
addpart=0;
}
else	//some extra predictors
{
if(backpart==-9999)
{printf("Error, the %d partition lists contain only %d of the %d predictors; you should either use \"--background NO\" to continue regardless, or use \"--background YES\" to create an extra partition containing the remaining %d predictors\n\n", num_parts, count, length, length-count);exit(1);}

printf("%d of the %d predictors are in at least one partition\n", count, length);
if(backpart==1)	//adding extra partition
{
printf("Will add a background partition containing the remaining %d predictors\n\n", length-count);
addpart=1;
for(j=0;j<length;j++)
{
if(usedpreds[j]==0){pindexes[num_parts][j]=1;pweights[num_parts][j]=1;}
}
}
else	//not adding extra partition, will deal with outside
{
printf("The remaining %d predictors are assumed not to contribute heritability\n\n", length-count);
addpart=2;
}
}	//end of count<length
}	//end of parttype=1

free(usedpreds);
}	//end of non-simple case

return(addpart);
}	//end of get_her_model

///////////////////////////

