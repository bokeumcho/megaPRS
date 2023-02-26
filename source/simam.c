//library includes

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <zlib.h>
#include <limits.h>

//will export functions from fortran

extern void dgemm_();
extern void sgemm_();
extern void dgemv_();
extern void dpotrf_();
extern void spotrf_();
extern void dpotri_();
extern void spotri_();
extern void dpotrs_();
extern void spotrs_();
extern void dsytrf_();
extern void ssytrf_();
extern void dsytri_();
extern void ssytri_();
extern void dsytrs_();
extern void ssytrs_();
extern void dsyev_();
extern void ssyev_();
extern void dsyevx_();
extern void dgesvd_();
extern void dgesdd_();

struct sorting_double{double value;int index;};

int compare_sorting_double (const void *a, const void *b)
{
  const struct sorting_double fa = *(const struct sorting_double *) a;
  const struct sorting_double fb = *(const struct sorting_double *) b;
  return (fa.value > fb.value) - (fa.value < fb.value);
}

#include "norm.c"

///////////////////////////

int main (int argc, const char * argv[])
{
//this line makes the buffer flush
setlinebuf(stdout);

//set random seed
srand((unsigned)clock()+(unsigned)time(NULL)+(unsigned)getpid());

//set seed for normal generator
zigset(rand());

//declarations

int nsamps=-9999, nsnps=-9999;
double minmaf=-9999, maxmaf=-9999;
char freqfile[500]="blank";
double power=-9999, her=-9999;
int gener=-9999, dist=-9999, pos=0;
double prop=-9999;
char outfile[500]="blank";

int i, j, g, gen, gen2, count, count2, found, one=1;
double sum, sum2, sumsq, sumsq2, sumsq3, mean, mean2, var, var2, var3;
double value, unifrand, alpha, beta;

double *mafs, *effects, *breed1, *breed2, *phens;
char *data;

int *ranks, low, high, pick, g1, g2;
int *p1, *p2;
char *data2;
struct sorting_double *sptrs;

char filename[500]="blank";
FILE *input, *output;


//get values

if(argc%2!=1)
{printf("Error, there must be an even number of arguments (not %d); arguments should be provided in pairs\n\n", argc-1);exit(1);}

count=1;
while(count<argc)
{
if(argv[count][0]!='-'||argv[count][1]!='-')
{printf("Error, Argument %d is incorrect (%s); all odd arguments should begin with \"--\"\n\n", count, argv[count]);exit(1);}
if(strcmp(argv[count+1],"blank")==0)
{printf("Error, none of the even arguments can be the word \"blank\"; the universe will probably implode or something\n\n");exit(1);}
count+=2;
}

for(count=1;count<argc;count+=2)
{
for(count2=count+2;count2<argc;count2+=2)
{
if(strcmp(argv[count],argv[count2])==0)
{printf("Error, Arguments %d and %d are the same (%s)\n\n", count, count2, argv[count]);exit(1);}
}
}

if(argc==3){printf("There is one pair of arguments:\n");}
else{printf("There are %d pairs of arguments:\n", (argc-1)/2);}
count=1;
while(count<argc)
{printf("%s %s\n", argv[count], argv[count+1]);count+=2;}
printf("\n");

count=1;
while(count<argc)
{
//see which argument specified, and get the corresponding number/file
found=0;

if(strcmp(argv[count],"--nsamps")==0)
{nsamps=atoi(argv[count+1]);found=1;}

if(strcmp(argv[count],"--nsnps")==0)
{nsnps=atoi(argv[count+1]);found=1;}

if(strcmp(argv[count],"--minmaf")==0)
{minmaf=atof(argv[count+1]);found=1;}

if(strcmp(argv[count],"--maxmaf")==0)
{maxmaf=atof(argv[count+1]);found=1;}

if(strcmp(argv[count],"--frequencies")==0)
{strcpy(freqfile,argv[count+1]);found=1;}

if(strcmp(argv[count],"--power")==0)
{power=atof(argv[count+1]);found=1;}

if(strcmp(argv[count],"--her")==0)
{her=atof(argv[count+1]);found=1;}

if(strcmp(argv[count],"--gener")==0)
{gener=atoi(argv[count+1]);found=1;}

if(strcmp(argv[count],"--proportion")==0)
{prop=atof(argv[count+1]);found=1;}

if(strcmp(argv[count],"--output")==0)
{strcpy(outfile,argv[count+1]);found=1;}

if(strcmp(argv[count],"--positive")==0)
{pos=atoi(argv[count+1]);found=1;}

if(found==0){printf("Error, %s is not a recognised argument\n",argv[count]);exit(1);}
count+=2;
}

if(nsamps==-9999){printf("Error, you must provide nsamps\n\n");exit(1);}
if(nsnps==-9999){printf("Error, you must provide nsnps\n\n");exit(1);}

if(strcmp(freqfile,"blank")==0)
{
if(minmaf==-9999||maxmaf==-9999){printf("Error, you must provide minmaf and maxmaf\n\n");exit(1);}
}
else
{
if(minmaf!=-9999||maxmaf!=-9999){printf("Error, there is no need to provide minmaf or maxmaf when using frequencies\n\n");exit(1);}
}

if(power==-9999){printf("Error, you must provide power\n\n");exit(1);}
if(her==-9999){printf("Error, you must provide her\n\n");exit(1);}
if(gener==-9999){printf("Error, you must provide gener\n\n");exit(1);}
if(prop==-9999){printf("Error, you must provide proportion\n\n");exit(1);}
if(strcmp(outfile,"blank")==0){printf("Error, you must provide output\n\n");exit(1);}

dist=prop*nsamps;

////////

value=(double)2*nsamps/1024*nsnps/1024/1024;
printf("Warning, to perform the analysis will require approximately %.1f Gb\n\n", value);

//get mafs
mafs=malloc(sizeof(double)*nsnps);

if(strcmp(freqfile,"blank")!=0)	//read from file
{
if((input=fopen(freqfile,"r"))==NULL)
{printf("Error opening %s\n\n", freqfile);exit(1);}
for(j=0;j<nsnps;j++)
{
if(fscanf(input, "%lf ", mafs+j)!=1)
{printf("Error reading Element %d from %s\n\n", j+1, freqfile);exit(1);}
if(mafs[j]>1||mafs[j]<0){printf("Error, MAF %d is %f\n", j+1, mafs[j]);exit(1);}
if(mafs[j]>.5||mafs[j]*nsamps<5){printf("Warning, MAF %d is %f\n", j+1, mafs[j]);}
if(mafs[j]>.5){mafs[j]=1-mafs[j];}
}
fclose(input);
}
else
{
for(j=0;j<nsnps;j++){mafs[j]=(double)rand()/RAND_MAX*(maxmaf-minmaf)+minmaf;}
}

//generate snp data
data=malloc(sizeof(char)*nsamps*nsnps);

for(j=0;j<nsnps;j++)
{
if(j%10000==0){printf("Making snp %d out of %d\n", j+1, nsnps);}

for(i=0;i<nsamps;i++)
{
unifrand=(double)rand()/RAND_MAX;
data[i+j*nsamps]=2;		//genotype 2
if(unifrand<1-pow(mafs[j],2)){data[i+j*nsamps]=1;}	//genotype 1
if(unifrand<pow(1-mafs[j],2)){data[i+j*nsamps]=0;}	//genotype 0
}
}

//get effect sizes
effects=malloc(sizeof(double)*nsnps);
for(j=0;j<nsnps;j++)
{
effects[j]=rand_safe()*pow(2*mafs[j]*(1-mafs[j]),power/2);
if(pos==1&&effects[j]<0){effects[j]=-effects[j];}
}

//get breeding values
breed1=malloc(sizeof(double)*nsamps);
breed2=malloc(sizeof(double)*nsamps);

for(i=0;i<nsamps;i++){breed1[i]=0;breed2[i]=0;}
for(j=0;j<nsnps;j++)
{
if(mafs[j]<0.01)
{
for(i=0;i<nsamps;i++){breed1[i]+=effects[j]*data[i+j*nsamps];}
}
else
{
for(i=0;i<nsamps;i++){breed2[i]+=effects[j]*data[i+j*nsamps];}
}
}

//get correlation
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(i=0;i<nsamps;i++)
{
sum+=breed1[i];sum2+=breed2[i];
sumsq+=pow(breed1[i],2);
sumsq2+=pow(breed2[i],2);
sumsq3+=breed1[i]*breed2[i];
}
mean=sum/nsamps;
mean2=sum2/nsamps;
var=sumsq/nsamps-mean*mean;
var2=sumsq2/nsamps-mean2*mean2;
var3=sumsq3/nsamps-mean*mean2;

printf("breeding cors %f and relative hers %f %f\n", var3*pow(var*var2,-.5), var, var2);

//now add to get phens
phens=malloc(sizeof(double)*nsamps);
for(i=0;i<nsamps;i++){phens[i]=breed1[i]+breed2[i];}

//get variance
sum=0;sumsq=0;
for(i=0;i<nsamps;i++){sum+=phens[i];sumsq+=pow(phens[i],2);}
mean=sum/nsamps;
var=sumsq/nsamps-mean*mean;

//scale to have variance her
value=pow(her/var,.5);
for(j=0;j<nsnps;j++){effects[j]*=value;}
for(i=0;i<nsamps;i++){phens[i]*=value;breed1[i]*=value;breed2[i]*=value;}

//add on noise
value=pow(1-her,.5);
for(i=0;i<nsamps;i++){phens[i]+=rand_safe()*value;}

//get variance and heritability
sum=0;sumsq=0;
for(i=0;i<nsamps;i++){sum+=phens[i];sumsq+=pow(phens[i],2);}
mean=sum/nsamps;
var2=sumsq/nsamps-mean*mean;

sum=0;sumsq=0;
for(i=0;i<nsamps;i++){sum+=phens[i]-breed1[i]-breed2[i];sumsq+=pow(phens[i]-breed1[i]-breed2[i],2);}
mean=sum/nsamps;
var=sumsq/nsamps-mean*mean;

printf("Starting var and her are %f and %f\n", var2, (var2-var)/var2);

//save data
printf("Saving starting data\n\n");

sprintf(filename, "%s.start.bed", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
gen=108;fwrite(&gen, sizeof(unsigned char), 1, output);
gen=27;fwrite(&gen, sizeof(unsigned char), 1, output);
gen=1;fwrite(&gen, sizeof(unsigned char), 1, output);

for(j=0;j<nsnps;j++)
{
gen=0;
for(i=0;i<nsamps;i++)
{
gen2=0;					//genotype 2
if(data[i+j*nsamps]==1){gen2=2;}	//genotype 1
if(data[i+j*nsamps]==0){gen2=3;}	//genotype 0
gen+=gen2<<(2*(i%4));
if(i%4==3||i==nsamps-1){fwrite(&gen, sizeof(unsigned char), 1, output);gen=0;}
}
}
fclose(output);

sprintf(filename, "%s.start.bim", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(j=0;j<nsnps;j++){fprintf(output,"1\tSNP%d\t0\t%d\tA\tB\n", j+1, 1000*(j+1));}
fclose(output);

sprintf(filename, "%s.start.fam", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(i=0;i<nsamps;i++){fprintf(output,"IND%d IND%d 0 0 0 0\n", i+1, i+1);}
fclose(output);

sprintf(filename, "%s.start.pheno", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(i=0;i<nsamps;i++){fprintf(output,"IND%d IND%d %.6f\n", i+1, i+1, phens[i]);}
fclose(output);

////////

ranks=malloc(sizeof(int)*nsamps);
p1=malloc(sizeof(int)*nsamps);
p2=malloc(sizeof(int)*nsamps);

data2=malloc(sizeof(double)*nsamps*nsnps);

sptrs=malloc(sizeof(struct sorting_double)*nsamps);

for(g=0;g<gener;g++)
{
printf("Generation %d of %d\n", g+1, gener);
//order the phenotypes and get ranks
for(i=0;i<nsamps;i++){sptrs[i].value=phens[i];sptrs[i].index=i;}
qsort(sptrs, nsamps, sizeof(struct sorting_double), compare_sorting_double);
for(i=0;i<nsamps;i++){ranks[sptrs[i].index]=i;}

//get nsamps pairs - pick first at random, second based on similarity
for(i=0;i<nsamps;i++)
{
p1[i]=(double)rand()/RAND_MAX*nsamps;

low=ranks[p1[i]]-dist;if(low<0){low=0;}
high=ranks[p1[i]]+dist;if(high>=nsamps){high=nsamps;}

p2[i]=p1[i];
while(p2[i]==p1[i])
{
pick=low+(double)rand()/RAND_MAX*(high-low);
p2[i]=sptrs[pick].index;
}
}

//get correlation
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(i=0;i<nsamps;i++)
{
sum+=phens[p1[i]];sum2+=phens[p2[i]];
sumsq+=pow(phens[p1[i]],2);
sumsq2+=pow(phens[p2[i]],2);
sumsq3+=phens[p1[i]]*phens[p2[i]];
}
mean=sum/nsamps;
mean2=sum2/nsamps;
var=sumsq/nsamps-mean*mean;
var2=sumsq2/nsamps-mean2*mean2;
var3=sumsq3/nsamps-mean*mean2;

printf("cors %d is %f\n", g+1, var3*pow(var*var2,-.5));

//make new data
#pragma omp parallel for private(j,i,g1,g2,unifrand) schedule(static)
for(j=0;j<nsnps;j++)
{
for(i=0;i<nsamps;i++)
{
g1=data[p1[i]+j*nsamps];
g2=data[p2[i]+j*nsamps];

data2[i+j*nsamps]=0;

if(g1==1)
{
unifrand=(double)rand()/RAND_MAX;
if(unifrand>.5){data2[i+j*nsamps]++;}
}
if(g1==2){data2[i+j*nsamps]++;}

if(g2==1)
{
unifrand=(double)rand()/RAND_MAX;
if(unifrand>.5){data2[i+j*nsamps]++;}
}
if(g2==2){data2[i+j*nsamps]++;}
}}

#pragma omp parallel for private(j,i) schedule(static)
for(j=0;j<nsnps;j++)
{
for(i=0;i<nsamps;i++){data[i+j*nsamps]=data2[i+j*nsamps];}
}

//now new breeding values
for(i=0;i<nsamps;i++){breed1[i]=0;breed2[i]=0;}
for(j=0;j<nsnps;j++)
{
if(mafs[j]<0.01)
{
for(i=0;i<nsamps;i++){breed1[i]+=effects[j]*data[i+j*nsamps];}
}
else
{
for(i=0;i<nsamps;i++){breed2[i]+=effects[j]*data[i+j*nsamps];}
}
}

//get correlation
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(i=0;i<nsamps;i++)
{
sum+=breed1[i];sum2+=breed2[i];
sumsq+=pow(breed1[i],2);
sumsq2+=pow(breed2[i],2);
sumsq3+=breed1[i]*breed2[i];
}
mean=sum/nsamps;
mean2=sum2/nsamps;
var=sumsq/nsamps-mean*mean;
var2=sumsq2/nsamps-mean2*mean2;
var3=sumsq3/nsamps-mean*mean2;

printf("breeding cors %f and relative hers %f %f\n", var3*pow(var*var2,-.5), var, var2);

//combine breeding values and add on noise
value=pow(1-her,.5);
for(i=0;i<nsamps;i++){phens[i]=breed1[i]+breed2[i]+rand_safe()*value;}

//get variance and heritability
sum=0;sumsq=0;
for(i=0;i<nsamps;i++){sum+=phens[i];sumsq+=pow(phens[i],2);}
mean=sum/nsamps;
var2=sumsq/nsamps-mean*mean;

sum=0;sumsq=0;
for(i=0;i<nsamps;i++){sum+=phens[i]-breed1[i]-breed2[i];sumsq+=pow(phens[i]-breed1[i]-breed2[i],2);}
mean=sum/nsamps;
var=sumsq/nsamps-mean*mean;

printf("Generation %d var and her are %f and %f\n", g+1, var2, (var2-var)/var2);
}	//end of g loop

//save data
printf("Saving final data\n\n");

sprintf(filename, "%s.end.bed", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
gen=108;fwrite(&gen, sizeof(unsigned char), 1, output);
gen=27;fwrite(&gen, sizeof(unsigned char), 1, output);
gen=1;fwrite(&gen, sizeof(unsigned char), 1, output);

for(j=0;j<nsnps;j++)
{
gen=0;
for(i=0;i<nsamps;i++)
{
gen2=0;					//genotype 2
if(data[i+j*nsamps]==1){gen2=2;}	//genotype 1
if(data[i+j*nsamps]==0){gen2=3;}	//genotype 0
gen+=gen2<<(2*(i%4));
if(i%4==3||i==nsamps-1){fwrite(&gen, sizeof(unsigned char), 1, output);gen=0;}
}
}
fclose(output);

sprintf(filename, "%s.end.bim", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(j=0;j<nsnps;j++){fprintf(output,"1\tSNP%d\t0\t%d\tA\tB\n", j+1, 1000*(j+1));}
fclose(output);

sprintf(filename, "%s.end.fam", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(i=0;i<nsamps;i++){fprintf(output,"IND%d IND%d 0 0 0 0\n", i+1, i+1);}
fclose(output);

sprintf(filename, "%s.end.pheno", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(i=0;i<nsamps;i++){fprintf(output,"IND%d IND%d %.6f\n", i+1, i+1, phens[i]);}
fclose(output);

return(0);
}	//end of main

///////////////////////////

