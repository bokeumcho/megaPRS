//I previously made randoms when constructing kinships
//Here are the related scripts

//some declarations in declare.c
double *randnorms, *ksamps;

//code from kins.c
randnorms=malloc(sizeof(double)*bitsize*num_vects);
ksamps=malloc(sizeof(double)*num_samples_use*num_vects);
for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){ksamps[i+g*num_samples_use]=0;}
}

if(num_vects>0)	//deal with random samples
{
for(g=0;g<num_vects;g++)
{
for(j=0;j<bitlength;j++){randnorms[j+g*bitlength]=rand_safe();}
}
alpha=1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_vects, &bitlength, &alpha, data, &num_samples_use, randnorms, &bitlength, &beta, ksamps, &num_samples_use);
}

//code from filekins.c

//in write_kins)
fprintf(output3, "Num_Preds %d\nSum_Weights %.4f\nDenominator %.4f\nOff_Diag_Variance %e\nNum_Random_Samplings %d", count, value, denom, 2*sumsq/ns/(ns-1)-pow(2*sum/ns/(ns-1),2), num_vects);

if(num_vects>0)	//write samplings
{
sprintf(filename5,"%s.grm.randoms", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}

value=pow(denom,-.5);
for(i=0;i<ns;i++)
{
for(g=0;g<num_vects;g++){fprintf(output5, "%.6f ", ksamps[i+g*ns]*value);}
fprintf(output5, "\n");
}
fclose(output5);
}

////////

int read_randoms(char *kinstem, double *R, int ns, int num_vects, char **ids3, int type)
//type=0 - normal, type=1 - quiet
{
int i, g, count, count2, found;
int *indexer, *indexer2;

char readchar, **wantids;

char filename[500];
FILE *input;


if(type==0){printf("Reading random samplings for kinship matrix with stem %s\n",kinstem);}

//first get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count2!=ns){printf("Error find indexes read_randoms, please tell Doug\n\n");exit(1);}

//open randoms
sprintf(filename, "%s.grm.randoms", kinstem);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

found=0;
for(i=0;i<count;i++)
{
if(i==indexer[found])	//will be using
{
for(g=0;g<num_vects;g++)
{
if(fscanf(input, "%lf%c", R+indexer2[found]+g*ns, &readchar)!=2)
{printf("Error reading Element %d of Row %d of %s\n\n", g+1, i+1, filename);exit(1);}
}
//get to end of row
while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
found++;
if(found==ns){break;}
}
else	//skip row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}	//end of i loop
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);

return(0);
}

////////

//some checks in consistent.c

if((priora!=-9999||priorb!=-9999||priorc!=-9999)&&mode!=128)
{printf("Error, you can only use \"--priora\", \"--priorb\" or \"--priorc\" with \"--calc-genes-reml\"\n\n");exit(1);}

/////////

//from kinfuns.c

//allocate samplings and set to zero
ksamps=malloc(sizeof(double)*ns*num_vects);
for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){ksamps[i+g*ns]=0;}
}

for(k=0;k<num_kins;k++)	//read through, adding or subtracting
{
value=1;if(scales!=NULL){value=scales[k];}
if(type==2&&k>0){value=-value;}
if(value!=0)
{
sprintf(filename, "%s.grm.randoms", kinstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
for(i=0;i<ns;i++)
{
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row %d of %s\n\n", i+1, filename);exit(1);}
if(strcmp(readstring,ids1[i])!=0||strcmp(readstring2,ids2[i])!=0)
{printf("Error, the sample IDs on Row %d of %s (%s and %s) do not match those on Row %d of %s.grm.id (%s and %s)\n\n", i+1, filename, readstring, readstring2, i+1, filename, ids1[i], ids2[i]);exit(1);}

for(g=0;g<num_vects;g++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", g+3, i+1, filename);exit(1);}
ksamps[i+g*ns]+=value*atof(readstring);
}
}
fclose(input);
}
}


///////////////////////////

