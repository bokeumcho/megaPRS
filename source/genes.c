/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Gene-based analyses

///////////////////////////

//count number of gene in this partition and set genemax to size of longest gene
count=0;genemax=0;
for(g=0;g<num_genes;g++)
{
if(gparts[g]==partition)
{
count++;
if(gends[g]-gstarts[g]>genemax){genemax=gends[g]-gstarts[g];}
}
}
printf("The longest gene/chunk contains %d predictors\n\n", genemax);

//allocate variables
data_warn(num_samples_use, genemax);
data=malloc(sizeof(double)*num_samples_use*genemax);

if(mode==137)
{
anal_warn(num_samples_use, num_samples_use);
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
}
else
{
if(num_kins==0){anal_warn(num_samples_use+2*genemax, genemax);}
else{anal_warn(5*num_samples_use+3*genemax, genemax);}

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
X=malloc(sizeof(double)*num_samples_use*genemax);

if(strcmp(sumsfile,"blank")!=0)
{
Xnss=malloc(sizeof(double)*genemax);
Xrhos=malloc(sizeof(double)*genemax);
Xsqs=malloc(sizeof(double)*genemax);
}

retain=malloc(sizeof(int)*genemax);
stats=malloc(sizeof(double)*8);
if(gene_perms>0){permlike=malloc(sizeof(double)*gene_perms);}

if(num_kins==0)
{
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTY=malloc(sizeof(double)*num_fixed);
ZTZZTY=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);

ZTX=malloc(sizeof(double)*num_fixed*genemax);
ZTZZTX=malloc(sizeof(double)*num_fixed*genemax);
XTCX=malloc(sizeof(double)*genemax*genemax);
E2=malloc(sizeof(double)*genemax);
}
else	//most allocations done within remladv.c
{vstarts=malloc(sizeof(double)*3);}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

if(mode==138)	//fill Y and Z, and if using, get factor - then stuff for nk=0 and nk=1
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

if(prev!=-9999){factor=get_factor(Y, num_samples_use, prev, ascer, NULL);}

if(num_kins==0)	//deal with covariates
{
//get ZTZ and inverse ZTZ
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &num_samples_use, &alpha, Z, &num_samples_use, Z, &num_samples_use, &beta, ZTZ, &num_fixed);
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

//get ZTY and ZTZZTY = inverse ZTZ ZTY
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_fixed, &alpha, Z, &num_samples_use, Y, &one, &beta, ZTY, &one);
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, ZTZZTY, &one);

//Yadj contains residuals from regressing Y on Z
for(i=0;i<num_samples_use;i++){Yadj[i]=Y[i];}
alpha=-1.0;beta=1.0;
dgemv_("N", &num_samples_use, &num_fixed, &alpha, Z, &num_samples_use, ZTZZTY, &one, &beta, Yadj, &one);
}
else	//get starting heritabilities, varnull and likenull (saved in vstarts)
{
printf("Solving Null Model\n");
(void)adv_reml(num_samples_use, num_fixed, 0, Y, Z, U, E, kintraces, NULL, -9999, -9999, NULL, vstarts, tol, maxiter);
printf("\n");
}
}

//deal with progress and on-the-fly files
sprintf(filename,"%sprogress.%d", folder, partition);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

if(mode==138)
{
sprintf(filename2,"%sremls.%d", folder, partition);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Gene_Number Gene_Name Num_Predictors Heritability SD Null_Likelihood Alt_Likelihood LRT_Stat LRT_P_Raw Her_Liability SD_Liability");
for(p=0;p<gene_perms;p++){fprintf(output2, " Perm_%d", p+1);}
fprintf(output2, "\n");
}

////////

count2=0;
for(g=0;g<num_genes;g++)
{
total=gends[g]-gstarts[g];

if(gparts[g]==partition)
{
if(mode==137)
{
if(count2%200==0)
{
printf("Calculating kinships for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating kinships for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
}
}
else
{
if((num_kins==0&&count2%200==0)||(num_kins==1&&count2%20==0))
{
printf("Performing REML for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Performing REML for Gene/Chunk %d of %d in Partition %d\n", count2+1, count, partition);
fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
}
}

shuffle=0;	//this is how many predictors we already have data for (can be more than total)
for(j=gstarts[g];j<end;j++)
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(j-start)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, gstarts[g]+shuffle, gends[g], keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data+(size_t)shuffle*num_samples_use, centres+gstarts[g]+shuffle, mults+gstarts[g]+shuffle, sqdevs+gstarts[g]+shuffle, num_samples_use, total-shuffle, missingvalue, power, 0, hwestand, weights+gstarts[g]+shuffle, 1, preds+gstarts[g]+shuffle);

//ready to do work
if(mode==137)	//calc-kins-genes
{
//find number of non-trivial predictors
count3=0;for(j=gstarts[g];j<gends[g];j++){count3+=(mults[j]!=-9999&&weights[j]>0);}

if(count3==0)	//gene trivial
{
printf("Warning, all %d predictors in Gene/Chunk %d are trivial", gends[g]-gstarts[g], g+1);
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n");
}
else	//non-trivial
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &total, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);

sprintf(outfile,"%sgeneships.%d", folder, g+1);
write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, preds+gstarts[g], keeppreds_use+gstarts[g], centres+gstarts[g], mults+gstarts[g], weights+gstarts[g], al1+gstarts[g], al2+gstarts[g], total, datafile, power, kingz, kinraw, 0);
}
}
else	//calc-genes-reml - note that will have summary statistics for all predictors
{
//find number of non-trivial predictors (indexed by retain)
count3=prune_gene(retain, gprune, data, num_samples_use, total, mults+gstarts[g], weights+gstarts[g], X, XTCX);

if(count3==0)	//gene trivial
{
printf("Warning, all %d predictors in Gene/Chunk %d are trivial", gends[g]-gstarts[g], g+1);
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n");

fprintf(output2, "%d %s 0 NA NA NA NA NA NA NA NA", g+1, gnames[g]);
for(p=0;p<gene_perms;p++){fprintf(output2, " NA");}
fprintf(output2, "\n");
}
else	//non-trivial
{
//fill X and, if using, summary statistics
if(strcmp(sumsfile,"blank")==0)
{sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, NULL, retain, data, NULL, NULL);}
else
{sumsq=fill_gene(X, Xnss, Xrhos, num_samples_use, total, NULL, retain, data, nss+gstarts[g], rhos+gstarts[g]);}

if(num_kins==0)	//compute and decompose XTCX
{
//compute ZTX and ZTZZTX = inverse ZTZ ZTX
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &count3, &num_samples_use, &alpha, Z, &num_samples_use, X, &num_samples_use, &beta, ZTX, &num_fixed);
dgemm_("N","N", &num_fixed, &count3, &num_fixed, &alpha, ZTZ, &num_fixed, ZTX, &num_fixed, &beta, ZTZZTX, &num_fixed);

//compute XTCX = XTX - XTZ ZTZ ZTX
alpha=1.0;beta=0.0;
dgemm_("T", "N", &count3, &count3, &num_samples_use, &alpha, X, &num_samples_use, X, &num_samples_use, &beta, XTCX, &count3);
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &count3, &count3, &num_fixed, &alpha, ZTX, &num_fixed, ZTZZTX, &num_fixed, &beta, XTCX, &count3);

if(strcmp(sumsfile,"blank")!=0)	//extract diagonals of XTCX (will not have covariates, so XTCX=XTX)
{
for(j=0;j<count3;j++){Xsqs[j]=XTCX[j+j*count3];}
}

if(shrink<1)	//deflate off diagonal terms
{
for(j=0;j<count3;j++)
{
for(j2=j+1;j2<count3;j2++){XTCX[j+j2*count3]*=shrink;XTCX[j2+j*count3]*=shrink;}
}
}
if(shrink>1)	//inflate diagonal terms
{
for(j=0;j<count3;j++){XTCX[j+j*count3]*=shrink;}
}

//decompose XTCX
lwork=-1;
dsyev_("V", "U", &count3, XTCX, &count3, E2, &wkopt, &lwork, &info);
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &count3, XTCX, &count3, E2, work, &lwork, &info);
free(work);

/*
if(strip>0)	//set end proportion of E2 to zero
{
sum=0;for(j=0;j<Xtotal;j++){sum+=E2[j];}
sum2=0;
for(j=Xtotal-1;j>=0;j--)
{
sum2+=E2[j];
if(sum2/sum>1-strip){E2[j]=0;}
}
}
*/
}

//analyse real data
if(strcmp(sumsfile,"blank")==0)	//have phenotypes
{
if(num_kins==0)
{(void)gene_reml(num_samples_use, num_fixed, Yadj, detZTZ, X, count3, XTCX, E2, sumsq, NULL, NULL, NULL, stats, tol, maxiter);}
else
{(void)adv_reml(num_samples_use, num_fixed, 1, Y, Z, U, E, kintraces, X, count3, sumsq, stats, vstarts, tol, maxiter);}
}
else	//using summary statistics
{(void)gene_reml(num_samples_use, -9999, Yadj, detZTZ, X, count3, XTCX, E2, sumsq, Xnss, Xrhos, Xsqs, stats, tol, maxiter);}

for(p=0;p<gene_perms;p++)
{
//when permuting, always use individual-level data
permute_int(order,num_samples_use);
sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, order, retain, data, NULL, NULL);

if(num_kins==0&&num_fixed>1)	//must re-calculate and decompose XTCX
{
//compute ZTX and ZTZZTX
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &count3, &num_samples_use, &alpha, Z, &num_samples_use, X, &num_samples_use, &beta, ZTX, &num_fixed);
dgemm_("N","N", &num_fixed, &count3, &num_fixed, &alpha, ZTZ, &num_fixed, ZTX, &num_fixed, &beta, ZTZZTX, &num_fixed);

//compute XTCX - not using summaries, so no need to extract diagonals
alpha=1.0;beta=0.0;
dgemm_("T", "N", &count3, &count3, &num_samples_use, &alpha, X, &num_samples_use, X, &num_samples_use, &beta, XTCX, &count3);
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &count3, &count3, &num_fixed, &alpha, ZTX, &num_fixed, ZTZZTX, &num_fixed, &beta, XTCX, &count3);

if(shrink<1)	//deflate off diagonal terms
{
for(j=0;j<count3;j++)
{
for(j2=j+1;j2<count3;j2++){XTCX[j+j2*count3]*=shrink;XTCX[j2+j*count3]*=shrink;}
}
}
if(shrink>1)	//inflate diagonal terms
{
for(j=0;j<count3;j++){XTCX[j+j*count3]*=shrink;}
}

//decompose XTCX
lwork=-1;
dsyev_("V", "U", &count3, XTCX, &count3, E2, &wkopt, &lwork, &info);
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &count3, XTCX, &count3, E2, work, &lwork, &info);
free(work);
}	//end of recalculating and decomposing XTCX

if(num_kins==0)
{permlike[p]=gene_reml(num_samples_use, num_fixed, Yadj, detZTZ, X, count3, XTCX, E2, sumsq, NULL, NULL, NULL, NULL, tol, maxiter);}
else
{permlike[p]=adv_reml(num_samples_use, num_fixed, 1, Y, Z, U, E, kintraces, X, count3, sumsq, NULL, vstarts, tol, maxiter);}
}	//end of p loop

//print results
fprintf(output2, "%d %s %d ", g+1, gnames[g], count3);
if(stats[1]!=-9999){fprintf(output2, "%.6f %.6f ", stats[0], stats[1]);}
else{fprintf(output2, "%.6f NA ", stats[0]);}
fprintf(output2, "%.4f %.4f %.4f %.4e ", stats[2], stats[3],  stats[4], stats[5]);
if(prev!=-9999)
{
if(stats[1]!=-9999){fprintf(output2,"%.6f %.6f ", stats[0]*factor, stats[1]*factor);}
else{fprintf(output2, "%.6f NA", stats[0]*factor);}
}
else{fprintf(output2, "NA NA");}
for(p=0;p<gene_perms;p++){fprintf(output2, " %.4f", permlike[p]);}
fprintf(output2, "\n");
}	//end of not trivial
}	//end of mode=138

count2++;
start=gstarts[g];if(gends[g]>end){end=gends[g];}
}}	//end of in right partition and g loop
printf("\n");

fclose(output);
if(mode==138){fclose(output2);}

if(mode==137){printf("Kinship matrices saved with stem %sgeneships\n\n", folder);}
else{printf("REML estimates saved in %sremls.%d\n\n", folder, partition);}

free(data);
if(mode==137){free(kins);}
else
{
free(order);free(Y);free(Z);free(X);
if(strcmp(sumsfile,"blank")!=0){free(Xnss);free(Xrhos);}
free(retain);free(stats);
if(gene_perms>0){free(permlike);}
if(num_kins==0)
{free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTY);free(ZTZZTY);free(Yadj);
free(ZTX);free(ZTZZTX);free(XTCX);free(E2);}
if(num_kins==1){free(vstarts);}
}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

