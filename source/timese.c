/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Do multiplications for pseudo jackknives

///////////////////////////

//set neff, the assumed number of individuals
sum=0;for(j=0;j<data_length;j++){sum+=nss[j];}
neff=sum/data_length;

if(bitsize>data_length){bitsize=data_length;}

//allocate variables
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

anal_warn(data_length,num_blocks);
randnorms=malloc(sizeof(double)*num_samples_use*num_blocks);
datarands=malloc(sizeof(double)*bitsize*num_blocks);
Mrhos=malloc(sizeof(double)*data_length*num_blocks);

//generate randnorms
for(k=0;k<num_blocks;k++)
{
for(i=0;i<num_samples_use;i++){randnorms[i+k*num_samples_use]=rnorm_safe();}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

////////

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Sampling summaries for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Sampling summaries for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

//read data and standardize
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1, preds+bitstart);

//get t(data) * randnorms * root(1/num_samples_use*(1-subprop)/subprop)
alpha=pow(1.0/num_samples_use*(1-subprop)/subprop,.5);beta=0.0;
dgemm_("T", "N", &bitlength, &num_blocks, &num_samples_use, &alpha, data, &num_samples_use, randnorms, &num_samples_use, &beta, datarands, &bitlength);

//Mrhos is rhos plus datarands/root(nss)
for(k=0;k<num_blocks;k++)
{
for(j=bitstart;j<bitend;j++)
{Mrhos[(size_t)k*data_length+j]=rhos[j]+datarands[(j-bitstart)+k*bitlength]*pow(nss[j],-.5);}
}
}	//end of bit loop
printf("\n");

//count number of trivial predictors (these will have Mrhos zero)
count=0;for(j=0;j<data_length;j++){count+=(mults[j]!=-9999);}
if(count==0){printf("Error, all %d predictors trivial\n\n", data_length);exit(1);}
if(count>0){printf("Warning, %d predictors trivial\n\n", count);}

////////

//save summaries
sprintf(filename2,"%s.jacks.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

for(k=0;k<num_blocks;k++)
{
fwrite(Mrhos+(size_t)k*data_length, sizeof(double), data_length, output2);
}
fclose(output2);

//save predictors
sprintf(filename3,"%s.jacks.predictors",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(j=0;j<data_length;j++){fprintf(output3,"%s %.1f\n", preds[j], nss[j]);}
fclose(output3);

//save root
sprintf(filename4,"%s.jacks.root", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Summaries %s\nTraining_Proportion %.4f\nNum_Blocks %d\n", sumsfile, subprop, num_blocks);
fprintf(output4,"Num_Predictors %d\nNum_Predictors_Used %d\n", num_preds, data_length);
fclose(output4);

printf("Jackknife samplings saved with prefix %s.jacks\n\n", outfile);

free(data);
free(Mrhos);
free(randnorms);free(datarands);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

