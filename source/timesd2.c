/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Do multiplications for pseudo summary statistics (when using pre-computed noise)

///////////////////////////

//allocate variables
datarands=malloc(sizeof(double)*data_length);
subrhos=malloc(sizeof(double)*data_length);
restrhos=malloc(sizeof(double)*data_length);

//read datarands
sprintf(filename,"%s.noise.values", noisestem);
read_values(filename,datarands,data_length,keeppreds_use,1,0,0);

//subrhos is rhos plus datarands*root((1-subprop)/subprop/nss)
value=pow((1-subprop)/subprop,.5);
for(j=0;j<data_length;j++){subrhos[j]=rhos[j]+datarands[j]*value*pow(nss[j],-.5);}

//restrhos is complement
for(j=0;j<data_length;j++){restrhos[j]=(rhos[j]-subprop*subrhos[j])/(1-subprop);}

//save summaries
sprintf(filename2,"%s.train.summaries",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor A1 A2 Direction Stat n\n");

sprintf(filename3,"%s.test.summaries",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Direction Stat n\n");

for(j=0;j<data_length;j++)
{
fprintf(output2, "%s %c %c %d %.4f %.1f\n", preds[j], al1[j], al2[j], (subrhos[j]>=0)-(subrhos[j]<0), subprop*nss[j]*pow(subrhos[j],2)/(1-pow(subrhos[j],2)), subprop*nss[j]);
fprintf(output3, "%s %c %c %d %.4f %.1f\n", preds[j], al1[j], al2[j], (restrhos[j]>=0)-(restrhos[j]<0), (1-subprop)*nss[j]*pow(restrhos[j],2)/(1-pow(restrhos[j],2)), (1-subprop)*nss[j]);
}

fclose(output2);
fclose(output3);

printf("New summary statistics saved in %s and %s\n\n", filename2, filename3);

free(datarands);free(subrhos);free(restrhos);

///////////////////////////

