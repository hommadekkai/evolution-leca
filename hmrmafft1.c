#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
int main(void)
{
	FILE *fr, *fl;
	char rstr[122], empty[250]="", tempstr[21], onumbstr[21];
	char command[200], gname[20000][21];
	char fnamel[100]="/Users/khomma/work/ensembl/hmr/hmrseq1list";
	char fnamer[100]="/Users/khomma/work/ensembl/hmr/hmrseq1/ENSP00000354687";
        int  v, w, x, y, z, gmax=-1;
	double gw, gx, gy, gz;
	
	for(x=0; x<20000; x++)	strncpy(gname[x],empty,21);
	if((fl=fopen(fnamel,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamel);      exit(1);        }
	strncpy(rstr,empty,122);
	while((fgets(rstr,121,fl))!=NULL){
		++gmax;
		if(gmax>=20000){	printf("The number of ortholog sequence files exceeded 20,000!\n");	exit(1);	}
		w = strlen(rstr);
		if(w>=21){	printf("The line length exceeded 21 in %s",rstr);	exit(1);	}
		strncpy(gname[gmax],rstr,w-1);
		strncpy(rstr,empty,122);
	}
	fclose(fl);
	for(w=0; w<=gmax; w++){
		if((w+1)%250==0) printf("Processing ortho gp. %d out of %d...\n",w+1,gmax+1);
		strncpy(onumbstr,empty,21);
		sprintf(onumbstr,"%d",w);
        	strncpy(command,empty,200);
        	strcpy(command,"/usr/local/bin/mafft --auto --quiet ");		/* The quite option is used to avoid progress report */
		strcat(command,"/Users/khomma/work/ensembl/hmr/hmrseq1/");
        	strcat(command,gname[w]);
        	strcat(command," > ");
		strcat(command,"/Users/khomma/work/ensembl/hmr/hmrmafft1/");
		strcat(command,gname[w]);
		system(command);
	}
	return 0;
}
