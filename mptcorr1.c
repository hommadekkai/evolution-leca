/* mptcorr1.c	This program outputs the correspondence between mouse protein and transcript IDs		210512
	basic.c -> mptcorr1.c: A new program [R2105 II].		  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
int main(void)
{
	FILE *fr, *fw;
	char *ps;
	char rstr[1001], empty[1010]="", tempstr[21], ptemp[21], ttemp[21];
	char fnamer[122]="/Users/khomma/work/ensembl/mouse/Mus_musculus.GRCm39.103.chr.gff3";
	char fnamew[122]="/Users/khomma/work/ensembl/mouse/mptcorr1";
	char ptcorr[100000][2][21];	/* 0: Protein ID (ENSMUSP), 1: Transcript ID (ENSMUST) */
        int w, x, y, z, ptmax=-1;
	double gw, gx, gy, gz;
	for(x=0; x<100000; x++){
		for(y=0; y<2; y++)	strncpy(ptcorr[x][y],empty,21);
	}
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
        strncpy(rstr,empty,1001);
        while((fgets(rstr,1000,fr))!=NULL){
		if(strstr(rstr,"ID=CDS:ENSMUSP")!=NULL && strstr(rstr,"Parent=transcript:ENSMUST")!=NULL){
			ps = strstr(rstr,"ENSMUSP");
			strncpy(ptemp,empty,21);
			for(x=0; x<20; x++){
				if(x>=7 && isdigit(*ps)==0)	break; 
				ptemp[x] = *ps;
				++ps;
			}
			ps = strstr(rstr,"ENSMUST");
			strncpy(ttemp,empty,21);
			for(x=0; x<20; x++){
				if(x>=7 && isdigit(*ps)==0)	break; 
				ttemp[x] = *ps;
				++ps;
			}
			if(ptmax>=0){	/* Check if the correspondence has been saved */
				for(x=0; x<=ptmax; x++){
					if(strcmp(ptcorr[x][0],ptemp)==0){	/* Already registered */
						if(strcmp(ptcorr[x][1],ttemp)==0)	goto nextr;	/* Transcript ID matched */
						else{
							printf("New correspondence for %s\n",ptemp);
							goto newcorr;
						}
					}
				}
			}
newcorr:		/* Register the new correspondence */
			++ptmax;
			if((ptmax+1)%2000==0)	printf("Identifying correspondence #%d...\n",ptmax+1);
			if(ptmax>=100000){	printf("The # of correspondence exceeded 100,000!\n");	exit(1);	}
			strcpy(ptcorr[ptmax][0],ptemp);
			strcpy(ptcorr[ptmax][1],ttemp);
		}
                
nextr:		strncpy(rstr,empty,1001);
        }
	/* Output the correspondences: */
	for(x=0; x<=ptmax; x++){
		fputs(ptcorr[x][0],fw);
		fputs("\t",fw);
		fputs(ptcorr[x][1],fw);
		fputs("\n",fw);
	}
	printf("%d protein-transcript correspondences were identified\n",ptmax+1);
        fclose(fr);        fclose(fw);
	return 0;
}
