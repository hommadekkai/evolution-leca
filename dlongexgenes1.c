#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fl;
char fnamer[122]="/Users/khomma/work/ensembl/dmel/dalluniqexon3A";
char fnamew[122]="/Users/khomma/work/ensembl/dmel/dlongexgenes1";
char rstr[501], empty[510]="", tempstr[21];
char longlist[100000][21];
char binstr[10][20]={"All","-120","-240","-360","-480","-600","-720","-840","-960","961bp-"};   /* #### */
int q, w, x, y, z, offset, longmax=-1, counter=0, hit=0, mlenbin;
int lendist[10], mlendist[10];
double gw, gx, gy, gz;
struct variant{
        char chr[21];
        char gname[21];         /* Gene name */
        char tname[21];         /* Transcript name */
        int  cdslen;
        int  start;
        int  end;
        char sense[2];          /* + or - */
        int  resstart;          /* The start res # of the exon */
        int  resend;            /* The end res # of the exon */
};
struct variant vlist;
void tread(char rstr[501]){    /* Parse the line to get the start and end of the terminal exon */
        int catcount=0;         /* Category counter */
        offset = 0;
        while(1){
                ++catcount;
                strncpy(tempstr,empty,21);
                for(x=0; x<20; x++){
                        if(isspace(rstr[x+offset])!=0)          break;
                        tempstr[x] = rstr[x+offset];
                }
                offset = x + offset;
                if(catcount==1)		strcpy(vlist.chr,tempstr);        /* Chromosome name */
                else if(catcount==2)    strcpy(vlist.gname,tempstr);       /* Gene name */
                else if(catcount==3)    strcpy(vlist.tname,tempstr);       /* Transcript name */
                else if(catcount==4)    vlist.cdslen = atoi(tempstr);       /* CDS length */
                else if(catcount==5)    vlist.start = atoi(tempstr);       /* Start nt */
                else if(catcount==6)    vlist.end = atoi(tempstr); 	   /* End nt */
                else if(catcount==7)    strcpy(vlist.sense,tempstr);       /* Sense */
                else if(catcount==8)    vlist.resstart = atoi(tempstr);    /* #### The start res. of the exon */
                else if(catcount==9){ 
                        vlist.resend = atoi(tempstr);      /* The end res. of the exon */
                        return;
                }

nonreg:         for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;          /* ### */
                offset = x + offset;
        }
        return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(x=0; x<100000; x++)	strncpy(longlist[x],empty,21);
	for(x=0; x<10; x++){
		lendist[x] = 0;
		mlendist[x] = 0;	/* The number of GENES with max exon length falling in each range */
	}
	for(mlenbin=9; mlenbin>=1; mlenbin--){	/* Max length bins */
		printf("Now processing %s...\n",binstr[mlenbin]);
		fseek(fr,0L,SEEK_SET);
        	strncpy(rstr,empty,501);
        	while((fgets(rstr,500,fr))!=NULL){
 	               if(rstr[0]=='>' || strstr(rstr,"ERROR")!=NULL || strstr(rstr,"Noncoding")!=NULL)        goto nextr;
			tread(rstr);
			w = abs(vlist.end - vlist.start) + 1;
			q = (w-1)/120 + 1;         /* Simple exon length bin */
                	if(q>9)     q = 9;

			if(mlenbin==9){
				++lendist[q];
				++lendist[0];
			}
			if(q==mlenbin){	/* A long internal exon */
				++hit;
				if(longmax>=0){		/* Check if the gene ID has already been registered */
					for(z=0; z<=longmax; z++)	if(strcmp(longlist[z],vlist.gname)==0)		goto nextr;
					/* Already registered, i.e., the gene has multiple long exons */
				}
				++longmax;
				if(longmax>=100000){ printf("The # of genes with long internal exons exceeded 100,000\n"); exit(1);	}
				strcpy(longlist[longmax],vlist.gname);
				++mlendist[mlenbin];
				++mlendist[0];
			}
nextr:         	 	strncpy(rstr,empty,501);
        	}
		if(mlenbin==7){		/* max exon length >= 721 bp */
			printf("%d genes had a total of %d long internal exons\n",longmax+1, hit);
			for(z=0; z<=longmax; z++){
				fputs(longlist[z],fw);		fputs("\n",fw);
			}
		}
		mlenbin = mlenbin;
	}
	printf("A total of %d genes were classified\n",mlendist[0]);
	printf("Exon len\tExon #s\n");
	for(q=0; q<10; q++){
		printf("%s\t%d\n",binstr[q],lendist[q]);
	}
	printf("Max exon len\tGene #s\n");
	for(q=0; q<10; q++){
		printf("%s\t%d\n",binstr[q],mlendist[q]);
	}
        fclose(fr);        fclose(fw);
	return 0;
}
