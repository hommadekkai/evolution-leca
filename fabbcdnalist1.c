#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fe;
char estr[591], rstr[5001], empty[5010]="", tempstr[21];
char fnamer[122]="/Users/khomma/work/ensembl/spomb/fcdnaorflist1";
char fnamew[122]="/Users/khomma/work/ensembl/spomb/fabbcdnalist1";
char fnamee[100]="/Users/khomma/work/ensembl/spomb/falluniqexon4AB";
char trlist[200000][21];
int w, x, y, z, offset, trmax=-1, counter=0, okflag=0, hitcount=0, hitflag=0;
double gw, gx, gy, gz;
struct variant{
        char chr[21];
        char gname[21];         /* Gene name */
        char tname[21];         /* Transcript name */
        int  cdslen;
        int  start;
        int  end;
        char sense[2];          /* + or - */
        char class[2];          /* F, T, or S (used for halluniqexon3B only) */
        int  resstart;          /* The start res # of the exon */
        int  resend;            /* The end res # of the exon */
};
struct variant vlist;
void sread(char estr[501]){    /* Parse the line to get the start and end the exon */
        int catcount=0;         /* Category counter */
        offset = 0;
        while(1){
                ++catcount;
                strncpy(tempstr,empty,21);
                for(x=0; x<20; x++){
                        if(isspace(estr[x+offset])!=0)          break;
                        tempstr[x] = estr[x+offset];
                }
                offset = x + offset;
                if(catcount==1)         strcpy(vlist.chr,tempstr);              /* Chromosome name */
                else if(catcount==2)    strcpy(vlist.gname,tempstr);            /* Gene name */
                else if(catcount==3)    strcpy(vlist.tname,tempstr);            /* Transcript name */
                else if(catcount==4)    vlist.cdslen = atoi(tempstr);           /* ##&& CDS length (nt) */
                else if(catcount==5)    vlist.start = atoi(tempstr);            /* Start genome address */
                else if(catcount==6)    vlist.end = atoi(tempstr);              /* End genome address */
                else if(catcount==7)    strcpy(vlist.sense,tempstr);            /* Sense */
                else if(catcount==8)    strcpy(vlist.class,tempstr);            /* Class: I, F, T, or S */
                else if(catcount==9)    vlist.resstart = atoi(tempstr);         /* ##&& The start CDS res. (nt) of the exon */
                else if(catcount==10){
                        vlist.resend = atoi(tempstr);                           /* ##&& The end CDS res. (nt) of the exon */
                        return;
                }
                for(x=0; x<50; x++)     if(isspace(estr[x+offset])==0)  break;
                offset = x + offset;
        }
        return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
        if((fe=fopen(fnamee,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamee);      exit(1);        }
	for(x=0; x<200000; x++)	strncpy(trlist[x],empty,21);
	/* Read in the list of uniq exons and save the transcript IDs */
       	strncpy(estr,empty,501);
        while((fgets(estr,500,fe))!=NULL){
		if(estr[0]!='>'){
			++counter;
			if(counter%10000==0)	printf("Now processing exon list line %d\n",counter);
			sread(estr);
			if(trmax>=0){	/* Check if the transcript ID has been registered */
				for(x=0; x<=trmax; x++)	if(strcmp(trlist[x],vlist.tname)==0)	goto nexte;	/* Registered */
			}
			++trmax;
			if(trmax>=200000){	printf("The # of unique transcript IDs exceeded 200,000\n");	exit(1);	}
			strcpy(trlist[trmax],vlist.tname);
		}
nexte:		strncpy(estr,empty,501);
        }
	printf("%d unique transcript IDs were found\n",trmax+1);
       	strncpy(rstr,empty,5001);	counter=0;
        while((fgets(rstr,5000,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%5000==0)	printf("Now processing cDNA %d...\n",counter);
			hitflag = 0;
			for(x=0; x<=trmax; x++){
				if(strstr(rstr,trlist[x])!=NULL){
					hitflag = 1;	/* OK */
					++hitcount;
					fputs(rstr,fw);
				}
			}
		}
		else if(hitflag==1)	fputs(rstr,fw);		/* Seq. line */

                strncpy(rstr,empty,5001);
        }
	printf("%d out of %d cDNAs were selected and saved\n",hitcount,counter);
        fclose(fr);        fclose(fw);
        fclose(fe);
	return 0;
}
