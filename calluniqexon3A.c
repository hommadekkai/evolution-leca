/* calluniqexon3A.c		This program outputs only genuine internal exons of nonredundnat genes		220202
	basic.c -> huniqexon1.c: A new program [R2011 I].
	 -> huniqexon2.c: Revised to output terminal and single exons instead of internal exons [R2011 I].
	 -> halluniqexon2.c: Revised to process ALL exons [R2104 I].
	 -> halluniqexon2A.c: Revised to process halluniqexon2 to ouput internal exons only [R2104 I].
	 -> halluniqexon3A.c: Revised to process data of exons by the conventional def. (####: changed) [R2202 I].
	 -> calluniqexon3A.c: Modified to process C. eleganse instead of human data [R2202 III].	  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw;
char rstr[122], tempstr[21], empty[500]="";
char gname[51], tname[1000][51];
char chname[50][21];    /* I, II, ..., X, Y, MtDNA */           /* ### */
int  chmax=-1;          /* ### */
int  w, x, y, z, offset, nex=-1, limv=-1, hit=0, gcount=0;	/* #### exmax -> nex */
double gw, gx, gy, gz;
struct variant{
        char chr[10];
        char gname[51];         /* Gene name */
        char tname[51];         /* Transcript name */
        int  cdslen;
        int  start;
        int  end;
        char sense[2];          /* + or - */
	int  exnumb;		/* Exon number */
	int  totex;		/* #### Total # of exons */
	int  resstart;		/* The start res # of the exon */
	int  resend;		/* The end res # of the exon */
};
struct variant vlist[10000];
void tread(char rstr[122]){    /* Parse the line to get the start and end of the terminal exon */
        int catcount=0;         /* Category counter */
        offset = 0;
        while(1){
                ++catcount;
                strncpy(tempstr,empty,21);
                for(x=0; x<50; x++){
                        if(isspace(rstr[x+offset])!=0)          break;
                        tempstr[x] = rstr[x+offset];
                }
                offset = x + offset;
                if(catcount==1){        /* ### */
                        strcpy(vlist[nex].chr,tempstr);        /* Chromosome name */
                        if(chmax>=0){
                                for(x=0; x<=chmax; x++)         if(strcmp(chname[x],vlist[nex].chr)==0)        goto nonreg;
                                        /* The name of the chromosome has been registered */
                        }
                        ++chmax;
                        if(chmax>=50){  printf("The # of chromosome names exceeded 50\n");      exit(1);        }
                        strcpy(chname[chmax],vlist[nex].chr);
                }       /* ### */
                else if(catcount==2)	strcpy(vlist[nex].gname,tempstr);	/* Gene name */
                else if(catcount==3)	strcpy(vlist[nex].tname,tempstr);	/* Transcript name */
                else if(catcount==4) 	vlist[nex].cdslen = atoi(tempstr);       /* CDS length */
                else if(catcount==5) 	vlist[nex].start = atoi(tempstr);       /* Start nt */
                else if(catcount==6)	vlist[nex].end = atoi(tempstr);	/* End nt */
                else if(catcount==7)	strcpy(vlist[nex].sense,tempstr);	/* Sense */
                else if(catcount==8)	vlist[nex].exnumb = atoi(tempstr);	/* Exon number */	
                else if(catcount==9)	vlist[nex].totex = atoi(tempstr);	/* #### Total # of exons */
                else if(catcount==10)	vlist[nex].resstart = atoi(tempstr);	/* #### The start res. of the exon */
                else if(catcount==11){			/* #### */
			vlist[nex].resend = atoi(tempstr);	/* The end res. of the exon */
			return;
		}
		
nonreg:         for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;          /* ### */
                offset = x + offset;
        }
        return;
}
void infoout( void ){   /* Output the regions of all exons */
        for(x=0; x<=nex; x++){        /* Examine all registered exons #### */
		if(vlist[x].exnumb==1 || vlist[x].exnumb==vlist[0].totex)	continue;	/* Skip the first and last exons #### */
		else if(vlist[x].start<0 || vlist[x].end<0)  continue;       /* ### The exon has been eliminated */
                ++hit;
                fputs(vlist[x].chr,fw);    fputs("\t",fw);         /* Chromosome name */
		fputs(vlist[x].gname,fw);	fputs("\t",fw);		/* Gene name */
		fputs(vlist[x].tname,fw);	fputs("\t",fw);		/* Transcript name */
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].cdslen);
                fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].start);       /* The start of the terminal or single exon */
                fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].end);       /* The end of the terminal or single exon */
                fputs(tempstr,fw);      fputs("\t",fw);
		fputs(vlist[x].sense,fw);	fputs("\t",fw);
			/* ### No need to output the exon # or the total # of exons */
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].resstart); 	/* The start res. of the exon */
                fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].resend); 		/* The end res. of the exon */
                fputs(tempstr,fw);      fputs("\n",fw);
                x = x;
        }
        return;
}
int main(void)
{
	char fnamer[122]="/Users/khomma/work/ensembl/cele/calluniqexon3";		/* List of nonredundant exons */
	char fnamew[122]="/Users/khomma/work/ensembl/cele/calluniqexon3A";		/* List of internal exons */

        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs(">Chromosome\tGene\tTranscript\tLength\tIntExStart\tIntExEnd\tSense\tRes start\tRes end\n",fw);
        while((fgets(rstr,121,fr))!=NULL){
		if(rstr[0]=='>' || strstr(rstr,"ERROR")!=NULL || strstr(rstr,"Noncoding")!=NULL)	goto nextr;
		else if(nex<0){                /* ### The first gene */
			++nex;
			if(nex>=10000){	printf("The # of variants exceeded 10,000 in %s",rstr);	exit(1);	}
			if(nex>limv)	limv = nex;
			tread(rstr);
		}	
		else if(strstr(rstr,vlist[0].tname)!=NULL){	/* The same transcript */
			++nex;
			if(nex>=10000){	printf("The # of variants exceeded 10,000 in %s",rstr);	exit(1);	}
			if(nex>limv)	limv = nex;
			tread(rstr);
		}
		else{		/* A new transcript (possibly of a new gene) is encountered. Process the last gene first */
			if(vlist[0].totex>=3){	/* #### Only genes with >=3 exons are output. nex is not the total # of exons */
				infoout();	/* Only genes with >=3 exons are output */
				++gcount;
				if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
			}
			nex = 0;
			tread(rstr);
		}
nextr:		strncpy(rstr,empty,122);
        }
	if(vlist[0].totex>=3){ /* ### Process the last gene */
		infoout();
		++gcount;
		if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
	}
        fclose(fr);        fclose(fw);
	printf("%d nonredundant exons in %d genes have been output, with the max # of exons = %d\n",hit,gcount,limv+1);
        /* ### */
        printf("The list of chromosome names: ");
        for(x=0; x<=chmax; x++){
                printf("%s",chname[x]);
                if(x<chmax)     printf(", ");
                else            printf("\n");
        }
        /* ### */
	return 0;
}
