/* calluniqexon3.c		This program removes redundancy in ALL exons				220222
	basic.c -> huniqexon1.c: A new program [R2011 I].
	 -> huniqexon2.c: Revised to output terminal and single exons instead of internal exons [R2011 I].
	 -> halluniqexon2.c: Revised to process ALL exons [R2104 I].
	 -> halluniqexon3.c: Revised to process all genuine exons [R2202 I].
	 -> calluniqexon3.c: Modified to process C. elegans instead of human data [R2202 III].		  */
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
int  w, x, y, z, offset, vmax=-1, limv=-1, hit=0, gcount=0;
double gw, gx, gy, gz;
struct variant{
        char chr[10];
        char gname[21];         /* #### Gene name */
        char tname[21];         /* #### Transcript name */
        int  cdslen;
        int  start;
        int  end;
        char sense[2];          /* + or - */
	int  exnumb;		/* Exon number */
	int  totex;		/* #### Total # of exons */
	char resstart[21];	/* ### The start res # of the exon, can be 'noncoding' */
	char resend[21];		/* ### The end res # of the exon, can be 'noncoding' */
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
                        strcpy(vlist[vmax].chr,tempstr);        /* Chromosome name */
                        if(chmax>=0){
                                for(x=0; x<=chmax; x++)         if(strcmp(chname[x],vlist[vmax].chr)==0)        goto nonreg;
                                        /* The name of the chromosome has been registered */
                        }
                        ++chmax;
                        if(chmax>=50){  printf("The # of chromosome names exceeded 50\n");      exit(1);        }
                        strcpy(chname[chmax],vlist[vmax].chr);
                }       /* ### */
                else if(catcount==2)	strcpy(vlist[vmax].gname,tempstr);	/* Gene name */
                else if(catcount==3)	strcpy(vlist[vmax].tname,tempstr);	/* Transcript name */
                else if(catcount==4) 	vlist[vmax].cdslen = atoi(tempstr);       /* CDS length */
                else if(catcount==5) 	vlist[vmax].start = atoi(tempstr);       /* Start nt */
                else if(catcount==6)	vlist[vmax].end = atoi(tempstr);	/* End nt */
                else if(catcount==7)	strcpy(vlist[vmax].sense,tempstr);	/* Sense */
                else if(catcount==8)	vlist[vmax].exnumb = atoi(tempstr);	/* Exon number */	
                else if(catcount==9)	vlist[vmax].totex = atoi(tempstr);	/* #### Total # of exons */	
                else if(catcount==10)	strcpy(vlist[vmax].resstart,tempstr);	/* #### The start res. of the exon */
                else if(catcount==11){
			strcpy(vlist[vmax].resend,tempstr);	/* #### The end res. of the exon */
			return;
		}
		
nonreg:         for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;          /* ### */
                offset = x + offset;
        }
        return;
}
void tunique(){ /* ### To leave the unique exons only */
        if(vmax==0)     return;         /* There is no need to remove redundandancy in this case */
        for(x=0; x<=vmax; x++){
                if(vlist[x].start<0 || vlist[x].end<0)          continue;       /* Eliminated due to redundancy */
                for(y=0; y<x; y++){     /* Compare with the previous (OKed) exons */
                        if(vlist[y].start<0 || vlist[y].end<0)          continue;       /* Eliminated due to redundancy */
                        else if(vlist[x].start==vlist[y].start && vlist[y].end==vlist[x].end){	/* #### */
                                /* Identical exons */
                                if(vlist[x].cdslen<vlist[y].cdslen){       /* Eliminate exon x, leave exon y */
                                        vlist[x].start = -1;            vlist[x].end = -1;
                                }
                                else{   /* Leave exon x, eliminate exon y */
                                        vlist[y].start = -1;            vlist[y].end = -1;
                                }
                        }
                }
        }
        return;
}
void infoout( void ){   /* Output the regions of all exons */
        for(x=0; x<=vmax; x++){        /* After elimination of redundant exons */
		if(vlist[x].start<0 || vlist[x].end<0)  continue;       /* ### The exon has been eliminated */
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

                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].exnumb);       /* Exon # */
                fputs(tempstr,fw);      fputs("\t",fw);

		/* #### */
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",vlist[x].totex);       /* Total # of exons */
                fputs(tempstr,fw);      fputs("\t",fw);
		/* #### */

                fputs(vlist[x].resstart,fw);    fputs("\t",fw);		/* #### The start res. of the exon, can be 'noncoding' */
                fputs(vlist[x].resend,fw);      fputs("\n",fw);		/* #### The end res. of the exon, can be 'noncoding' */
                x = x;
        }
        return;
}
int main(void)
{
	char fnamer[122]="/Users/khomma/work/ensembl/cele/callexon3";		/* List of all genuine exons */
	char fnamew[122]="/Users/khomma/work/ensembl/cele/calluniqexon3";

        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs(">Chromosome\tGene\tTranscript\tLength\tIntExStart\tIntExEnd\tSense\tExon #\tTotEx#\tRes start\tRes end\n",fw);
        while((fgets(rstr,121,fr))!=NULL){
		if(rstr[0]=='>' || strstr(rstr,"ERROR")!=NULL)	goto nextr; 
		else if(strstr(rstr,"Noncoding")!=NULL)		goto nextr;	/* Neglecting noncoding exons */
		else if(vmax<0){                /* ### The first gene */
			++gcount;
			if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
			++vmax;
			if(vmax>=10000){	printf("The # of variants exceeded 10,000 in %s",rstr);	exit(1);	}
			if(vmax>limv)	limv = vmax;
			tread(rstr);
		}	
		else if(strstr(rstr,vlist[0].gname)!=NULL){	/* A transcript of the same gene */
			++vmax;
			if(vmax>=10000){	printf("The # of variants exceeded 10,000 in %s",rstr);	exit(1);	}
			if(vmax>limv)	limv = vmax;
			tread(rstr);
		}
		else{		/* A new gene is encountered. Process the last gene first */
			tunique();	/* Elminate redundancies, if any */
			infoout();
			++gcount;
			if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
			vmax = 0;
			tread(rstr);
		}
nextr:		strncpy(rstr,empty,122);
        }
	if(vmax>=0){	/* Process the last gene */
		tunique();	/* Elminate redundancies, if any */
		infoout();
	}
        fclose(fr);        fclose(fw);
	printf("%d nonredundant exons in %d genes have been output, with the max # of variants = %d\n",hit,gcount,limv+1);
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
