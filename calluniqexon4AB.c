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
	int  totex;             /* #### Total # of exons */
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
                else if(catcount==6)	vlist[nex].end = atoi(tempstr);		/* End nt */
                else if(catcount==7)	strcpy(vlist[nex].sense,tempstr);	/* Sense */
                else if(catcount==8)	vlist[nex].exnumb = atoi(tempstr);	/* Exon number */	
		else if(catcount==9)    vlist[nex].totex = atoi(tempstr);       /* #### Total # of exons */
                else if(catcount==10)	vlist[nex].resstart = atoi(tempstr);	/* #### The start res. of the exon */
                else if(catcount==11){	/* #### */
			vlist[nex].resend = atoi(tempstr);	/* The end res. of the exon */
			return;
		}
		
nonreg:         for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;          /* ### */
                offset = x + offset;
        }
        return;
}
void infoout(void){   /* Output the regions of internal exons */
	if(nex<0)	return;		/* #### */	
        for(x=0; x<=nex; x++){ 		/* Examine all registered exons #### */
		if(vlist[x].start<0 || vlist[x].end<0)  continue;       /* ##&& The exon has been eliminated */
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

		if(vlist[x].exnumb==1)				fputs("F\t",fw);	/* #### */
		else if(vlist[x].exnumb==vlist[0].totex)	fputs("T\t",fw);	/* #### */
		else 						fputs("I\t",fw);	/* ##&& */

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
void singleout(void){   /* Output the region of a signle-exon gene  */
	if(nex<0)	return;		/* #### */
        for(x=0; x<=nex; x++){	/* Examine all registered exons #### nex=0 */ 
		if(vlist[x].exnumb!=1)	continue;	/* #### */
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

		fputs("S\t",fw);
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
	char fnamer[122]="/Users/khomma/work/ensembl/cele/calluniqexon4";	/* ##&& List of nonredundant exons */
	char fnamew[122]="/Users/khomma/work/ensembl/cele/calluniqexon4AB";	/* ##&& List of int, term, and single exons */
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs(">Chromosome\tGene\tTranscript\tLen(nt)\tIntExStart\tIntExEnd\tSense\tClass\tStart(nt)\tEnd(nt)\n",fw); 	/* ##&& */
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
			if(vlist[0].totex>=2){		/* #### Genes with 2 or more exons. Output the terminal exons */	
				infoout();
				++gcount;
				if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
			}
			else if(vlist[0].totex==1){	/* #### Single gene */
				if(vlist[0].cdslen != vlist[0].resend) goto skip;	/* This is wrong and should be eliminated */
				singleout();
				++gcount;
				if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
			}
skip:			nex = 0;
			tread(rstr);
		}
nextr:		strncpy(rstr,empty,122);
	}
	if(vlist[0].totex>=2){ /* #### Process the last gene */
		infoout();
		++gcount;
		if(gcount%1000==0)	printf("Now processing gene %d...\n",gcount);
	}
	else if(vlist[0].totex==1){     /* #### Single gene */
		singleout();
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
