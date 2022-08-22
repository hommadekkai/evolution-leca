/* mcdnaorflist1.c:	This program finds and outputs the longest ORF of each cDNA seq.		220218
	basic.c -> acdnaorflist1.c: A new program [R2203 II].
	-> mcdnaorflist1.c: Modified to process mouse instead of A thaliana data [R2203 VII].		  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw;
char rstr[1001], header[1001], empty[1010]="", tempstr[21];
char fnamer[122]="/Users/khomma/work/ensembl/mouse/Mus_musculus.GRCm39.cds.all.fa";
char fnamew[122]="/Users/khomma/work/ensembl/mouse/mcdnaorflist1";
char seq[200000], cseq[200000];
int  w, x, y, z, offset, numbseq=-1, counter=0, okcount=0;
int  clen, cmax;
double gw, gx, gy, gz;
void orffind()
{
	char *ps, *pt;
	int kx;
	cmax=0;		ps = seq;
	while(1){
		if(strstr(ps,"ATG")==NULL)	break;						/* No more start codons */
		ps = strstr(ps,"ATG");
		clen = strlen(seq)+1;
		/* Look for an inframe TAA */
		pt = ps;
		while(1){
			if(strstr(pt,"TAA")==NULL)	break;
                        pt = strstr(pt,"TAA");
			if(pt>ps && (pt-ps)%3==0){
				if((pt-ps)<clen)  	clen = pt-ps;
				break;
			}
			++pt;
		}
		/* Look for an inframe TAG */
		pt = ps;
		while(1){
			if(strstr(pt,"TAG")==NULL)	break;
                        pt = strstr(pt,"TAG");
			if(pt>ps && (pt-ps)%3==0){
				if((pt-ps)<clen)  	clen = pt-ps;
				break;
			}
			++pt;
		}
		/* Look for an inframe TGA */
		pt = ps;
		while(1){
			if(strstr(pt,"TGA")==NULL)	break;
                        pt = strstr(pt,"TGA");
			if(pt>ps && (pt-ps)%3==0){
				if((pt-ps)<clen)  	clen = pt-ps;
				break;
			}
			++pt;
		}
		if(clen<strlen(seq)+1 && clen>cmax)	cmax = clen;
		++ps;
	}
	if(cmax<=0)	return;
	/* Now that cmax (max ORF len in bp -3 for the stop codon) has been identified, find the corresponding ORF */
	ps = seq;
	while(1){
		if(strstr(ps,"ATG")==NULL)	break;						/* No more start codons */
		ps = strstr(ps,"ATG");
		clen = strlen(seq)+1;
		/* Look for an inframe TAA */
		pt = ps;
		while(1){
			if(strstr(pt,"TAA")==NULL)	break;
                        pt = strstr(pt,"TAA");
			if(pt>ps && (pt-ps)%3==0){
				if((pt-ps)<clen){
				  	clen = pt-ps;
					if(clen==cmax){	/* Found */
						for(kx=0; kx<=cmax+2; kx++){	/* cmax,cmax+1,cmax+2: Stop codon */
							cseq[kx] = *ps;
							++ps;
						}
						return;
					}
				}
				else	break;
			}
			++pt;
		}
		/* Look for an inframe TAG */
		pt = ps;
		while(1){
			if(strstr(pt,"TAG")==NULL)	break;
                        pt = strstr(pt,"TAG");
			if(pt>ps && (pt-ps)%3==0){
				if((pt-ps)<clen){
				  	clen = pt-ps;
					if(clen==cmax){	/* Found */
						for(kx=0; kx<=cmax+2; kx++){	/* cmax,cmax+1,cmax+2: Stop codon */
							cseq[kx] = *ps;
							++ps;
						}
						return;
					}
				}
				else	break;
			}
			++pt;
		}
		/* Look for an inframe TGA */
		pt = ps;
		while(1){
			if(strstr(pt,"TGA")==NULL)	break;
                        pt = strstr(pt,"TGA");
			if(pt>ps && (pt-ps)%3==0){
				if((pt-ps)<clen){
				  	clen = pt-ps;
					if(clen==cmax){	/* Found */
						for(kx=0; kx<=cmax+2; kx++){	/* cmax,cmax+1,cmax+2: Stop codon */
							cseq[kx] = *ps;
							++ps;
						}
						return;
					}
				}
				else	break;
			}
			++pt;
		}
		++ps;
	}
        return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
         strncpy(rstr,empty,1001);
        while((fgets(rstr,1000,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%5000==0)	printf("Now processing cDNA %d...\n",counter);

			/* Process the last seq. */
			if(counter>=2){
				orffind();	/* The longest ORF of the forward seq -> cseq */
				if(strlen(cseq)>0){
					/* Output the longest ORF */
					++okcount;
					fputs(header,fw);
					for(x=0; x<200000; x++){
						if(isprint(cseq[x])==0)		break;
						fputc(cseq[x],fw);
						if((x+1)%60==0)		fputs("\n",fw);
					}
					if(x%60!=0)	fputs("\n",fw);
				}	
			}

			/* Saving the header line */
			strncpy(header,empty,1001);		strcpy(header,rstr);

                        /* Reinitializations */
                        for(x=0; x<200000; x++){
			        seq[x] = empty[0];
				cseq[x] = empty[0];
			}
                        numbseq = -1;
		}
                else{   /* Read in the DNA seq. */
                        for(x=0; x<1000; x++){
                                if(isprint(rstr[x])==0)         break;
                                ++numbseq;
                                if(numbseq>=200000){    printf("The # of nucleotides exceeded 200,000\n");   exit(1);        }
                                seq[numbseq] = rstr[x];
                        }
                }

                strncpy(rstr,empty,1001);
        }
	/* Process the very last seq. */
	if(counter>=2){
		orffind();	/* The longest ORF of the forward seq -> cseq */
		if(strlen(cseq)>0){
			/* Output the longest ORF */
			++okcount;
			fputs(header,fw);
			for(x=0; x<200000; x++){
				if(isprint(cseq[x])==0)		break;
				fputc(cseq[x],fw);
				if((x+1)%60==0)		fputs("\n",fw);
			}
			if(x%60!=0)	fputs("\n",fw);
		}	
	}
	printf("%d cDNAs out of %d had their ORFs identified\n",okcount,counter);
        fclose(fr);        fclose(fw);
	return 0;
}
