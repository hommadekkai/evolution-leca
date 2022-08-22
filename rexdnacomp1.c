#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fe;
char rstr[1001], estr[501], empty[1010]="",tempstr[21];
char fnamer[122]="/Users/khomma/work/ensembl/rice/rabbcdnalist1";	/* Abbreviated cDNAs */
char fnamew[122]="/Users/khomma/work/ensembl/rice/rexdnacomp1";
char fnamee[100]="/Users/khomma/work/ensembl/rice/ralluniqexon4AB";
char seq[200000], trid[31], nucl[10]={"XATGC"};
char binstr[10][20]={"All","-120","-240","-360","-480","-600","-720","-840","-960","961bp-"};
char cstr[5][30]={"All","Internal exons","5'-term exons","3'-term exons","Single exons"};
char chistr[3][20]={"p<0.05","p<0.01","p<0.001"};
int p, q, r, s, t, u, v, w, x, y, z, offset, maxseq=0, numbseq=-1, counter=0, hitcount=0, hitflag, nclass;
double gu, gv, gw, gx, gy, gz;
double dna[5][10][5], obs11, obs12, obs21, obs22, obs1T, obs2T, obsT1, obsT2, obsTT,  exp11, exp12, exp21, exp22, chisq, obsTT;
double chiinv[3]={3.841, 6.636, 10.828};
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
                if(catcount==1)         strcpy(vlist.chr,tempstr);     		/* Chromosome name */
                else if(catcount==2)    strcpy(vlist.gname,tempstr);       	/* Gene name */
                else if(catcount==3)    strcpy(vlist.tname,tempstr);       	/* Transcript name */
                else if(catcount==4)    vlist.cdslen = atoi(tempstr);       	/* ##&& CDS length (nt) */
                else if(catcount==5)    vlist.start = atoi(tempstr);       	/* Start genome address */
                else if(catcount==6)    vlist.end = atoi(tempstr);         	/* End genome address */
                else if(catcount==7)    strcpy(vlist.sense,tempstr);       	/* Sense */
                else if(catcount==8)    strcpy(vlist.class,tempstr);       	/* Class: I, F, T, or S */
                else if(catcount==9)    vlist.resstart = atoi(tempstr);    	/* ##&& The start CDS res. (nt) of the exon */
                else if(catcount==10){
                        vlist.resend = atoi(tempstr);      			/* ##&& The end CDS res. (nt) of the exon */
                        return;
                }
                for(x=0; x<50; x++)     if(isspace(estr[x+offset])==0)  break;  
                offset = x + offset;
        }
        return;
}
void exsearch(void){	/* Look for the corresponding exon info. */
	fseek(fe,0L,SEEK_SET);
	strncpy(estr,empty,501);	hitflag = 0;
       	while((fgets(estr,500,fe))!=NULL){
		if(estr[0]!='>' && strstr(estr,trid)!=NULL){	/* Possible hit */
			sread(estr);
			if(strcmp(vlist.tname,trid)==0 && numbseq+1==vlist.cdslen){	/* Verified */
				/* Note that the CDS length must agree with the sequece length */
				if(vlist.class[0]=='I')		nclass = 1;
				else if(vlist.class[0]=='F')	nclass = 2;
				else if(vlist.class[0]=='T')	nclass = 3;
				else if(vlist.class[0]=='S')	nclass = 4;
				else{
					printf("What is this class?:%s from %s",vlist.class,estr);
					exit(1);
				}
	
				if(hitflag==0){
					++hitcount;
					hitflag = 1;
				}
				w = abs(vlist.resend - vlist.resstart) + 1; 	/* Exon length in bp */
				q = (w-1)/120 + 1;         /* Simple exon length bin */
				if(q>9)     q = 9;
				s = vlist.resstart -1;
				t = vlist.resend -1;
				for(x=s; x<=t; x++){
					for(z=1; z<=4; z++){
						if(toupper(seq[x])==nucl[z]){
							dna[nclass][q][z] += 1.;
							goto nextx;
						}	
					}
nextx:					x = x;
				}
			}
			else if(hitflag==1)	return; /* Unusual. The end of the exon lists of the trid */
		}	
		else if(estr[0]!='>' && strstr(estr,trid)==NULL && hitflag==1)	return;
			/* The end of the exon lists of the trid */
nexte:		strncpy(estr,empty,501);
	}
	return;
}
void statchi(int tx){	/* Chi sq statistics from bin tx to bin 9 */
	for(z=1; z<5; z++){	/* For each nucleotide */
		fputc(nucl[z],fw);	fputs("\t",fw);
		obsTT = dna[0][0][0];
		obs11 = 0.;
		for(y=tx; y<=9; y++)	obs11 += dna[1][y][z];
		obsT1 = dna[0][0][z];
		obs21 = obsT1 -obs11;
		obs1T = 0.;
		for(y=tx; y<=9; y++)	obs1T += dna[1][y][0];
		obs12 = obs1T - obs11;
		obsT2 = obsTT -obsT1;
		obs2T = obsTT -obs1T;
		obs22 = obs2T -obs21;
		if(abs(obs11+obs12+obs21+obs22-obsTT)>1.){  printf("Error at #1100\n");      exit(1);        }	/* ###&&& */
		else if(obsTT<1.){
			fputs("\n",fw);
			continue;
		}

		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",obs1T);	fputs(tempstr,fw);	fputs("\t",fw);
		if(obs1T<1.)	fputs("ND\t",fw);
		else{	
			gu = 100.*obs11/obs1T;          /* Freq. of ATGC occurence in internal exons length bin y */
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gu);	fputs(tempstr,fw);	fputs("\t",fw);
		}

		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",obs2T);	fputs(tempstr,fw);	fputs("\t",fw);
		if(obs2T<1.)      fputs("ND\t",fw);
		else{
			gv = 100.*obs21/obs2T;		/* Freq. of ATGC occurence in others */
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gv);	fputs(tempstr,fw);	fputs("\t",fw);
		}
		exp11 = obsTT*((obs11 + obs12)/obsTT)*((obs11 + obs21)/obsTT);
		exp21 = obsTT*((obs21 + obs22)/obsTT)*((obs11 + obs21)/obsTT);
		exp12 = obsTT*((obs11 + obs12)/obsTT)*((obs12 + obs22)/obsTT);
		exp22 = obsTT*((obs21 + obs22)/obsTT)*((obs12 + obs22)/obsTT);
		if(abs(exp11+exp12-obs11-obs12)>1.){   printf("Error at #1200\n");      exit(1);        }		 /* ###&&& */
		if(abs(exp11+exp21-obs11-obs21)>1.){   printf("Error at #1400\n");      exit(1);        }		 /* ###&&& */
		if(abs(exp12+exp22-obs12-obs22)>1.){   printf("Error at #1500\n");      exit(1);        }		 /* ###&&& */
		
                if(exp11>0. && exp12>0. && exp21>0. && exp22>0.){	/* ###&&& */
                        gw = (obs11 -exp11)*(obs11 -exp11)/exp11;
                        gx = (obs12 -exp12)*(obs12 -exp12)/exp12;
                        gy = (obs21 -exp21)*(obs21 -exp21)/exp21;
                        gz = (obs22 -exp22)*(obs22 -exp22)/exp22;
                        chisq = gw +gx +gy +gz;
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",chisq);	fputs(tempstr,fw);	fputs("\t",fw);
			for(p=2; p>=0; p--){	/* Diff. threshold of chi-sq values */
                               	if((p==2 && chisq>=chiinv[p]) || (p<=1 && chisq>=chiinv[p] && chisq<chiinv[p+1])){
					fputs("Significant at ",fw);
					fputs(chistr[p],fw);
					break;
				}
			}
			fputs("\n",fw);
		}
		else fputs("\n",fw);
	}
	return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fe=fopen(fnamee,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamee);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(x=0; x<5; x++){
		for(y=0; y<10; y++){
			for(z=0; z<5; z++)	dna[x][y][z] = 0.;
		}
	}
	for(x=0; x<200000; x++)		seq[x] = empty[0]; 
        strncpy(rstr,empty,1001);
        while((fgets(rstr,1000,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%5000==0)	printf("Now processing cDNA %d...\n",counter);

			/* Process the last variant */
			if(counter>=2){
				if(numbseq+1>=maxseq)	maxseq = numbseq + 1;
				exsearch();
			}

			/* Read in the transcript ID */
			offset = 1;
			strncpy(trid,empty,31);
			for(x=0; x<30; x++){
				if(isspace(rstr[x+offset])!=0)	break;		/* #### */
				trid[x] = rstr[x+offset];	
			}

			/* Reinitializations */
			for(x=0; x<200000; x++)		seq[x] = empty[0]; 
			numbseq = -1;
		}
		else{	/* Read in the DNA seq. */
			for(x=0; x<1000; x++){
				if(isprint(rstr[x])==0)		break;
				++numbseq;
				if(numbseq>=200000){	printf("The # of nucleotides exceeded 200,000 in %s\n",trid);	exit(1);	}
				seq[numbseq] = rstr[x];
			}
		}

                strncpy(rstr,empty,1001);
        }
	/* Process th very last variant */
	if(counter>=2){
		if(numbseq+1>=maxseq)	maxseq = numbseq + 1;
		exsearch();
	}
	printf("%d out of %d cDNAs were processed and the max # of nucleotides is %d\n",hitcount,counter,maxseq);

	/* Sum over ATGC */
        for(x=1; x<5; x++){
                for(y=1; y<10; y++){
                        for(z=1; z<5; z++)      dna[x][y][0] += dna[x][y][z];
                }
        }
	/* Sum over exon length bins */
        for(x=1; x<5; x++){
               	for(y=1; y<10; y++){
               		for(z=0; z<5; z++)	dna[x][0][z] += dna[x][y][z];
                }
        }
	/* Sum over classes */
        for(x=1; x<5; x++){
               	for(y=0; y<10; y++){
               		for(z=0; z<5; z++)	dna[0][y][z] += dna[x][y][z];
                }
        }

	/* chi-sq analysis */
	fputs("Nucleotide\t# in the bin\tFreq. in the bin\t# in others\tFreq. in others\tChi sq.\tSignificance\n",fw);
	for(y=1; y<10; y++){	/* For each exon length bins */
		fputs(">Internal exons bin ",fw);	
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",y);	fputs(tempstr,fw);	fputs(": ",fw);
		fputs(binstr[y],fw);
		if(y!=9)	fputs("bp",fw);
		fputs("\n",fw);
		for(z=1; z<5; z++){	/* For each nucleotide */
			fputc(nucl[z],fw);	fputs("\t",fw);
			obsTT = dna[0][0][0];
			obs11 = dna[1][y][z];
			obsT1 = dna[0][0][z];
			obs21 = obsT1 -obs11;
			obs1T = dna[1][y][0];
			obs12 = obs1T - obs11;
			obsT2 = obsTT -obsT1;
			obs2T = obsTT -obs1T;
			obs22 = obs2T -obs21;
			if(abs(obs11+obs12+obs21+obs22-obsTT)>0.1){  printf("Error at #100\n");      exit(1);        }
			else if(obsTT<1.){
				fputs("\n",fw);
				continue;
			}
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",obs1T);	fputs(tempstr,fw);	fputs("\t",fw);
			if(obs1T<1.)	fputs("ND\t",fw);
			else{	
				gu = 100.*obs11/obs1T;          /* Freq. of ATGC occurence in internal exons length bin y */
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gu);	fputs(tempstr,fw);	fputs("\t",fw);
			}

			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",obs2T);	fputs(tempstr,fw);	fputs("\t",fw);
			if(obs2T<1.)      fputs("ND\t",fw);
			else{
				gv = 100.*obs21/obs2T;		/* Freq. of ATGC occurence in others */
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gv);	fputs(tempstr,fw);	fputs("\t",fw);
			}
			exp11 = obsTT*((obs11 + obs12)/obsTT)*((obs11 + obs21)/obsTT);
			exp21 = obsTT*((obs21 + obs22)/obsTT)*((obs11 + obs21)/obsTT);
			exp12 = obsTT*((obs11 + obs12)/obsTT)*((obs12 + obs22)/obsTT);
			exp22 = obsTT*((obs21 + obs22)/obsTT)*((obs12 + obs22)/obsTT);
			if(abs(exp11+exp12-obs11-obs12)>0.1){   printf("Error at #200\n");      exit(1);        }
			if(abs(exp11+exp21-obs11-obs21)>0.1){   printf("Error at #400\n");      exit(1);        }
			if(abs(exp12+exp22-obs12-obs22)>0.1){   printf("Error at #500\n");      exit(1);        }
			
                        if(exp11>0.00001 && exp12>0.00001 && exp21>0.00001 && exp22>0.00001){
                                gw = (obs11 -exp11)*(obs11 -exp11)/exp11;
                                gx = (obs12 -exp12)*(obs12 -exp12)/exp12;
                                gy = (obs21 -exp21)*(obs21 -exp21)/exp21;
                                gz = (obs22 -exp22)*(obs22 -exp22)/exp22;
                                chisq = gw +gx +gy +gz;
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",chisq);	fputs(tempstr,fw);	fputs("\t",fw);
				for(p=2; p>=0; p--){	/* Diff. threshold of chi-sq values */
                                	if((p==2 && chisq>=chiinv[p]) || (p<=1 && chisq>=chiinv[p] && chisq<chiinv[p+1])){
						fputs("Significant at ",fw);
						fputs(chistr[p],fw);
						break;
					}
				}
				fputs("\n",fw);
			}
			else fputs("\n",fw);

		}
	}

	/* Bins 7-9 */
	fputs(">Internal exons bins 7-9: 721bp-\n",fw);	
	statchi(7);	/* Bins 7-9 */

	/* Bins 4-9 */
	fputs(">Internal exons bins 4-9: 361bp-\n",fw);	
	statchi(4);	/* Bins 4-9 */

        /* Output the overall freqs. of each amino acid */
        fputs("Nucleotide\tInt Ex Avg\tMulti Ex Avg\tAll Avg\n",fw);
        for(z=1; z<5; z++){    /* For each aa */
                fputc(nucl[z],fw);  fputs("\t",fw);
                gx = 100.*dna[1][0][z]/dna[1][0][0];          /* Internal exons */
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gx);     fputs(tempstr,fw);      fputs("\t",fw);
                gx = 100.*(dna[0][0][z]-dna[4][0][z])/(dna[0][0][0]-dna[4][0][0]);          /* Multiple exons */
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gx);     fputs(tempstr,fw);      fputs("\t",fw);
                gx = 100.*dna[0][0][z]/dna[0][0][0];;         /* All exons */
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gx);     fputs(tempstr,fw);      fputs("\n",fw);
        }

        fclose(fr);        fclose(fw);
        fclose(fe);
	return 0;
}
