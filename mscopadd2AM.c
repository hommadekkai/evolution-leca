/* mscopadd2AM.c		This program determines the max internal exon length downstream of SCOP domains		220811
	 basic.c -> hscopexon1.c: A new program [R2102 I].
	 -> hscopexbc1.c: Revised to identify the correspondence between SCOP and exon boundaries [R2102 I].
	 -> hscopinexbc1.c: Revised to limit analyses to internal exons, excluding terminal and single exons[R2102 I].
	 -> hscopinexbc1M.c: Revised to use the corrected list of exons AND to use not 80 but 40 aa bins[R2102 I].
	 -> hexonscop1.c: Completely revised to analyze the freq. of SCOP domains split to multiple exons [R2106 I].
	 -> hexonscopavg1.c: Revised to use the AVG expectation freq. [R2106 I].
	 -> hexonscopavg3.c: Revised to use BOTH right and left exon boundaries [R2106 I].
	 -> hrexonscopavg3.c: Revised to use a list of internal exons by a corrected definition R2204 I]. (DEBUGGED)
	 -> hscopnumb1.c: Revised to find the dependence of SCOP # on the longest internal exon of each gene (###) [R2208 I].
	 -> hscopnumb3.c: Revised to limit analysis to SCOP domains that are 3' to each exon and to calculate exp values (#+#) [R2208 I].
	 -> hscopnumb5.c: Revised to count SCOP domains that are encoded by the longest int exon or those 3' to it (#-#) [R2208 I].
	 -> hscopnumb5A.c: Modified to calculate Obs/ExpAll with t values and Obs/SimpleExpAll with t values [R2208 I].
	 -> hmultiscop1.c: Revised to compute the max internal exon legth between SCOP domains and the exp values [R2208 I].
	 -> hscopadd1.c: Revised to compute the dists. of max int exon len. for those with and without downstream SCOP [R2208 I].
	 -> mscopadd1.c: Revised to process mouse instead of human data [R2208 VII].
	 -> mscopadd2.c: Modified not to include genes without SCOP domains in Group B (#-#) [R2208 VII].
	 -> mscopadd2A.c: Modified to relax the requirement for the acceptable internal exons (#+#) [R2208 VII].
	 -> mscopadd2AM.c: Modified to use a redundant instead of nonredundant exon list (#++#) [R2208 VII]. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fl;
char rstr[122], empty[500]="",  tempstr[21];
char gname[100][21], lstr[201];
char fnamer[122]="/Users/khomma/work/ensembl/mouse/mscoplist12";
char fnamel[100]="/Users/khomma/work/ensembl/mouse/mallexon3";     	/* #++# Redundant list of exons, By new def. */
char fnamew[122]="/Users/khomma/work/ensembl/mouse/mscopadd2AM";		/* #+++# */
int t, u, v, w, x, y, z, offset, length, scmax=-1, exmax, hitflag, counter=0, hitcount=0, nbin, gmax;
int scop[10000][2];	/* 0: Start, 1: End of SCOP alignment */
int exon[10000][4];	/* 0: Start, 1: End of the exon-encoded residues, 2: Classification of each exon, 3: Exon length bin:
			 Classification =1: 5'-term exon, =2: interna; exon, =3: 3'-term exon, =0: Single exon */
int exlen, excdslen, binnumb, maxintex, sep, instart, inend, nsampl, totlen;
double gt, gu, gv, gw, gx, gy, gz, gs1, gs2, gs12, scaling=1000.;
double lobs0[10], lobs1[10], lobs2[10];	/* Tot. # of samples, sum of max int lengths, sum of max int len^2 in each bin */
double lnoscp0[10], lnoscp1[10], lnoscp2[10], startsep;
void scopparse(void)
{
	int i;
	char listread[21];
        offset = 15;
        strncpy(listread,empty,21);
	for(i=0; i<20; i++){
		if(rstr[i+offset]=='-')		break;
		listread[i] = rstr[i+offset];
	}
	offset = i + offset + 2;
	scop[scmax][0] = atoi(listread);

        strncpy(listread,empty,21);
	for(i=0; i<20; i++){
		if(rstr[i+offset]=='|')		break;
		listread[i] = rstr[i+offset];
	}
	scop[scmax][1] = atoi(listread);
        return;
}
int exonparse(int k)
{
	int i, j=0, l;
	char tempexon[51];
	offset = 0;
	while(1){
		++j;
		if(isprint(lstr[offset])==0)	break;

		strncpy(tempexon,empty,51);
		for(i=0; i<=50; i++){
			if(isspace(lstr[i+offset])!=0)	break;
			else	tempexon[i] = lstr[i+offset];
		}
		offset = i + offset;
		if(j==k){
			l = atoi(tempexon);
			return l;
		}
		for(i=0; i<50; i++)	if(isspace(lstr[i+offset])==0)  break;
		offset = i + offset;
	}
	return -1;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fl=fopen(fnamel,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamel);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs("Max internal exon length (nt) with (A) /without (B) downstream SCOPs\n",fw);
	for(z=0; z<10; z++){	/* Bins of separation between SCOP domains; 0: All, 1:1-40, 2:41-80,..., 8:280-320, 9:321- aa */
		/* The dist of cases with at least one SCOP domain donwstream */
		lobs0[z] = 0.;	lobs1[z] = 0.;	lobs2[z] = 0.;
		/* The dist of cases without SCOP domains donwstream */
		lnoscp0[z] = 0.;	lnoscp1[z] = 0.; 	lnoscp2[z] = 0.;
	}
        length = 0;     scmax = -1;
        for(x=0; x<10000; x++){
                for(y=0; y<2; y++)      scop[x][y] = -1;
        }
        gmax = -1;
        for(x=0; x<100; x++)    strncpy(gname[x],empty,21);
        strncpy(rstr,empty,122);
        while((fgets(rstr,121,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%1000==0)	printf("Now processing variant %d\n",counter);
		}
                else if(strncmp(rstr,"GeneName",8)==0){ /* Get gname: ENSMUSG00000091476 for instance */
                        offset = 10;
                        gmax = -1;
                        for(x=0; x<100; x++)    strncpy(gname[x],empty,21);
                        while(1){
                                ++gmax;
                                if(gmax>=100){  printf("The # of gene names exceeded 100 in %s",rstr);          exit(1);        }
                                for(x=0; x<20; x++){
                                        if(isprint(rstr[x+offset])==0)          goto nextr;
                                        else if(isspace(rstr[x+offset])!=0)     break;
                                        gname[gmax][x] = rstr[x+offset];
                                }
                                offset = x + offset + 1;
                        }
                }
		else if(strncmp(rstr,"LENGTH",6)==0){
			strncpy(tempstr,empty,21);
			offset = 10;
			for(x=0; x<20; x++){
				if(isprint(rstr[x+offset])==0)	break;
				tempstr[x] = rstr[x+offset];
			}
			length = atoi(tempstr);
		}
		else if(strncmp(rstr,"RP:SCP:REP",10)==0){
			++scmax;
			if(scmax>=10000){	printf("The # of SCOP alignments in %s exceeded 10,000\n",gname[0]);	exit(1);	}
			scopparse();
		}
		else if(strncmp(rstr,"//",2)==0){	/* Analyze the correspondence */
			if(gmax<0 || length<0)		goto reinit;
			/* Find the exon data */
			exmax = -1;	hitflag = 0;
			for(x=0; x<10000; x++){
				for(y=0; y<4; y++)	exon[x][y] = 0;
			}
			fseek(fl,0L,SEEK_SET);
			strncpy(lstr,empty,201);
                        while((fgets(lstr,200,fl))!=NULL){
                                if(hitflag==0){
                                        for(x=0; x<=gmax; x++){
                                                if(strstr(lstr,gname[x])!=NULL){        /* Possible hit. Check the total length */
                                                        w = exonparse(4);       /* Total CDS length (aa) */
                                                        if(length==w){  /* HIT! */
                                                                hitflag = 1;
                                                                ++exmax;
                                                                if(exmax>=10000){
                                                                        printf("The # of exon of %s exceeded 10,000\n",gname[x]);
                                                                        exit(1);
                                                                }
                                                                exon[exmax][0] = exonparse(10);         /* Exon res start */
                                                                exon[exmax][1] = exonparse(11);         /* Exon res end */
                                                                t = exonparse(8);               /* Exon number */
                                                                u = exonparse(9);               /* Total exon number */
                                                                if(u==1)        exon[exmax][2] = 0;     /* Single exon */
                                                                else if(t==1)   exon[exmax][2] = 1;     /* 5'-term exon */
                                                                else if(t==u)   exon[exmax][2] = 3;     /* 3'-term exon */
                                                                else            exon[exmax][2] = 2;     /* Internal exon */

                                                                /* Determine the exon length bin and register it as exon[][3] */
                                                                t = exonparse(5);               /* Exon start (nt) */
                                                                u = exonparse(6);               /* Exon end (nt) */
                                                                exlen = (abs)(u - t) + 1;
                                                                exon[exmax][3] = exlen;         /* ### */
                                                        }
                                                }
                                        }
                                }
                                else if(hitflag==1){
                                        hitflag = -1;
                                        for(x=0; x<=gmax; x++){
                                                if(strstr(lstr,gname[x])!=NULL){        /* Possible hit. Check the total length */
                                                        w = exonparse(4);       /* Total CDS length (aa) */
                                                        if(length==w){  /* HIT! */
                                                                hitflag = 1;
                                                                ++exmax;
                                                                if(exmax>=10000){
                                                                        printf("The # of exon of %s exceeded 10,000\n",gname[x]);
                                                                        exit(1);
                                                                }
                                                                exon[exmax][0] = exonparse(10);         /* Exon res start */
                                                                exon[exmax][1] = exonparse(11);         /* Exon res end */
                                                                t = exonparse(8);               /* Exon number */
                                                                u = exonparse(9);               /* Total exon number */
                                                                if(u==1)        exon[exmax][2] = 0;     /* Single exon */
                                                                else if(t==1)   exon[exmax][2] = 1;     /* 5'-term exon */
                                                                else if(t==u)   exon[exmax][2] = 3;     /* 3'-term exon */
                                                                else            exon[exmax][2] = 2;     /* Internal exon */

                                                                /* Determine the exon length bin and register it as exon[][3] */
                                                                t = exonparse(5);               /* Exon start (nt) */
                                                                u = exonparse(6);               /* Exon end (nt) */
                                                                exlen = (abs)(u - t) + 1;
                                                                exon[exmax][3] = exlen;         /* ### */
                                                        }
                                                }
                                        }
                                }
                                if(hitflag==-1) break;
                                /* The end of a gene */

                                strncpy(lstr,empty,201);
                        }

			if(hitflag==0 || exmax<0)	goto reinit;
			else if(exon[exmax][2]<=0)  	goto reinit;    /* Neglecting single-exon genes */
			else if(exmax<=1)		goto reinit;   /*  Neglecting two-exon genes */


			/* Identify the longest internal exon length between two neighbouring SCOP domains, maxintex */
			if(scmax<0)     goto reinit;    /* #-# Neglect genes without SCOP domains */

			++hitcount;
			if(scmax>=1){
				for(x=0; x<=scmax-1; x++){	/* One pair of contiguous SCOP domains at a time */
					sep = scop[x+1][0] - scop[x][1] +1;
					if(sep<=0)	continue;	/* sep=0 -> Qualified exons do not exist */
					else 		binnumb = (sep-1)/40 + 1;	/* 40 aa bins */
					if(binnumb>=10)	binnumb = 9;
					maxintex = -1;		/* maxintex becomes the max length (nt) of internal exons */	
					for(y=0; y<=exmax; y++){
						if(exon[y][2]!=2)	continue;	/* Must be an internal exon */
						else if(exon[y][0]<=scop[x][1] && scop[x][1]<exon[y][1]){
							/* #+# Case A */
							w = exon[y][3];		/* Exon len in nt */
							if(w>maxintex)	maxintex = w;
						}
						else if(scop[x][1]<exon[y][0] && exon[y][0]<=scop[x+1][0]){	/* Case BO */
							w = exon[y][3];		/* Exon len in nt */
							if(w>maxintex)	maxintex = w;
						}
					}
					if(maxintex>=0){
						gx = (double)maxintex/scaling;
						lobs0[binnumb] +=1.;
						lobs1[binnumb] +=gx;
						lobs2[binnumb] +=gx*gx;
					}
				}
			}
                        /* #-# Find the max int exon length downstream of the last SCOP
                          or from the start of internal exons in the absence of SCOP */
			/* Find the range of internal exons (instar-inend) first */
			instart = -1; 	inend = -1;
			for(y=0; y<=exmax; y++){
				if(exon[y][2]!=2)	continue;	/* Must be an internal exon */

				if(instart<0)			instart = exon[y][0];
				else if(exon[y][0]<instart)	instart = exon[y][0];

				if(inend<0)			inend = exon[y][1];
				else if(exon[y][1]>inend)	inend = exon[y][1];
			}

			if(scmax<0)			startsep = instart;
			else if(scop[scmax][1]<=inend)	startsep = scop[scmax][1];
			else	goto reinit;		/* No regions */

			if(inend-startsep-1<=0)	goto reinit;
			else			sep = inend -startsep -1;

			binnumb = (sep-1)/40 + 1;	/* 40 aa bins */
			if(binnumb>=10)	binnumb = 9;
			maxintex = -1;
			for(y=0; y<=exmax; y++){	/* Examine the internal exons */
				if(exon[y][2]!=2)	continue;	/* Must be an internal exon */
				else if(startsep<exon[y][1]){
					w = exon[y][3];		/* Exon len in nt */
					if(w>maxintex)	maxintex = w;
				}
			}
			if(maxintex>0){
				gx = (double)maxintex/scaling;
				lnoscp0[binnumb] +=1.;
				lnoscp1[binnumb] +=gx;
				lnoscp2[binnumb] +=gx*gx;
			}

reinit:			/* Reinitializations */
                        length = 0;
                        for(x=0; x<100; x++)    strncpy(gname[x],empty,21);
                        gmax = -1;
                        scmax = -1;
                        for(x=0; x<10000; x++){
                                for(y=0; y<2; y++)      scop[x][y] = -1;
                        }
		}
nextr:		strncpy(rstr,empty,122);
        }
	/* Sum over all max exon length bins */
	for(z=1; z<10; z++){
		lobs0[0] += lobs0[z];
		lobs1[0] += lobs1[z];
		lobs2[0] += lobs2[z];
		lnoscp0[0] += lnoscp0[z];
		lnoscp1[0] += lnoscp1[z];
		lnoscp2[0] += lnoscp2[z];
	}

	fputs("Separation\tA #\tA avg\tA SEM\tB #\tB avg\tB SEM\tA/B\tt\n",fw);
	for(z=0; z<10; z++){
		if(z==0)	fputs("All\t",fw);
		else{
			strncpy(tempstr,empty,21);
			sprintf(tempstr,"%d",40*z-39);
			fputs(tempstr,fw);
			if(z==9)	fputs(" aa~\t",fw);
			else{ 
				fputs("~",fw);
				strncpy(tempstr,empty,21);
				sprintf(tempstr,"%d",40*z);
				fputs(tempstr,fw);
				fputs(" aa\t",fw);
			}
		}
		/* A # */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",lobs0[z]);
		fputs(tempstr,fw);		fputs("\t",fw);

		if(lobs0[z]<0.1){
			fputs("\n",fw);
			continue;
		}

		/* A avg: With downstream SCOP */
		strncpy(tempstr,empty,21);      sprintf(tempstr,"%.1f",scaling*lobs1[z]/lobs0[z]);
		fputs(tempstr,fw);		fputs("\t",fw);

		/* A SEM */
		gx = lobs2[z] - (lobs1[z]*lobs1[z]/lobs0[z]);
		if(gx<0.){	printf("gx = %f at #100\n",gx);		exit(1);	}
		gy = sqrt(gx)/lobs0[z];		/* SEM */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",scaling*gy);
                fputs(tempstr,fw);              fputs("\t",fw);      
		gs1 = gy*sqrt(lobs0[z]);

		/* B # */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",lnoscp0[z]);
		fputs(tempstr,fw);		fputs("\t",fw);

		if(lnoscp0[z]<0.1){
			fputs("\n",fw);
			continue;
		}

		/* B avg: Without downstream SCOP */
		strncpy(tempstr,empty,21);      sprintf(tempstr,"%.1f",scaling*lnoscp1[z]/lnoscp0[z]);
		fputs(tempstr,fw);		fputs("\t",fw);

		/* B SEM */
		gx = lnoscp2[z] - (lnoscp1[z]*lnoscp1[z]/lnoscp0[z]);
		if(gx<0.){	printf("gx = %f at #100\n",gx);		exit(1);	}
		gy = sqrt(gx)/lnoscp0[z];		/* SEM */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",scaling*gy);
       	        fputs(tempstr,fw);              fputs("\t",fw);
		gs2 = gy*sqrt(lnoscp0[z]);

		/* A/B */
		if(lnoscp1[z]>=0.){
			gx = (lobs1[z]/lobs0[z])/(lnoscp1[z]/lnoscp0[z]); 
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gx);
			fputs(tempstr,fw);		fputs("\t",fw);
		}
		else	fputs("\t",fw);

		/* t value between A and B */
		gu = (lobs0[z]-1.)*gs1*gs1;
		gv = (lnoscp0[z]-1.)*gs2*gs2;
		gw = (gu + gv)/(lobs0[z] +lnoscp0[z] -2.);
		if(gw<0.){	printf("gw=%f at #300\n",gw);	exit(1);	}
		gs12 = sqrt(gw);
		if(gs1<=gs2 && gs12<gs1)	printf("S1=%f, S2=%f, S12=%f\n",gs1,gs2,gs12);
		else if(gs1>gs2 && gs12<gs2)	printf("S1=%f, S2=%f, S12=%f\n",gs1,gs2,gs12);
		if(gs1<=gs2 && gs2<gs12)	printf("S1=%f, S2=%f, S12=%f\n",gs1,gs2,gs12);
		else if(gs1>gs2 && gs12<gs12)	printf("S1=%f, S2=%f, S12=%f\n",gs1,gs2,gs12);

		gx = (lobs1[z]/lobs0[z]) - (lnoscp1[z]/lnoscp0[z]);	/* x1 - x2 */
		gy = (1./lobs0[z]) + (1./lnoscp0[z]);		/* 1/n1 + 1/n2 */
		gt = gx/sqrt(gy)/gs12;
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gt);
		fputs(tempstr,fw);		fputs("\n",fw);
	}
	printf("%d out of %d variants were analyzed\n",hitcount,counter);
        fclose(fr);	fclose(fl);        fclose(fw);
	return 0;
}
