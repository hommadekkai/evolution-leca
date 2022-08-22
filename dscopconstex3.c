#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fl;
char rstr[122], empty[500]="",  tempstr[21];
char idstr[51], lstr[201];
char fnamer[122]="/Users/khomma/work/ensembl/dmel/dscoplist12";
char fnamel[100]="/Users/khomma/work/ensembl/dmel/dalluniqexon3M";      /* Non-redundant list of exons, with int CDS length */
char fnamew[122]="/Users/khomma/work/ensembl/dmel/dscopconstex3";      /* #+# */
int t, u, v, w, x, y, z, offset, length, scmax=-1, exmax, hitflag, counter=0, hitcount=0, nbin;
int scop[10000][2];	/* 0: Start, 1: End of SCOP alignment */
int exon[10000][4];	/* 0: Start, 1: End of the exon-encoded residues, 2: Classification of each exon, 3: Exon length bin:
			 Classification =1: 5'-term exon, =2: interna; exon, =3: 3'-term exon, =0: Single exon */
int exlen, excdslen, intcds, intstart, intend, oldintcds, scnumb, sclen, maxscop=0;	/* #"# */
int totnume;	/* #++# */
double gw, gx, gy, gz, scaling=1000., exonprob;
double nobs[10], nexp[10], numb[10];
	/* 0:All exon len., 1: 1~120, 2: ~240, 3: ~360, 4: ~480, 5: ~600, 6: ~720, 7: ~840, 8: ~960, 9: 961~ bp */
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
	for(z=0; z<10; z++){	/* 0:All ex len, 1:1~120, 2:~240, 3:~360, 4:~480, 5:~600, 6:~720, 7:~840, 8:~960, 9:961~ bp */
		numb[z] = 0.;
		nobs[z] = 0.;
		nexp[z] = 0.;
	}

        strncpy(rstr,empty,122);
        while((fgets(rstr,121,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%1000==0)	printf("Now processing variant %d\n",counter);

                        /* ### DEBUGGED to convert ENSG31208_272into FBgn0031208 for example */
                        strncpy(tempstr,empty,21);
                        for(x=0; x<50; x++)     if(isdigit(rstr[x])!=0)         break;
                        offset = x;

                        /* Get 31208 from >ENSG31208_272, for example */
                        for(x=0; x<20; x++){
                                if(isdigit(rstr[x+offset])==0)  break;
                                tempstr[x] = rstr[x+offset];
                        }
                        w = strlen(tempstr);
                        strncpy(idstr,empty,51);        strcpy(idstr,"FBgn");
                        if(w<=6){
                                for(x=0; x<=6-w; x++)  strcat(idstr,"0");
                        }
                        strcat(idstr,tempstr);
                        /* ### */
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
			if(scmax>=10000){	printf("The # of SCOP alignments in %s exceeded 10,000\n",idstr);	exit(1);	}
			scopparse();
		}
		else if(strncmp(rstr,"//",2)==0){	/* Analyze the correspondence */
			/* ### DO NOT NEGLECT proteins without SCOP assingments */
			if(strlen(idstr)<2 || length<0)		goto reinit;
			/* Find the exon data */
			exmax = -1;	hitflag = 0;
			for(x=0; x<10000; x++){
				for(y=0; y<4; y++)	exon[x][y] = 0;
			}
			fseek(fl,0L,SEEK_SET);
			strncpy(lstr,empty,201);
        		while((fgets(lstr,200,fl))!=NULL){
				if(strstr(lstr,idstr)!=NULL){	/* Possible hit. Check the total length */
					w = exonparse(4);	/* Total CDS length (aa) */
					if(length-1==w){	/* $$$ HIT! */
						++exmax;
						if(exmax>=10000){	
							printf("The # of exon of %s (length=%d) exceeded 10,000\n",idstr,length);
							exit(1);
						}
						exon[exmax][0] = exonparse(10);		/* Exon res start */
						exon[exmax][1] = exonparse(11);		/* Exon res end */
						t = exonparse(8);		/* Exon number */
						u = exonparse(9);		/* Total exon number */
						if(u==1)	exon[exmax][2] = 0;	/* Single exon */
						else if(t==1)	exon[exmax][2] = 1;	/* 5'-term exon */
						else if(t==u)	exon[exmax][2] = 3;	/* 3'-term exon */
						else		exon[exmax][2] = 2;	/* Internal exon */

						/* Determine the exon length bin and register it as exon[][3] */
						t = exonparse(5);		/* Exon start (nt) */
						u = exonparse(6);		/* Exon end (nt) */
						exlen = (abs)(u - t) + 1;
						nbin = (exlen-1)/120 +1;
						if(nbin>=10)	nbin = 9;
						exon[exmax][3] = nbin;

						intstart = exonparse(12);	/* #"# */
						intend = exonparse(13);		/* #"# */
						intcds = intend - intstart +1;	/* #"# */
						if(hitflag==0)	oldintcds = intcds;	/* The first transcript */
						else if(intcds!=oldintcds)	break;	/* The end of the transcript */

						hitflag = 1;		/* Moved */
					}
					else if(length-1!=w && hitflag==1)	break;	/* $$$ The end of the transcript */
				}
				else if(hitflag==1 && strstr(lstr,idstr)==NULL) break;	/* The end of a gene */

				strncpy(lstr,empty,201);
			}
			if(hitflag==0 || exmax<0)	goto reinit;
			else if(intcds<=0)		goto reinit;	/* No CDS encoded by interna exons */
			/* ### DEBUGGED ### */

			/* Examine the observed # of exon boundaries that are inside SCOP domains and the expected #s */
			if(exon[exmax][2]<=0)    goto reinit;    /* Neglecting single-exon genes */

			/* Eliminate the SCOP domains that are not entirely encoded by internal exons */
			scnumb = scmax +1;	/* The number of acceptable SCOP domains */
			if(scnumb==0)   goto shunt;     /* ### */
                        /* #+# Those partially encoded by internal exons are partially left */
                        for(x=0; x<=scmax; x++){
                                if(intstart<=scop[x][0] && scop[x][0]<=intend){ /* The left end of the SCOP is OK */
                                        if(intstart<=scop[x][1] && scop[x][1]<=intend){ /* The right end of the SCOP is also OK */
                                                /* OK, both ends of the scop domain are encoded by internal exons ! */
                                                if(scop[x][1]-scop[x][0]+1 > maxscop)   maxscop = scop[x][1]-scop[x][0]+1;
                                        }
                                        else{   /* Only the left end of the SCOP domain is acceptable */
                                                scop[x][1] = intend;    /* The right end is shortened */
                                        }
                                }
                                else if(intstart<=scop[x][1] && scop[x][1]<=intend){    /* ONLY the right end of the SCOP is OK */
                                        scop[x][0] = intstart;  /* The left end is shortened */
                                }
                                else{ /* Both ends of the SCOP domain is not OK.  Eliminate the SCOP alignment */
                                        --scnumb;
                                        scop[x][0] = 0;         scop[x][1] = 0;
                                }
                        }
                        /* #+# */
shunt:                  /* ### DO NOT skip even if no acceptable SCOP domains exist */
			++hitcount;

                        /* #-# Calculate the expected and obs freqs in the internal exon renge only (no okflag used) */
                        for(y=0; y<=exmax; y++){        /* Examine the right AND left boundaries of all exons */
                                if(exon[y][2]!=2)       continue;       /* Only internal exons are examined */
                                nbin = exon[y][3];
                                excdslen = abs(exon[y][1] - exon[y][0])+ 1;     /* Exon length (aa) */
                                /* ### */
                                if(scmax<0){
                                        exonprob = 0.;
                                        goto skip;
                                }
                                /* ### */
				totnume = 0;	/* #++# */
                                for(x=0; x<=scmax; x++){        /* One SCOP alignment at a time */
                                        if(scop[x][0]==0 || scop[x][1]==0)      continue;       /* Eliminated SCOP domain */

                                        /* #+# Small corrections for the partially analyzed SCOP doamins */
                                        v = 0;          w = 0;
                                        if(scop[x][0]==intstart)        v = 1;
                                        if(scop[x][1]==intend)          w = 1;
                                        sclen = scop[x][1] - scop[x][0] + 1;    /* SCOP length */
                                        if(v==1)        ++sclen;
                                        if(w==1)        ++sclen;

                                        t = sclen - excdslen -1;
                                        u = intcds -excdslen +1;
					if(t<=0)	continue;	/* #++# The SCOP domain is too short */
                                        else if(t>u)	continue;	/* An aberrent case */
					else		totnume += t;	/* #++# */

                                        /* Check if SCOP domain x contains exon y */
                                        if(scop[x][0]-v<exon[y][0] && exon[y][0]<scop[x][1]+w){
                                                /* HIT! The left exon boundary is inside SCOP domain x */
                                                /* Examine the right boundary of the exon */

                                                if(scop[x][0]-v<exon[y][1] && exon[y][1]<scop[x][1]+w){
                                                        /* HIT! The right exon boundary is inside THE SAME SCOP domain x */
                                                        nobs[nbin] += 1.;
                                                }
                                        }
                                        /* #+# */
                                }
				/* #++# */
				u = intcds -excdslen +1;	
				if(totnume>0 && u>0)	exonprob = (double)totnume/(double)u;
				/* #++# */
				
skip:                           numb[nbin] += 1.;       /* ### Count the # of samples here */
                                nexp[nbin] += exonprob;
                        }      
                        /* #-# */

reinit:			/* Reinitializations */
			length = 0;	strncpy(idstr,empty,51);	scmax = -1;
			for(x=0; x<10000; x++){
				for(y=0; y<2; y++)	scop[x][y] = -1;
			}
		}
                strncpy(rstr,empty,122);
        }
	/* Sum over all the exon length bins */
	for(x=1; x<10; x++){
		numb[0] += numb[x];
		nobs[0] += nobs[x];
		nexp[0] += nexp[x];
	}
	/* Calculate and output the fractions */
	fputs("Exon len.\tNumber\tObs.freq.(%)\tExp.freq.(%)\tChi sq.\tObs/Exp\n",fw);
	for(z=0; z<10; z++){
                if(z==0)        fputs("All\t",fw);
                else{
                        strncpy(tempstr,empty,21);
                        sprintf(tempstr,"%d",120*z-119);
                        fputs(tempstr,fw);
                        if(z==9)        fputs(" bp~\t",fw);
                        else{
                                fputs("~",fw);
                                strncpy(tempstr,empty,21);
                                sprintf(tempstr,"%d",120*z);
                                fputs(tempstr,fw);
                                fputs(" bp\t",fw);
                        }
                }
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",numb[z]);	fputs(tempstr,fw);	fputs("\t",fw);
		if(numb[z]<=0){		fputs("\n",fw);		continue;	}	

                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.3f",100.*nobs[z]/numb[z]);      /* Observed (%) */
                fputs(tempstr,fw);      fputs("\t",fw);

                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.3f",100.*nexp[z]/numb[z]);    /* Expected (%) */
                fputs(tempstr,fw);      fputs("\t",fw);

                if(nexp[z]<0.0001){   	fputs("-\n",fw);	continue;	}
                gx = (nobs[z] - nexp[z])*(nobs[z] - nexp[z])/nexp[z];
                gy = (nexp[z]-nobs[z])*(nexp[z]-nobs[z])/(numb[z]-nexp[z]);
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.3f",gx+gy);                  /* Chi sq */
                fputs(tempstr,fw);      fputs("\t",fw);

                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.3f",nobs[z]/nexp[z]);    /* Obs/Exp */
                fputs(tempstr,fw);      fputs("\n",fw);
	}
	printf("%d out of %d variants were analyzed. The max SCOP length is %d aa\n",hitcount,counter,maxscop);
        fclose(fr);	fclose(fl);        fclose(fw);
	return 0;
}
