#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fl;
char rstr[122], empty[500]="",  tempstr[21];
char gname[100][21], lstr[201];
char fnamer[122]="/Users/khomma/work/ensembl/rat/nscoplist1";
char fnamel[100]="/Users/khomma/work/ensembl/rat/nalluniqexon3";        /* #### Non-redundant list of exons, By new def. */
char fnamew[122]="/Users/khomma/work/ensembl/rat/nscopnumb2";       /* #### */
char gpstr[5][20]={"All","~1/4","~1/2","~3/4","~1"};		/* ### */
int t, u, v, w, x, y, z, offset, length, scmax=-1, exmax, hitflag, counter=0, hitcount=0, nbin;
int scop[10000][2];	/* 0: Start, 1: End of SCOP alignment */
int exon[10000][4];	/* 0: Start, 1: End of the exon-encoded residues, 2: Classification of each exon, 3: Exon length bin:
			 Classification =1: 5'-term exon, =2: interna; exon, =3: 3'-term exon, =0: Single exon */
int exlen, excdslen, binmax, cdslen[1000000], cmax=-1, gmax;
int gpmax[4]={0,295,474,757}, gpnumb;	/* ### Max protein length in each length group */
double gt, gu, gv, gw, gx, gy, gz;
double sumscop0[5][10], sumscop1[5][10], sumscop2[5][10];	/* ### Tot. # of samples, sum of SCOP #s sum of scop #^2 in each bin */
double sumn[5], sumx1[5], sumx2[5], sumy1[5], sumy2[5], sumxy[5], corr, core, cort, scaling=1000.; /* ### */
double sumlen1[5][10], sumlen2[5][10];        /* ### */
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
int intcmp(const void *a, const void *b){
        return *(int *)a - *(int *)b;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fl=fopen(fnamel,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamel);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(y=0; y<5; y++){	/* ### 0: All, 1: 1/4, 2: 1/2, 3: 3/4; 4: 1 */
		sumn[y] =0.;	sumx1[y] = 0.;	 sumx2[y] = 0.;		sumy1[y] = 0.;   sumy2[y] = 0.;		sumxy[y] =0.;
		for(z=0; z<10; z++){	/*  0:All ex len, 1:1~120, 2:~240, 3:~360, 4:~480, 5:~600, 6:~720, 7:~840, 8:~960, 9:961~ bp */
			sumscop0[y][z] = 0.;	sumscop1[y][z] = 0.;	sumscop2[y][z] = 0.;
			sumlen1[y][z] = 0.;     sumlen2[y][z] = 0.;
		}
	}
	for(x=0; x<100000; x++) cdslen[x] = 0;

        strncpy(rstr,empty,122);
        while((fgets(rstr,121,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%1000==0)	printf("Now processing variant %d\n",counter);
		}
		else if(strncmp(rstr,"GeneName",8)==0){ /* Get gname: ENSRNOG00000047964 for instance */
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
			if(gmax<0 || length<0)          goto reinit;
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
			else if(exon[exmax][2]<=0)    goto reinit;    /* Neglecting single-exon genes */

			/* Determine the max exon length bin (1~9) */
			w = 0;	/* Max length of internal exons in nt */
			binmax = -1;
			for(x=0; x<=exmax; x++){
				if(exon[x][2]!=2)	continue;
				else if(exon[x][3]>w)	w = exon[x][3];
			}
			if(w<=0)	goto reinit;		/* No internal exons */
			++hitcount;

			/* ### Dermine the protein length group */
			gpnumb = 4;	/* Default gp number */
			for(x=1; x<=3; x++){
				if(length<=gpmax[x]){
					gpnumb = x;
					break;
				}
			}

			binmax = (w-1)/120 +1;
			if(binmax>=10)	binmax = 9;
			z = scmax + 1;	/* The number of SCOP domains */
			/* ### */
			sumscop0[gpnumb][binmax] += 1.; 
			sumscop1[gpnumb][binmax] += (double)z; 
			sumscop2[gpnumb][binmax] += (double)z *(double)z; 

			sumn[gpnumb] += 1.;
			sumx1[gpnumb] += (double)w/scaling;
			sumx2[gpnumb] += ((double)w/scaling) *((double)w/scaling);
			sumy1[gpnumb] += (double)z/scaling;
			sumy2[gpnumb] += ((double)z/scaling) *((double)z/scaling);
			sumxy[gpnumb] += ((double)w/scaling) * ((double)z/scaling);

                        sumlen1[gpnumb][binmax] += (double)length/scaling;
                        sumlen2[gpnumb][binmax] += ((double)length/scaling)*((double)length/scaling);
                        ++cmax;
                        if(cmax>=100000){       printf("The number of proteins exceeded 100,000\n");    exit(1);        }
                        cdslen[cmax] = length;

reinit:			/* Reinitializations */
			length = 0;
			for(x=0; x<100; x++)    strncpy(gname[x],empty,21);
			gmax = -1;
			scmax = -1;
			for(x=0; x<10000; x++){
				for(y=0; y<2; y++)	scop[x][y] = -1;
			}
		}
nextr:		strncpy(rstr,empty,122);
        }
        /* Sort the length array and output 1/4, 1/2, and 3/4 boundaries */
        qsort(cdslen,cmax+1,sizeof(int),intcmp);
        x = (cmax -3)/4;        /* 1/4 boundary */
        y = (cmax -1)/2;        /* 1/2 boundary */
        z = (3*cmax -1)/4;      /* 3/4 boundary */
	fputs("CDS len 1/4 boundary\t",fw);
        strncpy(tempstr,empty,21);      sprintf(tempstr,"%d",cdslen[x]);        fputs(tempstr,fw);      fputs("\t",fw);
        fputs("1/2 boundary\t",fw);
        strncpy(tempstr,empty,21);      sprintf(tempstr,"%d",cdslen[y]);        fputs(tempstr,fw);      fputs("\t",fw);
        fputs("3/4 boundary\t",fw);
        strncpy(tempstr,empty,21);      sprintf(tempstr,"%d",cdslen[z]);        fputs(tempstr,fw);      fputs("\n",fw);

	/* ### Sum over all max exon length bins */
	for(y=1; y<5; y++){
		for(z=1; z<10; z++){
			sumscop0[y][0] += sumscop0[y][z];
			sumscop1[y][0] += sumscop1[y][z];
			sumscop2[y][0] += sumscop2[y][z];
                	sumlen1[y][0] += sumlen1[y][z];
                	sumlen2[y][0] += sumlen2[y][z];
		}
	}
	/* ### Sum over all the protein length groups */
	for(y=1; y<5; y++){
		sumn[0] += sumn[y];	sumx1[0] += sumx1[y];	sumx2[0] += sumx2[y];
		sumy1[0] += sumy1[y];	sumy2[0] += sumy2[y];	sumxy[0] += sumxy[y];
		for(z=0; z<10; z++){
			sumscop0[0][z] += sumscop0[y][z];
			sumscop1[0][z] += sumscop1[y][z];
			sumscop2[0][z] += sumscop2[y][z];
                	sumlen1[0][z] += sumlen1[y][z];
                	sumlen2[0][z] += sumlen2[y][z];
		}
	}

	/* Correlation of max int. exon len and the # of SCOP domains */
	for(gpnumb=0; gpnumb<5; gpnumb++){
		fputs("Prot len gp: ",fw);
		fputs(gpstr[gpnumb],fw);	fputs("\n",fw);
		gt = sumx2[gpnumb] - (sumx1[gpnumb]*sumx1[gpnumb]/sumn[gpnumb]);
		gu = sumy2[gpnumb] - (sumy1[gpnumb]*sumy1[gpnumb]/sumn[gpnumb]);
		if(gt<0.){	printf("gt=%f\n",gt);	exit(1);	}
		if(gu<0.){	printf("gu=%f\n",gu);	exit(1);	}
		corr = (sumxy[gpnumb] - (sumx1[gpnumb]*sumy1[gpnumb]/sumn[gpnumb]))/sqrt(gt * gu);
		core = (1. - corr*corr)/sqrt(sumn[gpnumb]);
		cort = corr * sqrt(sumn[gpnumb] -2.)/(1. - corr*corr);

		fputs("r =\t",fw);
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.4f",corr);	fputs(tempstr,fw);
		fputs("\t+/-\t",fw);
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.4f",core);	fputs(tempstr,fw);
		fputs("\tt=\t",fw);
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",cort);	fputs(tempstr,fw);
		fputs("\n",fw);

		fputs("Max int. ex len.\t#\tAvg SCOP#\tSEM\tAvg len\tSEM\n",fw);
		for(z=0; z<10; z++){
			if(z==0)	fputs("All\t",fw);
			else{
				strncpy(tempstr,empty,21);
				sprintf(tempstr,"%d",120*z-119);
				fputs(tempstr,fw);
				if(z==9)	fputs(" bp~\t",fw);
				else{ 
					fputs("~",fw);
					strncpy(tempstr,empty,21);
					sprintf(tempstr,"%d",120*z);
					fputs(tempstr,fw);
					fputs(" bp\t",fw);
				}
			}
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",sumscop0[gpnumb][z]);
			fputs(tempstr,fw);		fputs("\t",fw);

			if(sumscop0[gpnumb][z]<0.1){
				fputs("\n",fw);
				continue;
			}
			strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",sumscop1[gpnumb][z]/sumscop0[gpnumb][z]);
			fputs(tempstr,fw);		fputs("\t",fw);

			gx = sumscop2[gpnumb][z] - (sumscop1[gpnumb][z]*sumscop1[gpnumb][z]/sumscop0[gpnumb][z]);
			if(gx<0.){	printf("gx = %f at #100\n",gx);		exit(1);	}
			gy = sqrt(gx)/sumscop0[gpnumb][z];		/* SEM */
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gy);
                	fputs(tempstr,fw);              fputs("\t",fw);

                	strncpy(tempstr,empty,21);      sprintf(tempstr,"%.2f",scaling*sumlen1[gpnumb][z]/sumscop0[gpnumb][z]);
                	fputs(tempstr,fw);              fputs("\t",fw);

                	gx = sumlen2[gpnumb][z] - (sumlen1[gpnumb][z]*sumlen1[gpnumb][z]/sumscop0[gpnumb][z]);
                	if(gx<0.){      printf("gx = %f at #200\n",gx);         exit(1);        }
                	gy = scaling*sqrt(gx)/sumscop0[gpnumb][z];              /* SEM */
                	strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gy);
                	fputs(tempstr,fw);              fputs("\n",fw);
		}
		fputs("\n",fw);
	}
	printf("%d out of %d variants were analyzed\n",hitcount,counter);
        fclose(fr);	fclose(fl);        fclose(fw);
	return 0;
}
