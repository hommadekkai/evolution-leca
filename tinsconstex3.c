#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fi, *fr, *fw, *fl;
char rstr[122], empty[500]="",  tempstr[21], tempexon[51];
char gname[100][21], lstr[201];
char fnamei[122]="/Users/khomma/work/ensembl/hmr/trinsexon1G";          /* #-# Classified list of inserted segments */
char fnamer[122]="/Users/khomma/work/ensembl/rat/nscoplist1";
char fnamel[100]="/Users/khomma/work/ensembl/rat/nalluniqexon3M";       /* Non-redundant list of exons, with the CDS range by int exons */
char fnamew[122]="/Users/khomma/work/ensembl/hmr/tinsconstex3";         /* #-# */
char insid[5000][21], trid[21];		/* #### */
char tpstr[4][20]={"All","5'-term","3'-term","Inside"};
char binstr[10][20]={"All","-120","-240","-360","-480","-600","-720","-840","-960","961-nt"};   /* #-# */
int t, u, v, w, x, y, z, offset, length, scmax=-1, exmax, hitflag, counter=0, hitcount=0, gmax;
int scop[10000][2];	/* 0: Start, 1: End of SCOP alignment */
int exon[10000][4];	/* 0: Start, 1: End of the exon-encoded residues, 2: Classification of each exon
			Classification =1: 5'-term exon, =2: interna; exon, =3: 3'-term exon, =0: Single exon
			constrained exon=1, otherwise 0 */
int excdslen, type;
int insseg[5000][3], insmax=-1;		/* 3: type of insertion segment (1,2,or3) */
int inscut[3]={0,5,10};	/* Ins. length cutoffs in a.a. residue */
int insstart, insend, exbin;	/* #-# */
double gexp[10]={31.298,34.291,29.685,22.073,17.793,15.575,17.303,19.379,19.721,15.690};	/* From R2204 nscopconstex3 */
double gv, gw, gx, gy, gz, scaling=1000.;	/* #-# */
double nobs[3][10][4], nins[3][10][4], gobs, sumint[10], sumconst[10];	/* #-# */ 
	/* nobs: # of ins segments in constraied exons, nins: # of ins seg segments 
	  10 length bins
	  3 ins length cutoffs; 3 +1(all) types of insertions are differentiated */
void insparse(void)
{
	int i, j;
	offset = 0;
	++insmax;
	if(insmax>=5000){ printf("The # of inserted segments exceeded 5,000!\n");	exit(1);	}	/* #### */
	insseg[insmax][2] = type;	/* Registering the type of the insertion seg */
	while(1){
		++j;
		strncpy(tempstr,empty,21);
		for(i=0; i<20; i++){
			if(isspace(rstr[i+offset])!=0)	break;
			else	tempstr[i] = rstr[i+offset];
		}
		offset = i + offset;
		if(j==2)	strcpy(insid[insmax],tempstr);		/* Transcript ID */
		else if(j==5)	insseg[insmax][0] = atoi(tempstr);
		else if(j==6){
			insseg[insmax][1] = atoi(tempstr);
			return;
		}
		for(i=0; i<20; i++)	if(isspace(rstr[i+offset])==0)  break;
		offset = i + offset;
	}
	return;
}
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
void exonparse(int k)
{
	int i, j=0, l;
	offset = 0;
	while(1){
		++j;
		strncpy(tempexon,empty,51);
		for(i=0; i<=50; i++){
			if(isspace(lstr[i+offset])!=0)	break;
			else	tempexon[i] = lstr[i+offset];
		}
		offset = i + offset;
		if(j==k)	return;
		else if(j==13)	return;

		for(i=0; i<50; i++)	if(isspace(lstr[i+offset])==0)  break;
		offset = i + offset;
	}
	return;
}
int main(void)
{
        if((fi=fopen(fnamei,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamei);	exit(1);	}
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fl=fopen(fnamel,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamel);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(z=0; z<3; z++){	/* 3 insertion length cutoffs */
		for(u=0; u<10; u++){    /* #-# Exon length bins */
			for(t=0; t<4; t++){	/* Types of ins seg */
				nobs[z][u][t] = 0.;
				nins[z][u][t] = 0.;
			}
		}
	}
	for(x=0; x<5000; x++){		/* #### */
		strncpy(insid[x],empty,21);
		for(y=0; y<3; y++)	insseg[x][y] = 0.;	/* y=0: start, y=1: end, y=2: type of ins seg */
	}
	for(u=0; u<10; u++){    /* #-# Exon length bins */
		sumint[u] = 0.;
		sumconst[u] = 0.;
	}

	/* Read in the classified list of inserted segments */
        strncpy(rstr,empty,122);
        while((fgets(rstr,121,fi))!=NULL){
		if(strstr(rstr,"IntEx5")!=NULL){
			type = 1;
			insparse();
		}
		else if(strstr(rstr,"IntEx3")!=NULL){
			type = 2;
			insparse();
		}
		else if(strstr(rstr,"IntExIn")!=NULL){
			type = 3;
			insparse();
		}
        	strncpy(rstr,empty,122);
	}
	fclose(fi);
	printf("There are %d insertion segments inside internal exons\n",insmax+1);	/* #### */

	/* Read in the SCOP alignment list */
        strncpy(rstr,empty,122);
        while((fgets(rstr,121,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%1000==0)	printf("Now processing variant %d\n",counter);
		}
		else if(strncmp(rstr,"GeneName",8)==0){ /* Get gname: ENSRNOG00000091476 for instance */
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
				for(y=0; y<4; y++)      exon[x][y] = 0; 
			}
			fseek(fl,0L,SEEK_SET);
			strncpy(lstr,empty,201);
                        while((fgets(lstr,200,fl))!=NULL){
                                if(hitflag==0){
                                        for(x=0; x<=gmax; x++){
                                                if(strstr(lstr,gname[x])!=NULL){        /* Possible hit. Check the total length */
							exonparse(4);
                                                        if(length==atoi(tempexon)){  /* HIT! */
                                                                hitflag = 1;
                                                                ++exmax;
                                                                if(exmax>=10000){
                                                                        printf("Exon # of %s(length=%d)exceeded 10,000\n",gname[x],length);
                                                                        exit(1);
                                                                }
								exonparse(3);
								strncpy(trid,empty,21);
								strcpy(trid,tempexon);	/* Transcript ID */
								exonparse(8);
								t = atoi(tempexon);	/* Exon number */
								exonparse(9);
								u = atoi(tempexon);     /* Total exon number */
								if(u==1)        exon[exmax][2] = 0;     /* Single exon */
								else if(t==1)   exon[exmax][2] = 1;     /* 5'-term exon */
								else if(t==u)   exon[exmax][2] = 3;     /* 3'-term exon */
								else            exon[exmax][2] = 2;     /* Internal exon */

								exonparse(10);
								exon[exmax][0] = atoi(tempexon);         /* Exon res start */
								exonparse(11);
								exon[exmax][1] = atoi(tempexon);         /* Exon res end */
								exonparse(12);
								insstart = atoi(tempexon);
								exonparse(13);
								insend = atoi(tempexon); 
                                                        }
                                                }
                                        }
                                }
                                else if(hitflag==1){
                                        hitflag = -1;
                                        for(x=0; x<=gmax; x++){
                                                if(strstr(lstr,gname[x])!=NULL){        /* Possible hit. Check the total length */
							exonparse(4);
                                                        if(length==atoi(tempexon)){  /* HIT! */
                                                                hitflag = 1;
                                                                ++exmax;
                                                                if(exmax>=10000){
                                                                        printf("Exon # of %s(length=%d)exceeded 10,000\n",gname[x],length);
                                                                        exit(1);
                                                                }
								exonparse(3);
								strncpy(trid,empty,21);
								strcpy(trid,tempexon);	/* Transcript ID */
								exonparse(8);
								t = atoi(tempexon);     /* Exon number */
								exonparse(9);
								u = atoi(tempexon);     /* Total exon number */
								if(u==1)        exon[exmax][2] = 0;     /* Single exon */
								else if(t==1)   exon[exmax][2] = 1;     /* 5'-term exon */
								else if(t==u)   exon[exmax][2] = 3;     /* 3'-term exon */
								else            exon[exmax][2] = 2;     /* Internal exon */
								exonparse(10);
								exon[exmax][0] = atoi(tempexon);         /* Exon res start */
								exonparse(11);
								exon[exmax][1] = atoi(tempexon);         /* Exon res end */
								exonparse(12);
								insstart = atoi(tempexon);
								exonparse(13);
								insend = atoi(tempexon); 
                                                        }
                                                }
                                        }
                                        if(hitflag==-1) break;
                                        /* The end of a gene */
                                }

                                strncpy(lstr,empty,201);
                        }
			if(hitflag==0 || exmax<0)	goto reinit;

			/* Examine the observed # of exon boundaries that are inside SCOP domains and the expected #s */
			if(exon[0][2]<=0)    goto reinit;    /* Neglecting single-exon genes */
			++hitcount;

			/* First, identify constrained exons and set exon[x][3] = 1 */
			/* #+# DO NOT NEGLECT genes without SCOP domains */
			for(x=0; x<=exmax; x++){
				v = exon[x][1] - exon[x][0] +1;
				exbin = ((v-1)/120) + 1;        /* #-#  Exon length bin */
				if(exbin>9) exbin = 9;		/* #-# */
			
				sumint[exbin] += (double)w;
				if(scmax<0)	continue;	/* #+# */
				for(y=0; y<=scmax; y++){
					if(scop[y][0]<exon[x][0] && exon[x][1]<scop[y][1]){	/* Constrained exon */
						exon[x][3] = 1;
						sumconst[exbin] += (double)w;
						break;
					}
				}
			}

			/* Look for insertion segment(s) */
			for(y=0; y<=insmax; y++){	/* One inssegment at a time */
				for(z=0; z<3; z++){
					if(insseg[y][1]-insseg[y][0]+1<inscut[z])	continue;
					if(strcmp(trid,insid[y])==0){	/* FOUND! */
						type = insseg[y][2];
						/* #+# Find if the residues in the ins segment is inside a constrained exon or not */
						for(w=insseg[y][0]; w<=insseg[y][1]; w++){	/* One residue at a time */
							for(x=0; x<=exmax; x++){
								if(exon[x][2]!=2)	continue;	/* Only internal exons accepted */
								else if(exon[x][0]<=w && w<=exon[x][1]){	/* Hit */
									v = exon[x][1] - exon[x][0] +1;
									exbin = ((v-1)/120) + 1;        /* #-#  Exon length bin */
									if(exbin>9) exbin = 9;		/* #-# */
									nins[z][exbin][type] += 1.;	/* #-# */
									if(exon[x][3]==1)	nobs[z][exbin][type] += 1.;	/* #-# */
										/* Hit in a constrained exon! */
									break;
								}
							}
						}
					}
				}
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
	/* Sum over all the types of ins seg */
	for(t=1; t<4; t++){
		for(z=0; z<3; z++){	/* 3 insertion length cutoffs */
			for(u=1; u<10; u++){    /* #-# Exon length bins */
				nins[z][u][0] += nins[z][u][t];
				nobs[z][u][0] += nobs[z][u][t];
			}
		}
	}

	/* Sum over all exon length bins #-# */
	for(z=0; z<3; z++){	/* 3 insertion length cutoffs */
		for(u=1; u<10; u++){    /* #-# Exon length bins */
			for(t=0; t<4; t++){
				nins[z][0][t] += nins[z][u][t];
				nobs[z][0][t] += nobs[z][u][t];
			}
		}
	}
	for(u=1; u<10; u++){    /* #-# Exon length bins */
		sumconst[0] += sumconst[u];
		sumint[0] += sumint[u];
	}

	/* All exon ranges */
	fputs("Ins seg in constrained exons\n",fw);
	for(z=0; z<3; z++){
		fputs("Insertion segment cutoff=",fw);
		strncpy(tempstr,empty,21);
		sprintf(tempstr,"%d",inscut[z]);
		fputs(tempstr,fw);	fputs(" aa\n",fw);

		fputs("Ins type\tExon len\tNumber\tObs.freq.(%)\tExp.freq.(%)\tChi sq.\tObs/Exp\n",fw);		/* #-# */
		for(t=0; t<4; t++){	/* Types of ins seg */
			for(u=0; u<10; u++){    /* #-# Exon length bins */
				if(u==0)	fputs(tpstr[t],fw);	/* #-# */	
				fputs("\t",fw);
				fputs(binstr[u],fw);	fputs("\t",fw);		/* #-# */

				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",nins[z][u][t]);	/* #-# */
				fputs(tempstr,fw);	fputs("\t",fw);

				if(nins[z][u][t]<=0.1){		/* #-# */
					fputs("\n",fw);
					continue;
				}
				gobs = nobs[z][u][t]/nins[z][u][t];	/* #-# Observed fraction */
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",100.*gobs);
				fputs(tempstr,fw);		fputs("\t",fw);

				if(sumint[u]<=0){
					fputs("\n",fw);
					continue;
				}
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gexp[u]);
				fputs(tempstr,fw);		fputs("\t",fw);

				/* Chi-sq */
				if(nins[z][u][t]<=0.1){		/* #-# */	
					fputs("\n",fw);
					continue;
				} 
				gx = nins[z][u][t] * gexp[u]/100.;	/* #-# */
				if(gx<0.000001){	/* #-# */
					fputs("\n",fw);
					continue;
				}
				gy = (nobs[z][u][t] -gx)*(nobs[z][u][t] -gx)/gx;	/* #-# */

				gv = nins[z][u][t]*nins[z][u][t]*(1.-gexp[u]/100.)*(1.-gexp[u]/100.);	/* #-# */
				if(abs(gv)<0.000001){	/* #-# */
					fputs("\n",fw);
					continue;
				}
				gz = (nobs[z][u][t] -gx)*(nobs[z][u][t] -gx)/gv;   /* #-# */
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gy+gz);
				fputs(tempstr,fw);	fputs("\t",fw);

				/* Obs/Exp */
				gw = nobs[z][u][t]/gx;		/* #-# */
				strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gw);
				fputs(tempstr,fw);	fputs("\n",fw);
			}
		}
	}
	printf("%d out of %d variants were analyzed\n",hitcount,counter);
        fclose(fr);	fclose(fl);        fclose(fw);
	return 0;
}
