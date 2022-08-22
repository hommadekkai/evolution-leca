/* mrinsexon1G.c		This program divides mouse-specific segments into exon classes 			220517
	basic.c -> minsexon1.c: A new program [R2105 II].
	-> minsexon1A.c: Modified to count IntEx+MultiEX, classify into 10 instead of 20 a.a. bins, and find the avg. # of exons.
	-> minsexon1B.c: Modified to classify into 40 a.a. bins and to output IntEx or MultiEx cases with >=201 a.a. as minsexon1B.
	-> minsexon1C.c: Modified to deifine MultiEx as multiple INTERNAL exons.
	-> minsexon1D.c: Modified to EXCLUDE insertion segments whose both ends coincide with exon boundaries +/- 5 a.a.
			and to classify into 20 a.a. bins.
	-> mrinsexon1E.c: Modified to make use of the exon list by a new exon definition (#&#) [R2205 I].
	-> mrinsexon1F.c: Modified to classify ins. to internal exons (IntExNC) into 5'-end ins (IntEx5), 3'-end ins (IntEx3),
			and ins to inside (InsExIn) [R2205 I].
	-> mrinsexon1G.c: Modified to change the allowance of exon coincidence from 5 to 1 a.a. (#+#) [R2205 I].    */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fe, *fw;
char *ps;
char fnamer[122]="/Users/khomma/work/ensembl/hmr/minsseg1";	/* List of mouse-specific regions */
char fnamee[122]="/Users/khomma/work/ensembl/mouse/mallexon3";	/* #&# List of all mouse exons */
char fnamew[122]="/Users/khomma/work/ensembl/hmr/mrinsexon1G";	/* #+# List of long insertions */
char rstr[1001], sstr[5001], empty[5010]="", tempstr[21], trid[21];
char tchr[10], tgname[51], oldgname[51], ttname[51], oldtname[51], tsense[2];
char clsstr[10][20]={"All","IntEx5","IntEx3","IntExIn","MultiExNC","IntExCI","MultiExCI","5tEx","3tEx","SingEx"}; 
	   /* #&# NC: No coninc, CI: coinc */
char enspid[21];		
int  v, w, x, y, z, offset, counter=0;
int  inscutoff[9]={1,11,21,31,41,51,61,71,81};   /* #&# The min length of inserted segment */
int  tcdslen, oldcdslen, tstart, tend, texnumb, tresstart, tresend;
int  insseg[1000][2], insmax, tlen, hitflag, varmax, hitexon[10000], nexon, hitcount;
int  class[9][10];	/* ##### The number of cases in each class (clsstr) with the five diff. inscutoff values */
int  cclass;		/* Which inscutoff values does the ins. length qualify: 0~12 */
int  firstex, lastex, numbex0=0, numbex1=0, numbex2=0, totres[9], excoinc[3], exlen;	/* #&# */
int  allowance=1;		/* #+# */
double gw, gx, gy, gz;
void segread(void){    /* Parse the mouse-specific segment string get the start and end of each segment */
	strncpy(trid,empty,21);
	for(x=0; x<1000; x++){
		for(y=0; y<2; y++)	insseg[x][y] = 0;	/* y=0: start, y=1: end */
	}
	insmax = -1;
	tlen = 0;	

	strncpy(enspid,empty,21);
	for(x=0; x<20; x++){
		if(isspace(sstr[x])!=0)	break;
		enspid[x] = sstr[x];

	}

	ps = strstr(sstr,"ENSMUST");
	for(x=0; x<20; x++){
		if(isspace(*ps)!=0)	break;
		trid[x] = *ps;
		++ps;
	}

	strncpy(tempstr,empty,21);
	ps = strstr(sstr,"LEN=");
	ps = ps + 4;
	for(x=0; x<20; x++){
		if(isdigit(*ps)==0)	break;
		tempstr[x] = *ps;
		++ps;
	}
	tlen = atoi(tempstr);
	
	for(x=0; x<100; x++)	if(sstr[x]==']')	break;
	offset = x + 2;
	while(1){
		++insmax;
		if(insmax>=1000){	printf("The # of insertion segments exceded 1000 in %s",sstr);	exit(1);	}
		strncpy(tempstr,empty,21);
		for(x=0; x<20; x++){
			if(isdigit(sstr[x+offset])==0)		break;
			tempstr[x] = sstr[x+offset];
		}
		offset = x + offset + 2;

		insseg[insmax][0] = atoi(tempstr);
		strncpy(tempstr,empty,21);
		for(x=0; x<20; x++){
			if(isdigit(sstr[x+offset])==0)		break;
			tempstr[x] = sstr[x+offset];
		}
		offset = x + offset + 1;
		insseg[insmax][1] = atoi(tempstr);
		if(isprint(sstr[offset])==0)	break;		/* End of line */
	}
	return;
}
struct variant{
        char chr[51];           /* Chromosome name */
        char gname[51];         /* Gene name */
        char tname[51];         /* Transcript name */
        int  cdslen;
        int  start;             /* Gene start address */
        int  end;               /* Gene end address */
        char sense[2];          /* + or - */
        int  exnumb;            /* Exon number */
	int  totex;             /* #&# Total exon number */
        int  resstart;          /* The start res # of the exon */
        int  resend;            /* The end res # of the exon */
};
struct variant vlist[10000];
void tread(void){    /* Parse the line to get the start and end of the inner exon */
        int catcount=0;         /* Category vhit */
        offset = 0;
        strncpy(vlist[varmax].chr,empty,10);         	strncpy(vlist[varmax].gname,empty,51);
       	strncpy(vlist[varmax].tname,empty,51);		strncpy(vlist[varmax].sense,empty,2);
        vlist[varmax].start = 0;     	vlist[varmax].end = 0;       	vlist[varmax].exnumb=0;      
	vlist[varmax].totex = 0;        vlist[varmax].resstart = 0;     vlist[varmax].resend = 0;       /* #### */
        while(1){
                ++catcount;
                strncpy(tempstr,empty,21);
                for(x=0; x<50; x++){
                        if(isspace(rstr[x+offset])!=0)          break;
                        tempstr[x] = rstr[x+offset];
                }
                offset = x + offset;
                if(catcount==1)         strcpy(vlist[varmax].chr,tempstr);              /* Chromosome name */
                else if(catcount==2)    strcpy(vlist[varmax].gname,tempstr);            /* Gene name */
                else if(catcount==3)    strcpy(vlist[varmax].tname,tempstr);            /* Transcript name */
                else if(catcount==4)    vlist[varmax].cdslen = atoi(tempstr);           /* CDS length */
                else if(catcount==5)    vlist[varmax].start = atoi(tempstr);            /* Start */
                else if(catcount==6)    vlist[varmax].end = atoi(tempstr);              /* End */
                else if(catcount==7)    strcpy(vlist[varmax].sense,tempstr);            /* Sense */
                else if(catcount==8)    vlist[varmax].exnumb = atoi(tempstr);           /* Exon number */
                else if(catcount==9)    vlist[varmax].totex = atoi(tempstr);            /* Total exon number */
                else if(catcount==10)    vlist[varmax].resstart = atoi(tempstr);        /* The start res. of the exon */
                else if(catcount==11){     
                        vlist[varmax].resend = atoi(tempstr);                           /* The end res. of the exon */
                        return;
                }

                for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;
                offset = x + offset;
        }
        return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fe=fopen(fnamee,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamee);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs("HumanID\tID\tLength(aa)\tClass\tSeg.start\tSeg.end\n",fw);
	for(y=0; y<9; y++){	
		totres[y] = 0;		/* Total # of ins. residues of IntEx or MultEx in each range */
		for(x=0; x<10; x++)	class[y][x] = 0; 
	}
	strncpy(sstr,empty,5001);
        while((fgets(sstr,5000,fr))!=NULL){
		if(strstr(sstr,"ENSMUST")==NULL)	goto nexts;
		++counter;
		if(counter%500==0)	printf("Now processing line %d...\n",counter);
		segread();

		/* Look for the exon assignments */
		fseek(fe,0L,SEEK_SET);
		strncpy(rstr,empty,1001);
		hitflag = 0;
		varmax = -1;
        	while((fgets(rstr,1000,fe))!=NULL){
			if(strstr(rstr,trid)!=NULL && strstr(rstr,"ERROR")==NULL && strstr(rstr,"Noncoding")==NULL){    /* ###### */
				hitflag = 1;
				++varmax;
				if(varmax>=10000){
					printf("The number of exons in %s exceeded 10,000!\n",trid);
					exit(1);
				}
				tread();
                                if(vlist[varmax].cdslen!=tlen && vlist[varmax].cdslen!=tlen-1){
                                        /* Allowing one residue less than the registered protein length */
					printf("Lengths of %s disagreed: %d vs. %d\n",trid,vlist[varmax].cdslen,tlen);
					goto nexts;	/* Lengths disagreed */
				}
			}
			else if(hitflag==1){
				/* The end of the exon data for this transcript */
				goto aexon;
			}
			strncpy(rstr,empty,1001);
		}
		if(hitflag==0)	goto nexts;

aexon:		/* Analyze the exon info */
		for(x=0; x<=insmax; x++){	
			exlen = insseg[x][1] - insseg[x][0] + 1;

			cclass = 8;	
			for(y=0; y<=7; y++){
				if(exlen<inscutoff[y+1]){
					cclass = y;
					break;
				}
			}
	
			for(y=0; y<10000; y++)	hitexon[y] = 0;	  
			/* 1: One end of the inserted segment exists in exon y, 2: The segment is entirely in exon y  */
			hitcount = 0;
			/* #&# Count the # of internal exon boundaries that coincide with the insertion segment +/- 1 a.a. #+# */
			for(z=0; z<3; z++)	excoinc[z] = 0;	
				/* Coincidence (1) or noconincidence(0) at 5'-end (z=1), 3'-end (z=2), and both ends (z=0) */
			/* ### Debugged ### */
			if(varmax>=0){
				for(y=0; y<=varmax; y++){
					if(vlist[y].exnumb!=1){	/* Neither the single exon nor the 5'-term exon */
						if(abs(insseg[x][0]-vlist[y].resstart)<=allowance){	/* #+# */
						 	++excoinc[1];	/* The 5' end agrees */
							++excoinc[0];
						}
					}
					if(vlist[y].exnumb!=vlist[y].totex){	/* Neither the single exon nor the 3'-term exon */
						if(abs(insseg[x][1]-vlist[y].resend)<=allowance){	/* #+# */
							++excoinc[2];	/* The 3' end agrees */
							++excoinc[0];
						}
					}
				}
			}
			/* #&# */

			for(y=0; y<=varmax; y++){
				if(vlist[y].resstart<=insseg[x][0] && insseg[x][0]<=vlist[y].resend){
					/* The start of the inserted segment is in exon y */
					++hitexon[y];
					++hitcount;
				}
                                else if(y==varmax && insseg[x][0]==vlist[y].resend+1){
                                        /* Allowing one residue less than the registered protein length */
                                        ++hitexon[y];
                                        ++hitcount;
                                }

				if(vlist[y].resstart<=insseg[x][1] && insseg[x][1]<=vlist[y].resend){
					/* The end of the inserted segment is in exon y */
					++hitexon[y];
					++hitcount;
				}
                                else if(y==varmax && insseg[x][1]==vlist[y].resend+1){
                                        /* Allowing one residue less than the registered protein length */
                                        ++hitexon[y];
                                        ++hitcount;
                                }
			}
			if(hitcount!=2){
			 	printf("hitcount=%d for the segnment %d-%d in %s\n",hitcount,insseg[x][0],insseg[x][1],trid);
				goto nextx;	/* An abnormal case */
			}

			/* Determnie the # of exons encoding the ins. segement */
			firstex = -1; 	lastex = -1;
			for(y=0; y<=varmax; y++){
				if(hitexon[y]==2){
					firstex = y;
					lastex = y;
				}
				else if(hitexon[y]==1){
					firstex = y;
					break;
				}
			}
			if(firstex>=0 && lastex<0){
				for(y=firstex+1; y<=varmax; y++){
					if(hitexon[y]==1){
						lastex = y;
						break;
					}
				}
			}

			nexon = 0;	
			if(vlist[0].totex==1){		/* #&# */	
				if(hitexon[0]==2)		nexon = 9;	/* #&# Single exon */
			}
			else{
				/* Check for the existence of an exon in which the inserted segment x wholly resides */
				for(y=0; y<=varmax; y++){
					if(hitexon[y]==2){	/* One exon */
						if(vlist[y].exnumb==1)  	nexon = 7;      /* #&# 5'-term exon */
						else if(vlist[y].exnumb==vlist[y].totex) 	nexon = 8;      /* #&# 3'-term exon */
						else if(excoinc[0]<=1){		/* Internal exon w/o coinc. */
							if(excoinc[1]>=1)	nexon = 1;	/* 5' ins in an int exon  */
							else if(excoinc[2]>=1)	nexon = 2;	/* 3' ins in an int exon */
							else			nexon = 3;	/* Ins inside an int exon */
						}
						else				nexon = 5;	/* #&# An internal exon w/ coinc. */
						break;	
					}
				}
				if(nexon==0){	/* ### Debugged. Encoded by multiple exons */
					if(excoinc[0]<=1)		nexon = 4;	/* ##### Multiple internal exons w/o coinc. */
					else				nexon = 6;	/* ##### Multiple internal exons w/ coinc. */
				}
			}

			/* #&# Limit the exon # stat. to IntExNC and MultExNC */
			if(firstex>=0 && lastex>=0 && lastex>=firstex && nexon<=4){	/* #&# */
				v = lastex - firstex + 1;	/* The number of exons that correspond to the inserted segment */
				++numbex0;
				numbex1 += v;
				numbex2 += v * v;
			}
	
			if(nexon<=4)	totres[cclass] += exlen;	/* exlen is the len. of the ins. seg. */
			
			/* Output the results */
			if(nexon==0){		
			 	printf("No exon hits for the segnment %d-%d in %s\n",insseg[x][0],insseg[x][1],trid);
				goto nextx;
			}
			else{
				++class[cclass][nexon];
				if(nexon<=4){
					/* Both ends of the ins. seg DO NOT coincide with exon boundaries */
					fputs(enspid,fw);		fputs("\t",fw);
					fputs(trid,fw);			fputs("\t",fw);
					strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",tlen);
					fputs(tempstr,fw);		fputs("\t",fw);
					fputs(clsstr[nexon],fw);	fputs("\t",fw);
					strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",insseg[x][0]);
					fputs(tempstr,fw);		fputs("\t",fw);
					strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",insseg[x][1]);
					fputs(tempstr,fw);		fputs("\n",fw);
				}
			}
nextx:			x = x;
               	} 
nexts:		strncpy(sstr,empty,5001);
        }
	z = 0;
	for(x=0; x<9; x++)	z += totres[x];
	printf("The total res. # of inserted segments is %d\n",z);

	gx = (double)numbex1/(double)numbex0;	/* Avg. # of exons */
	gy = (double)numbex2 - (double)(numbex1)*(numbex1)/(double)numbex0;
	if(gy<-0.0000001){	printf("gy=%f\n",gy);	exit(1);	}
	gz = sqrt(gy)/(double)numbex0;	/* SEM */
	printf("# of exons encoding ins. segments without exon boundary coincidence =%.3f +/- %.3f (n=%d)\n",gx,gz,numbex0);

	for(x=0; x<9; x++){
		for(y=1; y<10; y++)	class[x][0] += class[x][y];
	}
        printf("Class\t");
        for(x=0; x<9; x++){
                printf("#%d",x+1);
                if(x==8)        printf("\n");
                else            printf("\t");
        }

	printf("Ins seg\t");
        for(x=0; x<9; x++){
                printf("%daa~",inscutoff[x]);
                if(x==8)        printf("\n");
                else            printf("%d\t",inscutoff[x+1]-1);
        }

	/* The # of residues in each lenth range */
	printf("# of res.\t");
        for(x=0; x<9; x++){
		printf("%d",totres[x]);
		if(x==8)	printf("\n");
		else		printf("\t");
	}

	/* # of segments */
	for(y=0; y<10; y++){
		printf("%s\t",clsstr[y]);	
		for(x=0; x<9; x++){
			printf("%d",class[x][y]);
			if(x==8)	printf("\n");
			else		printf("\t");
		}
	}
        fclose(fe);	fclose(fr);        fclose(fw);
	return 0;
}
