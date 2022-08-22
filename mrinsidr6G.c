/* mrinsidr6G.c	This program determines the %IDR in inserted segments and that in other regions 	220519
	basic.c -> minsidr4.c: A new program [R2105 III].
	-> mainsidr4.c: Modified to include MultiEx as well as IntEx for analyses (#&#: modif. line)[R2105 III].
	-> mrinsidr4MA.c: Modified to use an internal exon list by new def. [R2205 II].
	-> mrinsidr6MA.c: Modified to make use of DISOPRED3 instead of POODLE-L results [R2205 II].
	-> mrinsidr6G.c: Modified to use a new list of inserted segments [R2205 II].  	  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fc, *fe, *fi;
char rstr[5001], empty[5010]="", cstr[122], estr[501], istr[501], tempstr[51], trid[51], protid[51];
char fnamer[122]="/Users/khomma/work/disopred/mdisopr1";        /* PARTIAL DISOPRED3 results */
char fnamec[122]="/Users/khomma/work/ensembl/mouse/mptcorr1";	/* Correspondence between protein(ENSMUSP) and transcr.(ENSMUST) IDs */
char fnamee[122]="/Users/khomma/work/ensembl/mouse/mallexon3";	/* Revised exon list */
char fnamei[122]="/Users/khomma/work/ensembl/hmr/mrinsexon1G";	/* NEW list of inserted segments */
char fnamew[122]="/Users/khomma/work/ensembl/hmr/mrinsidr6G";
char instr[10000][51];
int w, x, y, z, offset, plen, imax, hitflag, counter=0, insflag, sflag;
int trmax, insmax, insseg[10000][2], insmin=1;		/* The min length of acceptable inserted segments */
int sdmax=-1, sdseg[5000][2];	/* IDRs in a transcript */
int intstart, intend;		/* The start and end res. numbers of internal exons */
double gc, gw, gx, gy, gz, idr[4], tot[4];;
void insread(void){    /* Parse the line to get the start and end of the inserted segment */
        int catcount=0;         /* Category */
	int tstart, tend;	/* The start and end of the inserted segment */
	char ttranscr[51];
	strncpy(ttranscr,empty,51);
        offset = 0;
        while(1){	/* catcount=1: Human ID, 2: Mouse transcr., 3: Length, 4: Class, 5: Seg.start, 6; Seg.end */ 
                ++catcount;
                strncpy(tempstr,empty,51);
                for(x=0; x<50; x++){
                        if(isspace(istr[x+offset])!=0)          break;
                        tempstr[x] = istr[x+offset];
                }
                offset = x + offset;

		if(catcount==2)		strcpy(ttranscr,tempstr);
                else if(catcount==4 && (strstr(tempstr,"IntEx")==NULL && strstr(tempstr,"MultiEx")==NULL)){     /* #&# */
                        return; /* Terminal exon only or Single exon */
                }
                else if(catcount==5)    tstart = atoi(tempstr);       /* Start of the inserted segment */
                else if(catcount==6){
			tend = atoi(tempstr);
			if(tend-tstart+1 < insmin)	return;
			else{
				if(insmax+1>=10000){
					printf("The # of inserted segments exceeded 10,000!\n!");
					exit(1);
				}
				++insmax;
				strcpy(instr[insmax],ttranscr);
				insseg[insmax][0] = tstart;
				insseg[insmax][1] = tend;
				return;
			}	
                }

                for(x=0; x<50; x++)     if(isspace(istr[x+offset])==0)  break;
                offset = x + offset;
        }
        return;
}
struct variant{
        char tname[31];         /* Transcript name */
        int  cdslen;
        char sense[2];          /* + or - */
        int  exnumb;         	/* Exon number */
        int  resstart;          /* The start res # of the internal exons */
        int  resend;            /* The end res # of the internal exons */
};
struct variant trlist[500000];
void tread(void){    /* Parse the line to get the start and end of the internal exons */
        int catcount=0;         /* Category */
        offset = 0;
        strncpy(trlist[trmax].tname,empty,31);          strncpy(trlist[trmax].sense,empty,2);
	trlist[trmax].exnumb = 0;
        trlist[trmax].resstart = 0;                     trlist[trmax].resend = 0;
        while(1){
                ++catcount;
                strncpy(tempstr,empty,51);
                for(x=0; x<30; x++){
                        if(isspace(estr[x+offset])!=0)          break;
                        tempstr[x] = estr[x+offset];
                }
                offset = x + offset;

		/* ### */
                if(catcount==3)    strcpy(trlist[trmax].tname,tempstr);      	/* Transcript name */
                else if(catcount==4)    trlist[trmax].cdslen = atoi(tempstr);   /* CDS length */
                else if(catcount==7)    strcpy(trlist[trmax].sense,tempstr);    /* Sense */
                else if(catcount==8)    trlist[trmax].exnumb = atoi(tempstr);   /* Exon number */
                else if(catcount==9){		/* Total exon number */
			if(atoi(tempstr)==1){
				--trmax;	/* Single exon */
				return;
			}
			else if(trlist[trmax].exnumb==1){
				--trmax;	/* 5'-term exon */
				return;
			}
			else if(atoi(tempstr)==trlist[trmax].exnumb){
				--trmax;	/* 3'-term exon*/
				return;
			}
		}
                else if(catcount==10)    trlist[trmax].resstart = atoi(tempstr);   /* The start res. of the exon */
                else if(catcount==11){
                        trlist[trmax].resend = atoi(tempstr);     /* The end res. of the exon */
                        return;
                }
		/* ### */

                for(x=0; x<30; x++)     if(isspace(estr[x+offset])==0)  break;
                offset = x + offset;
        }
        return;
}
void sdread(void){    /* Parse the line to get the SDs */
        offset = 15;
	while(1){
		++sdmax;
		if(sdmax>=5000){	printf("The # of SDs in %s exceeded 5,000\n",protid);	exit(1);	}
		strncpy(tempstr,empty,51);
		for(x=0; x<50; x++){
			if(isdigit(rstr[x+offset])==0)	break;
			tempstr[x] = rstr[x+offset];
		}
		offset = x + offset + 2;
		sdseg[sdmax][0] = atoi(tempstr);

		strncpy(tempstr,empty,51);
		for(x=0; x<50; x++){
			if(isdigit(rstr[x+offset])==0)	break;
			tempstr[x] = rstr[x+offset];
		}
		offset = x + offset + 1;
		sdseg[sdmax][1] = atoi(tempstr);
		if(isprint(rstr[offset])==0)	break;
	}
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fc=fopen(fnamec,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamec);	exit(1);	}
        if((fe=fopen(fnamee,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamee);	exit(1);	}
        if((fi=fopen(fnamei,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamei);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(x=0; x<4; x++){	/* 0: All, 1: Inserted segments, 2: Other residues in the proteins, 3: Proteins without insertion */	
		idr[x] = 0.;	/* Total # of residues in IDR */
		tot[x] = 0.;	/* Total # of residues */
	}
	for(x=0; x<5000; x++){
		for(y=0; y<2; y++)	sdseg[x][y] = 0;
	}

	/* Read in the list of inserted segments */
	insmax = -1;
	for(x=0; x<1000; x++){
		strncpy(instr[x],empty,51);
		for(y=0; y<2; y++)	insseg[x][y] = 0;	/* y=0: start, 1: end of an ins. seg. */
	}
	strncpy(istr,empty,501);
       	while((fgets(istr,500,fi))!=NULL){
		if(strstr(istr,"ENSMUST")!=NULL)	insread();
		strncpy(istr,empty,501);
	}
	fclose(fi);
	printf("%d inserted segments in internal exons exist\n",insmax+1);
	
	/* Read in the list of transcripts with the start and end of the internal exons */
	strncpy(estr,empty,501);
	trmax = -1;
       	while((fgets(estr,500,fe))!=NULL){
		if(estr[0]!='>' && strstr(estr,"ERROR")==NULL && strstr(estr,"Noncoding")==NULL){
			++trmax;
			if(trmax>=500000){	printf("The # of tr. with internal exons exceeded 500,000!\n");	exit(1);	}
			tread();
		}
		strncpy(estr,empty,501);
	}
	fclose(fe);
	printf("There are %d transcripts with internal exons\n",trmax+1);

        strncpy(rstr,empty,5001);
	strncpy(protid,empty,51);	strncpy(trid,empty,51);		plen = -1; 	imax = -1;	hitflag = 0;
        while((fgets(rstr,5000,fr))!=NULL){
		if(rstr[0]=='>'){
			++counter;
			if(counter%2000==0)	printf("Now processing entry %d out of 51,229...\n",counter);
			hitflag = 1;
			offset = 1;
			for(x=0; x<50; x++){
				if(rstr[x+offset]=='_' || isspace(rstr[x+offset])!=0)	break;
				protid[x] = rstr[x+offset];
			}
			/* Look for the corresponding transcript ID */
			fseek(fc,0L,SEEK_SET);
			strncpy(cstr,empty,122);
        		while((fgets(cstr,121,fc))!=NULL){
				if(strstr(cstr,protid)!=NULL){		/* Hit */
					hitflag = 2;
					for(x=0; x<50; x++)	if(isspace(cstr[x])!=0)		break;
					offset = x; 
					for(x=0; x<50; x++)	if(isspace(cstr[x+offset])==0)	break;
					offset = x + offset;
					for(x=0; x<50; x++){
						if(isprint(cstr[x+offset])==0)	break;
						trid[x] = cstr[x+offset];
					}
				}	
				strncpy(cstr,empty,122);
			}
		}
		else if(strncmp(rstr,"LENGTH",6)==0 && hitflag==2){
			hitflag = 3;
			strncpy(tempstr,empty,51);
			if(strncmp(rstr,"LENGTH=",7)==0)	offset = 7;
			else	offset = 15;
			for(x=0; x<50; x++){
				if(isprint(rstr[x+offset])==0)	break;
				tempstr[x] = rstr[x+offset];
			}
			plen = atoi(tempstr);

			/* Look for the start and end of internal exons, if any */
		 	intstart = 0;
			intend = 0;
			for(x=0; x<=trmax; x++){
				if(strcmp(trid,trlist[x].tname)==0){	/* Transcript names agreed */
					if(trlist[x].cdslen==plen || trlist[x].cdslen==plen-1){
						/* Protein lengths agreed witnin the allowance of 1 residue */
						intstart = trlist[x].resstart;
						intend = trlist[x].resend;
						hitflag = 4;
						goto nextr;
					}
				}
			}
			/* Not found: protein lengths diaagreement and S, D cases */
			hitflag = 0;
		}
		else if(hitflag==4 && strncmp(rstr,"STRDOM",6)==0){	/* SD line */
			sdread();
		}
		else if(strncmp(rstr,"//",2)==0){
			if(hitflag!=4)		goto reinit;
			else if(intstart<=0 || intend<=0)	goto reinit;

			/* Look for inserted segments, if any */
			insflag = 0;
			for(x=0; x<=insmax; x++){
				if(strcmp(instr[x],trid)==0){
					insflag = 1;		/* There is at least one ins. segment in the transcript */
					break;
				}
			}

			/* Examine one residue in the internal exon segnment  at a time */
			for(x=intstart; x<=intend; x++){
				/* Check if x is in an SD */
				sflag = 0;		/* sflag becomes 1 for SD */
				if(sdmax>=0){	/* ### Corrected */
					for(y=0; y<=sdmax; y++){	/* ### Corrected */	
						if(sdseg[y][0]<=x && x<=sdseg[y][1]){	/* x is in SD */	
							sflag = 1;
							break;
						}
					}
				}

				if(insflag==1){
					for(y=0; y<=insmax; y++){
						if(strcmp(instr[y],trid)==0 && insseg[y][0]<=x && x<=insseg[y][1]){ /* x is in an ins.seg */
							tot[1] += 1.;
							if(sflag==0)	idr[1] += 1.;		/* ### */
							goto nextx;
						}
					}
					tot[2] += 1.;
					if(sflag==0)	idr[2] += 1.;		/* ### */
				}
				else{	/* The transript does not have an inserted segment */
					tot[3] += 1.;
					if(sflag==0)	idr[3] += 1.;		/* ### */
				}
nextx:				x = x;
			}
	
reinit:			/* Reninitializations */
			hitflag = 0;
			strncpy(protid,empty,51);	strncpy(trid,empty,51);		plen = -1; 	imax = -1;
			sdmax = -1;
			for(x=0; x<5000; x++){
				for(y=0; y<2; y++)	sdseg[x][y] = 0;
			}
			intstart = 0; 		intend = 0;
		}
       
nextr:	        strncpy(rstr,empty,5001);
        }

	for(x=1; x<4; x++){
		tot[0] += tot[x];
		idr[0] += idr[x];
	}

	/* Output fractions of IDR */
	fputs("DISOPRED3\n",fw);
	fputs("Class\t1:Ins.seg\t2:Non-ins.seg.\t1+2:Prot.w/ins.seg.\t3:Prot.w/o ins.seg\t0:All internal exons\n",fw);
	fputs("%IDR\t",fw);
	for(x=1; x<=2; x++){
		if(tot[x]>0.){
			gw = 100.*idr[x]/tot[x];
			strncpy(tempstr,empty,51);
			sprintf(tempstr,"%.2f",gw);
			fputs(tempstr,fw);	fputs("\t",fw);
		}
		else	fputs("-\t",fw);
	}
	/* Classes 1+2 */
	if(tot[1]+tot[2]>0.){
		gw = 100.*(idr[1]+idr[2])/(tot[1]+tot[2]);
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);
	}
	else	fputs("-\t",fw);

	for(x=3; x<=3; x++){
		if(tot[x]>0.){
			gw = 100.*idr[x]/tot[x];
			strncpy(tempstr,empty,51);
			sprintf(tempstr,"%.2f",gw);
			fputs(tempstr,fw);	fputs("\t",fw);
		}
		else	fputs("-\t",fw);
	}
	if(tot[0]>0.){
		gw = 100.*idr[0]/tot[0];
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gw);
		fputs(tempstr,fw);	fputs("\n",fw);
	}
	else	fputs("-\n",fw);
		
	/* ## Output the #s: Part 1 ## */
	fputs("\nClass\t1:Ins.Seg\t\t\t2:Non-ins. seg.\t\t\t1+2:Prot.w/ins seg.\t\n",fw);
	fputs("\tObs\tExp\tChiSq\tObs\tChiSq\tExp\tObs\tChiSq\n",fw);

	fputs("IDR\t",fw);
	gc = 0.;	/* Total sum of Chi sq. */
	gx = 0.;	/* Sum of Chi sq. in IDR */
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.0f",idr[x]);			/* Obs */
		fputs(tempstr,fw);	fputs("\t",fw);

		strncpy(tempstr,empty,51);
		gw = tot[x]*(idr[1]+idr[2])/(tot[1]+tot[2]);	/* Expected */
		sprintf(tempstr,"%.0f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);

		gy = idr[x] - gw;	/* Obs - Exp */
		gz = gy*gy/gw;		/* Chi-sq */
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);	fputs("\t",fw);
		gx += gz;
		gc += gz;
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",idr[1]+idr[2]);		/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gx);			/* Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	fputs("SD\t",fw);
	gx = 0.;	/* Sum of Chi sq. in SD */
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.0f",tot[x]-idr[x]);		/* Obs */
		fputs(tempstr,fw);	fputs("\t",fw);

		strncpy(tempstr,empty,51);
		gw = tot[x]*(tot[1]-idr[1]+tot[2]-idr[2])/(tot[1]+tot[2]);	/* Expected */
		sprintf(tempstr,"%.0f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);

		gy = tot[x]-idr[x] - gw;	/* Obs - Exp */
		gz = gy*gy/gw;		/* Chi-sq */
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);	fputs("\t",fw);
		gx += gz;
		gc += gz;
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",tot[1]-idr[1]+tot[2]-idr[2]);	/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gx);			/* Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	fputs("IDR+SD\t",fw);
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.0f",tot[x]);			/* Obs */
		fputs(tempstr,fw);	fputs("\t\t\t",fw);
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",tot[1]+tot[2]);			/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gc);			/* Total Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	/* ## Output the #s: Part 2 ## */
	fputs("\nClass\t1:Ins.Seg\t\t\t2+3:Rest of int.exons\t\t\t1+2+3:All int.exons\t\n",fw);
	fputs("\tObs\tExp\tChiSq\tObs\tChiSq\tExp\tObs\tChiSq\n",fw);
	fputs("IDR\t",fw);
	gc = 0.;	/* Total sum of Chi sq. */
	gx = 0.;	/* Sum of Chi sq. in IDR */
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		if(x==1)	sprintf(tempstr,"%.0f",idr[x]);			/* Obs */
		else		sprintf(tempstr,"%.0f",idr[2]+idr[3]);
		fputs(tempstr,fw);	fputs("\t",fw);

		strncpy(tempstr,empty,51);
		if(x==1)	gw = tot[x]*idr[0]/tot[0];	/* Expected */
		else		gw = (tot[2]+tot[3])*idr[0]/tot[0];
		sprintf(tempstr,"%.0f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);

		if(x==1)	gy = idr[x] - gw;	/* Obs - Exp */
		else		gy = idr[2]+idr[3] - gw;
		gz = gy*gy/gw;		/* Chi-sq */
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);	fputs("\t",fw);
		gx += gz;
		gc += gz;
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",idr[0]);			/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gx);			/* Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	fputs("SD\t",fw);
	gx = 0.;	/* Sum of Chi sq. in SD */
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		if(x==1)	sprintf(tempstr,"%.0f",tot[x]-idr[x]);		/* Obs */
		else		sprintf(tempstr,"%.0f",tot[2]-idr[2]+tot[3]-idr[3]);
		fputs(tempstr,fw);	fputs("\t",fw);

		strncpy(tempstr,empty,51);
		if(x==1)	gw = tot[x]*(tot[0]-idr[0])/tot[0];	/* Expected */
		else		gw = (tot[2]+tot[3])*(tot[0]-idr[0])/tot[0];
		sprintf(tempstr,"%.0f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);

		if(x==1)	gy = tot[x]-idr[x] - gw;	/* Obs - Exp */
		else		gy = tot[2]-idr[2]+tot[3]-idr[3] -gw;
		gz = gy*gy/gw;		/* Chi-sq */
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);	fputs("\t",fw);
		gx += gz;
		gc += gz;
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",tot[0]-idr[0]);	/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gx);			/* Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	fputs("IDR+SD\t",fw);
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		if(x==1)	sprintf(tempstr,"%.0f",tot[x]);			/* Obs */
		else		sprintf(tempstr,"%.0f",tot[2]+tot[3]);
		fputs(tempstr,fw);	fputs("\t\t\t",fw);
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",tot[0]);			/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gc);			/* Total Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);
	
	/* ## Output the #s: Part 3 ## */
	fputs("\nClass\t1+2:Prot.w/ins.seg\t\t\t3:Prot.w/o ins.seg\t\t\t1+2+3:All int.exons\t\n",fw);
	fputs("\tObs\tExp\tChiSq\tObs\tChiSq\tExp\tObs\tChiSq\n",fw);
	fputs("IDR\t",fw);
	gc = 0.;	/* Total sum of Chi sq. */
	gx = 0.;	/* Sum of Chi sq. in IDR */
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		if(x==1)	sprintf(tempstr,"%.0f",idr[1]+idr[2]);			/* Obs */
		else		sprintf(tempstr,"%.0f",idr[3]);
		fputs(tempstr,fw);	fputs("\t",fw);

		strncpy(tempstr,empty,51);
		if(x==1)	gw = (tot[1]+tot[2])*idr[0]/tot[0];	/* Expected */
		else		gw = tot[3]*idr[0]/tot[0];
		sprintf(tempstr,"%.0f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);

		if(x==1)	gy = idr[1]+idr[2] - gw;	/* Obs - Exp */
		else		gy = idr[3] - gw;
		gz = gy*gy/gw;		/* Chi-sq */
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);	fputs("\t",fw);
		gx += gz;
		gc += gz;
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",idr[0]);			/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gx);			/* Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	fputs("SD\t",fw);
	gx = 0.;	/* Sum of Chi sq. in SD */
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		if(x==1)	sprintf(tempstr,"%.0f",tot[1]-idr[1]+tot[2]-idr[2]);		/* Obs */
		else		sprintf(tempstr,"%.0f",tot[3]-idr[3]);
		fputs(tempstr,fw);	fputs("\t",fw);

		strncpy(tempstr,empty,51);
		if(x==1)	gw = (tot[1]+tot[2])*(tot[0]-idr[0])/tot[0];	/* Expected */
		else		gw = tot[3]*(tot[0]-idr[0])/tot[0];
		sprintf(tempstr,"%.0f",gw);
		fputs(tempstr,fw);	fputs("\t",fw);

		if(x==1)	gy = tot[1]-idr[1] +tot[2]-idr[2] - gw;	/* Obs - Exp */
		else		gy = tot[3]-idr[3] -gw;
		gz = gy*gy/gw;		/* Chi-sq */
		strncpy(tempstr,empty,51);
		sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);	fputs("\t",fw);
		gx += gz;
		gc += gz;
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",tot[0]-idr[0]);	/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gx);			/* Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);

	fputs("IDR+SD\t",fw);
	for(x=1; x<=2; x++){
		strncpy(tempstr,empty,51);
		if(x==1)	sprintf(tempstr,"%.0f",tot[1]+tot[2]);			/* Obs */
		else		sprintf(tempstr,"%.0f",tot[3]);
		fputs(tempstr,fw);	fputs("\t\t\t",fw);
	}
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.0f",tot[0]);			/* Obs */
	fputs(tempstr,fw);	fputs("\t",fw);
	strncpy(tempstr,empty,51);
	sprintf(tempstr,"%.2f",gc);			/* Total Chi-sq */
	fputs(tempstr,fw);	fputs("\n",fw);
	
        fclose(fr);    	fclose(fw);	fclose(fc);
	return 0;
}
