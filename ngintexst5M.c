#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fi, *fw, *fl, *fc;		/* ##### */
char fnamei[122]="/Users/khomma/work/dichot/ndichot1";                  /* DICHOT results */
char fnamer[122]="/Users/khomma/work/ensembl/rat/nalluniqexon3A";	/* #### Corrected list of ALL internal exons */
char fnamew[122]="/Users/khomma/work/ensembl/rat/ngintexst5M";		/* ### Avg and median values */
char fnamel[122]="/Users/khomma/work/ensembl/rat/ngintexstlist5M";	/* ### List of data */
char fnamec[122]="/Users/khomma/work/ensembl/rat/tptcorr1";             /* ### ENSRNOP-ENSRNOT correspondence */
char binstr[10][20]={"All","-120","-240","-360","-480","-600","-720","-840","-960","961-nt"};	/* #### */
char rstr[5001], istr[5001], empty[5010]="", tempstr[21], entryname[51], outstr[500], listread[21];	/* #### No abbentname */
char tchr[10], tgname[51], oldgname[51], ttname[51], oldtname[51], tsense[2];
char idcorr[70000][2][21];       /* ##### ENSRNOP-ENSRNOT correspondence */
int  p, q, t, u, v, w, x, y, z, offset, excount=0, vcounter=0, vhit=0, hitcount=0, cmax=-1;	/* ##### */
int  dmax=1000, sd[1001][2], sdmax=-1, ihit;
int  nmin=20;	/* nmin is the min # for fr. SD calculation */	
int  exlenbinmax=10, exnumb, smax, vmax=-1;	/* #### */
int  slist[500000], smax;	/* #### slist is the array for sorting */
int  tcdslen, oldcdslen, tstart, tend, tresstart, tresend;
int  numbsd, length;		/* The # of residues in SD in each exon and the length of each exon */
double gs, gt, gu, gv, gw, gx, gy, gz, gf, scaling=1000.;	/* #&# */
double pnumb[10], ptsd1[10], ptsd2[10], frsd1[10], frsd2[10], ptavsd[10], ptsemsd[10];	/* #### For statistics of SD contents */
	/* 0:All exon len., 1: 1-120 nt, 2: -240,..., 9: 961 nt- #### */
double totpdis1[10], totpdis2[10];  /* #### Total lengths of exons in each length bin, to calculate the avg. exon length */
double sn=0., sx1=0., sx2=0., sy1=0., sy2=0., sxy=0., corrr, corrt, corra, corrb;       /* Corr. of internal exon len. with %IDR #&# */
struct variant{
        char chr[51];		/* Chromosome name */
        char gname[51];         /* Gene name */
        char tname[51];         /* Transcript name */
        int  cdslen;
        int  start;		/* Gene start address */
        int  end;		/* Gene end address */
        char sense[2];          /* + or - */
        int  resstart;          /* The start res # of the exon */
        int  resend;            /* The end res # of the exon */
};
struct variant vlist[10000];
void tread(char rstr[122]){    /* Parse the line to get the start and end of the terminal exon */
        int catcount=0;         /* Category vhit */
        offset = 0;
	strncpy(tchr,empty,10);		strncpy(tgname,empty,51);	strncpy(ttname,empty,51);	strncpy(tsense,empty,2);
	tstart = 0;	tend = 0;	tresstart = 0;		tresend = 0;
        while(1){
                ++catcount;
                strncpy(tempstr,empty,21);
                for(x=0; x<50; x++){
                        if(isspace(rstr[x+offset])!=0)          break;
                        tempstr[x] = rstr[x+offset];
                }
                offset = x + offset;
                if(catcount==1)         strcpy(tchr,tempstr);        /* Chromosome name */
                else if(catcount==2)    strcpy(tgname,tempstr);      /* Gene name */
                else if(catcount==3)    strcpy(ttname,tempstr);      /* Transcript name */
                else if(catcount==4)    tcdslen = atoi(tempstr);       /* CDS length */
                else if(catcount==5)    tstart = atoi(tempstr);       /* Start */
                else if(catcount==6)    tend = atoi(tempstr);        /* End */
                else if(catcount==7)    strcpy(tsense,tempstr);      /* Sense */
                else if(catcount==8)    tresstart = atoi(tempstr);   /* The start res. of the exon */
                else if(catcount==9){
                        tresend = atoi(tempstr);     /* The end res. of the exon */
                        return;
                }

                for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;
                offset = x + offset;
        }
        return;
}
void assign(void)	/* Find the SD assignments, if any */
{
	int aw, ax, ay, az;
	fseek(fi,0L,SEEK_SET);
	sdmax = -1;
	for(ax=0; ax<1000; ax++){
		for(ay=0; ay<2; ay++)	sd[ax][ay] = 0;		/* start and end of each SD */
	}
	strncpy(istr,empty,5001);
       	while((fgets(istr,5000,fi))!=NULL){
		if(istr[0]=='>' && strstr(istr,entryname)!=NULL){	/* Possible hit */
			ihit = 1;
		}
		else if(ihit==1 && strncmp(istr,"LENGTH",5)==0){	/* length line */
			offset = 15;
			strncpy(tempstr,empty,21);
			for(ax=0; ax<20; ax++){
				if(isprint(istr[ax+offset])==0)	break;
				tempstr[ax] = istr[ax+offset];
			}
			aw = atoi(tempstr);
			if(aw!=oldcdslen)	ihit = 0;	/* False hit (a different variant) */
			else			ihit = 2;	/* Real hit */
		}
		else if(ihit==2 && (strncmp(istr,"KNWDOM",6)==0 || strncmp(istr,"NEWDOM",6)==0)){       /* SD line */
			offset = 15;
			if(sdmax>=1000){
				printf("The # of SDs in %s exceeded 1,000\n",entryname);
				exit(1);
			}
                       	while(1){
				++sdmax;
                       		strncpy(tempstr,empty,21);
                       		for(ax=0; ax<20; ax++){
                       			if(istr[ax+offset]=='-')         break;
                               			tempstr[ax] = istr[ax+offset];
                       		}
				sd[sdmax][0] = atoi(tempstr);	/* Start */
                       		offset = ax + offset + 2;
                       		strncpy(tempstr,empty,21);
                       		for(ax=0; ax<20; ax++){
                       			if(isprint(istr[ax+offset])==0)  break;
                       			else if(istr[ax+offset]=='|')    break;
                       			tempstr[ax] = istr[ax+offset];
                    		}
				sd[sdmax][1] = atoi(tempstr);	/* End */
                       		if(isprint(istr[ax+offset])==0)  break;
				else     offset = ax + offset + 1;
                       	}
		}
		else if(ihit==2 && strncmp(istr,"//",2)==0){
			ihit = 3;
			return;
		}
		strncpy(istr,empty,5001);
	}
	return;
}
void sdfraction(void)	/* Determine and output the fraction of SD in each exon */
{
	int kx, ky, kz;
	for(ky=0; ky<=vmax; ky++){		/* One (internal) exon at a time */
		++hitcount;	/* The total number of variants evaluated */
		length = abs(vlist[ky].end - vlist[ky].start) + 1;	/* #### Length of exon y in nt */
		q = (length-1)/120 + 1;		/* #### Exon length bin */
		if(q>exlenbinmax-1)	q = exlenbinmax -1;

		numbsd = 0.;	/* The # of residues in SD in each exon */
		if(sdmax>=0){
			for(kz=vlist[ky].resstart; kz<=vlist[ky].resend; kz++){	
				/* Examine all the residues (kz) in the exon */
				for(kx=0; kx<=sdmax; kx++){
					if(sd[kx][0]<=kz && kz<=sd[kx][1]){	/* Hit! The residue is in SD */
						++numbsd;
						break;
					}
				}
			}
		}
		pnumb[q] += 1.;
		gx = (double)(vlist[ky].resend - vlist[ky].resstart + 1);	/* #### CDS (in aa) */
		gy =  (double)numbsd;		/* #### */

		totpdis1[q] += gx/scaling;	/* Sum of exon len. in each bin */
		totpdis2[q] += (gx/scaling)*(gx/scaling);	/* Sum of exon len. in each bin */
		ptsd1[q] +=  gy/scaling;		/* Sum of residues in SD */
		ptsd2[q] += (gy/scaling)*(gy/scaling);

		frsd1[q] += gy/gx/scaling;		/* Fraction of SD */
		frsd2[q] += (gy/gx/scaling)*(gy/gx/scaling);

		/* #&# */
                gz = 1.-(gy/gx);        /* Fraction of IDR */
                sn += 1.;
                sx1 += gx/scaling;
                sx2 += (gx/scaling)*(gx/scaling);
                sy1 += gz;
                sy2 += gz*gz;
                sxy += gz*gx/scaling;
		/* #&# */

		strncpy(outstr,empty,500);
		strcpy(outstr,vlist[ky].gname);	strcat(outstr,"\t");
		strcat(outstr,vlist[ky].tname);	strcat(outstr,"\t");

		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",vlist[ky].cdslen);	/* Length of the entire protein */
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",vlist[ky].start);	/* Exon start address */	
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",vlist[ky].end);	/* Exon end address */
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strcat(outstr,vlist[ky].sense); strcat(outstr,"\t");		/* Sense */	
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",vlist[ky].resstart);	 /* Starting coding res. */
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",vlist[ky].resend); 	/* Ending coding res. */
		strcat(outstr,tempstr);		strcat(outstr,"\t");

		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",q);        /* Exon length bin */	
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",numbsd);   /* Residue # in SD */
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",length);   /* Exon length in nt */	
		strcat(outstr,tempstr);		strcat(outstr,"\t");
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",100.*gy/gx);     /* %SD */
		strcat(outstr,tempstr);		strcat(outstr,"\n");
		fputs(outstr,fl);
		ky = ky;
	}
	return;
}
void kokan(int a, int b)
{
        if(a==b)        return;
        t = slist[a];
        slist[a] = slist[b];
        slist[b] = t;
	return;
}
void quick_sort(int left, int right)
{
        int last, i;
        if(left >= right)       return;
        kokan(left, (left + right)/2);
        last = left;
        for(i=left + 1; i<= right; i++){
                if(slist[left]>slist[i]){
                        ++last;
                        kokan(last, i);
                }
        }
        kokan(left, last);
        quick_sort(left, last-1);
        quick_sort(last+1, right);
	return;
}
void listparse(int i)
{
	int j=0, k;
	offset = 0;
	strncpy(listread,empty,21);
	while(j<i){
		j++;
		strncpy(listread,empty,21);
		for(k=0; k<20; k++){
			if(isspace(rstr[k+offset])!=0)	break;
			listread[k] = rstr[k+offset];
		}
		offset = k + offset;
		for(k=0; k<50; k++)	if(isspace(rstr[k+offset])==0)	break;
		offset = k + offset;
	}
	return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
	if((fi=fopen(fnamei,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamei);      exit(1);        }
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
        if((fl=fopen(fnamel,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamel);     exit(1);        }

        /* ##### */
        if((fc=fopen(fnamec,"r")) == NULL){     printf("Cannot open the output file %s\n", fnamec);     exit(1);        }
        for(x=0; x<70000; x++){
                for(y=0; y<2; y++)      /* 0: ENSRNOP, 1: ENSRNOT */    strncpy(idcorr[x][y],empty,21);
        }
        strncpy(rstr,empty,5001);
        while((fgets(rstr,5000,fc))!=NULL){
                ++cmax;
                if(cmax>=70000){        printf("The # of ID correspondences exceeded 70,000\n");        exit(1);        }
                offset = 0;
                for(x=0; x<20; x++){
                        if(isspace(rstr[x+offset])!=0)  break;
                        idcorr[cmax][0][x] = rstr[x+offset];
                }
                offset = x + offset;
                for(x=0; x<50; x++)     if(isspace(rstr[x+offset])==0)  break;
                offset = x + offset;
                for(x=0; x<20; x++){
                        if(isspace(rstr[x+offset])!=0)  break;
                        idcorr[cmax][1][x] = rstr[x+offset];
                }
                strncpy(rstr,empty,5001);
        }
        fclose(fc);
        printf("%d correspondences of IDs were found\n",cmax);
        /* ##### */

	fputs(">Gene name\tVar name\tLen\tIntExStart\tIntExEnd\tSense\tRes start\tRes end\tExon len bin\tSD res\tExon len(nt)\t%SD\n",fl);
	vhit = 0; 		hitcount = 0;
	for(x=0; x<exlenbinmax; x++){		/* Exon length bins */
		pnumb[x] = 0.;		/* Total # of samples */
		ptsd1[x] = 0.;		/* Sum of total residues in SDs */
		ptsd2[x] = 0.;		/* Sum of (tot. res. in SDs)^2 */
		frsd1[x] = 0.;		/* Sum of fractional SD */
		frsd2[x] = 0.;		/* Sum of fr. SD^2 */
		ptavsd[x] = 0.;		/* Avg of total residues in SDs */
		ptsemsd[x] = 0.;	/* SEM of total residues in SDs */
		totpdis1[x] = 0.;	/* Total exon lengths */
		totpdis2[x] = 0.;	/* Total exon lengths^2 */
	}
	strncpy(entryname,empty,51);
        vmax = -1;
	strncpy(rstr,empty,5001);
        while((fgets(rstr,5000,fr))!=NULL){
		if(rstr[0]=='>')	goto nextu;
		else if(strncmp(rstr,"MT",2)==0){		/* Exclude non-nuclear coded genes */
			printf("Excluding %s",rstr);
			goto nextu;	
		}
		else{
			++excount;
			if(excount%10000==0)	printf("Now processing exon %d out of 210,853...\n",excount);	 /* #### */
			strncpy(oldgname,empty,51);	strncpy(oldtname,empty,51);
			strcpy(oldgname,tgname);	strcpy(oldtname,ttname);
			oldcdslen = tcdslen;
			tread(rstr);
			if(strcmp(ttname,oldtname)==0){	
				/* The exon is in the same variant as the previous one. Add data to the variant structure: */
				++vmax;
				if(vmax>=10000){ 	printf("The # of variants in %s exceeded 10,000\n",tgname);	exit(1);	}
				strcpy(vlist[vmax].chr,tchr);
				strcpy(vlist[vmax].gname,tgname);
				strcpy(vlist[vmax].tname,ttname);
				vlist[vmax].cdslen = tcdslen;
				vlist[vmax].start = tstart;
				vlist[vmax].end = tend;
				strcpy(vlist[vmax].sense,tsense);
				vlist[vmax].resstart = tresstart;
				vlist[vmax].resend =  tresend;
				goto nextu;
			}
			else{
				/* The current exon belongs to a new variant. Process data of the OLD variant */
				++vcounter;
				if(strlen(oldgname)<1)		goto update;	/* The first line, for which oldgname is empty */

                                /* ##### Find the ENSRNOP ID that corresponds to the ENSRNOT ID (oldtname) */
                                strncpy(entryname,empty,51);
                                for(x=0; x<=cmax; x++){
                                        if(strcmp(oldtname,idcorr[x][1])==0){   /* Found the correspondence */
                                                strcpy(entryname,idcorr[x][0]);         /* ENSRNOP ID */
                                                goto okcorr;
                                        }
                                }
                                goto update;    /* No correspondence found */
                                /* ##### */

okcorr:				/* ##### Find the IDR and SD assignments: */
				ihit = 0;	/* ihit: 1 -> possible hit, ihit: 2 -> real hit (i.e., the lengths agree) */
				assign();
				if(ihit!=3)	goto update;	/* ihit becomes 3 if assignments are found */

				/* Calcualte the avg. # of residues in SD and length of each exon of the protein: */
				++vhit;
				sdfraction();

update:				/* Update the data of the new variant: */
				vmax = -1;
				++vmax;		/* Actually 0 */
				strcpy(vlist[vmax].chr,tchr);
				strcpy(vlist[vmax].gname,tgname);
				strcpy(vlist[vmax].tname,ttname);
				vlist[vmax].cdslen = tcdslen;
				vlist[vmax].start = tstart;
				vlist[vmax].end = tend;
				strcpy(vlist[vmax].sense,tsense);
				vlist[vmax].resstart = tresstart;
				vlist[vmax].resend =  tresend;
			}
		}
                
nextu:		strncpy(rstr,empty,5001);
       	}
	/* Process the last variant */
	++vcounter;
        /* ##### Find the ENSRNOP ID that corresponds to the ENSRNOT ID (oldtname) */
        strncpy(entryname,empty,51);
        for(x=0; x<=cmax; x++){
                if(strcmp(oldtname,idcorr[x][1])==0){   /* Found the correspondence */
                        strcpy(entryname,idcorr[x][0]);         /* ENSRNOP ID */
                        goto fndcor;
                }
        }
        goto ctot;      /* No correspondence found */
        /* ##### */

fndcor: /* ###  Find the IDR and SD assignments: */
	ihit = 0;	/* ihit: 1 -> possible hit, ihit: 2 -> real hit (i.e., the lengths agree) */
	assign();
	if(ihit!=3)	goto ctot;	/* ihit becomes 3 if assignments are found */
	/* Calcualte the avg. # of residues in SD and length of each exon of the protein: */
	++vhit;
	sdfraction();
	
ctot:	/* Calculate the total of all exon length bins */
	for(q=1; q<exlenbinmax; q++){
		pnumb[0] += pnumb[q];
		totpdis1[0] += totpdis1[q];
		totpdis2[0] += totpdis2[q];
		ptsd1[0] += ptsd1[q];
		ptsd2[0] += ptsd2[q];
		frsd1[0] += frsd1[q];
		frsd2[0] += frsd2[q];
	}
	/* Calculate and output the avg and SEM of the number of SD residues in each exon length bin */
	fputs("DICHOT\tLen. range\tNumber\tAvg len\tSEM\tAvg SD res.\t SEM\tAvg %SD\tSEM\n",fw);
	for(q=0; q<exlenbinmax; q++){
		if(q==0)	fputs("All bins\t",fw);
		else{
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",q);
			fputs("Bin ",fw);		fputs(tempstr,fw);	fputs("\t",fw);
		}
		fputs(binstr[q],fw);		fputs("\t",fw);

		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.0f",pnumb[q]);
		fputs(tempstr,fw);		fputs("\t",fw);
		if(pnumb[q]<nmin){
			fputs("\n",fw);
			continue;
		}
		gx = scaling*totpdis1[q]/pnumb[q];	/* Avg. exon length in each exon length bin */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.1f",gx);
		fputs(tempstr,fw);		fputs("\t",fw);

		gy = totpdis2[q] - (totpdis1[q]*totpdis1[q]/pnumb[q]);
		if(gy<0.){	printf("gy=%f at q=%d\n",gy,q);		exit(1);	}
		gz = scaling*sqrt(gy)/pnumb[q];
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);		fputs("\t",fw);

		gx = scaling*ptsd1[q]/pnumb[q];	/* Avg. # of residues in SD in each exon length bin */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.1f",gx);
		fputs(tempstr,fw);		fputs("\t",fw);

		gy = ptsd2[q] - (ptsd1[q]*ptsd1[q]/pnumb[q]);
		if(gy<0.){	printf("gy=%f at q=%d\n",gy,q);		exit(1);	}
		gz = scaling*sqrt(gy)/pnumb[q];
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gz);
		fputs(tempstr,fw);		fputs("\t",fw);

		gx = 100.*scaling*frsd1[q]/pnumb[q];	/* Avg. %SD in each exon length bin */
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.2f",gx);
		fputs(tempstr,fw);		fputs("\t",fw);

		gy = frsd2[q] - (frsd1[q]*frsd1[q]/pnumb[q]);
		if(gy<0.){	printf("gy=%f at q=%d\n",gy,q);		exit(1);	}
		gz = 100.*scaling*sqrt(gy)/pnumb[q];
		strncpy(tempstr,empty,21);	sprintf(tempstr,"%.3f",gz);
		fputs(tempstr,fw);		fputs("\n",fw);
	}
	fclose(fr);	fclose(fl);
	printf("%d exons of %d variants out of a total of %d exons in %d variants were analyzed\n",hitcount,vhit,excount,vcounter);

	printf("Now calculating medians...\n");
	fputs("\nExon len bin\tLen range\tNumber\tMedian SD res #\tMedian tot len\tMedian %SD\n",fw);
	if((fl=fopen(fnamel,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamel);      exit(1);        }
	for(w=1; w<exlenbinmax; w++){
		if(w>1)		fseek(fl,0L,SEEK_SET);
		strncpy(tempstr,empty,21);
		sprintf(tempstr,"%d",w);
		fputs("Bin ",fw);
		fputs(tempstr,fw);	fputs("\t",fw);
		fputs(binstr[w],fw);	fputs("\t",fw);

		/* Read in the saved temporary files and detemine the median value of SD residues */
		smax = -1;
		for(x=0; x<500000; x++)	slist[x] = -1;		/* #### */
		strncpy(rstr,empty,5001);
        	while((fgets(rstr,5000,fl))!=NULL){
			if(rstr[0]=='>')	goto nextl;
			listparse(9);
			if(w!=atoi(listread))	goto nextl;	/* Wrong bin length bin */	
			++smax;
			if(smax>=500000){	printf("The # of elements in %s exceeded 500,000\n",fnamer);	exit(1);	}
				/* #### */
			listparse(10);		/* # of residues in SD */
			slist[smax] = atoi(listread);
nextl:			strncpy(rstr,empty,5001);
		}
		strncpy(tempstr,empty,21);
		sprintf(tempstr,"%d",smax+1);
		fputs(tempstr,fw);		fputs("\t",fw);
		if(smax<0){	fputs("\n",fw);		goto nextw;	}
		quick_sort(0, smax);
		if(smax%2==0){
			x = smax/2;
			gx = slist[x];
		}
		else{
			x = (smax-1)/2;
			gx = (slist[x]+slist[x+1])/2.;
		}
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.1f",gx);
                fputs(tempstr,fw);      fputs("\t",fw);

		/* Read in the saved list and detemine the median value of total length */
		fseek(fl,0L,SEEK_SET);
		smax = -1;
		for(x=0; x<500000; x++)	slist[x] = -1;		/* #### */
		strncpy(rstr,empty,5001);
        	while((fgets(rstr,5000,fl))!=NULL){
			if(rstr[0]=='>')	goto nextd;
			listparse(9);
			if(w!=atoi(listread))	goto nextd;	/* Wrong bin length bin */	
			++smax;
			if(smax>=500000){	printf("The # of elements in %s exceeded 500,000\n",fnamer);	exit(1);	}
				/* #### */
			listparse(11);		/* Exon length (a.a.) */
			slist[smax] = atoi(listread);
nextd:			strncpy(rstr,empty,5001);
		}
		quick_sort(0, smax);
		if(smax%2==0){
			x = smax/2;
			gx = slist[x];
		}
		else{
			x = (smax-1)/2;
			gx = (slist[x]+slist[x+1])/2.;
		}
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.1f",gx);
                fputs(tempstr,fw);      fputs("\t",fw);

		/* Read in the saved temporary files and detemine the median value of per cent SD */
		fseek(fl,0L,SEEK_SET);
		smax = -1;
		for(x=0; x<500000; x++)	slist[x] = -1;		/* #### */		

		strncpy(rstr,empty,5001);
        	while((fgets(rstr,5000,fl))!=NULL){
			if(rstr[0]=='>')	goto nextr;
			listparse(9);
			if(w!=atoi(listread))	goto nextr;	/* Wrong bin length bin */	

			++smax;
			if(smax>=500000){	printf("The # of elements in %s exceeded 500,000\n",fnamer);	exit(1);	}
				/* #### */
			listparse(12);		/* Per cent SD */
			gx = scaling*atof(listread);
			u = gx;
			slist[smax] = u;
nextr:			strncpy(rstr,empty,5001);
		}
		quick_sort(0, smax);
		if(smax%2==0){
			x = smax/2;
			gx = slist[x];
		}
		else{
			x = (smax-1)/2;
			gx = (slist[x]+slist[x+1])/2.;
		}
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%.2f",gx/scaling);
                fputs(tempstr,fw);      fputs("\n",fw);
nextw:		w = w;
	}
        /* #&# Calculate and output the correlation of internal exon lengths (x) with %IDR (y) */
        fputs("\nCorr. of each internal exon length with the %IDR it encodes\n",fw);
        fputs("#\tr\tt\tSlope a\tIntercept b\n",fw);
        strncpy(tempstr,empty,21);      sprintf(tempstr,"%.0f",sn);   fputs(tempstr,fw);
        if(sn>=3.){
                fputs("\t",fw);
                gs = sxy -(sx1*sy1/sn);
                gw = sx2 -(sx1*sx1/sn);
                gu = sy2 -(sy1*sy1/sn);
                if(gw*gu<0.){   printf("gt=%f, gu=%f at #100\n",gw,gu);         exit(1);        }
                corrr = gs/sqrt(gw*gu);    /* Correlation coeff. */
                corrt = corrr *sqrt(sn-2.)/sqrt(1.- corrr*corrr);
                corra = gs/gw;          /* Slope of the least sq. fit */
                corrb = sy1/sn -corra*sx1/sn;       /* y intercept of the least sq. fit */
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.4f",corrr);  fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",corrt);  fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.4f",100.*corra/scaling);       fputs(tempstr,fw);      fputs("\t",fw);
                        /* if y is expressed in per cent, a should be multiplied by 100 */
                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.4f",100.*corrb);    fputs(tempstr,fw);      fputs("\n",fw);
        }
        else    fputs("\n",fw);

	fclose(fl);	fclose(fw);
	return 0;
}
