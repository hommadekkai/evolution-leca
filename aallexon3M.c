#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw;
char empty[1002]="", gstr[1001], trname[51],  tempstr[21], class[21];
char *ps, *pt;
char chname[50][21];	/* I, II, ..., X, Y, MtDNA */
int exon[40000][6], exmax=-1, cds[40000][4], cdsmax=-1;		/* &&& */
        /* Exon; 0:+1 for pos strand,-1 for neg strand, 1:start; 2:end; 3: Chromosome; 4: start(a.a.); 5: end(a.a.) */ 
int v, w, x, y, z, offset, emax=-1, gflag=0, hit=0, counter=0, chmax=-1, cdslen=0;	/* &&& */
int tcount, tmax=-1;	/* The # of terminal exon regions per gene and its max. */
int ocount=0;           /* #&# The number of outputted exons */
int intstart, intend;   /* #### CDS length encoded by internal exons */
double gw, gx, gy, gz;
void infoout( void ){   /* ### Output the regions of ALL exons with exon # */
        int tx, ty, tz;
        for(ty=0; ty<=exmax; ty++){
                if(exon[0][0]==1)       tx = ty;                /* Positive strand */
                else                    tx = exmax - ty;        /* Negative strand */

                ++ocount;       /* #&# */
                ++tcount;
                w = exon[tx][3];
                fputs(chname[w],fw);    fputs("\t",fw);         /* Chromosome name */
                if(strlen(gstr)<1)      fputs("XXX\t",fw);
                else{
                        fputs(gstr,fw);         fputs("\t",fw);         /* Gene name */
                }
                if(strlen(trname)<1)    fputs("YYY\t",fw);
                else{
                        fputs(trname,fw);       fputs("\t",fw);         /* Transcript name */
                }

                /* &&& Output the total length */
                if(cdslen<0)    fputs("ERROR\t",fw);
                else{
                        strncpy(tempstr,empty,21);
                        sprintf(tempstr,"%d",cdslen);
                        fputs(tempstr,fw);      fputs("\t",fw);
                }
                /* &&& */

                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",exon[tx][1]);      /* The start of the last exon */
                fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",exon[tx][2]);      /* The end of the last exon */
                fputs(tempstr,fw);      fputs("\t",fw);
                if(exon[tx][0]==1)      fputs("+\t",fw);        /* &&& */
                else                    fputs("-\t",fw);        /* &&& */

                /* &&& */
                /* Output the exon number ###  */
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",ty+1);
                fputs(tempstr,fw);      fputs("\t",fw);

                /* #### Output the total # of exons */
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",exmax+1);
                fputs(tempstr,fw);      fputs("\t",fw);

                /* Output the start and end residues */
                if(exon[tx][4]>0){
                        strncpy(tempstr,empty,21);
                        sprintf(tempstr,"%d",exon[tx][4]);       /* The start res. # of the exon */
                        fputs(tempstr,fw);      fputs("\t",fw);
                }
                else    fputs("Noncoding\t",fw);

                if(exon[tx][5]>0){
                        tz = exon[tx][5];
                        if(exon[tx][4]>exon[tx][5])     tz = exon[tx][4];       /* Ad hoc correction */
                        strncpy(tempstr,empty,21);
                        sprintf(tempstr,"%d",tz);               /* The end res. # of the exon */
			fputs(tempstr,fw);      fputs("\t",fw); /* #"# */
                }
		else    fputs("Noncoding\t",fw);                /* #"# */
                /* &&& */

                /* #### Output the range of CDS encoded by internal exons, if any */
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",intstart);
                fputs(tempstr,fw);      fputs("\t",fw);
                strncpy(tempstr,empty,21);
                sprintf(tempstr,"%d",intend);
                fputs(tempstr,fw);      fputs("\n",fw);
                /* #### */

                ty = ty;
        }
        return;
}
void ecread(char rstr[1001]){	/* Parse the exon line to get the sense, start, end, and chromosome # */
	int excdsflag;
	strncpy(tempstr,empty,21);
	for(x=0; x<20; x++){
		if(isspace(rstr[x])!=0)	break;
		tempstr[x] = rstr[x];		/* Chr1, 2, .., 5, C, M */
	}
	w = -1;
	if(chmax>=0){
		for(z=0; z<=chmax; z++){
			if(strcmp(chname[z],tempstr)==0){	/* The chromosome name has been already registered as chname[x] */
				w = z;
				break;
			}
		}
	}
	if(w<0){
		++chmax;
		if(chmax>=50){	printf("The # of chromosome names exceeded 50\n");	exit(1);	}
		w = chmax;
		strcpy(chname[chmax],tempstr);
	}

	offset = x + 1;
	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])==0)	break;
	offset = x + offset;
	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])!=0) break;		/* The end of TAIR10 */
	offset = x + offset;
	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])==0)	break;
	offset = x + offset;
	strncpy(tempstr,empty,21);	/* Either exon or CDS */
	for(x=0; x<20; x++){
		if(isspace(rstr[x+offset])!=0)	break;
		tempstr[x] = rstr[x+offset];
	}
	offset = x + offset;
        excdsflag = 0;
        if(strcmp(tempstr,"exon")==0){
                excdsflag = 1;
                ++exmax;
                if(exmax>=40000){       printf("The # of exons exceeded 40,000 in %s\n",gstr);          exit(1);        }
                exon[exmax][3] = w;     /* Chromosome number */
        }
        else if(strcmp(tempstr,"CDS")==0){
                excdsflag = 2;
                ++cdsmax;
                if(cdsmax>=40000){      printf("The # of CDSs exceeded 40,000 in %s\n",gstr);           exit(1);        }
                cds[cdsmax][3] = w;     /* Chromosome number */
        }
        else    return;

	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])==0)	break;
	offset = x + offset;
	strncpy(tempstr,empty,21);	/* Start */
	for(x=0; x<20; x++){
		if(isspace(rstr[x+offset])!=0)	break;
		tempstr[x] = rstr[x+offset];
	}
	offset = x + offset;
        if(excdsflag==1)        exon[exmax][1] = atoi(tempstr);
        else if(excdsflag==2)   cds[cdsmax][1] = atoi(tempstr);
	
	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])==0)	break;
	offset = x + offset;
	strncpy(tempstr,empty,21);	/* End */
	for(x=0; x<20; x++){
		if(isspace(rstr[x+offset])!=0)	break;
		tempstr[x] = rstr[x+offset];
	}
	offset = x + offset;
        if(excdsflag==1)        exon[exmax][2] = atoi(tempstr);
        else if(excdsflag==2)   cds[cdsmax][2] = atoi(tempstr);

	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])==0)	break;
	offset = x + offset;	/* The beginning of a dot or something */
	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])!=0)	break;
	offset = x + offset;	/* The beginning of a dot or something */

	for(x=0; x<50; x++)	if(isspace(rstr[x+offset])==0)	break;
	offset = x + offset;
	if(rstr[offset]=='+'){
                if(excdsflag==1)        exon[exmax][0] = 1;
                else if(excdsflag==2)   cds[cdsmax][0] = 1;
	}
	else if(rstr[offset]=='-'){
                if(excdsflag==1)        exon[exmax][0] = -1;
                else if(excdsflag==2)   cds[cdsmax][0] = -1;
	}
	else{	printf("What is this character %c in %s",rstr[offset],rstr);	exit(1);	}
	return;
}
void calclen( void ){   /* Calculate the CDS length in a.a. */
        cdslen = 0;
        if(cdsmax<0)    return;
        for(y=0; y<=cdsmax; y++)        cdslen += cds[y][2] - cds[y][1] + 1;
        if(cdslen%3!=0) cdslen = -1;
        else                    cdslen = (cdslen/3) - 1;
        return;
}
void calcintcds( void ){        /* #### Calculate the CDS range encoded by internal exons in aa */
        int tx;
        intstart=0, intend=0;
        if(exmax<=1)    return;         /* One or two exons, No internal exons exist */
        if(cds[0][0]>0){        /* Positive strand */
                for(tx=1; tx<=exmax-1; tx++){
                        if(exon[tx][4]>0){
                                intstart = exon[tx][4];
                                break;
                        }
                }
                for(tx=exmax-1; tx>=0; tx--){
                        if(exon[tx][5]>0){
                                intend = exon[tx][5];
                                break;
                        }
                }
        }
        else{   /* Negative strand */
                for(tx=exmax-1; tx>=0; tx--){
                        if(exon[tx][4]>0){
                                intstart = exon[tx][4];
                                break;
                        }
                }
                for(tx=1; tx<=exmax-1; tx++){
                        if(exon[tx][5]>0){
                                intend = exon[tx][5];
                                break;
                        }
                }
        }
        return;
}       /* #### */
void exonrange( void ){         /* &&&  Determine the range of a.a. residues each exon encodes, if any */
        int cumres=0;
        v = 1;          /* The start res. # of the first exon */
        if(cdsmax<0)    return;
        if(cds[0][0]>0){        /* Positive strand */
                for(y=0; y<=cdsmax; y++){
                        cumres += cds[y][2] - cds[y][1] + 1;    /* The end nucleotide number of the CDS */
                        w = (cumres/3) + 1;     /* The end res. # of the current exon */
                        if(y==cdsmax)   --w;    /* Terminal codon */
                        for(x=0; x<=exmax; x++){        /* Identify the corresponding exon */
                                if(exon[x][1]<=cds[y][1] && cds[y][2]<=exon[x][2]){     /* Hit! */
					/* &&& */
					exon[x][4] = v;         /* The start res. number */
                                     	if(y==cdsmax)	exon[x][5] = cdslen;	/* -1 stop codon */
					else		exon[x][5] = w;         /* The end res. number */
					/* &&& */
                                        break;
                                }
                        }
                        v = w + 1;      /* The start res. of the next exon */
                }
        }
        else{           /* Negative strand */
                for(y=cdsmax; y>=0; y--){
                        cumres += cds[y][2] - cds[y][1] + 1;    /* The end nucleotide number of the CDS */
                        w = (cumres/3) + 1;     /* The end res. # of the current exon */
                        if(y==0)        --w;    /* Terminal codon */
                        for(x=0; x<=exmax; x++){        /* Identify the corresponding exon */
                                if(exon[x][1]<=cds[y][1] && cds[y][2]<=exon[x][2]){     /* Hit! */
					/* &&& */
					exon[x][4] = v;         /* The start res. number */
                                     	if(y==0)	exon[x][5] = cdslen;	/* -1 stop codon */
					else		exon[x][5] = w;         /* The end res. number */
					/* &&& */
                                        v = w + 1;      /* The start res. of the next exon */
                                        break;
                                }
                        }
                }
        }
        return;
}	/* &&& */
int main(void)
{
	char rstr[1001];
        char fnamer[122]="/Users/khomma/work/ensembl/athal0/TAIR10_GFF3_genes.gff";
        char fnamew[122]="/Users/khomma/work/ensembl/athal0/aallexon3M";		/* #"# */

	for(x=0; x<50; x++)	strncpy(chname[x],empty,21);
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs(">Chromosome\tGene\tTranscript\tLength\tExStart\tExEnd\tSense\tExon#\tTotEx#\tStart\tEnd\tIntStart\tIntEnd\n",fw); /* #### */
	strncpy(gstr,empty,1001);
	gflag = 0;	/* Gene flag */
	exmax = -1;	cdsmax = 1;		/* &&& */
	for(x=0; x<40000; x++){
		for(y=0; y<6; y++)	exon[x][y]=0;		/* &&& */
		for(y=0; y<4; y++)	cds[x][y]=0;		/* &&& */
	}
        while((fgets(rstr,1000,fr))!=NULL){
                if(strstr(rstr,"TAIR10")!=NULL && strstr(rstr,"gene")!=NULL && strstr(rstr,"ID=AT")!=NULL){
                        /* Verify that the new gene has actually been shown */
                        ps = strstr(rstr,"TAIR10");
                        pt = strstr(rstr,"gene");
                        /* printf("pt-ps = %ld in %s",pt-ps, rstr);     */
                        if(pt-ps>7){    /* Was a false hit */
                                if(strstr(rstr,"mRNA")!=NULL)   goto transcr;
                                else if(strstr(rstr,"exon")!=NULL || strstr(rstr,"CDS")!=NULL)  goto ecd;
                                else    goto nextr;
                        }

                        if(gflag==1 && exmax>=0 && cdsmax>=0){  /* Output the last variant */
                                ++hit;
                                if(exmax>emax)  emax = exmax;
				calclen();      /* Calculate the length of the CDS */
				/* #### DO NOT eliminae exons that have no CDS region */
                                exonrange();    /* &&& Assign the start and end res. # of encoded protein */
				calcintcds();   /* #### Calculate the range of the CDS encoded by internal exons */
                                if(exmax>=0){
                                        infoout();      /* Output the exon and CDS info. */
                                        if(tcount>tmax)         tmax = tcount;
                                }
                        }
                        gflag = 1;      /* Gene flag */
                        ++counter;
                        if(counter%1000==0)     printf("Now processing gene %d...\n",counter);

                        /* Read in the gene name (e.g., "AT1G01010") */
                        strncpy(gstr,empty,1001);
                        ps = strstr(rstr,"ID=AT");
                        ps += 3;
                        for(x=0; x<1000; x++){
                                if(*ps==';' || isprint(*ps)==0)     break;
                                gstr[x] = *ps;
                                ++ps;
                        }

                        /* Reinitializations */
                        tcount = 0;     /* The # of termnal exon regions per gene */
                        exmax = -1;     cdsmax = -1;
                        for(x=0; x<40000; x++){
				for(y=0; y<6; y++)	exon[x][y]=0;		/* &&& */
				for(y=0; y<4; y++)	cds[x][y]=0;		/* &&& */
                        }
                }
                else if(gflag==1 && strstr(rstr,"TAIR10")!=NULL && strstr(rstr,"mRNA")!=NULL){
transcr:                /* Possibly a new mRNA, i.e., a new set of exons and CDS must be initiated */
                        /* Verify */
                        ps = strstr(rstr,"TAIR10");
                        pt = strstr(rstr,"mRNA");
                        /*      printf("pt-ps = %ld in %s",pt-ps, rstr);        */
                        if(pt-ps>7){
                                if(strstr(rstr,"exon")!=NULL || strstr(rstr,"CDS")!=NULL)       goto ecd;
                                else goto nextr;        /* Was a false hit */
                        }
                        if(gflag==1 && exmax>=0 && cdsmax>=0){  /* Output the last variant */
                                ++hit;
                                if(exmax>emax)  emax = exmax;
				calclen();      /* Calculate the length of the CDS */
				/* #### DO NOT eliminate exons that have no CDS region */
                                exonrange();    /* &&& Assign the start and end res. # of encoded protein */
				calcintcds();   /* #### Calculate the range of the CDS encoded by internal exons */
                                if(exmax>=0){
                                        infoout();      /* Output the exon and CDS info. */
                                        if(tcount>tmax)         tmax = tcount;
                                }
                        }

                        /* Read in the transcript name (e.g., "AT1G01020.2;" */
                        if(strstr(rstr,"ID=AT")!=NULL){
                                ps = strstr(rstr,"ID=AT");
                                ps += 3;
                                for(x=0; x<50; x++){
                                        if(*ps==';' || isprint(*ps)==0)     break;
                                        trname[x] = *ps;
                                        ++ps;
                                }
                        }

                        /* Reinitializations */
                        gflag = gflag;          /* Gene flag unchanged */
                        exmax = -1;     cdsmax = -1;
                        for(x=0; x<40000; x++){
				for(y=0; y<6; y++)	exon[x][y]=0;		/* &&& */
                                for(y=0; y<4; y++)	cds[x][y]=0;
                        }
                }
                else if(gflag==1 && strstr(rstr,"TAIR10")!=NULL && (strstr(rstr,"exon")!=NULL || strstr(rstr,"CDS")!=NULL)){
ecd:                    /* Process the exon or CDS line */
                        /* Verify the hit first: */
                        ps = strstr(rstr,"TAIR10");
                        if(strstr(rstr,"exon")!=NULL)   pt = strstr(rstr,"exon");
                        else                            pt = strstr(rstr,"CDS");
                        /* printf("pt-ps = %ld in %s",pt-ps, rstr);     */
                        if(pt-ps>7)     goto nextr;     /* Was a false hit */
                        ecread(rstr);
                }

nextr:		strncpy(rstr,empty,1001);
        }
	/* Outputting the last gene, if any: */
	if(gflag==1 && exmax>=0 && cdsmax>=0){
		++hit;
		if(exmax>emax)	emax = exmax;
		calclen();      /* Calculate the length of the CDS */
		/* #### DO NOT eliminate exons that have no CDS region */
		exonrange();    /* &&& Assign the start and end res. # of encoded protein */
		calcintcds();   /* #### Calculate the range of the CDS encoded by internal exons */
		if(exmax>=0)	infoout();
	}
	if(tcount>tmax)		tmax = tcount;

	printf("In total %d exons in %d genes have been output\n",ocount,hit);      /* ### */
	printf("The max # of exons is %d, while the max # of exon regions per gene is %d\n",emax+1,tmax);	/* ### */
        fclose(fr);        fclose(fw);
	return 0;
}
