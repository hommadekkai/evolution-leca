/* minseg1.c	This program identifies and outupts protein segments unique to mouse among human-mouse-rat orthologs	210808
	basic.c -> minsseg1.c: A new program [R2105 II].		  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
int main(void)
{
	FILE *fr, *fw, *fl, *fc;	/* ### */
	char *ps;			/* ### */
	char rstr[122], lstr[122], empty[500]="", tempstr[21], schstr[21], cstr[122];	/* ### */
        char fnamel[100]="/Users/khomma/work/ensembl/hmr/hmrmafft1list";
	char fnamec[122]="/Users/khomma/work/ensembl/mouse/mptcorr1";	 /* ### Correspondece between protein IDs and transcript IDs */
	char fnamer[122]="/Users/khomma/work/ensembl/hmr/hmrmafft1/ENSP00000354687";
	char fnamew[122]="/Users/khomma/work/ensembl/hmr/minsseg1";
	char protname[3][51];
        int w, x, y, z, offset;
	int maffalign[3][40000];	/* 0: Human, 1: Mouse, 2: Rat; is 0 for an a.a., is 1 for a gap */
	int spaa[40000], spseg[5000][2];	/* Mouse-specific a.a.; Mouse-specific segement 0: srart, 1: end */
	int spp, maffmax[3], resnumb, spsmax, counter=0, hitcount=0, uplimit=-1;
	int first, last, totseg=0;	/* &&& */
	double gw, gx, gy, gz;

        if((fl=fopen(fnamel,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamel);      exit(1);        }
        if((fc=fopen(fnamec,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamec);      exit(1);        }	/* ### */
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	fputs(">Human protein\tMouse protein\tMouse transcr.\tProt.len.\tIns.seg\n",fw);
        strncpy(lstr,empty,122);
        while((fgets(lstr,121,fl))!=NULL){
		++counter;
		if(counter%1000==0)	printf("Now processing ortholog group %d...\n",counter);
		strncpy(tempstr,empty,21);
		w = strlen(lstr);
		strncpy(tempstr,lstr,w-1);
		strncpy(fnamer,empty,122);
		strcpy(fnamer,"/Users/khomma/work/ensembl/hmr/hmrmafft1/");
		strcat(fnamer,tempstr);
        	if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
		spp = -1;
		for(x=0; x<3; x++){
			strncpy(protname[x],empty,51);
			maffmax[x] = -1;
			for(y=0; y<40000; y++)	maffalign[x][y] = -1;
		}
		for(y=0; y<40000; y++)	spaa[y] = -1;		/* Mouse specific a.a., signified by 1 */
        	while((fgets(rstr,121,fr))!=NULL){
			if(rstr[0]=='>'){
				++spp;
				if(spp>=3){	printf("More than three seqs in %s\n",fnamer);	exit(1);	}
				offset = 1;
				for(x=0; x<50; x++){
					if(isspace(rstr[x+offset])!=0)		break;
					protname[spp][x] = rstr[x+offset];
				}
			}
			else{	/* Seq. line */
				for(x=0; x<100; x++){
					if(isspace(rstr[x])!=0)		break;
					++maffmax[spp];
					if(maffmax[spp]>=40000){	
						printf("The # of a.a. in sp.%d in %s exceeded 40,000\n",spp,fnamer);
						exit(1);
					}
					w = maffmax[spp];
					if(rstr[x]=='-')	maffalign[spp][w] = 1;	/* Gap */
					else			maffalign[spp][w] = 0;	
				}
			}

                	strncpy(rstr,empty,122);
        	}
        	fclose(fr);        
		if(maffmax[0]!=maffmax[1] || maffmax[0]!=maffmax[2]){	/* Lengths diagree */
			printf("In %s Human %d, Mouse %d, Rat %d residues\n",fnamer,maffmax[0]+1,maffmax[1]+1,maffmax[2]+1);
			goto nextl;
		}

		/* Determine the mouse-specific regions, if any */
		resnumb = 0;		/* Mouse residue number */

		/* &&&  Determine the first and last array # (x) for the commonly aligned segment */
		first = -1; 	last = -1;
		for(x=0; x<=maffmax[0]; x++){
			if(maffalign[0][x]==0 && maffalign[1][x]==0 && maffalign[2][x]==0){
				first = x;
				break;
			}
		}
		for(x=maffmax[0]; x>=0; x--){
			if(maffalign[0][x]==0 && maffalign[1][x]==0 && maffalign[2][x]==0){
				last = x;
				break;
			}
		}
		if(first<0 || last<0)	goto nextl;	/* No commonly aligned segment */
		/* &&& */

		for(x=0; x<=maffmax[1]; x++){
			if(maffalign[1][x]==1)	continue;	/* Gap in mouse, neglected */
			else if(maffalign[1][x]==0){
				++resnumb;
				if(resnumb>=40000){	printf("The # of mouse residues exceeded 40,000 at #100\n");	exit(1); 	}
				if(x>=first && x<=last && maffalign[0][x]==1 && maffalign[2][x]==1){	/* &&& */
					/* HIT! This residue is a mouse specific residue */
					spaa[resnumb] = 1;
				}
				else	spaa[resnumb] = 0;
			}
		}
		spsmax = -1;
		for(x=0; x<5000; x++){
			for(y=0; y<2; y++)	spseg[x][y] =0; 
		}
		for(x=1; x<=resnumb; x++){
			if(spaa[x]==1){	/* Find the endpoint of the mouse-specific region */
				for(y=x+1; y<=resnumb; y++)	if(spaa[y]==0)	break;
				++spsmax; 
				if(spsmax>=5000){	printf("The # of mouse-specific segments in %s exceeded 5,000\n",fnamer); exit(1); }
				spseg[spsmax][0] = x;
				spseg[spsmax][1] = y-1;
				x = y;		/* Shifting the res. number to examine */
				if(x>resnumb)	break;
			}
		}
		if(spsmax>=0){	/* Output the mouse specific regions */
			if(spsmax>uplimit)	uplimit = spsmax;
			++hitcount;
			fputs(protname[0],fw);	fputc(' ',fw);		/* Human protein ID */	
			fputs(protname[1],fw);	fputc(' ',fw);		/* Mouse protein ID */

			/* ### Find the transript ID that corresponds to the protein ID ### */
			fseek(fc,0L,SEEK_SET);
			strncpy(cstr,empty,122);
			strncpy(schstr,empty,21);
			for(x=0; x<20; x++){
				if(protname[1][x]=='_')	break;
				else if(isspace(protname[1][x])!=0)	break;
				else schstr[x] = protname[1][x];
			}
        		while((fgets(cstr,121,fc))!=NULL){
				if(strstr(cstr,schstr)!=NULL && strstr(cstr,"ENSMUST")!=NULL){	/* Hit */
					strncpy(tempstr,empty,21);
					ps = strstr(cstr,"ENSMUST");
					for(x=0; x<20; x++){
						if(isspace(*ps)!=0)	break;
						tempstr[x] = *ps;
						++ps;
					}
					fputs(tempstr,fw);
					fputs(" ",fw);
					goto lout;
				}
				strncpy(cstr,empty,122);
			}
lout:			/* ### Output the length ### */

			fputs("[LEN=",fw);
			strncpy(tempstr,empty,21);
			sprintf(tempstr,"%d",resnumb);
			fputs(tempstr,fw);
			fputs("] ",fw);
			for(x=0; x<=spsmax; x++){
				++totseg;	/* ### */
				strncpy(tempstr,empty,21);
				sprintf(tempstr,"%d",spseg[x][0]);
				fputs(tempstr,fw);
				fputs("->",fw);
				strncpy(tempstr,empty,21);
				sprintf(tempstr,"%d",spseg[x][1]);
				fputs(tempstr,fw);
				if(x<spsmax)	fputc('|',fw);
				else		fputs("\n",fw);
			}

		}

nextl:		strncpy(lstr,empty,122);
	}
        fclose(fl);	fclose(fc);	fclose(fw);
	printf("%d proteins out of %d had mouse-specific regions, with the max. of %d regions in a protein\n",hitcount,counter,uplimit+1);
	printf("A total of %d inserted segments exist\n",totseg);
	return 0;
}
