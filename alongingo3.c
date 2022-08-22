#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw, *fl;
char lstr[5001], empty[5010]="", tempstr[21];
char fnamer[122]="/Users/khomma/work/ensembl/athal0/alongexgenes4";
char fnamel[100]="/Users/khomma/work/uniprot/uniprot_sprot_211117.dat";
char fnamew[122]="/Users/khomma/work/ensembl/athal0/alongingo3";
char binstr[11][20]={"All","-120","-240","-360","-480","-600","-720","-840","-960","961bp-","Single-exon"};
char listread[2001];
char godescr[30000][101], fnstr[101];
char chistr[3][20]={"p<0.05","p<0.01","p<0.001"};
int  p, q, v, w, x, y, z, offset, maxp=-1, bincount=0, counter=0;
int  gonumb[30000], gomax=-1, binnumb=-1, gocount[30000][11], hitcount[11];
int  signifcount[2][3];		/* The # of GOs that passed the threshold, 0: Higher freq. 1: Lower freq. */
int  oksp=0;
double gu, gv, gw, gx, gy, gz, obs11, obs12, obs21, obs22, exp11, exp12, exp21, exp22, chisq, totcount;
double chiinv[3]={3.841, 6.636, 10.828};
struct protlist{
        char geneid[21];
	int  binnumb;
};
struct protlist list[30000];
void listparse(int i)
{
        int j=0, k;
        offset = 0;
        strncpy(listread,empty,2001);
        while(j<i){
                j++;
                strncpy(listread,empty,2001);
                for(k=0; k<20000; k++){
                        if(isspace(lstr[k+offset])!=0)  break;
                        listread[k] = lstr[k+offset];
                }
                offset = k + offset;
                for(k=0; k<50; k++)     if(isspace(lstr[k+offset])==0)  break;
                offset = k + offset;
        }
        return;
}
int main(void)
{
        if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
        if((fl=fopen(fnamel,"r")) == NULL){     printf("Cannot open the input file %s\n", fnamel);      exit(1);        }
        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(x=0; x<30000; x++){
		strncpy(godescr[x],empty,101);
		gonumb[x] = -1;
		for(y=0; y<11; y++)	gocount[x][y] = 0;	/* # of GO in each length bin y */
	}
	for(y=0; y<11; y++)	hitcount[y] = 0;	/* # of proteins in each length bin y */
	for(x=0; x<2; x++){
		for(y=0; y<3; y++)	signifcount[x][y] = 0;
	}

	/* Read in the list of proteins classifed by exon length */
	bincount = 0;
	strncpy(lstr,empty,5001);
        while((fgets(lstr,121,fr))!=NULL){
                if(strstr(lstr,"Gene name")!=NULL)      goto nextr;
                else if(lstr[0]=='>'){
                        strncpy(tempstr,empty,21);
                        offset = 4;
                        for(x=0; x<20; x++){
                                if(isdigit(lstr[x+offset])==0)  break;
                                tempstr[x] = lstr[x+offset];
                        }
                        bincount = atoi(tempstr);
                        goto nextr;
                }
                ++maxp;
                if(maxp>=35000){        printf("The # of proteins exceeded 35,000\n");  exit(1);        }
                listparse(1);
                strcpy(list[maxp].geneid,listread);
                list[maxp].binnumb = bincount;

nextr:          strncpy(lstr,empty,5001);
        }
	printf("%d proteins classified\n",maxp+1);

	/* Read in the GO # annotations to proteins */
	strncpy(lstr,empty,5001);	oksp = 0;
        while((fgets(lstr,5000,fl))!=NULL){
		if(strncmp(lstr,"ID",2)==0 && strstr(lstr,"_ARATH")!=NULL){
			oksp = 1;
			++counter;
			if(counter%2000==0)	printf("Now processing protein %d...\n",counter);
		}
		else if(oksp==1 && strstr(lstr,"DR   EnsemblPlants;")!=NULL && binnumb<0){
			/* Try to identify the int. max exon length bin # */
			for(x=0; x<=maxp; x++){
				if(strstr(lstr,list[x].geneid)!=NULL){	/* Hit */
					binnumb = list[x].binnumb;
					++hitcount[binnumb];
					break;
				}
			}
		}
		else if(oksp==1 && strstr(lstr,"DR   GO;")!=NULL && binnumb>0){
			offset = 12;
			strncpy(tempstr,empty,21);
			for(x=0; x<20; x++){
				if(lstr[x+offset]==';')		break;
				tempstr[x] = lstr[x+offset];
			}
			offset = x + offset + 2;
			strncpy(fnstr,empty,101);
			for(x=0; x<100; x++){
				if(lstr[x+offset]==';' || isprint(lstr[x+offset])==0)	break;
				fnstr[x] = lstr[x+offset];
			}
			w = atoi(tempstr);
			v = -1;
			if(gomax>=0){
				for(x=0; x<=gomax; x++){
					if(gonumb[x]==w){
						v = x;		/* GO # already registered as x */
						++gocount[x][binnumb];
						break;
					}
				}
			}
			if(v<0){	/* Not registered yet */
				++gomax;
				if(gomax>=20000){	printf("The GO numbers used exceeded 20,000!\n");	exit(1);	}
				gonumb[gomax] = w;
				strcpy(godescr[gomax],fnstr);
			}
		}
		else if(strncmp(lstr,"//",2)==0){
			oksp = 0;
			binnumb = -1;
		}
		strncpy(lstr,empty,5001);
	}
	/* Sum over exon length bins */
	for(x=0; x<=gomax; x++){
		for(y=1; y<11; y++)	gocount[x][0] += gocount[x][y];
	}
	for(y=1; y<11; y++)	hitcount[0] += hitcount[y];
        printf("%d out of %d had their GO #s classified. %d different GO #s exist\n",hitcount[0],counter,gomax+1);
        for(q=0; q<2; q++){     /* 0: Higher freq. in bins 7-9; 1: Lower freq. in bins 7-9 */
                for(p=2; p>=0; p--){    /* Diff. threshold of chi-sq values */
                        if(q==0)        fputs(">Higher freq. in bins 7-9\t",fw);
                        else            fputs(">Lower freq. in bins 7-9\t",fw);
                        fputs(chistr[p],fw);    fputs("\n",fw);
                        fputs("GO #\tChiSq\tFreq in bins 7-9\tFreq in other bins\tRatio\tGO description\n",fw);
                        for(x=0; x<=gomax; x++){
                                obs11 = (double)(gocount[x][7]+gocount[x][8]+gocount[x][9]);
                                obs21 = (double)(hitcount[7]+hitcount[8]+hitcount[9]) - obs11;
                                obs12 = (double)gocount[x][0] - obs11;
                                obs22 = (double)(hitcount[0] -gocount[x][0] +obs11 -(hitcount[7]+hitcount[8]+hitcount[9]));
                                totcount = (double)hitcount[0];
                                if(abs(obs11+obs12+obs21+obs22-totcount)>0.1){  printf("Error at #100\n");      exit(1);        }

                                gu = 100.*obs11/(obs11+obs21);          /* Freq. of occurence in bins 7-9 */
                                gv = 100.*obs12/(obs12+obs22);          /* Freq. of occurence in the other bins */
                                if(q==0 && gu<gv)       goto nextx;
                                else if(q==1 && gu>gv)  goto nextx;

                                exp11 = totcount*((obs11 + obs12)/totcount)*((obs11 + obs21)/totcount);
                                exp21 = totcount*((obs21 + obs22)/totcount)*((obs11 + obs21)/totcount);
                                exp12 = totcount*((obs11 + obs12)/totcount)*((obs12 + obs22)/totcount);
                                exp22 = totcount*((obs21 + obs22)/totcount)*((obs12 + obs22)/totcount);
                                if(abs(exp11+exp12-obs11-obs12)>0.1){   printf("Error at #200\n");      exit(1);        }
                                if(abs(exp21+exp22-obs21-obs22)>0.1){   printf("Error at #300\n");      exit(1);        }
                                if(abs(exp11+exp21-obs11-obs21)>0.1){   printf("Error at #400\n");      exit(1);        }
                                if(abs(exp12+exp22-obs12-obs22)>0.1){   printf("Error at #500\n");      exit(1);        }

                                if(exp11>0.00001 && exp12>0.00001 && exp21>0.00001 && exp22>0.00001){
                                        gw = (obs11 -exp11)*(obs11 -exp11)/exp11;
                                        gx = (obs12 -exp12)*(obs12 -exp12)/exp12;
                                        gy = (obs21 -exp21)*(obs21 -exp21)/exp21;
                                        gz = (obs22 -exp22)*(obs22 -exp22)/exp22;
                                        chisq = gw +gx +gy +gz;
                                        if((p==2 && chisq>=chiinv[p]) || (p<=1 && chisq>=chiinv[p] && chisq<chiinv[p+1])){
                                                /* Passed the threshold */
                                                ++signifcount[q][p];

                                                strncpy(tempstr,empty,21);      sprintf(tempstr,"%7d",gonumb[x]);
                                                for(z=0; z<7; z++){     /* Convert 14059 to 0014059 for example */
                                                        if(isdigit(tempstr[z])==0)      tempstr[z] = '0';
                                                }
                                                fputs(tempstr,fw);              fputs("\t",fw);

                                                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.2f",chisq);
                                                fputs(tempstr,fw);              fputs("\t",fw);

                                                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gu);
                                                fputs(tempstr,fw);              fputs("\t",fw);

                                                strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gv);
                                                fputs(tempstr,fw);              fputs("\t",fw);

                                                if(gv<=0.0001)  fputs("9999.99\t",fw);
                                                else{
                                                        strncpy(tempstr,empty,21);      sprintf(tempstr,"%.3f",gu/gv);
                                                        fputs(tempstr,fw);              fputs("\t",fw);
                                                }

                                                fputs(godescr[x],fw);           fputs("\n",fw);
                                        }
                                }
nextx:                          x = x;
                        }
                }
        }
        printf("p threshold\tHigher freq. in bins 7-9\tLower freq. in bins 7-9\n");
        for(p=2; p>=0; p--){
                printf("%s\t",chistr[p]);
                for(q=0; q<2; q++){
                        if(q==0)        printf("%d\t",signifcount[q][p]);
                        else            printf("%d\n",signifcount[q][p]);
                }
        }
        fclose(fr);        fclose(fw);
        fclose(fl);
	return 0;
}
