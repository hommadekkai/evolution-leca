#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
FILE *fr, *fw;
char rstr[201], empty[510]="", tempstr[21], listread[101];
char fnamer[122]="/Users/khomma/work/ensembl/human/hlongingo4";
char fnamew[122]="/Users/khomma/work/ensembl/total/commonlongingo4";
char spstr[4][50]={"human/h","mouse/m","athal0/a","spomb/f"};
int  w, x, y, z, offset, repeat, classnumb, highmax=-1, lowmax=-1, highhit=0, lowhit=0;
double gw, gx, gy, gz;
struct ingo{
        int  gonumb;	/* GO # */
        char godescr[101];	/* GO description */
	int  spnumb;	/* The # of appearances in the four spp */
};
struct ingo highlist[10000];
struct ingo lowlist[10000];
void listparse(int i)
{
        int j=0, k;
        offset = 0;
        strncpy(listread,empty,101);
        while(j<i){
                j++;
                strncpy(listread,empty,101);
                for(k=0; k<100; k++){
                        if(i!=6 && isspace(rstr[k+offset])!=0)  break;
			else if(isprint(rstr[k+offset])==0)	break;
                        listread[k] = rstr[k+offset];
                }
                offset = k + offset;
                for(k=0; k<50; k++){
		    	if(isspace(rstr[k+offset])==0)  break;
		}
                offset = k + offset;
        }
        return;
}
int main(void)
{

        if((fw=fopen(fnamew,"w")) == NULL){	printf("Cannot open the output file %s\n", fnamew);     exit(1);        }
	for(repeat=0; repeat<4; repeat++){
		strncpy(fnamer,empty,122);
		strcpy(fnamer,"/Users/khomma/work/ensembl/");
		strcat(fnamer,spstr[repeat]);
		strcat(fnamer,"longingo4");
        	if((fr=fopen(fnamer,"r")) == NULL){	printf("Cannot open the input file %s\n", fnamer);	exit(1);	}
		classnumb = -1;
        	strncpy(rstr,empty,201);
        	while((fgets(rstr,200,fr))!=NULL){
			if(rstr[0]=='>')	++classnumb;	/* 0~2: Higher freq, 3~5: Lower freq */
			else if(strncmp(rstr,"GO #",4)!=0){
				listparse(1);		/* Get the GO number */
				w = atoi(listread);

				if(classnumb<=2){
					if(highmax>=0){		/* Check if the GO number has been registered */
						for(x=0; x<=highmax; x++){
							if(highlist[x].gonumb == w){	/* Aready registered */
								++highlist[x].spnumb;
								goto nextr;
							}
						}
					}
					/* GO number has not been registered */
					++highmax;
					if(highmax>=10000){  printf("The # of GO #s in the high category exceeded 10,000\n");   exit(1); } 
					highlist[highmax].gonumb = w;
					listparse(6);
					strcpy(highlist[highmax].godescr,listread);	
					++highlist[highmax].spnumb;
				}
				else{
					if(lowmax>=0){		/* Check if the GO number has been registered */
						for(x=0; x<=lowmax; x++){
							if(lowlist[x].gonumb == w){	/* Aready registered */
								++lowlist[x].spnumb;
								goto nextr;
							}
						}
					}
					/* GO number has not been registered */
					++lowmax;
					if(lowmax>=10000){  printf("The # of GO #s in the low category exceeded 10,000\n");   exit(1); } 
					lowlist[lowmax].gonumb = w;
					listparse(6);
					strcpy(lowlist[lowmax].godescr,listread);	
					++lowlist[lowmax].spnumb;
				}

			}
                	
nextr:			strncpy(rstr,empty,201);
        	}
        	fclose(fr);
	}
	/* Process the two categories */
	printf("The high and low categories respectively have %d and %d different GO numbers\n",highmax+1,lowmax+1);
	fputs(">Commonly occurred in the high freq. category\n",fw);
	fputs("GO #\tGO description\n",fw);
	for(x=0; x<=highmax; x++){
		if(highlist[x].spnumb>=4){	/* The GO # appeared in all four lists */
			++highhit;
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",highlist[x].gonumb);
			fputs(tempstr,fw);		fputs("\t",fw);	
			fputs(highlist[x].godescr,fw);	fputs("\n",fw);
		}
	}

	fputs(">Commonly occurred in the low freq. category\n",fw);
	fputs("GO #\tGO description\n",fw);
	for(x=0; x<=lowmax; x++){
		if(lowlist[x].spnumb>=4){	/* The GO # appeared in all four lists */
			++lowhit;
			strncpy(tempstr,empty,21);	sprintf(tempstr,"%d",lowlist[x].gonumb);
			fputs(tempstr,fw);		fputs("\t",fw);	
			fputs(lowlist[x].godescr,fw);	fputs("\n",fw);
		}
	}
	printf("%d and %d GO numbered commonly occurred in the high and low categories, respectively\n",highhit,lowhit);
        fclose(fw);
	return 0;
}
