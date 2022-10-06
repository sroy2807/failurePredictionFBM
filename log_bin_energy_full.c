#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/*#define BIN (sqrt(2.0))*/
#define BIN 1.1
#define LOGBIN (log(BIN))
#define MAX 1000 /*number of log bins*/
#define MAX1 100000



double count[MAX],ind[MAX],mass[MAX];
char infile[100];

int main()
{
	int i,x,y,k,j,total;
int x1;
    double a,b,min,max,sum;
    double x2;
	FILE *fp;
	void lat_init();
	

	sprintf(infile,"energy_full_spectrum.dat");
    fp=fopen(infile,"r");
    for (i=0;i<MAX;i++){count[i]=0;}

max=0.0;
min=1000000; 

	while(fscanf(fp,"%d%lf",&x1,&x2)!=EOF){

	if(x2<min && x2!=0){min=x2;}a=(min*(pow(BIN,0.5)));
	if(x2>0){k=floor((0.5)+((log(x2/a))/(LOGBIN)));if(k<MAX){count[k]++;}}
	
    	}
	fclose(fp);

	total=0;
	for(k=0;k<MAX;k++){total+=count[k];}
	for(k=0;k<MAX;k++){ind[k]=((exp(log(a)+((k+0.5)*LOGBIN))))-((exp(log(a)+((k-0.5)*LOGBIN))));}
	//for(k=0;k<MAX;k++){printf("%d\t%.8e\n",k,ind[k]);}
	
	fp=fopen("logbin_energy_energy_full_spectrum.dat","w");
	for(k=0;k<MAX;k++){if(count[k]>0){fprintf(fp,"%.8e\t%.8e\n",(a*(pow(BIN,k))),(count[k]/(ind[k]*total)));}}
	fclose(fp);
}

