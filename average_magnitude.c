#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "mt19937ar.c"

#define Ns 1000
#define filename "average_magnitude_record_beta3.0_L1e3.dat"
#define sample_max 1000

FILE *fp;

double stress[Ns],Dif[Ns],R;
double stresstmp[Ns],X;
double threshold[Ns],Min=0,Deviding_Factor,S,delta,beta,D;
int broken[Ns],broken_prev[Ns],BF,NNR,NNL,D_NNR,D_NNL,D_TOT,P3;
int nnl[Ns], nnr[Ns],Pmax,NN, Done[Ns];

double threshold_distribution_powerlaw(double beta){
	R=(2*genrand_real2()-1);
	return pow(10,(R*beta));
	}

double break_minimum_stress(int *index){
    	int i=0;
    	double min=1000;
    	for(i=0;i<Ns;i++){
		Dif[i]=threshold[i]-stress[i];
        	if((Dif[i]<min) && (broken[i]!=1)){
            		min=Dif[i];
            		*index=i;
        		}
    		}
    	return(Dif[*index]);
	}

void take_input()
	{
    	int stime;
    	long ltime;
    	ltime=time(NULL);
    	stime=(unsigned)ltime/2;
    	long seedval=stime;
    	init_genrand(seedval);
	}

int main(){
	float Value;
    	int i,j,T_Max,Max=0,T_Fail,DeltaTau,Avalanche,Avalanche1,Number_prev,flag1,flag2;
	double Tmax_Tot=0,Tfail_Tot=0,DeltaTau_Tot=0,P[Ns],Energy1,Energy,Energy_prev,flag1_tot,flag2_tot;
    	double stress_increment,absolute_stress,absolute_stress_tot,Nsurvivors_tot,Number_Tot,amplitude_tot[Ns]; 
    	int breaking_fiber,sample,flag,amplitude[Ns],Q[Ns];
    	int Nsurvivors,Nsurvivors_prev=0,N1,Nsurvivors_prev_stress,Number,Number1;
    	take_input();
	Pmax=1;
	S=(1.0/(2*1.0*Pmax));
	beta=3.0;
	flag1_tot=0.0;
	flag2_tot=0.0;
	for(Number=1;Number<Ns;Number++){amplitude_tot[Number]=0;Q[Number]=0;}
	for(sample=0;sample<sample_max;sample++){
		for(Number=1;Number<Ns;Number++){amplitude[Number]=0;}
		for(i=0;i<Ns;i++){
		        nnl[i]=i-1;
		        nnr[i]=i+1;
		    	}
		nnl[0]=Ns-1;
		nnr[Ns-1]=0;
            	for(i=0;i<Ns;i++){
			threshold[i]=threshold_distribution_powerlaw(beta);
                	broken[i]=0;
			stress[i]=0;
			Done[i]=0;
            		}
		Nsurvivors=Ns;
		Energy=0;
		Energy1=0;
		Avalanche=0;
		Number=0;
		flag1=0;
		flag2=0;
		Energy_prev=0;
		flag=0;
		while(Nsurvivors>0){
			Energy=0.0;
			Number++;
			Q[Number]++;
			Nsurvivors_prev_stress=Nsurvivors;
			stress_increment=break_minimum_stress(&breaking_fiber);
			for(i=0;i<Ns;i++){stress[i]+=stress_increment;}
			for(i=0;i<Ns;i++){
				if((threshold[i]<=stress[i]) && (broken[i]!=1)){
					broken[i]=1;
					Nsurvivors--;
					Energy+=0.5*pow(threshold[i],2);
					if(Nsurvivors==0){break;}
					nnl[nnr[i]]=nnl[i];
					nnr[nnl[i]]=nnr[i];
					}
				}
			while(Nsurvivors_prev!=Nsurvivors && Nsurvivors>0){
				Nsurvivors_prev=Nsurvivors;
				for(i=0;i<Ns;i++){
					if(broken[i]==1 && Done[i]!=1){
						j=i;		
						while(broken[j]==1){j=nnl[j];}
						stress[j]+=S*stress[i];
						j=i;
						while(broken[j]==1){j=nnr[j];}
						stress[j]+=S*stress[i];
						Done[i]=1;
						}
					}
				for(i=0;i<Ns;i++){
                        		if((threshold[i]<=stress[i]) && (broken[i]!=1)){
                            			broken[i]=1;
                            			Nsurvivors--;
						Energy+=0.5*pow(threshold[i],2);
                            			if(Nsurvivors==0){break;}
						nnl[nnr[i]]=nnl[i];
						nnr[nnl[i]]=nnr[i];
						}
					}
				}
			if(Energy>Energy_prev){
				amplitude[Number]=Energy;
				Energy_prev=Energy;
				}
			}
		for(Number=1;Number<Ns;Number++){amplitude_tot[Number]+=amplitude[Number];}
		}
	for(Number=1;Number<Ns;Number++){amplitude_tot[Number]=((1.0*amplitude_tot[Number])/(1.0*sample_max));}
	fp=fopen(filename,"a");
	for(Number=1;Number<Ns;Number++){fprintf(fp,"%d\t%f\n",Number,amplitude_tot[Number]);}
	fclose(fp);
	return 0;
	}
