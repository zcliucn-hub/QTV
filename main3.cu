// Cardiac model ToR-ORd

// New adding I_NaL in HRd Model
// times or  multiple PARAs
// CaL, Ks, Kr and NaL
// Paras Plane of x_D and y_D (you can set as you want) 

/*
        Cardiac model ToR-ORd
%     Copyright (C) 2019 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>
//    TorORd model
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector> 
#include <fstream>
#include <random>
#include <iostream>
#include <algorithm>

#pragma mark PROBLEM CONSTANTS
////////////////////////////////////////////////////////////////////////////////
// Problem constants
////////////////////////////////////////////////////////////////////////////////

#define 		Device_NUM		0					//Select GPU or not

// Paras Set
//#define 	 	Long_QT			3
#define 	 	Para_short		20.0				//2.0 for LQT1 and 20.0 for LQT2
#define 	 	Para_long		0.0
#define 	 	v_th			0.005
#define 	 	step_adp		10

// Euler constants
#define 		dt            	0.01
#define 		STEP          	932165 			
#define 		TEMPO         	0

//Tissue Settings
#define 		Length       	101					//The Number of the Total Cells in Tissue (Length-2)
#define 		EAD_length		101					//The Number of the Long APD Cells (EAD_length-1)
#define 		NDIM      		46					//The dimension of single cell dynamics (NDIM variables)
#define 		ETM      		100					//The Max Number of EAD
#define 		N_gradient   	0					//Gradient from long to short 

// GPU blocks and threads
#define 		size_x			1
#define 		size_y			1
#define  		Blocknum     	(size_x*size_y)  	//The number of Paralleled Processes

// Parameter Domain
#define  		y_min          	1.0
#define  		y_max          	1.0
//#define  		y_per          	((y_max-y_min)/(size_y-1.0))
#define  		y_per          	1.0
#define  		x_D 			Para1
#define  		y_D           	(y_min+(Num_j/size_x)*y_per)
#define  		y_D1           	(y_min+(idy/size_x)*y_per)

// Format of Data Recording
#define 		Vshift       	10.0
#define 		tshift       	10.0
#define 		Lfactor      	8

// Simulation options
#define 		M_STIM        	53.0
#define 		STIM_OCT      	1.0 		
#define 		vdrg          	-72.0 		
#define 		vdrg2          	-50.0 
#define 		BCL           	4000.0
#define 		N_PREPACE     	20
#define 		Interval      	1500.0

// Diffusion constant
#define 		DIFF            0.0005
#define 		DX              0.015

#define 		L				Length 
//#define 			v1				v

//000
//#define         PCL1            932.165
//#define         RR_SD          65.36

//002
//#define         PCL1            654.868
//#define         RR_SD          67.93

//003
//#define         PCL1            668.11
//#define         RR_SD          90.73

//005
//#define         PCL1            766.186
//#define         RR_SD          54.64

//413
//#define         PCL1            721.618
//#define         RR_SD          91.125

//414
//#define         PCL1            691.891
//#define         RR_SD          80.592

//415
//#define         PCL1            743.407
//#define         RR_SD          123.624

//4000
//#define         PCL1            420.063
//#define         RR_SD          17.554

//4002
#define         PCL1            474.227
#define         RR_SD          36.856

#define         PCL2            500.0

#define         ye              0.2
#define         xe              1.8
#define         ECG_th          0.01

#define         burn            500
#define         record          1000
#define         rr_total        (burn+record)



#include "kernel.cu"
//#include "ECG.cu"

////////////////////////////////////////////////////////////////////////////////
//! update the variables
////////////////////////////////////////////////////////////////////////////////
__global__ void update(double *state_var,double *state_deri,double *vnew,
		double Fsti,int stim_site, int clamp, double T_time,
		double *APD, int *APD_gate,
		int *EAD_gate, int *EAD_times, double *EAD_amp,double *EAD_min,double *EAD_max,
		int *EAD_count,double *v_takeoff, double *EAD_period, int *EAD_times_pause, double *ECG, double Para1, double Para3)
{     	
		int i,Num,Long_QT;	
		double Itot,GKS_C,GCAL_C,GKR_C,GNAL_opt,v1;
		int idx=threadIdx.x, idy=blockIdx.x; 
		double ECGn,dVdx,r2,r3;

        Long_QT=Para3;
	
		if(idy<Blocknum)
		{
				if(idx>0&&idx<Length-1)
				{ 
						GCAL_C=1.0;
						GNAL_opt=0.0;
						GKS_C=1.0;
						GKR_C=1.0;
						if(Long_QT==1)
						{
								GKS_C=0.0;					//LQT 1
								if(idx<EAD_length) 
								{
										GKR_C=Para1;
								}
								else
								{
										GKR_C=Para_short; 
								}
								GNAL_opt=0.0;
						}
						else if(Long_QT==2)
						{
								if(idx<EAD_length) 
								{
										GKS_C=Para1;
								}
								else
								{
										GKS_C=Para_short;
								}
								GKR_C=0.3;					//LQT 2
								GNAL_opt=0.0;
						}
						else if(Long_QT==3)
						{
								if(idx<EAD_length) 
								{
										GKS_C=Para1;
								}
								else
								{
										GKS_C=Para_short;
								}
								GKR_C=1.0;	
								GNAL_opt=10.0;				//LQT 3
						}
						else if(Long_QT==4)
						{
								if(idx<EAD_length) 
								{
										GKS_C=Para_long;
								}
								else
								{
										GKS_C=Para1;
								}
								GKR_C=0.0;					//LQT 2
								GNAL_opt=0.0;
						}
						else if(Long_QT==5)
						{
								GKS_C=0.0;					//LQT 1
								if(idx<EAD_length) 
								{
										GKR_C=Para_long;
								}
								else
								{
										GKR_C=Para1;
								}
								GNAL_opt=0.0;
						}
						
						v1=v;
						
						if(cass<v_th)
						{
								Itot=comp_current(state_var,state_deri,idx,idy,GCAL_C,GKS_C,GKR_C,GNAL_opt);
						
								if(idx>=0&&idx<=2)
										dv = -Itot+Fsti;						
								else
										dv = -Itot;
								//Stimulus
								
								if(idx==1)
										vnew[idx+Length*idy]=state_var[pos(1,idx,idy)]+ dv*dt +
												DIFF*dt*(state_var[pos(1,idx+1,idy)]+state_var[pos(1,idx,idy)]-2.0*state_var[pos(1,idx,idy)])/(DX*DX); 
								else if(idx==Length-2)
										vnew[idx+Length*idy]=state_var[pos(1,idx,idy)]+ dv*dt +
												DIFF*dt*(state_var[pos(1,idx,idy)]+state_var[pos(1,idx-1,idy)]-2.0*state_var[pos(1,idx,idy)])/(DX*DX); 
								else
										vnew[idx+Length*idy]=state_var[pos(1,idx,idy)]+ dv*dt +
												DIFF*dt*(state_var[pos(1,idx+1,idy)]+state_var[pos(1,idx-1,idy)]-2.0*state_var[pos(1,idx,idy)])/(DX*DX); 
								//Couplings
							
								for (i = 2; i <NDIM; i++)
										state_var[pos(i,idx,idy)] = state_var[pos(i,idx,idy)]+ state_deri[pos(i,idx,idy)]*dt;
								//Update
						}
						else
						{		for(Num=0;Num<step_adp;Num++)
								{
										Itot=comp_current(state_var,state_deri,idx,idy,GCAL_C,GKS_C,GKR_C,GNAL_opt);
								
										if(idx>=0&&idx<=Length-2)
												dv = -Itot+Fsti;						
										else
												dv = -Itot;
										//Stimulus
										
										if(idx==1)
												vnew[idx+Length*idy]=state_var[pos(1,idx,idy)]+ dv*(dt/step_adp) +
														DIFF*(dt/step_adp)*(state_var[pos(1,idx+1,idy)]+state_var[pos(1,idx,idy)]-2.0*state_var[pos(1,idx,idy)])/(DX*DX); 
										else if(idx==Length-2)
												vnew[idx+Length*idy]=state_var[pos(1,idx,idy)]+ dv*(dt/step_adp) +
														DIFF*(dt/step_adp)*(state_var[pos(1,idx,idy)]+state_var[pos(1,idx-1,idy)]-2.0*state_var[pos(1,idx,idy)])/(DX*DX); 
										else
												vnew[idx+Length*idy]=state_var[pos(1,idx,idy)]+ dv*(dt/step_adp) +
														DIFF*(dt/step_adp)*(state_var[pos(1,idx+1,idy)]+state_var[pos(1,idx-1,idy)]-2.0*state_var[pos(1,idx,idy)])/(DX*DX); 
										//Couplings
										
										for (i = 2; i <NDIM; i++)
												state_var[pos(i,idx,idy)] = state_var[pos(i,idx,idy)]+ state_deri[pos(i,idx,idy)]*(dt/step_adp);
										
										v = vnew[idx+Length*idy];
										//Update
								}
						}
				}
								
				//The APDs and PVCs
				//Recording APDs 
				if(APD_gate[idx+Length*idy]==0&&v1<vdrg2&&vnew[idx+Length*idy]>vdrg2)  	//UP
				{
						APD_gate[idx+Length*idy]=1;
						if(idx==Length-2) EAD_times_pause[idy]=EAD_times_pause[idy]+1; //count EAD times
				}
				if(APD_gate[idx+Length*idy]==1&&v1>vdrg&&vnew[idx+Length*idy]<vdrg) APD_gate[idx+Length*idy]=0; //DOWN
				if(APD_gate[idx+Length*idy]==1) APD[idx+Length*idy]=APD[idx+Length*idy]+dt; 
				//vdrg=-62.0 //or -72.0
				
				//----------------------------------------------------------------------------------- 

				if(T_time>0&&fmod(T_time,1.5)<1e-6)
				{
			    ECGn=0;
				for (int cell_j1=1;cell_j1<Length-2;cell_j1=cell_j1+1)
				{
					dVdx=(state_var[pos(1,cell_j1-1,idy)]-state_var[pos(1,cell_j1+1,idy)])/(2.0*DX);
					r2=(cell_j1*DX-xe)*(cell_j1*DX-xe)+ye*ye;
					r3=pow(r2,1.5);
					ECGn=ECGn+DIFF*(dVdx*(xe-cell_j1*DX)/r3)*DX;
				}

                ECG[idy]=ECGn;
				}

				__syncthreads();
				v=vnew[idx+Length*idy];
		}
}


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int pace_count=0;
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) 
{
		double const Para1=atof(argv[1]);
		int const Para2=atoi(argv[2]);
        int const Para3=atoi(argv[3]);
		if(Device_NUM==1)
		{
				// ------------------- Input arguments, selecting the GPU ----------------   
				int CudaDevice=0;	
				if(argc>=3) 
				{
					int Device=atoi(argv[3]);	// Argument #1: cuda device (default: 0)
					if(Device>=0)
					CudaDevice=Device;	
				}
				cudaSetDevice(CudaDevice);
				//-------------------------------------------------------------------------
		}
		
		FILE *fileopen1,*fileopen2,*fileopen3;
		char fileSpec1[100],fileSpec2[100],fileSpec3[100];

		double rr;
        std::vector<double> rr_host;

		if(Para2==0)
		{
			snprintf(fileSpec1, 100, "Tissue_V_sim_mr0600_SD030_beta000_1.dat");
		    snprintf(fileSpec2, 100, "ECG_sim_mr0600_SD030_beta000_1.dat");
		    snprintf(fileSpec3, 100, "QT_sim_mr0600_SD030_beta000_1.dat");
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("sim_mr0600_SD030_beta000_test_1.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}
		else if(Para2==1)
		{
			snprintf(fileSpec1, 100, "Tissue_V_sim_mr0600_SD030_beta000_2.dat");
		    snprintf(fileSpec2, 100, "ECG_sim_mr0600_SD030_beta000_2.dat");
		    snprintf(fileSpec3, 100, "QT_sim_mr0600_SD030_beta000_2.dat");
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("sim_mr0600_SD030_beta000_test_2.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}
		else if(Para2==2)
		{
			snprintf(fileSpec1, 100, "Tissue_V_sim_mr0600_SD030_beta000_3.dat");
		    snprintf(fileSpec2, 100, "ECG_sim_mr0600_SD030_beta000_3.dat");
		    snprintf(fileSpec3, 100, "QT_sim_mr0600_SD030_beta000_3.dat");
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("sim_mr0600_SD030_beta000_test_3.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}
		else if(Para2==3)
		{
			snprintf(fileSpec1, 100, "Tissue_V_sim_mr0600_SD030_beta000_4.dat");
		    snprintf(fileSpec2, 100, "ECG_sim_mr0600_SD030_beta000_4.dat");
		    snprintf(fileSpec3, 100, "QT_sim_mr0600_SD030_beta000_4.dat");
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("sim_mr0600_SD030_beta000_test_4.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}
		else if(Para2==4)
		{
			snprintf(fileSpec1, 100, "Tissue_V_sim_mr0600_SD030_beta000_0.dat");
		    snprintf(fileSpec2, 100, "ECG_sim_mr0600_SD030_beta000_0.dat");
		    snprintf(fileSpec3, 100, "QT_sim_mr0600_SD030_beta000_0.dat");
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("sim_mr0600_SD030_beta000_test_0.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}

		/*if(Para2==0)
		{
			snprintf(fileSpec1, 100, "Tissue_V_real_4035_l%d.dat",Para3);
		    snprintf(fileSpec2, 100, "ECG_real_4035_l%d.dat",Para3);
		    snprintf(fileSpec3, 100, "QT_real_4035_l%d.dat",Para3);
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("RR_real_4035.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}
		else if(Para2==1)
		{
			snprintf(fileSpec1, 100, "Tissue_V_shuf_4035_l%d.dat",Para3);
		    snprintf(fileSpec2, 100, "ECG_shuf_4035_l%d.dat",Para3);
		    snprintf(fileSpec3, 100, "QT_shuf_4035_l%d.dat",Para3);
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("RR_shuf_4035.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }

			}
		else if(Para2==2)
		{
			snprintf(fileSpec1, 100, "Tissue_V_simu_4035_l%d.dat",Para3);
		    snprintf(fileSpec2, 100, "ECG_simu_4035_l%d.dat",Para3);
		    snprintf(fileSpec3, 100, "QT_simu_4035_l%d.dat",Para3);
		    fileopen1=fopen(fileSpec1,"w");
		    fileopen2=fopen(fileSpec2,"w");
		    fileopen3=fopen(fileSpec3,"w");

		    std::ifstream fin("RR_simu_4035.dat");
			while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}
		else
		{
		snprintf(fileSpec1, 100, "Tissue_V_gaus_4035_l%d.dat",Para3);
		snprintf(fileSpec2, 100, "ECG_gaus_4035_l%d.dat",Para3);
		snprintf(fileSpec3, 100, "QT_gaus_4035_l%d.dat",Para3);
		fileopen1=fopen(fileSpec1,"w");
		fileopen2=fopen(fileSpec2,"w");
		fileopen3=fopen(fileSpec3,"w");

		std::ifstream fin("RR_gaus_4035.dat");
		while (fin >> rr) {
            rr_host.push_back(rr);
        }
		}*/

		/*if(Para2==1)
		{
			std::mt19937 rng(12345); 
			std::shuffle(rr_host.begin(), rr_host.end(), rng);
		}
		else if(Para2==2)
		{
			std::vector<double> rr_hostr;
			std::uniform_int_distribution<size_t> dist(0, record - 1); 
			std::mt19937 rng(12345); 
			for (size_t ii = 0; ii < record; ii++) 
			{ 
				rr = rr_host[dist(rng)]; 
			    rr_hostr.push_back(rr); 
			}
			for (size_t ii = 0; ii < record; ii++) 
			{ 
				rr_host[ii]=rr_hostr[ii]; 
			}
			//rr_host=rr_hostr;
		}
		else if(Para2==3)
		{
			std::mt19937 rng(12345); 
			std::normal_distribution<double> gauss_dist(PCL1, RR_SD); 
			//std::vector<double> rr_gauss(10000); 
			//std::vector<double> rr_host(10000);
			for (size_t ii = 0; ii < rr_host.size(); ii++) 
			{ 
				rr = gauss_dist(rng); 
			    rr_host[ii]=rr; 
			}
		}*/

		double T_time,Fsti,Ts;
		int n,stim_gate,Num_j,j1,k,stim_site,clamp;
		double Fctrl(double T_time, int stim_gate, double Ts, double PCL);		

		//Host and Device Variables Settings 	
		size_t ArraySize = NDIM*Blocknum*Length*sizeof(double);	
		double *h_state_var,*h_state_deri,*h_vnew;
		double *d_state_var,*d_state_deri,*d_vnew;
		
		h_state_var=(double *)malloc(ArraySize);		
		h_state_deri=(double *)malloc(ArraySize);
		h_vnew=(double *)malloc(Blocknum*Length*sizeof(double));  
		cudaMalloc( (void**)&d_state_var,ArraySize );  
		cudaMalloc( (void**)&d_state_deri,ArraySize );
		cudaMalloc( (void**)&d_vnew,Blocknum*Length*sizeof(double) );
		// Get the initial state.
		get_initial_state(h_state_var,h_state_deri); 
		for(j1=0;j1<Blocknum;j1++)  
		{ 
				for(Num_j=0;Num_j<Length;Num_j++) 
				{
						h_vnew[Num_j+Length*j1]=h_state_var[pos(38,Num_j,j1)]; 	
				}
		}
		cudaMemcpy(d_state_var,h_state_var,ArraySize,cudaMemcpyHostToDevice);
		cudaMemcpy(d_state_deri,h_state_deri,ArraySize,cudaMemcpyHostToDevice);
		cudaMemcpy(d_vnew,h_vnew,Blocknum*Length*sizeof(double),cudaMemcpyHostToDevice);
		
		//APD and EAD Variables
		size_t ArraySize_doub = Blocknum*Length*sizeof(double);
		size_t ArraySize_int = Blocknum*Length*sizeof(int); 
		double *EAD_amp_max,*v_takeoff_max;
		double *APD,*EAD_amp,*EAD_min,*EAD_max,*v_takeoff,*EAD_period;
		int *APD_gate,*EAD_gate,*EAD_times,*EAD_count;
		double *dev_APD,*dev_EAD_amp,*dev_EAD_min,*dev_EAD_max,*dev_v_takeoff,*dev_EAD_period;
		int *dev_APD_gate,*dev_EAD_gate,*dev_EAD_times,*dev_EAD_count;

		double *ECG,*dev_ECG,*ECGo;
		double PCL,QT_m;
		int aa;
		int *ECG_up_reco;
		double current_stimuli,next_stimuli,RR_record;

        //aa=0;

		ECG=(double *)malloc(Blocknum*sizeof(double));
		ECGo=(double *)malloc(Blocknum*sizeof(double));
		//QT=(double *)malloc(ETM*sizeof(double));
		ECG_up_reco=(int *)malloc(Blocknum*sizeof(int));
		cudaMalloc( (void**)&dev_ECG,Blocknum*sizeof(double) );

		EAD_amp_max=(double *)malloc(ArraySize_doub);
		v_takeoff_max=(double *)malloc(ArraySize_doub); //MAX
		APD=(double *)malloc(ArraySize_doub);
		APD_gate=(int *)malloc(ArraySize_int);
		EAD_gate=(int *)malloc(ArraySize_int);
		EAD_times=(int *)malloc(ArraySize_int);
		EAD_count=(int *)malloc(ArraySize_int);	
		EAD_min=(double *)malloc(ArraySize_doub);
		EAD_max=(double *)malloc(ArraySize_doub);
		EAD_amp=(double *)malloc(ETM*ArraySize_doub);
		EAD_period=(double *)malloc(ETM*ArraySize_doub);
		v_takeoff=(double *)malloc(ETM*ArraySize_doub);				
		cudaMalloc((void**)&dev_APD,ArraySize_doub);
		cudaMalloc((void**)&dev_APD_gate,ArraySize_int);
		cudaMalloc((void**)&dev_EAD_gate,ArraySize_int);
		cudaMalloc((void**)&dev_EAD_times,ArraySize_int);
		cudaMalloc((void**)&dev_EAD_count,ArraySize_int);
		cudaMalloc((void**)&dev_EAD_max,ArraySize_doub);
		cudaMalloc((void**)&dev_EAD_min,ArraySize_doub);
		cudaMalloc((void**)&dev_EAD_amp,ETM*ArraySize_doub);
		cudaMalloc((void**)&dev_EAD_period,ETM*ArraySize_doub);
		cudaMalloc((void**)&dev_v_takeoff,ETM*ArraySize_doub);	
		for(Num_j=0;Num_j<Blocknum;Num_j++)
		{
			ECG[Num_j]=0;ECGo[Num_j]=0;ECG_up_reco[Num_j]=0;
				for(k=0;k<Length;k++)
				{		
						EAD_amp_max[k+Length*Num_j]=0.0;v_takeoff_max[k+Length*Num_j]=-100.0;
						APD_gate[k+Length*Num_j]=0;EAD_gate[k+Length*Num_j]=0;
						EAD_times[k+Length*Num_j]=0;EAD_count[k+Length*Num_j]=0;APD[k+Length*Num_j]=0.0;
						for(j1=0;j1<ETM;j1++)
						{
								EAD_amp[j1+ETM*(k+Length*Num_j)]=-10.0;v_takeoff[j1+ETM*(k+Length*Num_j)]=-100.0;
								EAD_period[j1+ETM*(k+Length*Num_j)]=0.0;
						}
				}
		}
		cudaMemcpy( dev_APD,APD,ArraySize_doub,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_APD_gate,APD_gate,ArraySize_int,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_gate,EAD_gate,ArraySize_int,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_times,EAD_times,ArraySize_int,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_count,EAD_count,ArraySize_int,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_max,EAD_max,ArraySize_doub,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_min,EAD_min,ArraySize_doub,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_amp,EAD_amp,ETM*ArraySize_doub,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_EAD_period,EAD_period,ETM*ArraySize_doub,cudaMemcpyHostToDevice );
		cudaMemcpy( dev_v_takeoff,v_takeoff,ETM*ArraySize_doub,cudaMemcpyHostToDevice );

		cudaMemcpy( dev_ECG,ECG,ArraySize_doub,cudaMemcpyHostToDevice );
		//APDs and EADs
								
		int *EAD_times_pause,*dev_EAD_times_pause;
		EAD_times_pause=(int *)malloc(Blocknum*sizeof(int));
		cudaMalloc((void**)&dev_EAD_times_pause,Blocknum*sizeof(int));
		for(j1=0;j1<Blocknum;j1++)
		{ 
				EAD_times_pause[j1]=0;
		}
		cudaMemcpy( dev_EAD_times_pause,EAD_times_pause,Blocknum*sizeof(int),cudaMemcpyHostToDevice );
		  
		stim_gate=1;Ts=100.0;stim_site=1;clamp=0;RR_record=0;next_stimuli=0;aa=0;
		for(n=0;RR_record<=rr_total;n++) 
		{
				T_time=n*dt;
				if(stim_gate==0&&T_time>10.0) 
				{
						stim_gate=1;Ts=T_time;
				} 

				/*if(RR_record>=burn)
				{PCL=rr_host[aa];}
				else
				{PCL=PCL1;}*/

				PCL=rr_host[aa];

				if(T_time>=next_stimuli&&T_time-next_stimuli<=STIM_OCT)
				{
					Fsti=M_STIM;
				}
				else if(T_time>next_stimuli+STIM_OCT)
				{
					Fsti=0;
					RR_record=RR_record+1;
					if(RR_record>=burn)
					{
						aa++;
					}
					current_stimuli=next_stimuli;
				    //rand_stimuli = normal(gen);
					//rand_stimuli=0;
					next_stimuli=PCL+current_stimuli;
				}
				else
				{
					Fsti=0;
				}

				//Fsti=Fctrl(T_time,stim_gate,Ts,PCL);
				update<<<Blocknum,Length>>>(d_state_var,d_state_deri,d_vnew,Fsti,stim_site,clamp,T_time,
											dev_APD,dev_APD_gate,
											dev_EAD_gate,dev_EAD_times,dev_EAD_amp,dev_EAD_min,dev_EAD_max,
											dev_EAD_count,dev_v_takeoff,dev_EAD_period,dev_EAD_times_pause,dev_ECG,Para1,Para3);
				if(n>=TEMPO&&fmod(n,150.0)<1e-6) 
				{
						cudaMemcpy(h_state_var,d_state_var,ArraySize,cudaMemcpyDeviceToHost);
						cudaMemcpy(ECG,dev_ECG,Blocknum*sizeof(double),cudaMemcpyDeviceToHost);

						if(RR_record>burn-100&&RR_record<burn+100)
						{

						fprintf(fileopen1,"%20.12f",T_time);
						double c_1=0.0e0;
						for(int cell_j1=2;cell_j1<Length-1;cell_j1=cell_j1+Lfactor)
						{
								fprintf(fileopen1,"%20.12f",h_state_var[pos(1,cell_j1,0)]-Vshift*c_1);c_1=c_1+1.0e0;
						}
						fprintf(fileopen1,"\n");
						//FILE 1 of Time Series
						}

						fprintf(fileopen2,"%20.12f %20.12f %20.12f %20.12f %20.12f\n",0.001*T_time,ECG[0],current_stimuli,RR_record,PCL);

						for(int trail=0;trail<Blocknum-1;trail++)
						{
							//if(fabs(ECGo[trail])<ECG_th&&fabs(ECG[trail]>ECG_th)&&ECG_up_reco[trail]<2)
							//{
							//	ECG_up_reco[trail]=ECG_up_reco[trail]+1;
							//}
							//dECG=fabs(ECG[trail]-ECGo[trail]);
							if(fabs(ECG[trail])<1e-3&&T_time-current_stimuli>250&&ECG_up_reco[trail]-RR_record<0)
							{
								QT_m=T_time-current_stimuli;
								ECG_up_reco[trail]=RR_record;
								fprintf(fileopen3,"%20.12f",T_time);
			                    fprintf(fileopen3,"%20.12f",QT_m);
								fprintf(fileopen3,"%20.12f",h_state_var[pos(2,0,0)]);
								fprintf(fileopen3,"%20.12f",h_state_var[pos(4,0,0)]);
								fprintf(fileopen3,"%20.12f",h_state_var[pos(6,0,0)]);
			                    fprintf(fileopen3,"\n");
							}
							ECGo[trail]=ECG[trail];
							//aa=aa+1;
							//if(bb==1&&aa<ETM)
							//{
							//	QT[aa]=QT_m;
							//	aa=aa+1;
							//	bb=0;
							//}
						}
				}
		}
				
		/*cudaMemcpy(EAD_times_pause,dev_EAD_times_pause,Blocknum*sizeof(int),cudaMemcpyDeviceToHost); 
		for(Num_j=0;Num_j<Blocknum;Num_j++)
		{
				if(EAD_times_pause[Num_j]>1)
						fprintf(fileopen2,"%20.12f%20.12f%5d\n",x_D,y_D,2);
			    else 
						fprintf(fileopen2,"%20.12f%20.12f%5d\n",x_D,y_D,1);  
				//---- PVC map --------------------------------------------------------------------------------------------------------------------------------------------------------
		}*/
		  
		free(h_state_var);free(h_vnew);
		cudaFree(d_state_var);cudaFree(d_vnew);
						
		free(EAD_amp);free(EAD_min);free(EAD_max);free(v_takeoff);free(EAD_amp_max);
		free(EAD_period);free(v_takeoff_max);
		free(EAD_gate);free(EAD_times);free(EAD_count);
			
		cudaFree(dev_EAD_amp);cudaFree(dev_EAD_min);cudaFree(dev_EAD_max);
		cudaFree(dev_v_takeoff);cudaFree(dev_EAD_period);
		cudaFree(dev_EAD_gate);cudaFree(dev_EAD_times);cudaFree(dev_EAD_count);cudaFree(dev_ECG);
		return(0);
} 

double Fctrl(double T_time, int stim_gate, double Ts, double PCL)
{
		double FF=0.0;
	
		if(stim_gate==1&&fmod(T_time,PCL)<=STIM_OCT) //stim_gate==1&&(T_time-Ts)-floor((T_time-Ts)/BCL)*BCL<=STIM_OCT
				FF=M_STIM;
		else if(stim_gate==2&&T_time-(Ts+Interval)<=STIM_OCT&&T_time-(Ts+Interval)>=0.0) //T_time>=Ts&&T_time-Ts<=STIM_OCT
				FF=M_STIM;
		else
				FF=0.0;
		return(FF);
}
