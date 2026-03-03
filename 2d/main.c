#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/OrderingMethods>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include <fftw3.h>
#include <omp.h>
#include <ctime>

#define On  (1)
#define Off (0)
#define Save_Fig (On)
#define Save_Data (Off)

#include "var_initialization.h"

//for coupling
double xinterv1,xinterv2,yinterv1,yinterv2;
double lcfd,lcfd_y;
double lenOZin,lenOZiny,lenOZtr,lenOZtry;
double ozstart[2];
double ozend[2];
complex** couple_etahat;
complex** couple_phihat;
complex** couple_uhat;
complex** couple_vhat;
double 	  alpha_couple;
double**  Chi_couple_H;
double**  Chi_couple_Hy;
int 	  iter_t=0;
int id_Xinterv1,id_Xinterv2;
int id_Yinterv1,id_Yinterv2;
int	nLin,nLiny,nLtr,nLtry;  
int nCFD,nCFDy;

//for kinematic
char    filename[256];
double* GQ;
double  xscale,tscale,potscale,dx_rec,dy_rec;
int 	Nxrec,Nyrec,Nzrec=20;
double  *xrec,*yrec;
double** depthrec;
int ip_Xinterv1,ip_Xinterv2;
int ip_Yinterv1,ip_Yinterv2;
int NtoutIP,NxoutIP,NyoutIP,NzoutIP,Ntrest,NxyoutIP,ntpart;
double  **toutIP,*trest,*xoutIP,*youtIP,*zoutIP,maxEta,minEta,minD,maxD;
FILE 	*dateta,*datvel,*datP,*datacc,*datPhi,*dat_record;
char var_elev[20],var_vel[20],var_P[20],var_acc[20],var_phi[20];
double *dy_hat,*y_hat;
double **Velx,**Vely,**Velz;
double **phi_prev,**depth_IP;
double **Af;
typedef Eigen::SparseMatrix<double> SpMat;
SpMat A(1,1);
typedef Eigen::LeastSquaresConjugateGradient<SpMat> solver;
//typedef Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int>> solver;
solver lscg;

#include "functions1.c"
#include "functions.c"
#include "breaking.c"
#include "UDW2p_HAWASSI2D.c"
#include "kinematic_modul.c"
#include "fun_odesolver.c"
			
int main(int argc, char* argv[])
{	
	//Checking argument
	if (argc != 4) {
        fprintf(stderr,"Usage (4 argv): Hawassi2D working_dir outname (integer)\n");
        return 1;
    }
	arg1 = argv[1];
	sprintf(arg2, "%s%s/", arg1,argv[2]);
	ncpu = atoi(argv[3]);	
	
	clock_t launch = clock();
	
	printf("======================================================\n"); 
	printf("HAWASSI-2D Software\n");
	printf("June 2021\n");
	printf("LabMath-Indonesia\n");
	printf("e-mail: info@labmath-indonesia.org\n");
	printf("======================================================\n"); 
	
	//read input file-----------------------------------------------
	char params[256];
	sprintf(params, "%s/input.txt", arg1);
	FILE* fpinput = fopen(params, "r");	
    if(fpinput==NULL){printf("Can not open input.txt!!\n");exit(0);}    
    read_input(fpinput);
    fclose(fpinput);
    
    //setup FFTW in multi-threads dan make general plan
    int ncpu_fft = 8;
    if (ncpu>8) {ncpu_fft = ncpu;}
    omp_set_num_threads(ncpu_fft);
    fftw_init_threads();
    fftw_plan_with_nthreads(ncpu_fft);

    //spatial setup-------------------------------------------------
    fun_param_init();
    tprev = dominfo.t[0];
    printf("Reading parameter done\n");

    //initial or influx setup----------------------------------------------
    #pragma omp parallel sections 
    {
		#pragma omp section
		waveinit=fun_generation_init(dominfo.Nx,dominfo.Ny);

		#pragma omp section
		if (strcmp(wavename,"zero")!=0){
			fun_generation_influx(dominfo.Nt);
		}
		else {
			meandepthinf = -dominfo.bathy.data[0][0]/2;
		}
	}
    
    //sea characteristics
	fun_calcseacharacteristic(wave_influx);
	
	//add ramp to original signal
	fun_ramp2d(ramp,taper);
	
	//area generation 	
	if (strcmp(wavename,"zero")!=0){
		fun_area_generation();	
	}
	
	//operator setup
	Oprt = fun_operator_setup(dominfo);	
	printf("Wave period = %.2f \n",seainfo.Tp);
	printf("Wave length = %.2f \n",seainfo.lambda);
	printf("Cpeak = %.2f \n",Oprt.Cpeak);
	
	//initiate coupling module
	couple_etahat = declare_2darray_complex(dominfo.Nx,dominfo.Ny);
	couple_uhat = declare_2darray_complex(dominfo.Nx,dominfo.Ny);
	couple_vhat = declare_2darray_complex(dominfo.Nx,dominfo.Ny);
	set_matrix_complex_val(couple_etahat,dominfo.Nx,dominfo.Ny,0.0);
	set_matrix_complex_val(couple_uhat,dominfo.Nx,dominfo.Ny,0.0);
	set_matrix_complex_val(couple_vhat,dominfo.Nx,dominfo.Ny,0.0);
	
	if (strcmp(coupling,"yes")==0){		
		alpha_couple  = grav*(2*M_PI/seainfo.kp)/(Oprt.Cpeak*lenOZin);
		printf("Decay rate = %.2f \n",alpha_couple);
		
		Chi_couple_H  = declare_2darray(dominfo.Nx,dominfo.Ny);
		set_matrix_val(Chi_couple_H,dominfo.Nx,dominfo.Ny,0.0);				
		define_coupling_char(dominfo.x,dominfo.y,dominfo.Nx,dominfo.Ny);		

		char* fout_name = new char[256];
		sprintf(fout_name,"%s/decay.dat",arg2);
		FILE* falpha=fopen(fout_name,"w");
		fprintf(falpha,"%lf",alpha_couple);
		delete[] fout_name;
		fclose(falpha);
	}
	
	//setup if breaking
	Sb_xhat = declare_2darray_complex(dominfo.Nx,dominfo.Ny);
	Sb_yhat = declare_2darray_complex(dominfo.Nx,dominfo.Ny);
	set_matrix_complex_val(Sb_xhat,dominfo.Nx,dominfo.Ny,0.0);
	set_matrix_complex_val(Sb_yhat,dominfo.Nx,dominfo.Ny,0.0);
	if(strcmp(breaking,"yes")==0){		
		parambreak(dominfo.Nx,dominfo.Ny);
	}
	
	//setup if kinematic yes
	if(strcmp(kinematic,"yes")==0){	
		phi_prev=declare_2darray(dominfo.Nx,dominfo.Ny);	
		set_matrix_val(phi_prev,dominfo.Nx,dominfo.Ny,0.0);
		read_kinematic(dominfo.bathy.pos,dominfo.Nx,dominfo.Ny,dominfo.Nt);
	}
	
	//ODE
	double eps_rel=1e-2;
	double eps_abs=eps_rel/10.0;
	size_t Ndat=6*dominfo.Nx*dominfo.Ny;
	y_hat  = new double[Ndat];
	dy_hat = new double[Ndat];
	set_init_ode(waveinit,dominfo.Nx,dominfo.Ny);
	set_init_double(dy_hat,Ndat);
	
	gsl_odeiv2_system sys;
	printf("Calculating time evolution...\n"); 
	gettimeofday(&start, NULL);
	
	if (runup_id==0){//no runup
		if (strcmp(evolinfo.model,"HS1")==0){
			if(strcmp(dominfo.bathy.Id,"flat")==0){
				sys = {rhs_AB1_2DF_uv, NULL, Ndat, NULL};				 
			}
			else{
				sys = {rhs_AB1_2DB_uv, NULL, Ndat, NULL};	
			}		 
		}
		else if (strcmp(evolinfo.model,"HS2")==0){
			if(strcmp(dominfo.bathy.Id,"flat")==0){
				sys = {rhs_AB2_2DF_uv, NULL, Ndat, NULL};				 
			}
			else{
				sys = {rhs_AB2_2DB_uv, NULL, Ndat, NULL};	
			}	
		}
		else if (strcmp(evolinfo.model,"HS3")==0){			
			if(strcmp(dominfo.bathy.Id,"flat")==0){
				sys = {rhs_AB3_2DF_uv, NULL, Ndat, NULL};				 
			}
			else{
				sys = {rhs_AB3_2DB_uv, NULL, Ndat, NULL};	
			}	
		}		
	}
	else{
		sys = {rhs_shore2D_uv, NULL, Ndat, NULL};	
	}	
	
	ode_solver_uv(sys,y_hat,dominfo.t,Ndat,dominfo.Nx,dominfo.Ny,dominfo.Nt,eps_abs,eps_rel); 
	delete[] y_hat;

	fftw_cleanup_threads();
	delete[] xgauge;
	delete[] ygauge;
	delete[] index_gauge[0];
	delete[] index_gauge[1];
	free_2darray(cfSA,dominfo.Nx,dominfo.Ny);
	
	gettimeofday(&end, NULL);	
	int n = dominfo.Nt;
	double rel = (((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6)/(dominfo.t[n-1]-dominfo.t[0]);
	printf("100%%\nCompleted.\n");
	printf("Relative comp. time=%g[percent]\n",rel);
	  
    return 0;
}
