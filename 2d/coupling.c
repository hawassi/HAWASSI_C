#include "coupling.h"

void HAWASSI_getElev_from_CFD(char* dir_force,double time,double dt,double** Elev,int Nx,int Ny,
int nt_partition)
{	
	double* xdata;
	double* ydata;
	double N,M;
	int    Nxdata,Nydata;
	
	//find this time at which time partition
	int ipart = 1+(time/dt)/nt_partition;
	int iter  = (time/dt)-(ipart-1)*nt_partition;
	
	//open file with constants data from the folder dir_force, otherwise sleep
	char* name_constant = new char[256];
	sprintf(name_constant,"%sConstant_%02d.dat",dir_force,ipart);	
	FILE* fc = fopen(name_constant,"rb");
	while (!fc){
		usleep(1);
		fc = fopen(name_constant,"rb");
	}
	fread(&N,sizeof(double),1,fc);Nxdata=(int) N;
	xdata = new double[Nxdata];
	fread(xdata,sizeof(double),Nxdata,fc);
	fread(&M,sizeof(double),1,fc);Nydata=(int) M;
	ydata = new double[Nydata];
	fread(ydata,sizeof(double),Nydata,fc);
	fclose(fc);	
	
	//open file with elevation data from the folder dir_force, otherwise sleep
	char* name_file = new char[256];
	sprintf(name_file,"%sCFD_eta_%02d.dat",dir_force,ipart);	
	FILE* feta = fopen(name_file,"rb");
	while (!feta){
		usleep(1);
		feta = fopen(name_file,"rb");
	}
	
	//reading data when it's ready, otherwise sleep
	int size;
	fseek(feta, 0, SEEK_END);
	size = ftell(feta);
	while (size<(iter+1)*Nxdata*Nydata*sizeof(double)){
		fflush(feta);
		fseek(feta, 0, SEEK_END);
		size = ftell(feta);		
	}
	fseek(feta, 0, SEEK_SET);
	fseek(feta,iter*Nxdata*Nydata*sizeof(double),SEEK_SET);	
	
	//temporary pointer
	double* Elev_temp = new double[Nxdata*Nydata];	
	for (int i=0;i<Nydata;i++){
	  for(int j=0;j<Nxdata;j++){
		  fread(&Elev_temp[i*Nxdata+j],sizeof(double),1,feta);
	  }
    }

    fclose(feta);
	delete[] name_file;  
	
	//interpolation
	double xmin = fun_min(xdata,Nxdata);
	double xmax = fun_max(xdata,Nxdata);
	double ymin = fun_min(ydata,Nxdata);
	double ymax = fun_max(ydata,Nxdata);
	
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, Nxdata, Nydata);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();		
	gsl_spline2d_init(spline, xdata, ydata, Elev_temp, Nxdata, Nydata);
	
	//filling to HW grid
	set_matrix_val(Elev,Nx,Ny,0.0);
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  if ((dominfo.x[j]>=xmin)&&(dominfo.x[j]<=xmax)&&(dominfo.y[i]>=ymin)&&(dominfo.y[i]<=ymax)){
			  Elev[i][j] = gsl_spline2d_eval(spline, dominfo.x[j],dominfo.y[i],xacc, yacc);
		  }
	  }
    }
    
    gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc); 
	delete[] Elev_temp;
	delete[] xdata;
	delete[] ydata;
}

void HAWASSI_getPot_from_CFD(char* dir_force,double time,double dt,double** Pot,int Nx,int Ny,
int nt_partition)
{	
	double* xdata;
	double* ydata;
	double N,M;
	int    Nxdata,Nydata;
	
	//find this time at which time partition
	int ipart = 1+(time/dt)/nt_partition;
	int iter  = (time/dt)-(ipart-1)*nt_partition;
	
	//open file with constants data from the folder dir_force, otherwise sleep
	char* name_constant = new char[256];
	sprintf(name_constant,"%sConstant_%02d.dat",dir_force,ipart);	
	FILE* fc = fopen(name_constant,"rb");
	while (!fc){
		usleep(1);
		fc = fopen(name_constant,"rb");
	}
	fread(&N,sizeof(double),1,fc);Nxdata=(int) N;
	xdata = new double[Nxdata];
	fread(xdata,sizeof(double),Nxdata,fc);
	fread(&M,sizeof(double),1,fc);Nydata=(int) M;
	ydata = new double[Nydata];
	fread(ydata,sizeof(double),Nydata,fc);
	fclose(fc);
	
	//open file with potential data from the folder dir_force, otherwise sleep
	char* name_file = new char[256];
	sprintf(name_file,"%sCFD_phi_%02d.dat",dir_force,ipart);	
	FILE* fphi = fopen(name_file,"rb");
	while (!fphi){
		usleep(1);
		fphi = fopen(name_file,"rb");
	}
	
	//reading data when it's ready, otherwise sleep
	int size;
	fseek(fphi, 0, SEEK_END);
	size = ftell(fphi);
	while (size<(iter+1)*Nxdata*Nydata*sizeof(double)){
		fflush(fphi);
		fseek(fphi, 0, SEEK_END);
		size = ftell(fphi);		
	}
	fseek(fphi, 0, SEEK_SET);
	fseek(fphi,iter*Nxdata*Nydata*sizeof(double),SEEK_SET);	
	
	//temporary pointer
	double* Pot_temp = new double[Nxdata*Nydata];	
	for (int i=0;i<Nydata;i++){
	  for(int j=0;j<Nxdata;j++){
		  fread(&Pot_temp[i*Nxdata+j],sizeof(double),1,fphi);
	  }
    }

    fclose(fphi);
	delete[] name_file;  
	
	//interpolation
	double xmin = fun_min(xdata,Nxdata);
	double xmax = fun_max(xdata,Nxdata);
	double ymin = fun_min(ydata,Nxdata);
	double ymax = fun_max(ydata,Nxdata);
	
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, Nxdata, Nydata);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();		
	gsl_spline2d_init(spline, xdata, ydata, Pot_temp, Nxdata, Nydata);
	
	//filling to HW grid
	set_matrix_val(Pot,Nx,Ny,0.0);
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  if ((dominfo.x[j]>=xmin)&&(dominfo.x[j]<=xmax)&&(dominfo.y[i]>=ymin)&&(dominfo.y[i]<=ymax)){
			  Pot[i][j] = gsl_spline2d_eval(spline, dominfo.x[j],dominfo.y[i],xacc, yacc);
		  }
	  }
    }
    
    gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc); 
	delete[] Pot_temp;
	delete[] xdata;
	delete[] ydata;
}

void HAWASSI_coupling_term(char* dir_force,double time,double dt,double** etaf,double** phif,
int nt_partition,int Nx,int Ny)
{
	double**  temp_eta 	= declare_2darray(Nx,Ny);
	double**  temp_phi 	= declare_2darray(Nx,Ny);
	double**  couple_eta= declare_2darray(Nx,Ny);
	double**  couple_phi= declare_2darray(Nx,Ny);
	
	printf("%.2f \n",time);
	
	HAWASSI_getElev_from_CFD(dir_force,time,dt,couple_eta,Nx,Ny,nt_partition);
	HAWASSI_getPot_from_CFD(dir_force,time,dt,couple_phi,Nx,Ny,nt_partition);
	
	//coupling term in real space
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		temp_eta[i][j] 	= alpha_couple*Chi_couple_H[i][j]*(couple_eta[i][j]-etaf[i][j]);
		temp_phi[i][j] 	= alpha_couple*Chi_couple_H[i][j]*(couple_phi[i][j]-phif[i][j]);
	  }
	}
	complex** temp_etahat = fftw_2d_r2c(temp_eta,Nx,Ny);
	complex** temp_phihat = fftw_2d_r2c(temp_phi,Nx,Ny);

	//filling
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		couple_etahat[i][j].Re = temp_etahat[i][j].Re;
		couple_phihat[i][j].Re = temp_phihat[i][j].Re;
		couple_etahat[i][j].Im = temp_etahat[i][j].Im;
		couple_phihat[i][j].Im = temp_phihat[i][j].Im;
	  }
	}
	
	free_2darray(temp_eta,Nx,Ny);
	free_2darray(temp_phi,Nx,Ny);
	free_2darray(couple_eta,Nx,Ny);
	free_2darray(couple_phi,Nx,Ny);
	free_2darray_complex(temp_etahat,Nx,Ny);
	free_2darray_complex(temp_phihat,Nx,Ny);	
}
