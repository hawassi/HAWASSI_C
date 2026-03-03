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
	double ymin = fun_min(ydata,Nydata);
	double ymax = fun_max(ydata,Nydata);
	
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

void fun_div(fftw_complex* resultdiv,fftw_complex* part_x,fftw_complex* part_y,double koef,int Nx,int Ny)
{
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			resultdiv[i*Nx+j][0] = koef*((-dominfo.kx[j]*part_x[i*Nx+j][1])+(-dominfo.ky[i]*part_y[i*Nx+j][1]));
			resultdiv[i*Nx+j][1] = koef*((dominfo.kx[j]*part_x[i*Nx+j][0])+(dominfo.ky[i]*part_y[i*Nx+j][0]));
		}
	}	
}

/*void HAWASSI_getPot_from_CFD(char* dir_force,double time,double dt,double** Pot,int Nx,int Ny,
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
	double ymin = fun_min(ydata,Nydata);
	double ymax = fun_max(ydata,Nydata);
	
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
}*/

void HAWASSI_getu_from_CFD(char* dir_force,double time,double dt,double** u,int Nx,int Ny,
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
	
	//open file with u data from the folder dir_force, otherwise sleep
	char* name_file = new char[256];
	sprintf(name_file,"%sCFD_u_%02d.dat",dir_force,ipart);	
	FILE* fu = fopen(name_file,"rb");
	while (!fu){
		usleep(1);
		fu = fopen(name_file,"rb");
	}
	
	//reading data when it's ready, otherwise sleep
	int size;
	fseek(fu, 0, SEEK_END);
	size = ftell(fu);
	while (size<(iter+1)*Nxdata*Nydata*sizeof(double)){
		fflush(fu);
		fseek(fu, 0, SEEK_END);
		size = ftell(fu);		
	}
	fseek(fu, 0, SEEK_SET);
	fseek(fu,iter*Nxdata*Nydata*sizeof(double),SEEK_SET);	
	
	//temporary pointer
	double* u_temp = new double[Nxdata*Nydata];	
	for (int i=0;i<Nydata;i++){
	  for(int j=0;j<Nxdata;j++){
		  fread(&u_temp[i*Nxdata+j],sizeof(double),1,fu);
	  }
    }

    fclose(fu);
	delete[] name_file;  
	
	//interpolation
	double xmin = fun_min(xdata,Nxdata);
	double xmax = fun_max(xdata,Nxdata);
	double ymin = fun_min(ydata,Nydata);
	double ymax = fun_max(ydata,Nydata);
	
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, Nxdata, Nydata);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();		
	gsl_spline2d_init(spline, xdata, ydata, u_temp, Nxdata, Nydata);
	
	//filling to HW grid
	set_matrix_val(u,Nx,Ny,0.0);
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  if ((dominfo.x[j]>=xmin)&&(dominfo.x[j]<=xmax)&&(dominfo.y[i]>=ymin)&&(dominfo.y[i]<=ymax)){
			  u[i][j] = gsl_spline2d_eval(spline, dominfo.x[j],dominfo.y[i],xacc, yacc);
		  }
	  }
    }
    
    gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc); 
	delete[] u_temp;
	delete[] xdata;
	delete[] ydata;
}

void HAWASSI_getv_from_CFD(char* dir_force,double time,double dt,double** v,int Nx,int Ny,
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
	
	//open file with u data from the folder dir_force, otherwise sleep
	char* name_file = new char[256];
	sprintf(name_file,"%sCFD_v_%02d.dat",dir_force,ipart);	
	FILE* fu = fopen(name_file,"rb");
	while (!fu){
		usleep(1);
		fu = fopen(name_file,"rb");
	}
	
	//reading data when it's ready, otherwise sleep
	int size;
	fseek(fu, 0, SEEK_END);
	size = ftell(fu);
	while (size<(iter+1)*Nxdata*Nydata*sizeof(double)){
		fflush(fu);
		fseek(fu, 0, SEEK_END);
		size = ftell(fu);		
	}
	fseek(fu, 0, SEEK_SET);
	fseek(fu,iter*Nxdata*Nydata*sizeof(double),SEEK_SET);	
	
	//temporary pointer
	double* v_temp = new double[Nxdata*Nydata];	
	for (int i=0;i<Nydata;i++){
	  for(int j=0;j<Nxdata;j++){
		  fread(&v_temp[i*Nxdata+j],sizeof(double),1,fu);
	  }
    }

    fclose(fu);
	delete[] name_file;  
	
	//interpolation
	double xmin = fun_min(xdata,Nxdata);
	double xmax = fun_max(xdata,Nxdata);
	double ymin = fun_min(ydata,Nydata);
	double ymax = fun_max(ydata,Nydata);
	
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, Nxdata, Nydata);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();		
	gsl_spline2d_init(spline, xdata, ydata, v_temp, Nxdata, Nydata);
	
	//filling to HW grid
	set_matrix_val(v,Nx,Ny,0.0);
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  if ((dominfo.x[j]>=xmin)&&(dominfo.x[j]<=xmax)&&(dominfo.y[i]>=ymin)&&(dominfo.y[i]<=ymax)){
			  v[i][j] = gsl_spline2d_eval(spline, dominfo.x[j],dominfo.y[i],xacc, yacc);
		  }
	  }
    }
    
    gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc); 
	delete[] v_temp;
	delete[] xdata;
	delete[] ydata;
}


void HAWASSI_coupling_term_uv(char* dir_force,double time,double dt,double** etaf,double** uf,double** vf,
int nt_partition,int Nx,int Ny)
{
	double**  temp_eta 	= declare_2darray(Nx,Ny);
	double**  temp_u 	= declare_2darray(Nx,Ny);
	double**  temp_v 	= declare_2darray(Nx,Ny);
	double**  couple_eta= declare_2darray(Nx,Ny);
	double**  couple_u= declare_2darray(Nx,Ny);
	double**  couple_v= declare_2darray(Nx,Ny);		
	
	HAWASSI_getElev_from_CFD(dir_force,time,dt,couple_eta,Nx,Ny,nt_partition);
	HAWASSI_getu_from_CFD(dir_force,time,dt,couple_u,Nx,Ny,nt_partition);
	HAWASSI_getv_from_CFD(dir_force,time,dt,couple_v,Nx,Ny,nt_partition);
		
	//coupling term in real space
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		temp_eta[i][j] 	= alpha_couple*(Chi_couple_H[i][j])*(couple_eta[i][j]-etaf[i][j]);
		temp_u[i][j] 	= alpha_couple*(Chi_couple_H[i][j])*(couple_u[i][j]-uf[i][j]);
		temp_v[i][j] 	= alpha_couple*(Chi_couple_H[i][j])*(couple_v[i][j]-vf[i][j]);
	  }
	}
	complex** temp_etahat = declare_2darray_complex(Nx,Ny);
	complex** temp_uhat = declare_2darray_complex(Nx,Ny);
	complex** temp_vhat = declare_2darray_complex(Nx,Ny);
	fftw_2d_r2c(temp_etahat,temp_eta,Nx,Ny);
	fftw_2d_r2c(temp_uhat,temp_u,Nx,Ny);
	fftw_2d_r2c(temp_vhat,temp_v,Nx,Ny);
	
	//filling
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		couple_etahat[i][j].Re = temp_etahat[i][j].Re;
		couple_etahat[i][j].Im = temp_etahat[i][j].Im;
		couple_uhat[i][j].Re = temp_uhat[i][j].Re;
		couple_uhat[i][j].Im = temp_uhat[i][j].Im;
		couple_vhat[i][j].Re = temp_vhat[i][j].Re;
		couple_vhat[i][j].Im = temp_vhat[i][j].Im;
	  }
	}
	
	free_2darray(temp_eta,Nx,Ny);
	free_2darray(temp_u,Nx,Ny);
	free_2darray(temp_v,Nx,Ny);
	free_2darray(couple_eta,Nx,Ny);
	free_2darray(couple_u,Nx,Ny);
	free_2darray(couple_v,Nx,Ny);
	free_2darray_complex(temp_etahat,Nx,Ny);
	free_2darray_complex(temp_uhat,Nx,Ny);
	free_2darray_complex(temp_vhat,Nx,Ny);
	
}

void print_data_gauges_uv(int i,double t,double** eta,double** u,double** v,FILE* fgauge1,FILE* fgauge2,FILE* fgauge3)
{
	int ii,jj;
	if (i==0){
		//first row
		fprintf(fgauge1,"0 ");
		fprintf(fgauge2,"0 ");
		fprintf(fgauge3,"0 ");
		for (int j=0;j<ngauge;j++){
			jj = index_gauge[0][j];//xgauge
			fprintf(fgauge1,"%f ",dominfo.x[jj]);
			fprintf(fgauge2,"%f ",dominfo.x[jj]);
			fprintf(fgauge3,"%f ",dominfo.x[jj]);
		}
		fprintf(fgauge1,"\n");
		fprintf(fgauge2,"\n");
		fprintf(fgauge3,"\n");
		//second row
		fprintf(fgauge1,"0 ");
		fprintf(fgauge2,"0 ");
		fprintf(fgauge3,"0 ");
		for (int j=0;j<ngauge;j++){
			ii = index_gauge[1][j];//ygauge
			fprintf(fgauge1,"%f ",dominfo.y[ii]);
			fprintf(fgauge2,"%f ",dominfo.y[ii]);
			fprintf(fgauge3,"%f ",dominfo.y[ii]);
		}
		fprintf(fgauge1,"\n");
		fprintf(fgauge2,"\n");
		fprintf(fgauge3,"\n");
	}
	fprintf(fgauge1,"%f ",t);
	fprintf(fgauge2,"%f ",t);
	fprintf(fgauge3,"%f ",t);
	
	for (int j=0;j<ngauge;j++) {
		ii = index_gauge[1][j];//ygauge
		jj = index_gauge[0][j];//xgaugef
		fprintf(fgauge1,"%f ",eta[ii][jj]);
		fprintf(fgauge2,"%f ",u[ii][jj]);
		fprintf(fgauge3,"%f ",v[ii][jj]);
	}
	fprintf(fgauge1,"\n");
	fprintf(fgauge2,"\n");
	fprintf(fgauge3,"\n");	
	
	fflush(fgauge1);
	fflush(fgauge2);
	fflush(fgauge3);
}


void get_eta(fftw_complex* eta, double t, const double* Zhat, int Nx, int Ny)
{	
	//Define eta
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
	// Create plans 
	fftw_plan plan_iftetahat   = fftw_plan_dft_2d(Ny,Nx, etahat, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	
		}
	}		
	fftwcomplex_2d_sym(etahat,Nx,Ny);
	fftw_execute(plan_iftetahat);
	fftw_destroy_plan(plan_iftetahat);	
	fftw_free(etahat);
}

void HSSgen_2depth(fftw_complex* phi, fftw_complex* phihat, double** Oprt_min,double** Oprt_plus,fftw_complex* result,int Nx,int Ny,int factIFFT)
{
	fftw_complex *gammin_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gamplus_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammin_L0minphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gamplus_L0plusphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0minphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0plusphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0minphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0plusphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));		
	
	// Create plans 
	fftw_plan plan_iftL0minphihat   = fftw_plan_dft_2d(Ny,Nx, L0minphihat,  L0minphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftL0plusphihat  = fftw_plan_dft_2d(Ny,Nx, L0plusphihat, L0plusphi,FFTW_BACKWARD, FFTW_ESTIMATE);
    
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		L0minphihat[i*Nx+j][0]=Oprt_min[i][j]*phihat[i*Nx+j][0];
		L0minphihat[i*Nx+j][1]=Oprt_min[i][j]*phihat[i*Nx+j][1];
		L0plusphihat[i*Nx+j][0]=Oprt_plus[i][j]*phihat[i*Nx+j][0];
		L0plusphihat[i*Nx+j][1]=Oprt_plus[i][j]*phihat[i*Nx+j][1];
	  }
	}
	
	fftwcomplex_2d_sym(L0minphihat,Nx,Ny);
	fftwcomplex_2d_sym(L0plusphihat,Nx,Ny);
    fftw_execute(plan_iftL0minphihat);
    fftw_execute(plan_iftL0plusphihat);
    
	fftw_destroy_plan(plan_iftL0minphihat);
	fftw_destroy_plan(plan_iftL0plusphihat);   
	
	fftw_complex* gammin_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* gamplus_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammin_L0minphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* gamplus_L0plusphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
    // Create plans 
    fftw_plan plan_fftgammin_phi   	   = fftw_plan_dft_2d(Ny,Nx, gammin_phi , gammin_phihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgamplus_phi      = fftw_plan_dft_2d(Ny,Nx, gamplus_phi, gamplus_phihat, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammin_L0minphi  = fftw_plan_dft_2d(Ny,Nx, gammin_L0minphi , gammin_L0minphihat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftgamplus_L0plusphi= fftw_plan_dft_2d(Ny,Nx, gamplus_L0plusphi, gamplus_L0plusphihat, FFTW_FORWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		gammin_phi[i*Nx+j][0]=Oprt.InterpD.Gam_min[i][j]*phi[i*Nx+j][0]/factIFFT;
		gammin_phi[i*Nx+j][1]=0;
		gamplus_phi[i*Nx+j][0]=Oprt.InterpD.Gam_plus[i][j]*phi[i*Nx+j][0]/factIFFT;
		gamplus_phi[i*Nx+j][1]=0;
		gammin_L0minphi[i*Nx+j][0]=Oprt.InterpD.Gam_min[i][j]*L0minphi[i*Nx+j][0]/(Nx*Ny);
		gammin_L0minphi[i*Nx+j][1]=0;
		gamplus_L0plusphi[i*Nx+j][0]=Oprt.InterpD.Gam_plus[i][j]*L0plusphi[i*Nx+j][0]/(Nx*Ny);
		gamplus_L0plusphi[i*Nx+j][1]=0;		
	  } 
	}
    
	fftw_execute(plan_fftgammin_phi);
    fftw_execute(plan_fftgamplus_phi);
	fftw_execute(plan_fftgammin_L0minphi);
	fftw_execute(plan_fftgamplus_L0plusphi);	
	
	fftw_destroy_plan(plan_fftgammin_phi);
    fftw_destroy_plan(plan_fftgamplus_phi);
	fftw_destroy_plan(plan_fftgammin_L0minphi);
	fftw_destroy_plan(plan_fftgamplus_L0plusphi);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			result[i*Nx+j][0] =(gammin_L0minphihat[i*Nx+j][0]+Oprt_min[i][j]*gammin_phihat[i*Nx+j][0]
								+gamplus_L0plusphihat[i*Nx+j][0]+Oprt_plus[i][j]*gamplus_phihat[i*Nx+j][0])/2;
			result[i*Nx+j][1] =(gammin_L0minphihat[i*Nx+j][1]+Oprt_min[i][j]*gammin_phihat[i*Nx+j][1]
								+gamplus_L0plusphihat[i*Nx+j][1]+Oprt_plus[i][j]*gamplus_phihat[i*Nx+j][1])/2;
		}
    }

	fftw_free(L0minphihat);fftw_free(L0plusphihat);
    fftw_free(L0minphi);fftw_free(L0plusphi);
	fftw_free(gammin_phi);fftw_free(gamplus_phi);
	fftw_free(gammin_L0minphi);fftw_free(gamplus_L0plusphi);
	fftw_free(gammin_phihat);fftw_free(gamplus_phihat);
	fftw_free(gammin_L0minphihat);fftw_free(gamplus_L0plusphihat);
    
}


void HSSgen_3depth(fftw_complex* phi, fftw_complex* phihat, double** Oprt_min,double** Oprt_mid,double** Oprt_plus,
fftw_complex* result,int Nx,int Ny,int factIFFT)
{
	fftw_complex *gammin_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammid_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gamplus_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammin_L0minphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammid_L0midphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gamplus_L0plusphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0minphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0midphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0plusphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0minphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0midphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0plusphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));		
	
	// Create plans 
	fftw_plan plan_iftL0minphihat   = fftw_plan_dft_2d(Ny,Nx, L0minphihat,  L0minphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftL0midphihat   = fftw_plan_dft_2d(Ny,Nx, L0midphihat,  L0midphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftL0plusphihat  = fftw_plan_dft_2d(Ny,Nx, L0plusphihat, L0plusphi,FFTW_BACKWARD, FFTW_ESTIMATE);
    
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		L0minphihat[i*Nx+j][0]=Oprt_min[i][j]*phihat[i*Nx+j][0];
		L0minphihat[i*Nx+j][1]=Oprt_min[i][j]*phihat[i*Nx+j][1];
		L0midphihat[i*Nx+j][0]=Oprt_mid[i][j]*phihat[i*Nx+j][0];
		L0midphihat[i*Nx+j][1]=Oprt_mid[i][j]*phihat[i*Nx+j][1];
		L0plusphihat[i*Nx+j][0]=Oprt_plus[i][j]*phihat[i*Nx+j][0];
		L0plusphihat[i*Nx+j][1]=Oprt_plus[i][j]*phihat[i*Nx+j][1];
	  }
	}
	
	fftwcomplex_2d_sym(L0minphihat,Nx,Ny);
	fftwcomplex_2d_sym(L0midphihat,Nx,Ny);
	fftwcomplex_2d_sym(L0plusphihat,Nx,Ny);
    fftw_execute(plan_iftL0minphihat);
    fftw_execute(plan_iftL0midphihat);
    fftw_execute(plan_iftL0plusphihat);
    
	fftw_destroy_plan(plan_iftL0minphihat);
	fftw_destroy_plan(plan_iftL0midphihat);
	fftw_destroy_plan(plan_iftL0plusphihat);   
	
	fftw_complex* gammin_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammid_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* gamplus_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammin_L0minphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammid_L0midphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* gamplus_L0plusphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
    // Create plans 
    fftw_plan plan_fftgammin_phi   	   = fftw_plan_dft_2d(Ny,Nx, gammin_phi , gammin_phihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammid_phi   	   = fftw_plan_dft_2d(Ny,Nx, gammid_phi , gammid_phihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgamplus_phi      = fftw_plan_dft_2d(Ny,Nx, gamplus_phi, gamplus_phihat, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammin_L0minphi  = fftw_plan_dft_2d(Ny,Nx, gammin_L0minphi , gammin_L0minphihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammid_L0midphi  = fftw_plan_dft_2d(Ny,Nx, gammid_L0midphi , gammid_L0midphihat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftgamplus_L0plusphi= fftw_plan_dft_2d(Ny,Nx, gamplus_L0plusphi, gamplus_L0plusphihat, FFTW_FORWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		gammin_phi[i*Nx+j][0]=Oprt.InterpD.Gam_min[i][j]*phi[i*Nx+j][0]/factIFFT;
		gammin_phi[i*Nx+j][1]=0;
		gammid_phi[i*Nx+j][0]=Oprt.InterpD.Gam_mid[i][j]*phi[i*Nx+j][0]/factIFFT;
		gammid_phi[i*Nx+j][1]=0;
		gamplus_phi[i*Nx+j][0]=Oprt.InterpD.Gam_plus[i][j]*phi[i*Nx+j][0]/factIFFT;
		gamplus_phi[i*Nx+j][1]=0;
		gammin_L0minphi[i*Nx+j][0]=Oprt.InterpD.Gam_min[i][j]*L0minphi[i*Nx+j][0]/(Nx*Ny);
		gammin_L0minphi[i*Nx+j][1]=0;
		gammid_L0midphi[i*Nx+j][0]=Oprt.InterpD.Gam_mid[i][j]*L0midphi[i*Nx+j][0]/(Nx*Ny);
		gammid_L0midphi[i*Nx+j][1]=0;
		gamplus_L0plusphi[i*Nx+j][0]=Oprt.InterpD.Gam_plus[i][j]*L0plusphi[i*Nx+j][0]/(Nx*Ny);
		gamplus_L0plusphi[i*Nx+j][1]=0;		
	  } 
	}
    
	fftw_execute(plan_fftgammin_phi);
	fftw_execute(plan_fftgammid_phi);
    fftw_execute(plan_fftgamplus_phi);
	fftw_execute(plan_fftgammin_L0minphi);
	fftw_execute(plan_fftgammid_L0midphi);
	fftw_execute(plan_fftgamplus_L0plusphi);	
	
	fftw_destroy_plan(plan_fftgammin_phi);
	fftw_destroy_plan(plan_fftgammid_phi);
    fftw_destroy_plan(plan_fftgamplus_phi);
	fftw_destroy_plan(plan_fftgammin_L0minphi);
	fftw_destroy_plan(plan_fftgammid_L0midphi);
	fftw_destroy_plan(plan_fftgamplus_L0plusphi);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			result[i*Nx+j][0] =(gammin_L0minphihat[i*Nx+j][0]+Oprt_min[i][j]*gammin_phihat[i*Nx+j][0]
									+gammid_L0midphihat[i*Nx+j][0]+Oprt_mid[i][j]*gammid_phihat[i*Nx+j][0]
									+gamplus_L0plusphihat[i*Nx+j][0]+Oprt_plus[i][j]*gamplus_phihat[i*Nx+j][0])/2;
			result[i*Nx+j][1] =(gammin_L0minphihat[i*Nx+j][1]+Oprt_min[i][j]*gammin_phihat[i*Nx+j][1]
									+gammid_L0midphihat[i*Nx+j][1]+Oprt_mid[i][j]*gammid_phihat[i*Nx+j][1]
									+gamplus_L0plusphihat[i*Nx+j][1]+Oprt_plus[i][j]*gamplus_phihat[i*Nx+j][1])/2;
		}
    }

	fftw_free(L0minphihat);fftw_free(L0midphihat);fftw_free(L0plusphihat);
    fftw_free(L0minphi);fftw_free(L0midphi);fftw_free(L0plusphi);
	fftw_free(gammin_phi);fftw_free(gammid_phi);fftw_free(gamplus_phi);
	fftw_free(gammin_L0minphi);fftw_free(gammid_L0midphi);fftw_free(gamplus_L0plusphi);
	fftw_free(gammin_phihat);fftw_free(gammid_phihat);fftw_free(gamplus_phihat);
	fftw_free(gammin_L0minphihat);fftw_free(gammid_L0midphihat);fftw_free(gamplus_L0plusphihat);    
}

void HSS_3depth(fftw_complex* phi, fftw_complex* phihat,double** Oprt_min,double** Oprt_plus,double** Oprt_mid,
double** Gam_min,double** Gam_plus,double** Gam_mid,fftw_complex* L0phihath,int Nx,int Ny,int factIFFT)
{
	fftw_complex *gammin_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammid_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gamplus_phihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammin_L0minphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gammid_L0midphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *gamplus_L0plusphihat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0minphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0midphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0plusphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0minphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0midphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* L0plusphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));		
	
	// Create plans 
	fftw_plan plan_iftL0minphihat   = fftw_plan_dft_2d(Ny,Nx, L0minphihat,  L0minphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftL0midphihat   = fftw_plan_dft_2d(Ny,Nx, L0midphihat,  L0midphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftL0plusphihat  = fftw_plan_dft_2d(Ny,Nx, L0plusphihat, L0plusphi,FFTW_BACKWARD, FFTW_ESTIMATE);
    
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		L0minphihat[i*Nx+j][0]=Oprt_min[i][j]*phihat[i*Nx+j][0];
		L0minphihat[i*Nx+j][1]=Oprt_min[i][j]*phihat[i*Nx+j][1];
		L0midphihat[i*Nx+j][0]=Oprt_mid[i][j]*phihat[i*Nx+j][0];
		L0midphihat[i*Nx+j][1]=Oprt_mid[i][j]*phihat[i*Nx+j][1];
		L0plusphihat[i*Nx+j][0]=Oprt_plus[i][j]*phihat[i*Nx+j][0];
		L0plusphihat[i*Nx+j][1]=Oprt_plus[i][j]*phihat[i*Nx+j][1];
	  }
	}
	
	fftwcomplex_2d_sym(L0minphihat,Nx,Ny);
	fftwcomplex_2d_sym(L0midphihat,Nx,Ny);
	fftwcomplex_2d_sym(L0plusphihat,Nx,Ny);
    fftw_execute(plan_iftL0minphihat);
    fftw_execute(plan_iftL0midphihat);
    fftw_execute(plan_iftL0plusphihat);
    
	fftw_destroy_plan(plan_iftL0minphihat);
	fftw_destroy_plan(plan_iftL0midphihat);
	fftw_destroy_plan(plan_iftL0plusphihat);   
	
	fftw_complex* gammin_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammid_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* gamplus_phi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammin_L0minphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gammid_L0midphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* gamplus_L0plusphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
    // Create plans 
    fftw_plan plan_fftgammin_phi   	   = fftw_plan_dft_2d(Ny,Nx, gammin_phi , gammin_phihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammid_phi   	   = fftw_plan_dft_2d(Ny,Nx, gammid_phi , gammid_phihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgamplus_phi      = fftw_plan_dft_2d(Ny,Nx, gamplus_phi, gamplus_phihat, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammin_L0minphi  = fftw_plan_dft_2d(Ny,Nx, gammin_L0minphi , gammin_L0minphihat,FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_fftgammid_L0midphi  = fftw_plan_dft_2d(Ny,Nx, gammid_L0midphi , gammid_L0midphihat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftgamplus_L0plusphi= fftw_plan_dft_2d(Ny,Nx, gamplus_L0plusphi, gamplus_L0plusphihat, FFTW_FORWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		gammin_phi[i*Nx+j][0]=Gam_min[i][j]*phi[i*Nx+j][0]/factIFFT;
		gammin_phi[i*Nx+j][1]=0;
		gammid_phi[i*Nx+j][0]=Gam_mid[i][j]*phi[i*Nx+j][0]/factIFFT;
		gammid_phi[i*Nx+j][1]=0;
		gamplus_phi[i*Nx+j][0]=Gam_plus[i][j]*phi[i*Nx+j][0]/factIFFT;
		gamplus_phi[i*Nx+j][1]=0;
		gammin_L0minphi[i*Nx+j][0]=Gam_min[i][j]*L0minphi[i*Nx+j][0]/(Nx*Ny);
		gammin_L0minphi[i*Nx+j][1]=0;
		gammid_L0midphi[i*Nx+j][0]=Gam_mid[i][j]*L0midphi[i*Nx+j][0]/(Nx*Ny);
		gammid_L0midphi[i*Nx+j][1]=0;
		gamplus_L0plusphi[i*Nx+j][0]=Gam_plus[i][j]*L0plusphi[i*Nx+j][0]/(Nx*Ny);
		gamplus_L0plusphi[i*Nx+j][1]=0;		
	  } 
	}
    
	fftw_execute(plan_fftgammin_phi);
	fftw_execute(plan_fftgammid_phi);
    fftw_execute(plan_fftgamplus_phi);
	fftw_execute(plan_fftgammin_L0minphi);
	fftw_execute(plan_fftgammid_L0midphi);
	fftw_execute(plan_fftgamplus_L0plusphi);	
	
	fftw_destroy_plan(plan_fftgammin_phi);
	fftw_destroy_plan(plan_fftgammid_phi);
    fftw_destroy_plan(plan_fftgamplus_phi);
	fftw_destroy_plan(plan_fftgammin_L0minphi);
	fftw_destroy_plan(plan_fftgammid_L0midphi);
	fftw_destroy_plan(plan_fftgamplus_L0plusphi);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			L0phihath[i*Nx+j][0] =(gammin_L0minphihat[i*Nx+j][0]+Oprt_min[i][j]*gammin_phihat[i*Nx+j][0]
						+gammid_L0midphihat[i*Nx+j][0]+Oprt_mid[i][j]*gammid_phihat[i*Nx+j][0]
						+gamplus_L0plusphihat[i*Nx+j][0]+Oprt_plus[i][j]*gamplus_phihat[i*Nx+j][0])/2;
			L0phihath[i*Nx+j][1] =(gammin_L0minphihat[i*Nx+j][1]+Oprt_min[i][j]*gammin_phihat[i*Nx+j][1]
						+gammid_L0midphihat[i*Nx+j][1]+Oprt_mid[i][j]*gammid_phihat[i*Nx+j][1]
						+gamplus_L0plusphihat[i*Nx+j][1]+Oprt_plus[i][j]*gamplus_phihat[i*Nx+j][1])/2;
		}
    }

	fftw_free(L0minphihat);fftw_free(L0midphihat);fftw_free(L0plusphihat);
    fftw_free(L0minphi);fftw_free(L0midphi);fftw_free(L0plusphi);
	fftw_free(gammin_phi);fftw_free(gammid_phi);fftw_free(gamplus_phi);
	fftw_free(gammin_L0minphi);fftw_free(gammid_L0midphi);fftw_free(gamplus_L0plusphi);
	fftw_free(gammin_phihat);fftw_free(gammid_phihat);fftw_free(gamplus_phihat);
	fftw_free(gammin_L0minphihat);fftw_free(gammid_L0midphihat);fftw_free(gamplus_L0plusphihat);    
}

void fun_adjust(fftw_complex* AChi_hat,fftw_complex* Ahat, int Nx, int Ny)
{ 
	fftw_complex* A    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* AChi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan ift_Ahat = fftw_plan_dft_2d(Ny,Nx,Ahat,A,FFTW_BACKWARD, FFTW_ESTIMATE);	
	fftw_plan fft_AChi = fftw_plan_dft_2d(Ny,Nx,AChi,AChi_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftwcomplex_2d_sym(Ahat,Nx,Ny);
	fftw_execute(ift_Ahat);
	fftw_destroy_plan(ift_Ahat);
	
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			AChi[i*Nx+j][0] =(A[i*Nx+j][0]*ChiAdj[i][j]/(Nx*Ny));
			AChi[i*Nx+j][1] =(A[i*Nx+j][1]*ChiAdj[i][j]/(Nx*Ny));
		}
    }
	
	fftw_execute(fft_AChi);
	fftw_destroy_plan(fft_AChi);
	
	fftw_free(A);
	fftw_free(AChi);	
}

void fun_sfriction(double t, int Nx,int Ny,int fact,fftw_complex* eta,fftw_complex* u,fftw_complex* v)
{	
	double** Sf_x   = declare_2darray(Nx,Ny); 
	double** Sf_y   = declare_2darray(Nx,Ny); 
	double** absuv  = declare_2darray(Nx,Ny);
	int factIFFT=Nx*Ny;
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			absuv[i][j] = sqrt(pow(u[i*Nx+j][0]/factIFFT,2)+pow(v[i*Nx+j][0]/factIFFT,2));			
		}
	}		
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Sf_x[i][j] = -fric_coef[i][j]*pow((1/(fric_depth[i][j]+eta[i*Nx+j][0]/factIFFT)),4.0/3.0)*u[i*Nx+j][0]/factIFFT*(absuv[i][j]);
			Sf_y[i][j] = -fric_coef[i][j]*pow((1/(fric_depth[i][j]+eta[i*Nx+j][0]/factIFFT)),4.0/3.0)*v[i*Nx+j][0]/factIFFT*(absuv[i][j]);
		}
	}
	fftw_2d_r2c(Sf_xhat,Sf_x,Nx,Ny);
	fftw_2d_r2c(Sf_yhat,Sf_y,Nx,Ny);
	
	free_2darray(Sf_x,Nx,Ny);
	free_2darray(Sf_y,Nx,Ny);
	free_2darray(absuv,Nx,Ny);
}

void fun_gradphi(int Nx,int Ny,const double* Zhat,fftw_complex* gradphi_xh,fftw_complex* gradphi_yh,fftw_complex* gradphi_xhath,fftw_complex* gradphi_yhath)
{	
	fftw_plan plan_iftgradphi_xhat = fftw_plan_dft_2d(Ny,Nx, gradphi_xhath, gradphi_xh, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_iftgradphi_yhat = fftw_plan_dft_2d(Ny,Nx, gradphi_yhath, gradphi_yh, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			gradphi_xhath[i*Nx+j][0]=-dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j];
			gradphi_xhath[i*Nx+j][1]=dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j];
			gradphi_yhath[i*Nx+j][0]=-dominfo.ky[i]*Zhat[3*Nx*Ny+i*Nx+j];
			gradphi_yhath[i*Nx+j][1]=dominfo.ky[i]*Zhat[2*Nx*Ny+i*Nx+j];
		}
	}
	
	fftwcomplex_2d_sym(gradphi_xhath,Nx,Ny);
	fftwcomplex_2d_sym(gradphi_yhath,Nx,Ny);
	fftw_execute(plan_iftgradphi_xhat);
	fftw_execute(plan_iftgradphi_yhat);
	fftw_destroy_plan(plan_iftgradphi_xhat);
	fftw_destroy_plan(plan_iftgradphi_yhat);
}

void fun_source(double t,int Nx,int Ny)
{	
	if (strcmp(wavename,"zero")!=0){
		complex** SkewSource_hat;
		int Nxy;
		if (strcmp(orientation,"vertical")==0){
			Nxy = Ny;
		}
		else {
			Nxy = Nx;
		}
		
		double* INsig_t 	 = new double[Nxy];
		double* INsig_skew_t = new double[Nxy];
		for (int j=0;j<Nxy;j++){
			INsig_t[j] 		= gsl_spline_eval(Spline_Signal[j].spline,t,Spline_Signal[j].acc);
			INsig_skew_t[j] = gsl_spline_eval(Spline_SkewSignal[j].spline,t,Spline_SkewSignal[j].acc);
		}
		
		
		//compute source and skew_source
		double**  Source 	 = declare_2darray(Nx,Ny);
		double**  SkewSource = declare_2darray(Nx,Ny);
		if (strcmp(orientation,"vertical")==0){
			if (strcmp(propagation,"Uni+")==0){
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						Source[i][j] 	= INsig_t[i]*newGam[j][i];
						SkewSource[i][j]= -1*INsig_skew_t[i]*SkewGam[j][i];
					}
				}
			}
			else if (strcmp(propagation,"Uni-")==0){
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						Source[i][j] 	= INsig_t[i]*newGam[j][i];
						SkewSource[i][j]= 1*INsig_skew_t[i]*SkewGam[j][i];
					}
				}
			}
			else {
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						Source[i][j] 	= 2*INsig_t[i]*newGam[j][i];
						SkewSource[i][j]= 0;
					}
				}
			}
		}
		else{
			if (strcmp(propagation,"Uni+")==0){
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						Source[i][j] 	= INsig_t[j]*newGam[j][i];
						SkewSource[i][j]= -1*INsig_skew_t[j]*SkewGam[j][i];
					}
				}
			}
			else if (strcmp(propagation,"Uni-")==0){
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						Source[i][j] 	= INsig_t[j]*newGam[j][i];
						SkewSource[i][j]= 1*INsig_skew_t[j]*SkewGam[j][i];
					}
				}
			}
			else {
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						Source[i][j] 	= 2*INsig_t[j]*newGam[j][i];
						SkewSource[i][j]= 0;
					}
				}
			}
		}

		//for ODE		
		complex** Source0_hat = declare_2darray_complex(Nx,Ny);
		SkewSource_hat = declare_2darray_complex(Nx,Ny);

		fftw_2d_r2c(Source0_hat,Source,Nx,Ny);
		fftw_2d_r2c(SkewSource_hat,SkewSource,Nx,Ny);	
		for (int i=0;i<Ny;i++){
			for (int j=0;j<Nx;j++){
				Source_hat[i][j].Re = Source0_hat[i][j].Re+SkewSource_hat[i][j].Re;
				Source_hat[i][j].Im = Source0_hat[i][j].Im+SkewSource_hat[i][j].Im;
			}
		}

			
		delete[] INsig_t;
		delete[] INsig_skew_t;
		free_2darray(Source,Nx,Ny);
		free_2darray(SkewSource,Nx,Ny);
		free_2darray_complex(Source0_hat,Nx,Ny);
		free_2darray_complex(SkewSource_hat,Nx,Ny);				
	}
}

void fun_sbreaking(double t,fftw_complex* eta,complex** delphiH_hat,int Nx,int Ny)
{
	double**  delphi   = declare_2darray(Nx,Ny);
	double**  dxFlux   = declare_2darray(Nx,Ny);
	double**  dyFlux   = declare_2darray(Nx,Ny);
	double**  Flux     = declare_2darray(Nx,Ny);
	int fact=Nx*Ny;

	ifftw_2d_c2r(delphi,delphiH_hat,Nx,Ny);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Flux[i][j] = -B[i][j]*pow(pb.delb,2)*(eta[i*Nx+j][0]/fact-dominfo.bathy.data[i][j])*pow(delphi[i][j],2);
		}
	}	
	fun_gradient(dxFlux,dyFlux,Flux,Nx,Ny,dominfo.dx,dominfo.dy);	
	
	double**  Sb_x   = declare_2darray(Nx,Ny);
	double**  Sb_y   = declare_2darray(Nx,Ny);
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Sb_x[i][j] = dxFlux[i][j]/(eta[i*Nx+j][0]/fact-dominfo.bathy.data[i][j]);
			Sb_y[i][j] = dyFlux[i][j]/(eta[i*Nx+j][0]/fact-dominfo.bathy.data[i][j]);
		}
	}
	fftw_2d_r2c(Sb_xhat,Sb_x,Nx,Ny);
	fftw_2d_r2c(Sb_yhat,Sb_y,Nx,Ny);
	
	free_2darray(Flux,Nx,Ny);
	free_2darray(dxFlux,Nx,Ny);free_2darray(dyFlux,Nx,Ny);
	free_2darray(Sb_x,Nx,Ny);free_2darray(Sb_y,Nx,Ny);
	free_2darray(delphi,Nx,Ny);			
}

void fun_swall(double t,double dx,double dy,int Nx,int Ny,const double* Zhat)
{
	fftw_complex* eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	get_eta(eta,t,Zhat,Nx,Ny);	
	
	if ((t>tprevW)){	
		for (int ii=0;ii<n_idxwall;ii++){
			Swall_ts[ii]=(refl_coef*eta[(int)(idxWall[ii].y*Nx+idxWall[ii].x)][0]/(Nx*Ny));
			Swt_skewline[ii]=Swt_skewline[ii]+Swall_ts[ii]*(t-tprevW);
		}		
		iterW	= iterW+1;
		tprevW	= t;
	}
	
	//Swt_skew*SkewGam
	double**  temp = declare_2darray(Nx,Ny);
	complex** temp_hat;
	set_matrix_val(temp,Nx,Ny,0);
	for (int iw=0;iw<4;iw++){
		if (iw==0){//skew line Ny
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					temp[i][j] = temp[i][j]+2*Swt_skewline[i]*(dxGam1[i][j]*dy)/(4*pow(M_PI,2));			
				}
			}
		}
		else if(iw==1){
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					temp[i][j] = temp[i][j]+2*Swt_skewline[i+Ny]*(dxGam2[i][j]*dy)/(4*pow(M_PI,2));			
				}
			}
		}
		else if(iw==2){
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					//temp[i][j] = temp[i][j]+2*Swt_skewline[j+2*Ny]*(dyGam1[j][i]*dx)/(4*pow(M_PI,2));			
				}
			}
		}
		else if(iw==3){
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					//temp[i][j] = temp[i][j]+2*Swt_skewline[j+2*Ny+Nx]*(dyGam2[j][i]*dx)/(4*pow(M_PI,2));			
				}
			}
		}
	}
	temp_hat = declare_2darray_complex(Nx,Ny);
	fftw_2d_r2c(temp_hat,temp,Nx,Ny);
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Sw_hat[i][j].Re = temp_hat[i][j].Re;
			Sw_hat[i][j].Im = temp_hat[i][j].Im;
		}
	}

	//update iter
	if (t>tprevW){
		
	}

	free_2darray(temp,Nx,Ny);
	free_2darray_complex(temp_hat,Nx,Ny);
	fftw_free(eta);
}

void fun_damping(double** dampcharH,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* dampetahat,
fftw_complex* dampuhat,fftw_complex* dampvhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;

	fftw_complex* dampeta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
	// Create plans 
	fftw_plan plan_fftdampeta = fftw_plan_dft_2d(Ny,Nx,dampeta,dampetahat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampu   = fftw_plan_dft_2d(Ny,Nx,dampu,dampuhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampv   = fftw_plan_dft_2d(Ny,Nx,dampv,dampvhat, FFTW_FORWARD, FFTW_ESTIMATE);

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dampeta[i*Nx+j][0]=Oprt.Cpeak*dampcharH[i][j]*eta[i*Nx+j][0]/factIFFT;
			dampeta[i*Nx+j][1]=0;
			dampu[i*Nx+j][0]=Oprt.Cpeak*dampcharH[i][j]*u[i*Nx+j][0]/factIFFT;
			dampu[i*Nx+j][1]=0;
			dampv[i*Nx+j][0]=Oprt.Cpeak*dampcharH[i][j]*v[i*Nx+j][0]/factIFFT;
			dampv[i*Nx+j][1]=0;
		}
	}
	
	fftw_execute(plan_fftdampeta);
	fftw_execute(plan_fftdampu); 
	fftw_execute(plan_fftdampv);      
	
    fftw_destroy_plan(plan_fftdampeta);
    fftw_destroy_plan(plan_fftdampu);
    fftw_destroy_plan(plan_fftdampv);

    fftw_free(dampeta);
    fftw_free(dampu);
    fftw_free(dampv);
}

void set_init_ode(wavestruc waveinit,int Nx, int Ny)
{
	/*if (runup_id==0){//runup
		//total initial eta0
		double 	 Hhere;		
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				Hhere = waveinit.profile[i][j]-dominfo.bathy.data[i][j];
				if (Hhere<H_minshore){
					waveinit.profile[i][j] = dominfo.bathy.data[i][j];
				}
			}
		}
		//initial phi0
	}	*/
	
	complex** eta_hat = declare_2darray_complex(Nx,Ny);
	complex** u_hat = declare_2darray_complex(Nx,Ny);
	complex** v_hat = declare_2darray_complex(Nx,Ny);
	fftw_2d_r2c(eta_hat,waveinit.profile,Nx,Ny);	
	fftw_2d_r2c(u_hat,waveinit.u,Nx,Ny);	
	fftw_2d_r2c(v_hat,waveinit.v,Nx,Ny);
		
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			y_hat[i*Nx+j]			=eta_hat[i][j].Re;
			y_hat[Nx*Ny+i*Nx+j]		=eta_hat[i][j].Im;
			y_hat[2*Nx*Ny+i*Nx+j]	=u_hat[i][j].Re;
			y_hat[3*Nx*Ny+i*Nx+j]	=u_hat[i][j].Im;
			y_hat[4*Nx*Ny+i*Nx+j]	=v_hat[i][j].Re;
			y_hat[5*Nx*Ny+i*Nx+j]	=v_hat[i][j].Im;
		}
	}	

	free_2darray_complex(eta_hat,Nx,Ny);
	free_2darray_complex(u_hat,Nx,Ny);
	free_2darray_complex(v_hat,Nx,Ny);
	
	//this is if runup exist
	/*if (runup == 0){//runup yes		
		//IF deep
		if (inf_category==1){
			complex* uhat 	  = new complex[nx];
			double*  ift_uhat = new double[nx];					
			for (int j=0;j<nx;j++){
				uhat[j].Re = y[j+2*nx];
				uhat[j].Im = y[j+3*nx];
			}
			ifft_complex_real_nx(uhat,ift_uhat);
						
			double* utot   = new double[nx];
			complex* utemp = new complex[nx];
			double  ur;
			double* grad_depth = new double[nx];
			gradient(depth,dx,nx,grad_depth);
			
			for (int j=0;j<nx;j++){
				Hhere = ift_uhat[j]+depth[j];
				if (Hhere<(Hmin_interp/5)){
					ur = -grad_depth[j];
				}
				else {
					ur = 0;
				}
				utot[j] = ift_uhat[j]+ur;
			}
			fft_real_complex_nx(utot,utemp);
		
			for (int j=0;j<nx;j++){
				y[j+2*nx] = utemp[j].Re;
				y[j+3*nx] = utemp[j].Im;
			}
			delete[] uhat;delete[] ift_uhat;delete[] utot;delete[] utemp;	
		}
	}*/			
	
}

/*void fun_Lphihat_if1(double t,int Nx,int Ny,const double* Zhat,fftw_complex* L0phihat)
{
	//allocation memory
	fftw_complex* gradphi_x=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gradphi_y=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gradphi_xhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gradphi_yhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* C2u_xhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* C2u_yhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* C2u_x 	   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* C2u_y 	   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	double** Csqr_min = declare_2darray(Nx,Ny);
	double** Csqr_plus= declare_2darray(Nx,Ny);
	double** Csqr_mid = declare_2darray(Nx,Ny);
	
	int factIFFT=Nx*Ny;

	fun_gradphi(Nx,Ny,Zhat,gradphi_x,gradphi_y,gradphi_xhat,gradphi_yhat);	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Csqr_min[i][j]  = pow(Oprt.InterpD.C2d_min[i][j],2);
			Csqr_plus[i][j] = pow(Oprt.InterpD.C2d_plus[i][j],2);
		}
	}
	
	if (ninterp==2){
		HSSgen_2depth(gradphi_x,gradphi_xhat,Csqr_min,Csqr_plus,C2u_xhat,Nx,Ny,factIFFT);	
		HSSgen_2depth(gradphi_y,gradphi_yhat,Csqr_min,Csqr_plus,C2u_yhat,Nx,Ny,factIFFT);		
	}
	else {
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				Csqr_mid[i][j]  = pow(Oprt.InterpD.C2d_mid[i][j],2);
			}
		}
		HSSgen_3depth(gradphi_x,gradphi_xhat,Csqr_min,Csqr_mid,Csqr_plus,C2u_xhat,Nx,Ny,factIFFT);
		HSSgen_3depth(gradphi_y,gradphi_yhat,Csqr_min,Csqr_mid,Csqr_plus,C2u_yhat,Nx,Ny,factIFFT);
	}	
	
	//plan and destroy
	fftw_plan plan_iftC2ux_hat = fftw_plan_dft_2d(Ny,Nx, C2u_xhat, C2u_x, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftC2uy_hat = fftw_plan_dft_2d(Ny,Nx, C2u_yhat, C2u_y, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftwcomplex_2d_sym(C2u_xhat,Nx,Ny);
	fftwcomplex_2d_sym(C2u_yhat,Nx,Ny);
	fftw_execute(plan_iftC2ux_hat);
    fftw_execute(plan_iftC2uy_hat);    
    fftw_destroy_plan(plan_iftC2ux_hat);
    fftw_destroy_plan(plan_iftC2uy_hat);
    
    complex **divC2u_hat  = declare_2darray_complex(Nx,Ny);
	double **divC2u 	  = declare_2darray(Nx,Ny);
	double **divC2u_wall  = declare_2darray(Nx,Ny);
    double **Flux_grad    = declare_2darray(Nx,Ny);
    
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			divC2u_hat[i][j].Re = -(-dominfo.kx[j]*C2u_xhat[i*Nx+j][1]-dominfo.ky[i]*C2u_yhat[i*Nx+j][1]);
			divC2u_hat[i][j].Im = -(dominfo.kx[j]*C2u_xhat[i*Nx+j][0]+dominfo.ky[i]*C2u_yhat[i*Nx+j][0]);
		}
	}
	ifftw_2d_c2r(divC2u,divC2u_hat,Nx,Ny);			
			
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			divC2u_wall[i][j] = divC2u[i][j]*dominfo.wallchar[i][j];
			Flux_grad[i][j]   = C2u_x[i*Nx+j][0]/factIFFT*grad_wallchar_x[i][j]+C2u_y[i*Nx+j][0]/factIFFT*grad_wallchar_y[i][j];
		}
	}
	complex** divC2u_wall_hat = fftw_2d_r2c(divC2u_wall,Nx,Ny);
	complex** Flux_grad_hat	  = fftw_2d_r2c(Flux_grad,Nx,Ny);	
			
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			L0phihat[i*Nx+j][0] = 1/grav*(divC2u_wall_hat[i][j].Re-Flux_grad_hat[i][j].Re);
			L0phihat[i*Nx+j][1] = 1/grav*(divC2u_wall_hat[i][j].Im-Flux_grad_hat[i][j].Im);
		}
	}
	
	free_2darray_complex(Flux_grad_hat,Nx,Ny);
	free_2darray_complex(divC2u_wall_hat,Nx,Ny);
	free_2darray_complex(divC2u_hat,Nx,Ny);	
	free_2darray(divC2u,Nx,Ny);
	free_2darray(divC2u_wall,Nx,Ny);
	free_2darray(Flux_grad,Nx,Ny);
	free_2darray(Csqr_min,Nx,Ny);
	free_2darray(Csqr_plus,Nx,Ny);
	free_2darray(Csqr_mid,Nx,Ny);
	
	fftw_free(gradphi_x);
	fftw_free(gradphi_y);
	fftw_free(gradphi_xhat);
	fftw_free(gradphi_yhat);
	fftw_free(C2u_xhat);
	fftw_free(C2u_yhat);
	fftw_free(C2u_x);
	fftw_free(C2u_y);
}

void fun_Lphihat(double t,int Nx,int Ny,const double* Zhat,fftw_complex* L0phihat)
{
	//allocation memory
	fftw_complex* gradphi_x=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gradphi_y=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gradphi_xhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* gradphi_yhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* Cu_xhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_yhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_x 	   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_y 	   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_xWall = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_yWall = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_xWall_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cu_yWall_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* C2u_xhat_w=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* C2u_yhat_w=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	
	int factIFFT=Nx*Ny;

	fun_gradphi(Nx,Ny,Zhat,gradphi_x,gradphi_y,gradphi_xhat,gradphi_yhat);	
		
	if (ninterp==2){
		HSSgen_2depth(gradphi_x,gradphi_xhat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_plus,Cu_xhat,Nx,Ny,factIFFT);	
		HSSgen_2depth(gradphi_y,gradphi_yhat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_plus,Cu_yhat,Nx,Ny,factIFFT);		
	}
	else {
		HSSgen_3depth(gradphi_x,gradphi_xhat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_mid,Oprt.InterpD.C2d_plus,Cu_xhat,Nx,Ny,factIFFT);
		HSSgen_3depth(gradphi_y,gradphi_yhat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_mid,Oprt.InterpD.C2d_plus,Cu_yhat,Nx,Ny,factIFFT);
	}

	//plan and destroy
	fftw_plan plan_iftCux_hat = fftw_plan_dft_2d(Ny,Nx, Cu_xhat, Cu_x, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftCuy_hat = fftw_plan_dft_2d(Ny,Nx, Cu_yhat, Cu_y, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftwcomplex_2d_sym(Cu_xhat,Nx,Ny);
	fftwcomplex_2d_sym(Cu_yhat,Nx,Ny);
	fftw_execute(plan_iftCux_hat);
    fftw_execute(plan_iftCuy_hat);    
    fftw_destroy_plan(plan_iftCux_hat);
    fftw_destroy_plan(plan_iftCuy_hat);
		
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Cu_xWall[i*Nx+j][0] = Cu_x[i*Nx+j][0]*dominfo.wallchar[i][j]/factIFFT;
			Cu_xWall[i*Nx+j][1] = 0;
			Cu_yWall[i*Nx+j][0] = Cu_y[i*Nx+j][0]*dominfo.wallchar[i][j]/factIFFT;
			Cu_yWall[i*Nx+j][1] = 0;
		}
	}

	//plan and destroy
	fftw_plan plan_ftCu_xWall = fftw_plan_dft_2d(Ny,Nx,Cu_xWall,Cu_xWall_hat,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan plan_ftCu_yWall = fftw_plan_dft_2d(Ny,Nx,Cu_yWall,Cu_yWall_hat,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(plan_ftCu_xWall);
    fftw_execute(plan_ftCu_yWall);    
    fftw_destroy_plan(plan_ftCu_xWall);
    fftw_destroy_plan(plan_ftCu_yWall);
		
	if (ninterp==2){
		HSSgen_2depth(Cu_xWall,Cu_xWall_hat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_plus,C2u_xhat_w,Nx,Ny,1);	
		HSSgen_2depth(Cu_yWall,Cu_yWall_hat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_plus,C2u_yhat_w,Nx,Ny,1);		
	}
	else {
		HSSgen_3depth(Cu_xWall,Cu_xWall_hat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_mid,Oprt.InterpD.C2d_plus,C2u_xhat_w,Nx,Ny,1);
		HSSgen_3depth(Cu_yWall,Cu_yWall_hat,Oprt.InterpD.C2d_min,Oprt.InterpD.C2d_mid,Oprt.InterpD.C2d_plus,C2u_yhat_w,Nx,Ny,1);
	}
		
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			L0phihat[i*Nx+j][0]= -1/grav*(-dominfo.kx[j]*C2u_xhat_w[i*Nx+j][1]-dominfo.ky[i]*C2u_yhat_w[i*Nx+j][1]);
			L0phihat[i*Nx+j][1]= -1/grav*(dominfo.kx[j]*C2u_xhat_w[i*Nx+j][0]+dominfo.ky[i]*C2u_yhat_w[i*Nx+j][0]);
		}
	}

	if (t>10){
		double** tes=declare_2darray(Nx,Ny);
		double* xx= new double[Nx];
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				tes[i][j] = L0phihat[i*Nx+j][0];
				xx[j] = j;
			}
		}
		FILE* gp1=popen("gnuplot -persistent","w");
		fun_plot_2d_trarray(gp1,dominfo.y,xx,tes,dominfo.Ny,dominfo.Nx);
		fflush(gp1);exit(0);
		pclose(gp1);
		free_2darray(tes,Nx,Ny);
	}	
		
	fftw_free(gradphi_x);
	fftw_free(gradphi_y);
	fftw_free(gradphi_xhat);
	fftw_free(gradphi_yhat);
	fftw_free(Cu_xhat);
	fftw_free(Cu_yhat);
	fftw_free(Cu_x);
	fftw_free(Cu_y);
	fftw_free(Cu_xWall);
	fftw_free(Cu_yWall);
	fftw_free(Cu_xWall_hat);
	fftw_free(Cu_yWall_hat);
	fftw_free(C2u_xhat_w);
	fftw_free(C2u_yhat_w);
}*/

void get_u(fftw_complex* phi, double t, const double* Zhat, int Nx, int Ny)
{	
	//Define phi
	fftw_complex* phihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
	// Create plans 
	fftw_plan plan_iftphihat   = fftw_plan_dft_2d(Ny,Nx, phihat, phi, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			phihat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			phihat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];	
		}
	}	
	fftwcomplex_2d_sym(phihat,Nx,Ny);
	fftw_execute(plan_iftphihat);
	fftw_destroy_plan(plan_iftphihat);
	fftw_free(phihat);
}

void  get_v(fftw_complex* phi, double t, const double* Zhat, int Nx, int Ny)
{	
	//Define phi
	fftw_complex* phihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
	// Create plans 
	fftw_plan plan_iftphihat   = fftw_plan_dft_2d(Ny,Nx, phihat, phi, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			phihat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			phihat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];	
		}
	}	
	fftwcomplex_2d_sym(phihat,Nx,Ny);
	fftw_execute(plan_iftphihat);
	fftw_destroy_plan(plan_iftphihat);
	fftw_free(phihat);
}

void fun_termshat_shore_uv_HSdirect(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* deluH1_hat,
fftw_complex* delvH1_hat,fftw_complex* deletaH1_hat,fftw_complex* deluH2_hat,fftw_complex* delvH2_hat,fftw_complex* deletaH2_hat,
domvar dominfo,const double* Zhat,double** WaveChar,double** H)
{
	int Nx = dominfo.Nx;
	int Ny = dominfo.Ny;
	int factIFFT=Nx*Ny;	
	
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){	
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];	
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];	
		}
	}
	fftw_complex* Cp2u_hat1= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2u_hat2= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2u     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2u_Xhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_Cp2u    		= fftw_plan_dft_2d(Ny,Nx,Cp2u_hat,Cp2u,FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_Cp2u_Xhat    = fftw_plan_dft_2d(Ny,Nx,Cp2u,Cp2u_Xhat,FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_complex* Cp2v_hat1= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2v_hat2= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2v_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2v     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cp2v_Yhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_Cp2v    		= fftw_plan_dft_2d(Ny,Nx,Cp2v_hat,Cp2v,FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_Cp2v_Yhat    = fftw_plan_dft_2d(Ny,Nx,Cp2v,Cp2v_Yhat,FFTW_FORWARD, FFTW_ESTIMATE);
	
	//prepare gamma
	double** gam_min1  = declare_2darray(Nx,Ny);
	double** gam_mid1  = declare_2darray(Nx,Ny);
	double** gam_plus1 = declare_2darray(Nx,Ny);
	double** gam_min2  = declare_2darray(Nx,Ny);
	double** gam_plus2 = declare_2darray(Nx,Ny);
	double** gam_mid2  = declare_2darray(Nx,Ny);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			gam_min1[i][j]  = gsl_spline_eval(sg_m1,H[i][j],acc_m1);
			gam_mid1[i][j]  = gsl_spline_eval(sg_c1,H[i][j],acc_c1);
			gam_plus1[i][j] = gsl_spline_eval(sg_p1,H[i][j],acc_p1);
			gam_min2[i][j]  = gsl_spline_eval(sg_m2,H[i][j],acc_m2);
			gam_mid2[i][j]  = gsl_spline_eval(sg_c2,H[i][j],acc_c2);
			gam_plus2[i][j] = gsl_spline_eval(sg_p2,H[i][j],acc_p2);
		}
	}	
	HSS_3depth(u,uhat,C2m1,C2p1,C2c1,gam_min1,gam_plus1,gam_mid1,Cp2u_hat1,Nx,Ny,factIFFT);
	HSS_3depth(u,uhat,C2m2,C2p2,C2c2,gam_min2,gam_plus2,gam_mid2,Cp2u_hat2,Nx,Ny,factIFFT);
	
	HSS_3depth(v,vhat,C2m1,C2p1,C2c1,gam_min1,gam_plus1,gam_mid1,Cp2v_hat1,Nx,Ny,factIFFT);
	HSS_3depth(v,vhat,C2m2,C2p2,C2c2,gam_min2,gam_plus2,gam_mid2,Cp2v_hat2,Nx,Ny,factIFFT);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Cp2u_hat[i*Nx+j][0] = Cp2u_hat1[i*Nx+j][0]+Cp2u_hat2[i*Nx+j][0];
			Cp2u_hat[i*Nx+j][1] = Cp2u_hat1[i*Nx+j][1]+Cp2u_hat2[i*Nx+j][1];
			
			Cp2v_hat[i*Nx+j][0] = Cp2v_hat1[i*Nx+j][0]+Cp2v_hat2[i*Nx+j][0];
			Cp2v_hat[i*Nx+j][1] = Cp2v_hat1[i*Nx+j][1]+Cp2v_hat2[i*Nx+j][1];
		}
	}	
	
	fftw_execute(plan_Cp2u);
	fftw_destroy_plan(plan_Cp2u);	
	fftw_execute(plan_Cp2v);
	fftw_destroy_plan(plan_Cp2v);	

	//update
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (H[i][j]<H_minShore){
				Cp2u[i*Nx+j][0] = grav*H[i][j]*u[i*Nx+j][0];
				Cp2v[i*Nx+j][0] = grav*H[i][j]*v[i*Nx+j][0];
			}
		}
	}
	fftw_execute(plan_Cp2u_Xhat);
	fftw_destroy_plan(plan_Cp2u_Xhat);
	fftw_execute(plan_Cp2v_Yhat);
	fftw_destroy_plan(plan_Cp2v_Yhat);

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			deluH1_hat[i*Nx+j][0] = -(-dominfo.kx[j]*Cp2u_Xhat[i*Nx+j][1])/grav;
			deluH1_hat[i*Nx+j][1] = -(dominfo.kx[j]*Cp2u_Xhat[i*Nx+j][0])/grav;
			delvH1_hat[i*Nx+j][0] = -(-dominfo.ky[i]*Cp2v_Yhat[i*Nx+j][1])/grav;
			delvH1_hat[i*Nx+j][1] = -(dominfo.ky[i]*Cp2v_Yhat[i*Nx+j][0])/grav;
			deletaH1_hat[i*Nx+j][0]	= -(etahat[i*Nx+j][0])*grav;
			deletaH1_hat[i*Nx+j][1]	= -(etahat[i*Nx+j][1])*grav;
		}
	}
	
	fftw_free(etahat);
	fftw_free(uhat);fftw_free(vhat);
	fftw_free(Cp2u_hat1);
	fftw_free(Cp2u_hat2);
	fftw_free(Cp2u_hat);
	fftw_free(Cp2u);
	fftw_free(Cp2u_Xhat);
	
	fftw_free(Cp2v_hat1);
	fftw_free(Cp2v_hat2);
	fftw_free(Cp2v_hat);
	fftw_free(Cp2v);
	fftw_free(Cp2v_Yhat);
	free_2darray(gam_min1,Nx,Ny);
	free_2darray(gam_mid1,Nx,Ny);
	free_2darray(gam_plus1,Nx,Ny);
	free_2darray(gam_min2,Nx,Ny);
	free_2darray(gam_mid2,Nx,Ny);
	free_2darray(gam_plus2,Nx,Ny);
}

void fun_Oprt_u(fftw_complex* u,fftw_complex* uhat,fftw_complex* result,fftw_complex* result_hat,
double** Oprt_m1,double** Oprt_p1,double** Oprt_c1,
double** gam_min1,double** gam_plus1,double** gam_mid1,
double** Oprt_m2,double** Oprt_p2,double** Oprt_c2,
double** gam_min2,double** gam_plus2,double** gam_mid2,int Nx,int Ny,int factIFFT)
{
	fftw_complex* Oprt_u_hat1= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Oprt_u_hat2= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	HSS_3depth(u,uhat,Oprt_m1,Oprt_p1,Oprt_c1,gam_min1,gam_plus1,gam_mid1,Oprt_u_hat1,Nx,Ny,factIFFT);
	HSS_3depth(u,uhat,Oprt_m2,Oprt_p2,Oprt_c2,gam_min2,gam_plus2,gam_mid2,Oprt_u_hat2,Nx,Ny,factIFFT);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			result_hat[i*Nx+j][0] = Oprt_u_hat1[i*Nx+j][0]+Oprt_u_hat2[i*Nx+j][0];
			result_hat[i*Nx+j][1] = Oprt_u_hat1[i*Nx+j][1]+Oprt_u_hat2[i*Nx+j][1];
		}
	}	
	fftw_plan plan_result = fftw_plan_dft_2d(Ny,Nx,result_hat,result,FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_result);
	fftw_destroy_plan(plan_result);	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			result[i*Nx+j][0]=result[i*Nx+j][0]/(factIFFT);
		}
	}
	
	fftw_free(Oprt_u_hat1);	
	fftw_free(Oprt_u_hat2);
}

void fun_termshat_shore_uv(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* delphiH_hat,
fftw_complex* ikxdeletaH23adj_hat,fftw_complex* ikydeletaH23adj_hat,domvar dominfo,const double* Zhat,double** WaveChar,double** H)
{
	int Nx = dominfo.Nx;
	int Ny = dominfo.Ny;
	int factIFFT=Nx*Ny;	

	//allocation memory
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//CuChar and CvChar
	fftw_complex* Cuhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* Cvhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* Cu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* Cv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* CuCharhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* CvCharhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* CuChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* CvChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* C2uCharhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* C2vCharhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* C2uChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* C2vChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//delphiH1 
	fftw_complex* delphiH1_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//C2u and C2v
	fftw_complex* C2uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* C2vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* C2u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));   
	fftw_complex* C2v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* M0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaM0uvChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* etaM0uvChar_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//LetaM0uvChar_hat        = funOprt_L2d_runup(g,Om2dsqInterp,etaM0uvChar_hat,etaM0uvChar);
	fftw_complex* LetaM0uvChar_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
	fftw_complex* LetaM0uvChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//etagradphi_hat.x_hat   = fft2(eta.*gradphi.x);
    //etagradphi_hat.y_hat   = fft2(eta.*gradphi.y);
	fftw_complex* etauhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
	fftw_complex* etau=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* etavhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
	fftw_complex* etav=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//divetagradphi_hat      = funOprt_div2d(dom.Kx,dom.Ky,etagradphi_hat);
	fftw_complex* divetagradphi_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//delphiH2_hat           = -(LetaM0uvChar_hat+divetagradphi_hat);
    //deletaH2_hat           = 0.5.*fft2((absuv2-M0uv.^2).*HeavChar);
	fftw_complex* delphiH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* deletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* deletaH2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* delphiH3_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
	fftw_complex* deletaH3_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	//nonlinear adjustment
	fftw_complex* delphiH23adj_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* delphiH23_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikxdeletaH23_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH23_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	//for first order
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){	
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];	
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];	
		}
	}	
	fun_Oprt_u(u,uhat,Cu,Cuhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	fun_Oprt_u(v,vhat,Cv,Cvhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			CuChar[i*Nx+j][0] = Cu[i*Nx+j][0]*dominfo.wallchar[i][j];
			CuChar[i*Nx+j][1] = 0;
			CvChar[i*Nx+j][0] = Cv[i*Nx+j][0]*dominfo.wallchar[i][j];
			CvChar[i*Nx+j][1] = 0;
		}
	}	
	fftw_plan plan_CuCharhat = fftw_plan_dft_2d(Ny,Nx,CuChar,CuCharhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_CvCharhat = fftw_plan_dft_2d(Ny,Nx,CvChar,CvCharhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_CuCharhat);
	fftw_execute(plan_CvCharhat);
	fftw_destroy_plan(plan_CuCharhat);	
	fftw_destroy_plan(plan_CvCharhat);	

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			CuChar[i*Nx+j][0] = CuChar[i*Nx+j][0]*factIFFT;//for operator
			CvChar[i*Nx+j][0] = CvChar[i*Nx+j][0]*factIFFT;//for operator
		}
	}	
	fun_Oprt_u(CuChar,CuCharhat,C2uChar,C2uCharhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	fun_Oprt_u(CvChar,CvCharhat,C2vChar,C2vCharhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);	
	fun_div(delphiH1_hat,C2uCharhat,C2vCharhat,-1.0/grav,Nx,Ny);	

	//for second order
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if(-dominfo.bathy.data[i][j]<Dmin_input){
				eta[i*Nx+j][0] = H[i][j]*factIFFT;
			}
		}
	}
	
	fun_Oprt_u(u,uhat,C2u,C2uhat,C2m1,C2p1,C2c1,gam_m1,gam_p1,gam_c1,C2m2,C2p2,C2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	fun_Oprt_u(v,vhat,C2v,C2vhat,C2m1,C2p1,C2c1,gam_m1,gam_p1,gam_c1,C2m2,C2p2,C2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);	  
	fun_div(M0uv_hat,C2uhat,C2vhat,-1.0/grav,Nx,Ny);

	fftw_plan plan_M0uv = fftw_plan_dft_2d(Ny,Nx,M0uv_hat,M0uv,FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_M0uv);
	fftw_destroy_plan(plan_M0uv);

	fftw_complex* etaM0uvChar_in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etaM0uvChar[i*Nx+j][0] = eta[i*Nx+j][0]/(factIFFT)*M0uv[i*Nx+j][0]/(factIFFT)*dominfo.wallchar[i][j];//normalized 
			etaM0uvChar[i*Nx+j][1] = 0;
			etaM0uvChar_in[i*Nx+j][0] = etaM0uvChar[i*Nx+j][0]*factIFFT;//for input operator 
			etaM0uvChar_in[i*Nx+j][1] = 0;
		}
	}	
	fftw_plan plan_etaM0uvChar_hat = fftw_plan_dft_2d(Ny,Nx,etaM0uvChar,etaM0uvChar_hat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_etaM0uvChar_hat);
	fftw_destroy_plan(plan_etaM0uvChar_hat);

	fun_Oprt_u(etaM0uvChar_in,etaM0uvChar_hat,LetaM0uvChar,LetaM0uvChar_hat,Om2m1,Om2p1,Om2c1,gam_m1,gam_p1,gam_c1,Om2m2,Om2p2,Om2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);//L need to be normalized by grav
	fftw_free(etaM0uvChar_in);	

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			LetaM0uvChar[i*Nx+j][0] = LetaM0uvChar[i*Nx+j][0]/grav;
			LetaM0uvChar[i*Nx+j][1] = LetaM0uvChar[i*Nx+j][1]/grav;
			LetaM0uvChar_hat[i*Nx+j][0] = LetaM0uvChar_hat[i*Nx+j][0]/grav;
			LetaM0uvChar_hat[i*Nx+j][1] = LetaM0uvChar_hat[i*Nx+j][1]/grav;
		}
	}	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etau[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*u[i*Nx+j][0]/factIFFT;
			etau[i*Nx+j][1] = 0;
			etav[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*v[i*Nx+j][0]/factIFFT;
			etav[i*Nx+j][1] = 0;
		}
	}	
	fftw_plan plan_etauhat = fftw_plan_dft_2d(Ny,Nx,etau,etauhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_etauhat);
	fftw_destroy_plan(plan_etauhat);
	fftw_plan plan_etavhat = fftw_plan_dft_2d(Ny,Nx,etav,etavhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_etavhat);
	fftw_destroy_plan(plan_etavhat);
	
	fun_div(divetagradphi_hat,etauhat,etavhat,1.0,Nx,Ny);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH2_hat[i*Nx+j][0] = -(LetaM0uvChar_hat[i*Nx+j][0]+divetagradphi_hat[i*Nx+j][0]);
			delphiH2_hat[i*Nx+j][1] = -(LetaM0uvChar_hat[i*Nx+j][1]+divetagradphi_hat[i*Nx+j][1]);
			deletaH2[i*Nx+j][0] = 0.5*(pow(u[i*Nx+j][0]/factIFFT,2)+pow(v[i*Nx+j][0]/factIFFT,2)-pow(M0uv[i*Nx+j][0]/factIFFT,2))*dominfo.wallchar[i][j];
			deletaH2[i*Nx+j][1] = 0;
		}
	}	
	fftw_plan plan_deletaH2_hat = fftw_plan_dft_2d(Ny,Nx,deletaH2,deletaH2_hat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_deletaH2_hat);
	fftw_destroy_plan(plan_deletaH2_hat);

	//for third order	
	if(strcmp(evolinfo.model,"HS3")!=0){
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				delphiH3_hat[i*Nx+j][0] = 0;
				delphiH3_hat[i*Nx+j][1] = 0;
				deletaH3_hat[i*Nx+j][0] = 0;
				deletaH3_hat[i*Nx+j][1] = 0;
			}
		}
	}
	else{
		//etaM0uv=eta.*M0uv;
		fftw_complex* etaM0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		fftw_complex* etaM0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));  		  
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				etaM0uv[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*M0uv[i*Nx+j][0]/factIFFT;
				etaM0uv[i*Nx+j][1] = 0;
			}
		}
		fftw_plan plan_etaM0uv_hat = fftw_plan_dft_2d(Ny,Nx,etaM0uv,etaM0uv_hat,FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_etaM0uv_hat);
		fftw_destroy_plan(plan_etaM0uv_hat);
		
		//LetaM0uv_hat        =funOprt_L2d_runup(g,Om2dsqInterp,fft2(etaM0uv),etaM0uv);
        //LetaM0uv            =funC_ifft2(LetaM0uv_hat);
        fftw_complex* LetaM0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
		fftw_complex* LetaM0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				etaM0uv[i*Nx+j][0] = etaM0uv[i*Nx+j][0]*factIFFT; //for operator
			}
		}
		fun_Oprt_u(etaM0uv,etaM0uv_hat,LetaM0uv,LetaM0uv_hat,Om2m1,Om2p1,Om2c1,gam_m1,gam_p1,gam_c1,Om2m2,Om2p2,Om2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
		
		//etaCharLetaM0uv     =eta.*HeavChar.*LetaM0uv;
		fftw_complex* etaCharLetaM0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		fftw_complex* etaCharLetaM0uvhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				LetaM0uv[i*Nx+j][0] = LetaM0uv[i*Nx+j][0]/grav;//normalized by grav
				
				LetaM0uv_hat[i*Nx+j][0]=LetaM0uv_hat[i*Nx+j][0]/grav;//normalized by grav
				LetaM0uv_hat[i*Nx+j][1]=LetaM0uv_hat[i*Nx+j][1]/grav;
				
				etaCharLetaM0uv[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*LetaM0uv[i*Nx+j][0]/(factIFFT)*dominfo.wallchar[i][j];
				etaCharLetaM0uv[i*Nx+j][1] = 0;
			}
		}
		fftw_plan plan_etaCharLetaM0uvhat = fftw_plan_dft_2d(Ny,Nx,etaCharLetaM0uv,etaCharLetaM0uvhat,FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_etaCharLetaM0uvhat);
		fftw_destroy_plan(plan_etaCharLetaM0uvhat);
		
		//LetaCharLetaM0uv_hat=funOprt_L2d_runup(g,Om2dsqInterp,fft2(etaCharLetaM0uv),etaCharLetaM0uv);
		fftw_complex* LetaCharLetaM0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
		fftw_complex* LetaCharLetaM0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				etaCharLetaM0uv[i*Nx+j][0] = etaCharLetaM0uv[i*Nx+j][0]*factIFFT;//for operator
			}
		}
		fun_Oprt_u(etaCharLetaM0uv,etaCharLetaM0uvhat,LetaCharLetaM0uv,LetaCharLetaM0uv_hat,Om2m1,Om2p1,Om2c1,gam_m1,gam_p1,gam_c1,Om2m2,Om2p2,Om2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
		
		//etaLetaM0uvChar     =eta.*funC_ifft2(LetaM0uvChar_hat);
		fftw_complex* etaLetaM0uvChar=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		fftw_complex* etaLetaM0uvCharhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				etaLetaM0uvChar[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*LetaM0uvChar[i*Nx+j][0]/factIFFT;
				etaLetaM0uvChar[i*Nx+j][1] = 0;
			}
		}
		fftw_plan plan_etaLetaM0uvCharhat = fftw_plan_dft_2d(Ny,Nx,etaLetaM0uvChar,etaLetaM0uvCharhat,FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_etaLetaM0uvCharhat);
		fftw_destroy_plan(plan_etaLetaM0uvCharhat);
		
		//LetaLetaCharM0uv_hat=funOprt_L2d_runup(g,Om2dsqInterp,fft2(etaLetaM0uvChar),etaLetaM0uvChar);
		fftw_complex* LetaLetaCharM0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
		fftw_complex* LetaLetaCharM0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				etaLetaM0uvChar[i*Nx+j][0] = etaLetaM0uvChar[i*Nx+j][0]*factIFFT;//for operator
			}
		}
		fun_Oprt_u(etaLetaM0uvChar,etaLetaM0uvCharhat,LetaLetaCharM0uv,LetaLetaCharM0uv_hat,Om2m1,Om2p1,Om2c1,gam_m1,gam_p1,gam_c1,
		Om2m2,Om2p2,Om2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
		
		//Chareta2M0uv_hat=fft2(Wallchar.*eta.^2.*M0uv);
        //gradChareta2M0uv_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,Chareta2M0uv_hat);
		fftw_complex* Chareta2M0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
		fftw_complex* Chareta2M0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				Chareta2M0uv[i*Nx+j][0] = dominfo.wallchar[i][j]*pow(eta[i*Nx+j][0]/factIFFT,2)*M0uv[i*Nx+j][0]/factIFFT;
				Chareta2M0uv[i*Nx+j][1] = 0;
			}
		}
		fftw_plan plan_Chareta2M0uv_hat = fftw_plan_dft_2d(Ny,Nx,Chareta2M0uv,Chareta2M0uv_hat,FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_Chareta2M0uv_hat);
		fftw_destroy_plan(plan_Chareta2M0uv_hat);

		//DivGradChareta2M0uv_hat=funOprt_div2d(dom.Kx,dom.Ky,gradChareta2M0uv_hat);
		fftw_complex* DivGradChareta2M0uv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				DivGradChareta2M0uv_hat[i*Nx+j][0]=-dominfo.kx[j]*dominfo.kx[j]*Chareta2M0uv_hat[i*Nx+j][0]-dominfo.ky[i]*dominfo.ky[i]*Chareta2M0uv_hat[i*Nx+j][0];
				DivGradChareta2M0uv_hat[i*Nx+j][1]=-dominfo.kx[j]*dominfo.kx[j]*Chareta2M0uv_hat[i*Nx+j][1]-dominfo.ky[i]*dominfo.ky[i]*Chareta2M0uv_hat[i*Nx+j][1];
			}
		}
		
        //DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,gradphi_hat);
        //DivGradphi  =funC_ifft2(DivGradphi_hat);    
		fftw_complex* DivGradphi_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		fftw_complex* DivGradphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 		   
	
		fun_div(DivGradphi_hat,uhat,vhat,1.0,Nx,Ny);
		
		fftw_plan plan_DivGradphi = fftw_plan_dft_2d(Ny,Nx,DivGradphi_hat,DivGradphi,FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan_DivGradphi);
		fftw_destroy_plan(plan_DivGradphi);	
	
		//eta2CharDivGradphi=eta.^2.*HeavChar.*DivGradphi;        
		fftw_complex* eta2CharDivGradphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		fftw_complex* eta2CharDivGradphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny)); 
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				eta2CharDivGradphi[i*Nx+j][0] = pow(eta[i*Nx+j][0]/factIFFT,2)*dominfo.wallchar[i][j]*DivGradphi[i*Nx+j][0]/factIFFT;
				eta2CharDivGradphi[i*Nx+j][1] = 0;
			}
		}
		fftw_plan plan_eta2CharDivGradphihat = fftw_plan_dft_2d(Ny,Nx,eta2CharDivGradphi,eta2CharDivGradphihat,FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_eta2CharDivGradphihat);
		fftw_destroy_plan(plan_eta2CharDivGradphihat);
		
		//Leta2CharDivGradphi_hat=funOprt_L2d_runup(g,Om2dsqInterp,fft2(eta2CharDivGradphi),eta2CharDivGradphi);
		fftw_complex* Leta2CharDivGradphi_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));     
		fftw_complex* Leta2CharDivGradphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				eta2CharDivGradphi[i*Nx+j][0] = eta2CharDivGradphi[i*Nx+j][0]*factIFFT;//for operator
			}
		}
		fun_Oprt_u(eta2CharDivGradphi,eta2CharDivGradphihat,Leta2CharDivGradphi,Leta2CharDivGradphi_hat,Om2m1,Om2p1,Om2c1,
		gam_m1,gam_p1,gam_c1,Om2m2,Om2p2,Om2c2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
		
		//delphiH3_hat  =0.5*(LetaCharLetaM0uv_hat+LetaLetaCharM0uv_hat+Leta2CharDivGradphi_hat+DivGradChareta2M0uv_hat);
        //deletaH3_hat  =fft2(M0uv.*(eta.*DivGradphi.*HeavChar+0.5*funC_ifft2(LetaM0uvChar_hat)+0.5*HeavChar.*LetaM0uv));
		fftw_complex* deletaH3=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
		fftw_plan plan_deletaH3_hat = fftw_plan_dft_2d(Ny,Nx,deletaH3,deletaH3_hat,FFTW_FORWARD, FFTW_ESTIMATE);
		//#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				delphiH3_hat[i*Nx+j][0] = 0.5*(LetaLetaCharM0uv_hat[i*Nx+j][0]/grav+DivGradChareta2M0uv_hat[i*Nx+j][0]+LetaCharLetaM0uv_hat[i*Nx+j][0]/grav
												+Leta2CharDivGradphi_hat[i*Nx+j][0]/grav);//term ini blow up
				delphiH3_hat[i*Nx+j][1] = 0.5*(LetaLetaCharM0uv_hat[i*Nx+j][1]/grav+DivGradChareta2M0uv_hat[i*Nx+j][1]+LetaCharLetaM0uv_hat[i*Nx+j][1]/grav
												+Leta2CharDivGradphi_hat[i*Nx+j][1]/grav);//term ini blow up
		
				deletaH3[i*Nx+j][0] = M0uv[i*Nx+j][0]/factIFFT*(eta[i*Nx+j][0]/factIFFT*DivGradphi[i*Nx+j][0]/factIFFT*dominfo.wallchar[i][j]
										+0.5*LetaM0uvChar[i*Nx+j][0]/factIFFT+0.5*dominfo.wallchar[i][j]*LetaM0uv[i*Nx+j][0]/factIFFT);
				deletaH3[i*Nx+j][1] = 0;
			}
		}
		
		fftw_execute(plan_deletaH3_hat);
		fftw_destroy_plan(plan_deletaH3_hat);		
		
		fftw_free(deletaH3);
		
		fftw_free(eta2CharDivGradphi);
		fftw_free(eta2CharDivGradphihat);
		fftw_free(Leta2CharDivGradphi_hat);
		fftw_free(Leta2CharDivGradphi);
		
		fftw_free(DivGradChareta2M0uv_hat);
		fftw_free(DivGradphi_hat);
		fftw_free(DivGradphi);
		
		fftw_free(Chareta2M0uv);
		fftw_free(Chareta2M0uv_hat);
		
		fftw_free(etaLetaM0uvChar);
		fftw_free(etaLetaM0uvCharhat);
		fftw_free(LetaLetaCharM0uv_hat);
		fftw_free(LetaLetaCharM0uv);
		
        fftw_free(LetaCharLetaM0uv_hat);
        fftw_free(LetaCharLetaM0uv);
        fftw_free(etaCharLetaM0uv);
        fftw_free(etaCharLetaM0uvhat);
        fftw_free(etaM0uv);
        fftw_free(etaM0uv_hat);
        fftw_free(LetaM0uv_hat);
        fftw_free(LetaM0uv);        
	}	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH23_hat[i*Nx+j][0] = delphiH2_hat[i*Nx+j][0]+delphiH3_hat[i*Nx+j][0];
			delphiH23_hat[i*Nx+j][1] = delphiH2_hat[i*Nx+j][1]+delphiH3_hat[i*Nx+j][1];
			ikxdeletaH23_hat[i*Nx+j][0] = -dominfo.kx[j]*(deletaH2_hat[i*Nx+j][1]+deletaH3_hat[i*Nx+j][1]);
			ikxdeletaH23_hat[i*Nx+j][1] = dominfo.kx[j]*(deletaH2_hat[i*Nx+j][0]+deletaH3_hat[i*Nx+j][0]);
			ikydeletaH23_hat[i*Nx+j][0] = -dominfo.ky[i]*(deletaH2_hat[i*Nx+j][1]+deletaH3_hat[i*Nx+j][1]);
			ikydeletaH23_hat[i*Nx+j][1] = dominfo.ky[i]*(deletaH2_hat[i*Nx+j][0]+deletaH3_hat[i*Nx+j][0]);
		}
	}
	fun_adjust(delphiH23adj_hat,delphiH23_hat,Nx,Ny);
	fun_adjust(ikxdeletaH23adj_hat,ikxdeletaH23_hat,Nx,Ny);
	fun_adjust(ikydeletaH23adj_hat,ikydeletaH23_hat,Nx,Ny);

	//for outputting to RHS
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH_hat[i*Nx+j][0] = delphiH1_hat[i*Nx+j][0]+delphiH23adj_hat[i*Nx+j][0];
			delphiH_hat[i*Nx+j][1] = delphiH1_hat[i*Nx+j][1]+delphiH23adj_hat[i*Nx+j][1];
		}
	}

	/*if (t>30){
		double* forline=new double[Ny];
		for (int i=0;i<Ny;i++){
			forline[i]=delphiH2_hat[i*Nx][0];
		}
		fun_plot_line(dominfo.y,forline,Ny);
		for (int i=0;i<Ny;i++){
			forline[i]=delphiH2_hat[i*Nx][1];
		}
		fun_plot_line(dominfo.y,forline,Ny);
		for (int i=0;i<Ny;i++){
			forline[i]=deletaH2[i*Nx][0];
		}
		fun_plot_line(dominfo.y,forline,Ny);
		delete[] forline;
	}*/
	
	fftw_free(delphiH23adj_hat);
	fftw_free(delphiH23_hat);
	fftw_free(ikxdeletaH23_hat);
	fftw_free(ikydeletaH23_hat);
	
	fftw_free(delphiH3_hat);
	fftw_free(deletaH3_hat);
	
	fftw_free(delphiH2_hat);
	fftw_free(deletaH2_hat);
	fftw_free(deletaH2);	
	fftw_free(divetagradphi_hat);

	fftw_free(delphiH1_hat);
	
	fftw_free(etauhat);
	fftw_free(etau);
	fftw_free(etavhat);
	fftw_free(etav);	
	
	fftw_free(LetaM0uvChar_hat);
	fftw_free(LetaM0uvChar);
	fftw_free(etaM0uvChar);
	fftw_free(etaM0uvChar_hat);
	
	fftw_free(M0uv);
	fftw_free(M0uv_hat);
	
	fftw_free(Cuhat);
	fftw_free(Cvhat);
	fftw_free(Cu);
	fftw_free(Cv);
	fftw_free(CuCharhat);
	fftw_free(CvCharhat);
	fftw_free(CuChar);
	fftw_free(CvChar);

	fftw_free(C2u);
	fftw_free(C2v);
	fftw_free(C2uhat);
	fftw_free(C2vhat);
	fftw_free(C2uCharhat);
	fftw_free(C2vCharhat);
	fftw_free(C2uChar);
	fftw_free(C2vChar);
	
	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(etahat);
}

void fun_delphiH1hat(double t,int Nx,int Ny,const double* Zhat,fftw_complex* delphiH1hat)
{
	fftw_complex *deluH1hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH1hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	if (strcmp(wall,"no")!=0){//wall exist
		if ((refl_coef==1)){//&&(strcmp(evolinfo.model,"HS1")==0)){
			complex **C2u_xhat = declare_2darray_complex(Nx,Ny);
			complex **C2u_yhat = declare_2darray_complex(Nx,Ny);
			double **Flux_x 	  = declare_2darray(Nx,Ny);
			double **Flux_y 	  = declare_2darray(Nx,Ny);
			complex **divC2u_hat  = declare_2darray_complex(Nx,Ny);
			double **divC2u 	  = declare_2darray(Nx,Ny);
			double **divC2u_wall  = declare_2darray(Nx,Ny);
			double **Flux_grad    = declare_2darray(Nx,Ny);
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					C2u_xhat[i][j].Re = pow(Oprt.C2d[i][j],2)*Zhat[2*Nx*Ny+i*Nx+j]/grav;
					C2u_xhat[i][j].Im = pow(Oprt.C2d[i][j],2)*Zhat[3*Nx*Ny+i*Nx+j]/grav;
					C2u_yhat[i][j].Re = pow(Oprt.C2d[i][j],2)*Zhat[4*Nx*Ny+i*Nx+j]/grav;
					C2u_yhat[i][j].Im = pow(Oprt.C2d[i][j],2)*Zhat[5*Nx*Ny+i*Nx+j]/grav;
				}
			}
			ifftw_2d_c2r(Flux_x,C2u_xhat,Nx,Ny);
			ifftw_2d_c2r(Flux_y,C2u_yhat,Nx,Ny);
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					divC2u_hat[i][j].Re = -(-dominfo.kx[j]*C2u_xhat[i][j].Im-dominfo.ky[i]*C2u_yhat[i][j].Im);
					divC2u_hat[i][j].Im = -(dominfo.kx[j]*C2u_xhat[i][j].Re+dominfo.ky[i]*C2u_yhat[i][j].Re);
				}
			}
			ifftw_2d_c2r(divC2u,divC2u_hat,Nx,Ny);			
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					divC2u_wall[i][j] = divC2u[i][j]*dominfo.wallchar[i][j];
					Flux_grad[i][j]   = Flux_x[i][j]*grad_wallchar_x[i][j]+Flux_y[i][j]*grad_wallchar_y[i][j];
				}
			}
			complex** divC2u_wall_hat = declare_2darray_complex(Nx,Ny);
			fftw_2d_r2c(divC2u_wall_hat,divC2u_wall,Nx,Ny);
			complex **Flux_grad_hat	  = declare_2darray_complex(Nx,Ny);
			fftw_2d_r2c(Flux_grad_hat,Flux_grad,Nx,Ny);	
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					delphiH1hat[i*Nx+j][0] = divC2u_wall_hat[i][j].Re-Flux_grad_hat[i][j].Re;
					delphiH1hat[i*Nx+j][1] = divC2u_wall_hat[i][j].Im-Flux_grad_hat[i][j].Im;				
				}
			}
			
			free_2darray_complex(Flux_grad_hat,Nx,Ny);
			free_2darray_complex(divC2u_wall_hat,Nx,Ny);
			free_2darray_complex(C2u_xhat,Nx,Ny);
			free_2darray_complex(C2u_yhat,Nx,Ny);
			free_2darray(Flux_x,Nx,Ny);
			free_2darray(Flux_y,Nx,Ny);	
			free_2darray_complex(divC2u_hat,Nx,Ny);	
			free_2darray(divC2u,Nx,Ny);
			free_2darray(divC2u_wall,Nx,Ny);
			free_2darray(Flux_grad,Nx,Ny);	
		}
		else {
			complex **Cu_xhat = declare_2darray_complex(Nx,Ny);
			complex **Cu_yhat = declare_2darray_complex(Nx,Ny);
			double **Cu_x 	  = declare_2darray(Nx,Ny);
			double **Cu_y 	  = declare_2darray(Nx,Ny);
			double **Cu_xWall = declare_2darray(Nx,Ny);
			double **Cu_yWall = declare_2darray(Nx,Ny);
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					Cu_xhat[i][j].Re = (Oprt.C2d[i][j])*Zhat[2*Nx*Ny+i*Nx+j];
					Cu_xhat[i][j].Im = (Oprt.C2d[i][j])*Zhat[3*Nx*Ny+i*Nx+j];
					Cu_yhat[i][j].Re = (Oprt.C2d[i][j])*Zhat[4*Nx*Ny+i*Nx+j];
					Cu_yhat[i][j].Im = (Oprt.C2d[i][j])*Zhat[5*Nx*Ny+i*Nx+j];
				}
			}
			ifftw_2d_c2r(Cu_x,Cu_xhat,Nx,Ny);
			ifftw_2d_c2r(Cu_y,Cu_yhat,Nx,Ny);
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					Cu_xWall[i][j] = Cu_x[i][j]*dominfo.wallchar[i][j];
					Cu_yWall[i][j] = Cu_y[i][j]*dominfo.wallchar[i][j];
				}
			}
			complex **Cu_xhat_w = declare_2darray_complex(Nx,Ny);
			complex **Cu_yhat_w = declare_2darray_complex(Nx,Ny);
			fftw_2d_r2c(Cu_xhat_w,Cu_xWall,Nx,Ny);
			fftw_2d_r2c(Cu_yhat_w,Cu_yWall,Nx,Ny);
			
			#pragma omp parallel for
			for (int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					deluH1hat[i*Nx+j][0] = (Oprt.C2d[i][j])/grav*Cu_xhat_w[i][j].Re;				
					deluH1hat[i*Nx+j][1] = (Oprt.C2d[i][j])/grav*Cu_xhat_w[i][j].Im;						
					
					delvH1hat[i*Nx+j][0] = (Oprt.C2d[i][j])/grav*Cu_yhat_w[i][j].Re;				
					delvH1hat[i*Nx+j][1] = (Oprt.C2d[i][j])/grav*Cu_yhat_w[i][j].Im;
					
					delphiH1hat[i*Nx+j][0] = dominfo.kx[j]*deluH1hat[i*Nx+j][1]+dominfo.ky[i]*delvH1hat[i*Nx+j][1];				
					delphiH1hat[i*Nx+j][1] = -dominfo.kx[j]*deluH1hat[i*Nx+j][0]-dominfo.ky[i]*delvH1hat[i*Nx+j][0];	
				}
			}
			
			free_2darray_complex(Cu_xhat,Nx,Ny);
			free_2darray_complex(Cu_yhat,Nx,Ny);
			free_2darray_complex(Cu_xhat_w,Nx,Ny);
			free_2darray_complex(Cu_yhat_w,Nx,Ny);
			free_2darray(Cu_x,Nx,Ny);
			free_2darray(Cu_y,Nx,Ny);
			free_2darray(Cu_xWall,Nx,Ny);
			free_2darray(Cu_yWall,Nx,Ny);	
		}
	}
	else {
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				deluH1hat[i*Nx+j][0] = pow(Oprt.C2d[i][j],2)/grav*Zhat[2*Nx*Ny+i*Nx+j];				
				deluH1hat[i*Nx+j][1] = pow(Oprt.C2d[i][j],2)/grav*Zhat[3*Nx*Ny+i*Nx+j];
				
				delvH1hat[i*Nx+j][0] = pow(Oprt.C2d[i][j],2)/grav*Zhat[4*Nx*Ny+i*Nx+j];
				delvH1hat[i*Nx+j][1] = pow(Oprt.C2d[i][j],2)/grav*Zhat[5*Nx*Ny+i*Nx+j];
				
				delphiH1hat[i*Nx+j][0] = dominfo.kx[j]*deluH1hat[i*Nx+j][1]+dominfo.ky[i]*delvH1hat[i*Nx+j][1];				
				delphiH1hat[i*Nx+j][1] = -dominfo.kx[j]*deluH1hat[i*Nx+j][0]-dominfo.ky[i]*delvH1hat[i*Nx+j][0];
			}
		}
	}
	fftw_free(deluH1hat);
	fftw_free(delvH1hat);
}

/*void fun_delphiH1hat_runup(fftw_complex* uh,fftw_complex* vh,fftw_complex* delphiH1hat,domvar dominfo,const double* Zhat)
{
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* u_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* v_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
    
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){	 
		u_hat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
		u_hat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
		v_hat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
		v_hat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];		
	  }
	}

	//interpolation 1
	fftw_complex* Cu = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Cv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* C_uhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* C_vhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fun_Oprt_u(uh,u_hat,Cu,C_uhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	fun_Oprt_u(vh,v_hat,Cv,C_vhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);

	//if (t>0){
		double** forplot=declare_2darray(Nx,Ny);
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				forplot[i][j]=Cv[i*Nx+j][0];
			}
		}
		FILE* gph=popen("gnuplot -persistent","w");
		fun_plot_2d_trarray(gph,dominfo.y,dominfo.x,forplot,Ny,Nx);free_2darray(forplot,Nx,Ny);
		fflush(gph);
		pclose(gph);getchar();
	//}
    
    //interpolation 2
    fftw_complex* Csqr_u = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Csqr_v = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* Csqr_uhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Csqr_vhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fun_Oprt_u(Cu,C_uhat,Csqr_u,Csqr_uhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	fun_Oprt_u(Cv,C_vhat,Csqr_v,Csqr_vhat,Cm1,Cp1,Cc1,gam_m1,gam_p1,gam_c1,Cm2,Cp2,Cc2,gam_m2,gam_p2,gam_c2,Nx,Ny,factIFFT);
	
	fun_div(delphiH1hat,Csqr_uhat,Csqr_vhat,-1.0/grav);

	fftw_free(Cu);fftw_free(Cv);fftw_free(C_uhat);fftw_free(C_vhat);
	fftw_free(Csqr_uhat);fftw_free(Csqr_vhat);fftw_free(Csqr_u);fftw_free(Csqr_v);
	fftw_free(u_hat);fftw_free(v_hat);
}*/

int rhs_shore2D_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{
	
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	fftw_complex* etaprev = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etafric = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
		
	//Getting eta not normalized yet by factIFFT
	fftw_complex* eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* v   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	get_eta(eta,t,Zhat,Nx,Ny);
	get_u(u,t,Zhat,Nx,Ny);
	get_v(v,t,Zhat,Nx,Ny);

	//Define active domain by a characteristic function
	double** geta=declare_2darray(Nx,Ny);
	complex** getahat=declare_2darray_complex(Nx,Ny);
	double** WaveChar  = declare_2darray(Nx,Ny);
	double** dampcharH = declare_2darray(Nx,Ny);
	double** 	H      = declare_2darray(Nx,Ny);
	set_matrix_val(WaveChar,Nx,Ny,1);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){//eta no scaling
			if ((eta[i*Nx+j][0]/factIFFT-(dominfo.bathy.min[i][j]+H_minShore))<0){
				eta[i*Nx+j][0] = (dominfo.bathy.min[i][j]+H_minShore)*factIFFT;
			}
			etaprev[i*Nx+j][0]= eta[i*Nx+j][0];
			dampcharH[i][j]   = 1-cfSA[i][j];
			geta[i][j]        = grav*(eta[i*Nx+j][0]/factIFFT+dominfo.bathy.plus[i][j]);
			H[i][j]  		  = ChiAdj[i][j]*eta[i*Nx+j][0]/factIFFT-dominfo.bathy.min[i][j];
			etafric[i*Nx+j][0]= ChiAdj[i][j]*eta[i*Nx+j][0];
			if (H[i][j]<H_minShore){
				WaveChar[i][j] = 0;
			}			
		}
	}
	fftw_2d_r2c(getahat,geta,Nx,Ny);

	//get friction term
    fun_sfriction(t,Nx,Ny,factIFFT,etafric,u,v);
	
	//damping and source
	fftw_complex *dampetahat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fun_source(t,Nx,Ny);
	fun_damping(dampcharH,etaprev,u,v,dampetahat,dampuhat,dampvhat);
	
	//start dynamic
	fftw_complex* ikxdeletaH23adj_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* ikydeletaH23adj_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex *delphiH_hat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	
	fun_termshat_shore_uv(t,eta,u,v,delphiH_hat,ikxdeletaH23adj_hat,ikydeletaH23adj_hat,dominfo,Zhat,WaveChar,H);//eta changed here

	//get B for breaking
	if (strcmp(breaking,"yes")==0){
		complex** delphiHhat=declare_2darray_complex(Nx,Ny);
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				delphiHhat[i][j].Re=delphiH_hat[i*Nx+j][0];
				delphiHhat[i][j].Im=delphiH_hat[i*Nx+j][1];
			}
		}
		set_matrix_val(B,Nx,Ny,0.0);
		breaking_process(t,Zhat,etaprev,u,v,factIFFT);// filling the global matrix B
		
		fun_sbreaking(t,eta,delphiHhat,Nx,Ny);
		free_2darray_complex(delphiHhat,Nx,Ny);
	}

	//time stepping ODE
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){		  
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re
										+delphiH_hat[i*Nx+j][0]
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im
										+delphiH_hat[i*Nx+j][1]
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (dominfo.kx[j]*getahat[i][j].Im-ikxdeletaH23adj_hat[i*Nx+j][0]+Sf_xhat[i][j].Re+Sb_xhat[i][j].Re
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-dominfo.kx[j]*getahat[i][j].Re-ikxdeletaH23adj_hat[i*Nx+j][1]+Sf_xhat[i][j].Im+Sb_xhat[i][j].Im
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (dominfo.ky[i]*getahat[i][j].Im-ikydeletaH23adj_hat[i*Nx+j][0]+Sf_yhat[i][j].Re+Sb_yhat[i][j].Re
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-dominfo.ky[i]*getahat[i][j].Re-ikydeletaH23adj_hat[i*Nx+j][1]+Sf_yhat[i][j].Im+Sb_yhat[i][j].Im
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];
		}
	}
	
	fftw_free(dampetahat);
	fftw_free(dampuhat);
	fftw_free(dampvhat);
	fftw_free(delphiH_hat);
	fftw_free(ikxdeletaH23adj_hat);
	fftw_free(ikydeletaH23adj_hat);
		
	//delete[] indexzero;
	free_2darray(geta,Nx,Ny);
	free_2darray_complex(getahat,Nx,Ny);
	fftw_free(eta);
	fftw_free(u);
	fftw_free(v);
	fftw_free(etaprev);fftw_free(etafric);
	
	free_2darray(WaveChar,Nx,Ny);
	free_2darray(dampcharH,Nx,Ny);
	free_2darray(H,Nx,Ny);

	return GSL_SUCCESS;
}

int rhs_shore2D_uv_HSdirect(double t, const double* Zhat, double* dt_Zhat,void *params)
{
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;

	fftw_complex* etaH = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	//Getting eta not normalized yet by factIFFT
	fftw_complex* eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* v   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	get_eta(eta,t,Zhat,Nx,Ny);
	get_u(u,t,Zhat,Nx,Ny);
	get_v(v,t,Zhat,Nx,Ny);
	
	
	//Define active domain by a characteristic function
	double** WaveChar  = declare_2darray(Nx,Ny);
	double** dampcharH = declare_2darray(Nx,Ny);
	double** 	H      = declare_2darray(Nx,Ny);
	set_matrix_val(WaveChar,Nx,Ny,1);

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			H[i][j] = ChiAdj[i][j]*eta[i*Nx+j][0]-dominfo.bathy.data[i][j];
			if (H[i][j]<H_minShore){
				WaveChar[i][j] = 0;
			}
		}
	}

	int* indexzero = new int[Ny];
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (WaveChar[i][j] == 0){
				indexzero[i] = j; break;
			}
		}
	}

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dampcharH[i][j]= dominfo.fbdy.charac[i][j];
			etaH[i][j]     = eta[i*Nx+j][0];
		}
	}

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=indexzero[i];j<Nx;j++){
			WaveChar[i][j] = 0;
			dampcharH[i][j]= 1;
			etaH[i][j]     = H[i][j];
		}
	}

	//dynamic
	fftw_complex *deletaH1_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex *deletaH2_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *deluH1_hat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *deluH2_hat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *delvH1_hat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH2_hat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *dampetahat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fun_termshat_shore_uv_HSdirect(t,eta,u,v,deluH1_hat,delvH1_hat,deletaH1_hat,deluH2_hat,delvH2_hat,deletaH2_hat,dominfo,Zhat,WaveChar,H);
	fun_source(t,Nx,Ny);
	fun_damping(dampcharH,etaH,u,v,dampetahat,dampuhat,dampvhat);
	
	//time stepping
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){					  
			dt_Zhat[i*Nx+j]  	   	= (deluH1_hat[i*Nx+j][0]+delvH1_hat[i*Nx+j][0]
										+deluH2_hat[i*Nx+j][0]+delvH2_hat[i*Nx+j][0]
										+Source_hat[i][j].Re-7*dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (deluH1_hat[i*Nx+j][1]+delvH1_hat[i*Nx+j][1]
										+deluH2_hat[i*Nx+j][1]+delvH2_hat[i*Nx+j][1]
										+Source_hat[i][j].Im-7*dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= ((-dominfo.kx[j]*deletaH1_hat[i*Nx+j][1])+Sf_xhat[i][j].Re+Sb_xhat[i][j].Re
										-7*dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= ((dominfo.kx[j]*deletaH1_hat[i*Nx+j][0])+Sf_xhat[i][j].Im+Sb_xhat[i][j].Im
										-7*dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= ((-dominfo.ky[i]*deletaH1_hat[i*Nx+j][1])+Sf_yhat[i][j].Re+Sb_yhat[i][j].Re
										-7*dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= ((dominfo.ky[i]*deletaH1_hat[i*Nx+j][0])+Sf_yhat[i][j].Im+Sb_yhat[i][j].Im
										-7*dampvhat[i*Nx+j][1])*dominfo.aal[i][j];
		}
	}	
	
	free_2darray(H,Nx,Ny);
	free_2darray(WaveChar,Nx,Ny);
	free_2darray(dampcharH,Nx,Ny);
	fftw_free(dampetahat);
	fftw_free(dampuhat);fftw_free(dampvhat);
	fftw_free(deluH1_hat);fftw_free(delvH1_hat);
	fftw_free(deletaH1_hat);
	fftw_free(deluH2_hat);fftw_free(delvH2_hat);
	fftw_free(deletaH2_hat);
	fftw_free(eta);fftw_free(u);fftw_free(v);
	fftw_free(etaH);
	delete[] indexzero;
	
	return GSL_SUCCESS;	
}

void print_xt(FILE* fpc,int ntpart,double* th,int Nx,int Ny)
{
	double NX = (double) Nx;
	double NT = (double) ntpart;
	double NY = (double) Ny;
	
	fwrite(&NX,sizeof(double),1,fpc);
	fwrite(dominfo.x,sizeof(double),Nx,fpc);
	fwrite(&NY,sizeof(double),1,fpc);
	fwrite(dominfo.y,sizeof(double),Ny,fpc);
	fwrite(&NT,sizeof(double),1,fpc);
	fwrite(th,sizeof(double),ntpart,fpc);
	fwrite(&seainfo.Hs,sizeof(double),1,fpc);
	fwrite(&seainfo.kp,sizeof(double),1,fpc);
	fwrite(&seainfo.wp,sizeof(double),1,fpc);	
	fwrite(dominfo.kx,sizeof(double),Nx,fpc);
	fwrite(dominfo.ky,sizeof(double),Ny,fpc);
	
	for (int kk=0;kk<Ny;kk++){
		fwrite(dominfo.bathy.data[kk],sizeof(double),Nx,fpc);
	}

	for (int kk=0;kk<Ny;kk++){
		fwrite(cfSA[kk],sizeof(double),Nx,fpc);
	}
	
	if (strcmp(wavename,"zero")!=0){
		for (int kk=0;kk<Ny;kk++){
			fwrite(ChiAdjPhi[kk],sizeof(double),Nx,fpc);
		}
	}
}

void ode_solver_uv(gsl_odeiv2_system sys,double* y_hat,double* tinterv,int N,int Nx,int Ny,int Nt,double eps_abs,double eps_rel)
{   	
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
	gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, N);
	gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new(eps_abs,eps_rel);
	gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (N);
	
	Source_hat 	= declare_2darray_complex(Nx,Ny);
	set_matrix_complex_val(Source_hat,Nx,Ny,0.0);
	
	float  float_eta,float_u,float_v;
	double t = tinterv[0]; 
	double dt= (tinterv[1]-tinterv[0]);
	double h = dt/5.0;   //starting step size;
	double ti,max_etahere;
	int    status;
	int ipart = 1;
	ntpart= (int) dominfo.Nt/npartition; 
	
	double*   tim=new double[Nt];//time for output
	complex** etahat=declare_2darray_complex(Nx,Ny);
	complex** uhat=declare_2darray_complex(Nx,Ny);
	complex** vhat=declare_2darray_complex(Nx,Ny);
	complex** dtetahat=declare_2darray_complex(Nx,Ny);
	double**  eta =declare_2darray(Nx,Ny);
	double**  u	=declare_2darray(Nx,Ny);
	double**  v	=declare_2darray(Nx,Ny);
	double**  dteta	=declare_2darray(Nx,Ny);
	for (int kk=0;kk<Ny;kk++){
		for(int ll=0;ll<Nx;ll++){
			eta[kk][ll] = waveinit.profile[kk][ll];
			u[kk][ll] = waveinit.u[kk][ll];
			v[kk][ll] = waveinit.v[kk][ll];
		}
	}
	
	char str_gauge1[256];
	char str_gauge2[256];
	char str_gauge3[256];	
	sprintf(str_gauge1, "%sHawassi_eta_gauge.dat",arg2);
	sprintf(str_gauge2, "%sHawassi_u_gauge.dat",arg2);
	sprintf(str_gauge3, "%sHawassi_v_gauge.dat",arg2);
	FILE* fgauge1;FILE* feta;
	FILE* fgauge2;FILE* fu;
	FILE* fgauge3;FILE* fv;
	if (ngauge > 0){
		fgauge1=fopen(str_gauge1,"wb");
		fgauge2=fopen(str_gauge2,"wb");
		fgauge3=fopen(str_gauge3,"wb");
	}
	
	double N1, N2x, N2y, N3, N4, N5;	
	int  npart = Nt/100.0;
	int  ip = 0;
	char str[20];
	int i25 = int(Nt/4);
	int i50 = int(Nt/2);
	int i75 = int(3*Nt/4);
	gsl_ieee_env_setup();
    for (int i=0;i<Nt;i++){
		if (i==i25){
			printf("25%%\n");
		}
		else if (i==i50){
			printf("50%%\n");
		}
		else if (i==i75){
			printf("75%%\n");
		}
		
		if (((i % npart)==0)&&(ip<100)){
			//printf("%d%%\n",ip);
			ip=ip+1;
		}
		ti=tinterv[i];
						
		while (t < ti){
			status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, ti, &h, y_hat);
			if (status != GSL_SUCCESS){
				break;
			}	
			tprev=t;		
		}
		tim[i]=t;
		
		#pragma omp parallel for
		for (int kk=0;kk<Ny;kk++){
			for(int ll=0;ll<Nx;ll++){
				etahat[kk][ll].Re=y_hat[kk*Nx+ll];
				etahat[kk][ll].Im=y_hat[Nx*Ny+kk*Nx+ll];
				uhat[kk][ll].Re=y_hat[2*Nx*Ny+kk*Nx+ll];
				uhat[kk][ll].Im=y_hat[3*Nx*Ny+kk*Nx+ll];
				vhat[kk][ll].Re=y_hat[4*Nx*Ny+kk*Nx+ll];
				vhat[kk][ll].Im=y_hat[5*Nx*Ny+kk*Nx+ll];
				dtetahat[kk][ll].Re=dy_hat[kk*Nx+ll];
				dtetahat[kk][ll].Im=dy_hat[Nx*Ny+kk*Nx+ll];
			}  
		}		
		ifftw_2d_c2r(eta,etahat,Nx,Ny);
		ifftw_2d_c2r(u,uhat,Nx,Ny);
		ifftw_2d_c2r(v,vhat,Nx,Ny);
		ifftw_2d_c2r(dteta,dtetahat,Nx,Ny);

		//update initial condition for runup case 
		if (runup_id==1){//runup yes			
			#pragma omp parallel for
			for (int kk=0;kk<Ny;kk++){
				for(int ll=0;ll<Nx;ll++){
					eta[kk][ll] = fmax(eta[kk][ll],dominfo.bathy.min[kk][ll]+H_minShore)*cfSA[kk][ll];					
				}
			}			
			fftw_2d_r2c(etahat,eta,Nx,Ny);
			#pragma omp parallel for
			for (int kk=0;kk<Ny;kk++){
				for(int ll=0;ll<Nx;ll++){
					y_hat[kk*Nx+ll]      =etahat[kk][ll].Re;
					y_hat[Nx*Ny+kk*Nx+ll]=etahat[kk][ll].Im;
				}
			}			
		}
		
		//opening file output 
		if ((i % ntpart)==0){
			char* str_eta=new char[256];
			char* str_u=new char[256];
			char* str_v=new char[256];	
			sprintf(str_eta, "%sHawassi_eta_%02d.dat",arg2,ipart);
			sprintf(str_u, "%sHawassi_u_%02d.dat",arg2,ipart);
			sprintf(str_v, "%sHawassi_v_%02d.dat",arg2,ipart);	
			feta=fopen(str_eta,"wb");
			fu=fopen(str_u,"wb");
			fv=fopen(str_v,"wb");
			
			//filling constants
			if (i==0){
				char* str_con=new char[256];
				sprintf(str_con, "%sConstant_%02d.dat",arg2,ipart);
				double* th = new double[ntpart];
				for (int j=0;j<ntpart;j++){
					th[j] = (i+j)*dominfo.dt;
				}
				FILE* fcon=fopen(str_con,"wb");
				print_xt(fcon,ntpart,th,Nx,Ny);
				fclose(fcon);
				delete[] th;
				delete[] str_con;
			}
			ipart++;
			delete[] str_eta;
			delete[] str_u;
			delete[] str_v;
		}	

		//writing output whole domain
		if (runup_id==1){//runup yes			
			#pragma omp parallel for
			for (int kk=0;kk<Ny;kk++){
				for(int ll=0;ll<Nx;ll++){
					eta[kk][ll] = eta[kk][ll]+dominfo.bathy.plus[kk][ll];					
				}
			}		
		}

		for (int kk=0;kk<Ny;kk++){
			for(int ll=0;ll<Nx;ll++){
				float_eta = float(eta[kk][ll]);
				fwrite(&float_eta,sizeof(float),1,feta);  
			}
		}
		fflush(feta);

		if ((strcmp(kinematic,"no")==0) && ((i % ntpart)==0)){
			for (int kk=0;kk<Ny;kk++){
				for(int ll=0;ll<Nx;ll++){
					float_u = float(u[kk][ll]);
					float_v = float(v[kk][ll]);
					fwrite(&float_u,sizeof(float),1,fu);
					fwrite(&float_v,sizeof(float),1,fv);  
				}
			}
				
			fflush(fu);
			fflush(fv);
		}
		
		//writing output at gauges
		if (ngauge>0){
			print_data_gauges_uv(i,t,eta,u,v,fgauge1,fgauge2,fgauge3);
		}		
		
		//IF kinematic yes
		if (strcmp(kinematic,"yes")==0){
			for (int kk=0;kk<Ny;kk++){
				for(int ll=0;ll<Nx;ll++){
					float_u = float(u[kk][ll]);
					float_v = float(v[kk][ll]);
					fwrite(&float_u,sizeof(float),1,fu);
					fwrite(&float_v,sizeof(float),1,fv);  
				}
			}
				
			fflush(fu);
			fflush(fv);
			
			if (strcmp(kinematic_type,"boundary")!=0){
				if ((i % ntpart)==0){
					sprintf(str, "%02d.dat", ipart);
					double Ntout;
					if (ipart<=npartition){
						Ntout = (double) ntpart;
					}
					else{
						Ntout = (double) Ntrest;
					}				
					N1  = (double) Ntout;
					N2x = (double) NxoutIP;
					N2y = (double) NyoutIP;
					N3  = (double) NzoutIP;
					N4  = N1*N2x*N2y;
					N5  = N1*N2x*N2y*N3;						
					
					//save header data for each partition
					save_P(var_P,str,N1,N2x,N2y,N3,N4,N5,ipart);			
					save_elev(var_elev,str,N1,N2x,N2y,N3,N4,N5,ipart);			
					save_vel(var_vel,str,N1,N2x,N2y,N3,N4,N5,ipart);								
					save_acc(var_acc,str,N1,N2x,N2y,N3,N4,N5,ipart);
					save_phi(var_phi,str);
					
					fflush(dateta);
					fflush(datP);
					fflush(datvel);
					fflush(datacc);
					fflush(datPhi);
				}
			}
			
			kinematic_modul(i,t,eta,u,v,dteta,dominfo.bathy.pos,Nx,Ny,dominfo.KK,seainfo.wp); 
			
			if (strcmp(kinematic_type,"boundary")!=0){		
				fflush(dateta);
				fflush(datP);
				fflush(datvel);
				fflush(datacc);
				fflush(datPhi);
			}
		}		
		
		//filling coupling term
		if (strcmp(coupling,"yes")==0){
			iter_t = i;
			HAWASSI_coupling_term_uv(dir_force,ti,dominfo.dt,eta,u,v,ntpart,Nx,Ny);
		}
		
		//time step failure checking
		max_etahere=fun_max2darray(eta,Nx,Ny);
		if (isnan(max_etahere)){
			printf("Time stepping failed at t=%fs.\n",t);
			exit(0);
		}
		
		//plot at halfway
		if (i==(1) || i==(dominfo.Nt/2.0)){
			FILE* gp=popen("gnuplot ","w");
			fprintf(gp,"set term png\n");
			fprintf(gp,"set output 'Eta.png'\n");
			fun_plot_2d_trarray(gp,dominfo.y,dominfo.x,eta,dominfo.Ny,dominfo.Nx);
			fflush(gp);
			pclose(gp);
		}	
	}
	//end of time integration
	
	fclose(feta);
	fclose(fu);
	fclose(fv);
	if (ngauge>0){
		fclose(fgauge1);
		fclose(fgauge2);
		fclose(fgauge3);
	}
	if (strcmp(breaking,"yes")==0){
		fclose(fbreak);
	}

	if ((strcmp(kinematic,"yes")==0)&&(strcmp(kinematic_type,"boundary")!=0)){
		fclose(dateta);
		fclose(datP);
		fclose(datvel);
		fclose(datacc);
		fclose(datPhi);
	}
	
	free_2darray_complex(etahat,Nx,Ny);
	free_2darray_complex(uhat,Nx,Ny);
	free_2darray_complex(vhat,Nx,Ny);
	free_2darray_complex(dtetahat,Nx,Ny);
	free_2darray(eta,Nx,Ny);
	free_2darray(u,Nx,Ny);
	free_2darray(v,Nx,Ny);
	free_2darray(dteta,Nx,Ny);
	
	delete[] tim;	
	gsl_odeiv2_evolve_free(e);
	gsl_odeiv2_control_free(c);
	gsl_odeiv2_step_free(s); 
	
	free_2darray_complex(Source_hat,Nx,Ny);
}

void fun_termshat_AB1_2DB_uv(fftw_complex* etah,fftw_complex* uh,fftw_complex* vh,fftw_complex* deluH1hat,fftw_complex* delvH1hat,
fftw_complex* dampetahat,fftw_complex* dampuhat,fftw_complex* dampvhat,domvar dominfo,Oprtvar Oprt,const double* Zhat)
{
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* eta_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* u_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* v_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	
	// Create plans 
	fftw_plan plan_iftetahat   	= fftw_plan_dft_2d(Ny,Nx, eta_hat, etah, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftuhat   	= fftw_plan_dft_2d(Ny,Nx, u_hat, uh, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftvhat   	= fftw_plan_dft_2d(Ny,Nx, v_hat, vh, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    #pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		eta_hat[i*Nx+j][0]=Zhat[i*Nx+j];
		eta_hat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
		u_hat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
		u_hat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
		v_hat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
		v_hat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];		
	  }
	}
	
	fftwcomplex_2d_sym(eta_hat,Nx,Ny);
	fftwcomplex_2d_sym(u_hat,Nx,Ny);
	fftwcomplex_2d_sym(v_hat,Nx,Ny);
	fftw_execute(plan_iftetahat);
    fftw_execute(plan_iftuhat);  
    fftw_execute(plan_iftvhat);    
    fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftuhat);
    fftw_destroy_plan(plan_iftvhat);
    
    //interpolation
    fftw_complex* Csqr_uhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Csqr_vhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	if (ninterp==2){
		HSSgen_2depth(uh,u_hat,Oprt.InterpD.Csqr2d_min,Oprt.InterpD.Csqr2d_plus,Csqr_uhat,Nx,Ny,factIFFT);
		HSSgen_2depth(vh,v_hat,Oprt.InterpD.Csqr2d_min,Oprt.InterpD.Csqr2d_plus,Csqr_vhat,Nx,Ny,factIFFT);
	}
	else {
		HSSgen_3depth(uh,u_hat,Oprt.InterpD.Csqr2d_min,Oprt.InterpD.Csqr2d_mid,Oprt.InterpD.Csqr2d_plus,Csqr_uhat,Nx,Ny,factIFFT);
		HSSgen_3depth(vh,v_hat,Oprt.InterpD.Csqr2d_min,Oprt.InterpD.Csqr2d_mid,Oprt.InterpD.Csqr2d_plus,Csqr_vhat,Nx,Ny,factIFFT);
	}   
	
	fftw_complex* dampeta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
    // Create plans 
	fftw_plan plan_fftdampeta  	= fftw_plan_dft_2d(Ny,Nx,dampeta,dampetahat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampu   	= fftw_plan_dft_2d(Ny,Nx,dampu,dampuhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampv   	= fftw_plan_dft_2d(Ny,Nx,dampv,dampvhat, FFTW_FORWARD, FFTW_ESTIMATE);		
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		dampeta[i*Nx+j][0]=7*Oprt.Cpeak*dominfo.fbdy.charac[i][j]*etah[i*Nx+j][0]/factIFFT;
		dampeta[i*Nx+j][1]=0;
		dampu[i*Nx+j][0]=7*Oprt.Cpeak*dominfo.fbdy.charac[i][j]*uh[i*Nx+j][0]/factIFFT;
		dampu[i*Nx+j][1]=0;	
		dampv[i*Nx+j][0]=7*Oprt.Cpeak*dominfo.fbdy.charac[i][j]*vh[i*Nx+j][0]/factIFFT;
		dampv[i*Nx+j][1]=0;		
		
		deluH1hat[i*Nx+j][0] = Csqr_uhat[i*Nx+j][1]/grav*dominfo.kx[j];
		deluH1hat[i*Nx+j][1] = -Csqr_uhat[i*Nx+j][0]/grav*dominfo.kx[j];
		delvH1hat[i*Nx+j][0] = Csqr_vhat[i*Nx+j][1]/grav*dominfo.ky[i];
		delvH1hat[i*Nx+j][1] = -Csqr_vhat[i*Nx+j][0]/grav*dominfo.ky[i];
		
	  } 
	}
    
	fftw_execute(plan_fftdampeta);
    fftw_execute(plan_fftdampu);
    fftw_execute(plan_fftdampv);	
	fftw_destroy_plan(plan_fftdampeta);
    fftw_destroy_plan(plan_fftdampu); 
    fftw_destroy_plan(plan_fftdampv);   
    
	fftw_free(Csqr_uhat);fftw_free(Csqr_vhat);
	fftw_free(dampeta);fftw_free(dampu);fftw_free(dampv);
	fftw_free(eta_hat);fftw_free(u_hat);fftw_free(v_hat);
}

void fun_termshat_AB1_2DF_uv(fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* dampetahat,fftw_complex* dampuhat,fftw_complex* dampvhat,domvar dominfo,double Cpeak,const double* Zhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex* dampeta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
	// Create plans 
	fftw_plan plan_iftetahat = fftw_plan_dft_2d(Ny,Nx, etahat, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftuhat   = fftw_plan_dft_2d(Ny,Nx, uhat, u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftvhat   = fftw_plan_dft_2d(Ny,Nx, vhat, v, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampeta= fftw_plan_dft_2d(Ny,Nx,dampeta,dampetahat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampu  = fftw_plan_dft_2d(Ny,Nx,dampu,dampuhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampv  = fftw_plan_dft_2d(Ny,Nx,dampv,dampvhat, FFTW_FORWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];	 
		}
	}
		
	fftwcomplex_2d_sym(etahat,Nx,Ny);
	fftwcomplex_2d_sym(uhat,Nx,Ny);
	fftwcomplex_2d_sym(vhat,Nx,Ny);	
	fftw_execute(plan_iftetahat);
	fftw_execute(plan_iftuhat);
	fftw_execute(plan_iftvhat); 
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dampeta[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*eta[i*Nx+j][0]/factIFFT;
			dampeta[i*Nx+j][1]=0;
			dampu[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*u[i*Nx+j][0]/factIFFT;
			dampu[i*Nx+j][1]=0;
			dampv[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*v[i*Nx+j][0]/factIFFT;
			dampv[i*Nx+j][1]=0;
		}
	}
	fftw_execute(plan_fftdampeta);
	fftw_execute(plan_fftdampu);
	fftw_execute(plan_fftdampv);
	
	fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftuhat);
    fftw_destroy_plan(plan_iftvhat);
    fftw_destroy_plan(plan_fftdampeta);
    fftw_destroy_plan(plan_fftdampu);
    fftw_destroy_plan(plan_fftdampv);
	
	fftw_free(dampeta);fftw_free(dampu);fftw_free(dampv);
	fftw_free(etahat);fftw_free(uhat);fftw_free(vhat);
}
/*
void fun_termshat_AB2_2DB_uv(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* delphiH1hat,fftw_complex* delphiH2hat,
fftw_complex* deletaH2hat,domvar dominfo,Oprtvar Oprt,const double* Zhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex* Lphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaM0uv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0uvhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

    // Create plans 	
	fftw_plan plan_iftLphihat  = fftw_plan_dft_2d(Ny,Nx, delphiH1hat, Lphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwcomplex_2d_sym(delphiH1hat,Nx,Ny);
    fftw_execute(plan_iftLphihat);
    fftw_destroy_plan(plan_iftLphihat);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];

			etaM0uv[i*Nx+j][0]=eta[i*Nx+j][0]*Lphi[i*Nx+j][0]/pow(factIFFT,2);
			etaM0uv[i*Nx+j][1]=0;			
			eta_u[i*Nx+j][0]=eta[i*Nx+j][0]*u[i*Nx+j][0]/pow(factIFFT,2);
			eta_u[i*Nx+j][1]=0;
			eta_v[i*Nx+j][0]=eta[i*Nx+j][0]*v[i*Nx+j][0]/pow(factIFFT,2);
			eta_v[i*Nx+j][1]=0;
		}
	}	
		
	// Create plans 
	fftw_plan plan_fftetaM0uv   = fftw_plan_dft_2d(Ny,Nx,etaM0uv,etaM0uvhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_u		= fftw_plan_dft_2d(Ny,Nx,eta_u,eta_uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_v		= fftw_plan_dft_2d(Ny,Nx,eta_v,eta_vhat,FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(plan_fftetaM0uv);
	fftw_execute(plan_ffteta_u);
	fftw_execute(plan_ffteta_v);

	fftw_destroy_plan(plan_fftetaM0uv);
	fftw_destroy_plan(plan_ffteta_u);
	fftw_destroy_plan(plan_ffteta_v);

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etaM0uv[i*Nx+j][0]=etaM0uv[i*Nx+j][0]*factIFFT;//for operator			
		}
	}	

	//interpolation
    fftw_complex* Csqr_etaM0uvhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* LetaM0uv_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	if (ninterp==2){
		HSSgen_2depth(etaM0uv,etaM0uvhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,Csqr_etaM0uvhat,Nx,Ny,factIFFT);
	}
	else {
		HSSgen_3depth(etaM0uv,etaM0uvhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,Csqr_etaM0uvhat,Nx,Ny,factIFFT);
	}

	fftw_complex* deletaH2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));    
	fftw_plan plan_fftdeletaH2	= fftw_plan_dft_2d(Ny,Nx,deletaH2,deletaH2hat,FFTW_FORWARD, FFTW_ESTIMATE);	

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){			
			LetaM0uv_hat[i*Nx+j][0]=Csqr_etaM0uvhat[i*Nx+j][0]/grav;
			LetaM0uv_hat[i*Nx+j][1]=Csqr_etaM0uvhat[i*Nx+j][1]/grav;
			delphiH2hat[i*Nx+j][0] =-(LetaM0uv_hat[i*Nx+j][0]-dominfo.kx[j]*eta_uhat[i*Nx+j][1]-dominfo.ky[i]*eta_vhat[i*Nx+j][1]);//
			delphiH2hat[i*Nx+j][1] =-(LetaM0uv_hat[i*Nx+j][1]+dominfo.kx[j]*eta_uhat[i*Nx+j][0]+dominfo.ky[i]*eta_vhat[i*Nx+j][0]);//

			deletaH2[i*Nx+j][0] = 0.5*(pow(u[i*Nx+j][0]/factIFFT,2)+pow(v[i*Nx+j][0]/factIFFT,2)-pow(Lphi[i*Nx+j][0]/factIFFT,2));
			deletaH2[i*Nx+j][1] = 0;
		}
	}
	
	fftw_execute(plan_fftdeletaH2);	
	fftw_destroy_plan(plan_fftdeletaH2);	
	
	fftw_free(Csqr_etaM0uvhat);fftw_free(LetaM0uv_hat);
	
	fftw_free(uhat);
	fftw_free(eta_u); fftw_free(eta_uhat); 
	fftw_free(etaM0uv);fftw_free(etaM0uvhat);

	fftw_free(vhat);
	fftw_free(eta_v); fftw_free(eta_vhat); 
	
	fftw_free(etahat);
	fftw_free(deletaH2);
	fftw_free(Lphi);
}
*/

void fun_termshat_AB2_2DF_uv(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,
fftw_complex* ikxdeletaH2_hatadj,fftw_complex* ikydeletaH2_hatadj,fftw_complex* deluH2hat_adj,fftw_complex* delvH2hat_adj,
fftw_complex* dampetahat,fftw_complex* dampuhat,fftw_complex* dampvhat,domvar dominfo,double Cpeak,const double* Zhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	// Create plans 
	fftw_plan plan_iftetahat = fftw_plan_dft_2d(Ny,Nx, etahat, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftuhat   = fftw_plan_dft_2d(Ny,Nx, uhat, u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftvhat   = fftw_plan_dft_2d(Ny,Nx, vhat, v, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftLphihat= fftw_plan_dft_2d(Ny,Nx, Lphihat, Lphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftM0uhat = fftw_plan_dft_2d(Ny,Nx, M0uhat, M0u, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_iftM0vhat = fftw_plan_dft_2d(Ny,Nx, M0vhat, M0v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];		
			
			M0uhat[i*Nx+j][0]=pow(Oprt.C2d[i][j],2)/grav*dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j];
			M0uhat[i*Nx+j][1]=-pow(Oprt.C2d[i][j],2)/grav*dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j];
			M0vhat[i*Nx+j][0]=pow(Oprt.C2d[i][j],2)/grav*dominfo.ky[i]*Zhat[5*Nx*Ny+i*Nx+j];
			M0vhat[i*Nx+j][1]=-pow(Oprt.C2d[i][j],2)/grav*dominfo.ky[i]*Zhat[4*Nx*Ny+i*Nx+j];
			
			Lphihat[i*Nx+j][0] = M0uhat[i*Nx+j][0]+M0vhat[i*Nx+j][0];
			Lphihat[i*Nx+j][1] = M0uhat[i*Nx+j][1]+M0vhat[i*Nx+j][1];			
		}
	}
	
	fftwcomplex_2d_sym(etahat,Nx,Ny);
	fftwcomplex_2d_sym(uhat,Nx,Ny);
	fftwcomplex_2d_sym(vhat,Nx,Ny);
	fftwcomplex_2d_sym(Lphihat,Nx,Ny);
	fftwcomplex_2d_sym(M0uhat,Nx,Ny);
	fftwcomplex_2d_sym(M0vhat,Nx,Ny);	

	fftw_execute(plan_iftetahat);
    fftw_execute(plan_iftuhat);
    fftw_execute(plan_iftvhat);
    fftw_execute(plan_iftLphihat);
    fftw_execute(plan_iftM0uhat);
	fftw_execute(plan_iftM0vhat);

	fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftuhat);
    fftw_destroy_plan(plan_iftvhat);
    fftw_destroy_plan(plan_iftLphihat);
	fftw_destroy_plan(plan_iftM0uhat);
	fftw_destroy_plan(plan_iftM0vhat);
	
	fftw_complex* dampeta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    fftw_complex* etaM0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    fftw_complex* deletaH2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* deletaH2hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    // Create plans 
    fftw_plan plan_fftdampeta   = fftw_plan_dft_2d(Ny,Nx,dampeta,dampetahat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampu     = fftw_plan_dft_2d(Ny,Nx,dampu,dampuhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampv     = fftw_plan_dft_2d(Ny,Nx,dampv,dampvhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0u    = fftw_plan_dft_2d(Ny,Nx,etaM0u,etaM0uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_u		= fftw_plan_dft_2d(Ny,Nx,eta_u,eta_uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0v    = fftw_plan_dft_2d(Ny,Nx,etaM0v,etaM0vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_v		= fftw_plan_dft_2d(Ny,Nx,eta_v,eta_vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdeletaH2	= fftw_plan_dft_2d(Ny,Nx,deletaH2,deletaH2hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dampeta[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*eta[i*Nx+j][0]/factIFFT;
			dampeta[i*Nx+j][1]=0;
			dampu[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*u[i*Nx+j][0]/factIFFT;
			dampu[i*Nx+j][1]=0;
			dampv[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*v[i*Nx+j][0]/factIFFT;
			dampv[i*Nx+j][1]=0;
			
			etaM0u[i*Nx+j][0]=eta[i*Nx+j][0]*M0u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0u[i*Nx+j][1]=0;
			eta_u[i*Nx+j][0]=eta[i*Nx+j][0]*u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_u[i*Nx+j][1]=0;
			
			etaM0v[i*Nx+j][0]=eta[i*Nx+j][0]*M0v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0v[i*Nx+j][1]=0;
			eta_v[i*Nx+j][0]=eta[i*Nx+j][0]*v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_v[i*Nx+j][1]=0;

			deletaH2[i*Nx+j][0]= 0.5*((pow(u[i*Nx+j][0],2)+pow(v[i*Nx+j][0],2))-pow(Lphi[i*Nx+j][0],2))/pow(factIFFT,2)*dominfo.wallchar[i][j];
			deletaH2[i*Nx+j][1]= 0;
		}
	}
	fftw_execute(plan_fftdampeta);
	fftw_execute(plan_fftdampu);
	fftw_execute(plan_fftdampv);
	fftw_execute(plan_fftetaM0u);
	fftw_execute(plan_ffteta_u);
	fftw_execute(plan_fftetaM0v);
	fftw_execute(plan_ffteta_v);
	fftw_execute(plan_fftdeletaH2);
	
	fftw_destroy_plan(plan_fftdampeta);
	fftw_destroy_plan(plan_fftdampu);
	fftw_destroy_plan(plan_fftdampv);
	fftw_destroy_plan(plan_fftetaM0u);
	fftw_destroy_plan(plan_ffteta_u);
	fftw_destroy_plan(plan_fftetaM0v);
	fftw_destroy_plan(plan_ffteta_v);
	fftw_destroy_plan(plan_fftdeletaH2);	
	
	fftw_complex* LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* deluH2_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0v_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* delvH2_hat 	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikxdeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			ikxdeletaH2_hat[i*Nx+j][0] = -dominfo.kx[j]*deletaH2hat[i*Nx+j][1];
			ikxdeletaH2_hat[i*Nx+j][1] = dominfo.kx[j]*deletaH2hat[i*Nx+j][0];
			
			ikydeletaH2_hat[i*Nx+j][0] = -dominfo.ky[i]*deletaH2hat[i*Nx+j][1];
			ikydeletaH2_hat[i*Nx+j][1] = dominfo.ky[i]*deletaH2hat[i*Nx+j][0];
			
			LetaM0u_hat[i*Nx+j][0]=Oprt.L2d[i][j]*etaM0uhat[i*Nx+j][0];
			LetaM0u_hat[i*Nx+j][1]=Oprt.L2d[i][j]*etaM0uhat[i*Nx+j][1];			

			LetaM0v_hat[i*Nx+j][0]=Oprt.L2d[i][j]*etaM0vhat[i*Nx+j][0];
			LetaM0v_hat[i*Nx+j][1]=Oprt.L2d[i][j]*etaM0vhat[i*Nx+j][1];
			
			deluH2_hat[i*Nx+j][0] = -LetaM0u_hat[i*Nx+j][0]+dominfo.kx[j]*eta_uhat[i*Nx+j][1];
			deluH2_hat[i*Nx+j][1] = -LetaM0u_hat[i*Nx+j][1]-dominfo.kx[j]*eta_uhat[i*Nx+j][0];

			delvH2_hat[i*Nx+j][0] = -LetaM0v_hat[i*Nx+j][0]+dominfo.ky[i]*eta_vhat[i*Nx+j][1];
			delvH2_hat[i*Nx+j][1] = -LetaM0v_hat[i*Nx+j][1]-dominfo.ky[i]*eta_vhat[i*Nx+j][0];
		}
	}
	fun_adjust(deluH2hat_adj,deluH2_hat,Nx,Ny);
	fun_adjust(delvH2hat_adj,delvH2_hat,Nx,Ny);
	fun_adjust(ikxdeletaH2_hatadj,ikxdeletaH2_hat,Nx,Ny);	
	fun_adjust(ikydeletaH2_hatadj,ikydeletaH2_hat,Nx,Ny);
	
	fftw_free(ikxdeletaH2_hat);
	fftw_free(ikydeletaH2_hat);
	
	fftw_free(etahat);
	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(Lphi);fftw_free(Lphihat);	
	fftw_free(dampeta);fftw_free(dampu);fftw_free(dampv);
	fftw_free(M0uhat);fftw_free(M0u);
	fftw_free(M0vhat);fftw_free(M0v);
	fftw_free(etaM0u);fftw_free(etaM0uhat);
	fftw_free(etaM0v);fftw_free(etaM0vhat);
	fftw_free(eta_u); fftw_free(eta_uhat); 
	fftw_free(eta_v); fftw_free(eta_vhat); 			
	fftw_free(LetaM0u_hat);
	fftw_free(LetaM0v_hat);	
	fftw_free(deluH2_hat);
	fftw_free(delvH2_hat);
	fftw_free(deletaH2);
	fftw_free(deletaH2hat);
}

void fun_termshat_AB2_2DB_uv(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* deluH1hat,fftw_complex* delvH1hat,
fftw_complex* ikxdeletaH2_hatadj,fftw_complex* ikydeletaH2_hatadj,fftw_complex* deluH2hat_adj,fftw_complex* delvH2hat_adj,domvar dominfo,double Cpeak,const double* Zhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	// Create plans 
	fftw_plan plan_iftetahat = fftw_plan_dft_2d(Ny,Nx, etahat, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftuhat   = fftw_plan_dft_2d(Ny,Nx, uhat, u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftvhat   = fftw_plan_dft_2d(Ny,Nx, vhat, v, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftLphihat= fftw_plan_dft_2d(Ny,Nx, Lphihat, Lphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftM0uhat = fftw_plan_dft_2d(Ny,Nx, M0uhat, M0u, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_iftM0vhat = fftw_plan_dft_2d(Ny,Nx, M0vhat, M0v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];		
			
			M0uhat[i*Nx+j][0]=deluH1hat[i*Nx+j][0];
			M0uhat[i*Nx+j][1]=deluH1hat[i*Nx+j][1];
			M0vhat[i*Nx+j][0]=delvH1hat[i*Nx+j][0];
			M0vhat[i*Nx+j][1]=delvH1hat[i*Nx+j][1];
			
			Lphihat[i*Nx+j][0] = deluH1hat[i*Nx+j][0]+delvH1hat[i*Nx+j][0];
			Lphihat[i*Nx+j][1] = deluH1hat[i*Nx+j][1]+delvH1hat[i*Nx+j][1];			
		}
	}
	
	fftwcomplex_2d_sym(etahat,Nx,Ny);
	fftwcomplex_2d_sym(uhat,Nx,Ny);
	fftwcomplex_2d_sym(vhat,Nx,Ny);
	fftwcomplex_2d_sym(Lphihat,Nx,Ny);
	fftwcomplex_2d_sym(M0uhat,Nx,Ny);
	fftwcomplex_2d_sym(M0vhat,Nx,Ny);	

	fftw_execute(plan_iftetahat);
    fftw_execute(plan_iftuhat);
    fftw_execute(plan_iftvhat);
    fftw_execute(plan_iftLphihat);
    fftw_execute(plan_iftM0uhat);
	fftw_execute(plan_iftM0vhat);

	fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftuhat);
    fftw_destroy_plan(plan_iftvhat);
    fftw_destroy_plan(plan_iftLphihat);
	fftw_destroy_plan(plan_iftM0uhat);
	fftw_destroy_plan(plan_iftM0vhat);
	
    fftw_complex* etaM0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    fftw_complex* deletaH2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* deletaH2hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    // Create plans 
	fftw_plan plan_fftetaM0u    = fftw_plan_dft_2d(Ny,Nx,etaM0u,etaM0uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_u		= fftw_plan_dft_2d(Ny,Nx,eta_u,eta_uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0v    = fftw_plan_dft_2d(Ny,Nx,etaM0v,etaM0vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_v		= fftw_plan_dft_2d(Ny,Nx,eta_v,eta_vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdeletaH2	= fftw_plan_dft_2d(Ny,Nx,deletaH2,deletaH2hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			
			etaM0u[i*Nx+j][0]=eta[i*Nx+j][0]*M0u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0u[i*Nx+j][1]=0;
			eta_u[i*Nx+j][0]=eta[i*Nx+j][0]*u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_u[i*Nx+j][1]=0;
			
			etaM0v[i*Nx+j][0]=eta[i*Nx+j][0]*M0v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0v[i*Nx+j][1]=0;
			eta_v[i*Nx+j][0]=eta[i*Nx+j][0]*v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_v[i*Nx+j][1]=0;

			deletaH2[i*Nx+j][0]= 0.5*((pow(u[i*Nx+j][0],2)+pow(v[i*Nx+j][0],2))-pow(Lphi[i*Nx+j][0],2))/pow(factIFFT,2)*dominfo.wallchar[i][j];
			deletaH2[i*Nx+j][1]= 0;
		}
	}
	fftw_execute(plan_fftetaM0u);
	fftw_execute(plan_ffteta_u);
	fftw_execute(plan_fftetaM0v);
	fftw_execute(plan_ffteta_v);
	fftw_execute(plan_fftdeletaH2);
	
	fftw_destroy_plan(plan_fftetaM0u);
	fftw_destroy_plan(plan_ffteta_u);
	fftw_destroy_plan(plan_fftetaM0v);
	fftw_destroy_plan(plan_ffteta_v);
	fftw_destroy_plan(plan_fftdeletaH2);	
	
	fftw_complex* LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* deluH2_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0v_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* delvH2_hat 	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex* ikxdeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	//interpolation
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etaM0u[i*Nx+j][0]=etaM0u[i*Nx+j][0]*factIFFT;//for operator	
			etaM0v[i*Nx+j][0]=etaM0v[i*Nx+j][0]*factIFFT;//for operator			
		}
	}		
	
	if (ninterp==2){
		HSSgen_2depth(etaM0u,etaM0uhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,LetaM0u_hat,Nx,Ny,factIFFT);
		HSSgen_2depth(etaM0v,etaM0vhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,LetaM0v_hat,Nx,Ny,factIFFT);
	}
	else {
		HSSgen_3depth(etaM0u,etaM0uhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,LetaM0u_hat,Nx,Ny,factIFFT);
		HSSgen_3depth(etaM0v,etaM0vhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,LetaM0v_hat,Nx,Ny,factIFFT);
	}
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			ikxdeletaH2_hat[i*Nx+j][0] = -dominfo.kx[j]*deletaH2hat[i*Nx+j][1];
			ikxdeletaH2_hat[i*Nx+j][1] = dominfo.kx[j]*deletaH2hat[i*Nx+j][0];
			
			ikydeletaH2_hat[i*Nx+j][0] = -dominfo.ky[i]*deletaH2hat[i*Nx+j][1];
			ikydeletaH2_hat[i*Nx+j][1] = dominfo.ky[i]*deletaH2hat[i*Nx+j][0];
			
			deluH2_hat[i*Nx+j][0] = dominfo.kx[j]*eta_uhat[i*Nx+j][1]-LetaM0u_hat[i*Nx+j][0]/grav;
			deluH2_hat[i*Nx+j][1] = -dominfo.kx[j]*eta_uhat[i*Nx+j][0]-LetaM0u_hat[i*Nx+j][1]/grav;

			delvH2_hat[i*Nx+j][0] = dominfo.ky[i]*eta_vhat[i*Nx+j][1]-LetaM0v_hat[i*Nx+j][0]/grav;
			delvH2_hat[i*Nx+j][1] = -dominfo.ky[i]*eta_vhat[i*Nx+j][0]-LetaM0v_hat[i*Nx+j][1]/grav;
		}
	}
	fun_adjust(deluH2hat_adj,deluH2_hat,Nx,Ny);
	fun_adjust(delvH2hat_adj,delvH2_hat,Nx,Ny);
	fun_adjust(ikxdeletaH2_hatadj,ikxdeletaH2_hat,Nx,Ny);	
	fun_adjust(ikydeletaH2_hatadj,ikydeletaH2_hat,Nx,Ny);
	
	fftw_free(ikxdeletaH2_hat);
	fftw_free(ikydeletaH2_hat);
	
	fftw_free(etahat);
	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(Lphi);fftw_free(Lphihat);	
	
	fftw_free(M0uhat);fftw_free(M0u);
	fftw_free(M0vhat);fftw_free(M0v);
	fftw_free(etaM0u);fftw_free(etaM0uhat);
	fftw_free(etaM0v);fftw_free(etaM0vhat);
	fftw_free(eta_u); fftw_free(eta_uhat); 
	fftw_free(eta_v); fftw_free(eta_vhat); 			
	fftw_free(LetaM0u_hat);
	fftw_free(LetaM0v_hat);	
	fftw_free(deluH2_hat);
	fftw_free(delvH2_hat);
	fftw_free(deletaH2);	
	fftw_free(deletaH2hat);
}


void fun_termshat_AB3_2DF_uv(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* ikxdeletaH3_hatadj,fftw_complex* ikydeletaH3_hatadj,
fftw_complex* deluH3hat_adj,fftw_complex* delvH3hat_adj,fftw_complex* ikxdeletaH2_hatadj,fftw_complex* ikydeletaH2_hatadj,
fftw_complex* deluH2hat_adj,fftw_complex* delvH2hat_adj,fftw_complex* dampetahat,fftw_complex* dampuhat,fftw_complex* dampvhat,domvar dominfo,double Cpeak,const double* Zhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	// Create plans 
	fftw_plan plan_iftetahat = fftw_plan_dft_2d(Ny,Nx, etahat, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftuhat   = fftw_plan_dft_2d(Ny,Nx, uhat, u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftvhat   = fftw_plan_dft_2d(Ny,Nx, vhat, v, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftLphihat= fftw_plan_dft_2d(Ny,Nx, Lphihat, Lphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftM0uhat = fftw_plan_dft_2d(Ny,Nx, M0uhat, M0u, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_iftM0vhat = fftw_plan_dft_2d(Ny,Nx, M0vhat, M0v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];		
			
			M0uhat[i*Nx+j][0]=pow(Oprt.C2d[i][j],2)/grav*dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j];
			M0uhat[i*Nx+j][1]=-pow(Oprt.C2d[i][j],2)/grav*dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j];
			M0vhat[i*Nx+j][0]=pow(Oprt.C2d[i][j],2)/grav*dominfo.ky[i]*Zhat[5*Nx*Ny+i*Nx+j];
			M0vhat[i*Nx+j][1]=-pow(Oprt.C2d[i][j],2)/grav*dominfo.ky[i]*Zhat[4*Nx*Ny+i*Nx+j];
			Lphihat[i*Nx+j][0] = M0uhat[i*Nx+j][0]+M0vhat[i*Nx+j][0];
			Lphihat[i*Nx+j][1] = M0uhat[i*Nx+j][1]+M0vhat[i*Nx+j][1];			
		}
	}
	fftwcomplex_2d_sym(etahat,Nx,Ny);
	fftwcomplex_2d_sym(uhat,Nx,Ny);
	fftwcomplex_2d_sym(vhat,Nx,Ny);
	fftwcomplex_2d_sym(Lphihat,Nx,Ny);
	fftwcomplex_2d_sym(M0uhat,Nx,Ny);
	fftwcomplex_2d_sym(M0vhat,Nx,Ny);	

	fftw_execute(plan_iftetahat);
    fftw_execute(plan_iftuhat);
    fftw_execute(plan_iftvhat);
    fftw_execute(plan_iftLphihat);
    fftw_execute(plan_iftM0uhat);
	fftw_execute(plan_iftM0vhat);

	fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftuhat);
    fftw_destroy_plan(plan_iftvhat);
    fftw_destroy_plan(plan_iftLphihat);
	fftw_destroy_plan(plan_iftM0uhat);
	fftw_destroy_plan(plan_iftM0vhat);
	
	fftw_complex* dampeta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    fftw_complex* etaM0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* deletaH2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* deletaH2hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    // Create plans 
    fftw_plan plan_fftdampeta   = fftw_plan_dft_2d(Ny,Nx,dampeta,dampetahat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampu     = fftw_plan_dft_2d(Ny,Nx,dampu,dampuhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampv     = fftw_plan_dft_2d(Ny,Nx,dampv,dampvhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0u    = fftw_plan_dft_2d(Ny,Nx,etaM0u,etaM0uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_u		= fftw_plan_dft_2d(Ny,Nx,eta_u,eta_uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0v    = fftw_plan_dft_2d(Ny,Nx,etaM0v,etaM0vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_v		= fftw_plan_dft_2d(Ny,Nx,eta_v,eta_vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdeletaH2	= fftw_plan_dft_2d(Ny,Nx,deletaH2,deletaH2hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dampeta[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*eta[i*Nx+j][0]/factIFFT;
			dampeta[i*Nx+j][1]=0;
			dampu[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*u[i*Nx+j][0]/factIFFT;
			dampu[i*Nx+j][1]=0;
			dampv[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*v[i*Nx+j][0]/factIFFT;
			dampv[i*Nx+j][1]=0;
			
			etaM0u[i*Nx+j][0]=eta[i*Nx+j][0]*M0u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0u[i*Nx+j][1]=0;
			eta_u[i*Nx+j][0]=eta[i*Nx+j][0]*u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_u[i*Nx+j][1]=0;
			
			etaM0v[i*Nx+j][0]=eta[i*Nx+j][0]*M0v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0v[i*Nx+j][1]=0;
			eta_v[i*Nx+j][0]=eta[i*Nx+j][0]*v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_v[i*Nx+j][1]=0;

			deletaH2[i*Nx+j][0]= 0.5*((pow(u[i*Nx+j][0],2)+pow(v[i*Nx+j][0],2))-pow(Lphi[i*Nx+j][0],2))/pow(factIFFT,2)
								*dominfo.wallchar[i][j];
			deletaH2[i*Nx+j][1]= 0;
		}
	}
	fftw_execute(plan_fftdampeta);
	fftw_execute(plan_fftdampu);
	fftw_execute(plan_fftdampv);
	fftw_execute(plan_fftetaM0u);
	fftw_execute(plan_ffteta_u);
	fftw_execute(plan_fftetaM0v);
	fftw_execute(plan_ffteta_v);
	fftw_execute(plan_fftdeletaH2);
	
	fftw_destroy_plan(plan_fftdampeta);
	fftw_destroy_plan(plan_fftdampu);
	fftw_destroy_plan(plan_fftdampv);
	fftw_destroy_plan(plan_fftetaM0u);
	fftw_destroy_plan(plan_ffteta_u);
	fftw_destroy_plan(plan_fftetaM0v);
	fftw_destroy_plan(plan_ffteta_v);
	fftw_destroy_plan(plan_fftdeletaH2);	
	
	fftw_complex* LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* deluH2_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0v_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* delvH2_hat 	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_plan plan_iftLetaM0u_hat  = fftw_plan_dft_2d(Ny,Nx, LetaM0u_hat, LetaM0u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftLetaM0v_hat  = fftw_plan_dft_2d(Ny,Nx, LetaM0v_hat, LetaM0v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* LetaLphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaLphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_iftLetaLphihat= fftw_plan_dft_2d(Ny,Nx, LetaLphihat, LetaLphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* divgradphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* divgradphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_plan plan_iftdivgradphihat= fftw_plan_dft_2d(Ny,Nx, divgradphihat, divgradphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* eta2M0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta2M0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_ffteta2M0u     = fftw_plan_dft_2d(Ny,Nx,eta2M0u,eta2M0uhat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftw_complex* eta2M0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta2M0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_ffteta2M0v     = fftw_plan_dft_2d(Ny,Nx,eta2M0v,eta2M0vhat,FFTW_FORWARD, FFTW_ESTIMATE);		
	
	fftw_complex* dxu_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* dxu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_iftdxu_hat= fftw_plan_dft_2d(Ny,Nx, dxu_hat, dxu, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* dyv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* dyv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_iftdyv_hat= fftw_plan_dft_2d(Ny,Nx, dyv_hat, dyv, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* ikxdeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			ikxdeletaH2_hat[i*Nx+j][0] = -dominfo.kx[j]*deletaH2hat[i*Nx+j][1];
			ikxdeletaH2_hat[i*Nx+j][1] = dominfo.kx[j]*deletaH2hat[i*Nx+j][0];
			
			ikydeletaH2_hat[i*Nx+j][0] = -dominfo.ky[i]*deletaH2hat[i*Nx+j][1];
			ikydeletaH2_hat[i*Nx+j][1] = dominfo.ky[i]*deletaH2hat[i*Nx+j][0];
			
			LetaM0u_hat[i*Nx+j][0]=Oprt.L2d[i][j]*etaM0uhat[i*Nx+j][0];
			LetaM0u_hat[i*Nx+j][1]=Oprt.L2d[i][j]*etaM0uhat[i*Nx+j][1];			

			LetaM0v_hat[i*Nx+j][0]=Oprt.L2d[i][j]*etaM0vhat[i*Nx+j][0];
			LetaM0v_hat[i*Nx+j][1]=Oprt.L2d[i][j]*etaM0vhat[i*Nx+j][1];
			
			LetaLphihat[i*Nx+j][0]=LetaM0u_hat[i*Nx+j][0]+LetaM0v_hat[i*Nx+j][0];//Oprt.L2d[i][j]*etaLphihat[i*Nx+j][0];
			LetaLphihat[i*Nx+j][1]=LetaM0u_hat[i*Nx+j][1]+LetaM0v_hat[i*Nx+j][1];//Oprt.L2d[i][j]*etaLphihat[i*Nx+j][1];
			
			deluH2_hat[i*Nx+j][0] = -LetaM0u_hat[i*Nx+j][0]+dominfo.kx[j]*eta_uhat[i*Nx+j][1];
			deluH2_hat[i*Nx+j][1] = -LetaM0u_hat[i*Nx+j][1]-dominfo.kx[j]*eta_uhat[i*Nx+j][0];

			delvH2_hat[i*Nx+j][0] = -LetaM0v_hat[i*Nx+j][0]+dominfo.ky[i]*eta_vhat[i*Nx+j][1];
			delvH2_hat[i*Nx+j][1] = -LetaM0v_hat[i*Nx+j][1]-dominfo.ky[i]*eta_vhat[i*Nx+j][0];
			
			divgradphihat[i*Nx+j][0]=-dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j]-dominfo.ky[i]*Zhat[5*Nx*Ny+i*Nx+j];
			divgradphihat[i*Nx+j][1]=dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j]+dominfo.ky[i]*Zhat[4*Nx*Ny+i*Nx+j];
			
			eta2M0u[i*Nx+j][0]=pow(eta[i*Nx+j][0],2)*M0u[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
			eta2M0u[i*Nx+j][1]=0;
			
			eta2M0v[i*Nx+j][0]=pow(eta[i*Nx+j][0],2)*M0v[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
			eta2M0v[i*Nx+j][1]=0;
			
			dxu_hat[i*Nx+j][0] = -dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j];
			dxu_hat[i*Nx+j][1] = dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j];
			
			dyv_hat[i*Nx+j][0] = -dominfo.ky[i]*Zhat[5*Nx*Ny+i*Nx+j];
			dyv_hat[i*Nx+j][1] = dominfo.ky[i]*Zhat[4*Nx*Ny+i*Nx+j];
		}
	}
	fun_adjust(deluH2hat_adj,deluH2_hat,Nx,Ny);
	fun_adjust(delvH2hat_adj,delvH2_hat,Nx,Ny);
	fun_adjust(ikxdeletaH2_hatadj,ikxdeletaH2_hat,Nx,Ny);	
	fun_adjust(ikydeletaH2_hatadj,ikydeletaH2_hat,Nx,Ny);	
		
	fftwcomplex_2d_sym(LetaM0u_hat,Nx,Ny);
	fftwcomplex_2d_sym(LetaM0v_hat,Nx,Ny);
	fftwcomplex_2d_sym(divgradphihat,Nx,Ny);
	fftwcomplex_2d_sym(LetaLphihat,Nx,Ny);
	fftwcomplex_2d_sym(dxu_hat,Nx,Ny);
	fftwcomplex_2d_sym(dyv_hat,Nx,Ny);
	
	fftw_execute(plan_iftLetaM0u_hat);
	fftw_execute(plan_iftLetaM0v_hat);	
	fftw_execute(plan_iftdivgradphihat);
	fftw_execute(plan_ffteta2M0u);
	fftw_execute(plan_ffteta2M0v);
	fftw_execute(plan_iftLetaLphihat);
	fftw_execute(plan_iftdxu_hat);
	fftw_execute(plan_iftdyv_hat);
	
	fftw_destroy_plan(plan_iftLetaM0u_hat);
	fftw_destroy_plan(plan_iftLetaM0v_hat);	
	fftw_destroy_plan(plan_iftdivgradphihat);
	fftw_destroy_plan(plan_ffteta2M0u);
	fftw_destroy_plan(plan_ffteta2M0v);
	fftw_destroy_plan(plan_iftLetaLphihat);
	fftw_destroy_plan(plan_iftdxu_hat);
	fftw_destroy_plan(plan_iftdyv_hat);	

	fftw_complex* etaLetaM0u 	=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaLetaM0u_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaLetaM0v	=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaLetaM0v_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_fftetaLetaM0u     = fftw_plan_dft_2d(Ny,Nx,etaLetaM0u,etaLetaM0u_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_plan plan_fftetaLetaM0v     = fftw_plan_dft_2d(Ny,Nx,etaLetaM0v,etaLetaM0v_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftw_complex* eta2dxu  	 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta2dxu_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta2dyv	 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta2dyv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_ffteta2dxu= fftw_plan_dft_2d(Ny,Nx,eta2dxu,eta2dxu_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_plan plan_ffteta2dyv= fftw_plan_dft_2d(Ny,Nx,eta2dyv,eta2dyv_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftw_complex* deletaH3=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* deletaH3hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_fftdeletaH3 = fftw_plan_dft_2d(Ny,Nx,deletaH3,deletaH3hat,FFTW_FORWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  
		etaLetaM0u[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*LetaM0u[i*Nx+j][0]/factIFFT;
		etaLetaM0u[i*Nx+j][1] = 0;
		
		etaLetaM0v[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*LetaM0v[i*Nx+j][0]/factIFFT;
		etaLetaM0v[i*Nx+j][1] = 0;
		
		eta2dxu[i*Nx+j][0] = pow(eta[i*Nx+j][0],2)*dxu[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
		eta2dxu[i*Nx+j][1] = 0;
		
		eta2dyv[i*Nx+j][0] = pow(eta[i*Nx+j][0],2)*dyv[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
		eta2dyv[i*Nx+j][1] = 0;
        
		deletaH3[i*Nx+j][0]=(Lphi[i*Nx+j][0]/factIFFT*eta[i*Nx+j][0]/factIFFT*divgradphi[i*Nx+j][0]/factIFFT
								+Lphi[i*Nx+j][0]/factIFFT*LetaLphi[i*Nx+j][0]/factIFFT);
		deletaH3[i*Nx+j][1]=0;
	  }
	}	
	
	fftw_execute(plan_fftetaLetaM0u);
	fftw_execute(plan_fftetaLetaM0v);
	fftw_execute(plan_ffteta2dxu);
	fftw_execute(plan_ffteta2dyv);
	fftw_execute(plan_fftdeletaH3);
	
	fftw_destroy_plan(plan_fftetaLetaM0u);
	fftw_destroy_plan(plan_fftetaLetaM0v);
	fftw_destroy_plan(plan_ffteta2dxu);
	fftw_destroy_plan(plan_ffteta2dyv);
	fftw_destroy_plan(plan_fftdeletaH3);
	
	fftw_complex* deluH3temp_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* delvH3temp_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	//last round
	fftw_complex* ikxdeletaH3_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH3_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  ikxdeletaH3_hat[i*Nx+j][0] = -dominfo.kx[j]*deletaH3hat[i*Nx+j][1];
		  ikxdeletaH3_hat[i*Nx+j][1] = dominfo.kx[j]*deletaH3hat[i*Nx+j][0];
			
		  ikydeletaH3_hat[i*Nx+j][0] = -dominfo.ky[i]*deletaH3hat[i*Nx+j][1];
		  ikydeletaH3_hat[i*Nx+j][1] = dominfo.ky[i]*deletaH3hat[i*Nx+j][0];
		  		  	
		  deluH3temp_hat[i*Nx+j][0] = Oprt.L2d[i][j]*etaLetaM0u_hat[i*Nx+j][0]
										+0.5*Oprt.L2d[i][j]*eta2dxu_hat[i*Nx+j][0]										
										-0.5*(pow(dominfo.kx[j],2)*eta2M0uhat[i*Nx+j][0]);										
		  deluH3temp_hat[i*Nx+j][1] = Oprt.L2d[i][j]*etaLetaM0u_hat[i*Nx+j][1]		  
										+0.5*Oprt.L2d[i][j]*eta2dxu_hat[i*Nx+j][1]										
										-0.5*(pow(dominfo.kx[j],2)*eta2M0uhat[i*Nx+j][1]);
										
		  delvH3temp_hat[i*Nx+j][0] = Oprt.L2d[i][j]*etaLetaM0v_hat[i*Nx+j][0]		  
										+0.5*Oprt.L2d[i][j]*eta2dyv_hat[i*Nx+j][0]										
										-0.5*(pow(dominfo.ky[i],2)*eta2M0vhat[i*Nx+j][0]);
		  delvH3temp_hat[i*Nx+j][1] = Oprt.L2d[i][j]*etaLetaM0v_hat[i*Nx+j][1]		  
										+0.5*Oprt.L2d[i][j]*eta2dyv_hat[i*Nx+j][1]										
										-0.5*(pow(dominfo.ky[i],2)*eta2M0vhat[i*Nx+j][1]);
	  }
	}
	
	/*if (t>95){
				double* forline=new double[dominfo.Nx];
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=Lphi[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=LetaLphi[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				delete[] forline;
	}*/
	
	//adding nonlinear adjustment
	fun_adjust(deluH3hat_adj,deluH3temp_hat,Nx,Ny); fftw_free(deluH3temp_hat);	
	fun_adjust(delvH3hat_adj,delvH3temp_hat,Nx,Ny); fftw_free(delvH3temp_hat);
	fun_adjust(ikxdeletaH3_hatadj,ikxdeletaH3_hat,Nx,Ny);	
	fun_adjust(ikydeletaH3_hatadj,ikydeletaH3_hat,Nx,Ny);		
	
	fftw_free(ikxdeletaH2_hat);
	fftw_free(ikydeletaH2_hat);
	fftw_free(ikxdeletaH3_hat);
	fftw_free(ikydeletaH3_hat);
	fftw_free(etahat);
	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(Lphi);fftw_free(Lphihat);	
	fftw_free(dampeta);fftw_free(dampu);fftw_free(dampv);
	fftw_free(M0uhat);fftw_free(M0u);
	fftw_free(M0vhat);fftw_free(M0v);
	fftw_free(etaM0u);fftw_free(etaM0uhat);
	fftw_free(etaM0v);fftw_free(etaM0vhat);
	fftw_free(eta_u); fftw_free(eta_uhat); 
	fftw_free(eta_v); fftw_free(eta_vhat); 			
	fftw_free(LetaM0u_hat);
	fftw_free(LetaM0v_hat);	
	fftw_free(deluH2_hat);
	fftw_free(delvH2_hat);
	fftw_free(deletaH2);fftw_free(deletaH3);
	fftw_free(divgradphihat);
	fftw_free(divgradphi);	
	fftw_free(LetaLphihat);fftw_free(LetaLphi);
	
	fftw_free(LetaM0u);
	fftw_free(LetaM0v);
	fftw_free(eta2M0u);
	fftw_free(eta2M0uhat);	
	fftw_free(eta2M0v);
	fftw_free(eta2M0vhat);
	fftw_free(dxu_hat);
	fftw_free(dxu);
	fftw_free(dyv_hat);
	fftw_free(dyv);
	fftw_free(etaLetaM0u);
	fftw_free(etaLetaM0u_hat);
	fftw_free(etaLetaM0v);
	fftw_free(etaLetaM0v_hat);
	fftw_free(eta2dxu);
	fftw_free(eta2dxu_hat);
	fftw_free(eta2dyv);
	fftw_free(eta2dyv_hat);
	fftw_free(deletaH2hat);
	fftw_free(deletaH3hat);
}

void fun_termshat_AB3_2DB_uv(double t,fftw_complex* eta,fftw_complex* u,fftw_complex* v,fftw_complex* deluH1hat,fftw_complex* delvH1hat,
fftw_complex* ikxdeletaH3_hatadj,fftw_complex* ikydeletaH3_hatadj,fftw_complex* deluH3hat_adj,fftw_complex* delvH3hat_adj,
fftw_complex* ikxdeletaH2_hatadj,fftw_complex* ikydeletaH2_hatadj,fftw_complex* deluH2hat_adj,fftw_complex* delvH2hat_adj,
fftw_complex* dampetahat,fftw_complex* dampuhat,fftw_complex* dampvhat,domvar dominfo,double Cpeak,const double* Zhat)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex* etahat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* Lphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* M0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	// Create plans 
	fftw_plan plan_iftetahat = fftw_plan_dft_2d(Ny,Nx, etahat, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftuhat   = fftw_plan_dft_2d(Ny,Nx, uhat, u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftvhat   = fftw_plan_dft_2d(Ny,Nx, vhat, v, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftLphihat= fftw_plan_dft_2d(Ny,Nx, Lphihat, Lphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftM0uhat = fftw_plan_dft_2d(Ny,Nx, M0uhat, M0u, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_iftM0vhat = fftw_plan_dft_2d(Ny,Nx, M0vhat, M0v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etahat[i*Nx+j][0]=Zhat[i*Nx+j];
			etahat[i*Nx+j][1]=Zhat[Nx*Ny+i*Nx+j];	 
			uhat[i*Nx+j][0]=Zhat[2*Nx*Ny+i*Nx+j];
			uhat[i*Nx+j][1]=Zhat[3*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][0]=Zhat[4*Nx*Ny+i*Nx+j];
			vhat[i*Nx+j][1]=Zhat[5*Nx*Ny+i*Nx+j];		
			
			M0uhat[i*Nx+j][0]=deluH1hat[i*Nx+j][0];
			M0uhat[i*Nx+j][1]=deluH1hat[i*Nx+j][1];
			M0vhat[i*Nx+j][0]=delvH1hat[i*Nx+j][0];
			M0vhat[i*Nx+j][1]=delvH1hat[i*Nx+j][1];
			
			Lphihat[i*Nx+j][0] = M0uhat[i*Nx+j][0]+M0vhat[i*Nx+j][0];
			Lphihat[i*Nx+j][1] = M0uhat[i*Nx+j][1]+M0vhat[i*Nx+j][1];			
		}
	}
	fftwcomplex_2d_sym(etahat,Nx,Ny);
	fftwcomplex_2d_sym(uhat,Nx,Ny);
	fftwcomplex_2d_sym(vhat,Nx,Ny);
	fftwcomplex_2d_sym(Lphihat,Nx,Ny);
	fftwcomplex_2d_sym(M0uhat,Nx,Ny);
	fftwcomplex_2d_sym(M0vhat,Nx,Ny);	

	fftw_execute(plan_iftetahat);
    fftw_execute(plan_iftuhat);
    fftw_execute(plan_iftvhat);
    fftw_execute(plan_iftLphihat);
    fftw_execute(plan_iftM0uhat);
	fftw_execute(plan_iftM0vhat);

	fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftuhat);
    fftw_destroy_plan(plan_iftvhat);
    fftw_destroy_plan(plan_iftLphihat);
	fftw_destroy_plan(plan_iftM0uhat);
	fftw_destroy_plan(plan_iftM0vhat);
	
	fftw_complex* dampeta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* dampv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    fftw_complex* etaM0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaM0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta_vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaLphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* etaLphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* deletaH2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* deletaH2hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    
    // Create plans 
    fftw_plan plan_fftdampeta   = fftw_plan_dft_2d(Ny,Nx,dampeta,dampetahat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampu     = fftw_plan_dft_2d(Ny,Nx,dampu,dampuhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdampv     = fftw_plan_dft_2d(Ny,Nx,dampv,dampvhat, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0u    = fftw_plan_dft_2d(Ny,Nx,etaM0u,etaM0uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_u		= fftw_plan_dft_2d(Ny,Nx,eta_u,eta_uhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaM0v    = fftw_plan_dft_2d(Ny,Nx,etaM0v,etaM0vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_ffteta_v		= fftw_plan_dft_2d(Ny,Nx,eta_v,eta_vhat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftetaLphi   = fftw_plan_dft_2d(Ny,Nx,etaLphi,etaLphihat,FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan plan_fftdeletaH2	= fftw_plan_dft_2d(Ny,Nx,deletaH2,deletaH2hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dampeta[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*eta[i*Nx+j][0]/factIFFT;
			dampeta[i*Nx+j][1]=0;
			dampu[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*u[i*Nx+j][0]/factIFFT;
			dampu[i*Nx+j][1]=0;
			dampv[i*Nx+j][0]=7*Cpeak*dominfo.fbdy.charac[i][j]*v[i*Nx+j][0]/factIFFT;
			dampv[i*Nx+j][1]=0;
			
			etaM0u[i*Nx+j][0]=eta[i*Nx+j][0]*M0u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0u[i*Nx+j][1]=0;
			eta_u[i*Nx+j][0]=eta[i*Nx+j][0]*u[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_u[i*Nx+j][1]=0;
			
			etaM0v[i*Nx+j][0]=eta[i*Nx+j][0]*M0v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			etaM0v[i*Nx+j][1]=0;
			eta_v[i*Nx+j][0]=eta[i*Nx+j][0]*v[i*Nx+j][0]/pow(factIFFT,2)*dominfo.wallchar[i][j];
			eta_v[i*Nx+j][1]=0;
			
			//etaLphi[i*Nx+j][0]=eta[i*Nx+j][0]*Lphi[i*Nx+j][0]/pow(factIFFT,2);
			//etaLphi[i*Nx+j][1]=0;

			deletaH2[i*Nx+j][0]= 0.5*((pow(u[i*Nx+j][0],2)+pow(v[i*Nx+j][0],2))-pow(Lphi[i*Nx+j][0],2))/pow(factIFFT,2)
								*dominfo.wallchar[i][j];
			deletaH2[i*Nx+j][1]= 0;
		}
	}
	fftw_execute(plan_fftdampeta);
	fftw_execute(plan_fftdampu);
	fftw_execute(plan_fftdampv);
	fftw_execute(plan_fftetaM0u);
	fftw_execute(plan_ffteta_u);
	fftw_execute(plan_fftetaM0v);
	fftw_execute(plan_ffteta_v);
	fftw_execute(plan_fftetaLphi);
	fftw_execute(plan_fftdeletaH2);
	
	fftw_destroy_plan(plan_fftdampeta);
	fftw_destroy_plan(plan_fftdampu);
	fftw_destroy_plan(plan_fftdampv);
	fftw_destroy_plan(plan_fftetaM0u);
	fftw_destroy_plan(plan_ffteta_u);
	fftw_destroy_plan(plan_fftetaM0v);
	fftw_destroy_plan(plan_ffteta_v);
	fftw_destroy_plan(plan_fftetaLphi);
	fftw_destroy_plan(plan_fftdeletaH2);	
	
	fftw_complex* LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* deluH2_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0v_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
	fftw_complex* delvH2_hat 	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaM0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_plan plan_iftLetaM0u_hat  = fftw_plan_dft_2d(Ny,Nx, LetaM0u_hat, LetaM0u, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_iftLetaM0v_hat  = fftw_plan_dft_2d(Ny,Nx, LetaM0v_hat, LetaM0v, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* LetaLphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaLphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_iftLetaLphihat= fftw_plan_dft_2d(Ny,Nx, LetaLphihat, LetaLphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* divgradphihat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* divgradphi=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_plan plan_iftdivgradphihat= fftw_plan_dft_2d(Ny,Nx, divgradphihat, divgradphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* eta2M0u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta2M0uhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_ffteta2M0u     = fftw_plan_dft_2d(Ny,Nx,eta2M0u,eta2M0uhat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftw_complex* eta2M0v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* eta2M0vhat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_ffteta2M0v     = fftw_plan_dft_2d(Ny,Nx,eta2M0v,eta2M0vhat,FFTW_FORWARD, FFTW_ESTIMATE);		
	
	fftw_complex* dxu_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* dxu=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_iftdxu_hat= fftw_plan_dft_2d(Ny,Nx, dxu_hat, dxu, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* dyv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* dyv=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_iftdyv_hat= fftw_plan_dft_2d(Ny,Nx, dyv_hat, dyv, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_complex* ikxdeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH2_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etaM0u[i*Nx+j][0]=etaM0u[i*Nx+j][0]*factIFFT;//for operator
			etaM0v[i*Nx+j][0]=etaM0v[i*Nx+j][0]*factIFFT;//for operator	
		}
	}	

	//interpolation
	if (ninterp==2){
		HSSgen_2depth(etaM0u,etaM0uhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,LetaM0u_hat,Nx,Ny,factIFFT);
		HSSgen_2depth(etaM0v,etaM0vhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,LetaM0v_hat,Nx,Ny,factIFFT);
	}
	else {
		HSSgen_3depth(etaM0u,etaM0uhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,LetaM0u_hat,Nx,Ny,factIFFT);
		HSSgen_3depth(etaM0v,etaM0vhat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,LetaM0v_hat,Nx,Ny,factIFFT);
	}
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			ikxdeletaH2_hat[i*Nx+j][0] = -dominfo.kx[j]*deletaH2hat[i*Nx+j][1];
			ikxdeletaH2_hat[i*Nx+j][1] = dominfo.kx[j]*deletaH2hat[i*Nx+j][0];
			
			ikydeletaH2_hat[i*Nx+j][0] = -dominfo.ky[i]*deletaH2hat[i*Nx+j][1];
			ikydeletaH2_hat[i*Nx+j][1] = dominfo.ky[i]*deletaH2hat[i*Nx+j][0];
			
			LetaM0u_hat[i*Nx+j][0]=LetaM0u_hat[i*Nx+j][0]/grav;
			LetaM0u_hat[i*Nx+j][1]=LetaM0u_hat[i*Nx+j][1]/grav;			

			LetaM0v_hat[i*Nx+j][0]=LetaM0v_hat[i*Nx+j][0]/grav;
			LetaM0v_hat[i*Nx+j][1]=LetaM0v_hat[i*Nx+j][1]/grav;
			
			LetaLphihat[i*Nx+j][0]=LetaM0u_hat[i*Nx+j][0]+LetaM0v_hat[i*Nx+j][0];//Oprt.L2d[i][j]*etaLphihat[i*Nx+j][0];
			LetaLphihat[i*Nx+j][1]=LetaM0u_hat[i*Nx+j][1]+LetaM0v_hat[i*Nx+j][1];//Oprt.L2d[i][j]*etaLphihat[i*Nx+j][1];
			
			deluH2_hat[i*Nx+j][0] = -LetaM0u_hat[i*Nx+j][0]+dominfo.kx[j]*eta_uhat[i*Nx+j][1];
			deluH2_hat[i*Nx+j][1] = -LetaM0u_hat[i*Nx+j][1]-dominfo.kx[j]*eta_uhat[i*Nx+j][0];

			delvH2_hat[i*Nx+j][0] = -LetaM0v_hat[i*Nx+j][0]+dominfo.ky[i]*eta_vhat[i*Nx+j][1];
			delvH2_hat[i*Nx+j][1] = -LetaM0v_hat[i*Nx+j][1]-dominfo.ky[i]*eta_vhat[i*Nx+j][0];
			
			divgradphihat[i*Nx+j][0]=-dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j]-dominfo.ky[i]*Zhat[5*Nx*Ny+i*Nx+j];
			divgradphihat[i*Nx+j][1]=dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j]+dominfo.ky[i]*Zhat[4*Nx*Ny+i*Nx+j];
			
			eta2M0u[i*Nx+j][0]=pow(eta[i*Nx+j][0],2)*M0u[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
			eta2M0u[i*Nx+j][1]=0;
			
			eta2M0v[i*Nx+j][0]=pow(eta[i*Nx+j][0],2)*M0v[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
			eta2M0v[i*Nx+j][1]=0;
			
			dxu_hat[i*Nx+j][0] = -dominfo.kx[j]*Zhat[3*Nx*Ny+i*Nx+j];
			dxu_hat[i*Nx+j][1] = dominfo.kx[j]*Zhat[2*Nx*Ny+i*Nx+j];
			
			dyv_hat[i*Nx+j][0] = -dominfo.ky[i]*Zhat[5*Nx*Ny+i*Nx+j];
			dyv_hat[i*Nx+j][1] = dominfo.ky[i]*Zhat[4*Nx*Ny+i*Nx+j];
		}
	}
	fun_adjust(deluH2hat_adj,deluH2_hat,Nx,Ny);
	fun_adjust(delvH2hat_adj,delvH2_hat,Nx,Ny);
	fun_adjust(ikxdeletaH2_hatadj,ikxdeletaH2_hat,Nx,Ny);	
	fun_adjust(ikydeletaH2_hatadj,ikydeletaH2_hat,Nx,Ny);	
		
	fftwcomplex_2d_sym(LetaM0u_hat,Nx,Ny);
	fftwcomplex_2d_sym(LetaM0v_hat,Nx,Ny);
	fftwcomplex_2d_sym(divgradphihat,Nx,Ny);
	fftwcomplex_2d_sym(LetaLphihat,Nx,Ny);
	fftwcomplex_2d_sym(dxu_hat,Nx,Ny);
	fftwcomplex_2d_sym(dyv_hat,Nx,Ny);
	
	fftw_execute(plan_iftLetaM0u_hat);
	fftw_execute(plan_iftLetaM0v_hat);	
	fftw_execute(plan_iftdivgradphihat);
	fftw_execute(plan_ffteta2M0u);
	fftw_execute(plan_ffteta2M0v);
	fftw_execute(plan_iftLetaLphihat);
	fftw_execute(plan_iftdxu_hat);
	fftw_execute(plan_iftdyv_hat);
	
	fftw_destroy_plan(plan_iftLetaM0u_hat);
	fftw_destroy_plan(plan_iftLetaM0v_hat);	
	fftw_destroy_plan(plan_iftdivgradphihat);
	fftw_destroy_plan(plan_ffteta2M0u);
	fftw_destroy_plan(plan_ffteta2M0v);
	fftw_destroy_plan(plan_iftLetaLphihat);
	fftw_destroy_plan(plan_iftdxu_hat);
	fftw_destroy_plan(plan_iftdyv_hat);	

	fftw_complex* etaLetaM0u 	=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaLetaM0u_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaLetaM0v	=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etaLetaM0v_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_fftetaLetaM0u     = fftw_plan_dft_2d(Ny,Nx,etaLetaM0u,etaLetaM0u_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_plan plan_fftetaLetaM0v     = fftw_plan_dft_2d(Ny,Nx,etaLetaM0v,etaLetaM0v_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftw_complex* eta2dxu  	 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta2dxu_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta2dyv	 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta2dyv_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_ffteta2dxu= fftw_plan_dft_2d(Ny,Nx,eta2dxu,eta2dxu_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_plan plan_ffteta2dyv= fftw_plan_dft_2d(Ny,Nx,eta2dyv,eta2dyv_hat,FFTW_FORWARD, FFTW_ESTIMATE);	
	
	fftw_complex* deletaH3=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* deletaH3hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_plan plan_fftdeletaH3 = fftw_plan_dft_2d(Ny,Nx,deletaH3,deletaH3hat,FFTW_FORWARD, FFTW_ESTIMATE);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  
		etaLetaM0u[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*LetaM0u[i*Nx+j][0]/factIFFT;
		etaLetaM0u[i*Nx+j][1] = 0;
		
		etaLetaM0v[i*Nx+j][0] = eta[i*Nx+j][0]/factIFFT*LetaM0v[i*Nx+j][0]/factIFFT;
		etaLetaM0v[i*Nx+j][1] = 0;
		
		eta2dxu[i*Nx+j][0] = pow(eta[i*Nx+j][0],2)*dxu[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
		eta2dxu[i*Nx+j][1] = 0;
		
		eta2dyv[i*Nx+j][0] = pow(eta[i*Nx+j][0],2)*dyv[i*Nx+j][0]/pow(factIFFT,3)*dominfo.wallchar[i][j];
		eta2dyv[i*Nx+j][1] = 0;
        
		deletaH3[i*Nx+j][0]=(Lphi[i*Nx+j][0]/factIFFT*eta[i*Nx+j][0]/factIFFT*divgradphi[i*Nx+j][0]/factIFFT
								+Lphi[i*Nx+j][0]/factIFFT*LetaLphi[i*Nx+j][0]/factIFFT);
		deletaH3[i*Nx+j][1]=0;
	  }
	}	
	
	fftw_execute(plan_fftetaLetaM0u);
	fftw_execute(plan_fftetaLetaM0v);
	fftw_execute(plan_ffteta2dxu);
	fftw_execute(plan_ffteta2dyv);
	fftw_execute(plan_fftdeletaH3);
	
	fftw_destroy_plan(plan_fftetaLetaM0u);
	fftw_destroy_plan(plan_fftetaLetaM0v);
	fftw_destroy_plan(plan_ffteta2dxu);
	fftw_destroy_plan(plan_ffteta2dyv);
	fftw_destroy_plan(plan_fftdeletaH3);
	
	fftw_complex* deluH3temp_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* delvH3temp_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	//last round
	fftw_complex* ikxdeletaH3_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* ikydeletaH3_hat=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaLetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* LetaLetaM0v_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* Leta2dxu_hat    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    fftw_complex* Leta2dyv_hat    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	//interpolation
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etaLetaM0u[i*Nx+j][0]  = etaLetaM0u[i*Nx+j][0]*factIFFT;//for operator
			etaLetaM0v[i*Nx+j][0]  = etaLetaM0v[i*Nx+j][0]*factIFFT;//for operator
			eta2dxu[i*Nx+j][0] = eta2dxu[i*Nx+j][0]*factIFFT;//for operator
			eta2dyv[i*Nx+j][0] = eta2dyv[i*Nx+j][0]*factIFFT;//for operator
		}
	}	
    
	if (ninterp==2){
		HSSgen_2depth(etaLetaM0u,etaLetaM0u_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,LetaLetaM0u_hat,Nx,Ny,factIFFT);
		HSSgen_2depth(etaLetaM0v,etaLetaM0v_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,LetaLetaM0v_hat,Nx,Ny,factIFFT);
		HSSgen_2depth(eta2dxu,eta2dxu_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,Leta2dxu_hat,Nx,Ny,factIFFT);
		HSSgen_2depth(eta2dyv,eta2dyv_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_plus,Leta2dyv_hat,Nx,Ny,factIFFT);
	}
	else {
		HSSgen_3depth(etaLetaM0u,etaLetaM0u_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,LetaLetaM0u_hat,Nx,Ny,factIFFT);
		HSSgen_3depth(etaLetaM0v,etaLetaM0v_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,LetaLetaM0v_hat,Nx,Ny,factIFFT);
		HSSgen_3depth(eta2dxu,eta2dxu_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,Leta2dxu_hat,Nx,Ny,factIFFT);
		HSSgen_3depth(eta2dyv,eta2dyv_hat,Oprt.InterpD.Om2dSq_min,Oprt.InterpD.Om2dSq_mid,Oprt.InterpD.Om2dSq_plus,Leta2dyv_hat,Nx,Ny,factIFFT);
	}
		
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  ikxdeletaH3_hat[i*Nx+j][0] = -dominfo.kx[j]*deletaH3hat[i*Nx+j][1];
		  ikxdeletaH3_hat[i*Nx+j][1] = dominfo.kx[j]*deletaH3hat[i*Nx+j][0];
			
		  ikydeletaH3_hat[i*Nx+j][0] = -dominfo.ky[i]*deletaH3hat[i*Nx+j][1];
		  ikydeletaH3_hat[i*Nx+j][1] = dominfo.ky[i]*deletaH3hat[i*Nx+j][0];
		  		  	
		  deluH3temp_hat[i*Nx+j][0] = LetaLetaM0u_hat[i*Nx+j][0]/grav
										+0.5*Leta2dxu_hat[i*Nx+j][0]/grav										
										-0.5*(pow(dominfo.kx[j],2)*eta2M0uhat[i*Nx+j][0]);										
		  deluH3temp_hat[i*Nx+j][1] = LetaLetaM0u_hat[i*Nx+j][1]/grav		  
										+0.5*Leta2dxu_hat[i*Nx+j][1]/grav										
										-0.5*(pow(dominfo.kx[j],2)*eta2M0uhat[i*Nx+j][1]);
										
		  delvH3temp_hat[i*Nx+j][0] = LetaLetaM0v_hat[i*Nx+j][0]/grav		  
										+0.5*Leta2dyv_hat[i*Nx+j][0]/grav										
										-0.5*(pow(dominfo.ky[i],2)*eta2M0vhat[i*Nx+j][0]);
		  delvH3temp_hat[i*Nx+j][1] = LetaLetaM0v_hat[i*Nx+j][1]/grav		  
										+0.5*Leta2dyv_hat[i*Nx+j][1]/grav										
										-0.5*(pow(dominfo.ky[i],2)*eta2M0vhat[i*Nx+j][1]);
	  }
	}
	
	/*if (t>95){
				double* forline=new double[dominfo.Nx];
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=Lphi[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=LetaLphi[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				delete[] forline;
	}*/
	
	//adding nonlinear adjustment
	fun_adjust(deluH3hat_adj,deluH3temp_hat,Nx,Ny); fftw_free(deluH3temp_hat);	
	fun_adjust(delvH3hat_adj,delvH3temp_hat,Nx,Ny); fftw_free(delvH3temp_hat);
	fun_adjust(ikxdeletaH3_hatadj,ikxdeletaH3_hat,Nx,Ny);	
	fun_adjust(ikydeletaH3_hatadj,ikydeletaH3_hat,Nx,Ny);		
	
	fftw_free(ikxdeletaH2_hat);
	fftw_free(ikydeletaH2_hat);
	fftw_free(ikxdeletaH3_hat);
	fftw_free(ikydeletaH3_hat);
	fftw_free(etahat);
	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(Lphi);fftw_free(Lphihat);	
	fftw_free(dampeta);fftw_free(dampu);fftw_free(dampv);
	fftw_free(M0uhat);fftw_free(M0u);
	fftw_free(M0vhat);fftw_free(M0v);
	fftw_free(etaM0u);fftw_free(etaM0uhat);
	fftw_free(etaM0v);fftw_free(etaM0vhat);
	fftw_free(eta_u); fftw_free(eta_uhat); 
	fftw_free(eta_v); fftw_free(eta_vhat); 			
	fftw_free(LetaM0u_hat);
	fftw_free(LetaM0v_hat);	
	fftw_free(deluH2_hat);
	fftw_free(delvH2_hat);
	fftw_free(deletaH2);fftw_free(deletaH3);
	fftw_free(etaLphi);fftw_free(etaLphihat);
	fftw_free(divgradphihat);
	fftw_free(divgradphi);	
	fftw_free(LetaLphihat);fftw_free(LetaLphi);
	
	fftw_free(LetaM0u);
	fftw_free(LetaM0v);
	fftw_free(eta2M0u);
	fftw_free(eta2M0uhat);	
	fftw_free(eta2M0v);
	fftw_free(eta2M0vhat);
	fftw_free(dxu_hat);
	fftw_free(dxu);
	fftw_free(dyv_hat);
	fftw_free(dyv);
	fftw_free(etaLetaM0u);
	fftw_free(etaLetaM0u_hat);
	fftw_free(etaLetaM0v);
	fftw_free(etaLetaM0v_hat);
	fftw_free(eta2dxu);
	fftw_free(eta2dxu_hat);
	fftw_free(eta2dyv);
	fftw_free(eta2dyv_hat);
	fftw_free(deletaH2hat);
	fftw_free(deletaH3hat);
	
	fftw_free(LetaLetaM0u_hat);
	fftw_free(LetaLetaM0v_hat);
	fftw_free(Leta2dxu_hat);
	fftw_free(Leta2dyv_hat);
}


int rhs_AB1_2DB_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex *dampetahat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *deluH1hat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH1hat= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* etah=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* uh=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* vh=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fun_termshat_AB1_2DB_uv(etah,uh,vh,deluH1hat,delvH1hat,dampetahat,dampuhat,dampvhat,dominfo,Oprt,Zhat);	
	fun_source(t,Nx,Ny);

	//get friction term
    if (strcmp(friction,"yes")==0){
		fun_sfriction(t,Nx,Ny,factIFFT,etah,uh,vh);
	}

	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){					  
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re+couple_etahat[i][j].Re
										+deluH1hat[i*Nx+j][0]+delvH1hat[i*Nx+j][0]
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im+couple_etahat[i][j].Im
										+deluH1hat[i*Nx+j][1]+delvH1hat[i*Nx+j][1]
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
			
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (grav*dominfo.kx[j]*Zhat[Nx*Ny+i*Nx+j]+Sb_xhat[i][j].Re+Sf_xhat[i][j].Re+couple_uhat[i][j].Re
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-grav*dominfo.kx[j]*Zhat[i*Nx+j]+Sb_xhat[i][j].Im+Sf_xhat[i][j].Im+couple_uhat[i][j].Im
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (grav*dominfo.ky[i]*Zhat[Nx*Ny+i*Nx+j]+Sb_yhat[i][j].Re+Sf_yhat[i][j].Re+couple_vhat[i][j].Re
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-grav*dominfo.ky[i]*Zhat[i*Nx+j]+Sb_yhat[i][j].Im+Sf_yhat[i][j].Im+couple_vhat[i][j].Im
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];
		}
	}

	fftw_free(etah);fftw_free(uh);fftw_free(vh);
	fftw_free(dampetahat);
	fftw_free(dampuhat);
	fftw_free(dampvhat);
	fftw_free(deluH1hat);
	fftw_free(delvH1hat);
	
	return GSL_SUCCESS;
}

int rhs_AB1_2DF_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{
	int Nx=dominfo.Nx; 
	int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex *dampetahat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delphiH1hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* eta=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* u=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex* v=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	fun_delphiH1hat(t,Nx,Ny,Zhat,delphiH1hat);
	fun_termshat_AB1_2DF_uv(eta,u,v,dampetahat,dampuhat,dampvhat,dominfo,Oprt.Cpeak,Zhat);
	fun_source(t,Nx,Ny);

	//get friction term
    if (strcmp(friction,"yes")==0){
		fun_sfriction(t,Nx,Ny,factIFFT,eta,u,v);
	}
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re+couple_etahat[i][j].Re
										+delphiH1hat[i*Nx+j][0]
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im+couple_etahat[i][j].Im
										+delphiH1hat[i*Nx+j][1]
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (grav*dominfo.kx[j]*Zhat[Nx*Ny+i*Nx+j]+Sf_xhat[i][j].Re+couple_uhat[i][j].Re
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-grav*dominfo.kx[j]*Zhat[i*Nx+j]+Sf_xhat[i][j].Im+couple_uhat[i][j].Im
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (grav*dominfo.ky[i]*Zhat[Nx*Ny+i*Nx+j]+Sf_yhat[i][j].Re+couple_vhat[i][j].Re
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-grav*dominfo.ky[i]*Zhat[i*Nx+j]+Sf_yhat[i][j].Im+couple_vhat[i][j].Im
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];			
		}
	}

	fftw_free(eta);fftw_free(u);fftw_free(v);
	fftw_free(dampetahat);
	fftw_free(dampuhat);
	fftw_free(dampvhat);	
	fftw_free(delphiH1hat);
	
	return GSL_SUCCESS;
}

int rhs_AB2_2DF_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{ 
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex *dampetahat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikxdeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikydeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	fftw_complex *deluH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *delphiH1hat 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *eta        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *u        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *v        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	fun_delphiH1hat(t,Nx,Ny,Zhat,delphiH1hat);
	fun_termshat_AB2_2DF_uv(t,eta,u,v,ikxdeletaH2_hatadj,ikydeletaH2_hatadj,deluH2hat_adj,delvH2hat_adj,dampetahat,dampuhat,dampvhat,dominfo,Oprt.Cpeak,Zhat);
	fun_source(t,Nx,Ny);
	
	complex** delphiH_hat = declare_2darray_complex(Nx,Ny);
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH_hat[i][j].Re = delphiH1hat[i*Nx+j][0]+deluH2hat_adj[i*Nx+j][0]+delvH2hat_adj[i*Nx+j][0];										
			delphiH_hat[i][j].Im = delphiH1hat[i*Nx+j][1]+deluH2hat_adj[i*Nx+j][1]+delvH2hat_adj[i*Nx+j][1];										
		}	
	}
	
	
	//get friction term
    if (strcmp(friction,"yes")==0){
		fun_sfriction(t,Nx,Ny,factIFFT,eta,u,v);
	}

	//get B for breaking
	if (strcmp(breaking,"yes")==0){
		set_matrix_val(B,Nx,Ny,0.0);
		if (tprev<t){
			breaking_process(t,Zhat,eta,u,v,factIFFT);// filling the global matrix B
			fun_sbreaking(t,eta,delphiH_hat,Nx,Ny);
		}
	}		
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){			
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re+couple_etahat[i][j].Re
										+delphiH_hat[i][j].Re
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im+couple_etahat[i][j].Im
										+delphiH_hat[i][j].Im
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (grav*dominfo.kx[j]*Zhat[Nx*Ny+i*Nx+j]+Sb_xhat[i][j].Re+Sf_xhat[i][j].Re+couple_uhat[i][j].Re
										-ikxdeletaH2_hatadj[i*Nx+j][0]
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-grav*dominfo.kx[j]*Zhat[i*Nx+j]+Sb_xhat[i][j].Im+Sf_xhat[i][j].Im+couple_uhat[i][j].Im
										-ikxdeletaH2_hatadj[i*Nx+j][1]
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (grav*dominfo.ky[i]*Zhat[Nx*Ny+i*Nx+j]+Sb_yhat[i][j].Re+Sf_yhat[i][j].Re+couple_vhat[i][j].Re
										-ikydeletaH2_hatadj[i*Nx+j][0]
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-grav*dominfo.ky[i]*Zhat[i*Nx+j]+Sb_yhat[i][j].Im+Sf_yhat[i][j].Im+couple_vhat[i][j].Im
										-ikydeletaH2_hatadj[i*Nx+j][1]
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];			
		}
	}
	
	free_2darray_complex(delphiH_hat,Nx,Ny);
	fftw_free(eta);fftw_free(u);fftw_free(v);
	fftw_free(dampetahat);fftw_free(dampuhat);fftw_free(dampvhat);	
	fftw_free(ikxdeletaH2_hatadj);fftw_free(ikydeletaH2_hatadj);	
	fftw_free(deluH2hat_adj);
	fftw_free(delvH2hat_adj);
	fftw_free(delphiH1hat);
	  
	return GSL_SUCCESS;
}

int rhs_AB2_2DB_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{ 
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex *dampetahat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *ikxdeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikydeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));	
	fftw_complex *deluH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *deluH1hat 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH1hat 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *eta        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *u        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *v        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	fun_termshat_AB1_2DB_uv(eta,u,v,deluH1hat,delvH1hat,dampetahat,dampuhat,dampvhat,dominfo,Oprt,Zhat);	
	fun_termshat_AB2_2DB_uv(t,eta,u,v,deluH1hat,delvH1hat,ikxdeletaH2_hatadj,ikydeletaH2_hatadj,deluH2hat_adj,delvH2hat_adj,dominfo,Oprt.Cpeak,Zhat);
	fun_source(t,Nx,Ny);
	
	complex** delphiH_hat = declare_2darray_complex(Nx,Ny);
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH_hat[i][j].Re = deluH1hat[i*Nx+j][0]+delvH1hat[i*Nx+j][0]+deluH2hat_adj[i*Nx+j][0]+delvH2hat_adj[i*Nx+j][0];										
			delphiH_hat[i][j].Im = deluH1hat[i*Nx+j][1]+delvH1hat[i*Nx+j][1]+deluH2hat_adj[i*Nx+j][1]+delvH2hat_adj[i*Nx+j][1];										
		}	
	}	
	
	//get friction term
    if (strcmp(friction,"yes")==0){
		fun_sfriction(t,Nx,Ny,factIFFT,eta,u,v);
	}

	//get B for breaking
	if (strcmp(breaking,"yes")==0){
		set_matrix_val(B,Nx,Ny,0.0);
		if (tprev<t){
			breaking_process(t,Zhat,eta,u,v,factIFFT);// filling the global matrix B
			fun_sbreaking(t,eta,delphiH_hat,Nx,Ny);
		}
	}		
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){			
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re+couple_etahat[i][j].Re
										+delphiH_hat[i][j].Re
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im+couple_etahat[i][j].Im
										+delphiH_hat[i][j].Im
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (grav*dominfo.kx[j]*Zhat[Nx*Ny+i*Nx+j]+Sb_xhat[i][j].Re+Sf_xhat[i][j].Re+couple_uhat[i][j].Re
										-ikxdeletaH2_hatadj[i*Nx+j][0]
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-grav*dominfo.kx[j]*Zhat[i*Nx+j]+Sb_xhat[i][j].Im+Sf_xhat[i][j].Im+couple_uhat[i][j].Im
										-ikxdeletaH2_hatadj[i*Nx+j][1]
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (grav*dominfo.ky[i]*Zhat[Nx*Ny+i*Nx+j]+Sb_yhat[i][j].Re+Sf_yhat[i][j].Re+couple_vhat[i][j].Re
										-ikydeletaH2_hatadj[i*Nx+j][0]
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-grav*dominfo.ky[i]*Zhat[i*Nx+j]+Sb_yhat[i][j].Im+Sf_yhat[i][j].Im+couple_vhat[i][j].Im
										-ikydeletaH2_hatadj[i*Nx+j][1]
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];			
		}
	}
	
	free_2darray_complex(delphiH_hat,Nx,Ny);
	fftw_free(eta);fftw_free(u);fftw_free(v);
	fftw_free(dampetahat);fftw_free(dampuhat);fftw_free(dampvhat);	
	fftw_free(ikxdeletaH2_hatadj);fftw_free(ikydeletaH2_hatadj);
	fftw_free(deluH1hat);fftw_free(delvH1hat);
	fftw_free(deluH2hat_adj);
	fftw_free(delvH2hat_adj);
	 
	return GSL_SUCCESS;
}



int rhs_AB3_2DF_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{ 
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex *dampetahat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikxdeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikydeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikxdeletaH3_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikydeletaH3_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *deluH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *deluH3hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH3hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delphiH1hat 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *eta        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *u        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *v        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	fun_delphiH1hat(t,Nx,Ny,Zhat,delphiH1hat);
	fun_termshat_AB3_2DF_uv(t,eta,u,v,ikxdeletaH3_hatadj,ikydeletaH3_hatadj,deluH3hat_adj,delvH3hat_adj,ikxdeletaH2_hatadj,ikydeletaH2_hatadj,deluH2hat_adj,delvH2hat_adj,dampetahat,dampuhat,dampvhat,dominfo,Oprt.Cpeak,Zhat);
	fun_source(t,Nx,Ny);
	
	/*if (t>80){
				double* forline=new double[dominfo.Nx];
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=eta[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=u[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				delete[] forline;
	}*/
	
	complex** delphiH_hat = declare_2darray_complex(Nx,Ny);
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH_hat[i][j].Re = delphiH1hat[i*Nx+j][0]
										+deluH2hat_adj[i*Nx+j][0]+delvH2hat_adj[i*Nx+j][0]
										+deluH3hat_adj[i*Nx+j][0]+delvH3hat_adj[i*Nx+j][0];
			delphiH_hat[i][j].Im = delphiH1hat[i*Nx+j][1]
										+deluH2hat_adj[i*Nx+j][1]+delvH2hat_adj[i*Nx+j][1]
										+deluH3hat_adj[i*Nx+j][1]+delvH3hat_adj[i*Nx+j][1];
		}
	}

	//get friction term
    if (strcmp(friction,"yes")==0){
		fun_sfriction(t,Nx,Ny,factIFFT,eta,u,v);
	}

	//get B for breaking
	if (strcmp(breaking,"yes")==0){
		set_matrix_val(B,Nx,Ny,0.0);
		if (tprev<t){
			breaking_process(t,Zhat,eta,u,v,factIFFT);// filling the global matrix B
			fun_sbreaking(t,eta,delphiH_hat,Nx,Ny);
		}
	}		
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){			
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re+couple_etahat[i][j].Re
										+delphiH_hat[i][j].Re
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im+couple_etahat[i][j].Im
										+delphiH_hat[i][j].Im
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (grav*dominfo.kx[j]*Zhat[Nx*Ny+i*Nx+j]+Sb_xhat[i][j].Re+Sf_xhat[i][j].Re+couple_uhat[i][j].Re
										-ikxdeletaH2_hatadj[i*Nx+j][0]
										-ikxdeletaH3_hatadj[i*Nx+j][0]
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-grav*dominfo.kx[j]*Zhat[i*Nx+j]+Sb_xhat[i][j].Im+Sf_xhat[i][j].Im+couple_uhat[i][j].Im
										-ikxdeletaH2_hatadj[i*Nx+j][1]
										-ikxdeletaH3_hatadj[i*Nx+j][1]
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (grav*dominfo.ky[i]*Zhat[Nx*Ny+i*Nx+j]+Sb_yhat[i][j].Re+Sf_yhat[i][j].Re+couple_vhat[i][j].Re
										-ikydeletaH2_hatadj[i*Nx+j][0]
										-ikydeletaH3_hatadj[i*Nx+j][0]
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-grav*dominfo.ky[i]*Zhat[i*Nx+j]+Sb_yhat[i][j].Im+Sf_yhat[i][j].Im+couple_vhat[i][j].Im
										-ikydeletaH2_hatadj[i*Nx+j][1]
										-ikydeletaH3_hatadj[i*Nx+j][1]
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];			
		}
	}
	
	free_2darray_complex(delphiH_hat,Nx,Ny);
	fftw_free(eta);fftw_free(u);fftw_free(v);
	fftw_free(dampetahat);fftw_free(dampuhat);fftw_free(dampvhat);	
	fftw_free(ikxdeletaH2_hatadj);fftw_free(ikydeletaH2_hatadj);
	fftw_free(ikxdeletaH3_hatadj);fftw_free(ikydeletaH3_hatadj);
	fftw_free(deluH2hat_adj);fftw_free(deluH3hat_adj);
	fftw_free(delvH2hat_adj);fftw_free(delvH3hat_adj);
	fftw_free(delphiH1hat);
	  
	return GSL_SUCCESS;
}

int rhs_AB3_2DB_uv(double t, const double* Zhat, double* dt_Zhat,void *params)
{ 
	int Nx=dominfo.Nx; int Ny=dominfo.Ny;
	int factIFFT=Nx*Ny;
	
	fftw_complex *dampetahat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampuhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *dampvhat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikxdeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikydeletaH2_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikxdeletaH3_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *ikydeletaH3_hatadj    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *deluH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *deluH3hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH2hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH3hat_adj      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	
	fftw_complex *deluH1hat 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *delvH1hat 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *eta        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *u        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	fftw_complex *v        	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	fun_termshat_AB1_2DB_uv(eta,u,v,deluH1hat,delvH1hat,dampetahat,dampuhat,dampvhat,dominfo,Oprt,Zhat);	
	fun_termshat_AB3_2DB_uv(t,eta,u,v,deluH1hat,delvH1hat,ikxdeletaH3_hatadj,ikydeletaH3_hatadj,deluH3hat_adj,delvH3hat_adj,
	ikxdeletaH2_hatadj,ikydeletaH2_hatadj,deluH2hat_adj,delvH2hat_adj,dampetahat,dampuhat,dampvhat,dominfo,Oprt.Cpeak,Zhat);
	fun_source(t,Nx,Ny);
	
	/*if (t>80){
				double* forline=new double[dominfo.Nx];
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=eta[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				for (int i=0;i<dominfo.Nx;i++){
					forline[i]=u[i][0]/(Nx*Ny);
				}
				fun_plot_line(dominfo.x,forline,dominfo.Nx);
				delete[] forline;
	}*/
	
	complex** delphiH_hat = declare_2darray_complex(Nx,Ny);
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			delphiH_hat[i][j].Re = deluH1hat[i*Nx+j][0]+delvH1hat[i*Nx+j][0]
										+deluH2hat_adj[i*Nx+j][0]+delvH2hat_adj[i*Nx+j][0]
										+deluH3hat_adj[i*Nx+j][0]+delvH3hat_adj[i*Nx+j][0];
			delphiH_hat[i][j].Im = deluH1hat[i*Nx+j][1]+delvH1hat[i*Nx+j][1]
										+deluH2hat_adj[i*Nx+j][1]+delvH2hat_adj[i*Nx+j][1]
										+deluH3hat_adj[i*Nx+j][1]+delvH3hat_adj[i*Nx+j][1];
		}
	}

	//get friction term
    if (strcmp(friction,"yes")==0){
		fun_sfriction(t,Nx,Ny,factIFFT,eta,u,v);
	}

	//get B for breaking
	if (strcmp(breaking,"yes")==0){
		set_matrix_val(B,Nx,Ny,0.0);
		if (tprev<t){
			breaking_process(t,Zhat,eta,u,v,factIFFT);// filling the global matrix B
			fun_sbreaking(t,eta,delphiH_hat,Nx,Ny);
		}
	}		
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){			
			dt_Zhat[i*Nx+j]  	   	= (Source_hat[i][j].Re+couple_etahat[i][j].Re
										+delphiH_hat[i][j].Re
										-dampetahat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[Nx*Ny+i*Nx+j]  	= (Source_hat[i][j].Im+couple_etahat[i][j].Im
										+delphiH_hat[i][j].Im
										-dampetahat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[2*Nx*Ny+i*Nx+j]	= (grav*dominfo.kx[j]*Zhat[Nx*Ny+i*Nx+j]+Sb_xhat[i][j].Re+Sf_xhat[i][j].Re+couple_uhat[i][j].Re
										-ikxdeletaH2_hatadj[i*Nx+j][0]
										-ikxdeletaH3_hatadj[i*Nx+j][0]
										-dampuhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[3*Nx*Ny+i*Nx+j]	= (-grav*dominfo.kx[j]*Zhat[i*Nx+j]+Sb_xhat[i][j].Im+Sf_xhat[i][j].Im+couple_uhat[i][j].Im
										-ikxdeletaH2_hatadj[i*Nx+j][1]
										-ikxdeletaH3_hatadj[i*Nx+j][1]
										-dampuhat[i*Nx+j][1])*dominfo.aal[i][j];
										
			dt_Zhat[4*Nx*Ny+i*Nx+j]	= (grav*dominfo.ky[i]*Zhat[Nx*Ny+i*Nx+j]+Sb_yhat[i][j].Re+Sf_yhat[i][j].Re+couple_vhat[i][j].Re
										-ikydeletaH2_hatadj[i*Nx+j][0]
										-ikydeletaH3_hatadj[i*Nx+j][0]
										-dampvhat[i*Nx+j][0])*dominfo.aal[i][j];
			dt_Zhat[5*Nx*Ny+i*Nx+j]	= (-grav*dominfo.ky[i]*Zhat[i*Nx+j]+Sb_yhat[i][j].Im+Sf_yhat[i][j].Im+couple_vhat[i][j].Im
										-ikydeletaH2_hatadj[i*Nx+j][1]
										-ikydeletaH3_hatadj[i*Nx+j][1]
										-dampvhat[i*Nx+j][1])*dominfo.aal[i][j];			
		}
	}
	
	free_2darray_complex(delphiH_hat,Nx,Ny);
	fftw_free(eta);fftw_free(u);fftw_free(v);
	fftw_free(dampetahat);fftw_free(dampuhat);fftw_free(dampvhat);	
	fftw_free(ikxdeletaH2_hatadj);fftw_free(ikydeletaH2_hatadj);
	fftw_free(ikxdeletaH3_hatadj);fftw_free(ikydeletaH3_hatadj);
	fftw_free(deluH2hat_adj);fftw_free(deluH3hat_adj);
	fftw_free(delvH2hat_adj);fftw_free(delvH3hat_adj);
	fftw_free(deluH1hat);fftw_free(delvH1hat);
	  
	return GSL_SUCCESS;
}

