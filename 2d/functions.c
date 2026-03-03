#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))

double phase_velocity_i(double k_i,double d_i)
{
	double Up_i;
	
	if(k_i==0)Up_i=sqrt(grav*d_i);
	else Up_i=(k_i<0?-1:1)*sqrt(grav*k_i*tanh(d_i*k_i))/k_i;
	
	return Up_i;
}

double** funOprt_Ug2d(double* k,double* kt,double depth,double ej[2],int N,int Nt)
{
	double K;
	double** Ug2d = declare_2darray(N,Nt);
	
	for (int i=0;i<N;i++){
		for (int j=0;j<Nt;j++) {
			K = sqrt(pow(k[i],2)+pow(kt[j],2));
			if (K==0){
				Ug2d[j][i] = sqrt(grav*depth);
			}
			else {
				Ug2d[j][i] = fun_sign(k[i]*ej[0]+kt[j]*ej[1])*sqrt(grav)/2./pow((K*tanh(depth*K)),(0.5))*(tanh(depth*K)+K*(1-pow(tanh(depth*K),2))*depth);	
			}
		}
	}	
	return Ug2d;
}

double Up_exact(double k_i,double d_i,double g)
{
	double Up_i;
	
	if (k_i==0) {
		Up_i=sqrt(g*d_i);
	}
	else {
		Up_i=(k_i<0?-1:1)*sqrt(g*k_i*tanh(d_i*k_i))/k_i;
	}	
	
	return Up_i;
}

double** funOprt_Ug2d_i(double** Kx,double depth,int N,int Nt)
{
	double** Ug2d = declare_2darray(N,Nt);
	double K;
	for (int i=0;i<N;i++){
		for (int j=0;j<Nt;j++) {
			K = Kx[j][i];
			if (K==0){
				Ug2d[j][i] = sqrt(grav*depth);
			}
			else {
				Ug2d[j][i] = fun_sign(K)*sqrt(grav)/2./pow((K*tanh(depth*K)),(0.5))*(tanh(depth*K)+K*(1-pow(tanh(depth*K),2))*depth);	
			}
		}
	}	
	return Ug2d;
}

double fun_exact_disp1d_val(double k, double D,double U)
{
	double Om=fun_sign(k)*sqrt(grav*k*tanh(k*D))+k*U;
	return Om ;
}

double fun_exact_Cp1d_val(double k,double D,double U)
{
	double Cp;
	if (k==0){
		Cp=sqrt(grav*D)+U;
	}
	else{
		Cp=(fun_sign(k)*sqrt(grav*k*tanh(k*D))+k*U)/k;
	}
	return Cp;
}

double fun_exact_Ug_val(double k_i,double d_i,double U)
{
	double Ug_i;	
	if(k_i==0) {
		Ug_i = sqrt(grav*d_i)+U;
	}
	else {
		Ug_i=((k_i<0?-1:1)*sqrt(grav)/2/sqrt(k_i*tanh(d_i*k_i))*(tanh(d_i*k_i)+k_i*(1-pow(tanh(d_i*k_i),2))*d_i ))+U;
	}	
	return Ug_i;
}

void fun_exact_Ug(double* Ug,double* k_i,double d_i,int N,double U)
{
	for (int j=0;j<N;j++){
		Ug[j] = fun_exact_Ug_val(k_i[j],fabs(d_i),U);
	}
}

double fun_exact_Cp_val(double kx, double ky, double depth, Vec2d u)
{
	double normK=sqrt(pow(kx,2)+pow(ky,2));
	double Cp;
	if (normK==0){
		Cp=sqrt(grav*depth);
	}
	else{
		Cp=sqrt(grav*tanh(normK*depth)/normK)+(kx*u.x+ky*u.y)/normK;
	}
	return Cp;
}

double fun_exact_disp_val(double kx, double ky, double depth, double thetadir, double ux, double uy)
{
	double signO=fun_sign((double) (kx*cos(thetadir)+ky*sin(thetadir)));
	double normK=sqrt(pow(kx,2)+pow(ky,2));
	double Om=signO*sqrt(grav*normK*tanh(normK*depth))+kx*ux+ky*uy;
	return Om;
}

double* fun_exact_disp1d(double* k, double D,int N,double U)
{
	double* Om=new double[N];	
	for(int i=0;i<N;i++){
		Om[i]=fun_sign(k[i])*sqrt(grav*k[i]*tanh(k[i]*D))+k[i]*U;
	}
	return Om;
}

void read_input(FILE* fparams)
{
	//initialization reading data
   	printf("Reading input parameters and files...\n");
   	struct stat st = {0};
    if (stat(arg2, &st) == -1) {
		mkdir(arg2, 0700);
	}
	
	//reading parameters line by line
	char line[256];
	double mdir;
	fgets(line, sizeof(line), fparams);	
	fscanf(fparams,"%*s %*s %s \n",wavename);
	fscanf(fparams,"%*s %*s %lf %lf \n",&ramp,&taper);
	fscanf(fparams,"%*s %*s %s \n",initial);
	fgets(line, sizeof(line), fparams);
	fscanf(fparams,"%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf\n",&dom_t.start,&dom_t.end,&dom_t.delta);
	fscanf(fparams,"%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf \n %*s %*s %s\n",&seainfo.Tp,&seainfo.Hs,&seainfo.gamJS,seainfo.s);
	fscanf(fparams,"%*s %*s %lf %lf \n %*s %*s %lf\n",&seainfo.current.x,&seainfo.current.y,&mdir);
	fscanf(fparams,"%*s %*s %s %s\n",orientation,propagation);
	fgets(line, sizeof(line), fparams);
	fscanf(fparams,"%*s %*s %lf %lf\n %*s %*s %lf %lf\n",&dom_x.start,&dom_y.start,&dom_x.end,&dom_y.end);
	fscanf(fparams,"%*s %*s %lf %lf %lf %lf\n",&lineinflux.xstart,&lineinflux.ystart,&lineinflux.xend,&lineinflux.yend);
	fscanf(fparams,"%*s %*s %lf %lf\n %*s %*s %lf %lf %lf %lf\n ",&dom_x.delta,&dom_y.delta,&dominfo.fbdy.Lleft,
			&dominfo.fbdy.Lright,&dominfo.fbdy.Lbottom,&dominfo.fbdy.Ltop);
	fscanf(fparams,"%*s %*s %s \n %*s %*s %lf \n %*s %*s %lf \n",dominfo.bathy.Id,&dominfo.bathy.depth,&Dmin_input); 
	fscanf(fparams,"%*s %*s %lf \n",&dominfo.Nadj);
	fgets(line, sizeof(line), fparams);	
	fscanf(fparams,"%*s %*s %s \n",evolinfo.model);
	fscanf(fparams,"%*s %*s %s \n",breaking);
	fscanf(fparams,"%*s %*s %s \n",friction);	
	fscanf(fparams,"%*s %*s %s %s \n",kinematic,kinematic_type);
	fscanf(fparams,"%*s %*s %d \n",&ngauge);
	fscanf(fparams,"%*s %*s %d \n",&npartition);
	fgets(line, sizeof(line), fparams);
	fscanf(fparams,"%*s %*s %lf \n %*s %*s %lf %lf \n",&pb.KBC,&pb.TC,&pb.Tchar);
	fgets(line, sizeof(line), fparams);
	fscanf(fparams,"%*s %*s %lf %lf\n %*s %*s %lf %lf\n",&xinterv1,&yinterv1,&xinterv2,&yinterv2);
	fscanf(fparams,"%*s %*s %lf %lf\n",&dx_rec,&dy_rec);	
	fscanf(fparams,"%*s %*s %d \n %*s %*s %d \n",&nsurf,&ndeep);
	fgets(line, sizeof(line), fparams);
	fscanf(fparams,"%*s %*s %s \n",cut_temp);
	fscanf(fparams,"%*s %*s %d %s\n",&ninterp,mid_temp);
	fscanf(fparams,"%*s %*s %s \n",min_z);
	fscanf(fparams,"%*s %*s %s \n",max_z);
	fscanf(fparams,"%*s %*s %s ",wall);
	if ((strcmp(wall,"no")!=0)||(strcmp(wall,"user_defined")!=0)){
		nwall = atoi(wall);
		Xwall[0] = new double[nwall];
		Xwall[1] = new double[nwall];
		Ywall[0] = new double[nwall];
		Ywall[1] = new double[nwall];
		for (int j=0;j<nwall;j++){
			fscanf(fparams,"%lf %lf %lf %lf",&Xwall[0][j],&Ywall[0][j],&Xwall[1][j],&Ywall[1][j]);
		}
	}
	fscanf(fparams,"\n");
	fscanf(fparams,"%*s %*s %d %lf %s\n",&walltype,&refl_coef,wallmethod);
	fgets(line, sizeof(line), fparams);
	fscanf(fparams,"%*s %*s %s \n",coupling);
	double OZin,OZiny,OZtr,OZtry;
	fscanf(fparams,"%*s %*s %lf %lf\n",&ozstart[0],&ozstart[1]);
	fscanf(fparams,"%*s %*s %lf %lf\n",&ozend[0],&ozend[1]);	
	fscanf(fparams,"%*s %*s %lf %lf\n",&OZin,&OZiny);
	fscanf(fparams,"%*s %*s %lf %lf\n",&OZtr,&OZtry);
	fscanf(fparams,"%*s %*s %s \n",dir_force);
	
	lenOZin  = OZin;
	lenOZiny = OZiny;
	lenOZtr  = OZtr;
	lenOZtry = OZtry;
	
	pb.delb = 1.5;
	seainfo.mdir = mdir*M_PI/180.0;
	dominfo.U 	 = sqrt(pow(seainfo.current.x,2)+pow(seainfo.current.y,2));
	strcpy(bath,dominfo.bathy.Id);
	
	if (strcmp(cut_temp,"default")==0){
		if (strcmp(evolinfo.model,"HS1")==0){
			evolinfo.cutfrac = 2;
		}
		else if (strcmp(evolinfo.model,"HS2")==0){
			evolinfo.cutfrac = 4;
		}
		else if (strcmp(evolinfo.model,"HS3")==0){
			evolinfo.cutfrac = 6;
		}		
	} 
	else {
		evolinfo.cutfrac = atoi(cut_temp);
	}
}

double fun_cutoff(double R)
{
	double a,b,c,coff;
	if (strcmp(evolinfo.model,"HS1")==0){
		a=0.8086;
		b=0.2413;
		c=0.0044;
		if (R==1){
			coff=1;
		}
		else if(R==0){
			coff=0;
		}
		else {
			coff=(-b+sqrt(pow(b,2)-4*a*(c-R))/(2*a));
		}
	}
	else {
		a=0.7289;
		b=0.3229;
		c=-0.00058;
		if(R==1){
			coff=1;
		}
		else if(R==0){
			coff=0;
		}
		else {
			coff=(-b+sqrt(pow(b,2)-4*a*(c-R))/(2*a));
		}
	}
	
	return coff;
}

void fun_wallchar(double* x, double* y,int Nx,int Ny,double dx,double dy)
{
	cfSAwall = declare_2darray(Nx,Ny);
	set_matrix_val(cfSAwall,Nx,Ny,1);
	if (strcmp(wall,"user_defined")==0){
		double cutoff = fun_cutoff(refl_coef);
		char* wallinp 	 = new char[256];
		sprintf(wallinp,"%s/wall.dat",arg1);
		FILE* fwall = fopen(wallinp,"r");
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (!fscanf(fwall,"%lf ", &dominfo.wallchar[i][j])){
					printf("The data size of wall.dat is incorrect!");
					break;
				}
				cfSAwall[i][j] = 1-dominfo.wallchar[i][j];
				if (cfSAwall[i][j]<1){
					cfSAwall[i][j] = 0;//directly damp
				}
				dominfo.wallchar[i][j] = 1-dominfo.wallchar[i][j]*cutoff;
			}
		}
		fclose(fwall);
		
	}
	else if (strcmp(wall,"no")!=0){
		double** Heav  = declare_2darray(Nx,Ny);
		int* iwx[2];
		int* iwy[2];
		iwx[0] = new int[nwall];
		iwx[1] = new int[nwall];
		iwy[0] = new int[nwall];
		iwy[1] = new int[nwall];		
		
		set_matrix_val(Heav,Nx,Ny,1);
		for (int iw=0;iw<nwall;iw++){			
			iwx[0][iw] = round((Xwall[0][iw]-x[0])/dx);//position wall
			iwx[1][iw] = round((Xwall[1][iw]-x[0])/dx);
			iwy[0][iw] = round((Ywall[0][iw]-y[0])/dy);
			iwy[1][iw] = round((Ywall[1][iw]-y[0])/dy);
			
			double** heavx = fun_Heav2d(x,Nx,Xwall[0][iw]-x[0],x[Nx-1]-Xwall[1][iw],Ny);
			double** heavy = fun_Heav2d(y,Ny,Ywall[0][iw]-y[0],y[Ny-1]-Ywall[1][iw],Nx);			
			for(int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					Heav[i][j] = Heav[i][j]*(1-heavx[i][j]*heavy[j][i]);
				}
			}
			free_2darray(heavx,Nx,Ny);	
			free_2darray(heavy,Ny,Nx);
		}				
		
		set_matrix_val(dominfo.wallchar,Nx,Ny,1);
		if (strcmp(wallmethod,"influxing")!=0){//this if for energy truncation
			double cutoff = fun_cutoff(refl_coef);
			for(int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					if (Heav[i][j]==0){						
						dominfo.wallchar[i][j] = 1-cutoff;
					}
				}
			}			
		}
		
		//damping inside wall
		int Nxw = iwx[1][0]-iwx[0][0]+1;
		int Nyw = iwy[1][0]-iwy[0][0]+1;
		
		double* xw = new double[Nxw];		
		double* yw = new double[Nyw];
		for (int i=0;i<Nxw;i++){
			xw[i] = x[(iwx[0][0])]+i*dx;
		}
		for (int i=0;i<Nyw;i++){
			yw[i] = y[(iwy[0][0])]+i*dy;
		}
		
		for(int i=iwy[0][0];i<=iwy[1][0];i++){
			for(int j=iwx[0][0];j<=iwx[1][0];j++){
				cfSAwall[i][j] = 0;
			}
		}
			
		delete[] xw;
		delete[] yw;

		//define index wall
		idxWall = new Vec2d[Nx*Ny];
		int id  = 0;
		for (int iw=0;iw<nwall;iw++){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					if (((i==iwx[0][iw])||(i==iwx[1][iw]))&&((j>=iwy[0][iw])&&(j<=iwy[1][iw]))){
						idxWall[id].x = i;
						idxWall[id].y = j;
						id++;
					}
				}
			}
			for(int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					if (((i==iwy[0][iw])||(i==iwy[1][iw]))&&((j>=iwx[0][iw])&&(j<=iwx[1][iw]))){
						idxWall[id].x = j;
						idxWall[id].y = i;
						id++;
					}
				}
			}
		}
		n_idxwall = id;				

		//define directory boundary
		dirWall = new int[n_idxwall];
		for (int i=0;i<n_idxwall;i++){
			for (int iw=0;iw<nwall;iw++){
				if (idxWall[i].x==iwx[0][iw]){
					dirWall[i] = -1;
				}
				else if (idxWall[i].x==iwx[1][iw]){
					dirWall[i] = 1;
				}
				else if (idxWall[i].y==iwy[0][iw]){
					dirWall[i] = -1;
				}
				else if (idxWall[i].y==iwy[1][iw]){
					dirWall[i] = 1;
				}
			}
		}
	}
}

void fun_gradient(double** grad_wallchar_x,double** grad_wallchar_y,double** wallchar,
int Nx,int Ny,double dx,double dy)
{
	double** Fgrad_wallchar_x=declare_2darray(Nx,Ny);
	double** Fgrad_wallchar_y=declare_2darray(Nx,Ny);
	double** Bgrad_wallchar_x=declare_2darray(Nx,Ny);
	double** Bgrad_wallchar_y=declare_2darray(Nx,Ny);
	
	//forward difference
	#pragma omp parallel for
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx-1;j++){
			Fgrad_wallchar_x[i][j]= (wallchar[i][j+1]-wallchar[i][j])/dx;
		}
		Fgrad_wallchar_x[i][Nx-1] = (wallchar[i][0]-wallchar[i][Nx-1])/dx;
	}
	#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny-1;j++){
			Fgrad_wallchar_y[j][i]= (wallchar[j+1][i]-wallchar[j][i])/dy;
		}
		Fgrad_wallchar_y[Ny-1][i] = (wallchar[0][i]-wallchar[Ny-1][i])/dy;
	}
	
	//backward difference
	#pragma omp parallel for
	for(int i=0;i<Ny;i++){
		for(int j=1;j<Nx;j++){
			Bgrad_wallchar_x[i][j]= (wallchar[i][j]-wallchar[i][j-1])/dx;
		}
		Bgrad_wallchar_x[i][0] = (wallchar[i][0]-wallchar[i][Nx-1])/dx;
	}
	#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		for(int j=1;j<Ny;j++){
			Bgrad_wallchar_y[j][i]= (wallchar[j][i]-wallchar[j-1][i])/dy;
		}
		Bgrad_wallchar_y[0][i] = (wallchar[0][i]-wallchar[Ny-1][i])/dy;
	}
	
	#pragma omp parallel for
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (Fgrad_wallchar_x[i][j]>0){
				Fgrad_wallchar_x[i][j] = 0;
			}
			if (Fgrad_wallchar_y[i][j]>0){
				Fgrad_wallchar_y[i][j] = 0;
			}
			if (Bgrad_wallchar_x[i][j]<0){
				Bgrad_wallchar_x[i][j] = 0;
			}
			if (Bgrad_wallchar_y[i][j]<0){
				Bgrad_wallchar_y[i][j] = 0;
			}
		}
	}
	
	#pragma omp parallel for
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			grad_wallchar_x[i][j] = Fgrad_wallchar_x[i][j]+Bgrad_wallchar_x[i][j];
			grad_wallchar_y[i][j] = Fgrad_wallchar_y[i][j]+Bgrad_wallchar_y[i][j];
		}
	}
	free_2darray(Fgrad_wallchar_x,Nx,Ny);
	free_2darray(Fgrad_wallchar_y,Nx,Ny);
	free_2darray(Bgrad_wallchar_x,Nx,Ny);
	free_2darray(Bgrad_wallchar_y,Nx,Ny);
}

void fun_wall_influxchar(double** wall_gam,double** wall_skewgam,double* x,double* y,double dx,double dy,int Nx,int Ny)
{/*
	set_matrix_val(wall_gam,Nx,Ny,0.0);
    set_matrix_val(wall_skewgam,Nx,Ny,0.0);
	if (strcmp(wall,"rectangle")==0){		
		int iwx[2],iwy[2];
		iwx[0] = round((Xwall[0]-x[0])/dx);
		iwx[1] = round((Xwall[1]-x[0])/dx);
		iwy[0] = round((Ywall[0]-y[0])/dy);
		iwy[1] = round((Ywall[1]-y[0])/dy);
		double dist,om;
		double fact = pow((2*M_PI),2)/(dx*dy);		
        
        //for wall vertical boundary 
        double*  Ugx	 = new double[Nx];
		double*  Ugx_hat = new double[Nx];
        complex* temp	 = new complex[Nx];
        complex* temp2   = new complex[Nx];
        double** dxGam   = declare_2darray(Nx,Ny);     
        dxGam1  = declare_2darray(Nx,Ny);
        dxGam2  = declare_2darray(Nx,Ny); 
        
        set_matrix_val(dxGam,Nx,Ny,0.0);
        set_matrix_val(dxGam1,Nx,Ny,0.0);
        set_matrix_val(dxGam2,Nx,Ny,0.0);
        for (int xx=0;xx<2;xx++){
			if (xx==0){
				dist = x[iwx[xx]]-x[0];
			}
			else{
				dist = x[iwx[xx]]-x[0];
			}
			if ((x[iwx[xx]]==x[0])||(x[iwx[xx]]==x[Nx-1])){break;}
			
			for(int j=0;j<Ny;j++){	
				if ((j>=iwy[0])&&(j<=iwy[1])){		
					fun_exact_Ug(Ugx_hat,dominfo.kx,dominfo.bathy.data[j][(iwx[xx])],Nx,dominfo.U);
					ifft_1d_real_real(Ugx,Ugx_hat,Nx);//size Nx
					circshift(Ugx,round(dist/dx),Nx);
					fft_1d(temp,Ugx,Nx);
					for (int i=0;i<Nx;i++){
						om = omega_i(dominfo.kx[i],-dominfo.bathy.data[j][i]);						
						wall_gam[j][i]  += fact*Ugx[i];
						temp2[i].Re = -fact*om*temp[i].Im;
						temp2[i].Im = fact*om*temp[i].Re;
					}		
					ifft_1d_com_re(temp2,dxGam[j],Nx);//size Nx tempSkew
					if (xx==0){
						for (int i=0;i<Nx;i++){
							dxGam1[j][i] += dxGam[j][i];
						} 
					}
					else {
						for (int i=0;i<Nx;i++){
							dxGam2[j][i] += -dxGam[j][i];
						}
					}					
				}     
			}
        }
        delete[] Ugx;
        delete[] Ugx_hat;
        delete[] temp;
        delete[] temp2;       
        free_2darray(dxGam,Nx,Ny);  
		
		//for horizontal wall boundary
		double*  Ugy	 = new double[Ny];
		double*  Ugy_hat = new double[Ny];
        complex* tempy   = new complex[Ny];
        complex* tempy2  = new complex[Ny];
        double** dyGam    = declare_2darray(Ny,Nx);        
        dyGam1   = declare_2darray(Ny,Nx);
        dyGam2   = declare_2darray(Ny,Nx);        
        set_matrix_val(dyGam,Ny,Nx,0.0);        
        set_matrix_val(dyGam1,Ny,Nx,0.0);
        set_matrix_val(dyGam2,Ny,Nx,0.0);
		for (int yy=0;yy<2;yy++){
			if (yy==0){
				dist = y[iwy[yy]]-y[0];
			}
			else{
				dist = y[iwy[yy]]-y[0];
			}
			if ((y[iwy[yy]]==y[0])||(y[iwy[yy]]==y[Ny-1])){break;}
			
			for(int j=0;j<Nx;j++){	
				if ((j>=iwx[0])&&(j<=iwx[1])){		
					fun_exact_Ug(Ugy_hat,dominfo.ky,dominfo.bathy.data[(iwy[yy])][j],Ny,dominfo.U);
					ifft_1d_real_real(Ugy,Ugy_hat,Ny);//size Ny
					circshift(Ugy,round(dist/dy),Ny);
					fft_1d(tempy,Ugy,Ny);
					for (int i=0;i<Ny;i++){
						om = omega_i(dominfo.ky[i],-dominfo.bathy.data[i][j]);
						wall_gam[i][j]  += fact*Ugy[i];
						tempy2[i].Re = -fact*om*tempy[i].Im;
						tempy2[i].Im = fact*om*tempy[i].Re;
					}	
					ifft_1d_com_re(tempy2,dyGam[j],Ny);//size Ny   
					if (yy==0){
						for (int i=0;i<Ny;i++){
							dyGam1[j][i] += dyGam[j][i];
						} 
					}
					else {
						for (int i=0;i<Ny;i++){
							dyGam2[j][i] += -dyGam[j][i];
						}
					}	
				}     
			}
        }
        delete[] Ugy;
        delete[] Ugy_hat;
        delete[] tempy;
        delete[] tempy2;
        free_2darray(dyGam,Ny,Nx);
        
        //Gam*wallchar this is for wallchar in freq. dependent
        for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				wall_gam[i][j]    = wall_gam[i][j];//*dominfo.wallchar[i][j];
				wall_skewgam[i][j]= (dxGam1[i][j]+dxGam2[i][j]+dyGam1[j][i]+dyGam2[j][i]);//*dominfo.wallchar[i][j];
			}
		}
	}	*/
}

void fun_friction(double** coef,int Nx,int Ny)
{
	//initiation
	Sf_xhat = declare_2darray_complex(Nx,Ny);
	Sf_yhat = declare_2darray_complex(Nx,Ny);
	set_matrix_complex_val(Sf_xhat,Nx,Ny,0.0);
	set_matrix_complex_val(Sf_yhat,Nx,Ny,0.0);
		
	fric_depth= declare_2darray(Nx,Ny);
	set_matrix_val(coef,Nx,Ny,0.0);

	if (runup_id==1){
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				fric_depth[i][j]=-dominfo.bathy.min[i][j];
			}
		}
	}
	else{
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				fric_depth[i][j]=-dominfo.bathy.data[i][j];
			}
		}
	}

	double coef_temp;

	if (runup_id==1){
		coef_temp=pow(0.02,2)*grav;
		set_matrix_val(coef,Nx,Ny,coef_temp);
	}
	
	if(strcmp(friction,"yes")==0){
		FILE* frinput;
		char* frfile = new char[256];
		sprintf(frfile,"%s/cf.dat",arg1);
		frinput = fopen(frfile, "r");
		if(frinput==NULL){
			printf("Can not open %s!!\n",frfile);
			exit(0);
		} 
		else {
			for(int i=0;i<Ny;i++){
				for(int j=0;j<Nx;j++){
					if (!fscanf(frinput,"%lf ", &coef_temp)){
						printf("The data size is incorrect!");
						break;
					}
					coef[i][j]=fmax(coef[i][j],coef_temp);
				}
			}
		}
		delete[] frfile;
		fclose(frinput);
    }
}

int getMaxPrimeFactor(int n) {
   int i, max = -1;
   while(n % 2 == 0) {
      max = 2;
      n = n/2; //reduce n by dividing this by 2
   }
   for(i = 3; i <= sqrt(n); i=i+2){ //i will increase by 2, to get only odd numbers
      while(n % i == 0) {
         max = i;
         n = n/i;
      }
   }
   if(n > 2) {
      max = n;
   }
   return max;
}

void fun_param_init()
{	
	//spatial setup
	int Nx=calNpoint(dom_x.start, dom_x.end, dom_x.delta);	
	
	//checking prime factor of Nx
	int primeNx 	= getMaxPrimeFactor(Nx);
	if (primeNx>13){
		printf("Number of points in x axis has a largest prime factor>13.\nPlease adjust dx!\n");
		exit(0);
	}	
	dominfo.x=new double[Nx];
	fun_intervaltovector(dominfo.x,dom_x.start,dom_x.end,Nx);
	
	int Ny=calNpoint(dom_y.start, dom_y.end, dom_y.delta);
	//checking prime factor of Ny
	int primeNy 	= getMaxPrimeFactor(Ny);
	if (primeNy>13){
		printf("Number of points in y axis has a largest prime factor>13.\nPlease adjust dy!\n");
		exit(0);
	}
	dominfo.y=new double[Ny];
	fun_intervaltovector(dominfo.y,dom_y.start,dom_y.end,Ny);	

	dominfo.Nx=Nx;
	dominfo.Ny=Ny;
	dominfo.dx=dominfo.x[1]-dominfo.x[0];
	dominfo.dy=dominfo.y[1]-dominfo.y[0];	

	//fourier space
	dominfo.kx=new double[Nx];
	freqspace(dominfo.kx,(double) dominfo.x[Nx-1]-dominfo.x[0],Nx);
	dominfo.ky=new double[Ny];
	freqspace(dominfo.ky,(double) dominfo.y[Ny-1]-dominfo.y[0],Ny);

	dominfo.KK=declare_2darray(dominfo.Nx,dominfo.Ny);
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dominfo.KK[i][j]=sqrt(pow(dominfo.kx[j],2)+pow(dominfo.ky[i],2));
		}
	}

	//bathymetry
	dominfo.bathy.data=declare_2darray(Nx,Ny);
	dominfo.bathy.pos=declare_2darray(Nx,Ny);
	dominfo.bathy.plus=declare_2darray(Nx,Ny);
	dominfo.bathy.min=declare_2darray(Nx,Ny);
	dominfo.bathy.Char=declare_2darray(Nx,Ny);
	fun_bathy2d(dominfo.bathy.data,Nx,Ny,dominfo.bathy);

	//aliasing
	dominfo.aal=declare_2darray(Nx,Ny);	
	fun_aal2D(dominfo.aal,evolinfo.cutfrac,dominfo.kx,dominfo.ky,dominfo.KK,dominfo.Nx,dominfo.Ny);

	//temporal setup
	int Nt=calNpoint(dom_t.start, dom_t.end, dom_t.delta);
	if ((Nt%2)!=0) Nt=Nt-1;
	dominfo.t=new double[Nt];
	fun_intervaltovector_ds(dominfo.t,dom_t.start,dom_t.delta,Nt); 
	dominfo.Nt=Nt;
	dominfo.dt=dominfo.t[2]-dominfo.t[1];	

	//wall setting
	dominfo.wallchar = declare_2darray(Nx,Ny);
	set_matrix_val(dominfo.wallchar,Nx,Ny,1);
	Sw_hat = declare_2darray_complex(Nx,Ny);
	set_matrix_complex_val(Sw_hat,Nx,Ny,0);
	if (strcmp(wall,"no")!=0){
		wall_gam 	 = declare_2darray(Nx,Ny);
		wall_skewgam = declare_2darray(Nx,Ny);
		grad_wallchar_x = declare_2darray(Nx,Ny);
		grad_wallchar_y = declare_2darray(Nx,Ny);
		fun_wallchar(dominfo.x,dominfo.y,Nx,Ny,dominfo.dx,dominfo.dy);
		fun_wall_influxchar(wall_gam,wall_skewgam,dominfo.x,dominfo.y,dominfo.dx,dominfo.dy,Nx,Ny);
		fun_gradient(grad_wallchar_x,grad_wallchar_y,dominfo.wallchar,Nx,Ny,dominfo.dx,dominfo.dy);
		iterW  = 0;
		tprevW = dominfo.t[0];
		Swall_ts   = new double[n_idxwall];
		Swall_prev = new double[n_idxwall];
		set_array_val(Swall_ts,n_idxwall,0);
		Swt_skew = declare_2darray(Nx,Ny);
		Swt_skewline = new double[n_idxwall];
		set_init_double(Swt_skewline,n_idxwall);
	}

	//damping
	dominfo.fbdy.charac=declare_2darray(Nx,Ny);
	fun_fbdy2D(dominfo.fbdy.charac,dominfo.x,dominfo.y,dominfo.Nx,dominfo.Ny,dominfo.fbdy);

	//friction coefficient
	fric_coef = declare_2darray(Nx,Ny);
	fun_friction(fric_coef,Nx,Ny);

	//gauges index
	if (ngauge>0){
		FILE* gaugeinput;
		char* gaugefile = new char[256];
		sprintf(gaugefile,"%s/gauges.dat",arg1);
		gaugeinput = fopen(gaugefile, "r");
		if(gaugeinput==NULL){
			printf("Can not open %s!!\n",gaugefile);
			exit(0);
		}
		xgauge = new double[ngauge];
		ygauge = new double[ngauge];
		index_gauge[0] = new int[ngauge];
		index_gauge[1] = new int[ngauge];
		for (int j=0;j<ngauge;j++){
			fscanf(gaugeinput,"%lf %lf\n",&xgauge[j],&ygauge[j]);
			index_gauge[0][j] = round((xgauge[j]-dominfo.x[0])/dominfo.dx);
			index_gauge[1][j] = round((ygauge[j]-dominfo.y[0])/dominfo.dy);
			if ((index_gauge[0][j]>=Nx) || (index_gauge[1][j]>=Ny)){
				printf("Gauges are outside of the domain!\n");exit(0);
			}
		}
		fclose(gaugeinput);
		delete[] gaugefile;
	}	
}

void cfv(double** ChiAdj,int id,double* xi,double Xi,double del,int Nx)
{	
	for (int j=0;j<Nx;j++) {
		ChiAdj[id][j] = fun_heav(xi[j]-Xi)*GSL_MAX(fun_sign(xi[j]-(Xi+del)),(1-cos((xi[j]-Xi)*M_PI/del))/2.0);
	}
}

void cfh(double** ChiAdj,int id,double* xi,double Xi,double del,int Nx)
{	
	for (int j=0;j<Nx;j++) {
		ChiAdj[j][id] = (fun_heav(xi[j]-Xi))*GSL_MAX(fun_sign(xi[j]-(Xi+del)),(1-cos((xi[j]-Xi)*M_PI/del))/2.0);
	}
}

void fun_bathy2d(double** bathy,int Nx,int Ny,bathyvar bathdat)
{	
	FILE* btyinput;
	char* btyfile = new char[256];
	if(strcmp(bathdat.Id,"user_defined")==0){
		sprintf(btyfile,"%s/bath.dat",arg1);
		btyinput = fopen(btyfile, "r");
		if(btyinput==NULL){
			printf("Can not open %s!!\n",btyfile);
			exit(0);
		} 
    }
		
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (strcmp(bathdat.Id,"flat")==0){
				bathy[i][j]=-dominfo.bathy.depth;
			}
			else if(strcmp(bathdat.Id,"user_defined")==0){
				if (!fscanf(btyinput,"%lf ", &bathy[i][j])){
					printf("The data size is incorrect!");
					break;
				}
				bathy[i][j] = -bathy[i][j];
			}
		}
	}
    delete[] btyfile;
    
	if(strcmp(bathdat.Id,"user_defined")==0){
		fclose(btyinput);
	}

	//bathy characteristic adn runup
	set_matrix_val(dominfo.bathy.Char,Nx,Ny,1.0);
	double max_bath;
	max_bath = fun_max2darray(bathy,Nx,Ny);
	if (max_bath>0){
		runup_id = 1;//runup exist along Y
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (dominfo.bathy.data[i][j]>0){
					dominfo.bathy.Char[i][j]=0;
				}
			}
		}
	}
	else{
		runup_id = 0;
	}

	//bathyplus and min
	if (runup_id==1){
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (dominfo.bathy.data[i][j]<0){
					dominfo.bathy.min[i][j]=dominfo.bathy.data[i][j];
					dominfo.bathy.plus[i][j]=0;
				}
				else{
					dominfo.bathy.min[i][j]=0;
					dominfo.bathy.plus[i][j]=dominfo.bathy.data[i][j];
				}
			}
		}
	}

	//bathy.pos
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			dominfo.bathy.pos[i][j]=-dominfo.bathy.data[i][j];
		}
	}
}

void fun_calcseacharacteristic(waveinf wave_influx)
{
	if (strcmp(wavename,"zero")==0){
		seainfo.wp = 2*M_PI/seainfo.Tp;
		seainfo.kp = invers_omega(seainfo.wp,meandepthinf);
		seainfo.lambda=2*M_PI/seainfo.kp;
		ChiAdj = declare_2darray(dominfo.Nx,dominfo.Ny);
		set_matrix_val(ChiAdj,dominfo.Nx,dominfo.Ny,1.0);
	}
	else {
		Insig_hat  = declare_2darray_complex(wave_influx.Nxy,dominfo.Nt);
		fftw_2d_r2c(Insig_hat,wave_influx.eta,wave_influx.Nxy,dominfo.Nt);
		
		if (strcmp(wavename,"user_defined")==0){
			seainfo.Hs = 4.0*sqrt(fun_var2darray(wave_influx.eta,wave_influx.Nxy,dominfo.Nt));
		}
		printf("Hs=%.2f\n",seainfo.Hs);
		double** abs_eta = fun_abs_complex(Insig_hat,wave_influx.Nxy,dominfo.Nt);
		int* idx_max = new int[2];
		idx_max_matrix(idx_max,abs_eta,wave_influx.Nxy,dominfo.Nt);
		int id = idx_max[1];	
		seainfo.wp = fabs(wave_influx.ww[id]);
		seainfo.Tp = 2*M_PI/seainfo.wp;
		seainfo.kp = invers_omega(seainfo.wp,meandepthinf);
		seainfo.lambda=2*M_PI/seainfo.kp;
		
		delete[] idx_max;
		free_2darray(abs_eta,wave_influx.Nxy,dominfo.Nt);

		//nonlinear adjustment
		ChiAdj    = declare_2darray(dominfo.Nx,dominfo.Ny);
		//set_matrix_val(ChiAdj,dominfo.Nx,dominfo.Ny,1.0);
		if (strcmp(evolinfo.model,"HS1")!=0){		
			if (dominfo.Nadj==0.0){
				set_matrix_val(ChiAdj,dominfo.Nx,dominfo.Ny,1.0);
			}
			else {
				set_matrix_val(ChiAdj,dominfo.Nx,dominfo.Ny,1.0);
				int indx,indy;
				double** ChiAdj1 = declare_2darray(dominfo.Nx,dominfo.Ny);
				double** ChiAdj2 = declare_2darray(dominfo.Nx,dominfo.Ny);
				if (strcmp(orientation,"vertical")==0){
					if  ((wave_influx.y[0]==dominfo.y[0])&&(wave_influx.y[wave_influx.nlines-1]==dominfo.y[dominfo.Ny-1])){	
						for (int j=0;j<wave_influx.nlines;j++){
							indx = round((wave_influx.x[j]-dominfo.x[0])/dominfo.dx);
							indy = round((wave_influx.y[j]-dominfo.y[0])/dominfo.dy);
							//cfv(ChiAdj,indy,dominfo.x,dominfo.x[indx],dominfo.Nadj*seainfo.lambda,dominfo.Nx);
							cfv(ChiAdj1,indy,dominfo.x,dominfo.x[indx],dominfo.Nadj*seainfo.lambda,dominfo.Nx);
							cfv(ChiAdj2,indy,dominfo.x,dominfo.x[indx]-dominfo.Nadj*seainfo.lambda,dominfo.Nadj*seainfo.lambda,dominfo.Nx);
							for (int i=0;i<dominfo.Nx;i++)	{
								ChiAdj[indy][i] = 1-ChiAdj2[indy][i]+ChiAdj1[indy][i];
							}
						}
					}
					else{
						fbdyvar ftemp;
						double dista = dominfo.Nadj*seainfo.lambda;
						ftemp.Lleft=dista;
						ftemp.Lright=dista;
						ftemp.Ltop=dista;
						ftemp.Lbottom=dista;
						double xtemp1,xtemp2,ytemp1,ytemp2;
						xtemp1=wave_influx.x[0]-dista;
						xtemp2=wave_influx.x[wave_influx.nlines-1]+dista;
						ytemp1=wave_influx.y[0]-dista;
						ytemp2=wave_influx.y[wave_influx.nlines-1]+dista;
						int nxtemp=int((xtemp2-xtemp1)/dominfo.dx)+1;
						int nytemp=int((ytemp2-ytemp1)/dominfo.dy)+1;
						double* xtemp=new double[nxtemp];
						double* ytemp=new double[nytemp];
						for (int i=0;i<nxtemp;i++){
							xtemp[i]=xtemp1+i*dominfo.dx;
						}
						for (int i=0;i<nytemp;i++){
							ytemp[i]=ytemp1+i*dominfo.dy;
						}
						double** ChiAdjtemp=declare_2darray(nxtemp,nytemp);						
						fun_cfSA2D(ChiAdjtemp,xtemp,ytemp,nxtemp,nytemp,ftemp);
						fun_interp2(nxtemp,nytemp,xtemp,ytemp,ChiAdjtemp,dominfo.Nx,dominfo.Ny,dominfo.x,dominfo.y,ChiAdj);
						for (int i=0;i<dominfo.Ny;i++){
							for(int j=0;j<dominfo.Nx;j++){								
								if ((ChiAdj[i][j])==(-9999)){
									ChiAdj[i][j]=1;
								}
								else{
									ChiAdj[i][j]=1-ChiAdj[i][j];
									if(abs(ChiAdj[i][j])<0.001){
										ChiAdj[i][j]=0;
									}
								}								
							}							
						}
						
						delete[] xtemp;
						delete[] ytemp;
						free_2darray(ChiAdjtemp,nxtemp,nytemp);
					}
				}
				else {		
					if  ((wave_influx.x[0]==dominfo.x[0])&&(wave_influx.x[wave_influx.nlines-1]==dominfo.x[dominfo.Nx-1])){			
						for (int j=0;j<wave_influx.nlines;j++){
							indx = round((wave_influx.x[j]-dominfo.x[0])/dominfo.dx);
							indy = round((wave_influx.y[j]-dominfo.y[0])/dominfo.dy);
							//cfh(ChiAdj,indx,dominfo.y,dominfo.y[indy],dominfo.Nadj*seainfo.lambda,dominfo.Ny);
							cfh(ChiAdj1,indx,dominfo.y,dominfo.y[indy],dominfo.Nadj*seainfo.lambda,dominfo.Ny);
							cfh(ChiAdj2,indx,dominfo.y,dominfo.y[indy]-dominfo.Nadj*seainfo.lambda,dominfo.Nadj*seainfo.lambda,dominfo.Ny);
							for (int i=0;i<dominfo.Ny;i++)	{
								ChiAdj[i][indx] = 1-ChiAdj2[i][indx]+ChiAdj1[i][indx];
							}						
						}
					}
					else{
						fbdyvar ftemp;
						double dista = dominfo.Nadj*seainfo.lambda;
						ftemp.Lleft=dista;
						ftemp.Lright=dista;
						ftemp.Ltop=dista;
						ftemp.Lbottom=dista;
						double xtemp1,xtemp2,ytemp1,ytemp2;
						xtemp1=wave_influx.x[0]-dista;
						xtemp2=wave_influx.x[wave_influx.nlines-1]+dista;
						ytemp1=wave_influx.y[0]-dista;
						ytemp2=wave_influx.y[wave_influx.nlines-1]+dista;
						int nxtemp=int((xtemp2-xtemp1)/dominfo.dx)+1;
						int nytemp=int((ytemp2-ytemp1)/dominfo.dy)+1;
						double* xtemp=new double[nxtemp];
						double* ytemp=new double[nytemp];
						for (int i=0;i<nxtemp;i++){
							xtemp[i]=xtemp1+i*dominfo.dx;
						}
						for (int i=0;i<nytemp;i++){
							ytemp[i]=ytemp1+i*dominfo.dy;
						}
						double** ChiAdjtemp=declare_2darray(nxtemp,nytemp);						
						fun_cfSA2D(ChiAdjtemp,xtemp,ytemp,nxtemp,nytemp,ftemp);
						fun_interp2(nxtemp,nytemp,xtemp,ytemp,ChiAdjtemp,dominfo.Nx,dominfo.Ny,dominfo.x,dominfo.y,ChiAdj);
						for (int i=0;i<dominfo.Ny;i++){
							for(int j=0;j<dominfo.Nx;j++){								
								if ((ChiAdj[i][j])==(-9999)){
									ChiAdj[i][j]=1;
								}
								else{
									ChiAdj[i][j]=1-ChiAdj[i][j];
									if(abs(ChiAdj[i][j])<0.001){
										ChiAdj[i][j]=0;
									}
								}								
							}							
						}
						
						delete[] xtemp;
						delete[] ytemp;
						free_2darray(ChiAdjtemp,nxtemp,nytemp);
					}
				}
				free_2darray(ChiAdj1,dominfo.Nx,dominfo.Ny);
				free_2darray(ChiAdj2,dominfo.Nx,dominfo.Ny);
			}
		}
		else{
			set_matrix_val(ChiAdj,dominfo.Nx,dominfo.Ny,1.0);
		}

		//nonlinear adjustment
		/*FILE* gp4=popen("gnuplot -persistent","w");
		fun_plot_2d_trarray(gp4,dominfo.y,dominfo.x,ChiAdj,dominfo.Ny,dominfo.Nx);
		fflush(gp4);
		pclose(gp4);*/
	
		//adjustment to smooth Phi
		double** ChiAdjPhi1 = declare_2darray(dominfo.Nx,dominfo.Ny);
		double** ChiAdjPhi2 = declare_2darray(dominfo.Nx,dominfo.Ny);
		ChiAdjPhi = declare_2darray(dominfo.Nx,dominfo.Ny);
		set_matrix_val(ChiAdjPhi,dominfo.Nx,dominfo.Ny,1.0);
		
		int ix,iy;
		if (strcmp(orientation,"vertical")==0){
			for (int j=0;j<wave_influx.nlines;j++){
				ix = round((wave_influx.x[j]-dominfo.x[0])/dominfo.dx);
				iy = round((wave_influx.y[j]-dominfo.y[0])/dominfo.dy);
				cfv(ChiAdjPhi1,iy,dominfo.x,dominfo.x[ix],seainfo.lambda,dominfo.Nx);
				cfv(ChiAdjPhi2,iy,dominfo.x,dominfo.x[ix]-seainfo.lambda,seainfo.lambda,dominfo.Nx);
				for (int i=0;i<dominfo.Nx;i++)	{
					ChiAdjPhi[iy][i] = 1-ChiAdjPhi2[iy][i]+ChiAdjPhi1[iy][i];
				}
			}
		}
		else {
			for (int j=0;j<wave_influx.nlines;j++){
				ix = round((wave_influx.x[j]-dominfo.x[0])/dominfo.dx);
				iy = round((wave_influx.y[j]-dominfo.y[0])/dominfo.dy);
				cfh(ChiAdjPhi1,ix,dominfo.y,dominfo.y[iy],seainfo.lambda,dominfo.Ny);
				cfh(ChiAdjPhi2,ix,dominfo.y,dominfo.y[iy]-seainfo.lambda,seainfo.lambda,dominfo.Ny);
				for (int i=0;i<dominfo.Ny;i++)	{
					ChiAdjPhi[i][ix] = 1-ChiAdjPhi2[i][ix]+ChiAdjPhi1[i][ix];
				}
			}				
		}
		free_2darray(ChiAdjPhi1,dominfo.Nx,dominfo.Ny);
		free_2darray(ChiAdjPhi2,dominfo.Nx,dominfo.Ny);		
		
		/*FILE* gp=popen("gnuplot -persistent","w");
		fun_plot_2d_trarray(gp,dominfo.y,dominfo.x,ChiAdjPhi,dominfo.Ny,dominfo.Nx);
		fflush(gp);
		pclose(gp);getchar();*/
	}	
	
}

void funOprt_runup2D(int Nx,int Ny) //if runup exist
{
	double maxetainit = 1.5*seainfo.Hs;
	
	//along X
	double H_mid,H_plus,H_min,k_cut,wcut;
	double nupeak;
	int id;

	//compute kmax
	if (strcmp(orientation,"vertical")==0){
		id      = (int)ceil(Nx/evolinfo.cutfrac);	
		k_cut   = fabs(dominfo.kx[id]);
	}
	else{
		id      = (int)ceil(Ny/evolinfo.cutfrac);
		k_cut   = fabs(dominfo.ky[id]);
	}

	H_plus= -fun_min2darray(dominfo.bathy.data,Nx,Ny)+maxetainit;
	nupeak= 2*M_PI/seainfo.Tp;
	H_minShore  = 0.01*seainfo.Hs;//pow((3.0*nupeak/4.0/k_cut),2)/grav;
	H_min 		= Dmin_input;
	
	double Dperlambda,KD,Hmintemp;
	Dperlambda=1.0/200;//this is for very longwaves i.e T=100 s, if 1/20 (SWE) this leads to D~~200 (too much)
    KD=2.0*M_PI*Dperlambda;
    Hmintemp=grav*(KD*tanh(KD))/(pow(nupeak,2));
    H_min=fmax(H_min,Hmintemp);
	
	if (strcmp(mid_temp,"default")==0){	
		H_mid = (H_plus+H_min)/2;
	}
	else{
		H_mid = atof(mid_temp);
	}
	
	int Npoints=100;
	double* Htot1 = new double[Npoints];
	double* Htot2 = new double[Npoints];
	double  dh1,dh2,H_min1,H_plus1,H_mid1,H_min2,H_plus2,H_mid2;

	dh1 =  (H_mid-H_min)/(Npoints-1);
	dh2 =  (H_plus-H_mid)/(Npoints-1);
	for (int j=0;j<Npoints;j++){
		Htot1[j] = H_min+j*dh1;
		Htot2[j] = H_mid+j*dh2;
	}
	H_min1	= H_min;
	H_plus1	= H_mid;
	H_mid1	= (H_min1+H_plus1)/2;
	
	H_min2	= H_mid;
	H_plus2	= H_plus;
	H_mid2	= H_min2+(H_plus2-H_min2)/2; //IF deep water not included yet!!!
	wcut	= omega_i(k_cut,H_plus);
	
	double* Cp2_H1 		= new double[Npoints];
	double* Cp2_min1 	= new double[Npoints];
	double* Cp2_plus1 	= new double[Npoints];
	double* Cp2_mid1 	= new double[Npoints];
	double* Cp2kL_H1 	= new double[Npoints];
	double* Cp2kL_min1 	= new double[Npoints];
	double* Cp2kL_plus1 = new double[Npoints];
	double* Cp2kL_mid1 	= new double[Npoints];
	double* Cp2IG_H1 	= new double[Npoints];
	double* Cp2IG_min1 	= new double[Npoints];
	double* Cp2IG_plus1 = new double[Npoints];
	double* Cp2IG_mid1 	= new double[Npoints];
	double* gam_min1 	= new double[Npoints];
	double* gam_plus1 	= new double[Npoints];
	double* gam_mid1 	= new double[Npoints];

	double* Cp2_H2 		= new double[Npoints];
	double* Cp2_min2 	= new double[Npoints];
	double* Cp2_plus2 	= new double[Npoints];
	double* Cp2_mid2 	= new double[Npoints];
	double* Cp2kL_H2 	= new double[Npoints];
	double* Cp2kL_min2 	= new double[Npoints];
	double* Cp2kL_plus2 = new double[Npoints];
	double* Cp2kL_mid2 	= new double[Npoints];
	double* Cp2IG_H2 	= new double[Npoints];
	double* Cp2IG_min2 	= new double[Npoints];
	double* Cp2IG_plus2 = new double[Npoints];
	double* Cp2IG_mid2 	= new double[Npoints];
	double* knu1        = new double[Npoints];
	double* knu2        = new double[Npoints];
	double* k_half 		= new double[int(Nx/2)];
	double* gam_min2 	= new double[Npoints];
	double* gam_plus2 = new double[Npoints];
	double* gam_mid2 	= new double[Npoints];

	for (int j=0;j<int(Nx/2);j++){
		k_half[j] = dominfo.kx[j];
	}

	double Hj1,Hj2,kappanu_H1,kappanuIG_H1,kappanuL_H1;
	double kappanu_H2,kappanuIG_H2,kappanuL_H2;
	double A1,B1,C1,D1,E1,F1,G1,H1,I1,det1;
	double A2,B2,C2,D2,E2,F2,G2,H2,I2,det2;
	
	for (int j=0;j<Npoints;j++){
		Hj1	= Htot1[j];
		Hj2	= Htot2[j];   
		kappanu_H1  	= invers_omega(nupeak,Hj1);
		kappanuIG_H1   	= invers_omega(nupeak/10.0,Hj1); 
		kappanuL_H1 	= invers_omega(wcut,Hj1);
		kappanu_H2  	= invers_omega(nupeak,Hj2);
		kappanuIG_H2   	= invers_omega(nupeak/10.0,Hj2);
		kappanuL_H2 	= invers_omega(wcut,Hj2);

		knu1[j] = kappanu_H1;
		knu2[j] = kappanu_H2;
		
		Cp2_min1[j] = pow(phase_velocity_i(kappanu_H1,H_min1),2);
		Cp2_plus1[j]= pow(phase_velocity_i(kappanu_H1,H_plus1),2);
		Cp2_mid1[j] = pow(phase_velocity_i(kappanu_H1,H_mid1),2);
		Cp2_H1[j]   = pow(phase_velocity_i(kappanu_H1,Hj1),2);

		Cp2kL_min1[j] = pow(phase_velocity_i(kappanuL_H1,H_min1),2);
		Cp2kL_plus1[j]= pow(phase_velocity_i(kappanuL_H1,H_plus1),2);
		Cp2kL_mid1[j] = pow(phase_velocity_i(kappanuL_H1,H_mid1),2);
		Cp2kL_H1[j]   = pow(phase_velocity_i(kappanuL_H1,Hj1),2);

		Cp2IG_min1[j] = pow(phase_velocity_i(kappanuIG_H1,H_min1),2);
		Cp2IG_plus1[j]= pow(phase_velocity_i(kappanuIG_H1,H_plus1),2);
		Cp2IG_mid1[j] = pow(phase_velocity_i(kappanuIG_H1,H_mid1),2);
		Cp2IG_H1[j]   = pow(phase_velocity_i(kappanuIG_H1,Hj1),2);

		Cp2_min2[j] = pow(phase_velocity_i(kappanu_H2,H_min2),2);
		Cp2_plus2[j]= pow(phase_velocity_i(kappanu_H2,H_plus2),2);
		Cp2_mid2[j] = pow(phase_velocity_i(kappanu_H2,H_mid2),2);
		Cp2_H2[j]   = pow(phase_velocity_i(kappanu_H2,Hj2),2);

		Cp2kL_min2[j] = pow(phase_velocity_i(kappanuL_H2,H_min2),2);
		Cp2kL_plus2[j]= pow(phase_velocity_i(kappanuL_H2,H_plus2),2);
		Cp2kL_mid2[j] = pow(phase_velocity_i(kappanuL_H2,H_mid2),2);
		Cp2kL_H2[j]   = pow(phase_velocity_i(kappanuL_H2,Hj2),2);

		Cp2IG_min2[j] = pow(phase_velocity_i(kappanuIG_H2,H_min2),2);
		Cp2IG_plus2[j]= pow(phase_velocity_i(kappanuIG_H2,H_plus2),2);
		Cp2IG_mid2[j] = pow(phase_velocity_i(kappanuIG_H2,H_mid2),2);
		Cp2IG_H2[j]   = pow(phase_velocity_i(kappanuIG_H2,Hj2),2);
	
		A1= (Cp2kL_mid1[j]*Cp2IG_plus1[j])-(Cp2kL_plus1[j]*Cp2IG_mid1[j]);
		B1=-(Cp2kL_min1[j]*Cp2IG_plus1[j])+(Cp2kL_plus1[j]*Cp2IG_min1[j]);
		C1= (Cp2kL_min1[j]*Cp2IG_mid1[j])-(Cp2kL_mid1[j]*Cp2IG_min1[j]);

		D1=-(Cp2_mid1[j]*Cp2IG_plus1[j])+(Cp2_plus1[j]*Cp2IG_mid1[j]);
		E1= (Cp2_min1[j]*Cp2IG_plus1[j])-(Cp2_plus1[j]*Cp2IG_min1[j]);
		F1=-(Cp2_min1[j]*Cp2IG_mid1[j])+(Cp2_mid1[j]*Cp2IG_min1[j]);

		G1= (Cp2_mid1[j]*Cp2kL_plus1[j])-(Cp2_plus1[j]*Cp2kL_mid1[j]);
		H1=-(Cp2_min1[j]*Cp2kL_plus1[j])+(Cp2_plus1[j]*Cp2kL_min1[j]);
		I1= (Cp2_min1[j]*Cp2kL_mid1[j])-(Cp2_mid1[j]*Cp2kL_min1[j]);

		det1=Cp2_min1[j]*A1+Cp2_mid1[j]*B1+Cp2_plus1[j]*C1;

		gam_min1[j]  = (A1*Cp2_H1[j]+D1*Cp2kL_H1[j]+G1*Cp2IG_H1[j])/det1;
		gam_mid1[j]  = (B1*Cp2_H1[j]+E1*Cp2kL_H1[j]+H1*Cp2IG_H1[j])/det1;
		gam_plus1[j] = (C1*Cp2_H1[j]+F1*Cp2kL_H1[j]+I1*Cp2IG_H1[j])/det1;

		A2= (Cp2kL_mid2[j]*Cp2IG_plus2[j])-(Cp2kL_plus2[j]*Cp2IG_mid2[j]);
		B2=-(Cp2kL_min2[j]*Cp2IG_plus2[j])+(Cp2kL_plus2[j]*Cp2IG_min2[j]);
		C2= (Cp2kL_min2[j]*Cp2IG_mid2[j])-(Cp2kL_mid2[j]*Cp2IG_min2[j]);

		D2=-(Cp2_mid2[j]*Cp2IG_plus2[j])+(Cp2_plus2[j]*Cp2IG_mid2[j]);
		E2= (Cp2_min2[j]*Cp2IG_plus2[j])-(Cp2_plus2[j]*Cp2IG_min2[j]);
		F2=-(Cp2_min2[j]*Cp2IG_mid2[j])+(Cp2_mid2[j]*Cp2IG_min2[j]);

		G2= (Cp2_mid2[j]*Cp2kL_plus2[j])-(Cp2_plus2[j]*Cp2kL_mid2[j]);
		H2=-(Cp2_min2[j]*Cp2kL_plus2[j])+(Cp2_plus2[j]*Cp2kL_min2[j]);
		I2= (Cp2_min2[j]*Cp2kL_mid2[j])-(Cp2_mid2[j]*Cp2kL_min2[j]);

		det2=Cp2_min2[j]*A2+Cp2_mid2[j]*B2+Cp2_plus2[j]*C2;

		gam_min2[j]  = ( A2*Cp2_H2[j]+D2*Cp2kL_H2[j]+G2*Cp2IG_H2[j])/det2;
		gam_mid2[j]  = ( B2*Cp2_H2[j]+E2*Cp2kL_H2[j]+H2*Cp2IG_H2[j])/det2;
		gam_plus2[j] = ( C2*Cp2_H2[j]+F2*Cp2kL_H2[j]+I2*Cp2IG_H2[j])/det2;
	}

	int 	Np=2*Npoints-1;
	double* HtotComb    = new double[Np];
	double* gam_min1tot = new double[Np];
	double* gam_mid1tot = new double[Np];
	double* gam_plus1tot= new double[Np];
	double* gam_min2tot = new double[Np];
	double* gam_mid2tot = new double[Np];
	double* gam_plus2tot= new double[Np];
	
	set_init_double(HtotComb,Np);
	set_init_double(gam_min1tot,Np);
	set_init_double(gam_mid1tot,Np);
	set_init_double(gam_plus1tot,Np);
	set_init_double(gam_min2tot,Np);
	set_init_double(gam_mid2tot,Np);
	set_init_double(gam_plus2tot,Np);	

	for (int j=0;j<(Np);j++){
		if (j<Npoints){
			HtotComb[j] 	= Htot1[j];
			gam_min1tot[j]  = gam_min1[j];
			gam_mid1tot[j]  = gam_mid1[j];
			gam_plus1tot[j] = gam_plus1[j];
			gam_min2tot[j]	= 0;
			gam_mid2tot[j]	= 0;
			gam_plus2tot[j]	= 0;
		}
		else{
			HtotComb[j] 	= Htot2[j+1-Npoints];
			gam_min1tot[j]  = 0;
			gam_mid1tot[j]  = 0;
			gam_plus1tot[j] = 0;
			gam_min2tot[j]	= gam_min2[j+1-Npoints];
			gam_mid2tot[j]	= gam_mid2[j+1-Npoints];
			gam_plus2tot[j]	= gam_plus2[j+1-Npoints];
		}
	}
	
	//allocation
	acc_p1  	= gsl_interp_accel_alloc ();
	sg_p1 		= gsl_spline_alloc (gsl_interp_cspline, Np);

	acc_c1  	= gsl_interp_accel_alloc ();
	sg_c1 		= gsl_spline_alloc (gsl_interp_cspline, Np);

	acc_m1 		= gsl_interp_accel_alloc ();
	sg_m1 		= gsl_spline_alloc (gsl_interp_cspline, Np);

	acc_p2 		= gsl_interp_accel_alloc ();
	sg_p2 		= gsl_spline_alloc (gsl_interp_cspline, Np);

	acc_c2 		= gsl_interp_accel_alloc ();
	sg_c2 		= gsl_spline_alloc (gsl_interp_cspline, Np);;

	acc_m2 		= gsl_interp_accel_alloc ();
	sg_m2 		= gsl_spline_alloc (gsl_interp_cspline, Np);    

	//interpolation    		  
	gsl_spline_init (sg_p1,HtotComb,gam_plus1tot,Np);
	gsl_spline_init (sg_c1,HtotComb,gam_mid1tot,Np);
	gsl_spline_init (sg_m1,HtotComb,gam_min1tot,Np);
	gsl_spline_init (sg_p2,HtotComb,gam_plus2tot,Np);
	gsl_spline_init (sg_c2,HtotComb,gam_mid2tot,Np);
	gsl_spline_init (sg_m2,HtotComb,gam_min2tot,Np);
	
	//For HS
	gam_m1  = declare_2darray(Nx,Ny);
	gam_c1  = declare_2darray(Nx,Ny);
	gam_p1  = declare_2darray(Nx,Ny);
	gam_m2  = declare_2darray(Nx,Ny);
	gam_c2  = declare_2darray(Nx,Ny);
	gam_p2  = declare_2darray(Nx,Ny);
	
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (-dominfo.bathy.data[i][j]<H_min){
				gam_m1[i][j]  = 0;
				gam_c1[i][j]  = 0;
				gam_p1[i][j]  = 0;
				gam_m2[i][j]  = 0;
				gam_c2[i][j]  = 0;
				gam_p2[i][j]  = 0;
			}
			else{
				gam_m1[i][j]  = gsl_spline_eval(sg_m1,-dominfo.bathy.data[i][j],acc_m1);
				gam_c1[i][j]  = gsl_spline_eval(sg_c1,-dominfo.bathy.data[i][j],acc_c1);
				gam_p1[i][j]  = gsl_spline_eval(sg_p1,-dominfo.bathy.data[i][j],acc_p1);
				gam_m2[i][j]  = gsl_spline_eval(sg_m2,-dominfo.bathy.data[i][j],acc_m2);
				gam_c2[i][j]  = gsl_spline_eval(sg_c2,-dominfo.bathy.data[i][j],acc_c2);
				gam_p2[i][j]  = gsl_spline_eval(sg_p2,-dominfo.bathy.data[i][j],acc_p2);
			}
		}
	}		
	
	delete[] Htot1; delete[] Htot2;
	delete[] Cp2_H1; 		
	delete[] Cp2_min1; 	
	delete[] Cp2_plus1; 	
	delete[] Cp2_mid1; 
	delete[] Cp2kL_H1; 
	delete[] Cp2kL_min1;
	delete[] Cp2kL_plus1;
	delete[] Cp2kL_mid1;
	delete[] Cp2IG_H1;
	delete[] Cp2IG_min1; 
	delete[] Cp2IG_plus1; 
	delete[] Cp2IG_mid1;
	delete[] gam_min1; 
	delete[] gam_plus1; 
	delete[] gam_mid1; 

	delete[] Cp2_H2;
	delete[] Cp2_min2;
	delete[] Cp2_plus2;
	delete[] Cp2_mid2;
	delete[] Cp2kL_H2;
	delete[] Cp2kL_min2;
	delete[] Cp2kL_plus2;
	delete[] Cp2kL_mid2;
	delete[] Cp2IG_H2;
	delete[] Cp2IG_min2;
	delete[] Cp2IG_plus2;
	delete[] Cp2IG_mid2;
	delete[] knu1;
	delete[] knu2;
	delete[] k_half;
	delete[] gam_min2;
	delete[] gam_plus2;
	delete[] gam_mid2;
	delete[] HtotComb;

	delete[] gam_min1tot;
	delete[] gam_mid1tot;
	delete[] gam_plus1tot;
	delete[] gam_min2tot;
	delete[] gam_mid2tot;
	delete[] gam_plus2tot;
	
	//operator C2
	C2m1 = declare_2darray(Nx,Ny);
	C2p1 = declare_2darray(Nx,Ny);
	C2c1 = declare_2darray(Nx,Ny);
	C2m2 = declare_2darray(Nx,Ny);
	C2p2 = declare_2darray(Nx,Ny);
	C2c2 = declare_2darray(Nx,Ny);
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Nx;j++){
			C2m1[i][j] = pow(fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_min1,seainfo.current),2);
			C2p1[i][j] = pow(fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_plus1,seainfo.current),2);
			C2c1[i][j] = pow(fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_mid1,seainfo.current),2);
			C2m2[i][j] = pow(fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_min2,seainfo.current),2);
			C2p2[i][j] = pow(fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_plus2,seainfo.current),2);
			C2c2[i][j] = pow(fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_mid2,seainfo.current),2);
		}
	}
	
	//operator C
	Cm1 = declare_2darray(Nx,Ny);
	Cp1 = declare_2darray(Nx,Ny);
	Cc1 = declare_2darray(Nx,Ny);
	Cm2 = declare_2darray(Nx,Ny);
	Cp2 = declare_2darray(Nx,Ny);
	Cc2 = declare_2darray(Nx,Ny);
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Nx;j++){
			Cm1[i][j] = (fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_min1,seainfo.current));
			Cp1[i][j] = (fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_plus1,seainfo.current));
			Cc1[i][j] = (fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_mid1,seainfo.current));
			Cm2[i][j] = (fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_min2,seainfo.current));
			Cp2[i][j] = (fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_plus2,seainfo.current));
			Cc2[i][j] = (fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],H_mid2,seainfo.current));
		}
	}
	
	//operator Om2
	Om2m1 = declare_2darray(Nx,Ny);
	Om2p1 = declare_2darray(Nx,Ny);
	Om2c1 = declare_2darray(Nx,Ny);
	Om2m2 = declare_2darray(Nx,Ny);
	Om2p2 = declare_2darray(Nx,Ny);
	Om2c2 = declare_2darray(Nx,Ny);
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Nx;j++){
			Om2m1[i][j] = pow(fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],H_min1,seainfo.current.x,seainfo.current.y,seainfo.mdir),2);
			Om2p1[i][j] = pow(fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],H_plus1,seainfo.current.x,seainfo.current.y,seainfo.mdir),2);
			Om2c1[i][j] = pow(fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],H_mid1,seainfo.current.x,seainfo.current.y,seainfo.mdir),2);
			Om2m2[i][j] = pow(fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],H_min2,seainfo.current.x,seainfo.current.y,seainfo.mdir),2);
			Om2p2[i][j] = pow(fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],H_plus2,seainfo.current.x,seainfo.current.y,seainfo.mdir),2);
			Om2c2[i][j] = pow(fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],H_mid2,seainfo.current.x,seainfo.current.y,seainfo.mdir),2); 
		}
	}	
	
	Oprt.HminDisp = H_min;
}

Oprtvar fun_operator_setup(domvar dominfo)
{
	Oprtvar Opr;
	int Nx=dominfo.Nx;
	int Ny=dominfo.Ny;
	
	Oprt.Om2d=declare_2darray(Nx,Ny);
	Oprt.L2d=declare_2darray(Nx,Ny);
	Oprt.Cpeak2d=declare_2darray(Nx,Ny);
	Oprt.Csqr2d=declare_2darray(Nx,Ny);
	Oprt.C2d=declare_2darray(Nx,Ny);

	if (strcmp(dominfo.bathy.Id,"flat")!=0){//VARYING BOTTOM 
		Oprt.Cpeak=fun_exact_Cp1d_val(seainfo.kp,meandepthinf,dominfo.U);
		if (runup_id==0){//no runup
			Oprt.InterpD.Gam_min =declare_2darray(Nx,Ny);	
			Oprt.InterpD.Gam_plus=declare_2darray(Nx,Ny);
			Oprt.InterpD.L2d_min =declare_2darray(Nx,Ny);		
			Oprt.InterpD.L2d_plus=declare_2darray(Nx,Ny);
			Oprt.InterpD.Om2d_min =declare_2darray(Nx,Ny);		
			Oprt.InterpD.Om2d_plus=declare_2darray(Nx,Ny);
			Oprt.InterpD.Om2dSq_min =declare_2darray(Nx,Ny);		
			Oprt.InterpD.Om2dSq_plus=declare_2darray(Nx,Ny);
			Oprt.InterpD.C2d_min  =declare_2darray(Nx,Ny);		
			Oprt.InterpD.C2d_plus =declare_2darray(Nx,Ny);
			Oprt.InterpD.Csqr2d_min  =declare_2darray(Nx,Ny);		
			Oprt.InterpD.Csqr2d_plus =declare_2darray(Nx,Ny);
			
			if (ninterp==2){
				funOprt_DispInterpolation_2p(Oprt,dominfo.bathy.data,dominfo.KK,dominfo.kx,dominfo.ky,
										   dominfo.Nx,dominfo.Ny,seainfo.wp,seainfo.current,seainfo.mdir);
										   
			}
			else if (ninterp==3){
				Oprt.InterpD.Gam_mid =declare_2darray(Nx,Ny);
				Oprt.InterpD.L2d_mid =declare_2darray(Nx,Ny);
				Oprt.InterpD.Om2d_mid=declare_2darray(Nx,Ny);
				Oprt.InterpD.Om2dSq_mid=declare_2darray(Nx,Ny);
				Oprt.InterpD.C2d_mid =declare_2darray(Nx,Ny);
				Oprt.InterpD.Csqr2d_mid =declare_2darray(Nx,Ny);
				funOprt_DispInterpolation_3p(Oprt,dominfo.bathy.data,dominfo.KK,dominfo.kx,dominfo.ky,
										   dominfo.Nx,dominfo.Ny,seainfo.wp,seainfo.current,seainfo.mdir);										   
			}
			else{
				printf("Error in parameter ninterp!!\n"); exit(0);
			}
		}
		else{//RUNUP operator
			funOprt_runup2D(dominfo.Nx,dominfo.Ny);
		}
	}
	else{//THIS IS FOR FLAT BOTTOM
		Oprt.Cpeak=fun_exact_Cp1d_val(seainfo.kp,dominfo.bathy.depth,dominfo.U);
		for(int i=0;i<dominfo.Ny;i++){
			for(int j=0;j<dominfo.Nx;j++){
				Oprt.Om2d[i][j]=fun_exact_disp_val(dominfo.kx[j],dominfo.ky[i],dominfo.bathy.depth,seainfo.mdir,seainfo.current.x,seainfo.current.y);														
				Oprt.L2d[i][j]=Oprt.Om2d[i][j]*Oprt.Om2d[i][j]/grav;
				Oprt.Cpeak2d[i][j]= Oprt.Cpeak;
				Oprt.C2d[i][j] = fun_exact_Cp_val(dominfo.kx[j],dominfo.ky[i],dominfo.bathy.depth,seainfo.current);
			}
		}
	}
	
	for(int i=0;i<dominfo.Ny;i++){
		for(int j=0;j<dominfo.Nx;j++){
			Oprt.Csqr2d[i][j] = pow(fun_exact_Cp1d_val(dominfo.KK[i][j],dominfo.bathy.depth,dominfo.U),2);
		}
	}
	
	return Oprt;
}

void funOprt_DispInterpolation_2p(Oprtvar Oprt,double** bathydat,double** KK,double* kx, double* ky,int Nx,int Ny,double nu,Vec2d u,double mdir)
{
	double D_min,D_max,KK_min,KK_max,k_cut;
	
	D_min=-fun_max2darray(bathydat,Nx,Ny);
	D_max =-fun_min2darray(bathydat,Nx,Ny);
	KK_min=0;
	KK_max=fun_max2darray(KK,Nx,Ny);
	
	int N=101;
	double* DDepth=new double[N];
	double* kk=new double[N];
	double dD=(D_max-D_min)/(N-1);
	double dkk=(KK_max-KK_min)/(N-1);
	for(int i=0;i<N;i++){
		DDepth[i]=D_min+i*dD;
		kk[i]=KK_min+i*dkk;
	}
	//DDepth[0]=DDepth[0]-0.5*dD;
	DDepth[N-1]=DDepth[N-1]+0.1*dD;
	
	double  U=sqrt(pow(u.x,2)+pow(u.y,2));
	double  Cp_min, Cp_plus, Cp_D, Om2_min, Om2_plus, Om2_D,det;
	double* kappanu=new double[N];
	double* Ggam_min=new double[N];
	double* Ggam_plus=new double[N];
	double* Om;
	
	for(int i=0;i<N;i++){
		Om			=fun_exact_disp1d(kk,DDepth[i],N,U);
		kappanu[i]	=fun_invOm(N,Om,kk,nu);
		Cp_D  		=fun_exact_Cp1d_val(kappanu[i], DDepth[i],U);
		Cp_min		=fun_exact_Cp1d_val(kappanu[i], D_min,U);
		Cp_plus		=fun_exact_Cp1d_val(kappanu[i], D_max,U);
		Om2_D		=pow((double) fun_exact_disp1d_val(kappanu[i],DDepth[i],U),2);
		Om2_min		=pow((double) fun_exact_disp1d_val(kappanu[i],D_min,U),2);
		Om2_plus	=pow((double) fun_exact_disp1d_val(kappanu[i],D_max,U),2);
		det			=Cp_min*Om2_plus -Cp_plus*Om2_min;
		Ggam_min[i] = ( Om2_plus*Cp_D-Cp_plus*Om2_D )/det;
		Ggam_plus[i]= ( -Om2_min*Cp_D+Cp_min*Om2_D)/det;
    }
    
	const gsl_interp_type *T = gsl_interp_cspline;
	gsl_interp_accel   *xacc = gsl_interp_accel_alloc();
	gsl_spline *spline_gam_min  = gsl_spline_alloc(T,N);
	gsl_spline *spline_gam_plus = gsl_spline_alloc(T,N);

	// initialize interpolation 
	gsl_spline_init(spline_gam_min,DDepth,Ggam_min, N);
	gsl_spline_init(spline_gam_plus,DDepth,Ggam_plus, N);
	
	double Om2d_min,Om2d_plus;
	double kappa;
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Oprt.InterpD.Gam_min[i][j] =gsl_spline_eval(spline_gam_min, -bathydat[i][j], xacc);
			Oprt.InterpD.Gam_plus[i][j]=gsl_spline_eval(spline_gam_plus, -bathydat[i][j], xacc);
			Om2d_min =fun_exact_disp_val(kx[j],ky[i],D_min,mdir,u.x,u.y);
			Om2d_plus=fun_exact_disp_val(kx[j],ky[i],D_max,mdir,u.x,u.y);

			Oprt.InterpD.L2d_min[i][j] =Om2d_min*Om2d_min/grav;
			Oprt.InterpD.L2d_plus[i][j]=Om2d_plus*Om2d_plus/grav;
			
			Oprt.InterpD.Om2d_min[i][j] =Om2d_min;
			Oprt.InterpD.Om2d_plus[i][j]=Om2d_plus;

			Oprt.InterpD.Om2dSq_min[i][j] =pow(Om2d_min,2);
			Oprt.InterpD.Om2dSq_plus[i][j]=pow(Om2d_plus,2);	
					
			Oprt.InterpD.C2d_min[i][j] =fun_exact_Cp_val(kx[j],ky[i],D_min,u);
			Oprt.InterpD.C2d_plus[i][j]=fun_exact_Cp_val(kx[j],ky[i],D_max,u);
			
			Oprt.InterpD.Csqr2d_min[i][j] = pow(Oprt.InterpD.C2d_min[i][j],2);
			Oprt.InterpD.Csqr2d_plus[i][j]= pow(Oprt.InterpD.C2d_plus[i][j],2);
			
			Om=fun_exact_disp1d(kk,-bathydat[i][j],N,U);
			kappa=fun_invOm(N,Om,kk,nu); 
			Oprt.Cpeak2d[i][j]=fun_exact_Cp1d_val(kappa,-bathydat[i][j],U);
		}
	}
	
	gsl_spline_free(spline_gam_min);
	gsl_spline_free(spline_gam_plus);
	gsl_interp_accel_free(xacc);
	
	delete[] DDepth;
	delete[] kappanu;
	delete[] kk;
	delete[] Om;
	delete[] Ggam_min;
	delete[] Ggam_plus;
}

void funOprt_DispInterpolation_3p(Oprtvar Oprt,double** bathydat,double** KK,double* kx, double* ky,int Nx,int Ny,double nu,Vec2d u,double mdir)
{
	double D_min,D_mid,D_max,KK_min,KK_max,k_cut;	
	
	D_min =-fun_max2darray(bathydat,Nx,Ny);
	D_max =-fun_min2darray(bathydat,Nx,Ny);
	KK_min=0;
	KK_max=fun_max2darray(KK,Nx,Ny); // must be larger than nu
	if (KK_max < nu){
			KK_max = 2*nu;
	}	
    
	if (strcmp(mid_temp,"default")==0){
		D_mid = (D_max+D_min)/2;
	}
	else {
		D_mid = atof(mid_temp);
	}
	
	int N=100;
	double* DDepth=new double[N];
	double* kk=new double[N];
	double dD=(D_max-D_min)/(N-1);
	double dkk=(KK_max-KK_min)/(N-1);
	for(int i=0;i<N;i++){
		DDepth[i]=D_min+i*dD;
		kk[i]=KK_min+i*dkk;
	}
	//DDepth[0]=DDepth[0]-0.5*dD;
	DDepth[N-1]=DDepth[N-1]+0.1*dD;

	double  U=sqrt(pow(u.x,2)+pow(u.y,2));
	double  Cp_min, Cp_mid, Cp_plus, Cp_D, Om2_min, Om2_mid, Om2_plus, Om2_D, det;
	double  C0_min, C0_mid, C0_plus, C0_D;
	double* kappanu=new double[N];
	double* Ggam_min=new double[N];
	double* Ggam_mid=new double[N];
	double* Ggam_plus=new double[N];
	double* Om;
	double  A,B,C,D,E,F,G,H,I;
	
	for(int i=0;i<N;i++){
		Om			=fun_exact_disp1d(kk,DDepth[i],N,U);
		kappanu[i]	=fun_invOm(N,Om,kk,nu);
		Cp_D  		=fun_exact_Cp1d_val(kappanu[i], DDepth[i],U);
		Cp_min		=fun_exact_Cp1d_val(kappanu[i], D_min,U);
		Cp_mid		=fun_exact_Cp1d_val(kappanu[i], D_mid,U);
		Cp_plus		=fun_exact_Cp1d_val(kappanu[i], D_max,U);
		Om2_D		=pow((double) fun_exact_disp1d_val(kappanu[i],DDepth[i],U),2);
		Om2_min		=pow((double) fun_exact_disp1d_val(kappanu[i],D_min,U),2);
		Om2_mid		=pow((double) fun_exact_disp1d_val(kappanu[i],D_mid,U),2);
		Om2_plus	=pow((double) fun_exact_disp1d_val(kappanu[i],D_max,U),2);
		C0_min   = sqrt(grav*D_min);
		C0_mid   = sqrt(grav*D_mid);
		C0_plus  = sqrt(grav*D_max);
		C0_D     = sqrt(grav*DDepth[i]);

		A= (Om2_mid*Cp_plus)-(Om2_plus*Cp_mid);
		B=-(Om2_min*Cp_plus)+(Om2_plus*Cp_min);
		C= (Om2_min*Cp_mid)-(Om2_mid*Cp_min);
		D=-(C0_mid*Cp_plus)+(C0_plus*Cp_mid);
		E= (C0_min*Cp_plus)-(C0_plus*Cp_min);

		F=-(C0_min*Cp_mid)+(C0_mid*Cp_min);
		G= (C0_mid*Om2_plus)-(C0_plus*Om2_mid);
		H=-(C0_min*Om2_plus)+(C0_plus*Om2_min);
		I= (C0_min*Om2_mid)-(C0_mid*Om2_min);

		det=C0_min*A+C0_mid*B+C0_plus*C;

		Ggam_min[i]  = ( A*C0_D+D*Om2_D+G*Cp_D )/det;
		Ggam_mid[i]  = ( B*C0_D+E*Om2_D+H*Cp_D )/det;
		Ggam_plus[i] = ( C*C0_D+F*Om2_D+I*Cp_D )/det;
    }	

	const gsl_interp_type *T = gsl_interp_cspline;
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_spline *spline_gam_min  = gsl_spline_alloc(T,N);
	gsl_spline *spline_gam_mid  = gsl_spline_alloc(T,N);
	gsl_spline *spline_gam_plus = gsl_spline_alloc(T,N);

	// initialize interpolation 
	gsl_spline_init(spline_gam_min,DDepth,Ggam_min, N);
	gsl_spline_init(spline_gam_mid,DDepth,Ggam_mid, N);
	gsl_spline_init(spline_gam_plus,DDepth,Ggam_plus, N);

	double Om2d_min,Om2d_plus,Om2d_mid;
	double kappa;
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			Oprt.InterpD.Gam_min[i][j] =gsl_spline_eval(spline_gam_min, -bathydat[i][j], xacc);
			Oprt.InterpD.Gam_mid[i][j] =gsl_spline_eval(spline_gam_mid, -bathydat[i][j], xacc);
			Oprt.InterpD.Gam_plus[i][j]=gsl_spline_eval(spline_gam_plus, -bathydat[i][j], xacc);
			
			Om2d_min =fun_exact_disp_val(kx[j],ky[i],D_min,mdir,u.x,u.y);
			Om2d_mid =fun_exact_disp_val(kx[j],ky[i],D_mid,mdir,u.x,u.y);
			Om2d_plus=fun_exact_disp_val(kx[j],ky[i],D_max,mdir,u.x,u.y);
			
			Oprt.InterpD.L2d_min[i][j] =Om2d_min*Om2d_min/grav;
			Oprt.InterpD.L2d_mid[i][j] =Om2d_mid*Om2d_mid/grav;
			Oprt.InterpD.L2d_plus[i][j]=Om2d_plus*Om2d_plus/grav;
			
			Oprt.InterpD.Om2d_min[i][j] =Om2d_min;
			Oprt.InterpD.Om2d_mid[i][j] =Om2d_mid;
			Oprt.InterpD.Om2d_plus[i][j]=Om2d_plus;

			Oprt.InterpD.Om2dSq_min[i][j] =pow(Om2d_min,2);
			Oprt.InterpD.Om2dSq_mid[i][j] =pow(Om2d_mid,2);
			Oprt.InterpD.Om2dSq_plus[i][j]=pow(Om2d_plus,2);
			
			Oprt.InterpD.C2d_min[i][j] =fun_exact_Cp_val(kx[j],ky[i],D_min,u);
			Oprt.InterpD.C2d_plus[i][j]=fun_exact_Cp_val(kx[j],ky[i],D_max,u);
			Oprt.InterpD.C2d_mid[i][j] =fun_exact_Cp_val(kx[j],ky[i],D_mid,u);
			
			Oprt.InterpD.Csqr2d_min[i][j] = pow(fun_exact_Cp1d_val(dominfo.KK[i][j],D_min,dominfo.U),2);
			Oprt.InterpD.Csqr2d_mid[i][j] = pow(fun_exact_Cp1d_val(dominfo.KK[i][j],D_mid,dominfo.U),2);
			Oprt.InterpD.Csqr2d_plus[i][j]= pow(fun_exact_Cp1d_val(dominfo.KK[i][j],D_max,dominfo.U),2);
			
			Om=fun_exact_disp1d(kk,-bathydat[i][j],N,U);
			kappa=fun_invOm(N,Om,kk,nu); 
			Oprt.Cpeak2d[i][j]=fun_exact_Cp1d_val(kappa,-bathydat[i][j],U);
		}
	}

	gsl_spline_free(spline_gam_min);
	gsl_spline_free(spline_gam_mid);
	gsl_spline_free(spline_gam_plus);
	gsl_interp_accel_free(xacc);
	delete[] DDepth;
	delete[] kappanu;
	delete[] kk;
	delete[] Om;
	delete[] Ggam_min;
	delete[] Ggam_mid;
	delete[] Ggam_plus;
}

double fun_invOm(int N,double* Om,double* kk,double nu)
{
	const gsl_interp_type *T = gsl_interp_cspline;
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(T,N);
	
	gsl_spline_init(spline,Om,kk,N);	
	
	double result = gsl_spline_eval(spline, nu, xacc);
	
	gsl_spline_free(spline);
	gsl_interp_accel_free(xacc); 
	
	return result;
}

void fun_invOm_vec(double* result,int N,double* Om,double* kk,double* nu,int Nnu)
{
	const gsl_interp_type *T = gsl_interp_cspline;
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(T,N);
	
	gsl_spline_init(spline,Om,kk,N);	

	#pragma omp parallel for
	for (int j=0;j<Nnu;j++){
		if (nu[j]==0){
			result[j] = 0;
		}
		else {
			result[j] = gsl_spline_eval(spline, nu[j], xacc);
		}		
	}
	
	gsl_spline_free(spline);
	gsl_interp_accel_free(xacc); 
	
}

complex fun_exp_ix(double var){
	complex result;
	result.Re=cos(var);
	result.Im=sin(var);
	return result;
}

void freqspace(double* freq_space,double len,int N)
{
	double delta_om_k	= 2*M_PI/len;
	double fmin 		= -0.5*2*M_PI/(len/(N-1));
	double fmax;
	if (((N/2)*2 == N)) {//even N
		fmax = -fmin-delta_om_k;
	}
	else {//odd N
		fmax = -fmin;
	}
	double d_update = (fmax-fmin)/(N-1);
	
	for(int i=0;i<int(N/2);i++){
		freq_space[i]	= fmin+(N/2+i)*d_update;
	}
	int N_i	= 0;
	for(int i=int(N/2);i<N;i++)	{
		freq_space[i]	= fmin+N_i*d_update;
		N_i++;
	}	
	freq_space[0] = 1E-10;
}

void fun_aal2D(double** aal,int cutfrac,double* kx,double* ky,double** kk,int Nx,int Ny)
{
	int idx = ((int)floor((double)Nx/(double)cutfrac)-1);
	int idy = ((int)floor((double)Ny/(double)cutfrac)-1);
	if (idx<0){
		idx=0;
	}
	if (idy<0){
		idy=0;
	}
	double kxmax=fabs(kx[idx]);
	double kymax=fabs(ky[idy]);
	double Kmax;
	
	set_matrix_val(aal,Nx,Ny,1);
	if (idx>0 && idy>0){
		Kmax=sqrt(pow(kxmax,2)+pow(kymax,2));
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (fabs(kk[i][j])>Kmax){
					aal[i][j]=0;	
				}
			}
		}
	}
	else {
		double alx,aly;
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (fabs(kx[j])>kxmax){
					alx=0;	
				}
				else{
					alx=1;
				}
				if (fabs(ky[i])>kymax){
					aly=0;	
				}
				else{
					aly=1;
				}
				aal[i][j]=alx*aly;
			}
		}
	}
}

/*void fun_fbdy2D_old(double** fbdychar,double* x,double* y,int Nx,int Ny,fbdyvar fbdy)
{
	cfSA = declare_2darray(Nx,Ny);
	fun_cfSA2D(cfSA,x,y,Nx,Ny,fbdy);

	//IF COUPLING
	int istart1,istart2;
	double lendamp1,lendamp2,lendamp3,lendamp4;
	if (strcmp(coupling,"yes")==0){
		//in x-axis	
		id_Xinterv1	= int((xinterv1-dominfo.x[0])/dominfo.dx);
		nLin		= int((lenOZ)/dominfo.dx);  
		nCFD 		= int((lcfd)/dominfo.dx);
		id_Xinterv2	= id_Xinterv1+3*nLin+nCFD;
		
		lendamp1 = lenOZ;//(lcfd)/2;
		lendamp2 = lenOZ;		
		istart1 = (xinterv1+2*lenOZ-dominfo.x[0])/dominfo.dx-1;
		istart2 = (xinterv1+4*lenOZ+lcfd-dominfo.x[0])/dominfo.dx-2;		
		
		for(int i=0;i<Ny;i++){//additional damping in the middle (CFD zone)
			for(int j=0;j<Nx;j++){
				if ((x[j]>=x[istart1]) && (x[j]<=(x[istart1]+lendamp1))){
					cfSA[i][j] = cfSA[i][j]*(1+cos(M_PI*(x[j]-x[istart1])/lendamp1))/2;
				}				
				else if ((x[j]>(x[istart2]-lendamp2)) && (x[j]<x[istart2])){
					cfSA[i][j] = cfSA[i][j]*(1-cos(M_PI*(x[j]-x[istart2]+lendamp2)/lendamp2))/2;
				}
				else if ((x[j]>(x[istart1]+lendamp1))&&(x[j]<=(x[istart2]-lendamp2))){
					cfSA[i][j] = 0;
				}
			}
		}
	}
	
	int indxl=fun_closest(x, Nx,(double) x[0]+fbdy.Lleft);
	int indxr=fun_closest(x, Nx,(double) x[Nx-1]-fbdy.Lright);
	int indyb=fun_closest(y, Ny,(double) y[0]+fbdy.Lbottom);
	int indyt=fun_closest(y, Ny,(double) y[Ny-1]-fbdy.Ltop);
	
	double L=7;
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (j<=indxl){
				L=fbdy.Lleft;
			}
			if(j>=indxr){
				L=fbdy.Lright;	
			}
			if(i<=indyb&&(j>indxl&&j<indxr)){
				L=fbdy.Lbottom;
			}
			if(i>=indyt&&(j>indxl&&j<indxr)){
				L=fbdy.Ltop;
			}
			if(L==0){
				L=7;
			}
			fbdychar[i][j]=(1-cfSA[i][j])/(0.9*L);
		}
	}

	//IF COUPLING
	if (strcmp(coupling,"yes")==0){
		for(int i=0;i<Ny;i++){//additional damping in the middle
			for(int j=0;j<Nx;j++){
				if ((x[j]>=(x[istart1]))&&(x[j]<=(x[istart2]-lendamp2))){
					fbdychar[i][j]=(1-cfSA[i][j])/(0.9*lendamp1);
				}
				else if ((x[j]<x[istart2])&&(x[j]>=(x[istart2]-lendamp2))){
					fbdychar[i][j]=(1-cfSA[i][j])/(0.9*lendamp2);
				}
			}
		}
	}	
}*/

void fun_fbdy2D(double** fbdychar,double* x,double* y,int Nx,int Ny,fbdyvar fbdy)
{
	cfSA = declare_2darray(Nx,Ny);
	fun_cfSA2D(cfSA,x,y,Nx,Ny,fbdy);
	
	//IF WALL then add strong damping inside wall
	if ((strcmp(wall,"no")!=0)){
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){			
				cfSA[i][j]=fmin(cfSA[i][j],cfSAwall[i][j]);			
			}
		}
	}
	
	//IF COUPLING
	int istart1,istart2,jstart1,jstart2;
	double lendamp1,lendamp2,lendamp3,lendamp4;
	if (strcmp(coupling,"yes")==0){
		//in x-axis	
		nLin		= (int)((lenOZin)/dominfo.dx);  
		nLtr		= (int)((lenOZtr)/dominfo.dx); 
		
		//in y-axis
		nLiny		= int((lenOZiny)/dominfo.dy);  
		nLtry		= int((lenOZtry)/dominfo.dy); 
		
		//id_Xinterv1	= (int)ceil((xinterv1-dominfo.x[0])/dominfo.dx);		 
		//nCFD 		= (int)((lcfd)/dominfo.dx);		
		//lendamp1 = nCFD*dominfo.dx/2;
		//lendamp2 = nCFD*dominfo.dx/2;		
		//istart1 = id_Xinterv1+2*nLin;
		//istart2 = id_Xinterv1+4*nLin+nCFD;			
		
		//id_Yinterv1	= (int)ceil((yinterv1-dominfo.y[0])/dominfo.dy);		 
		//nCFDy 		= int((lcfd_y)/dominfo.dy);
		//id_Yinterv2	= id_Yinterv1+nCFDy;
		
		/*if ((lenOZy!=0)){
			lendamp3 = nCFDy*dominfo.dy/2;
			lendamp4 = nCFDy*dominfo.dy/2;
		}
		else{
			id_Yinterv2 = 0;
			id_Yinterv2 = dominfo.Ny-1;
			lendamp3 = 0;
			lendamp4 = 0;
		}
		jstart1 = id_Yinterv1+2*nLiny;
		jstart2 = id_Yinterv1+4*nLiny+nCFDy;
		
		//additional damping in the middle (CFD zone)		
		int Nxx,Nyy;		
		Nxx = istart2-istart1+1;
		Nyy = jstart2-jstart1+1;
		double** cfSA_cfd = declare_2darray(Nxx,Nyy);
		double *xx = new double[Nxx];
		double *yy = new double[Nyy];
		fbdyvar fbdy_cfd;
		fbdy_cfd.Lleft   = lendamp1;
		fbdy_cfd.Lright  = lendamp2;
		fbdy_cfd.Lbottom = lendamp3;
		fbdy_cfd.Ltop 	 = lendamp4;
		
		for (int i=0;i<Nxx;i++){
			xx[i] = dominfo.x[istart1+i];
		}
		for (int i=0;i<Nyy;i++){
			yy[i] = dominfo.y[jstart1+i];
		}
		fun_cfSA2D(cfSA_cfd,xx,yy,Nxx,Nyy,fbdy_cfd);	
		
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if ((i>=jstart1)&&(i<=jstart2)&&(j>=istart1)&&(j<=istart2)){
					cfSA[i][j] = cfSA[i][j]*(1-cfSA_cfd[i-jstart1][j-istart1]);
				}
			}
		}	
		
		free_2darray(cfSA_cfd,Nxx,Nyy);
		delete[] xx;
		delete[] yy;*/
	}
	
	int indxl=fun_closest(x, Nx,(double) x[0]+fbdy.Lleft);
	int indxr=fun_closest(x, Nx,(double) x[Nx-1]-fbdy.Lright);
	int indyb=fun_closest(y, Ny,(double) y[0]+fbdy.Lbottom);
	int indyt=fun_closest(y, Ny,(double) y[Ny-1]-fbdy.Ltop);
	
	double L=7;
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if (j<=indxl){
				L=fbdy.Lleft;
			}
			if(j>=indxr){
				L=fbdy.Lright;	
			}
			if(i<indyb&&(j>indxl&&j<indxr)){
				L=fbdy.Lbottom;
			}
			if(i>indyt&&(j>indxl&&j<indxr)){
				L=fbdy.Ltop;
			}
			if(L==0){
				L=7;
			}
			fbdychar[i][j]=(1-cfSA[i][j])/(0.9*L);
		}
	}
	
	if ((strcmp(wall,"no")!=0)&&(strcmp(evolinfo.model,"HS1")!=0)){
		for(int i=0;i<Ny;i++){//additional stronger damping
			for(int j=0;j<Nx;j++){
				fbdychar[i][j]=10*fbdychar[i][j];
			}
		}
	}
	
	/*FILE* gp1=popen("gnuplot -persistent","w");
	fun_plot_2d_trarray(gp1,dominfo.y,dominfo.x,fbdychar,dominfo.Ny,dominfo.Nx);
	fflush(gp1);
	pclose(gp1);*/
	
}

void fun_cfSA2D(double** cfSA,double* x,double* y,int Nx, int Ny,fbdyvar fbdy)
{
	int indxl=fun_closest(x, Nx,(double) x[0]+fbdy.Lleft);
	int indxr=fun_closest(x, Nx,(double) x[Nx-1]-fbdy.Lright);
	int indyb=fun_closest(y, Ny,(double) y[0]+fbdy.Lbottom);
	int indyt=fun_closest(y, Ny,(double) y[Ny-1]-fbdy.Ltop);

	double cfX,cfXL, cfXR, A, B,C,D;
	double* heavXL;
	double* heavXR;
	heavXL=fun_heaviside(x,Nx,x[0]);
	heavXR=fun_heaviside(x,Nx,x[indxr]);
	
	double cfY,cfYB, cfYT, E, F,G,H;
	double* heavYB;
	double* heavYT;
	heavYB=fun_heaviside(y,Ny,y[0]);
	heavYT=fun_heaviside(y,Ny,y[indyt]);

	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			A=(double) fun_sign(x[j]-x[indxl]);
			B=(1-cos((x[j]-x[0])*M_PI/fbdy.Lleft))/2;
			cfXL=heavXL[j]*(A>B?A:B);
			C=(double) fun_sign(x[j]-x[Nx-1]);
			D=(1-cos((x[j]-x[indxr])*M_PI/fbdy.Lright))/2;
			cfXR=heavXR[j]*(C>D?C:D);
			if(fbdy.Lleft==0){cfXL=1;}
			if(fbdy.Lright==0){cfXR=0;}
			cfX=cfXL-cfXR;
			//printf("xi=%lf xl=%lf cf=%lf\n",x[j],x[indxl],cfX[j]);
			E=(double) fun_sign(y[i]-y[indyb]);
			F=(1-cos((y[i]-y[0])*M_PI/fbdy.Lbottom))/2;
			cfYB=heavYB[i]*(E>F?E:F);
			G=(double) fun_sign(y[i]-y[Ny-1]);
			H=(1-cos((y[i]-y[indyt])*M_PI/fbdy.Ltop))/2;
			cfYT=heavYT[i]*(G>H?G:H);
			if(fbdy.Lbottom==0){cfYB=1;}
			if(fbdy.Ltop==0){cfYT=0;}
			cfY=cfYB-cfYT;
			cfSA[i][j]=cfX*cfY;
		}
	}
   
	delete[] heavXL;
	delete[] heavYB;
	delete[] heavXR;
	delete[] heavYT;
}

wavestruc fun_generation_init( int Nx, int Ny)
{	
	wavestruc waveinit;
	waveinit.profile = declare_2darray(Nx,Ny);	
	waveinit.u     	 = declare_2darray(Nx,Ny);
	waveinit.v     	 = declare_2darray(Nx,Ny);	
	
	if (strcmp(initial,"user_defined")==0){
		FILE* finit;
		char* initfile = new char[256];
		sprintf(initfile,"%s/initial.dat",arg1);		
		finit = fopen(initfile, "r");
		if(finit==NULL){
			printf("Can not open %s!!\n",initfile);
			exit(0);
		}
		
		//get eta init
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (!fscanf(finit,"%lf ", &waveinit.profile[i][j])){
					printf("The data size is incorrect!");
					break;
				}				
			}
		}
		printf("Reading initial eta done\n");
		
		//get u init
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (!fscanf(finit,"%lf ", &waveinit.u[i][j])){
					waveinit.u[i][j] = 0.0;					
				}				
			}
		}
		
		//get v init
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				if (!fscanf(finit,"%lf ", &waveinit.v[i][j])){
					waveinit.v[i][j] = 0.0;					
				}				
			}
		}
		printf("Reading initial velocity done\n");
		
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				waveinit.profile[i][j] = waveinit.profile[i][j]*(cfSA[i][j]);
				waveinit.u[i][j]       = waveinit.u[i][j]*(cfSA[i][j]);
				waveinit.v[i][j]       = waveinit.v[i][j]*(cfSA[i][j]);
			}
		}
	}
	else{
		set_matrix_val(waveinit.profile,Nx,Ny,0.0);
		set_matrix_val(waveinit.u,Nx,Ny,0.0);
		set_matrix_val(waveinit.v,Nx,Ny,0.0);
	}
	
	return waveinit;
}



void harmonic_influx(int nlines,int indx0,int indx1,int indy0,int indy1)
{
	wave_influx.eta = declare_2darray(wave_influx.Nxy,dominfo.Nt);
	set_matrix_val(wave_influx.eta,wave_influx.Nxy,dominfo.Nt,0.0);
	double W0=2*M_PI/seainfo.Tp;
	double K0=invers_omega(W0,meandepthinf);
	if (indx0==indx1){//vertical
		#pragma omp parallel for
		for (int i=0;i<dominfo.Nt;i++){
			for (int ii=0;ii<wave_influx.nlines;ii++){
				wave_influx.eta[i][ii+wave_influx.idy0] = ramp2dori[i][ii]*seainfo.Hs/2*cos(K0*sin(seainfo.mdir)*dominfo.y[ii]-W0*dominfo.t[i]);
				//ampl.*cos(K0.*cos(theta0)*(dom.X(indxInfl)-xinit)+K0.*sin(theta0)*dom.Y(indyInfl)-W0.*timesig);
			}
		}
	}
	else {
		#pragma omp parallel for
		for (int i=0;i<dominfo.Nt;i++){
			for (int ii=0;ii<wave_influx.nlines;ii++){
				wave_influx.eta[i][ii+wave_influx.idx0] = ramp2dori[i][ii]*seainfo.Hs/2*cos(K0*cos(seainfo.mdir)*dominfo.x[ii]-W0*dominfo.t[i]);
				//ampl.*cos(K0.*cos(theta0)*dom.X(indxInfl)+K0.*sin(theta0)*(dom.Y(indyInfl)-yinit)-W0.*timesig);
			}
		}
	}
}

double* fun_jonswap(double* omg, double omgp, double* halfomsig, double Tp, double gamJS, double Hs,int Nt)
{
	int 	Nomg 	= Nt/2;                   
	double* sigma 	= (double*) malloc (Nomg*sizeof(double));
	
	for (int i=0;i<Nomg;i++){
		if (omg[i]<=omgp){
			sigma[i] 	= 0.07;
		}
		else {
			sigma[i]	= 0.09;
		}
	}
	double f,gampangkat;
	double* JS		= (double*) malloc(Nomg*sizeof(double));
	double* JPM 	= (double*) malloc(Nomg*sizeof(double));
	double* JS_NOR 	= (double*) malloc(Nomg*sizeof(double));
	double alpha    = 0.0081;
	double f_p      = omgp/(2*M_PI);	
	
	set_init_double(JS,Nomg);
	for (int i=1;i<Nomg;i++){
		f 			= omg[i]/(2*M_PI);
		gampangkat 	= exp(-0.5*pow(((f/f_p-1)/sigma[i]),2));
		JPM[i]      = alpha*pow(grav,2)*pow((2*M_PI),(-4))*pow(f,(-5))*exp(-(5.0/4.0)*pow((f/f_p),(-4))); //Pierson - Moskowitz
		JS[i]       = alpha*pow(grav,2)*pow((2*M_PI),(-4))*pow(f,(-5))*exp(-(5.0/4.0)*pow((f/f_p),(-4)))*pow(gamJS,(gampangkat)); //JONSWAP
	}
	
	// calculating Hs from non-normalized JONSWAP spectrum 
	double m0;
	m0      = trapz(omg,JS,Nomg);  
	
	// Now normalize the spectrum based on definition of zeroth-order moment 
	// and its relation with Hs : mo = int(omg,JS(omega)) and Hs = 4*sqrt(mo); so:
	
	for (int i=0;i<Nomg;i++){
		JS_NOR[i]  = (pow(Hs,2)/16)*JS[i]/(m0);
	}
	free(sigma);
	free(JS);
	free(JPM);	

	//correction or callibration of the spectrum
	double* JS_corr = new double[Nomg];
	char  input_corr[128];
	FILE* fcorr;
	sprintf(input_corr, "%s/SPcorrection.dat", arg1);	
	fcorr	= fopen(input_corr, "r");
	
	if(fcorr==NULL){
		printf("There is no SPcorrection.dat file!\nRunning without any spectrum correction.\n");
		for (int j=0;j<Nomg;j++){
			JS_corr[j] = JS_NOR[j];
		}
	}
	else {
		int nlines = countlines(fcorr);
		double* freq_corr = (double*) malloc(sizeof(double)*nlines);
		double* fact_corr = (double*) malloc(sizeof(double)*nlines);
		double* omg_corr  = (double*) malloc(sizeof(double)*nlines);
		for (int j=0;j<nlines;j++){
			fscanf(fcorr,"%lf %lf\n",&freq_corr[j],&fact_corr[j]);
			omg_corr[j] = freq_corr[j]*2*M_PI;
		}
		fclose(fcorr);

		gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
		gsl_spline *spline 	   = gsl_spline_alloc (gsl_interp_cspline, nlines);
		gsl_spline_init (spline, omg_corr, fact_corr, nlines);
		for (int j=0;j<Nomg;j++){
			if ((omg[j]>=omg_corr[0]) && (omg[j]<=omg_corr[nlines-1])){
				JS_corr[j] = gsl_spline_eval(spline,omg[j],acc)*JS_NOR[j];
			}
			else {
				JS_corr[j] = JS_NOR[j];
			}
		}
		free(freq_corr);
		free(omg_corr);
		free(fact_corr);
	}	
	return JS_NOR;
}

void jonswap_influx(double T0, double Tend, double Tp, double dt, int Nt,double* omsig,double* halfomsig,double gamJS,double Hs,int nlines,int indx0,int indy0)
{
	double   omgp = 2*M_PI/Tp;
	double 	 domg = wave_influx.ww[1]-wave_influx.ww[0];
	
	//get the Jonswap spectra
	double* JS_sp = fun_jonswap(omsig,omgp,halfomsig,Tp,gamJS,Hs,Nt);	
	
	//preparation of directional spreading
	int 	Nhalf = Nt/2;
	double* Rtheta = new double[Nhalf];
	double* cos2pdf= new double[Nhalf];
	double  ds 	   = (double) M_PI/(Nhalf-1);
	double* theta  = new double[Nhalf];
	fun_intervaltovector_ds(theta,-M_PI/2.0,ds,Nhalf);
	double  StdTheta;
	
	if (strcmp(seainfo.s,"inf")==0){//long crested wave
		set_array_val(Rtheta,Nhalf,0.0);
		StdTheta = 0.0;			
		int idx_theta0   = fun_closest(theta,Nhalf,0.0);
		for (int j=0;j<Nhalf;j++){
			if (j==idx_theta0){
				cos2pdf[j] = 1/(theta[1]-theta[0]);
			}
			else {
				cos2pdf[j] = 0;
			}			
		}		
	}
	else {//short crested wave
		double s 	= atof(seainfo.s);
		double coef = (tgamma(s+1.0)/(sqrt(M_PI)*tgamma(s+0.5)));
		for (int j=0;j<Nhalf;j++){
			cos2pdf[j] = coef*pow(cos(theta[j]),2*s);
		}
		
		double* cos2cdf= new double[Nhalf];
		double  total  = 0;
		for (int j=0;j<Nhalf;j++){
			total 		= total+cos2pdf[j];
			cos2cdf[j]	= total;
		}
		
		double* yadd = new double[Nhalf];
		for (int j=0;j<Nhalf;j++){
			yadd[j] = exp(-12.0)/M_PI*theta[j];
		}		
		for (int j=0;j<Nhalf;j++){
			cos2cdf[j]	= cos2cdf[j]/cos2cdf[Nhalf-1]+yadd[j];
		}
		
		//generate random number
		double* RN = new double[Nhalf];
		srand(time(0)); 
		for (int j=0;j<Nhalf;j++){
			RN[j] = (double) rand()/(RAND_MAX);
		}

		//interpolation for random phases
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_spline *spline	  = gsl_spline_alloc (gsl_interp_cspline, Nhalf);				  
		gsl_spline_init (spline,cos2cdf,theta,Nhalf);

		double* Rcos2s = new double[Nhalf];
		for (int j=0;j<Nhalf;j++){
			Rcos2s[j] = gsl_spline_eval(spline,RN[j],acc);
		}
		
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);

		//calculating standard deviation
		double mean_temp    = mean_array(Rcos2s,Nhalf);
		double mean_cos2pdf = mean_temp+seainfo.mdir;

		//returning outside IF
		StdTheta  = sqrt(var_array(Rcos2s,Nhalf,mean_temp));
		for (int j=0;j<Nhalf;j++){
			Rtheta[j] = Rcos2s[j];
		}

		delete[] Rcos2s;
		delete[] cos2cdf;
		delete[] RN;
		delete[] yadd;		
	}

	//building 2d spectrum
	double* a_bar = new double[Nhalf];
	for (int j=0;j<Nhalf;j++){
		a_bar[j] = sqrt(2*JS_sp[j]*domg);
	}

	//generate random phase
	double* RP = new double[Nhalf];
	srand(time(0)); 
	for (int j=0;j<Nhalf;j++){
		RP[j] = 2*M_PI*(double) rand()/(RAND_MAX);
	}
	
	//Vertical or horizontal orientation computing the signal
	double   eta_xy;
	double 	 KK,WW;
	wave_influx.eta = declare_2darray(wave_influx.Nxy,dominfo.Nt);
	set_matrix_val(wave_influx.eta,wave_influx.Nxy,dominfo.Nt,0.0);
	
	if (strcmp(orientation,"vertical")==0){
		for (int i=0;i<dominfo.Nt;i++){
			for (int ii=0;ii<wave_influx.nlines;ii++){
				eta_xy = 0;
				for (int iii=0;iii<Nhalf;iii++){
					KK      = wave_influx.kw[iii];
					WW      = wave_influx.ww[iii];					
					eta_xy += a_bar[iii]*cos(-WW*dominfo.t[i]+RP[iii]+KK*cos(Rtheta[iii]+seainfo.mdir)*(dominfo.x[indx0]-dominfo.x[indx0])
							  +KK*sin(Rtheta[iii]+seainfo.mdir)*(wave_influx.y[ii]-dominfo.y[indy0]));					
				}
				wave_influx.eta[i][ii+wave_influx.idy0]= eta_xy*ramp2dori[i][ii];
			}
		}
	}
	else {
		for (int i=0;i<dominfo.Nt;i++){
			for (int ii=0;ii<nlines;ii++){
				eta_xy = 0;
				for (int iii=0;iii<Nhalf;iii++){
					KK      = wave_influx.kw[iii];
					WW      = wave_influx.ww[iii];
					eta_xy += a_bar[iii]*cos(KK*sin(Rtheta[iii]+seainfo.mdir)*(dominfo.y[indy0]-dominfo.y[indy0])
							 +KK*cos(Rtheta[iii]+seainfo.mdir)*(wave_influx.x[ii]-dominfo.x[indx0])
							 -WW*dominfo.t[i]+RP[iii]);					
				}
				wave_influx.eta[i][ii+wave_influx.idx0]= eta_xy*ramp2dori[i][ii];
			}
		}
	}
}

void fun_tapered2d(const double* xy,double xmin,double xmax,double LrampT,double LrampXY,const double dxy)
{
	int idl = round((xmin+LrampXY-xy[0])/dxy);
	int idr = round((xmax-LrampXY-xy[0])/dxy);
	int idb = round((LrampT)/dominfo.dt);
	int idt = round((dominfo.t[dominfo.Nt-1]-LrampT-dominfo.t[0])/dominfo.dt);
	
	double cfXl,cfXr,cfTb,cfTt;
	double** cfX = declare_2darray(wave_influx.Nxy,dominfo.Nt);
	double** cfT = declare_2darray(wave_influx.Nxy,dominfo.Nt);

	if (taper==0){
		set_matrix_val(cfX,wave_influx.Nxy,dominfo.Nt,1.0);
	}
	else {	
		//#pragma omp parallel for
		for (int j=0;j<wave_influx.Nxy;j++){
			cfXl = fun_heav(xy[j]-xmin)*GSL_MAX(fun_sign(xy[j]-xy[idl]),0.5*(1.0-cos((xy[j]-xmin)*M_PI/LrampXY)));
			cfXr = fun_heav(xy[j]-xy[idr])*GSL_MAX(fun_sign(xy[j]-xmax),0.5*(1.0-cos((xy[j]-xy[idr])*M_PI/LrampXY)));
			for (int i=0;i<dominfo.Nt;i++){
				cfX[i][j]= cfXl-cfXr;
			}
		}
	}
	
	if (ramp==0){
		set_matrix_val(cfT,wave_influx.Nxy,dominfo.Nt,1.0);
	}
	else {	
		#pragma omp parallel for
		for (int i=0;i<dominfo.Nt;i++){
			cfTb = fun_heav(dominfo.t[i]-dominfo.t[0])*GSL_MAX(fun_sign(dominfo.t[i]-dominfo.t[idb]),0.5*(1-cos((dominfo.t[i]-dominfo.t[0])*M_PI/LrampT)));
			cfTt = fun_heav(dominfo.t[i]-dominfo.t[idt])*GSL_MAX(fun_sign(dominfo.t[i]-dominfo.t[dominfo.Nt-1]),0.5*(1-cos((dominfo.t[i]-dominfo.t[idt])*M_PI/LrampT)));
			for (int j=0;j<wave_influx.Nxy;j++){
				cfT[i][j]= cfTb-cfTt;
			}
		}
	}
	
	#pragma omp parallel for
	for (int i=0;i<dominfo.Nt;i++){
		for (int j=0;j<wave_influx.Nxy;j++){
			ramp2d[i][j]= cfX[i][j]*cfT[i][j];
		}
	}
	free_2darray(cfX,wave_influx.Nxy,dominfo.Nt);
	free_2darray(cfT,wave_influx.Nxy,dominfo.Nt);
}

void fun_ramp2d(double ramp_coef,double tapered_coef)
{
	ramp2d = declare_2darray(wave_influx.Nxy,dominfo.Nt);
	if ((ramp_coef==0.0) && (tapered_coef==0.0)){
		set_matrix_val(ramp2d,wave_influx.Nxy,dominfo.Nt,1.0);
	}
	else {	
		double LrampT = ramp_coef*(seainfo.Tp);		
		if (strcmp(orientation,"vertical")==0){
			double ymin=dominfo.y[0];
			double ymax=dominfo.y[dominfo.Ny-1];
			double Linf=ymax-ymin;
			double LrampXY=tapered_coef/2.0*Linf;
			fun_tapered2d(dominfo.y,ymin,ymax,LrampT,LrampXY,dominfo.dy);
		}
		else {
			double xmin=dominfo.x[0];
			double xmax=dominfo.x[dominfo.Nx-1];
			double Linf=xmax-xmin;
			double LrampXY=tapered_coef*Linf;
			fun_tapered2d(dominfo.x,xmin,xmax,LrampT,LrampXY,dominfo.dx);
		}
	}
}

void fun_tapered2dori(const double* xy,double xmin,double xmax,double LrampT,double LrampXY,const double dxy)
{
	int idl = round((xmin+LrampXY-xy[0])/dxy);
	int idr = round((xmax-LrampXY-xy[0])/dxy);
	int idb = round((LrampT)/dominfo.dt);
	int idt = round((dominfo.t[dominfo.Nt-1]-LrampT-dominfo.t[0])/dominfo.dt);
	
	double cfXl,cfXr,cfTb,cfTt;
	double** cfX = declare_2darray(wave_influx.nlines,dominfo.Nt);
	double** cfT = declare_2darray(wave_influx.nlines,dominfo.Nt);

	if (taper==0){
		set_matrix_val(cfX,wave_influx.nlines,dominfo.Nt,1.0);
	}
	else {	
		//#pragma omp parallel for
		for (int j=0;j<wave_influx.nlines;j++){
			cfXl = fun_heav(xy[j]-xmin)*GSL_MAX(fun_sign(xy[j]-xy[idl]),0.5*(1.0-cos((xy[j]-xmin)*M_PI/LrampXY)));
			cfXr = fun_heav(xy[j]-xy[idr])*GSL_MAX(fun_sign(xy[j]-xmax),0.5*(1.0-cos((xy[j]-xy[idr])*M_PI/LrampXY)));
			for (int i=0;i<dominfo.Nt;i++){
				cfX[i][j]= cfXl-cfXr;
			}
		}
	}
	
	if (ramp==0){
		set_matrix_val(cfT,wave_influx.nlines,dominfo.Nt,1.0);
	}
	else {	
		#pragma omp parallel for
		for (int i=0;i<dominfo.Nt;i++){
			cfTb = fun_heav(dominfo.t[i]-dominfo.t[0])*GSL_MAX(fun_sign(dominfo.t[i]-dominfo.t[idb]),0.5*(1-cos((dominfo.t[i]-dominfo.t[0])*M_PI/LrampT)));
			cfTt = fun_heav(dominfo.t[i]-dominfo.t[idt])*GSL_MAX(fun_sign(dominfo.t[i]-dominfo.t[dominfo.Nt-1]),0.5*(1-cos((dominfo.t[i]-dominfo.t[idt])*M_PI/LrampT)));
			for (int j=0;j<wave_influx.nlines;j++){
				cfT[i][j]= cfTb-cfTt;
			}
		}
	}
	
	#pragma omp parallel for
	for (int i=0;i<dominfo.Nt;i++){
		for (int j=0;j<wave_influx.nlines;j++){
			ramp2dori[i][j]= cfX[i][j]*cfT[i][j];
		}
	}
	free_2darray(cfX,wave_influx.nlines,dominfo.Nt);
	free_2darray(cfT,wave_influx.nlines,dominfo.Nt);
}


void fun_ramp2dori(double ramp_coef,double tapered_coef)
{
	ramp2dori = declare_2darray(wave_influx.nlines,dominfo.Nt);
	if ((ramp_coef==0.0) && (tapered_coef==0.0)){
		set_matrix_val(ramp2dori,wave_influx.nlines,dominfo.Nt,1.0);
	}
	else {	
		double LrampT = ramp_coef*(seainfo.Tp);		
		if (strcmp(orientation,"vertical")==0){
			double ymin=wave_influx.y[0];
			double ymax=wave_influx.y[wave_influx.nlines-1];
			double Linf=ymax-ymin;
			double LrampXY=tapered_coef/2.0*Linf;
			fun_tapered2dori(wave_influx.y,ymin,ymax,LrampT,LrampXY,dominfo.dy);
		}
		else {
			double xmin=wave_influx.x[0];
			double xmax=wave_influx.x[wave_influx.nlines-1];
			double Linf=xmax-xmin;
			double LrampXY=tapered_coef*Linf;
			fun_tapered2dori(wave_influx.x,xmin,xmax,LrampT,LrampXY,dominfo.dx);
		}
	}
}

int read_row(const char *file_name)
{
    FILE *myfile = fopen(file_name, "r");
    int newRows = 0; 
    int ch;

	while(!feof(myfile)){
		ch = fgetc(myfile);
		if(ch == '\n') {
			newRows++;
		}
	}    
	fclose(myfile);
	
    return newRows+1; 
}

int read_col(const char *file_name)
{
    FILE *myfile = fopen(file_name, "r");
    int newCols = 0; 
    int ch;

    while(!feof(myfile)){
		ch = fgetc(myfile);
		if ((ch == ' ')){
			newCols++;
		}
		else if(ch == '\n') {
			break;
		}
	}
	fclose(myfile);

	return newCols; 
}

void fun_generation_influx(int Nt)
{
	if (strcmp(orientation,"vertical")==0){
		wave_influx.Nxy = dominfo.Ny;
	}
	else{
		wave_influx.Nxy = dominfo.Nx;
	}
	
	//user_defined influxing
	if (strcmp(wavename,"user_defined")==0){
		FILE* finflux;
		char influxfile[256];
		sprintf(influxfile,"%s/influx.dat",arg1);		
		finflux = fopen(influxfile, "r");
		if(finflux==NULL){
			printf("Can not open %s!!\n",influxfile);
			exit(0);
		}
		
		//reading influx data
		int row,col;
		row=dominfo.Nt+2;
		col=read_col(influxfile);
		if (col>wave_influx.Nxy){
			col = wave_influx.Nxy;
		}
		
		int 	 nlines = col;//countlines(finflux);belum bener itung kolomnya
		wave_influx.nlines = nlines;
		double*  x 		= new double[nlines];
		double*  y 		= new double[nlines];
		double** eta	= declare_2darray(nlines,Nt);
		double trash;
		fscanf(finflux,"%lf ", &trash);
		for(int i=0;i<nlines;i++){
			fscanf(finflux,"%lf ", &x[i]);
		}
		
		fscanf(finflux,"%lf ", &trash);
		for(int i=0;i<nlines;i++){
			fscanf(finflux,"%lf ", &y[i]);
		}	
		
		for(int j=0;j<dominfo.Nt;j++){
			fscanf(finflux,"%lf ",&dominfo.t[j]);
			//printf("%f\n",dominfo.t[j]);
			for(int i=0;i<nlines;i++){
				fscanf(finflux,"%lf ", &eta[j][i]);
			}
		}
		fclose(finflux);		
	
		//set influxing based on orientation		
		wave_influx.eta = declare_2darray(wave_influx.Nxy,Nt);
		set_matrix_val(wave_influx.eta,wave_influx.Nxy,Nt,0.0);
		wave_influx.x   = new double[nlines];
		wave_influx.y   = new double[nlines];
		double *bath_inf= new double[nlines];

		int indx[nlines],indy[nlines];
		if (strcmp(orientation,"vertical")==0){		
			for (int j=0;j<nlines;j++){
				indx[j] = round((x[j]-dominfo.x[0])/dominfo.dx);
				indy[j] = round((y[j]-dominfo.y[0])/dominfo.dy);
				bath_inf[j]      = dominfo.bathy.data[indy[j]][indx[j]];	
				wave_influx.x[j] = dominfo.x[indx[j]];
				wave_influx.y[j] = dominfo.y[indy[j]];
			}
		}
		else {		
			for (int j=0;j<nlines;j++){
				indx[j] = round((x[j]-dominfo.x[0])/dominfo.dx);
				indy[j] = round((y[j]-dominfo.y[0])/dominfo.dy);
				bath_inf[j]      = dominfo.bathy.data[indy[j]][indx[j]];	
				wave_influx.x[j] = dominfo.x[indx[j]];
				wave_influx.y[j] = dominfo.y[indy[j]];
			}
		}
		
		//add ramp for influxline
		fun_ramp2dori(ramp,taper);
		
		if (strcmp(orientation,"vertical")==0){		
			for (int j=0;j<nlines;j++){
				for (int i=0;i<Nt;i++){
					wave_influx.eta[i][indy[j]] = eta[i][j]*ramp2dori[i][j];
				}
			}
		}
		else {		
			for (int j=0;j<nlines;j++){
				for (int i=0;i<Nt;i++){
					wave_influx.eta[i][indx[j]] = eta[i][j]*ramp2dori[i][j];
				}
			}
		}
		
		meandepthinf = -mean_array(bath_inf,nlines);		
		
		delete[] bath_inf;
		delete[] x;
		delete[] y;
		free_2darray(eta,nlines,Nt);		
	}
	
	//jonswap or harmonic
	int indx0,indx1,indy0,indy1,nlines;
	if ((strcmp(wavename,"jonswap")==0)||(strcmp(wavename,"harmonic")==0)){
		indx0  = round((lineinflux.xstart-dominfo.x[0])/dominfo.dx);
		indy0  = round((lineinflux.ystart-dominfo.y[0])/dominfo.dy);
		indx1  = round((lineinflux.xend-dominfo.x[0])/dominfo.dx);
		indy1  = round((lineinflux.yend-dominfo.y[0])/dominfo.dy);
		if (strcmp(orientation,"vertical")==0){
			nlines = (indy1-indy0)+1;
		}
		else {
			nlines = (indx1-indx0)+1;
		}
		
		double *bath_inf= new double[nlines];
		wave_influx.x   = new double[nlines];
		wave_influx.y 	= new double[nlines];
		wave_influx.nlines = nlines;
		wave_influx.idx0 = indx0;
		wave_influx.idy0 = indy0;

		double mgrad;
		int    idtemp;
		if (strcmp(orientation,"vertical")==0){
			if (indx1!=indx0){
				mgrad = (dominfo.y[indy1]-dominfo.y[indy0])/(dominfo.x[indx1]-dominfo.x[indx0]);
			}
			else {
				mgrad = 0.0;
			}
			wave_influx.Nxy = dominfo.Ny;
			for (int j=0;j<nlines;j++){
				idtemp = round(((dominfo.x[indx0]+j*mgrad/(dominfo.dy))-dominfo.x[0])/dominfo.dx);
				wave_influx.x[j] = dominfo.x[idtemp];
				wave_influx.y[j] = dominfo.y[indy0+j];
				bath_inf[j]  	 = dominfo.bathy.data[indy0+j][idtemp];
			}
		}
		else {
			if (indy1!=indy0){
				mgrad = (dominfo.y[indy1]-dominfo.y[indy0])/(dominfo.x[indx1]-dominfo.x[indx0]);
			}
			else{
				mgrad = 0;
			}
			wave_influx.Nxy = dominfo.Nx;
			for (int j=0;j<nlines;j++){
				idtemp = round(((dominfo.y[indy0]+j*mgrad*(dominfo.dx))-dominfo.y[0])/dominfo.dy);
				wave_influx.x[j] = dominfo.x[indx0+j];
				wave_influx.y[j] = dominfo.y[idtemp];
				bath_inf[j]  	 = dominfo.bathy.data[idtemp][indx0+j];
			}
		}		
		meandepthinf = -mean_array(bath_inf,nlines);
		seainfo.wp   = 2*M_PI/seainfo.Tp;
		seainfo.kp   = invers_omega(seainfo.wp,meandepthinf);
		fun_ramp2dori(ramp,taper);
	}
	
	//for all influx cases
	wave_influx.ww 		= new double[Nt];
	freqspace(wave_influx.ww,dominfo.t[Nt-1]-dominfo.t[0],Nt);
	
	double wmax 		= fun_max(wave_influx.ww,Nt);
	wave_influx.wmax 	= wmax;
	wave_influx.kmax 	= invers_omega(wave_influx.wmax,meandepthinf);
	wave_influx.kw      = new double[dominfo.Nt];
	if (strcmp(bath,"flat")==0){
		for (int j=0;j<dominfo.Nt;j++){
			wave_influx.kw[j]=invers_omega(wave_influx.ww[j],meandepthinf);
		}
	}
	else{		
		int ntemp = 10000;
		double dktemp  = (2.2*wave_influx.kmax)/(ntemp-1);
		double*  ktemp = new double[ntemp];
		fun_intervaltovector_ds(ktemp,-1.1*wave_influx.kmax,dktemp,ntemp);
		double*  wtemp = fun_exact_disp1d(ktemp,meandepthinf,ntemp,dominfo.U);		
		fun_invOm_vec(wave_influx.kw,ntemp,wtemp,ktemp,wave_influx.ww,dominfo.Nt);//k(omega) from time
		delete[] wtemp;
		delete[] ktemp;
	}
	
	if (strcmp(wavename,"jonswap")==0){
		double* halfomsig = (double*) malloc(Nt/2*sizeof(double)); 
		for (int i=1;i<Nt/2;i++){
			halfomsig[i] = wave_influx.ww[i];
		}
		jonswap_influx(dominfo.t[0],dominfo.t[Nt-1],seainfo.Tp,dominfo.dt,Nt,wave_influx.ww,halfomsig,seainfo.gamJS,seainfo.Hs,nlines,indx0,indy0);
		free(halfomsig);
		FILE* fin=fopen(strcat(arg1,"influx.dat"),"w");
		fprintf(fin,"0 ");
		for(int i=0;i<wave_influx.Nxy;i++){
			fprintf(fin,"%lf ", wave_influx.x[i]);
		}fprintf(fin,"\n");
		
		fprintf(fin,"0 ");
		for(int i=0;i<wave_influx.Nxy;i++){
			fprintf(fin,"%lf ", wave_influx.y[i]);
		}fprintf(fin,"\n");
		
		if (strcmp(orientation,"vertical")==0){
			for(int j=0;j<Nt;j++){
				fprintf(fin,"%lf ", dominfo.t[j]);
				for(int i=0;i<wave_influx.Nxy;i++){
					fprintf(fin,"%lf ",wave_influx.eta[j][i+wave_influx.idy0]);
				}
				if (j<Nt-1){fprintf(fin,"\n");}			
			}	
		}
		else{
			for(int j=0;j<Nt;j++){
				fprintf(fin,"%lf ", dominfo.t[j]);
				for(int i=0;i<wave_influx.Nxy;i++){
					fprintf(fin,"%lf ",wave_influx.eta[j][i+wave_influx.idx0]);
				}
				if (j<Nt-1){fprintf(fin,"\n");}			
			}	
		}
		fclose(fin);
	}
	else if (strcmp(wavename,"harmonic")==0){
		harmonic_influx(nlines,indx0,indx1,indy0,indy1);
		FILE* fin=fopen(strcat(arg1,"influx.dat"),"w");
		fprintf(fin,"0 ");
		for(int i=0;i<wave_influx.Nxy;i++){
			fprintf(fin,"%lf ", wave_influx.x[i]);
		}fprintf(fin,"\n");
		
		fprintf(fin,"0 ");
		for(int i=0;i<wave_influx.Nxy;i++){
			fprintf(fin,"%lf ", wave_influx.y[i]);
		}fprintf(fin,"\n");

		if (strcmp(orientation,"vertical")==0){
			for(int j=0;j<Nt;j++){
				fprintf(fin,"%lf ", dominfo.t[j]);
				for(int i=0;i<wave_influx.Nxy;i++){
					fprintf(fin,"%lf ",wave_influx.eta[j][i+wave_influx.idy0]);
				}
				if (j<Nt-1){fprintf(fin,"\n");}			
			}	
		}
		else{
			for(int j=0;j<Nt;j++){
				fprintf(fin,"%lf ", dominfo.t[j]);
				for(int i=0;i<wave_influx.Nxy;i++){
					fprintf(fin,"%lf ",wave_influx.eta[j][i+wave_influx.idx0]);
				}
				if (j<Nt-1){fprintf(fin,"\n");}			
			}	
		}
		fclose(fin);
	}	
	
}

void plotting_preproc(double** Sig,double** gam)
{
	#pragma omp parallel sections 
	{	//initial
		#pragma omp section
		if (strcmp(initial,"zero")!=0){
			FILE* gp00=popen("gnuplot ","w");
			fprintf(gp00,"set term png\n");
			fprintf(gp00,"set output 'Initial.png'\n");
			fun_plot_2d_trarray(gp00,dominfo.y,dominfo.x,waveinit.profile,dominfo.Ny,dominfo.Nx);
			fflush(gp00);
			pclose(gp00);
		}
		
		//damping
		#pragma omp section
		{
			FILE* gp0=popen("gnuplot ","w");
			fprintf(gp0,"set term png\n");
			fprintf(gp0,"set output 'Damping.png'\n");
			fun_plot_2d_trarray(gp0,dominfo.y,dominfo.x,dominfo.fbdy.charac,dominfo.Ny,dominfo.Nx);
			fflush(gp0);
			pclose(gp0);
		}
		
		//influx signal
		#pragma omp section
		{
			FILE* gp=popen("gnuplot ","w");
			fprintf(gp,"set term png\n");
			fprintf(gp,"set output 'Influx_mod.png'\n");
			if (strcmp(orientation,"vertical")==0){
				fun_plot_2d_trarray(gp,dominfo.y,dominfo.t,Sig,dominfo.Ny,dominfo.Nt);
			}
			else {
				fun_plot_2d_array(gp,dominfo.x,dominfo.t,Sig,dominfo.Nx,dominfo.Nt);
			}
			fflush(gp);
			pclose(gp);
		}
			
		//spatial domain
		#pragma omp section
		{
			double** gam_damp = declare_2darray(dominfo.Nx,dominfo.Ny);
			for (int i=0;i<dominfo.Ny;i++){
				for(int j=0;j<dominfo.Nx;j++){
					gam_damp[i][j]=gam[j][i];
				}
			}
			FILE* gp1=popen("gnuplot ","w");
			fprintf(gp1,"set term png\n");
			fprintf(gp1,"set output 'Gamma.png'\n");
			fun_plot_2d_trarray(gp1,dominfo.y,dominfo.x,gam_damp,dominfo.Ny,dominfo.Nx);
			fflush(gp1);
			pclose(gp1);
			free_2darray(gam_damp,dominfo.Nx,dominfo.Ny);
		}
		
		//if varying bottom
		#pragma omp section		
		if (strcmp(dominfo.bathy.Id,"user_defined")==0){
			FILE* gp2=popen("gnuplot -persistent","w");
			fprintf(gp2,"set term png\n");
			fprintf(gp2,"set output 'Bath.png'\n");
			fun_plot_2d_trarray(gp2,dominfo.y,dominfo.x,dominfo.bathy.data,dominfo.Ny,dominfo.Nx);
			fflush(gp2);
			pclose(gp2);
		}
		
		//if wall
		#pragma omp section
		if (strcmp(wall,"no")!=0){
			FILE* gp3=popen("gnuplot ","w");
			fprintf(gp3,"set term png\n");
			fprintf(gp3,"set output 'Wall.png'\n");
			fun_plot_2d_trarray(gp3,dominfo.y,dominfo.x,dominfo.wallchar,dominfo.Ny,dominfo.Nx);
			fflush(gp3);
			pclose(gp3);
		}
		
		//if friction
		#pragma omp section
		if (strcmp(friction,"no")!=0){
			FILE* gp4=popen("gnuplot ","w");
			fprintf(gp4,"set term png\n");
			fprintf(gp4,"set output 'Friction.png'\n");
			fun_plot_2d_trarray(gp4,dominfo.y,dominfo.x,fric_coef,dominfo.Ny,dominfo.Nx);
			fflush(gp4);
			pclose(gp4);
		}
	}
}

void fun_area_generation()
{
	complex** GamKK1_hat;
	double**  UgKK_hat;
	double 	  ffact;
	
	newGam = declare_2darray(dominfo.Ny,dominfo.Nx);
	SkewGam= declare_2darray(dominfo.Ny,dominfo.Nx);
	
	double** Gam     = declare_2darray(dominfo.Ny,dominfo.Nx);
	double** gamheav = declare_2darray(dominfo.Ny,dominfo.Nx);
	set_matrix_val(Gam,dominfo.Ny,dominfo.Nx,0.);
	
	double fact = pow((2*M_PI),2)/(dominfo.dx*dominfo.dy);
	
	//this is for invers omega interpolation
	double*  K1 = new double[dominfo.Nt];
	for (int j=0;j<dominfo.Nt;j++){
		K1[j] = wave_influx.kw[j];
	}
	
	int Nxy = wave_influx.Nxy;
	double dxy;
	double*  k1;
	double*  kk;
	double** KK;

	if (strcmp(orientation,"vertical")==0){
		dxy = dominfo.dy;		
		double     ymin = wave_influx.y[0];
		double 	   ymax = wave_influx.y[wave_influx.nlines-1];
		double** heav1	= fun_Heav2d(dominfo.y,dominfo.Ny,dominfo.fbdy.Lbottom,dominfo.fbdy.Ltop,dominfo.Nx);
		double   Yinfb 	= ymin-dominfo.y[0];
		double   Yinft  = dominfo.y[dominfo.Ny-1]-ymax;
		double** heav2	= fun_Heav2d(dominfo.y,dominfo.Ny,Yinfb,Yinft,dominfo.Nx);
		
		kk  = new double[dominfo.Nx];
		k1  = new double[Nxy];				
		KK  = declare_2darray(Nxy,dominfo.Nt);
		for (int i=0;i<Nxy;i++){
			k1[i] = dominfo.ky[i];
			for (int j=0;j<dominfo.Nt;j++){
				KK[j][i] = sqrt(pow(k1[i],2)+pow(K1[j],2));
			}
		}
		for (int j=0;j<dominfo.Nx;j++){
			kk[j] = dominfo.kx[j];
		}		
		
		int indyinf;
		double dist;	
		
		//prepare Gam and SkewGam (Ny*Nx)
		double*  Ugx	 = new double[dominfo.Nx];
		double*  Ugx_hat = new double[dominfo.Nx];
        complex* temp	 = new complex[dominfo.Nx];
        complex* temp2   = new complex[dominfo.Nx];
        double** dxGam   = declare_2darray(dominfo.Nx,dominfo.Ny);
        set_matrix_val(dxGam,dominfo.Nx,dominfo.Ny,0.0);
        int ind;
        
		for (int jj=0;jj<wave_influx.nlines;jj++){
			ind 	= round((wave_influx.y[jj]-dominfo.y[0])/dominfo.dy);	
			dist	= wave_influx.x[jj]-dominfo.x[0];
			
			fun_exact_Ug(Ugx_hat,kk,meandepthinf,dominfo.Nx,dominfo.U);
			
			ifft_1d_real_real(Ugx,Ugx_hat,dominfo.Nx);//size Nx
			circshift(Ugx,round((double) dist/dominfo.dx),dominfo.Nx);
			
			fft_1d(temp,Ugx,dominfo.Nx);
            for (int i=0;i<dominfo.Nx;i++){
				Gam[i][ind] = Ugx[i];
				temp2[i].Re = -kk[i]*temp[i].Im;
				temp2[i].Im = kk[i]*temp[i].Re;
			}		
			ifft_1d_com_re(temp2,dxGam[ind],dominfo.Nx);//size Nx
        }        
        delete[] temp; 
        delete[] temp2;
        delete[] Ugx;
        delete[] Ugx_hat;
		
		double** dxgamheav	= declare_2darray(dominfo.Ny,dominfo.Nx);
        for (int i=0;i<dominfo.Nx;i++){
			for (int j=0;j<dominfo.Ny;j++) {
				gamheav[i][j]   = fact*Gam[i][j]*heav1[i][j]*heav2[i][j];
				dxgamheav[i][j] = fact*dxGam[j][i]*heav1[i][j]*heav2[i][j];
			}
		}
		free_2darray(heav1,dominfo.Ny,dominfo.Nx);
		free_2darray(heav2,dominfo.Ny,dominfo.Nx);
		free_2darray(dxGam,dominfo.Nx,dominfo.Ny);
		
		complex** Gam_hat    = declare_2darray_complex(dominfo.Ny,dominfo.Nx);
		complex** dxGam_hat  = declare_2darray_complex(dominfo.Ny,dominfo.Nx);
		fftw_2d_r2c(Gam_hat,gamheav,dominfo.Ny,dominfo.Nx);
		fftw_2d_r2c(dxGam_hat,dxgamheav,dominfo.Ny,dominfo.Nx);
		
		ifftw_2d_c2r(newGam,Gam_hat,dominfo.Ny,dominfo.Nx);
		ifftw_2d_c2r(SkewGam,dxGam_hat,dominfo.Ny,dominfo.Nx);
		
		free_2darray(dxgamheav,dominfo.Ny,dominfo.Nx);
		free_2darray_complex(Gam_hat,dominfo.Ny,dominfo.Nx);
		free_2darray_complex(dxGam_hat,dominfo.Ny,dominfo.Nx);
	}
	else {//horizontal influx line
		dxy = dominfo.dx;		
		double     xmin = wave_influx.x[0];
		double 	   xmax = wave_influx.x[wave_influx.nlines-1];
		double** heav1	= fun_Heav2d(dominfo.x,dominfo.Nx,dominfo.fbdy.Lleft,dominfo.fbdy.Lright,dominfo.Ny);
		double   Xinfl 	= xmin-dominfo.x[0];
		double   Xinfr  = dominfo.x[dominfo.Nx-1]-xmax;
		double** heav2	= fun_Heav2d(dominfo.x,dominfo.Nx,Xinfl,Xinfr,dominfo.Ny);
		
		kk  = new double[dominfo.Ny];
		k1  = new double[Nxy];				
		KK  = declare_2darray(Nxy,dominfo.Nt);
		for (int i=0;i<Nxy;i++){
			k1[i] = dominfo.kx[i];
			for (int j=0;j<dominfo.Nt;j++){
				KK[j][i] = sqrt(pow(k1[i],2)+pow(K1[j],2));
			}
		}
		for (int j=0;j<dominfo.Ny;j++){
			kk[j] = dominfo.ky[j];
		}		
		
		int indxinf;
		double dist;

		//prepare Gam and SkewGam (Ny*Nx)
		double*  Ugy	 = new double[dominfo.Ny];
		double*  Ugy_hat = new double[dominfo.Ny];
        complex* temp	 = new complex[dominfo.Ny];
        complex* temp2   = new complex[dominfo.Ny];
        double** dyGam   = declare_2darray(dominfo.Ny,dominfo.Nx);
        set_matrix_val(dyGam,dominfo.Ny,dominfo.Nx,0.0);
        int ind;
        
		for (int jj=0;jj<wave_influx.nlines;jj++){
			ind 	= round((wave_influx.x[jj]-dominfo.x[0])/dxy);		
			dist	= wave_influx.y[jj]-dominfo.y[0];
			
			fun_exact_Ug(Ugy_hat,kk,meandepthinf,dominfo.Ny,dominfo.U);
			ifft_1d_real_real(Ugy,Ugy_hat,dominfo.Ny);	
			circshift(Ugy,round((double) dist/dominfo.dy),dominfo.Ny);
			
			fft_1d(temp,Ugy,dominfo.Ny);
            for (int i=0;i<dominfo.Ny;i++){
				Gam[ind][i] = Ugy[i];
				temp2[i].Re = -kk[i]*temp[i].Im;
				temp2[i].Im = kk[i]*temp[i].Re;
			}		
			ifft_1d_com_re(temp2,dyGam[ind],dominfo.Ny);
        }        
        delete[] temp; delete[] temp2;
        delete[] Ugy;
        delete[] Ugy_hat;
		
		double** dygamheav	= declare_2darray(dominfo.Ny,dominfo.Nx);
        for (int i=0;i<dominfo.Nx;i++){
			for (int j=0;j<dominfo.Ny;j++) {
				gamheav[i][j]   = fact*Gam[i][j]*heav1[j][i]*heav2[j][i];
				dygamheav[i][j] = fact*dyGam[i][j]*heav1[j][i]*heav2[j][i];
			}
		}
		free_2darray(heav1,dominfo.Nx,dominfo.Ny);
		free_2darray(heav2,dominfo.Nx,dominfo.Ny);
		free_2darray(dyGam,dominfo.Ny,dominfo.Nx);
		
		complex** Gam_hat    = declare_2darray_complex(dominfo.Ny,dominfo.Nx);
		complex** dyGam_hat  = declare_2darray_complex(dominfo.Ny,dominfo.Nx);
		fftw_2d_r2c(Gam_hat,gamheav,dominfo.Ny,dominfo.Nx);
		fftw_2d_r2c(dyGam_hat,dygamheav,dominfo.Ny,dominfo.Nx);
		
		ifftw_2d_c2r(newGam,Gam_hat,dominfo.Ny,dominfo.Nx);
		ifftw_2d_c2r(SkewGam,dyGam_hat,dominfo.Ny,dominfo.Nx);
		
		free_2darray(dygamheav,dominfo.Ny,dominfo.Nx);
		free_2darray_complex(Gam_hat,dominfo.Ny,dominfo.Nx);
		free_2darray_complex(dyGam_hat,dominfo.Ny,dominfo.Nx);
	}
		
	//prepare for ff and Skew_ff (Nxy*Nt)
    double ej[2]; ej[0]=1; ej[1]=0;
	double**  UgKK1 	= funOprt_Ug2d(k1,K1,meandepthinf,ej,Nxy,dominfo.Nt);
	
	complex** UgxyKK1 	= ifftw_2d_r2c(UgKK1,Nxy,dominfo.Nt);       
    GamKK1_hat = fftw_2d_c2c(UgxyKK1,Nxy,dominfo.Nt);   
        
	UgKK_hat   = funOprt_Ug2d(k1,K1,meandepthinf,ej,Nxy,dominfo.Nt);
	ffact      = dxy/(2*M_PI);
		
	double** WW1  = declare_2darray(Nxy,dominfo.Nt);
	double** Op_C = declare_2darray(Nxy,dominfo.Nt);
	double K;
	for (int i=0;i<Nxy;i++){
		for (int j=0;j<dominfo.Nt;j++){
			K = fun_sign(K1[j]+k1[i])*sqrt(pow(K1[j],2)+pow(k1[i],2));
			WW1[j][i]  = fun_exact_disp1d_val(K,meandepthinf,dominfo.U);
			Op_C[j][i] = WW1[j][i]/K1[j];
		}		
	}
	
	Insig_hat  = declare_2darray_complex(wave_influx.Nxy,dominfo.Nt);
	fftw_2d_r2c(Insig_hat,wave_influx.eta,wave_influx.Nxy,dominfo.Nt);//no scaling
	
	complex** 	ff_hat 		= declare_2darray_complex(Nxy,dominfo.Nt);
	complex** 	Skewf_hat	= declare_2darray_complex(Nxy,dominfo.Nt);
	complex tempAB;
	for (int i=0;i<Nxy;i++){
		for (int j=0;j<dominfo.Nt;j++) {
			tempAB = fun_complex_division(Insig_hat[j][i],GamKK1_hat[j][i]);
			ff_hat[j][i].Re = ffact*tempAB.Re/(2*M_PI)*abs(K1[j])/KK[j][i]*UgKK_hat[j][i];
			ff_hat[j][i].Im = ffact*tempAB.Im/(2*M_PI)*abs(K1[j])/KK[j][i]*UgKK_hat[j][i];
			
			Skewf_hat[j][i].Re = Op_C[j][i]*ff_hat[j][i].Im/(WW1[j][i]);
			Skewf_hat[j][i].Im = -Op_C[j][i]*ff_hat[j][i].Re/(WW1[j][i]);
		}
	}
	
	double** ff		= declare_2darray(Nxy,dominfo.Nt);
	double** Skewf	= declare_2darray(Nxy,dominfo.Nt);
	double**  Signal;
	double**  SkewSignal;
	
	ifftw_2d_c2r(ff,ff_hat,Nxy,dominfo.Nt);
	ifftw_2d_c2r(Skewf,Skewf_hat,Nxy,dominfo.Nt);
	
	//applying ramp and taper
	for (int i=0;i<wave_influx.Nxy;i++){
		for (int j=0;j<dominfo.Nt;j++){
			ff[j][i]    = ff[j][i]*ramp2d[j][i];
			Skewf[j][i] = Skewf[j][i]*ramp2d[j][i];
		}
	}	
	Signal 		= transpose(ff,Nxy,dominfo.Nt);
	SkewSignal	= transpose(Skewf,Nxy,dominfo.Nt);
	

		/*double* forline=new double[dominfo.Nt];
		for (int i=0;i<dominfo.Nt;i++){
			forline[i]=Skewf_hat[i][0].Re;
		}
		fun_plot_line(dominfo.t,forline,dominfo.Nt);
		delete[] forline;*/
		
	
	//signal and skewsignal interpolation
	Spline_Signal	  = new interp_sig[Nxy];
	Spline_SkewSignal = new interp_sig[Nxy];
	for (int j=0;j<Nxy;j++){
		Spline_Signal[j].acc 	= gsl_interp_accel_alloc ();
		Spline_Signal[j].spline = gsl_spline_alloc (gsl_interp_cspline,dominfo.Nt);
		Spline_Signal[j].init 	= gsl_spline_init(Spline_Signal[j].spline,dominfo.t,Signal[j],dominfo.Nt);
		
		Spline_SkewSignal[j].acc 	= gsl_interp_accel_alloc ();
		Spline_SkewSignal[j].spline = gsl_spline_alloc (gsl_interp_cspline,dominfo.Nt);
		Spline_SkewSignal[j].init 	= gsl_spline_init(Spline_SkewSignal[j].spline,dominfo.t,SkewSignal[j],dominfo.Nt);
	}	

	delete[] kk;
	delete[] k1;
	delete[] K1;
	
	//plotting
	plotting_preproc(Signal,newGam);
	
	free_2darray(KK,Nxy,dominfo.Nt);
	free_2darray(Gam,dominfo.Ny,dominfo.Nx);
	free_2darray(gamheav,dominfo.Ny,dominfo.Nx);
	free_2darray(UgKK1,Nxy,dominfo.Nt);
	free_2darray(UgKK_hat,Nxy,dominfo.Nt);
	free_2darray_complex(UgxyKK1,Nxy,dominfo.Nt);
	free_2darray_complex(GamKK1_hat,Nxy,dominfo.Nt);
	free_2darray(Op_C,Nxy,dominfo.Nt);
	free_2darray(WW1,Nxy,dominfo.Nt); 
	free_2darray_complex(ff_hat,Nxy,dominfo.Nt);
	free_2darray_complex(Skewf_hat,Nxy,dominfo.Nt);	
	
	free_2darray(Signal,dominfo.Nt,Nxy);
	free_2darray(SkewSignal,dominfo.Nt,Nxy);
	free_2darray(ff,Nxy,dominfo.Nt);
	free_2darray(Skewf,Nxy,dominfo.Nt);
} 

void fftw_to_complex(int Nx,int Ny,fftw_complex* A,complex** B)
{
	#pragma omp parallel for
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			B[i][j].Re = A[i*Nx+j][0];
			B[i][j].Im = A[i*Nx+j][1];
		}
	}	
}

void define_coupling_char(double* x,double* y,int Nx,int Ny)
{ 
	int id_h0 		= (int)ceil((ozstart[0]-dominfo.x[0])/dominfo.dx);//id_Xinterv1+nLin;
    int id_h1 		= id_h0+nLtr+nLin;    
    int id_h3 		= (int)ceil((ozend[0]-dominfo.x[0])/dominfo.dx);//id_h0+2*nLtr+nCFD;
    int id_h2 		= id_h3-nLtr-nLin;
    int id_peak1	= id_h0+nLtr;
	int id_peak2	= id_h3-nLtr;
        
    int id_h0y 		= (int)ceil((ozstart[1]-dominfo.y[0])/dominfo.dy);//id_Yinterv1+nLiny;
	int id_h1y 		= id_h0y+nLtry+nLiny;
	int id_h3y 		= (int)ceil((ozend[1]-dominfo.y[0])/dominfo.dy);//id_h0y+2*nLtry+nCFDy;
	int id_h2y 		= id_h3y-nLtry-nLiny;	
    int id_peak1y	= id_h0y+nLtry; 
	int id_peak2y	= id_h3y-nLtry;
	
	//x-direction
	int Nx1,Ny1;		
	Nx1 = id_h3-id_h0+1;
	Ny1 = id_h3y-id_h0y+1;	
	double** cfSA_x = declare_2darray(Nx1,Ny1);
	double *xx = new double[Nx1];
	double *yy = new double[Ny1];
	fbdyvar fbdy_cfd;
	fbdy_cfd.Lleft   = x[id_peak1]-x[id_h0];
	fbdy_cfd.Lright  = x[id_peak1]-x[id_h0];
	fbdy_cfd.Lbottom = y[id_peak1y]-y[id_h0y];
	fbdy_cfd.Ltop 	 = y[id_peak1y]-y[id_h0y];
	
	for (int i=0;i<Nx1;i++){
		xx[i] = dominfo.x[id_h0+i];
	}
	for (int i=0;i<Ny1;i++){
		yy[i] = dominfo.y[id_h0y+i];
	}
	fun_cfSA2D(cfSA_x,xx,yy,Nx1,Ny1,fbdy_cfd);
	
	//y-direction
	int Nx2,Ny2;		
	Nx2 = id_peak2-id_peak1+1;
	Ny2 = id_peak2y-id_peak1y+1;	
	double** cfSA_y = declare_2darray(Nx2,Ny2);
	double *xx2 = new double[Nx2];
	double *yy2 = new double[Ny2];
	fbdy_cfd.Lleft   = x[id_h1]-x[id_peak1];
	fbdy_cfd.Lright  = x[id_peak2]-x[id_h2];
	fbdy_cfd.Lbottom = y[id_h1y]-y[id_peak1y];
	fbdy_cfd.Ltop 	 = y[id_peak2y]-y[id_h2y];
		
	for (int i=0;i<Nx2;i++){
		xx2[i] = dominfo.x[id_peak1+i];
	}
	for (int i=0;i<Ny2;i++){
		yy2[i] = dominfo.y[id_peak1y+i];
	}
	fun_cfSA2D(cfSA_y,xx2,yy2,Nx2,Ny2,fbdy_cfd);
	
	//filling the characteristic function	
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if ((i>=id_h0y)&&(i<=id_h3y)&&(j>=id_h0)&&(j<=id_h3)){
				Chi_couple_H[i][j] = (cfSA_x[i-id_h0y][j-id_h0]);
			}
		}
	}	
	
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if ((i>=id_peak1y)&&(i<=id_peak2y)&&(j>=id_peak1)&&(j<=id_peak2)){
				Chi_couple_H[i][j] = (1-cfSA_y[i-id_peak1y][j-id_peak1]);
			}
		}
	}
	
	/*FILE* gp=popen("gnuplot -persistent","w");
	fun_plot_2d_trarray(gp,y,x,Chi_couple_H,Ny,Nx);
	fflush(gp);
	pclose(gp);*/

}

/*void logfile(char* argv2,char* argv1,int n)
{
	char data_log[128];
	sprintf(data_log, "%s/logfile.txt", argv2);
    FILE* flog   = fopen(data_log,"w");
    
    time_t tnow;   // not a primitive datatype
    time(&tnow);
    fprintf(flog,"HAWASSI 2D simulation has been run on %s\n",ctime(&tnow));
    fprintf(flog,"File input: %s/input.txt \n", argv1);
    fprintf(flog,"MODEL DESCRIPTION:\n"
    "Dynamic Model \t\t: %s\n"
    "Dispersion Model \t: exact\n"
    "Breaking \t\t: %s\n",evolinfo.model,breaking);
    if (strcmp(breaking,"yes")==0){
		fprintf(flog,"\t Initiation \t: %.2f\n"
					 "\t Termination \t: %.2f %.2f\n",parKBC,parTC,parTchar);
	}
	fprintf(flog,"Bottom friction \t: %s\n",friction);
	if (strcmp(friction,"yes")==0){
		fprintf(flog,"The friction coefficient data is from %s/cf.dat\n",argv2);
	}
	fprintf(flog,"Kinematics \t\t: %s\n",kinematic);
	if (strcmp(kinematic,"yes")==0){
		fprintf(flog,"\t Overlay zone \t: [%.2f,%.2f]\n",xinterv1,xinterv2);
		fprintf(flog,"\t Non-uniform vertical grid : [%.2f,%.2f]\n",zoutIP[0],zoutIP[NzoutIP-1]);
	}
	
	if (strcmp(coupling,"yes")==0){
		
	}
	
	fprintf(flog,"\n");
	fprintf(flog,"INFLUX DESCRIPTION :\n"
	"Signal type \t\t: %s\n",wavename);
	if ((strcmp(wavename,"jonswap")==0) || (strcmp(wavename,"harmonic")==0)){		
		fprintf(flog,"Significant wave Height (Hs) \t: %.2f\n"
			"Peak period (Tp) \t\t: %.2f\n",Hs_in,Tp);
	}
	if (strcmp(wavename,"zero")!=0){
		fprintf(flog,"Derived info:\n"
			"\t Peak frequency (nu) \t: %.2f[rad/s]\n"
			"\t Peak wave-number (kp) \t: %.2f\n"
			"\t Peak wave-length \t: %.2f[m]\n"
			"\t Steepness (kp*(Hs./2)) : %.2f\n",nu_peak,kp,lambda_p,kp*(Hs_in/2));
	}    
	
	fprintf(flog,"\n");
	fprintf(flog,"INITIAL WAVE CONDITIONS : \n"
		"Initial condition : %s\n",initial);
		
	fprintf(flog,"\n");
	fprintf(flog,"NUMERICAL SETTINGS : \n"
		"Spatial interval: (%.2f,%.2f)[m]\n"
		"Damping zone \t: %.2f[m]\n"
		"Number of Nodes : %d\n"
		"Grid size (dx) \t: %.2f[m]\n"
        "Cutfrac k \t: %d\n"
		"Time interval \t: (%.2f,%.2f)[s]\n"
        "Time step (dt) \t: %.2f[s]\n",x[0],x[Nx-1],damp,Nx,dx,cutfracwn,tim[0],tim[n-1],tim[1]-tim[0]);
    fprintf(flog,"Bathymetry \t: %s\n",bath);
    if (strcmp(bath,"flat")==0){
		fprintf(flog,"\t Depth \t: %.2f[m]\n",depthflat);
	}
	else {
		fprintf(flog,"The bathymetry data is from %s/bath.dat\n",argv2);
	}
	
	fprintf(flog,"\n");
	
	if (strcmp(wavename,"zero")!=0){
		fprintf(flog,"INFLUXING : \n"
			"Generation method : %s\n"
			"Direction \t: %s\n"
			"Influx position : %.2f[m]\n",influxing,propagation,Xinflux);
		if (strcmp(dynmodel,"HsF1")!=0){
			fprintf(flog,"Nonlinear adjustment: %g*lambda_peak\n",adjcoef);
		}
	}

	fprintf(flog,"\n");
	fprintf(flog,"OUTPUT SETTING: \n"
		//"ODE solver \t\t: ode45\n"
		"No of Partition : %d\n",npartition);
	
	fprintf(flog,"\n");
	fprintf(flog,"After simulation:\n"
		//"CPUtime for ODEs and write output: %.2f sec\n"
		"Relative time consuming: %g\n"
		"Data saved in the folder: %s",rel,argv2);  
	fclose(flog);
}*/
