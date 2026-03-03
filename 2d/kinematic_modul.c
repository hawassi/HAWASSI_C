void save_elev(char var[20],char str[20],double N1,double N2x,double N2y,double N3,double N4,double N5,int ip)
{
	char fdata[256];
	strcpy(fdata,arg2);
	strcat(fdata,"/");
	strcat(fdata,"Hawassi_");
	strcat(fdata,var);
	strcat(fdata,"_");
	strcat(fdata,str);
	dateta = fopen(fdata,"wb");
	
	fwrite(&N4,sizeof(double),1,dateta);
	fwrite(&N1,sizeof(double),1,dateta);
	fwrite(&N2x,sizeof(double),1,dateta);
	fwrite(&N2y,sizeof(double),1,dateta);
	
	if (ip<=npartition){
		fwrite(toutIP[ip-1],sizeof(double),NtoutIP,dateta);
	}
	else {
		fwrite(trest,sizeof(double),Ntrest,dateta);
	}

	fwrite(xoutIP,sizeof(double),NxoutIP,dateta);
	fwrite(youtIP,sizeof(double),NyoutIP,dateta);
	fwrite(zoutIP,sizeof(double),NzoutIP,dateta);
}

void save_vel(char var[20],char str[20],double N1,double N2x,double N2y,double N3,double N4,double N5,int ip)
{
	char fdata[256];
	strcpy(fdata,arg2);
	strcat(fdata,"/");
	strcat(fdata,"Hawassi_");
	strcat(fdata,var);
	strcat(fdata,"_");
	strcat(fdata,str);
	datvel = fopen(fdata,"wb");
	
	fwrite(&N5,sizeof(double),1,datvel);
	fwrite(&N1,sizeof(double),1,datvel);
	fwrite(&N2x,sizeof(double),1,datvel);
	fwrite(&N2y,sizeof(double),1,datvel);
	fwrite(&N3,sizeof(double),1,datvel);
	
	if (ip<=npartition){
		fwrite(toutIP[ip-1],sizeof(double),NtoutIP,datvel);
	}
	else {		
		fwrite(trest,sizeof(double),Ntrest,datvel);
	}
	fwrite(xoutIP,sizeof(double),NxoutIP,datvel);
	fwrite(youtIP,sizeof(double),NyoutIP,datvel);
	fwrite(zoutIP,sizeof(double),NzoutIP,datvel);
}

void save_P(char var[20],char str[20],double N1,double N2x,double N2y,double N3,double N4,double N5,int ip)
{
	char fdata[256];
	strcpy(fdata,arg2);
	strcat(fdata,"/");
	strcat(fdata,"Hawassi_");
	strcat(fdata,var);
	strcat(fdata,"_");
	strcat(fdata,str);
	datP = fopen(fdata,"wb");
	
	fwrite(&N5,sizeof(double),1,datP);
	fwrite(&N1,sizeof(double),1,datP);
	fwrite(&N2x,sizeof(double),1,datP);
	fwrite(&N2y,sizeof(double),1,datP);
	fwrite(&N3,sizeof(double),1,datP);
	
	if (ip<=npartition){
		for (int j=0;j<NtoutIP;j++){
			fwrite(&toutIP[ip-1][j],sizeof(double),1,datP);
		}
	}
	else {		
		fwrite(trest,sizeof(double),Ntrest,datP);
	}
	fwrite(xoutIP,sizeof(double),NxoutIP,datP);
	fwrite(youtIP,sizeof(double),NyoutIP,datP);
	fwrite(zoutIP,sizeof(double),NzoutIP,datP);
}

void save_acc(char var[20],char str[20],double N1,double N2x,double N2y,double N3,double N4,double N5,int ip)
{
	char fdata[256];
	strcpy(fdata,arg2);
	strcat(fdata,"/");
	strcat(fdata,"Hawassi_");
	strcat(fdata,var);
	strcat(fdata,"_");
	strcat(fdata,str);
	datacc = fopen(fdata,"wb");
	
	fwrite(&N5,sizeof(double),1,datacc);
	fwrite(&N1,sizeof(double),1,datacc);
	fwrite(&N2x,sizeof(double),1,datacc);
	fwrite(&N2y,sizeof(double),1,datacc);
	fwrite(&N3,sizeof(double),1,datacc);
	
	if (ip<=npartition){
		fwrite(toutIP[ip-1],sizeof(double),NtoutIP,datacc);
	}
	else {		
		fwrite(trest,sizeof(double),Ntrest,datacc);
	}
	fwrite(xoutIP,sizeof(double),NxoutIP,datacc);
	fwrite(youtIP,sizeof(double),NyoutIP,datacc);
	fwrite(zoutIP,sizeof(double),NzoutIP,datacc);
}

void save_phi(char var[20],char str[20])
{
	char fdata[256];
	strcpy(fdata,arg2);
	strcat(fdata,"/");
	strcat(fdata,"Hawassi_");
	strcat(fdata,var);
	strcat(fdata,"_");
	strcat(fdata,str);
	datPhi = fopen(fdata,"wb");
}


void fun_spmatrix(double** Af,int nx, int ny,double dx,double dy)
{ 
	//define memory for array as global variable
	int L;
	int* indx= new int[nx*ny];
	int* indy= new int[nx*ny];
	int* ind= new int[nx*ny];
	int** rind = declare_2darray_int(nx*ny,2); 
	int** cind = declare_2darray_int(nx*ny,2);  
	double** dfdx = declare_2darray(nx*ny,2); 
	double*  rhs = new double[2*nx*ny];
	int m = (nx-2)*ny;
		
	//do the leading edge in x, forward difference
	L=0;
	indx[0] = 0;
	#pragma parallel for
	for (int j=0;j<ny;j++){
		indy[j] 	= j;
		ind[j]  	= indy[j];
		rind[0][j]	= L+indy[j];
		rind[1][j]	= L+indy[j];
		cind[0][j]	= ind[j];
		cind[1][j]	= ind[j]+ny;
		dfdx[0][j]	= -1/dx;
		dfdx[1][j]	= 1/dx;
		Af[L+j][0]  = rind[0][j];
		Af[L+j][1]  = rind[1][j];
		Af[L+j][2]  = cind[0][j];
		Af[L+j][3]  = cind[1][j];
		Af[L+j][4]  = dfdx[0][j];
		Af[L+j][5]  = dfdx[1][j];
	}
	
	//interior partials in x, central difference
	L=L+ny;
	if (nx>2){
		int ix = 0;
		#pragma parallel for
		for (int j=0;j<m;j++){
			if ((j % ny)==0){
				ix++;
			}
			indx[j] = ix;
			indy[j] = (j % ny);
			ind[j]  = indy[j] + (indx[j])*ny;
			rind[0][j]	= L+j;
			rind[1][j]	= L+j;
			cind[0][j]	= ind[j]-ny;
			cind[1][j]	= ind[j]+ny;
			dfdx[0][j]	= -1/(2*dx);
			dfdx[1][j]	= 1/(2*dx);
			Af[L+j][0]  = rind[0][j];
			Af[L+j][1]  = rind[1][j];
			Af[L+j][2]  = cind[0][j];
			Af[L+j][3]  = cind[1][j];
			Af[L+j][4]  = dfdx[0][j];
			Af[L+j][5]  = dfdx[1][j];
		}		
	}
	
	//do the trailing edge in x, backward difference	
	L = L+m;
	indx[0] = nx-1;
	#pragma parallel for
	for (int j=0;j<ny;j++){
		indy[j] 	= j;
		ind[j]  	= indy[j] + (indx[0])*ny;
		rind[0][j]	= L+indy[j];
		rind[1][j]	= L+indy[j];
		cind[0][j]	= ind[j]-ny;
		cind[1][j]	= ind[j];		
		dfdx[0][j]	= -1/dx;
		dfdx[1][j]	= 1/dx;
		
		Af[L+j][0]  = rind[0][j];
		Af[L+j][1]  = rind[1][j];
		Af[L+j][2]  = cind[0][j];
		Af[L+j][3]  = cind[1][j];
		Af[L+j][4]  = dfdx[0][j];
		Af[L+j][5]  = dfdx[1][j];
	}
	
	//do the leading edge in y, forward difference
	L = L+ny;
	indy[0] = 0;
	#pragma parallel for
	for (int j=0;j<nx;j++){
		indx[j] 	= j;
		ind[j]  	= indy[0] + (indx[j])*ny;
		rind[0][j]	= L+indx[j];
		rind[1][j]	= L+indx[j];
		cind[0][j]	= ind[j];
		cind[1][j]	= ind[j]+1;		
		dfdx[0][j]	= -1/dy;
		dfdx[1][j]	= 1/dy;
		
		Af[L+j][0]  = rind[0][j];
		Af[L+j][1]  = rind[1][j];
		Af[L+j][2]  = cind[0][j];
		Af[L+j][3]  = cind[1][j];
		Af[L+j][4]  = dfdx[0][j];
		Af[L+j][5]  = dfdx[1][j];
	}

	//interior partials in y, use a central difference
	L = L+nx;	
	if (ny>2){
		m = nx*(ny-2);
		int ix = -1;
		#pragma parallel for
		for (int j=0;j<m;j++){
			if ((j % (ny-2))==0){
				ix++;
			}
			indy[j] = (j % (ny-2))+1;
			indx[j] = ix;
			ind[j]  = indy[j] + (indx[j])*ny;	
			rind[0][j]	= L+j;
			rind[1][j]	= L+j;
			cind[0][j]	= ind[j]-1;
			cind[1][j]	= ind[j]+1;
			dfdx[0][j]	= -1/(2*dy);
			dfdx[1][j]	= 1/(2*dy);
			Af[L+j][0]  = rind[0][j];
			Af[L+j][1]  = rind[1][j];
			Af[L+j][2]  = cind[0][j];
			Af[L+j][3]  = cind[1][j];
			Af[L+j][4]  = dfdx[0][j];
			Af[L+j][5]  = dfdx[1][j];
		}		
	}

	//do the trailing edge in y, backward diffeence
	L = L+m;
	indy[0] = ny-1;
	#pragma parallel for
	for (int j=0;j<nx;j++){
		indx[j] 	= j;
		ind[j]  	= indy[0] + (indx[j])*ny;
		rind[0][j]	= L+indx[j];
		rind[1][j]	= L+indx[j];
		cind[0][j]	= ind[j]-1;
		cind[1][j]	= ind[j];		
		dfdx[0][j]	= -1/dy;
		dfdx[1][j]	= 1/dy;
		
		Af[L+j][0]  = rind[0][j];
		Af[L+j][1]  = rind[1][j];
		Af[L+j][2]  = cind[0][j];
		Af[L+j][3]  = cind[1][j];
		Af[L+j][4]  = dfdx[0][j];
		Af[L+j][5]  = dfdx[1][j];
	}
	
	//filling sparse matrix shifting one column to the left
	Eigen::initParallel();
    Eigen::setNbThreads(ncpu);
	A.resize(2*nx*ny,nx*ny-1);
	A.reserve(3*2*nx*ny);
	#pragma parallel for
	for (int j=0;j<2*nx*ny;j++){	
		if (int (Af[j][2]) > 0)	{
			A.insert(int (Af[j][0]),int (Af[j][2]-1)) = Af[j][4];	
		}
		if (int (Af[j][3]) > 0){
			A.insert(int (Af[j][1]),int (Af[j][3]-1)) = Af[j][5];
		}
	}
	//A.makeCompressed();
	//ldlt.compute(A.transpose()*A);
	
	lscg.compute(A);
	
	delete[] indx;
	delete[] indy;
	delete[] ind;
	free_2darray_int(rind,nx*ny,2); 
	free_2darray_int(cind,nx*ny,2); 
	free_2darray(dfdx,nx*ny,2); 
	delete[] rhs;
}

void read_kinematic(double** depth, int Nx, int Ny,int Nt)
{	
	depth_IP=declare_2darray(Nx,Ny);
	for (int i=0;i<Ny;i++) {
		for (int j=0;j<Nx;j++) {
			depth_IP[i][j] = depth[i][j]*dominfo.bathy.Char[i][j];
		}
	}
	
	/*FILE* gph=popen("gnuplot -persistent","w");
	fun_plot_2d_trarray(gph,dominfo.y,dominfo.x,depth_IP,Ny,Nx);
	fflush(gph);
	pclose(gph);getchar();*/
			
	Af=declare_2darray(6,2*Nx*Ny);	
	fun_spmatrix(Af,Nx,Ny,dominfo.dx,dominfo.dy);
	
	ip_Xinterv1	= (int)floor((xinterv1-dominfo.x[0])/dominfo.dx);
	ip_Xinterv2	= (int)floor((xinterv2-dominfo.x[0])/dominfo.dx);
	ip_Yinterv1	= (int)floor((yinterv1-dominfo.y[0])/dominfo.dy);
	ip_Yinterv2	= (int)floor((yinterv2-dominfo.y[0])/dominfo.dy);
	
	//define max-min depth for interpolation
	maxD 	= fun_max2darray(depth,Nx,Ny);
	minD 	= fun_min2darray(depth,Nx,Ny);
	
	//grid-x for kinematic_modul
	NxoutIP = ip_Xinterv2-ip_Xinterv1+1;
	xoutIP 	= new double[NxoutIP];		
	for (int i=0;i<NxoutIP;i++){
		xoutIP[i] = dominfo.x[ip_Xinterv1+i];
	}
	
	//grid-y for kinematic_modul
	NyoutIP = ip_Yinterv2-ip_Yinterv1+1;
	youtIP 	= new double[NyoutIP];		
	for (int i=0;i<NyoutIP;i++){
		youtIP[i] = dominfo.y[ip_Yinterv1+i];
	}
	NxyoutIP= NxoutIP*NyoutIP;
	
	//grid-time
	NtoutIP = dominfo.Nt/npartition;
	Ntrest  = (Nt-npartition*NtoutIP);
	toutIP  = new double*[npartition];
	trest   = new double[Ntrest];
	for (int j=0;j<npartition;j++){
		toutIP[j] = new double[NtoutIP];
		for (int tt=0;tt<NtoutIP;tt++){
			toutIP[j][tt] = dominfo.t[tt+(j*NtoutIP)];
		}
	}
	for (int j=0;j<Ntrest;j++){
		trest[j] = dominfo.t[npartition*NtoutIP+j];
	}
	
	//estimate vertical grid
	double zsurf1,zsurf2;	
	if (strcmp(min_z,"default")==0){
		zsurf1 = -2.5*seainfo.Hs;
	}
	else {
		zsurf1 = atof(min_z);
	}
		
	if (strcmp(max_z,"default")==0){
		zsurf2 = 2.5*seainfo.Hs;
	}
	else {
		zsurf2 = atof(max_z);
	}
		
	//define zout, z at surface layer must be lower than maximum eta
	maxEta = zsurf2; //user
	minEta = zsurf1;
	
	//compute zoutIP and NzoutIP
	if (strcmp(kinematic_type,"boundary")==0){//using GQ points, z changes every time	
		int Nspg=20;	
		GQ = new double[Nspg];
		getGQPoints(GQ);
		xscale = maxD;
		tscale = sqrt(xscale/grav);
		potscale= xscale*xscale/tscale;
		
		//horizontal grid
		Nxrec = ((xoutIP[NxoutIP-1]-xoutIP[0])/dx_rec)+1;
		Nyrec = ((youtIP[NyoutIP-1]-youtIP[0])/dy_rec)+1;
		xrec  = new double[Nxrec];
		yrec  = new double[Nyrec];
		for (int i=0;i<Nxrec;i++) {
			xrec[i]=xoutIP[0]+i*dx_rec;
		}
		for (int i=0;i<Nyrec;i++) {
			yrec[i]=youtIP[0]+i*dy_rec;
		}		
		
		//vertical grid	
		double minH = minEta+minD;
		if (minH<0){
			minH=minD;
		}
		double dz1percent = 0.01*(minH);
		NzoutIP    = int ((zsurf2+maxD)/dz1percent)+1;
		dz1percent = (zsurf2+maxD)/(NzoutIP-1);
		zoutIP 	   = new double[NzoutIP];	 
		for (int i=0;i<NzoutIP;i++) {
			zoutIP[i] = -maxD+i*dz1percent;
		}
		
		//depth at record
		depthrec = declare_2darray(Nxrec,Nyrec);
		fun_interp2(Nx,Ny,dominfo.x,dominfo.y,depth,Nxrec,Nyrec,xrec,yrec,depthrec);//depth at rec domain	
		
		//writing header record
		sprintf(filename,"%sHAWASSI2D.rec",arg2);
		UDW2p_setUp_hw2D(Nxrec,Nyrec,Nspg,Nt,dx_rec,dy_rec,dominfo.dt,xscale,tscale,depthrec);	
		
	}
	else {		
		double* zsurf  = (double*) malloc (sizeof(double)*nsurf);
		double  dzsurf = (maxEta-minEta)/(nsurf-1);
		double 	kpeak  = seainfo.kp;
			
		for (int j=0;j<nsurf;j++){
			zsurf[j] = minEta+j*dzsurf;
		}
			
		//compute z at deeper layer using Riemann upper sum 
		double  Flux,tol;
		
		Flux = 1.0/(kpeak*sinh(kpeak*maxD))*(cosh(kpeak*(maxD+minEta))-1);
		tol  = Flux/(ndeep);	
		
		double* zdeep = (double*) malloc (sizeof(double)*ndeep);
		zdeep[0]  = -1*(minEta-dzsurf);
			
		for (int i=1;i<ndeep;i++) {			
			zdeep[i]  = zdeep[i-1]+(tol/(sinh(kpeak*(maxD-zdeep[i-1]))/sinh(kpeak*maxD)));
		}
			
		//combining surface layer and deeper layer
		NzoutIP = nsurf+ndeep;
		zoutIP 	= new double[NzoutIP];	 
		for (int i=0;i<nsurf;i++) {
			zoutIP[i] = zsurf[i];
		}
		for (int i=nsurf;i<NzoutIP;i++) {
			zoutIP[i] = -1*zdeep[i-nsurf];
		}
		sorting(zoutIP,NzoutIP);		
		free(zsurf);free(zdeep);
	}
	
	//define variable string
	sprintf(var_elev,"Elev");
	sprintf(var_vel,"Vel");
	sprintf(var_P,"P");
	sprintf(var_acc,"Acc");
	sprintf(var_phi,"phi");
	Velx = declare_2darray(NzoutIP,NxyoutIP);
	Vely = declare_2darray(NzoutIP,NxyoutIP);
	Velz = declare_2darray(NzoutIP,NxyoutIP);
	set_matrix_val(Velx,NzoutIP,NxyoutIP,0.0);
	set_matrix_val(Vely,NzoutIP,NxyoutIP,0.0);
	set_matrix_val(Velz,NzoutIP,NxyoutIP,0.0);	
}

void fun_lscgsolver(double **f,double** fx,double** fy,double dx,double dy,int nx, int ny,double f11)
{   
	//define memory for array as global variable
	int L;
	int indx[nx*ny];
	int indy[nx*ny];
	int ind[nx*ny];
	int m = (nx-2)*ny;
	double  rhs[2*nx*ny];
	
	//do the leading edge in x, forward difference
	L=0;
	#pragma parallel for
	for (int j=0;j<ny;j++){
		rhs[L+j] 	= fx[j][0];
	}
	
	//interior partials in x, central difference
	L=L+ny;
	if (nx>2){
		int ix = 0;
		#pragma parallel for
		for (int j=0;j<m;j++){
			if ((j % ny)==0){
				ix++;
			}
			indx[j] = ix;
			indy[j] = (j % ny);
			rhs[L+j]= fx[(indy[j])][(indx[j])];
		}		
	}
	
	//do the trailing edge in x, backward difference	
	L = L+m;
	#pragma parallel for
	for (int j=0;j<ny;j++){
		rhs[L+j] 	= fx[j][nx-1];
	}
	
	//do the leading edge in y, forward difference
	L = L+ny;
	#pragma parallel for
	for (int j=0;j<nx;j++){
		rhs[L+j] 	= fy[0][j];
	}

	//interior partials in y, use a central difference
	L = L+nx;	
	if (ny>2){
		m = nx*(ny-2);
		int ix = -1;
		#pragma parallel for
		for (int j=0;j<m;j++){
			if ((j % (ny-2))==0){
				ix++;
			}
			indy[j] = (j % (ny-2))+1;
			indx[j] = ix;
			rhs[L+j] 	= fy[(indy[j])][(indx[j])];
		}		
	}

	//do the trailing edge in y, backward diffeence
	L = L+m;
	#pragma parallel for
	for (int j=0;j<nx;j++){
		rhs[L+j] 	= fy[ny-1][j];
	}
	
	//allocate memory for Ax=b
	Eigen::VectorXd b(2*nx*ny);
	Eigen::VectorXd x(nx*ny-1);
	
	//filling sparse matrix shifting one column to the left
	#pragma parallel for
	for (int j=0;j<2*nx*ny;j++){		
		b(j) = rhs[j];
	}
	
	x = (lscg).solve(b);
	
	double xnew[nx*ny];
    xnew[0]=f11;
    for (int i=0;i<nx*ny-1;i++) {
		xnew[i+1]=x(i);
	}

	for (int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			f[i][j]=xnew[j*ny+i];
		}
	}
}

void sinhcosh2IP_2D(int sizD,double* D,double Dmin,double Dplus,double** kk,double nupeak,double* gam_min,double* gam_plus)
{
	double*  kappanu_D	=new double[sizD];
	double*  cosh_D		=new double[sizD];
	double*  cosh_min	=new double[sizD];
	double*  cosh_plus	=new double[sizD];
	double*  sinh_D		=new double[sizD];
	double*  sinh_min	=new double[sizD];
	double*  sinh_plus	=new double[sizD];	
	
	for (int j=0; j < sizD; j++) {
		if (D[j] < 0.0001) {
			kappanu_D[j] = kk[1][1];
		}		
		else {
			kappanu_D[j] = invers_omega(nupeak,D[j]);
		}
		//printf("%f \n",D[j]);
		cosh_min[j]	 = cosh(kappanu_D[j]*Dmin);
		cosh_plus[j] = cosh(kappanu_D[j]*Dplus);
		cosh_D[j]    = cosh(kappanu_D[j]*D[j]);
		sinh_min[j]  = sinh(kappanu_D[j]*Dmin);
		sinh_plus[j] = sinh(kappanu_D[j]*Dplus);
		sinh_D[j]    = sinh(kappanu_D[j]*D[j]);
	}
	
	double det[sizD],detA[sizD],detB[sizD],gminA[sizD],gminB[sizD],gplusA[sizD],gplusB[sizD];
	
	inner_product_array(detA,cosh_min,sinh_plus,sizD);
	inner_product_array(detB,cosh_plus,sinh_min,sizD);
	inner_product_array(gminA,sinh_plus,cosh_D,sizD);
	inner_product_array(gminB,cosh_plus,sinh_D,sizD);
	inner_product_array(gplusA,sinh_min,cosh_D,sizD);
	inner_product_array(gplusB,cosh_min,sinh_D,sizD);
	
	for (int i=0; i<sizD;i++) {
		det[i] 		   = detA[i]-detB[i];
		gam_min[i]     = (gminA[i]-gminB[i])/det[i];				
		gam_plus[i]    = (-gplusA[i]+gplusB[i])/det[i];
	}
	
	delete[] kappanu_D;
	delete[] cosh_D;
	delete[] cosh_min;
	delete[] cosh_plus;
	delete[] sinh_D;
	delete[] sinh_min;
	delete[] sinh_plus;	
}

void kinematic_modul(int it,double t, double** eta, double** u, double** v, double** dteta, double** depth,int Nx,int Ny,double** kk,double nupeak)
{			
	//compute phi using integral gradient
	double**  phi	=declare_2darray(Nx,Ny);
	double**  phitemp=declare_2darray(Nx,Ny);
	double**  dtphi	=declare_2darray(Nx,Ny);
	
	fun_lscgsolver(phitemp,u,v,dx,dy,Nx,Ny,0.0);
	
	for (int i=0;i<Ny;i++) {
		for (int j=0;j<Nx;j++) {
			phitemp[i][j]  = phitemp[i][j]*ChiAdjPhi[i][j];
			phi[i][j]      = phitemp[i][j]*cfSA[i][j];
			dtphi[i][j]    = (phi[i][j]-phi_prev[i][j])/dominfo.dt;
			phi_prev[i][j] = phi[i][j];
		}
	}

	complex** phihat=declare_2darray_complex(Nx,Ny);
	complex** dtphihat=declare_2darray_complex(Nx,Ny);
	fftw_2d_r2c(phihat,phi,Nx,Ny);
	fftw_2d_r2c(dtphihat,dtphi,Nx,Ny);

	/*if (it==(dominfo.Nt/10)){
		FILE* gph=popen("gnuplot -persistent","w");
		fun_plot_2d_trarray(gph,dominfo.y,dominfo.x,dtphi,Ny,Nx);
		fflush(gph);
		pclose(gph);
	}*/	
	
	//local header
	int 	 Nxy = Nx*Ny;
	int 	 indx1=ip_Xinterv1,indx2=ip_Xinterv2;
	int 	 indy1=ip_Yinterv1,indy2=ip_Yinterv2;
	double 	 Hplus,Hmin,dhtot,DZmin,DZplus;
	
	// set constant and variables 
	int 	 const NH = 100;
	double   Dshallow = 1.0/20*2*M_PI/seainfo.kp;
		
	typedef struct{
	double gam_Hmin0;
	double Htot;
	}Htot_gam_Hmin0;

	typedef struct{
	double gam_Hplus0;
	double Htot;
	}Htot_gam_Hplus0;

	typedef struct{
	double gam_DZmin0;
	double DZtot;
	}DZtot_gam_DZmin0;

	typedef struct{
	double gam_DZplus0;
	double DZtot;
	}DZtot_gam_DZplus0;
	
    // initialize and read input data at t		
	double*  gam_Hmin0  = new double[NH+1];
	double*  gam_Hplus0 = new double[NH+1];
	double*  gam_DZmin0 = new double[NH+1];
	double*  gam_DZplus0= new double[NH+1];	

	//define H
	double Htot[NH+1];
	double maxEtaH=fun_max2darray(eta,Nx,Ny);
	double mEta=fmax(maxEta,maxEtaH);
	Hplus 	= mEta + maxD;
	if ((runup_id==1)){//this is for runup (!!!)
		minD    = H_minShore;
		Hmin 	= H_minShore;	
	}
	else {		
		Hmin 	= fmax(minEta + minD,0);
	}	
	dhtot 	= (Hplus-Hmin)/NH;
	
	for (int i=0;i<(NH+1);i++) {
		Htot[i] = Hmin + i*dhtot;
	}	
	sinhcosh2IP_2D(NH+1,Htot,Hmin,Hplus,kk,nupeak,gam_Hmin0,gam_Hplus0);

	gsl_interp_accel *accHm0 = gsl_interp_accel_alloc ();
	gsl_spline *splineHm0    = gsl_spline_alloc (gsl_interp_cspline, NH+1);
    gsl_spline_init (splineHm0,Htot,gam_Hmin0,NH+1);
    
    gsl_interp_accel *accHp0  = gsl_interp_accel_alloc ();
    gsl_spline *splineHp0 	  = gsl_spline_alloc (gsl_interp_cspline,NH+1);
    gsl_spline_init (splineHp0,Htot,gam_Hplus0,NH+1); 
   
    //define for DZ
    double 	DZtot[NH+1];
    double 	cosh_kD[Nxy],sinh_kD[Nxy];
	double  gam_DZmin,gam_DZplus;   
	DZmin 	= fmax(0,minD);
	DZplus 	= maxD;
	gsl_interp_accel *accDZm0  = gsl_interp_accel_alloc ();
	gsl_spline *splineDZm0 	   = gsl_spline_alloc (gsl_interp_cspline,NH+1);
	gsl_interp_accel *accDZp0  = gsl_interp_accel_alloc ();
	gsl_spline *splineDZp0 	   = gsl_spline_alloc (gsl_interp_cspline,NH+1);						

	if (minD!=maxD) {
		double Dh;
		for (int i=0;i<(NH+1);i++) {
			DZtot[i] = DZmin+i*(DZplus-DZmin)/NH;
		}		
		sinhcosh2IP_2D(NH+1,DZtot,DZmin,DZplus,kk,nupeak,gam_DZmin0,gam_DZplus0);		
		gsl_spline_init (splineDZm0,DZtot,gam_DZmin0,NH+1);		
		gsl_spline_init (splineDZp0,DZtot,gam_DZplus0,NH+1);

		for (int i=0;i<Ny;i++) {
			for (int j=0;j<Nx;j++) {
				Dh=fmax(depth_IP[i][j],DZmin);
				if (depth_IP[i][j]>0){
					gam_DZmin 		= gsl_spline_eval(splineDZm0,Dh,accDZm0);
					gam_DZplus		= gsl_spline_eval(splineDZp0,Dh,accDZp0);			
					cosh_kD[i*Nx+j] = gam_DZmin*cosh(kk[i][j]*DZmin)+gam_DZplus*cosh(kk[i][j]*DZplus);
					sinh_kD[i*Nx+j] = gam_DZmin*sinh(kk[i][j]*DZmin)+gam_DZplus*sinh(kk[i][j]*DZplus);
				}
				else{
					cosh_kD[i*Nx+j] = 0.0;
					sinh_kD[i*Nx+j] = 0.0;
				}
			}
		}
	}
	
	//allocate memory for output variable in (xy,z) NxyoutIP=Nxout*Nyout
	double** (Phi)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dtPhi)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dxPhi)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dyPhi)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dzPhi)=declare_2darray(NzoutIP,NxyoutIP);
	double** (Ptot)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dt_dxPHI)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dt_dyPHI)=declare_2darray(NzoutIP,NxyoutIP);
	double** (dt_dzPHI)=declare_2darray(NzoutIP,NxyoutIP);           
	
	//calculate the interior properties at time i
	complex phihat_i[Nxy];
	double**  H=declare_2darray(Nx,Ny);
	double** gradx_H=declare_2darray(Nx,Ny);
	double** grady_H=declare_2darray(Nx,Ny);
	double** gradx_D=declare_2darray(Nx,Ny);
	double** grady_D=declare_2darray(Nx,Ny);
	double** gam_Hmin=declare_2darray(Nx,Ny);
	double** gam_Hplus=declare_2darray(Nx,Ny);
	double 	cosh_kH[Nxy],tanh_kH[Nxy];
	double  minEta_i;			
	for (int i=0;i<Ny;i++) {
		for (int j=0;j<Nx;j++) {
			H[i][j]			= fmax(eta[i][j]+depth_IP[i][j],Hmin);
			gam_Hmin[i][j]  = gsl_spline_eval(splineHm0,H[i][j],accHm0);
			gam_Hplus[i][j] = gsl_spline_eval(splineHp0,H[i][j],accHp0);
			cosh_kH[i*Nx+j]	= gam_Hmin[i][j]*cosh(kk[i][j]*Hmin)+gam_Hplus[i][j]*cosh(kk[i][j]*Hplus);
			tanh_kH[i*Nx+j]	= gam_Hmin[i][j]*tanh(kk[i][j]*Hmin)+gam_Hplus[i][j]*tanh(kk[i][j]*Hplus);
			phihat_i[i*Nx+j].Re = phihat[i][j].Re;
			phihat_i[i*Nx+j].Im = phihat[i][j].Im;
		}
	}	
	minEta_i 	= fun_min2darray(eta,Nx,Ny);
	gradient2D(H,dominfo.dx,dominfo.dy,Nx,Ny,gradx_H,grady_H);
	gradient2D(depth_IP,dominfo.dx,dominfo.dy,Nx,Ny,gradx_D,grady_D);
	
	//#pragma omp parallel
	{
		//define new variable for computing the interior properties at each zz	
		double  ZZ,ratio,SINH;
		double*  cosh_kDZ = new double[Nxy];
		double*  sinh_kDZ = new double[Nxy];	
		complex**	factor_nonifft    = declare_2darray_complex(Nx,Ny);
		complex**   tempdtPhi_nonifft = declare_2darray_complex(Nx,Ny);
		complex**   tempdxPhi_nonifft = declare_2darray_complex(Nx,Ny);
		complex** 	tempdyPhi_nonifft = declare_2darray_complex(Nx,Ny);
		complex**   tempdzPhi_nonifft = declare_2darray_complex(Nx,Ny);
		complex**   tempPhi_nonifft = declare_2darray_complex(Nx,Ny);

		double** factor = declare_2darray(Nx,Ny);
		double** tempdtPhi_ifft = declare_2darray(Nx,Ny);
		double** tempdxPhi_ifft = declare_2darray(Nx,Ny);
		double** tempdyPhi_ifft = declare_2darray(Nx,Ny);
		double** tempPhi_ifft = declare_2darray(Nx,Ny);

		double** tempdzPhi = declare_2darray(Nx,Ny);
		double** tempdtPhi = declare_2darray(Nx,Ny);
		double** tempdxPhi = declare_2darray(Nx,Ny);
		double** tempdyPhi = declare_2darray(Nx,Ny);
		double** tempPhi = declare_2darray(Nx,Ny);
		
		//#pragma omp for
		for (int zz=0;zz<NzoutIP;zz++) {
			if (minD==maxD) {//flat bottom
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						cosh_kDZ[i*Nx+j] 	= cosh(kk[i][j]*(zoutIP[zz]+depth[0][0]));
						sinh_kDZ[i*Nx+j] 	= sinh(kk[i][j]*(zoutIP[zz]+depth[0][0]));
					}
				}
			}
			else {//nonflat bottom
				for (int i=0;i<Ny;i++){
					for (int j=0;j<Nx;j++){
						cosh_kDZ[i*Nx+j]=sinh_kD[i*Nx+j]*sinh(kk[i][j]*zoutIP[zz])+cosh_kD[i*Nx+j]*cosh(kk[i][j]*zoutIP[zz]);
						sinh_kDZ[i*Nx+j]=sinh_kD[i*Nx+j]*cosh(kk[i][j]*zoutIP[zz])+cosh_kD[i*Nx+j]*sinh(kk[i][j]*zoutIP[zz]);
						if (cosh_kDZ[i*Nx+j] != cosh_kDZ[i*Nx+j]){cosh_kDZ[i*Nx+j]=0;}
						if (sinh_kDZ[i*Nx+j] != sinh_kDZ[i*Nx+j]){sinh_kDZ[i*Nx+j]=0;}					
					} 
				}	
			}

			for (int i=0;i<Ny;i++){			
				for (int j=0;j<Nx;j++){		
					ratio = cosh_kDZ[i*Nx+j]/cosh_kH[i*Nx+j];
					SINH  = sinh_kDZ[i*Nx+j]/cosh_kH[i*Nx+j];	
					if (runup_id==1 && depth[i][j]<Dshallow && depth[i][j]>0){
						ratio = 1;
						SINH  = kk[i][j]*(depth[i][j]+zoutIP[zz]);
					}
					factor_nonifft[i][j].Re 	= phihat_i[i*Nx+j].Re*ratio*(-kk[i][j]*tanh_kH[i*Nx+j]);
					factor_nonifft[i][j].Im 	= phihat_i[i*Nx+j].Im*ratio*(-kk[i][j]*tanh_kH[i*Nx+j]);
					tempdtPhi_nonifft[i][j].Re 	= dtphihat[i][j].Re*ratio;	
					tempdtPhi_nonifft[i][j].Im	= dtphihat[i][j].Im*ratio;
					tempdzPhi_nonifft[i][j].Re  = phihat_i[i*Nx+j].Re*kk[i][j]*SINH;			
					tempdzPhi_nonifft[i][j].Im  = phihat_i[i*Nx+j].Im*kk[i][j]*SINH;
					
					tempdxPhi_nonifft[i][j].Re	= -dominfo.kx[j]*phihat_i[i*Nx+j].Im*ratio;	
					tempdxPhi_nonifft[i][j].Im	= dominfo.kx[j]*phihat_i[i*Nx+j].Re*ratio;
					tempdyPhi_nonifft[i][j].Re	= -dominfo.ky[i]*phihat_i[i*Nx+j].Im*ratio;	
					tempdyPhi_nonifft[i][j].Im	= dominfo.ky[i]*phihat_i[i*Nx+j].Re*ratio;
					
					tempPhi_nonifft[i][j].Re	= phihat_i[i*Nx+j].Re*ratio;	
					tempPhi_nonifft[i][j].Im	= phihat_i[i*Nx+j].Im*ratio;
				}//endfor x(j)
			}

			ifftw_2d_c2r(factor,factor_nonifft,Nx,Ny);
			ifftw_2d_c2r(tempdtPhi_ifft,tempdtPhi_nonifft,Nx,Ny);
			ifftw_2d_c2r(tempdzPhi,tempdzPhi_nonifft,Nx,Ny);
			ifftw_2d_c2r(tempdxPhi_ifft,tempdxPhi_nonifft,Nx,Ny);
			ifftw_2d_c2r(tempdyPhi_ifft,tempdyPhi_nonifft,Nx,Ny);
			ifftw_2d_c2r(tempPhi_ifft,tempPhi_nonifft,Nx,Ny);
			
			for (int i=0;i<Ny;i++){					
				for (int j=0;j<Nx;j++) {
					if ((zoutIP[zz]>eta[i][j]) || (zoutIP[zz]<dominfo.bathy.data[i][j])) {
						tempdtPhi[i][j] = NAN;
						tempdzPhi[i][j] = NAN;
						tempdxPhi[i][j] = NAN;
						tempdyPhi[i][j] = NAN;	
						tempPhi[i][j]   = NAN;									
					} 
					else { 
						tempdtPhi[i][j] = tempdtPhi_ifft[i][j]+dteta[i][j]*factor[i][j];
						tempdxPhi[i][j] = tempdxPhi_ifft[i][j]+gradx_H[i][j]*factor[i][j]+gradx_D[i][j]*tempdzPhi[i][j];
						tempdyPhi[i][j] = tempdyPhi_ifft[i][j]+grady_H[i][j]*factor[i][j]+grady_D[i][j]*tempdzPhi[i][j];
						tempPhi[i][j]   = tempPhi_ifft[i][j];	
					}											
				}//end for j
			}
			
			int xo = 0;
			int yo = 0;
			for (int yo=0;yo<NyoutIP;yo++){					
				for (int xo=0;xo<NxoutIP;xo++){
					Phi[yo*NxoutIP+xo][zz]   = tempPhi[indy1+yo][indx1+xo];
					dtPhi[yo*NxoutIP+xo][zz] = tempdtPhi[indy1+yo][indx1+xo];
					dzPhi[yo*NxoutIP+xo][zz] = tempdzPhi[indy1+yo][indx1+xo];
					dxPhi[yo*NxoutIP+xo][zz] = tempdxPhi[indy1+yo][indx1+xo];
					dyPhi[yo*NxoutIP+xo][zz] = tempdyPhi[indy1+yo][indx1+xo];
					Ptot[yo*NxoutIP+xo][zz]  = -dtPhi[yo*NxoutIP+xo][zz]-grav*zoutIP[zz]-0.5*(pow(dxPhi[yo*NxoutIP+xo][zz],2)+pow(dyPhi[yo*NxoutIP+xo][zz],2)+2*dzPhi[yo*NxoutIP+xo][zz]);				
				}
			}		
		}//end zz
		
		free_2darray_complex(factor_nonifft,Nx,Ny);
		free_2darray_complex(tempdtPhi_nonifft,Nx,Ny);
		free_2darray_complex(tempdxPhi_nonifft,Nx,Ny);
		free_2darray_complex(tempdyPhi_nonifft,Nx,Ny);
		free_2darray_complex(tempdzPhi_nonifft,Nx,Ny);
		free_2darray_complex(tempPhi_nonifft,Nx,Ny);

		free_2darray(factor,Nx,Ny);
		free_2darray(tempdtPhi_ifft,Nx,Ny);
		free_2darray(tempdxPhi_ifft,Nx,Ny);
		free_2darray(tempdyPhi_ifft,Nx,Ny);
		free_2darray(tempPhi_ifft,Nx,Ny);
		
		free_2darray(tempdzPhi,Nx,Ny);
		free_2darray(tempdtPhi,Nx,Ny);
		free_2darray(tempdxPhi,Nx,Ny);
		free_2darray(tempdyPhi,Nx,Ny);
		free_2darray(tempPhi,Nx,Ny);	
		
		delete[] cosh_kDZ;
		delete[] sinh_kDZ;	
	}
	gsl_spline_free (splineHm0);
    gsl_interp_accel_free (accHm0);
    gsl_spline_free (splineHp0);
    gsl_interp_accel_free (accHp0);
    gsl_spline_free (splineDZm0); 
    gsl_interp_accel_free (accDZm0);
    gsl_spline_free (splineDZp0);
    gsl_interp_accel_free (accDZp0);
    
    // compute fluid acceleration	
    for (int i=0;i<NxyoutIP;i++) {
		for (int j=0; j<NzoutIP; j++) {
			dt_dxPHI[i][j] = (dxPhi[i][j]-Velx[i][j])/(dominfo.dt);
			dt_dyPHI[i][j] = (dyPhi[i][j]-Vely[i][j])/(dominfo.dt);
			dt_dzPHI[i][j] = (dzPhi[i][j]-Velz[i][j])/(dominfo.dt);
			Velx[i][j] = dxPhi[i][j];
			Vely[i][j] = dyPhi[i][j];
			Velz[i][j] = dzPhi[i][j];	
		}		
	}	
	
	// save elevation data at overlay zone
	double** eta_cut = declare_2darray(NxoutIP,NyoutIP);
	for (int iy=0;iy<NyoutIP;iy++){
		for (int ix=0;ix<NxoutIP;ix++){
			eta_cut[iy][ix] = eta[indy1+iy][indx1+ix];
		}
	}
	
	if (strcmp(kinematic_type,"boundary")==0){//interpolation for HAWASSI record
		//interpolation for surface elevation and surface potential
		double**  etar = declare_2darray(Nx,Ny);
		double**  phir = declare_2darray(Nx,Ny);
		double** etarec= declare_2darray(Nxrec,Nyrec);
		double** phirec= declare_2darray(Nxrec,Nyrec);
		
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				phir[i][j]=phi[i][j]/potscale;
			}
		}
		
		#pragma omp parallel for
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				etar[i][j]=eta[i][j]/xscale;
			}
		}
				
		fun_interp2(Nx,Ny,dominfo.x,dominfo.y,etar,Nxrec,Nyrec,xrec,yrec,etarec);
		fun_interp2(Nx,Ny,dominfo.x,dominfo.y,phir,Nxrec,Nyrec,xrec,yrec,phirec);
		
		//interpolation for 4 boundaries
		double**  Lphi= declare_2darray(NzoutIP,NyoutIP);
		double**  Rphi = declare_2darray(NzoutIP,NyoutIP);
		double**  Fphi = declare_2darray(NzoutIP,NxoutIP);
		double**  Bphi = declare_2darray(NzoutIP,NxoutIP);		
		
		#pragma omp parallel for
		for (int iz=0;iz<NzoutIP;iz++){//x=0 and x=end
			for (int iy=0;iy<NyoutIP;iy++){
				Lphi[iy][iz]=Phi[iy*NxoutIP][iz]/potscale;	
				Rphi[iy][iz]=Phi[iy*NxoutIP+(NxoutIP-1)][iz]/potscale;
			}
		} 
		
		#pragma omp parallel for
		for (int iz=0;iz<NzoutIP;iz++){// y=0 and y=end
			for (int ix=0;ix<NxoutIP;ix++){
				Fphi[ix][iz]=Phi[ix][iz]/potscale;
				Bphi[ix][iz]=Phi[(NyoutIP-1)*NxoutIP+ix][iz]/potscale;
			}
		} 		
			
		double**  Lphirec0  = declare_2darray(NzoutIP,Nyrec);
		double**  Rphirec0  = declare_2darray(NzoutIP,Nyrec);
		double**  Fphirec0  = declare_2darray(NzoutIP,Nxrec);
		double**  Bphirec0  = declare_2darray(NzoutIP,Nxrec);
		
		fun_interp2(NzoutIP,NyoutIP,zoutIP,youtIP,Lphi,NzoutIP,Nyrec,zoutIP,yrec,Lphirec0);
		fun_interp2(NzoutIP,NyoutIP,zoutIP,youtIP,Rphi,NzoutIP,Nyrec,zoutIP,yrec,Rphirec0);		
		fun_interp2(NzoutIP,NxoutIP,zoutIP,xoutIP,Fphi,NzoutIP,Nxrec,zoutIP,xrec,Fphirec0);
		fun_interp2(NzoutIP,NxoutIP,zoutIP,xoutIP,Bphi,NzoutIP,Nxrec,zoutIP,xrec,Bphirec0);
		
		int Nzrec=20;
		double*  zrec = new double[Nzrec];
		double** Lphirec  = declare_2darray(Nzrec,Nyrec);
		double** Rphirec  = declare_2darray(Nzrec,Nyrec);
		double** Fphirec  = declare_2darray(Nzrec,Nxrec);
		double** Bphirec  = declare_2darray(Nzrec,Nxrec);	
		
		//left and right boundary
		int Nzout_nonan;		
		for (int i=0;i<Nyrec;i++){
			Nzout_nonan=0;	
			while (!isnan(Lphirec0[i][Nzout_nonan])){
				Nzout_nonan++;
			}
			for (int ig=0;ig<Nzrec;ig++){//x=0
				zrec[ig] = -depthrec[i][0]+GQ[ig]*(etarec[i][0]+depthrec[i][0]);	
			}
			fun_interp1(Nzout_nonan,zoutIP,Lphirec0[i],Nzrec,zrec,Lphirec[i]);//left
			
			Nzout_nonan=0;
			while (!isnan(Rphirec0[i][Nzout_nonan])){
				Nzout_nonan++;
			}			
			for (int ig=0;ig<Nzrec;ig++){
				zrec[ig] = -depthrec[i][Nxrec-1]+GQ[ig]*(etarec[i][Nxrec-1]+depthrec[i][Nxrec-1]);	
			}
			fun_interp1(Nzout_nonan,zoutIP,Rphirec0[i],Nzrec,zrec,Rphirec[i]);//right								
		}
		
		//front and back boundary
		for (int j=0;j<Nxrec;j++){
			Nzout_nonan=0;
			while (!isnan(Fphirec0[j][Nzout_nonan])){
				Nzout_nonan++;
			}	
			for (int ig=0;ig<Nzrec;ig++){
				zrec[ig] = -depthrec[0][j]+GQ[ig]*(etarec[0][j]+depthrec[0][j]);
			}
			fun_interp1(Nzout_nonan,zoutIP,Fphirec0[j],Nzrec,zrec,Fphirec[j]);//front
			
			Nzout_nonan=0;
			while (!isnan(Bphirec0[j][Nzout_nonan])){
				Nzout_nonan++;
			}	
			for (int ig=0;ig<Nzrec;ig++){
				zrec[ig] = -depthrec[Nyrec-1][j]+GQ[ig]*(etarec[Nyrec-1][j]+depthrec[Nyrec-1][j]);
			}
			fun_interp1(Nzout_nonan,zoutIP,Bphirec0[j],Nzrec,zrec,Bphirec[j]);//back	
		}	
		
		//writing to Hawassi record
		UDW2p_write_hw2D(filename,Nxrec,Nyrec,Nzrec,etarec,phirec,Lphirec,Rphirec,Fphirec,Bphirec);
		
		if (it==200){
			
			/*FILE* gp0=popen("gnuplot -persistent","w");
			fun_plot_2d_trarray(gp0,yrec,xrec,etarec,Nyrec,Nxrec);
			fflush(gp0);
			pclose(gp0);
			
			FILE* gp1=popen("gnuplot -persistent","w");
			fun_plot_2d_trarray(gp1,yrec,xrec,phirec,Nyrec,Nxrec);
			fflush(gp1);
			pclose(gp1);
			
			FILE* gp2=popen("gnuplot -persistent","w");
			fun_plot_2d_array(gp2,yrec,zoutIP,Lphirec0,Nyrec,NzoutIP);
			fflush(gp2);
			pclose(gp2);
			
			FILE* gp3=popen("gnuplot -persistent","w");
			fun_plot_2d_array(gp3,yrec,zoutIP,Rphirec0,Nyrec,NzoutIP);
			fflush(gp3);
			pclose(gp3);
			
			FILE* gp4=popen("gnuplot -persistent","w");
			fun_plot_2d_array(gp4,xrec,zoutIP,Fphirec0,Nxrec,NzoutIP);
			fflush(gp4);
			pclose(gp4);
			
			FILE* gp5=popen("gnuplot -persistent","w");
			fun_plot_2d_array(gp5,xrec,zoutIP,Bphirec0,Nxrec,NzoutIP);
			fflush(gp5);
			pclose(gp5);*/
			
		}		
		
		free_2darray(etar,Nx,Ny);
		free_2darray(phir,Nx,Ny);
		free_2darray(etarec,Nxrec,Nyrec);
		free_2darray(phirec,Nxrec,Nyrec);
		
		free_2darray(Lphi,NzoutIP,NyoutIP);
		free_2darray(Rphi,NzoutIP,NyoutIP);
		free_2darray(Fphi,NzoutIP,NxoutIP);
		free_2darray(Bphi,NzoutIP,NxoutIP);	
		
		free_2darray(Lphirec0,NzoutIP,Nyrec);	
		free_2darray(Rphirec0,NzoutIP,Nyrec);	
		free_2darray(Fphirec0,NzoutIP,Nxrec);	
		free_2darray(Bphirec0,NzoutIP,Nxrec);
		
		delete[] zrec;	
		free_2darray(Lphirec,Nzrec,Nyrec);	
		free_2darray(Rphirec,Nzrec,Nyrec);	
		free_2darray(Fphirec,Nzrec,Nxrec);	
		free_2darray(Bphirec,Nzrec,Nxrec);
	}
	else {
	
		#pragma omp parallel sections num_threads(4)
		{
			#pragma omp section
			{
				for (int i=0;i<NyoutIP;i++){
					fwrite(eta_cut[i],sizeof(double),NxoutIP,dateta);//save elevation at overlay
				}
				
				for (int kk=0;kk<Ny;kk++){
					fwrite(phitemp[kk],sizeof(double),Nx,datPhi);//save surface potential at whole domain	
				}
			}
			#pragma omp section
			{		
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(Ptot[i*NxoutIP+j],sizeof(double),NzoutIP,datP);//save pressure
					}
				}			
				
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dtPhi[i*NxoutIP+j],sizeof(double),NzoutIP,datP);//save pressure
					}
				}
			}		
			
			#pragma omp section 		
			{
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dxPhi[i*NxoutIP+j],sizeof(double),NzoutIP,datvel);
					}
				}
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dyPhi[i*NxoutIP+j],sizeof(double),NzoutIP,datvel);
					}
				}
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dzPhi[i*NxoutIP+j],sizeof(double),NzoutIP,datvel);
					}
				}
			}
			
			#pragma omp section 
			{
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dt_dxPHI[i*NxoutIP+j],sizeof(double),NzoutIP,datacc);//save acceleration
					}
				}
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dt_dyPHI[i*NxoutIP+j],sizeof(double),NzoutIP,datacc);//save acceleration
					}
				}
				for (int i=0;i<NyoutIP;i++){
					for (int j=0;j<NxoutIP;j++){
						fwrite(dt_dzPHI[i*NxoutIP+j],sizeof(double),NzoutIP,datacc);//save acceleration
					}
				}
			}
		}
	}
	
	// freeing memory		
	delete[] (gam_Hmin0);
	delete[] (gam_Hplus0);
	delete[] (gam_DZmin0);
	delete[] (gam_DZplus0);
	
	free_2darray(H,Nx,Ny);
	free_2darray(gam_Hmin,Nx,Ny);
	free_2darray(gam_Hplus,Nx,Ny);
	free_2darray(gradx_H,Nx,Ny);
	free_2darray(grady_H,Nx,Ny);
	free_2darray(gradx_D,Nx,Ny);
	free_2darray(grady_D,Nx,Ny);
	free_2darray(eta_cut,NxoutIP,NyoutIP);	
	
	free_2darray(phi,Nx,Ny);
	free_2darray(phitemp,Nx,Ny);
	free_2darray(dtphi,Nx,Ny);
	free_2darray_complex(phihat,Nx,Ny);
	free_2darray_complex(dtphihat,Nx,Ny);
	
	free_2darray(Phi,NzoutIP,NxyoutIP);
	free_2darray(dtPhi,NzoutIP,NxyoutIP);
	free_2darray(dxPhi,NzoutIP,NxyoutIP);
	free_2darray(dyPhi,NzoutIP,NxyoutIP);
	free_2darray(dzPhi,NzoutIP,NxyoutIP);
	free_2darray(Ptot,NzoutIP,NxyoutIP);
	free_2darray(dt_dxPHI,NzoutIP,NxyoutIP);
	free_2darray(dt_dyPHI,NzoutIP,NxyoutIP);
	free_2darray(dt_dzPHI,NzoutIP,NxyoutIP);  
}
