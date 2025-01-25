void define_wall(complex* gamW_hat, complex* gamW_skew_hat) {
    Swall_hat 		 = new complex[Nx];
    
    if (strcmp(wall, "yes") == 0){
        int wallInfdir;		
        idx_wall 	 = closest(x, Xwall, Nx);
        depthXwall	 = depth[idx_wall];
        ITERdtW		 = 0;
        tprevW		 = ttrace[0];
		        
        if  (Xwall > Xinflux){
	    wallInfdir = -1;
        }
        else {
	    wallInfdir = 1;
        }	
        
        //only Area source
        double* vgx = new double[Nx];
        double  dista;
        
        ifft_real_real(Ug, vgx, Nx);
        
        dista       = x[idx_wall] - x[0];
        
        circshift(vgx, idx_wall, Nx);
        fft_real_complex(vgx, gamW_hat, Nx);
        scaling_influx(gamW_hat, 1/dx, dx, sqrt(g * depthinf));
        
        for (int j = 0; j < Nx; j++){
            gamW_skew_hat[j].Re = -omega_i(k[j], depthinf) * gamW_hat[j].Im;
            gamW_skew_hat[j].Im = omega_i(k[j], depthinf) * gamW_hat[j].Re;
        }

        if (wallInfdir == 1){
            for (int j = 0; j < Nx; j++){
                gamW_skew_hat[j].Re = -gamW_skew_hat[j].Re;
                gamW_skew_hat[j].Im = -gamW_skew_hat[j].Im;
            }
        }
        
        delete[] vgx;
    }
    else {
        set_init_complex(Swall_hat, Nx);
    }	
}

void wall_setup(double *ChiAdj_wall, double *ChiAdjL, double *ChiAdjR) {
    ChiAdjL = cf(x, Xwall, adjcoef * lambda_p, Nx);
    ChiAdjR = cf(x, Xwall - adjcoef * lambda_p, adjcoef * lambda_p, Nx);
    
    ChiAdj_wall  = new double[Nx];
    for (int j = 0; j < Nx; j++)	{
        ChiAdj_wall[j] = ChiAdj[j] * (1 - ChiAdjR[j] + ChiAdjL[j]);
    }

    delete[] ChiAdjL; 
    delete[] ChiAdjR; 			

    for (int j = 0; j < Nx; j++){
        ChiAdj[j] = ChiAdj_wall[j];
    }

    delete[] ChiAdj_wall;
}

void check_runup_dynmodel(int runup, char* dynmodel, size_t Ndim) {
    if (runup == 1){	
        if(strcmp(dynmodel, "HS1") == 0){ 
            gsl_odeiv_system sys = {rhs_HsF1, NULL, Ndim, NULL};
            ode_solver_Ham(sys);   
        }    
        else if(strcmp(dynmodel, "HS2") == 0){
            gsl_odeiv_system sys = {rhs_HsF2, NULL, Ndim, NULL};
            ode_solver_Ham(sys);  
        }
        else if(strcmp(dynmodel, "HS3") == 0){
            gsl_odeiv_system sys = {rhs_HsF3, NULL, Ndim, NULL};
            ode_solver_Ham(sys);      
        }  
        else if(strcmp(dynmodel, "HS4") == 0) {
            gsl_odeiv_system sys = {rhs_HsF4, NULL, Ndim, NULL};
            ode_solver_Ham(sys);
        }
    }
    else {
        gsl_odeiv_system sys = {rhs_shore, NULL, Ndim, NULL};
        ode_solver_Ham(sys);
    } 
}

int IsExpired(char* dueDay) {
    time_t currentTime;
    struct tm *timeInfo;
    char buffer[10];
    int value;

    time (&currentTime);
    timeInfo = localtime(&currentTime);

    strftime(buffer, sizeof(buffer), "%Y%m%d", timeInfo);

    value = strcmp(dueDay, buffer);

    if (value >= 0) return 1; // Not expired.
    else return -1; // Expired
}

bool isvalueinarray(int val, double *arr, int size) {
    int i;
    for (i = 0; i < size; i++) {
        if (arr[i] > val)
            return true; //yes
    }
    return false; //no
}

void readVector(int dim, char *path, double *vec) {
    FILE *fptr;
    fptr = fopen(path, "rb");
    
    if (fptr == NULL) {
        printf("Cannot open file \n");
        return;
    }
              
    for (int i = 0; i < dim; ++i){
        if (!fscanf(fptr, "%lf", &vec[i])) 
	break;
    }
    
    fclose(fptr);
}

void readVectorBin(int dim, char *path, double *vec) {
    FILE *infile;
    infile = fopen(path, "rb");
    
    if (infile == NULL) { 
        printf("Error opening file\n"); 
        return; 
    }       
      
    for (int i = 0; i < dim; i++) { 
	fread(&vec[i], sizeof(double), 1, infile);
    }
    
    fclose(infile);
}

double search_v(int n, double* array, search_type maxmin) {
    double temp = array[0];
    double vm;
    
    switch(maxmin) {
        case MAXIMUM:
            for (int i = 0; i < n; i++) {
                if(array[i] > temp) {
	            temp = array[i];
                }
            }
            break;
        case MINIMUM:
            for (int i = 0; i < n; i++) {
                if(array[i] < temp) {
	            temp = array[i];
                }
            }
            break;
        default:
            printf("Invalid calculation type\n");
            break;
    }
    
    return vm = temp;
}

double search_m(int row, int col, double** array, search_type maxmin) {
    double mm;
    double temp = array[0][0]; 
    
    switch(maxmin) {
        case MAXIMUM:
            for (int i = 0; i < row; i++) { 
                for (int j = 0; j < col; j++) { 
                    if (array[i][j] > temp) { 
                        temp = array[i][j];
                    } 
                } 
            }
            break;
        case MINIMUM:
            for (int i = 0; i < row; i++) { 
                for (int j = 0; j < col; j++) { 
                    if (array[i][j] < temp) { 
                        temp = array[i][j];
                    } 
                } 
            } 
            break;
        default:
            printf("Invalid calculation type\n");
            break;
    }
    
    return mm = temp;
}

int countlines(FILE *fp) { 
    char line[256];
    int i = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        i++;
    } 
    
    rewind(fp);
    
    return i; 
}

void set_memory() {
    k		= new double[Nx];
    dampchar	= new double[Nx];
    gamX	= new complex[Nx];
    spatgamX	= new double[Nx];
    Ug		= new double[Nx];
    aal		= new double[Nx];
    ww		= new double[n];
    spinflux	= new double[n];	
    insig 	= new double[n];
    insig_skew  = new double[n];
    insig_hat 	= new complex[n];
    kom  	= new double[n];
    GVom 	= new double[n];
    gamX_skew 	= new complex[Nx];	
    nu_peak     = 0;
}

void set_init_complex(complex* Rf_hat, int Nx) {
    for (int i = 0; i < Nx; i++){
        Rf_hat[i].Re = 0;
        Rf_hat[i].Im = 0;
    } 
}

void save_all_output(char* argv2, char* argv1) {
    logfile(argv2, argv1, n);	
}

void free_memory() {
    delete[] wave; 
    delete[] x;
    delete[] k; 
    delete[] dampchar; 
    delete[] gamX;
    delete[] spatgamX; 
    delete[] Ug; 
    delete[] ww;
    delete[] spinflux;  
    delete[] aal;	
}

double prod_normal(complex A, complex B, num_type type) {
    double prod;
    switch (type) {
        case REAL:
            prod = (A.Re * B.Re) - (A.Im * B.Im);
            break;
        case IMAG:
            prod = (A.Re * B.Im) + (A.Im * B.Re);
            break;
        default:
            printf("Invalid calculation type\n");
            break;
    }
    return prod;
}

double prod_fftw(complex A, double B_Re, double B_Im, num_type type){
    double prod;
    switch (type) {
        case REAL:
            prod = (A.Re * B_Re) - (A.Im * B_Im);
            break;
        case IMAG:
            prod = (A.Re * B_Im) + (A.Im * B_Re);
            break;
        default:
            printf("Invalid calculation type\n");
            break;
    }
    return prod;
}

double prod_1ik(double A_Im, double B_Re, double B_Im, num_type type){
    double prod;
    double A_Re = 0;
    switch (type) {
        case REAL:
            prod = (A_Re * B_Re) - (A_Im * B_Im);
            break;
        case IMAG:
            prod = (A_Re * B_Im) + (A_Im * B_Re);
            break;
        default:
            printf("Invalid calculation type\n");
            break;
    }
    return prod;
}

void set_initial_ode(double* y, int N_x) {
    int nx = N_x / 4;
    for (int i = 0; i < (2 * nx); i++) {
        y[i] = 0;//etahat.Re and etahat.Im are zero
        y[i + (2 * nx)] = 0;//uhat.Re and uhat.Im are zero
    }

    //if initial value exist
    if (strcmp(initial, "zero") != 0) {
        double*  initwave  = new double[nx];
        complex* init_hat  = new complex[nx];
        
        if (strcmp(initial, "user_defined") == 0) {			
            double*  init_u    = new double[nx];
            complex* init_uhat = new complex[nx];
            
            for (int j = 0; j < Nx; j++) {
	        init_u[j]   = temp_u[j] * (1 - dampchar[j]);
            }	
            
            fft_real_complex_nx(init_u, init_uhat);
            
            delete[] temp_u; 
            delete[] init_u;
            
            for (int j = 0; j < nx; j++) {
                y[j + 2*nx] = init_uhat[j].Re;
                y[j + 3*nx] = init_uhat[j].Im;
            }
            
            delete[] init_uhat;
        }
        
        for (int j = 0; j < Nx; j++) {
	    initwave[j] = temp_wave[j] * (1 - dampchar[j]);
        }	
        
        fft_real_complex_nx(initwave, init_hat);
        delete[] temp_wave; 
        delete[] initwave;
        
        for (int j = 0; j < nx; j++) {
            y[j]       = init_hat[j].Re;
            y[j + nx]  = init_hat[j].Im;
        }
        
        delete[] init_hat;
    }

    if (runup == 0) {
        complex* yhat 	  = new complex[nx];
        double*  ift_yhat = new double[nx];					
        for (int j = 0; j < nx; j++){
            yhat[j].Re = y[j];
            yhat[j].Im = y[j + nx];
        }
        
        ifft_complex_real_nx(yhat, ift_yhat);

        double Hhere, yr;
        double*  y_tot  = new double[nx];
        complex* y_temp = new complex[nx];

        for (int j = 0; j < nx; j++){
            Hhere = ift_yhat[j] + depth[j];
            if (Hhere < Hmin_interp){
	        yr = -depth[j];
            }
            else {
	        yr = 0;
            }
            y_tot[j] = ift_yhat[j] + yr;
        }
        
        fft_real_complex_nx(y_tot, y_temp);
        
        for (int j = 0; j < nx; j++){
            y[j]      = y_temp[j].Re;
            y[j + nx] = y_temp[j].Im;
        }	
        
        delete[] yhat;
        delete[] ift_yhat;
        delete[] y_tot;
        delete[] y_temp;
    }		
}

void save_elevation(char str[20], double N1, double N2, double N4) {
    char felev[256];
    
    strcpy(felev, argv2);
    strcat(felev, "/Hawassi_Elev_");
    strcat(felev, str);
    
    dateta = fopen(felev, "wb");

    fwrite(&N4, sizeof(double), 1, dateta);
    fwrite(&N1, sizeof(double), 1, dateta);
    fwrite(&N2, sizeof(double), 1, dateta);
    
    if (N1 == Ntrest){
        fwrite(trest, sizeof(double), N1, dateta);
    }
    else {
	fwrite(toutIP[ip - 1], sizeof(double), NtoutIP, dateta);
    }
    
    fwrite(xoutIP, sizeof(double), NxoutIP, dateta);
}

void save_pressure(char str[20], double N1, double N2, double N3, double N5) {
    char fp[256];
    
    strcpy(fp, argv2);
    strcat(fp, "/Hawassi_Pressure_");
    strcat(fp, str);
    
    datP = fopen(fp, "wb");
			    
    fwrite(&N5, sizeof(double), 1, datP);
    fwrite(&N1, sizeof(double), 1, datP);
    fwrite(&N2, sizeof(double), 1, datP);
    fwrite(&N3, sizeof(double), 1, datP);
    
    if (N1 == Ntrest){
        fwrite(trest, sizeof(double), N1, datP);
    }
    else {
	fwrite(toutIP[ip - 1], sizeof(double), NtoutIP, datP);
    }
    
    fwrite(xoutIP, sizeof(double), NxoutIP, datP);
    fwrite(zoutIP, sizeof(double), NzoutIP, datP); 
}

void save_velocity(char str[20], double N1, double N2, double N3, double N5) {
    char fvel[256];
    
    strcpy(fvel, argv2);
    strcat(fvel, "/Hawassi_Vel_");
    strcat(fvel, str);
    
    datvel = fopen(fvel, "wb");
    
    fwrite(&N5, sizeof(double), 1, datvel);
    fwrite(&N1, sizeof(double), 1, datvel);
    fwrite(&N2, sizeof(double), 1, datvel);
    fwrite(&N3, sizeof(double), 1, datvel);
  
    if (N1 == Ntrest){
        fwrite(trest, sizeof(double), N1, datvel);
    }
    else {
	fwrite(toutIP[ip - 1], sizeof(double), NtoutIP, datvel);
    }
    
    fwrite(xoutIP, sizeof(double), NxoutIP, datvel);
    fwrite(zoutIP, sizeof(double), NzoutIP, datvel);  
}

void save_acceleration(char str[20], double N1, double N2, double N3, double N5) {
    char facc[256];
    
    strcpy(facc, argv2);
    strcat(facc, "/Hawassi_Acc_");
    strcat(facc, str);
    
    datacc = fopen(facc, "wb");
    
    fwrite(&N5, sizeof(double), 1, datacc);
    fwrite(&N1, sizeof(double), 1, datacc);
    fwrite(&N2, sizeof(double), 1, datacc);
    fwrite(&N3, sizeof(double), 1, datacc);
    
    if (N1 == Ntrest){
        fwrite(trest, sizeof(double), N1, datacc);
    }
    else {
	fwrite(toutIP[ip - 1], sizeof(double), NtoutIP, datacc);
    }
    
    fwrite(xoutIP, sizeof(double), NxoutIP, datacc);
    fwrite(zoutIP, sizeof(double), NzoutIP, datacc);
}

void print_data_gauge(FILE* dateta_gauge, FILE* dateta_t_gauge, double t, double* etaf, double* dtetaf) {
    double etag, dtetag;
    
    gsl_interp_accel *acc1  = gsl_interp_accel_alloc();
    gsl_spline *spline1     = gsl_spline_alloc(gsl_interp_cspline, Nx);				
    gsl_spline *spline2     = gsl_spline_alloc(gsl_interp_cspline, Nx);
					      
    gsl_spline_init (spline1, x, etaf, Nx);
    gsl_spline_init (spline2, x, dtetaf, Nx);
		    
    fprintf(dateta_gauge, "%f ", t);
    fprintf(dateta_t_gauge, "%f ", t);
    
    for (int j = 0; j < ngauge; j++) {
        etag 	= gsl_spline_eval(spline1, gauge[j], acc1);
        dtetag 	= gsl_spline_eval(spline2, gauge[j], acc1);
        fprintf(dateta_gauge, "%f ", etag);
        fprintf(dateta_t_gauge, "%f ", dtetag);
    }
    
    fprintf(dateta_gauge, "\n");
    fprintf(dateta_t_gauge, "\n");	
			    
    gsl_spline_free(spline1);
    gsl_spline_free(spline2);
    gsl_interp_accel_free(acc1);
}

void ode_solver_Ham(gsl_odeiv_system sys) {
    double eps_abs = 1e-5;
    double eps_rel = 1e-4;  

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step *s            = gsl_odeiv_step_alloc (T, 4 * Nx);
    gsl_odeiv_control *c         = gsl_odeiv_control_y_new(eps_abs, eps_rel);
    gsl_odeiv_evolve *e          = gsl_odeiv_evolve_alloc (4 * Nx);

    double* y = new double[4 * Nx];
    set_initial_ode(y, 4 * Nx);

    double t = ttrace[0]; 
    double h = 1e-4;   //starting step size;

    gsl_ieee_env_setup();

    tim = new double[n];		
	
    for (int i = 0; i < n; i++) {	
        double  ti = wave[i].t;
        while (t < ti) {		
            int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y);
            if (status != GSL_SUCCESS) {printf ("Time integration is failed\n"); break;}        
            tprev = t;
        }

        tim[i] = t;

        //fwriting full domain data
        double* etaf 	     = new double[Nx];
        double* uf 	     = new double[Nx];
        double* dtetaf 	     = new double[Nx];
        double* dtphif 	     = new double[Nx];
        double* phif 	     = new double[Nx];
        
        complex* eta_hatf    = new complex[Nx];
        complex* u_hatf      = new complex[Nx];
        complex* dteta_hatf  = new complex[Nx];
        complex* dtphi_hatf  = new complex[Nx];
        complex* phi_hatf    = new complex[Nx];
        
        for (int j = 0; j < Nx; j++) {
            eta_hatf[j].Re 	= y[j];
            eta_hatf[j].Im 	= y[Nx + j];
            u_hatf[j].Re 	= y[2*Nx + j];
            u_hatf[j].Im 	= y[3*Nx + j];
            dteta_hatf[j].Re	= dy[j]; 
            dteta_hatf[j].Im	= dy[j + Nx];
            
            if (k[j] == 0){
                dtphi_hatf[j].Re = 0; 
                dtphi_hatf[j].Im = 0;
                phi_hatf[j].Re   = 0; 
                phi_hatf[j].Im   = 0;				
            }
            
            else {
                dtphi_hatf[j].Re = prod_1ik(-1/k[j], dy[j + 2*Nx], dy[j + 3*Nx], REAL); 
                dtphi_hatf[j].Im = prod_1ik(-1/k[j], dy[j + 2*Nx], dy[j + 3*Nx], IMAG);
                phi_hatf[j].Re   = prod_1ik(-1/k[j], u_hatf[j].Re, u_hatf[j].Im, REAL); 
                phi_hatf[j].Im   = prod_1ik(-1/k[j], u_hatf[j].Re, u_hatf[j].Im, IMAG);
            }
        }
        
        ifft_complex_real_nx(eta_hatf, etaf);
        ifft_complex_real_nx(u_hatf, uf);
        ifft_complex_real_nx(dteta_hatf, dtetaf);
        ifft_complex_real_nx(phi_hatf, phif);

        //IF WALL EXIST
        if (strcmp(wall, "yes") == 0){
	    etaf[idx_wall] = (1 + refl_coef) * etaf[idx_wall];
        }

        //in partition file
        if (((i % NtoutIP) == 0)){
            double persen= (ip + 0.0) / (npartition) * 100.0;
            printf("%g%%\n", persen);
            ip += 1;	
            
            char str[20];
            sprintf(str, "%02d.dat", ip);

            char f1[256];
            strcpy(f1, argv2);
            strcat(f1, "/Constants_");
            strcat(f1, str);
            datconst = fopen(f1, "wb");
            
            if (ip <= npartition){
	        print_xt(datconst, NtoutIP, Nx);
            }
            else {
	        print_xt(datconst, Ntrest, Nx);
            }
            
            fflush(datconst);

            char f2[256];
            strcpy(f2, argv2);
            strcat(f2, "/Hawassi_eta_");
            strcat(f2, str);
            dateta_full = fopen(f2, "wb");

            char f3[256];
            strcpy(f3, argv2);
            strcat(f3, "/Hawassi_u_");
            strcat(f3, str);
            datu_full = fopen(f3, "wb");
	        
            char f4[256];
            strcpy(f4, argv2);
            strcat(f4, "/Hawassi_phi_");
            strcat(f4, str);
            datphi_full = fopen(f4, "wb");		
        }
        
        fwrite(etaf, sizeof(double), Nx, dateta_full);
        fwrite(uf, sizeof(double), Nx, datu_full);
        fwrite(phif, sizeof(double), Nx, datphi_full);
        
        if  ((i % NtoutIP) == (NtoutIP - 1) || (i == n)){
            fflush(dateta_full);
            fflush(datu_full);
            fflush(datphi_full);
        }		

        if (strcmp(kinematic, "yes") == 0) {
            if ((i % NtoutIP) == 0){				
                char str[20];
                sprintf(str, "%02d.dat", ip);
                
                if (ip <= npartition){
                    double N1, N2, N3, N4, N5;
                    N1 = (double) NtoutIP;
                    N2 = (double) NxoutIP;
                    N3 = (double) NzoutIP;
                    N4 = N1 * N2;
                    N5 = N1 * N2 * N3;
                    
                    //save header data for first partition
                    save_elevation(str, N1, N2, N4);
                    save_pressure(str, N1, N2, N3, N5);
                    save_velocity(str, N1, N2, N3, N5);
                    save_acceleration(str, N1, N2, N3, N5);
                }
                else {
                    double N1, N2, N3, N4, N5;
                    N1 = (double) Ntrest;
                    N2 = (double) NxoutIP;
                    N3 = (double) NzoutIP;
                    N4 = N1 * N2;
                    N5 = N1 * N2 * N3;
                    
                    //save header data for first partition
                    save_elevation(str, N1, N2, N4);
                    save_pressure(str, N1, N2, N3, N5);
                    save_velocity(str, N1, N2, N3, N5);
                    save_acceleration(str, N1, N2, N3, N5);
                }		
	    }
	    
            kinematic_modul(t, etaf, u_hatf, dtetaf, dtphi_hatf);
            
            if  ((i % NtoutIP) == (NtoutIP - 1) || (i == n - 1)){
                fflush(dateta);
                fflush(datP);
                fflush(datvel);
                fflush(datacc);
            }		
        }

        if (ngauge > 0){
            if (i == 0){
                char fg1[256];
                strcpy(fg1, argv2);
                strcat(fg1, "/Hawassi_eta_gauge.txt");
                dateta_gauge = fopen(fg1, "w");	
	                
                char fg2[256];
                strcpy(fg2,argv2);
                strcat(fg2, "/Hawassi_dteta_gauge.txt");
                dateta_t_gauge = fopen(fg2, "w");
                
                char fg3[256];
                strcpy(fg3, argv2);
                strcat(fg3, "/Hawassi_dxPhi_gauge.txt");
                datphi_x_gauge = fopen(fg3, "w");
                
                char fg4[256];
                strcpy(fg4, argv2);
                strcat(fg4, "/Hawassi_dzPhi_gauge.txt");
                datphi_z_gauge = fopen(fg4, "w");
                
                fprintf(dateta_gauge, "0 ");
                fprintf(dateta_t_gauge, "0 ");
                fprintf(datphi_x_gauge, "0 ");
                fprintf(datphi_z_gauge, "0 ");
                
                for (int j = 0; j < ngauge; j++){
                    fprintf(dateta_gauge, "%f ", gauge[j]);
                    fprintf(dateta_t_gauge, "%f ", gauge[j]);
                    fprintf(datphi_x_gauge, "%f ", gauge[j]);
                    fprintf(datphi_z_gauge, "%f ", gauge[j]);
                }
                
                fprintf(dateta_gauge, "\n");
                fprintf(dateta_t_gauge, "\n");
                fprintf(datphi_x_gauge, "\n");
                fprintf(datphi_z_gauge, "\n");
                
                print_data_gauge(dateta_gauge, dateta_t_gauge, tim[i], etaf, dtetaf);
                surf_velocity_z(tim[i], datphi_x_gauge, datphi_z_gauge, phi_hatf, etaf, k, eta_hatf);				
            }		
            else {
                print_data_gauge(dateta_gauge, dateta_t_gauge, tim[i], etaf, dtetaf);
                surf_velocity_z(tim[i], datphi_x_gauge, datphi_z_gauge, phi_hatf, etaf, k, eta_hatf);							
                if (i == (n - 1)){
                    fclose(dateta_gauge);
                    fclose(dateta_t_gauge);
                    fclose(datphi_x_gauge);
                    fclose(datphi_z_gauge);
                }
            }
        }			

        if ((i == (n - 1)) && ((npartition == 1) || (Ntrest == 0))){
            double persen = 100.0;
            printf("%g%%\n", persen);
        }

        delete[] eta_hatf; 
        delete[] u_hatf; 
        delete[] phi_hatf;
        delete[] etaf; 
        delete[] uf; 
        delete[] phif;
        delete[] dteta_hatf; 
        delete[] dtphi_hatf;
        delete[] dtetaf; 
        delete[] dtphif;			
    }
    
    delete[] y; 
    delete[] dy;
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    gsl_spline_free (splineskew);
    gsl_interp_accel_free (accskew);
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);	
}

void print_xt(FILE* fpc, int N1, int Nx) {
    double NX = (double) Nx;
    double NT = (double) N1;
    
    fwrite(&NX, sizeof(double), 1, fpc);
    fwrite(x, sizeof(double), Nx, fpc);
    fwrite(&NT, sizeof(double), 1, fpc);
    
    if (N1 == Ntrest){
        fwrite(trest, sizeof(double), N1, fpc);
    }
    else {
	fwrite(toutIP[ip-1], sizeof(double), N1, fpc);
    }
    
    fwrite(depth, sizeof(double), Nx, fpc);
    fwrite(k, sizeof(double), Nx, fpc);
    fwrite(&nu_peak, sizeof(double), 1, fpc);
}

void define_bath(char* argv1) {
    if (strcmp(bath, "user_defined") == 0){
        char  input_bath[128];
        read_bath_file(argv1, input_bath); 

        depthinf = depth[idx_inf];
        
        //update bath if runup exist
        for (int j = 0; j < Nx; j++){
            if ((depth[j]) < 0 ){
                runup = 0;
                strcpy(dynmodel, "HS2");			
                break;
            }
        }		
    } else if (strcmp(bath, "flat") == 0){
        for (int i = 0; i < Nx; i++) {
            depth[i] = depthflat;
            depthinf = depthflat;			
        }
    }
    //define bathy
    bathy = new double[Nx];
    for (int i = 0; i < Nx; i++){
        bathy[i] = -depth[i];
    }	
}

void define_gauges(char* argv1) {
    if (ngauge > 0){
        char  input_gauges[128];
        FILE* fgauge;
        sprintf(input_gauges, "%s/gauges.dat", argv1);	
        fgauge	= fopen(input_gauges, "r");
        
        if(fgauge == NULL){
            printf("Can not open gauges.dat file (one column and ngauge rows data describing the gauge positions)!!\n");
            exit(0);
        }
        
        gauge = new double[ngauge];
        
        for (int j = 0; j < ngauge; j++){
	    fscanf(fgauge, "%lf \n", &gauge[j]);
        }
        
        fclose(fgauge);				
    }
}

void define_friction(char* argv1) {
    if (strcmp(friction, "yes") == 0){
        char  input_cf[128];
        FILE* fcf;
        sprintf(input_cf, "%s/cf.dat", argv1);	
        fcf = fopen(input_cf, "r");
        
        if(fcf == NULL){
            printf("Can not open cf.dat file!!\n");
            exit(0);
        }
        
        for (int j = 0; j < Nx; j++){
	    fscanf(fcf,"%lf \n",&Cf[j]);
        }
        
        fclose(fcf);
    }    
}

void define_wind(char* argv1) {
    if (strcmp(wind, "yes") == 0) {
        char  input_wd[128];
        FILE* fw;
        sprintf(input_wd, "%s/wind.dat", argv1);	
        fw = fopen(input_wd, "r");
        
        if(fw == NULL){
            printf("Can not open wind.dat file!!\n");
            exit(0);
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < Nx; j++){
	        fscanf(fw, "%lf ", &Cw[i*Nx + j]);
            }
            
            fscanf(fw, "\n");
        }
        
        fclose(fw);
    }	
}

void data_saving_break(FILE* fbr, double* dataBreak, int Nx) {	
    int j; 
    int ncb;
    ncb = int(dataBreak[1]);
    
    fprintf(fbr, "%lf %d ", dataBreak[0], ncb);
    for (j = 0; j < ncb; j++) {
        fprintf(fbr, "%d ", int(dataBreak[j + 2]));
    }
    
    fprintf(fbr, "\n");
}

void print_data(double* eta, double t, FILE* fp, int en) {
    fprintf(fp, "%g ", t);
    
    for(int i = 0; i < en; i++){
	fprintf(fp, "%g ", eta[i]);
    }
    
    fprintf(fp, "\n");
}

void print_influx(FILE* fin, marin* wave, int n) {
    for (int i = 0; i < n; i++){
        fprintf(fin, "%lf %lf\n", wave[i].t, wave[i].eta);
    }
}

void print_const_data(FILE* fpc, int n) {
    fprintf(fpc, "%d %d \n", Nx, n);
    fprintf(fpc, "%d %s %s\n", spaceinterpnum, wavename, dynmodel);
}

marin* read_data(FILE* fp, int* n) {
    Vmarin p_tmp; 
    double time;
    double signal;

    int step = 0;
    while (fscanf(fp, "%lf %lf", &time, &signal) != EOF)
    {
        marin wave_temp;
        wave_temp.t = time;
        wave_temp.eta = signal;
        p_tmp.push_back(wave_temp);
        step++;
    }
    
    (*n) = p_tmp.size(); // n <- total particle number   
 
    marin* orgwave = new marin[(*n)];
    for(int i = 0; i < (*n); i++){
        orgwave[i] = p_tmp[i];
    }
    
    p_tmp.clear();
    
    return orgwave;
}

double my_f(double x, void *params) {
    fsolve_p *p = (fsolve_p *) params;
    
    double om = p -> w;
    double depth1 = p -> d;  
    
    return om - omega_i(x, depth1);
}
     
double my_df (double x, void *params) {
    fsolve_p *p = (fsolve_p *) params;
    double depth1 = p -> d;  
    
    return -group_velocity_i(x, depth1);
}
     
void my_fdf (double x, void *params, double *f, double *df) {
    fsolve_p *p = (fsolve_p *) params;
    double om = p -> w;
    double depth1 = p -> d;
	    
    *f = om - omega_i(x, depth1);
    *df = -group_velocity_i(x, depth1);
}    
 
double invers_omega(double w1, double depth1) {
    int status;
    double value = 0.;
    int iter = 0, max_iter = 1000;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0, x = 0.; 
    gsl_function_fdf FDF;
      
    fsolve_p params = {w1, depth1};
     
    FDF.f = &my_f;
    FDF.df = &my_df;
    FDF.fdf = &my_fdf;
    FDF.params = &params;
     
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, &FDF, x);
     
    do {
        iter++;
        status = gsl_root_fdfsolver_iterate(s);
        x0 = x;
        x = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(x, x0, 0, 1e-3);
     
        if (status == GSL_SUCCESS){
            value = x;
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fdfsolver_free(s);
    
    return value;	
}
 
double mean(marin* wave, int n) {
    double signal = 0;
    
    for (int i = 0; i < n; i++){
        signal += wave[i].eta;
    }
    
    return (signal / n);
}

double mean_array(double* wave, int n) {
    double signal = 0;
    for (int i = 0; i < n; i++){
        signal += wave[i];
    }
    
    return (signal / n);
}

double variance(marin* wave, int n, double avrg) {
    double var = 0;
    for (int i = 0; i < n; i++){
        var += pow(wave[i].eta - avrg, 2);
    }
    
    return (var / n);
}

double var_array(double* wave, int n, double avrg) {
    double var = 0;
    for (int i = 0; i < n; i++){
        var += pow(wave[i] - avrg, 2);
    }
    
    return (var / n);
}

void ifft_1d(complex* ifft_result, complex* input_array, int N) {
    fftw_complex* a, *b;
    fftw_plan     plan_b;

    /* Allocate input & output array */
    a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    /* Create plans */
    plan_b = fftw_plan_dft_1d(N, a, b, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Populate input data */
    for (int i = 0; i < N; i++) {
        a[i][0] = input_array[i].Re;
        a[i][1] = input_array[i].Im;
    }

    /* Compute inverse DFT */
    fftw_execute(plan_b);
    for (int i = 0; i < N; i++) {
        ifft_result[i].Re = b[i][0] / N;
        ifft_result[i].Im = b[i][1] / N;
    }

    /* Free memory */
    fftw_destroy_plan(plan_b);
    fftw_free(a);
    fftw_free(b);
}

void fft_1d(complex* fft_result, complex* input_array, int N) {
    fftw_complex* a, *b;
    fftw_plan     plan_f;

    /* Allocate input & output array */
    a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    /* Create plans */
    plan_f = fftw_plan_dft_1d(N, a, b, FFTW_FORWARD,  FFTW_ESTIMATE);

    /* Populate input data */
    for (int i = 0; i < N; i++)	{
        a[i][0] = input_array[i].Re;
        a[i][1] = input_array[i].Im;
    }

    /* Forward DFT */
    fftw_execute(plan_f);  

    for (int i = 0; i < N; i++)	{
        fft_result[i].Re = b[i][0];
        fft_result[i].Im = b[i][1];
    }
      
    /* Free memory */
    fftw_destroy_plan(plan_f);
    fftw_free(a);
    fftw_free(b);
}

void absolute_complex(double* absolute, complex* input, int n) {
    for(int i = 0; i < n; i++) absolute[i] = sqrt(pow(input[i].Re, 2) + pow(input[i].Im, 2));
}

void abs_fft_signal(double* abs_shat, marin* wave, int n) {
    complex* signal = new complex[n];
    complex* shat   = new complex[n];
    
    for(int i = 0; i < n; i++){
        signal[i].Re = wave[i].eta;
        signal[i].Im = 0;
    }
    
    fft_1d(shat, signal, n);
    absolute_complex(abs_shat, shat,n);
    
    delete[] signal;
    delete[] shat;
}
 
double sum_array(double* array, int n) {
    double value = 0;
    for(int i = 0; i < n; i++) value += array[i];
    
    return value;
}

void inner_product_array(double* output, double* array1, double* array2, double* array3, double* array4, int n) {
    if(array3 == NULL && array4 == NULL) {
        for(int i = 0; i < n; i++) output[i] = array1[i] * array2[i];
    }
    else {
        for(int i = 0; i < n; i++) output[i] = array1[i] * array2[i] * array3[i] * array4[i];
    }
}

void inner_product_real_complex(complex* output, double* array1, complex* array2, int n) {
    for(int i = 0; i < n; i++){
        output[i].Re = array1[i] * array2[i].Re;
        output[i].Im = array1[i] * array2[i].Im;
    }
}
    
void nsqrsp_array(double* nsqrsp, double* abs_shat, int n) {
    double* shat_square = new double[n];
    inner_product_array(shat_square, abs_shat, abs_shat, NULL, NULL, n);

    double  sum = sum_array(shat_square, n);
    for(int i = 0; i < n; i++) nsqrsp[i] = sqrt(shat_square[i] / sum);
    
    delete[] shat_square;
}

void freqspace(double* freq_space, double len, int N) {
    double delta_om_k	= 2 * Pi / len;
    double fmin 	= -0.5 * 2 * Pi / (len / (N - 1));
    double fmax;
    
    if (((N / 2) * 2 == N)) {
        fmax = -fmin - delta_om_k;
    }
    else {
	fmax = -fmin;
    }
    
    double d_update = (fmax - fmin) / (N - 1);

    for(int i = 0; i < int(N / 2); i++){
	freq_space[i]	= fmin + (N/2 + i) * d_update;
    }
    
    int N_i	= 0;
    
    for(int i = int (N / 2); i < N; i++) {
        freq_space[i]	= fmin + N_i * d_update;
        N_i++;
    }
}

void filter(int ord, double *a, double *b, int np, double *x, double *y) {
    int i, j;
    y[0] = b[0] * x[0];
    for (i = 1; i < ord; i++) {
        y[i] = 0.0;
        for (j = 0; j < i + 1; j++) {
            y[i] = y[i] + b[j] * x[i - j];
        }
        for (j = 0; j < i; j++) {
            y[i] = y[i] - a[j + 1] * y[i - j - 1];
        }
    }
    for (i = ord; i < np; i++) {
        y[i] = 0.0;
        for (j = 0; j < ord; j++){
            y[i] = y[i] + b[j] * x[i - j];
        }
        for (j = 0; j < ord - 1; j++){
            y[i] = y[i] - a[j + 1] * y[i - j - 1];
        }
    }
} 

void set_coef_filter(double* a, double* b) {
    int spsmooth = 1;
    for(int i = 0; i < (spsmooth); i++){
        b[i] = 1. / (spsmooth);
        a[i] = 0;
    }
    a[0] = 1;
}

void circshift(double* array, int shift, int N) {
    int INDX = fabs(shift);
    double temp;
    if (shift > 0){
	for(int i = 0; i < INDX; i++) {
	    //CircularShiftRight
	    temp = array[N - 1];
            for(int j = N - 1; j > 0; j--){
                array[j] = array[j - 1];
            }
            array[0] = temp;
	} 
    }
    else {
	for(int i = 0; i < INDX; i++) {
	    //CircularShiftLeft
	    temp = array[0];
            for(int j = 0; j < N - 1; j++){
	        array[j] = array[j + 1];
            }
            array[N - 1] = temp;
	}
    }
}

void nsqrspec(marin* wave, double* ww, double* spinflux, double* varinflux, int n) {
    int spsmooth        = 1;
    double* abs_shat	= new double[n];
    double  gem	        = mean(wave, n);
    (*varinflux)	= variance(wave, n, gem);

    abs_fft_signal(abs_shat, wave, n);
    nsqrsp_array(spinflux, abs_shat, n);
    double length_time	= wave[n - 1].t - wave[0].t;
    freqspace(ww, length_time, n);

    delete[] abs_shat;
}


void search_max(double const* array, int n, double *max, int *indx_max) {
    double temp =  array[0];
    for(int i = 0; i < n / 2; i++){
        if(array[i] >= temp){
            temp = array[i];
            (*indx_max) = i;
        }
    }
    (*max) = temp;
}

int idx_max(double* array, int n) {
    double temp = array[0];
    int    idx  = 0;
    
    for(int i = 0; i < n; i++) {
        if(array[i] > temp){
            temp = array[i];
            idx  = i;
        }
    }
    
    return idx;
} 

double mean_periode(double* ww, double* spinflux, int n) {
    double* temp	= new double[n];
    double* temp1	= new double[n];
    double  normal	= 0;
    double  value	= 0;
    
    for(int i = 0; i < n; i++) {
        if(i > int(n / 2)) temp[i] = ww[i]; 
        else temp[i] = 0;
        temp1[i] = temp[i] * spinflux[i];
        value   += temp1[i] * spinflux[i];
        normal  += spinflux[i] * spinflux[i];
    }
    
    return fabs(2 * value / normal);
}

void interp_all(double* t, double* y, int n, gsl_interp_accel **acc_var, gsl_spline **spline_var) {
    *acc_var = gsl_interp_accel_alloc();
    *spline_var = gsl_spline_alloc(gsl_interp_cspline, n);
    
    gsl_spline_init(*spline_var, t, y, n);
}

void set_influx(marin* wave, int n, double* nu_peak, double* ww, double* spinflux) {	
    double  varinflux;
    double  max;
    int     index_max = 0;

    nsqrspec(wave, ww, spinflux, &varinflux, n); 
    search_max(spinflux, n, &max, &index_max);	

    (*nu_peak) = fabs(ww[index_max]);
}
         
void dampzone(double* dampchar, double* x, double Bleft, double Bright, double dampL, double dampR) {
    double  hhR, hhL, hR0, hR1, hR2, hL0, hL1, hL2;
    
    for(int i = 0; i < Nx; i++) {
        hR0	    = Bright - x[i] < 0? -1 : 1;
        hR1	    = 1 + ((Bright + dampR - x[i]) < 0? -1 : 1);
        hR2	    = (1 + cos(Pi * (x[i] - Bright) / dampR)) / 2;
        hhR	    = hR0 > (hR1 < hR2? hR1 : hR2)? hR0 : (hR1 < hR2? hR1 : hR2);
        hL0	    = 1 + ((x[i] - Bleft + dampL) < 0? -1 : 1);
        hL1	    = (x[i] - Bleft) < 0? -1 : 1;
        hL2	    = (1 - cos(Pi * (x[i] - Bleft + dampL) / dampL)) / 2;
        hhL	    = hL0 < (hL1 > hL2? hL1 : hL2)? hL0 : (hL1 > hL2? hL1 : hL2);
        dampchar[i] = 1 - (hhR * hhL);
    }
}

void group_velocity(double* Ug, double* k, int Nk, double* d) {
    if (strcmp(dispersion, "SWE") == 0){
        for(int i = 0; i < Nk; i++) {
	    Ug[i] = sqrt(g * d[i]);
        }		
    }
    else {
        for(int i = 0; i < Nk; i++) {
            if(k[i] == 0) Ug[i] = sqrt(g * d[i]);
            else Ug[i] = (k[i] < 0? -1 : 1) * sqrt(g) / 2 / sqrt(k[i] * tanh(d[i] * k[i])) * (tanh(d[i] * k[i]) + k[i] * (1 - pow(tanh(d[i] * k[i]), 2)) * d[i]);
        }
    }	
}

double group_velocity_i(double k_i, double d_i) {
    double Ug_i;

    if (strcmp(dispersion, "SWE") == 0){
        Ug_i = sqrt(g * d_i);	
    }	
    else {
        if(k_i == 0) {
	    Ug_i = sqrt(g * d_i);
        }
        else {
	    Ug_i = (k_i < 0? -1 : 1) * sqrt(g) / 2 / sqrt(k_i * tanh(d_i * k_i)) * (tanh(d_i * k_i) + k_i * (1 - pow(tanh(d_i * k_i), 2)) * d_i);
        }
    }

    return Ug_i;
}

void omega(double* om, double* k, double* d, int N) {
    if (strcmp(dispersion, "SWE") == 0){
        for(int i = 0; i < N; i++) {
	    om[i] = sqrt(g * d[i]) * k[i];
        }		
    }	
    else {
        for (int i = 0; i < N; i++){
            if (k[i] == 0){
	        om[i] = 0;
            }
            else {
	        om[i] = (k[i] < 0? -1 : 1) * sqrt(g * k[i] * tanh(d[i] * k[i]));	
            }
        }
    }
}

double omega_i(double k_i, double d_i) {
    double om;

    if (strcmp(dispersion, "SWE") == 0){
        om = sqrt(g * d_i) * k_i;		
    }
    else {
        if (k_i == 0){
	    om = 0;
        }
        else {
	    om = (k_i < 0? -1 : 1) * sqrt(g * k_i * tanh(d_i * k_i));	
        }
    }	

    return om;
}

double phase_velocity_i(double k_i, double d_i) {
    double Up_i;

    if (strcmp(dispersion, "SWE") == 0){
        Up_i = sqrt(g * d_i);		
    }
    else {
        if(k_i == 0) Up_i = sqrt(g * d_i);
        else Up_i = (k_i < 0? -1 : 1) * sqrt(g * k_i * tanh(d_i * k_i)) / k_i;
    }

    return Up_i;
}

double* fun_csqr(double* k, double d, int Nx) {
    double* cpsqr = new double[Nx];
    for (int i = 0; i < Nx; i++){
        cpsqr[i] = pow(phase_velocity_i(k[i], d), 2);
    }
    
    return cpsqr;
}

void ifft_real_real(double* input, double* output, int N)
{
    complex* signal = new complex[N];
    complex* shat   = new complex[N];
    
    for(int i = 0; i < N; i++){
        signal[i].Re = input[i];
        signal[i].Im = 0;
    }
    
    ifft_1d(shat, signal, N);
    
    for(int i = 0; i < N; i++){
	output[i] = shat[i].Re;
    }
    
    delete[] signal;
    delete[] shat;
}

void ifft_fft_Hsf2 (const double* z, double* L, fftw_complex* eta_u_hat, fftw_complex* LetaM0u_hat, fftw_complex* u2_M0u2_hat, fftw_complex* damping_eta, fftw_complex* damping_u, double* dampchar, int N) {
    fftw_complex* etahat_input, *eta;
    fftw_complex* uhat, *u;
    fftw_complex* M0uhat, *M0u;

    etahat_input   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    eta            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    uhat           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    u              = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0uhat         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0u            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for(int i = 0; i < N; i++){
        etahat_input[i][0] = z[i];
        etahat_input[i][1] = z[i + N];
        uhat[i][0]  	   = z[i + 2*N];
        uhat[i][1]    	   = z[i + 3*N];
        M0uhat[i][0]	   = prod_1ik(M0[i].Im, z[i + 2*N], z[i + 3*N], REAL);
        M0uhat[i][1]	   = prod_1ik(M0[i].Im, z[i + 2*N], z[i + 3*N], IMAG);
    }
    
    fftw_execute_dft(plan_ift_nx, etahat_input, eta);
    fftw_execute_dft(plan_ift_nx, uhat, u);
    fftw_execute_dft(plan_ift_nx, M0uhat, M0u);    
    
    fftw_complex  *eta_u      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta_M0u    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *etaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *u2_M0u2    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *damp_eta   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *damp_u     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
      
    for (int i = 0; i < N; i++){//real from ift
        eta_u[i][0]    = (eta[i][0] / N) * (u[i][0] / N);
        eta_u[i][1]    = 0;
        eta_M0u[i][0]  = (eta[i][0] / N) * (M0u[i][0] / N);
        eta_M0u[i][1]  = 0;
        u2_M0u2[i][0]  = (pow(u[i][0] / N, 2) - pow(M0u[i][0] / N, 2));
        u2_M0u2[i][1]  = 0;  
        damp_eta[i][0] = dampcoef[i] * (eta[i][0] / N) * dampchar[i];
        damp_eta[i][1] = 0;
        damp_u[i][0]   = dampcoef[i] * (u[i][0] / N) * dampchar[i];
        damp_u[i][1]   = 0;
    }
    
    fftw_execute_dft(plan_ft_nx, eta_u, eta_u_hat);  
    fftw_execute_dft(plan_ft_nx, eta_M0u, etaM0u_hat);  
    fftw_execute_dft(plan_ft_nx, u2_M0u2, u2_M0u2_hat);  
    fftw_execute_dft(plan_ft_nx, damp_eta, damping_eta);
    fftw_execute_dft(plan_ft_nx, damp_u, damping_u);

    for (int i = 0; i < N; i++){
        LetaM0u_hat[i][0] = L[i] * etaM0u_hat[i][0];
        LetaM0u_hat[i][1] = L[i] * etaM0u_hat[i][1];
    }	
     
    fftw_free(etahat_input);
    fftw_free(eta);
    fftw_free(uhat);
    fftw_free(u);
    fftw_free(M0uhat);
    fftw_free(M0u);
    fftw_free(eta_M0u); 
    fftw_free(etaM0u_hat);
    fftw_free(eta_u);
    fftw_free(u2_M0u2); 
    fftw_free(damp_eta);
    fftw_free(damp_u); 
}

void ifft_fft_Hsf1(const double* z, double* L, fftw_complex* damping_eta, fftw_complex* damping_u, double* dampchar, int N) {
    fftw_complex* etahat_input, *eta;
    fftw_complex* uhat, *u;

    etahat_input   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    eta            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    uhat           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    u 	           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for(int i = 0; i < N; i++){
        etahat_input[i][0] = z[i];
        etahat_input[i][1] = z[i + N];
        uhat[i][0]         = z[i + 2*Nx];
        uhat[i][1]         = z[i + 3*Nx];
    }
  
    fftw_execute_dft(plan_ift_nx, etahat_input, eta);
    fftw_execute_dft(plan_ift_nx, uhat, u);
	
    fftw_complex  *damp_eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *damp_u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++){//real
        damp_eta[i][0] = dampcoef[i] * (eta[i][0] / N) * dampchar[i];
        damp_eta[i][1] = 0;
        damp_u[i][0]   = dampcoef[i] * (u[i][0] / N) * dampchar[i];
        damp_u[i][1]   = 0;
    }

    fftw_execute_dft(plan_ft_nx, damp_eta, damping_eta);
    fftw_execute_dft(plan_ft_nx, damp_u, damping_u);
     
    fftw_free(etahat_input);
    fftw_free(eta);
    fftw_free(uhat); 
    fftw_free(u);
    fftw_free(damp_eta);
    fftw_free(damp_u); 
}

void ifft_fft_Hsf3(const double* z, double* L, fftw_complex* LetaM0u_hat, fftw_complex* deletaH3_hat, fftw_complex* deluH3_hat, int N) {	
    fftw_complex* etahat_input, *eta;
    fftw_complex* dxu_hat, *dxu;
    fftw_complex* M0u_hat, *M0u;
    fftw_complex* LetaM0u;
      
    // Allocate input & output array 
    etahat_input   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    eta            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxu_hat        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxu            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);	
    M0u_hat        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0u            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);	
    LetaM0u        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	  
    for(int i = 0; i < N; i++){//spectral
        etahat_input[i][0] = z[i];
        etahat_input[i][1] = z[i + N];
        dxu_hat[i][0]      = prod_1ik(k[i], z[i + 2*N], z[i + 3*N], REAL);
        dxu_hat[i][1]      = prod_1ik(k[i], z[i + 2*N], z[i + 3*N], IMAG);
        M0u_hat[i][0]      = prod_1ik(M0[i].Im, z[i + 2*N], z[i + 3*N], REAL);
        M0u_hat[i][1]      = prod_1ik(M0[i].Im, z[i + 2*N], z[i + 3*N], IMAG);
    }
        
    fftw_execute_dft(plan_ift_nx, etahat_input, eta);
    fftw_execute_dft(plan_ift_nx, dxu_hat, dxu);
    fftw_execute_dft(plan_ift_nx, M0u_hat, M0u);
    fftw_execute_dft(plan_ift_nx, LetaM0u_hat, LetaM0u);	

    fftw_complex  *etaM0dxetaM0u      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *etaM0dxetaM0u_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta2_dxu           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta2_dxu_hat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta2_M0u           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta2_M0u_hat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *deletaH3           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++){//real
        etaM0dxetaM0u[i][0]  = (eta[i][0] / N) * (LetaM0u[i][0] / N);
        etaM0dxetaM0u[i][1]  = 0;
        eta2_dxu[i][0]       = pow(eta[i][0] / N, 2) * (dxu[i][0] / N);
        eta2_dxu[i][1]       = 0;
        eta2_M0u[i][0]       = pow(eta[i][0] / N, 2) * (M0u[i][0] / N);
        eta2_M0u[i][1]       = 0;
        deletaH3[i][0]       = (M0u[i][0] / N) * ((eta[i][0] / N) * (dxu[i][0] / N) + (LetaM0u[i][0] / N)) ;
        deletaH3[i][1]       = 0;
    }

    fftw_execute_dft(plan_ft_nx, etaM0dxetaM0u, etaM0dxetaM0u_hat);  
    fftw_execute_dft(plan_ft_nx, eta2_dxu, eta2_dxu_hat);  
    fftw_execute_dft(plan_ft_nx, eta2_M0u, eta2_M0u_hat); 
    fftw_execute_dft(plan_ft_nx, deletaH3, deletaH3_hat); 
	
    for(int i = 0; i < N; i++){//spectral
        deluH3_hat[i][0] = L[i] * etaM0dxetaM0u_hat[i][0] + (1/2.) * L[i] * eta2_dxu_hat[i][0] - (1/2.) * pow(k[i], 2) * eta2_M0u_hat[i][0];
        deluH3_hat[i][1] = L[i] * etaM0dxetaM0u_hat[i][1] + (1/2.) * L[i] * eta2_dxu_hat[i][1] - (1/2.) * pow(k[i], 2) * eta2_M0u_hat[i][1];
    }
     
    fftw_free(etahat_input);
    fftw_free(eta);
    fftw_free(dxu);
    fftw_free(dxu_hat);
    fftw_free(M0u);
    fftw_free(M0u_hat);	 
    fftw_free(etaM0dxetaM0u_hat);
    fftw_free(etaM0dxetaM0u);
    fftw_free(eta2_dxu); 
    fftw_free(eta2_dxu_hat);
    fftw_free(eta2_M0u); 
    fftw_free(eta2_M0u_hat);
    fftw_free(deletaH3);
    fftw_free(LetaM0u);
}

void ifft_complex_real_nx(complex* input, double* output) {
    fftw_complex* a, *b;

    a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    for(int i = 0; i < Nx; i++){
        a[i][0] = input[i].Re;
        a[i][1] = input[i].Im;
    }
    
    fftw_execute_dft(plan_ift_nx, a, b);

    for (int i = 0; i < Nx; i++){
	output[i] = b[i][0] / Nx;
    }

    fftw_free(a);
    fftw_free(b);
} 

void ifft_complex_real(complex* input, double* output, int N) {
    fftw_complex* a, *b;
    fftw_plan     plan_b;

    /* Allocate input & output array */
    a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    /* Create plans */
    plan_b = fftw_plan_dft_1d(N, a, b, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(int i = 0; i < N; i++){
        a[i][0] = input[i].Re;
        a[i][1] = input[i].Im;
    }
    fftw_execute(plan_b);

    for (int i = 0; i < N; i++){
	output[i] = b[i][0] / N;
    }
    /* Free memory */
    fftw_destroy_plan(plan_b);

    fftw_free(a);
    fftw_free(b);
} 

void fft_real_complex(double* input, complex* output, int N) {
    fftw_complex  *a, *b;
    fftw_plan     plan_f;

    /* Allocate input & output array */
    a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    /* Create plans */
    plan_f = fftw_plan_dft_1d(N, a, b,FFTW_FORWARD, FFTW_ESTIMATE);

    /* Populate input data */
    for (int i = 0; i < N; i++)	{
        a[i][0] = input[i];
        a[i][1] = 0;    
    }

    /* Forward DFT */
    fftw_execute(plan_f);  

    for (int i = 0; i < N; i++){
        output[i].Re = b[i][0];
        output[i].Im = b[i][1];
    }
      
    /* Free memory */
    fftw_destroy_plan(plan_f);
    fftw_free(a);
    fftw_free(b);
}

void fft_real_complex_nx(double* input, complex* output) {
    fftw_complex  *a, *b;

    /* Allocate input & output array */
    a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    /* Populate input data */
    for (int i = 0; i < Nx; i++)	{
	    a[i][0] = input[i];
	    a[i][1] = 0;    
    }

    /* Forward DFT */
    fftw_execute_dft(plan_ft_nx, a, b);  

    for (int i = 0; i < Nx; i++){
        output[i].Re = b[i][0];
        output[i].Im = b[i][1];
    }
      
    /* Free memory */
    fftw_free(a);
    fftw_free(b);
}

void scaling_influx(complex* gamX, double fact, double dx, double cp) {
    for(int i = 0; i < Nx; i++){
        gamX[i].Re = gamX[i].Re * fact;
        gamX[i].Im = gamX[i].Im * fact;
    }
    
    double scale = gamX[0].Re * dx;
    
    for(int i = 0; i < Nx; i++){
        gamX[i].Re = gamX[i].Re * cp / scale;
        gamX[i].Im = gamX[i].Im * cp / scale;
    }
}

int closest(double* array, double ref, int N) {
    int     index = 0;
    double  min   = 10 * (array[N / 2] - array[(N / 2) - 1]);
    for(int i = 0; i < N; i++) {
        if (fabs(array[i] - ref) < min){
            min = fabs(array[i] - ref);
            index = i;
        }
        else continue;
    }
    return index;
}

void alias(double* al, double* k, double frac) {
    int idk  = int(Nx / frac) - 1;
    double k0 = k[idk];
    for(int i = 0; i < Nx; i++){
        if (i != idk){
            al[i] = ((k0 - fabs(k[i])) / fabs(k0 - fabs(k[i])) + 1) / 2;
        }
        else {
	        al[i] = 0.5;
        }	
    }
}

void set_operator(int Ndim) {
    dy	  = new double[Ndim]; set_init_double(dy,Ndim);
    Csqr  = new double[Nx]; set_init_double(Csqr,Nx);
    L	  = new double[Nx]; set_init_double(L,Nx);
    L0    = new double[Nx]; set_init_double(L0,Nx);
    M0    = new complex[Nx]; set_init_complex(M0,Nx);
    M1    = new complex[Nx]; set_init_complex(M1,Nx);
}

void operator_L(double* L, double* k, double d) {
    for (int i = 0; i < Nx; i++){
        L[i] = pow(k[i], 2) * pow(phase_velocity_i(k[i], d), 2) / g;
    }
}

void operator_L0(double* L0, double* k) {
    for (int i = 0; i < Nx; i++){
	L0[i] = k[i] * Csqr[i] * k[i] / g;
    }
}

void operator_M0(complex* M0, double* k) {
    for (int i = 0; i < Nx; i++){
        M0[i].Re = 0;
        M0[i].Im = -k[i] * Csqr[i] / g; 
    }
}

void operator_M1(complex* M1, double* k) {
    for (int i = 0; i < Nx; i++){
        M1[i].Re = 0;
        M1[i].Im = -Csqr[i] * k[i] / g; 
    }
}
	
double trapz(double* x, double* array, int N) {
    double sum = 0.0;
    for (int i = 0; i < N; i++){
        if (i == 0 || i == N - 1) 
            sum += array[i] / 2;
        else
            sum += array[i]; 
    }   
    return sum * (x[1] - x[0]); 
}

double trapz_non(double* x, double* array, int N) {
    double sum = 0.0;
    for (int i = 0; i < N - 1; i++){
        sum += (array[i] + array[i + 1]) * (x[i + 1] - x[i]) / 2;
    }
    return sum; 
}

void cumtrapz(double* in, int n, double dt, double* out) {
    double temp   = 0;
    out[0]        = temp;
    for (int i = 1; i < n; i++){
        temp 	= temp + 0.5 * dt * (in[i - 1] + in[i]);
        out[i]  = temp;
    }
}

void gradient(double* array, double d, int n, double* grad) {
    grad[0] 	= (array[1] - array[0]) / d;
    grad[n-1] 	= (array[n - 1] - array[n - 2]) / d;
    for (int i = 1; i < n - 1; i++) {
        grad[i] = (array[i + 1] - array[i - 1]) / (2 * d);
    }
}

int valueinarraylarger(double val, double* arr, int n) {
    int i;
    for(i = 0; i < n; i++) {
        if(arr[i] > val)
            return 1;
    }
    return 0;
}

int allzeros(int* arr, int n) {   
    int temp = 1;//all element is zero
    for(int i = 0; i < n; i++) {
        if(arr[i] != 0)
            temp = 0;
            break;
    }
    return temp;
}

int allvalueless(double val, double* arr, int n) {
    double temp;
    temp = search_v(n, arr, MAXIMUM);
    
    if (temp < val) {
        return 1;
    }
    else {
	return 0;
    }
}

void sinhcosh2IP(int sizD, double* D, double Dmin, double Dplus, int Nx, double* k, double nupeak, double* gam_min, double* gam_plus) {
    double*  kappanu_D	= new double[sizD];
    double*  cosh_D	= new double[sizD];
    double*  cosh_min	= new double[sizD];
    double*  cosh_plus	= new double[sizD];
    double*  sinh_D	= new double[sizD];
    double*  sinh_min	= new double[sizD];
    double*  sinh_plus	= new double[sizD];	
		
    for (int j = 0; j < sizD; j++) {
        if (D[j] == 0) {
	    kappanu_D[j] = k[2];
        }		
        else {
	    kappanu_D[j] = invers_omega(nupeak, D[j]);
        }
        
        cosh_min[j]	= cosh(kappanu_D[j] * Dmin);
        cosh_plus[j]    = cosh(kappanu_D[j] * Dplus);
        cosh_D[j]       = cosh(kappanu_D[j] * D[j]);
        sinh_min[j]     = sinh(kappanu_D[j] * Dmin);
        sinh_plus[j]    = sinh(kappanu_D[j] * Dplus);
        sinh_D[j]       = sinh(kappanu_D[j] * D[j]);
    }

    double det[sizD], detA[sizD], detB[sizD], gminA[sizD], gminB[sizD], gplusA[sizD], gplusB[sizD];

    inner_product_array(detA, cosh_min, sinh_plus, NULL, NULL, sizD);
    inner_product_array(detB, cosh_plus, sinh_min, NULL, NULL, sizD);
    inner_product_array(gminA, sinh_plus, cosh_D, NULL, NULL, sizD);
    inner_product_array(gminB, cosh_plus, sinh_D, NULL, NULL, sizD);
    inner_product_array(gplusA, sinh_min, cosh_D, NULL, NULL, sizD);
    inner_product_array(gplusB, cosh_min, sinh_D, NULL, NULL, sizD);
	
    for (int i = 0; i < sizD; i++) {
        det[i] 	    = detA[i] - detB[i];
        gam_min[i]  = (gminA[i] - gminB[i]) / det[i];				
        gam_plus[i] = (-gplusA[i] + gplusB[i]) / det[i];
    }
}

int findfirstzero(double *arr, int size) {
    int i;
    int idx = 0;
    for (i = 0; i < size; i++) {
        if (arr[i] != 0){
	    idx += 1;
        }
        else
            return idx;
    }
    return 0;
}
 
int findfirstnegative(double *arr, int size) {
    int i;
    int idx = 0;
    for (i = 0; i < size; i++) {
        if (arr[i] >= 0){
	    idx = idx + 1;
        }
        else
            return idx;
    }
    return 0;
}

void subarrayx(int Nx, double* array, int j1, int Nxnew, double* arraynew) {	
    for (int j = 0; j < Nxnew; j++){
        arraynew[j] = array[j1 + j];	
    }
}

void subvector(double *array, int j1, int Nxnew, double *arraynew) {
    for (int i = 0; i < Nxnew; i++){
	arraynew[i] = array[j1 + i];		
    }
}

void modified_influx(marin* wave, int n, double* nu_peak, double* ww, double* spinflux, double* insig, double* insig_skew, double* kom, double* GVom) {	
    for (int j = 0; j < n; j++){
        insig[j]  = wave[j].eta;
    }

    if (strcmp(wavename, "user_defined") == 0) { 
        double* ramp 	= new double[n];
        double* endramp = new double[n];
        double  A, B1, B2;
       
        for (int j = 0; j < n; j++) {
	    A 	= ((ttrace[j] - ttrace[0]) - nTramp * Tp_peak);
	    B1 	= A / fabs(A);
	    B2 	= (1 - cos((ttrace[j] - ttrace[0]) * Pi / nTramp / Tp_peak)) / 2;
	    if (B1 > B2) {
	        ramp[j] = B1;
	    } 
	    else {
		ramp[j] = B2;
	    }			
	}
	
	for (int j = 0; j < n; j++) {
	    endramp[j] = ramp[n - 1 - j];
	    insig[j]   = ramp[j] * insig[j] * endramp[j];
	}	
	
	delete[] ramp; 
	delete[] endramp;		
    }	

    cumtrapz(insig, n, dt, insig_skew);
    fft_real_complex(insig, insig_hat, n);

    for (int i = 0; i < n; i++) {
        if (ww[i] == 0) {
	    kom[i]  = 0;
        }
        else {
	    kom[i] = invers_omega(ww[i], depthinf);
        }
        GVom[i] = group_velocity_i(kom[i], depthinf);
    }	
		
    complex* modinsighat 	= new complex[n];
    complex* insig_skewhat	= new complex[n];	
    complex* modinsig_skewhat 	= new complex[n];
    double*  modinsig		= new double[n];
    double*  modinsig_skew 	= new double[n];
	
    if (strcmp(influxing, "Point") == 0){
        inner_product_real_complex(modinsighat, GVom, insig_hat, n);	
        ifft_complex_real(modinsighat, modinsig, n);	
	        
        fft_real_complex(insig_skew, insig_skewhat, n);
        inner_product_real_complex(modinsig_skewhat, GVom, insig_skewhat, n);
        ifft_complex_real(modinsig_skewhat, modinsig_skew, n);       
       
        for (int j = 0; j < n; j++){
            insig[j]      = modinsig[j]; //modified influx for point generation
            insig_skew[j] = modinsig_skew[j];
        }        
    }
    else if  (strcmp(influxing, "AreaShort") == 0){
        for (int j = 0; j < n; j++) {
	    GVom[j]  = GVom[j] * sqrt(g * depthinf) / GVom[0];   
        }   
    }
    
    delete[] modinsighat;
    delete[] modinsig_skewhat;
    delete[] modinsig; 
    delete[] modinsig_skew;
    delete[] insig_skewhat;
}

void fftw_to_real(fftw_complex* etafftw, double* eta, int Nx) {
    for (int j = 0; j < Nx; j++) {
        eta[j] = etafftw[j][0];
    }
}

void fftw_to_complex(fftw_complex* etafftw, complex* eta, int Nx) {
    for (int j = 0; j < Nx; j++) {
        eta[j].Re = etafftw[j][0];
        eta[j].Im = etafftw[j][1];
    }
}

void set_init_CB(CB* CrestB, int n) {
    for (int j = 0; j < n; j++){
        CrestB[j].Bindex    = 0;
        CrestB[j].Ucrest    = 0;
        CrestB[j].DirProp   = 0;
        CrestB[j].Bposition = 0;
        CrestB[j].Btime     = 0;
        CrestB[j].Bcount    = 0;
    }
}

void set_init_double(double* B, int n) {
    for (int j = 0; j < n; j++){
	B[j]   = 0;
    }
}

void set_init_int(int* B, int n) {
    for (int j = 0; j < n; j++){
	B[j]   = 0;
    }
}

void sorting(double* arr, int size) {
    double temp;
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if(arr[j] < arr[i]) {
                temp 	= arr[i];
                arr[i] 	= arr[j];
                arr[j] 	= temp;
            }
        }
    }    
}

void logfile(char* argv2, char* argv1, int n) {
    char data_log[128];
    sprintf(data_log, "%s/logfile.txt", argv2);
    FILE* flog   = fopen(data_log, "w");
 
    time_t tnow;
    time(&tnow);
    
    fprintf(flog, "HAWASSI v240719 simulation has been run on %s\n", ctime(&tnow));
    fprintf(flog, "File input: %s/input.txt \n", argv1);
    fprintf(flog, "MODEL DESCRIPTION:\n"
    "Dynamic Model \t\t: %s\n"
    "Dispersion Model \t: %s\n"
    "Breaking \t\t: %s\n", dynmodel, dispersion, breaking);
    
    if (strcmp(breaking, "yes") == 0){
        fprintf(flog,"\t Initiation \t: %.2f\n"
		     "\t Termination \t: %.2f %.2f\n", parKBC, parTC, parTchar);
    }
    fprintf(flog, "Bottom friction \t: %s\n", friction);
    
    if (strcmp(friction, "yes") == 0){
	fprintf(flog, "The friction coefficient data is from %s/cf.dat\n", argv2);
    }
    fprintf(flog, "Kinematics \t\t: %s\n", kinematic);
    
    if (strcmp(kinematic, "yes") == 0){
        fprintf(flog, "\t Overlay zone \t: [%.2f,%.2f]\n", xinterv1, xinterv2);
        fprintf(flog, "\t Non-uniform vertical grid : [%.2f,%.2f]\n", zoutIP[0], zoutIP[NzoutIP - 1]);
    }

    fprintf(flog, "\n");
    fprintf(flog, "INFLUX DESCRIPTION :\n"
    "Signal type \t\t: %s\n", wavename);
    if ((strcmp(wavename, "jonswap") == 0) || (strcmp(wavename, "harmonic") == 0)){		
        fprintf(flog, "Significant wave Height (Hs) \t: %.2f\n"
		    "Peak period (Tp) \t\t: %.2f\n", Hs, Tp);
    }
    if (strcmp(wavename, "zero") != 0){
        fprintf(flog, "Derived info:\n"
	        "\t Peak frequency (nu) \t: %.2f[rad/s]\n"
	        "\t Peak wave-number (kp) \t: %.2f\n"
	        "\t Peak wave-length \t: %.2f[m]\n"
	        "\t Steepness (kp*(Hs./2)) : %.2f\n", nu_peak, kp, lambda_p, kp * (Hs_in / 2));
    }    
	
    fprintf(flog, "\n");
    
    fprintf(flog, "INITIAL WAVE CONDITIONS : \n"
	    "Initial condition : %s\n", initial);
	    
    fprintf(flog, "\n");
    
    fprintf(flog, "NUMERICAL SETTINGS : \n"
	    "Spatial interval: (%.2f,%.2f)[m]\n"
	    "Damping zone \t: %.2f[m] and %.2f[m]\n"
	    "Number of Nodes : %d\n"
	    "Grid size (dx) \t: %.2f[m]\n"
            "Cutfrac k \t: %d\n"
	    "Time interval \t: (%.2f,%.2f)[s]\n"
            "Time step (dt) \t: %.2f[s]\n", x[0], x[Nx - 1], dampL, dampR, Nx, dx, cutfracwn, tim[0], tim[n - 1], tim[1] - tim[0]);
            
    fprintf(flog, "Bathymetry \t: %s\n", bath);
    
    if (strcmp(bath, "flat") == 0){
        fprintf(flog, "\t Depth \t: %.2f[m]\n", depthflat);
    }
    else {
	fprintf(flog, "The bathymetry data is from %s/bath.dat\n", argv2);
    }
	
    fprintf(flog, "\n");
	
    if (strcmp(wavename, "zero") != 0){
        fprintf(flog, "INFLUXING : \n"
	        "Generation method : %s\n"
	        "Direction \t: %s\n"
	        "Influx position : %.2f[m]\n", influxing, propagation, Xinflux);
        if (strcmp(dynmodel, "HsF1") != 0){
	    fprintf(flog, "Nonlinear adjustment: %g*lambda_peak\n", adjcoef);
        }
    }

    fprintf(flog, "\n");
    fprintf(flog, "OUTPUT SETTING: \n"
	         "No of Partition : %d\n", npartition);

    fprintf(flog, "\n");
    fprintf(flog, "After simulation:\n"
	         "Relative time consuming: %g\n"
	         "Data saved in the folder: %s", rel, argv2);  
}

void set_operator_bathy(double nu, double* mid) {
    if (mid[0] == 3){
        Up3IP(nu);
    }
    else {
	Up2IP(nu);
    }
}

void Up2IP(double nu) {
    double D_min   = search_v(Nx, depth, MINIMUM);
    double D_plus  = search_v(Nx, depth, MAXIMUM);
    int n          = 101;
    double dD 	   = (D_plus - D_min) / (n - 1);
    
    double Cp_min, Cp_plus, Cp_D, Om_min, Om_plus, Om_D, det, kappanu_D;
    
    double* Ddepth   = new double[n];
    double* Ggam_min = new double[n];
    double* Ggam_plus= new double[n];
		
    for (int i = 0; i < n; i++) {
        Ddepth[i]   = D_min + i * dD;
        kappanu_D   = invers_omega(nu, Ddepth[i]);
        Cp_min      = pow(phase_velocity_i(kappanu_D, D_min), 2);
        Cp_plus     = pow(phase_velocity_i(kappanu_D, D_plus), 2);
        Cp_D        = pow(phase_velocity_i(kappanu_D, Ddepth[i]), 2);
        Om_min      = omega_i(kappanu_D, D_min);
        Om_plus     = omega_i(kappanu_D, D_plus);
        Om_D        = omega_i(kappanu_D, Ddepth[i]);
        det 	    = Cp_min * Om_plus - Cp_plus * Om_min;
        Ggam_min[i] = (Om_plus * Cp_D - Cp_plus * Om_D) / det;
        Ggam_plus[i]= (-Om_min * Cp_D + Cp_min * Om_D) / det;
    }   
	
    gsl_interp_accel *acc_plus = gsl_interp_accel_alloc();
    gsl_spline *spline_plus    = gsl_spline_alloc(gsl_interp_cspline, n);
	    
    gsl_interp_accel *acc_min  = gsl_interp_accel_alloc();;
    gsl_spline *spline_min     = gsl_spline_alloc(gsl_interp_cspline, n);;	  
      		  
    gsl_spline_init(spline_plus, Ddepth, Ggam_plus, n);
    gsl_spline_init(spline_min, Ddepth, Ggam_min, n);
	    
    for (int i = 0; i < Nx; i++){
        gam_plus[i] = gsl_spline_eval(spline_plus, depth[i], acc_plus);
        gam_min[i]  = gsl_spline_eval(spline_min, depth[i], acc_min);
        Csq_min[i]  = pow(phase_velocity_i(k[i], D_min), 2);
        Csq_plus[i] = pow(phase_velocity_i(k[i], D_plus), 2);
    }
		
    gsl_spline_free(spline_plus);
    gsl_interp_accel_free(acc_plus);
    gsl_spline_free(spline_min);
    gsl_interp_accel_free(acc_min); 
    
    delete[] Ddepth; 
    delete[] Ggam_plus; 
    delete[] Ggam_min;				
}

void Up3IP(double nu) {
    //first
    double D_min   = search_v(Nx, depth, MINIMUM);
    double D_plus  = search_v(Nx, depth, MAXIMUM);
    double D_mid   = (D_plus + D_min) / 2;

    if (strcmp(mid_temp, "default") != 0){
        D_mid = mid[1];
    }	

    int n = 101;
    double dD = (D_plus - D_min) / (n - 1);
    
    double Cp_min, Cp_mid, Cp_plus, Cp_D, Om_min, Om_mid, Om_plus, Om_D, det, kappanu_D;
    double C0_min, C0_mid, C0_plus, C0_D;
    double A, B, C, D, E, F, G, H, J;

    double* Ddepth   = new double[n];
    double* Ggam_min = new double[n];
    double* Ggam_mid = new double[n];
    double* Ggam_plus= new double[n];
		
    for (int i = 0; i < n; i++) {
        Ddepth[i]   = D_min + i * dD;
        kappanu_D   = invers_omega(nu, Ddepth[i]);
        
        C0_min      = sqrt(g * D_min);
        C0_mid      = sqrt(g * D_mid);
        C0_plus     = sqrt(g * D_plus);
        C0_D        = sqrt(g * Ddepth[i]);		
        
        Cp_min      = pow(phase_velocity_i(kappanu_D, D_min), 2);
        Cp_mid      = pow(phase_velocity_i(kappanu_D, D_mid), 2);
        Cp_plus     = pow(phase_velocity_i(kappanu_D, D_plus), 2);
        Cp_D        = pow(phase_velocity_i(kappanu_D, Ddepth[i]), 2);
        
        Om_min      = omega_i(kappanu_D, D_min);
        Om_mid      = omega_i(kappanu_D, D_mid);
        Om_plus     = omega_i(kappanu_D, D_plus);
        Om_D        = omega_i(kappanu_D, Ddepth[i]);
        
        A           = (Om_mid * Cp_plus) - (Om_plus * Cp_mid);
        B           = -(Om_min * Cp_plus) + (Om_plus * Cp_min);
        C           = (Om_min * Cp_mid) - (Om_mid * Cp_min);
        D           = -(C0_mid * Cp_plus) + (C0_plus * Cp_mid);
        E           = (C0_min * Cp_plus) - (C0_plus * Cp_min);

        F           = -(C0_min * Cp_mid) + (C0_mid * Cp_min);
        G           = (C0_mid * Om_plus) - (C0_plus * Om_mid);
        H           = -(C0_min * Om_plus) + (C0_plus * Om_min);
        J           = (C0_min * Om_mid) - (C0_mid * Om_min);
        
        det         = C0_min * A + C0_mid * B + C0_plus * C;
        
        Ggam_min[i] = (A * C0_D + D * Om_D + G * Cp_D ) / det;
        Ggam_mid[i] = (B * C0_D + E * Om_D + H * Cp_D ) / det;
        Ggam_plus[i]= (C * C0_D + F * Om_D + J * Cp_D ) / det;
    }   
	
    gsl_interp_accel *acc_plus = gsl_interp_accel_alloc();
    gsl_spline *spline_plus    = gsl_spline_alloc(gsl_interp_cspline, n);
	    
    gsl_interp_accel *acc_min  = gsl_interp_accel_alloc();
    gsl_spline *spline_min     = gsl_spline_alloc(gsl_interp_cspline, n);  

    gsl_interp_accel *acc_mid  = gsl_interp_accel_alloc ();;
    gsl_spline *spline_mid     = gsl_spline_alloc (gsl_interp_cspline, n);;	
      		  
    gsl_spline_init(spline_plus, Ddepth, Ggam_plus, n);
    gsl_spline_init(spline_min, Ddepth, Ggam_min, n);
    gsl_spline_init(spline_mid, Ddepth, Ggam_mid, n);
	    
    for (int i = 0; i < Nx; i++){
        Csq_min[i] = pow(phase_velocity_i(k[i], D_min), 2);
        Csq_mid[i] = pow(phase_velocity_i(k[i], D_mid), 2);
        Csq_plus[i]= pow(phase_velocity_i(k[i], D_plus), 2);
        
        gam_min[i] = gsl_spline_eval(spline_min, depth[i], acc_min);
        gam_mid[i] = gsl_spline_eval(spline_mid, depth[i], acc_mid);
        gam_plus[i]= gsl_spline_eval(spline_plus, depth[i], acc_plus);		
    }

    gsl_spline_free(spline_min);
    gsl_interp_accel_free(acc_min); 
    gsl_spline_free(spline_mid);
    gsl_interp_accel_free(acc_mid); 
    gsl_spline_free(spline_plus);
    gsl_interp_accel_free(acc_plus); 

    delete[] Ddepth; 
    delete[] Ggam_plus; 
    delete[] Ggam_min;
    delete[] Ggam_mid;			
}

void HSS_runup(complex* uh, complex* C_result, double* Cm, double* Cc, double* Cp, double* gm, double* gc, double* gp) {
    complex* Cminuh  = (complex*) malloc(sizeof(complex) * Nx);
    complex* Cplusuh = (complex*) malloc(sizeof(complex) * Nx);
	    
    for (int i = 0; i < Nx; i++) {//spectral
        Cminuh[i].Re   = Cm[i] * uh[i].Re;
        Cminuh[i].Im   = Cm[i] * uh[i].Im;
        Cplusuh[i].Re  = Cp[i] * uh[i].Re;
        Cplusuh[i].Im  = Cp[i] * uh[i].Im;			
    }
    
    double* ift_uh          = (double*) malloc(sizeof(double) * Nx);
    double* ift_Cminuh      = (double*) malloc(sizeof(double) * Nx);
    double* ift_Cplusuh     = (double*) malloc(sizeof(double) * Nx);
    double* gam_ift_Cminuh  = (double*) malloc(sizeof(double) * Nx);
    double* gam_ift_Cplusuh = (double*) malloc(sizeof(double) * Nx);
    double* gammin_ift_uh   = (double*) malloc(sizeof(double) * Nx);
    double* gamplus_ift_uh  = (double*) malloc(sizeof(double) * Nx);
    
    complex* fft_gammin_uu  = (complex*) malloc(sizeof(complex) * Nx);
    complex* fft_gamplus_uu = (complex*) malloc(sizeof(complex) * Nx);
    complex* HSS1min        = (complex*) malloc(sizeof(complex) * Nx);
    complex* HSS1plus       = (complex*) malloc(sizeof(complex) * Nx);
    complex* HSS2min        = (complex*) malloc(sizeof(complex) * Nx);
    complex* HSS2plus       = (complex*) malloc(sizeof(complex) * Nx);

    ifft_complex_real_nx(uh, ift_uh);
    ifft_complex_real_nx(Cminuh, ift_Cminuh);
    free(Cminuh);
    ifft_complex_real_nx(Cplusuh, ift_Cplusuh); 
    free(Cplusuh);
	
    for (int i = 0; i < Nx; i++) {//spectral
        gammin_ift_uh[i]   = gm[i] * ift_uh[i];
        gamplus_ift_uh[i]  = gp[i] * ift_uh[i];
        gam_ift_Cminuh[i]  = gm[i] * ift_Cminuh[i];
        gam_ift_Cplusuh[i] = gp[i] * ift_Cplusuh[i];	
    }
    
    free(ift_Cminuh);
    free(ift_Cplusuh);
		
    fft_real_complex_nx(gam_ift_Cminuh, HSS1min);
    free(gam_ift_Cminuh);
    fft_real_complex_nx(gam_ift_Cplusuh, HSS1plus);
    free(gam_ift_Cplusuh);		
    fft_real_complex_nx(gammin_ift_uh, fft_gammin_uu);
    free(gammin_ift_uh);
    fft_real_complex_nx(gamplus_ift_uh, fft_gamplus_uu);
    free(gamplus_ift_uh);
	    
    for (int i = 0; i < Nx; i++) {//spectral
        HSS2min[i].Re   = Cm[i] * fft_gammin_uu[i].Re;
        HSS2min[i].Im   = Cm[i] * fft_gammin_uu[i].Im;
        HSS2plus[i].Re  = Cp[i] * fft_gamplus_uu[i].Re;
        HSS2plus[i].Im  = Cp[i] * fft_gamplus_uu[i].Im;	
    }
    
    free(fft_gammin_uu);
    free(fft_gamplus_uu);	
	
    for (int i = 0; i < Nx; i++){	
        C_result[i].Re 	= (HSS1min[i].Re + HSS1plus[i].Re + HSS2min[i].Re + HSS2plus[i].Re) / 2;
        C_result[i].Im 	= (HSS1min[i].Im + HSS1plus[i].Im + HSS2min[i].Im + HSS2plus[i].Im) / 2;           
    }
    
    free(HSS1min);
    free(HSS1plus);
    free(HSS2min);
    free(HSS2plus);
	
    if (mid[0] == 3) {
        complex* Cmiduh  = (complex*) malloc(sizeof(complex) * Nx);
	        
        for (int i = 0; i < Nx; i++) {//spectral
            Cmiduh[i].Re  = Cc[i] * uh[i].Re;
            Cmiduh[i].Im  = Cc[i] * uh[i].Im;		
        }

        double* ift_Cmiduh     = (double*) malloc(sizeof(double) * Nx);
        double* gam_ift_Cmiduh = (double*) malloc(sizeof(double) * Nx);
        double* gammid_ift_uh  = (double*) malloc(sizeof(double) * Nx);
        complex* fft_gammid_uu = (complex*) malloc(sizeof(complex) * Nx);
        complex* HSS1mid       = (complex*) malloc(sizeof(complex) * Nx);
        complex* HSS2mid       = (complex*) malloc(sizeof(complex) * Nx);
        
        ifft_complex_real_nx(Cmiduh, ift_Cmiduh);
        free(Cmiduh);
        
        for (int i = 0; i < Nx; i++) {//spectral
            gammid_ift_uh[i]   = gc[i] * ift_uh[i];
            gam_ift_Cmiduh[i]  = gc[i] * ift_Cmiduh[i];	
        }
        
        free(ift_Cmiduh);
		    
        fft_real_complex_nx(gam_ift_Cmiduh, HSS1mid);
        free(gam_ift_Cmiduh);	
        fft_real_complex_nx(gammid_ift_uh, fft_gammid_uu);
        free(gammid_ift_uh);
          
        for (int i = 0; i < Nx; i++) {//spectral
            HSS2mid[i].Re   = Cc[i] * fft_gammid_uu[i].Re;
            HSS2mid[i].Im   = Cc[i] * fft_gammid_uu[i].Im;
        }
        
        free(fft_gammid_uu);	
        
        //update Csqr_bathy
        for (int i = 0; i < Nx; i++){	
            C_result[i].Re  = C_result[i].Re + (HSS1mid[i].Re + HSS2mid[i].Re) / 2;
            C_result[i].Im  = C_result[i].Im + (HSS1mid[i].Im + HSS2mid[i].Im) / 2;           
        }
        
        free(HSS1mid); 	
        free(HSS2mid);
    }	

    free(ift_uh);	
}

double Cder_i(double k, double d) {
    double cder;
    if (k == 0){
        cder = (g / 2) / sqrt(g * d);
    }
    else {
	cder = g / 2 * sqrt(k / (g * tanh(k * d))) * (1 - pow(tanh(k * d), 2));
    }
    return cder;
}

void Up3IP_runup(double nupeak) {
    double max_k, k_cut, dh;
    double H_min, H_plus, H_mid;

    max_k = search_v(Nx, k, MAXIMUM);
    k_cut = max_k / cutfracwn;
    H_min = 0.0167 * Hs_in;

    if (depth[0] > depth[Nx - 1]){
        H_plus = depth[0] + MaxEtaInit;
    }
    else {
	H_plus = depth[Nx - 1] + MaxEtaInit;
    }
    
    H_mid = H_plus / 6;
    dh    = (H_plus - H_min) / (Ninterp - 1);

    if (strcmp(mid_temp, "default") != 0){
        H_mid = mid[1];
    }

    //allocate vector
    double Cp_min, Cp_mid, Cp_plus, Cp_H, Cder_min, Cder_mid, Cder_plus, Cder_H, det, knu;
    double C0_min, C0_mid, C0_plus, C0_H;
    double A, B, C, D, E, F, G, H, J;
    double *Hj      = new double[Ninterp];
    double *gm_min  = new double[Ninterp];
    double *gm_mid  = new double[Ninterp];
    double *gm_plus = new double[Ninterp];	
	
    for (int j = 0; j < Ninterp; j++){
        Hj[j]	    = H_min + j * dh; 
        knu   	    = invers_omega(nupeak, Hj[j]);
        
        Cp_min      = (phase_velocity_i(knu, H_min));
        Cp_mid      = (phase_velocity_i(knu, H_mid));
        Cp_plus     = (phase_velocity_i(knu, H_plus));
        Cp_H        = (phase_velocity_i(knu, Hj[j]));		
        
        C0_min      = sqrt(g * H_min);
        C0_mid      = sqrt(g * H_mid);
        C0_plus     = sqrt(g * H_plus);
        C0_H        = sqrt(g * Hj[j]);		
        
        Cder_min    = Cder_i(knu, H_min);
        Cder_mid    = Cder_i(knu, H_mid);
        Cder_plus   = Cder_i(knu, H_plus);
        Cder_H      = Cder_i(knu, Hj[j]);
        
        A           = (C0_mid * Cder_plus) - (C0_plus * Cder_mid);
        B           = -(C0_min * Cder_plus) + (C0_plus * Cder_min);
        C           = (C0_min * Cder_mid) - (C0_mid * Cder_min);
        D           = -(Cp_mid * Cder_plus) + (Cp_plus * Cder_mid);
        E           = (Cp_min * Cder_plus) - (Cp_plus * Cder_min);

        F           = -(Cp_min * Cder_mid) + (Cp_mid * Cder_min);
        G           = (Cp_mid * C0_plus) - (Cp_plus * C0_mid);
        H           = -(Cp_min * C0_plus) + (Cp_plus * C0_min);
        J           = (Cp_min * C0_mid) - (Cp_mid * C0_min);

        det         = Cp_min * A + Cp_mid * B + Cp_plus * C;
        
        gm_min[j]  = (A * Cp_H + D * C0_H + G * Cder_H) / det;
        gm_mid[j]  = (B * Cp_H + E * C0_H + H * Cder_H) / det;
        gm_plus[j] = (C * Cp_H + F * C0_H + J * Cder_H) / det;
    }

    //interpolation    		  
    gsl_spline_init(spline_sgp, Hj, gm_plus, Ninterp);
    gsl_spline_init(spline_sgc, Hj, gm_mid, Ninterp);
    gsl_spline_init(spline_sgm, Hj, gm_min, Ninterp);		

    //Operator for HS direct
    Cm 	  = new double[Nx];
    Cp    = new double[Nx];
    Cc    = new double[Nx];
    Cderm = new double[Nx];
    Cderp = new double[Nx];
    Cderc = new double[Nx];
    
    for (int i = 0; i < Nx; i++){
        Cm[i]    = phase_velocity_i(k[i], H_min);
        Cp[i]    = phase_velocity_i(k[i], H_plus);
        Cc[i]    = phase_velocity_i(k[i], H_mid);
        Cderm[i] = Cder_i(k[i], H_min);
        Cderp[i] = Cder_i(k[i], H_plus);
        Cderc[i] = Cder_i(k[i], H_mid);
    }	
    Hmin_interp = H_min;

    delete[] gm_min;
    delete[] gm_mid;
    delete[] gm_plus;
}

void Up3IP_runup_m2(double nupeak) {
    double H_mid, H_plus, H_min, k_cut, wcut, H_minShore;
    double max_k;

    H_plus      = GSL_MAX(depth[0], depth[Nx - 1]) + MaxEtaInit;
    max_k       = search_v(Nx, k, MAXIMUM);
    k_cut       = max_k / cutfracwn;
    H_minShore  = pow((3 * nupeak / 4 / k_cut), 2) / g;
    H_min       = 5 * H_minShore;
    H_mid       = H_plus / 6;
    
    if (strcmp(mid_temp, "default") != 0){
	    H_mid = mid[1];
    }

    int Npoints   = 100;
    double* Htot1 = new double[Npoints];
    double* Htot2 = new double[Npoints];
    double  dh1, dh2, H_min1, H_plus1, H_mid1, H_min2, H_plus2, H_mid2;

    dh1 =  (H_mid - H_min) / (Npoints - 1);
    dh2 =  (H_plus - H_mid) / (Npoints - 1);
    
    for (int j = 0; j < Npoints; j++){
        Htot1[j] = H_min + j * dh1;
        Htot2[j] = H_mid + j * dh2;
    }
    
    H_min1	= H_min;
    H_plus1	= H_mid;
    H_mid1	= (H_min1 + H_plus1) / 2;
    H_min2	= H_mid;
    H_plus2	= H_plus;
    H_mid2	= (H_min2 + H_plus2) / 2;
    k_cut	= k[int(Nx / cutfracwn)];
    wcut	= omega_i(k_cut, H_plus);
	
    double* Cp2_H1 	= new double[Npoints];
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

    double* Cp2_H2 	= new double[Npoints];
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
    double* k_half 	= new double[int(Nx/2)];
    double* gam_min2 	= new double[Npoints];
    double* gam_plus2   = new double[Npoints];
    double* gam_mid2 	= new double[Npoints];

    for (int j = 0; j < int(Nx / 2); j++){
        k_half[j] = k[j];
    }

    double Hj1, Hj2, kappanu_H1, kappanuIG_H1, kappanuL_H1;
    double kappanu_H2, kappanuIG_H2, kappanuL_H2;
    double A1, B1, C1, D1, E1, F1, G1, H1, I1, det1;
    double A2, B2, C2, D2, E2, F2, G2, H2, I2, det2;

    for (int j = 0; j < Npoints; j++){
        Hj1	        = Htot1[j];
        Hj2	        = Htot2[j];    
        kappanu_H1  	= invers_omega(nupeak, Hj1);
        kappanuIG_H1   	= invers_omega(nupeak / 10, Hj1);
        kappanuL_H1 	= invers_omega(wcut, Hj1);
        kappanu_H2  	= invers_omega(nupeak, Hj2);
        kappanuIG_H2   	= invers_omega(nupeak / 10, Hj2);
        kappanuL_H2 	= invers_omega(wcut, Hj2);

        knu1[j]         = kappanu_H1;
        knu2[j]         = kappanu_H2;
        
        Cp2_min1[j]     = pow(phase_velocity_i(kappanu_H1, H_min1), 2);
        Cp2_plus1[j]    = pow(phase_velocity_i(kappanu_H1, H_plus1), 2);
        Cp2_mid1[j]     = pow(phase_velocity_i(kappanu_H1, H_mid1), 2);
        Cp2_H1[j]       = pow(phase_velocity_i(kappanu_H1, Hj1), 2);

        Cp2kL_min1[j]   = pow(phase_velocity_i(kappanuL_H1, H_min1), 2);
        Cp2kL_plus1[j]  = pow(phase_velocity_i(kappanuL_H1, H_plus1), 2);
        Cp2kL_mid1[j]   = pow(phase_velocity_i(kappanuL_H1, H_mid1), 2);
        Cp2kL_H1[j]     = pow(phase_velocity_i(kappanuL_H1, Hj1), 2);

        Cp2IG_min1[j]   = pow(phase_velocity_i(kappanuIG_H1, H_min1), 2);
        Cp2IG_plus1[j]  = pow(phase_velocity_i(kappanuIG_H1, H_plus1), 2);
        Cp2IG_mid1[j]   = pow(phase_velocity_i(kappanuIG_H1, H_mid1), 2);
        Cp2IG_H1[j]     = pow(phase_velocity_i(kappanuIG_H1, Hj1), 2);

        Cp2_min2[j]     = pow(phase_velocity_i(kappanu_H2, H_min2), 2);
        Cp2_plus2[j]    = pow(phase_velocity_i(kappanu_H2, H_plus2), 2);
        Cp2_mid2[j]     = pow(phase_velocity_i(kappanu_H2, H_mid2), 2);
        Cp2_H2[j]       = pow(phase_velocity_i(kappanu_H2, Hj2), 2);

        Cp2kL_min2[j]   = pow(phase_velocity_i(kappanuL_H2, H_min2), 2);
        Cp2kL_plus2[j]  = pow(phase_velocity_i(kappanuL_H2, H_plus2), 2);
        Cp2kL_mid2[j]   = pow(phase_velocity_i(kappanuL_H2, H_mid2), 2);
        Cp2kL_H2[j]     = pow(phase_velocity_i(kappanuL_H2, Hj2), 2);

        Cp2IG_min2[j]   = pow(phase_velocity_i(kappanuIG_H2, H_min2), 2);
        Cp2IG_plus2[j]  = pow(phase_velocity_i(kappanuIG_H2, H_plus2), 2);
        Cp2IG_mid2[j]   = pow(phase_velocity_i(kappanuIG_H2, H_mid2), 2);
        Cp2IG_H2[j]     = pow(phase_velocity_i(kappanuIG_H2, Hj2), 2);

        A1              = (Cp2kL_mid1[j] * Cp2IG_plus1[j]) - (Cp2kL_plus1[j] * Cp2IG_mid1[j]);
        B1              = -(Cp2kL_min1[j] * Cp2IG_plus1[j]) + (Cp2kL_plus1[j] * Cp2IG_min1[j]);
        C1              = (Cp2kL_min1[j] * Cp2IG_mid1[j]) - (Cp2kL_mid1[j] * Cp2IG_min1[j]);

        D1              = -(Cp2_mid1[j] * Cp2IG_plus1[j]) + (Cp2_plus1[j] * Cp2IG_mid1[j]);
        E1              = (Cp2_min1[j] * Cp2IG_plus1[j]) - (Cp2_plus1[j] * Cp2IG_min1[j]);
        F1              = -(Cp2_min1[j] * Cp2IG_mid1[j]) + (Cp2_mid1[j] * Cp2IG_min1[j]);

        G1              = (Cp2_mid1[j] * Cp2kL_plus1[j]) - (Cp2_plus1[j] * Cp2kL_mid1[j]);
        H1              = -(Cp2_min1[j] * Cp2kL_plus1[j]) + (Cp2_plus1[j] * Cp2kL_min1[j]);
        I1              = (Cp2_min1[j] * Cp2kL_mid1[j]) - (Cp2_mid1[j] * Cp2kL_min1[j]);

        det1            = Cp2_min1[j] * A1 + Cp2_mid1[j] * B1 + Cp2_plus1[j] * C1;

        gam_min1[j]     = (A1 * Cp2_H1[j] + D1 * Cp2kL_H1[j] + G1 * Cp2IG_H1[j]) / det1;
        gam_mid1[j]     = (B1 * Cp2_H1[j] + E1 * Cp2kL_H1[j] + H1 * Cp2IG_H1[j]) / det1;
        gam_plus1[j]    = (C1 * Cp2_H1[j] + F1 * Cp2kL_H1[j] + I1 * Cp2IG_H1[j]) / det1;

        A2              = (Cp2kL_mid2[j] * Cp2IG_plus2[j]) - (Cp2kL_plus2[j] * Cp2IG_mid2[j]);
        B2              = -(Cp2kL_min2[j] * Cp2IG_plus2[j]) + (Cp2kL_plus2[j] * Cp2IG_min2[j]);
        C2              = (Cp2kL_min2[j] * Cp2IG_mid2[j]) - (Cp2kL_mid2[j] * Cp2IG_min2[j]);

        D2              = -(Cp2_mid2[j] * Cp2IG_plus2[j]) + (Cp2_plus2[j] * Cp2IG_mid2[j]);
        E2              = (Cp2_min2[j] * Cp2IG_plus2[j]) - (Cp2_plus2[j] * Cp2IG_min2[j]);
        F2              = -(Cp2_min2[j] * Cp2IG_mid2[j]) + (Cp2_mid2[j] * Cp2IG_min2[j]);

        G2              = (Cp2_mid2[j] * Cp2kL_plus2[j]) - (Cp2_plus2[j] * Cp2kL_mid2[j]);
        H2              = -(Cp2_min2[j] * Cp2kL_plus2[j]) + (Cp2_plus2[j] * Cp2kL_min2[j]);
        I2              = (Cp2_min2[j] * Cp2kL_mid2[j]) - (Cp2_mid2[j] * Cp2kL_min2[j]);

        det2            = Cp2_min2[j] * A2 + Cp2_mid2[j] * B2 + Cp2_plus2[j] * C2;

        gam_min2[j]     = (A2 * Cp2_H2[j] + D2 * Cp2kL_H2[j] + G2 * Cp2IG_H2[j]) / det2;
        gam_mid2[j]     = (B2 * Cp2_H2[j] + E2 * Cp2kL_H2[j] + H2 * Cp2IG_H2[j]) / det2;
        gam_plus2[j]    = (C2 * Cp2_H2[j] + F2 * Cp2kL_H2[j] + I2 * Cp2IG_H2[j]) / det2;
    }	
	
	int 	Np          = 2 * Npoints - 1;
	double* HtotComb    = new double[Np];
	double* gam_min1tot = new double[Np];
	double* gam_mid1tot = new double[Np];
	double* gam_plus1tot= new double[Np];
	double* gam_min2tot = new double[Np];
	double* gam_mid2tot = new double[Np];
	double* gam_plus2tot= new double[Np];
	
	for (int j = 0; j < (Npoints); j++){
	    HtotComb[j]     = Htot1[j];
	    gam_min1tot[j]  = gam_min1[j];
	    gam_mid1tot[j]  = gam_mid1[j];
	    gam_plus1tot[j] = gam_plus1[j];	
	}
	
	for (int j = Npoints; j < Np; j++){
	    HtotComb[j]     = Htot2[j + 1 - Npoints];
	    gam_min2tot[j]  = gam_min2[j + 1 - Npoints];
	    gam_mid2tot[j]  = gam_mid2[j + 1 - Npoints];
	    gam_plus2tot[j] = gam_plus2[j + 1 - Npoints];
	}

	//interpolation    		  
	gsl_spline_init (sg_p1, HtotComb, gam_plus1tot, Np);
	gsl_spline_init (sg_c1, HtotComb, gam_mid1tot, Np);
	gsl_spline_init (sg_m1, HtotComb, gam_min1tot, Np);
	gsl_spline_init (sg_p2, HtotComb, gam_plus2tot, Np);
	gsl_spline_init (sg_c2, HtotComb, gam_mid2tot, Np);
	gsl_spline_init (sg_m2, HtotComb, gam_min2tot, Np);

	delete[] Htot1; 
	delete[] Htot2;
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

	//operator
	C2m1 = new double[Nx];
	C2p1 = new double[Nx];
	C2c1 = new double[Nx];
	C2m2 = new double[Nx];
	C2p2 = new double[Nx];
	C2c2 = new double[Nx];
	
	for (int j = 0; j < Nx; j++){
	    C2m1[j] = pow(phase_velocity_i(k[j], H_min1), 2);
	    C2p1[j] = pow(phase_velocity_i(k[j], H_plus1), 2);
	    C2c1[j] = pow(phase_velocity_i(k[j], H_mid1), 2);
	    C2m2[j] = pow(phase_velocity_i(k[j], H_min2), 2);
	    C2p2[j] = pow(phase_velocity_i(k[j], H_plus2), 2);
	    C2c2[j] = pow(phase_velocity_i(k[j], H_mid2), 2);
	}
	
	Hmin_interp = H_min;
	Hplus_interp= H_plus;
}

void define_operator(int Ndim, double nu) {
    if (strcmp(bath, "flat") == 0){
        set_operator(Ndim);
        Csqr = fun_csqr(k, depthinf, Nx);    
        operator_L(L, k, depthinf);
        operator_L0(L0, k);
        operator_M0(M0, k);
        operator_M1(M1, k); 
    }
    else {	
        if (runup == 1){
            Csq_min = new double[Nx];
            Csq_mid = new double[Nx];
            Csq_plus= new double[Nx];
            gam_min = new double[Nx];
            gam_mid = new double[Nx];
            gam_plus= new double[Nx];
            dy	    = new double[Ndim];
            Csqr    = new double[Nx];
            Csqr    = fun_csqr(k, depthinf, Nx); 
            Csqr_v  = new complex[Nx];
            set_operator_bathy(nu, mid);	
        }
	else {
	    dy	 = new double[Ndim];
	    if (strcmp(wavename, "zero") != 0){
	        MaxEtaInit = 1.5 * Hs_in;
	    }
	    Up3IP_runup_m2(nu);
	}
    }
}

void HSS(complex* uh, complex* Csqr_bathy) {
    complex* Cminuh  = (complex*) malloc(sizeof(complex) * Nx);
    complex* Cplusuh = (complex*) malloc(sizeof(complex) * Nx);
	    
    for (int i = 0; i < Nx; i++) {//spectral
        Cminuh[i].Re   = Csq_min[i] * uh[i].Re;
        Cminuh[i].Im   = Csq_min[i] * uh[i].Im;
        Cplusuh[i].Re  = Csq_plus[i] * uh[i].Re;
        Cplusuh[i].Im  = Csq_plus[i] * uh[i].Im;			
    }
    
    double* ift_uh          = (double*) malloc(sizeof(double)*Nx);
    double* ift_Cminuh      = (double*) malloc(sizeof(double)*Nx);
    double* ift_Cplusuh     = (double*) malloc(sizeof(double)*Nx);
    double* gam_ift_Cminuh  = (double*) malloc(sizeof(double)*Nx);
    double* gam_ift_Cplusuh = (double*) malloc(sizeof(double)*Nx);
    double* gammin_ift_uh   = (double*) malloc(sizeof(double)*Nx);
    double* gamplus_ift_uh  = (double*) malloc(sizeof(double)*Nx);
    complex* fft_gammin_uu  = (complex*) malloc(sizeof(complex)*Nx);
    complex* fft_gamplus_uu = (complex*) malloc(sizeof(complex)*Nx);
    complex* HSS1min        = (complex*) malloc(sizeof(complex)*Nx);
    complex* HSS1plus       = (complex*) malloc(sizeof(complex)*Nx);
    complex* HSS2min        = (complex*) malloc(sizeof(complex)*Nx);
    complex* HSS2plus       = (complex*) malloc(sizeof(complex)*Nx);

    ifft_complex_real_nx(uh, ift_uh);
    ifft_complex_real_nx(Cminuh, ift_Cminuh);free(Cminuh);
    ifft_complex_real_nx(Cplusuh, ift_Cplusuh);free(Cplusuh);
	
    for (int i = 0; i < Nx; i++) {
        gammin_ift_uh[i]   = gam_min[i] * ift_uh[i];
        gamplus_ift_uh[i]  = gam_plus[i] * ift_uh[i];
        gam_ift_Cminuh[i]  = gam_min[i] * ift_Cminuh[i];
        gam_ift_Cplusuh[i] = gam_plus[i] * ift_Cplusuh[i];	
    }
    
    free(ift_Cminuh);
    free(ift_Cplusuh);
	    
    fft_real_complex_nx(gam_ift_Cminuh, HSS1min);
    free(gam_ift_Cminuh);
    fft_real_complex_nx(gam_ift_Cplusuh, HSS1plus);
    free(gam_ift_Cplusuh);		
    fft_real_complex_nx(gammin_ift_uh, fft_gammin_uu);
    free(gammin_ift_uh);
    fft_real_complex_nx(gamplus_ift_uh, fft_gamplus_uu);
    free(gamplus_ift_uh);
	    
    for (int i = 0; i < Nx; i++) {//spectral
        HSS2min[i].Re   = Csq_min[i] * fft_gammin_uu[i].Re;
        HSS2min[i].Im   = Csq_min[i] * fft_gammin_uu[i].Im;
        HSS2plus[i].Re  = Csq_plus[i] * fft_gamplus_uu[i].Re;
        HSS2plus[i].Im  = Csq_plus[i] * fft_gamplus_uu[i].Im;	
    }
    
    free(fft_gammin_uu);
    free(fft_gamplus_uu);	
	
    for (int i = 0; i < Nx; i++){	
        Csqr_bathy[i].Re = (HSS1min[i].Re + HSS1plus[i].Re + HSS2min[i].Re + HSS2plus[i].Re) / 2;
        Csqr_bathy[i].Im = (HSS1min[i].Im + HSS1plus[i].Im + HSS2min[i].Im + HSS2plus[i].Im) / 2;           
    }
    
    free(HSS1min);
    free(HSS1plus);
    free(HSS2min);
    free(HSS2plus);
	
    if (mid[0] == 3) {
        complex* Cmiduh  = (complex*) malloc(sizeof(complex) * Nx);
	        
        for (int i = 0; i < Nx; i++) {//spectral
            Cmiduh[i].Re  = Csq_mid[i] * uh[i].Re;
            Cmiduh[i].Im  = Csq_mid[i] * uh[i].Im;		
        }

        double* ift_Cmiduh     = (double*) malloc(sizeof(double) * Nx);
        double* gam_ift_Cmiduh = (double*) malloc(sizeof(double) * Nx);
        double* gammid_ift_uh  = (double*) malloc(sizeof(double) * Nx);
        complex* fft_gammid_uu = (complex*) malloc(sizeof(complex) * Nx);
        complex* HSS1mid       = (complex*) malloc(sizeof(complex) * Nx);
        complex* HSS2mid       = (complex*) malloc(sizeof(complex) * Nx);
	    
        ifft_complex_real_nx(Cmiduh, ift_Cmiduh);
        free(Cmiduh);
        
        for (int i = 0; i < Nx; i++) {//spectral
            gammid_ift_uh[i]   = gam_mid[i] * ift_uh[i];
            gam_ift_Cmiduh[i]  = gam_mid[i] * ift_Cmiduh[i];	
        }
        
        free(ift_Cmiduh);
		    
        fft_real_complex_nx(gam_ift_Cmiduh, HSS1mid);
        free(gam_ift_Cmiduh);	
        fft_real_complex_nx(gammid_ift_uh, fft_gammid_uu);
        free(gammid_ift_uh);
	        
        for (int i = 0; i < Nx; i++) {//spectral
            HSS2mid[i].Re   = Csq_mid[i] * fft_gammid_uu[i].Re;
            HSS2mid[i].Im   = Csq_mid[i] * fft_gammid_uu[i].Im;
        }
        
        free(fft_gammid_uu);	
        
        //update Csqr_bathy
        for (int i = 0; i < Nx; i++){	
            Csqr_bathy[i].Re = Csqr_bathy[i].Re + (HSS1mid[i].Re + HSS2mid[i].Re) / 2;
            Csqr_bathy[i].Im = Csqr_bathy[i].Im + (HSS1mid[i].Im + HSS2mid[i].Im) / 2;           
        }
        
        free(HSS1mid); 	
        free(HSS2mid);
    }	
	
    free(ift_uh);	
}

void Op_M0(fftw_complex* A_hat_fftw, fftw_complex* M0A_hat, int N) {
    complex* A_hat      = new complex[N];
    double*  A          = new double[N];
    complex* HSS_result = new complex[N];

    fftw_to_complex(A_hat_fftw, A_hat,N);
    ifft_complex_real(A_hat, A, N);
    HSS(A_hat,HSS_result);
    
    for (int j = 0; j < Nx; j++){
        M0A_hat[j][0] = -1 / g * prod_1ik(k[j], HSS_result[j].Re, HSS_result[j].Im, REAL);
        M0A_hat[j][1] = -1 / g * prod_1ik(k[j], HSS_result[j].Re, HSS_result[j].Im, IMAG);
    }
    
    delete[] A; 
    delete[] A_hat;
    delete[] HSS_result;
}

void Op_M1(fftw_complex* A_hat_fftw, fftw_complex* M1A_hat, int N) {
    complex* A_hat      = new complex[N];
    double*  A          = new double[N];
    complex* HSS_result = new complex[N];

    fftw_to_complex(A_hat_fftw, A_hat,N);
    ifft_complex_real(A_hat, A, N);
    HSS(A_hat, HSS_result);
    
    for (int j = 0; j < Nx; j++){
        M1A_hat[j][0] = -1 / g * prod_1ik(k[j], HSS_result[j].Re, HSS_result[j].Im, REAL);
        M1A_hat[j][1] = -1 / g * prod_1ik(k[j], HSS_result[j].Re, HSS_result[j].Im, IMAG);
    }
    
    delete[] A; 
    delete[] A_hat;
    delete[] HSS_result;
}

//HS4
void ifft_fft_Hsf4(const double* z, double* L, fftw_complex* LetaM0u_hat, fftw_complex* deletaH4_hat, fftw_complex* deluH4_hat, int N) {
    fftw_complex *etahat_input, *eta;
    fftw_complex* uhat, *u;
    fftw_complex *dxu_hat, *dxu;
    fftw_complex *M0u_hat, *M0u, *dxeta_hat;
    fftw_complex* LetaM0u;
    
    etahat_input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    eta          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    uhat         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    u 	         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxu_hat      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxu          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0u          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0u_hat      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxeta_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);  
    LetaM0u        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++){
        etahat_input[i][0] = z[i];
        etahat_input[i][1] = z[i + N];
        uhat[i][0]         = z[i + 2*N];
        uhat[i][1]         = z[i + 3*N];
        M0u_hat[i][0]	   = prod_1ik(M0[i].Im, z[i + 2*N], z[i + 3*N], REAL);
        M0u_hat[i][1]	   = prod_1ik(M0[i].Im, z[i + 2*N], z[i + 3*N], IMAG);
        dxu_hat[i][0]      = prod_1ik(k[i], z[i + 2*N], z[i + 3*N], REAL);
        dxu_hat[i][1]      = prod_1ik(k[i], z[i + 2*N], z[i + 3*N], IMAG);
        dxeta_hat[i][0]     = prod_1ik(k[i], z[i], z[i + N], REAL);
        dxeta_hat[i][1]     = prod_1ik(k[i], z[i], z[i + N], IMAG);
    }

    fftw_execute_dft(plan_ift_nx, etahat_input, eta);
    fftw_execute_dft(plan_ift_nx, uhat, u);
    fftw_execute_dft(plan_ift_nx, dxu_hat, dxu);
    fftw_execute_dft(plan_ift_nx, M0u_hat, M0u);
    fftw_execute_dft(plan_ift_nx, LetaM0u_hat, LetaM0u);

    fftw_complex *dx2eta_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *dx2eta       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *etaM0dxetaM0u      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *etaM0dxetaM0u_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        //dx2eta_hat      =1i.*k.*1i.*k.*z1
        dx2eta_hat[i][0] = prod_1ik(k[i], dxeta_hat[i][0], dxeta_hat[i][1], REAL);
        dx2eta_hat[i][1] = prod_1ik(k[i], dxeta_hat[i][0], dxeta_hat[i][1], IMAG);
        etaM0dxetaM0u[i][0]  = (eta[i][0] / N) * (LetaM0u[i][0] / N); //ada di Hsf3 juga
        etaM0dxetaM0u[i][1]  = 0;
    }
    
    //dx2eta   =Ifft(dx2eta_hat) & fft(eta.*M0dxetaM0u)
    fftw_execute_dft(plan_ift_nx, dx2eta_hat, dx2eta);
    fftw_execute_dft(plan_ft_nx, etaM0dxetaM0u, etaM0dxetaM0u_hat); 
    
    fftw_complex *dxetaM0dxetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M0dxetaM0dxetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M0dxetaM0dxetaM0u     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        //dxetaM0dxetaM0u_hat  =1i.*k.*fft(eta.*M0dxetaM0u);
        dxetaM0dxetaM0u_hat[i][0] = prod_1ik(k[i], etaM0dxetaM0u_hat[i][0], etaM0dxetaM0u_hat[i][1], REAL);
        dxetaM0dxetaM0u_hat[i][1] = prod_1ik(k[i], etaM0dxetaM0u_hat[i][0], etaM0dxetaM0u_hat[i][1], IMAG);
    } 

    //M0dxetaM0dxetaM0u_hat=M0.*dxetaM0dxetaM0u_hat;
    for (int i = 0; i < N; i++) {
        M0dxetaM0dxetaM0u_hat[i][0] = (M0[i].Im / N) * (dxetaM0dxetaM0u_hat[i][0] / N);
        M0dxetaM0dxetaM0u_hat[i][1] = (M0[i].Im / N) * (dxetaM0dxetaM0u_hat[i][1] / N);

    }
    //Op_M0(dxetaM0dxetaM0u_hat, M0dxetaM0dxetaM0u_hat, N); //problem ifft tidak bisa menggunakan Op_M0

    //M0dxetaM0dxetaM0u    =Ifft(M0dxetaM0dxetaM0u_hat);
    fftw_execute_dft(plan_ift_nx, M0dxetaM0dxetaM0u_hat, M0dxetaM0dxetaM0u);
 
    fftw_complex *eta2dxu     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2dxu_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        eta2dxu[i][0] = pow(eta[i][0] / N, 2) * (dxu[i][0] / N);
        eta2dxu[i][1] = 0;
    }

    //fft(eta.^2.*dxu)
    fftw_execute_dft(plan_ft_nx, eta2dxu, eta2dxu_hat);

    fftw_complex *dxeta2dxu 	    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *dxeta2dxu_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M0dxeta2dxu_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M0dxeta2dxu = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        //dxeta2dxu_hat     =1i.*k.*fft(eta.^2.*dxu)
        dxeta2dxu_hat[i][0] = prod_1ik(k[i], eta2dxu_hat[i][0], eta2dxu_hat[i][1], REAL);
        dxeta2dxu_hat[i][1] = prod_1ik(k[i], eta2dxu_hat[i][0], eta2dxu_hat[i][1], IMAG);
    }

    //M0dxeta2dxu_hat    =M0.*dxeta2dxu_hat;
    //Op_M0(dxeta2dxu_hat, M0dxeta2dxu_hat, N);
    for (int i = 0; i < N; i++) {
        M0dxeta2dxu_hat[i][0] = (M0[i].Im / N) * (dxeta2dxu_hat[i][0] / N);
        M0dxeta2dxu_hat[i][1] = (M0[i].Im / N) * (dxeta2dxu_hat[i][1] / N);
    }

    //M0dxeta2dxu        =Ifft(M0dxeta2dxu_hat);
    fftw_execute_dft(plan_ift_nx, M0dxeta2dxu_hat, M0dxeta2dxu);

    fftw_complex *etaM0dxetaM0dxetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *etaM0dxetaM0dxetaM0u     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M1etaM0dxetaM0dxetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); 
    
    for (int i = 0; i < N; i++) {
        etaM0dxetaM0dxetaM0u[i][0] = (eta[i][0] / N) * (M0dxetaM0dxetaM0u[i][0] / N);
        etaM0dxetaM0dxetaM0u[i][1] = 0;
    }
    
    //etaM0dxetaM0dxetaM0u_hat  =fft(eta.*M0dxetaM0dxetaM0u);
    fftw_execute_dft(plan_ft_nx, etaM0dxetaM0dxetaM0u, etaM0dxetaM0dxetaM0u_hat);

    //M1etaM0dxetaM0dxetaM0u_hat=M1.*etaM0dxetaM0dxetaM0u_hat;
    //Op_M1(etaM0dxetaM0dxetaM0u_hat, M1etaM0dxetaM0dxetaM0u_hat, N);
    for (int i = 0; i < N; i++) {
        M1etaM0dxetaM0dxetaM0u_hat[i][0] = (M1[i].Im / N) * (etaM0dxetaM0dxetaM0u_hat[i][0] / N);
        M1etaM0dxetaM0dxetaM0u_hat[i][1] = (M1[i].Im / N) * (etaM0dxetaM0dxetaM0u_hat[i][1] / N);
    }
    
    fftw_complex *etaM0dxeta2dxu_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *etaM0dxeta2dxu       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M1etaM0dxeta2dxu_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        etaM0dxeta2dxu[i][0] = (eta[i][0] / N) * (M0dxeta2dxu[i][0] / N);
        etaM0dxeta2dxu[i][1] = 0;
    }
    
    //etaM0dxeta2dxu_hat  =fft(eta.*M0dxeta2dxu);
    fftw_execute_dft(plan_ft_nx, etaM0dxeta2dxu, etaM0dxeta2dxu_hat);

    //M1etaM0dxeta2dxu_hat=M1.*etaM0dxeta2dxu_hat;
    //Op_M1(etaM0dxeta2dxu_hat, M1etaM0dxeta2dxu_hat, N);
    for (int i = 0; i < N; i++) {
        M1etaM0dxeta2dxu_hat[i][0] = (M1[i].Im / N) * (etaM0dxeta2dxu_hat[i][0] / N);
        M1etaM0dxeta2dxu_hat[i][1] = (M1[i].Im / N) * (etaM0dxeta2dxu_hat[i][1] / N);
    }

    fftw_complex *dxM0u_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *dxM0u       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        dxM0u_hat[i][0] = prod_1ik(k[i], M0u_hat[i][0], M0u_hat[i][1], REAL);
        dxM0u_hat[i][1] = prod_1ik(k[i], M0u_hat[i][0], M0u_hat[i][1], IMAG);
    }

    //Ifft(1i.*k.*M0u_hat))
    fftw_execute_dft(plan_ift_nx, dxM0u_hat, dxM0u);

    fftw_complex *eta3dxM0u = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta3dxM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        eta3dxM0u[i][0] = pow(eta[i][0] / N, 3) * (dxM0u[i][0] / N);
        eta3dxM0u[i][1] = 0;
    }
    
    //fft(eta.^3.*Ifft(1i.*k.*M0u_hat))
    fftw_execute_dft(plan_ft_nx, eta3dxM0u, eta3dxM0u_hat);

    fftw_complex *dxeta3dxM0u_hat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M1dxeta3dxM0u_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        //dxeta3dxM0u_hat     =1i.*k.*fft(eta.^3.*Ifft(1i.*k.*M0u_hat));
        dxeta3dxM0u_hat[i][0] = prod_1ik(k[i], eta3dxM0u_hat[i][0], eta3dxM0u_hat[i][1], REAL);
        dxeta3dxM0u_hat[i][1] = prod_1ik(k[i], eta3dxM0u_hat[i][0], eta3dxM0u_hat[i][1], IMAG);
    }

    //M1dxeta3dxM0u_hat   =M1.*dxeta3dxM0u_hat;
    //Op_M1(dxeta3dxM0u_hat, M1dxeta3dxM0u_hat, N);
    for (int i = 0; i < N; i++) {
        M1dxeta3dxM0u_hat[i][0] = (M1[i].Im / N) * (dxeta3dxM0u_hat[i][0] / N);
        M1dxeta3dxM0u_hat[i][1] = (M1[i].Im / N) * (dxeta3dxM0u_hat[i][1] / N);
    }

    fftw_complex *eta2dx2etaM0u_hat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2dx2etaM0u         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *M1eta2dx2etaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        eta2dx2etaM0u[i][0] = pow(eta[i][0] / N, 2) * (dx2eta[i][0] / N) * (M0u[i][0] / N);
        eta2dx2etaM0u[i][1] = 0;
    }

    //eta2dx2etaM0u_hat   =fft(eta.^2.*dx2eta.*M0u);
    fftw_execute_dft(plan_ft_nx, eta2dx2etaM0u, eta2dx2etaM0u_hat);

    //M1eta2dx2etaM0u_hat   =M1.*eta2dx2etaM0u_hat;
    //Op_M1(eta2dx2etaM0u_hat, M1eta2dx2etaM0u_hat, N);
    for (int i = 0; i < N; i++) {
        M1eta2dx2etaM0u_hat[i][0] = (M1[i].Im / N) * (eta2dx2etaM0u_hat[i][0] / N);
        M1eta2dx2etaM0u_hat[i][1] = (M1[i].Im / N) * (eta2dx2etaM0u_hat[i][1] / N);
    }

    //deleta
    fftw_complex *eta2dxu2_hat      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2dxu2          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2M0u2          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2M0u2_hat      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *dxeta2M0u2_hat    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        eta2dxu2[i][0]          = -(pow(eta[i][0] / N, 2) * pow(dxu[i][0] / N, 2));
        eta2dxu2[i][1]          = 0;
        eta2M0u2[i][0]          = (pow(eta[i][0] / N, 2)) * (pow(M0u[i][0] / N, 2));
        eta2M0u2[i][1]          = 0;
        dxeta2M0u2_hat[i][0]    = prod_1ik(k[i], eta2M0u2_hat[i][0], eta2M0u2_hat[i][1], REAL);
        dxeta2M0u2_hat[i][1]    = prod_1ik(k[i], eta2M0u2_hat[i][0], eta2M0u2_hat[i][1], IMAG);
    }

    fftw_execute_dft(plan_ft_nx, eta2dxu2, eta2dxu2_hat);
    fftw_execute_dft(plan_ft_nx, eta2M0u2, eta2M0u2_hat);

    fftw_complex *dx2eta2M0u2_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        dx2eta2M0u2_hat[i][0] = prod_1ik(k[i], dxeta2M0u2_hat[i][0], dxeta2M0u2_hat[i][1], REAL);
        dx2eta2M0u2_hat[i][1] = prod_1ik(k[i], dxeta2M0u2_hat[i][0], dxeta2M0u2_hat[i][1], IMAG);
    }
    
    for (int i = 0; i < N; i++) {
        deletaH4_hat[i][0] = (1/2.) * eta2dxu2_hat[i][0] - pow(LetaM0u[i][0], 2) - (2 * (M0u[i][0] * M0dxetaM0dxetaM0u[i][0])) - (2 * (eta[i][0] * dxu[i][0] * LetaM0u[i][0])) - (M0u[i][0] * M0dxeta2dxu[i][0]) + (pow(eta[i][0], 2) * pow(dxM0u[i][0], 2)) - (eta[i][0] * dx2eta[i][0] * pow(M0u[i][0], 2)) - ((1/4.) * dx2eta2M0u2_hat[i][0]);
        deletaH4_hat[i][1] = 0;
    }
    
    //deluh
    fftw_complex *eta3dxu          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta3dxu_hat      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2M0dxetaM0u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *eta2M0dxetaM0u_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        eta3dxu[i][0] = (pow(eta[i][0] / N, 3)) * (dxu[i][0] / N);
        eta3dxu[i][1] = (pow(eta[i][1] / N, 3)) * (dxu[i][1] / N);
        eta2M0dxetaM0u[i][0] = (pow(eta[i][0] / N, 2)) * (LetaM0u[i][0] / N);
        eta2M0dxetaM0u[i][1] = 0;
    }

    fftw_execute_dft(plan_ft_nx, eta3dxu, eta3dxu_hat);
    fftw_execute_dft(plan_ft_nx, eta2M0dxetaM0u, eta2M0dxetaM0u_hat);

    fftw_complex *dxeta3dxu_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *dxeta2M0dxetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    for (int i = 0; i < N; i++) {
        dxeta3dxu_hat[i][0] = prod_1ik(k[i], eta3dxu_hat[i][0], eta3dxu_hat[i][1], REAL);
        dxeta3dxu_hat[i][1] = prod_1ik(k[i], eta3dxu_hat[i][0], eta3dxu_hat[i][1], IMAG);
        dxeta2M0dxetaM0u_hat[i][0] = prod_1ik(k[i], eta2M0dxetaM0u_hat[i][0], eta2M0dxetaM0u_hat[i][0], REAL);
        dxeta2M0dxetaM0u_hat[i][1] = prod_1ik(k[i], eta2M0dxetaM0u_hat[i][0], eta2M0dxetaM0u_hat[i][1], IMAG);
    }

    for (int i = 0; i < N; i++) {
        deluH4_hat[i][0] = -((- (1/3.) * dxeta3dxu_hat[i][0]) 
                              - (M1etaM0dxetaM0dxetaM0u_hat[i][0]) 
                              - ((1/2.) * M1etaM0dxeta2dxu_hat[i][0]) 
                              - ((1/2.) * dxeta2M0dxetaM0u_hat[i][0])
                              - ((1/3.) * M1dxeta3dxM0u_hat[i][0])
                              - ((1/2.) * M1eta2dx2etaM0u_hat[i][0]));
        deluH4_hat[i][1] = -((- (1/3.) * dxeta3dxu_hat[i][1]) 
                              - (M1etaM0dxetaM0dxetaM0u_hat[i][1]) 
                              - ((1/2.) * M1etaM0dxeta2dxu_hat[i][1]) 
                              - ((1/2.) * dxeta2M0dxetaM0u_hat[i][1])
                              - ((1/3.) * M1dxeta3dxM0u_hat[i][1])
                              - ((1/2.) * M1eta2dx2etaM0u_hat[i][1]));
    }
    
    fftw_free(etahat_input); fftw_free(eta);
    fftw_free(uhat); fftw_free(u);
    fftw_free(dxu_hat); fftw_free(dxu);
    fftw_free(M0u_hat); fftw_free(M0u);
    fftw_free(dxeta_hat);
    fftw_free(LetaM0u);

    fftw_free(dx2eta_hat); fftw_free(dx2eta);
    fftw_free(etaM0dxetaM0u); fftw_free(etaM0dxetaM0u_hat);

    fftw_free(dxetaM0dxetaM0u_hat); 
    fftw_free(M0dxetaM0dxetaM0u_hat); fftw_free(M0dxetaM0dxetaM0u);

    fftw_free(eta2dxu); fftw_free(eta2dxu_hat);

    fftw_free(dxeta2dxu); fftw_free(dxeta2dxu_hat);
    fftw_free(M0dxeta2dxu_hat);

    fftw_free(etaM0dxetaM0dxetaM0u_hat); fftw_free(etaM0dxetaM0dxetaM0u);
    fftw_free(M1etaM0dxetaM0dxetaM0u_hat);

    fftw_free(etaM0dxeta2dxu_hat); fftw_free(etaM0dxeta2dxu);
    fftw_free(M1etaM0dxeta2dxu_hat); fftw_free(M0dxeta2dxu);

    fftw_free(eta3dxM0u);  fftw_free(eta3dxM0u_hat);
    fftw_free(dxeta3dxM0u_hat); fftw_free(M1dxeta3dxM0u_hat);

    fftw_free(eta2dx2etaM0u_hat); fftw_free(eta2dx2etaM0u);
    fftw_free(M1eta2dx2etaM0u_hat);

    fftw_free(eta2dxu2_hat); fftw_free(eta2dxu2);
    fftw_free(dxM0u_hat); fftw_free(dxM0u);
    fftw_free(eta2M0u2); fftw_free(eta2M0u2_hat);
    fftw_free(dxeta2M0u2_hat);

    fftw_free(dx2eta2M0u2_hat);
    fftw_free(eta3dxu); fftw_free(eta3dxu_hat);
    fftw_free(eta2M0dxetaM0u); fftw_free(eta2M0dxetaM0u_hat);

    fftw_free(dxeta3dxu_hat); fftw_free(dxeta2M0dxetaM0u_hat);
}

void ifft_fft_Hsf2_v(const double* z, double* L, fftw_complex* eta_u_hat, fftw_complex* LetaM0u_hat, fftw_complex* u2_M0u2_hat, fftw_complex* damping_eta, fftw_complex* damping_u, double* dampchar, int N) {
    //first loop
    fftw_complex* etahat_input, *eta;
    fftw_complex* uhat, *u;
    fftw_complex* M0uhat, *M0u;

    // Allocate input & output array 
    etahat_input   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    eta            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    uhat           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    u              = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0uhat         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0u            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
    for(int i = 0; i < N; i++){//spectral
        etahat_input[i][0] = z[i];
        etahat_input[i][1] = z[i + N];
        uhat[i][0]  	   = z[i + 2*N];
        uhat[i][1]    	   = z[i + 3*N];
    }
    
    Op_M0(uhat, M0uhat, N);

    //IFFT
    fftw_execute_dft(plan_ift_nx, etahat_input, eta);
    fftw_execute_dft(plan_ift_nx, uhat, u);
    fftw_execute_dft(plan_ift_nx, M0uhat, M0u);

    //second loop  
    fftw_complex  *eta_u      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta_M0u    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *etaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *u2_M0u2    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *damp_eta   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *damp_u     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++){
        eta_u[i][0]    = (eta[i][0] / N) * (u[i][0] / N);
        eta_u[i][1]    = 0;
        eta_M0u[i][0]  = (eta[i][0] / N) * (M0u[i][0] / N);
        eta_M0u[i][1]  = 0;
        u2_M0u2[i][0]  = (pow(u[i][0] / N, 2) - pow(M0u[i][0] / N, 2));
        u2_M0u2[i][1]  = 0;  
        damp_eta[i][0] = dampcoef[i] * (eta[i][0] / N) * dampchar[i];
        damp_eta[i][1] = 0;
        damp_u[i][0]   = dampcoef[i] * (u[i][0] / N) * dampchar[i];
        damp_u[i][1]   = 0;
    }

    // Forward DFT
    fftw_execute_dft(plan_ft_nx, eta_u, eta_u_hat);
    fftw_execute_dft(plan_ft_nx, eta_M0u, etaM0u_hat);
    fftw_execute_dft(plan_ft_nx, u2_M0u2, u2_M0u2_hat);
    fftw_execute_dft(plan_ft_nx, damp_eta, damping_eta);
    fftw_execute_dft(plan_ft_nx, damp_u, damping_u);

    //third loop
    fftw_complex* M1etaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    Op_M1(etaM0u_hat, M1etaM0u_hat, N);
	    
    for (int i = 0; i < N; i++){//spectral
        LetaM0u_hat[i][0] = prod_1ik(k[i], M1etaM0u_hat[i][0], M1etaM0u_hat[i][1], REAL);
        LetaM0u_hat[i][1] = prod_1ik(k[i], M1etaM0u_hat[i][0], M1etaM0u_hat[i][1], IMAG);
    }
	
    fftw_free(M1etaM0u_hat);   
    fftw_free(etahat_input);
    fftw_free(eta);
    fftw_free(uhat);
    fftw_free(u);
    fftw_free(M0uhat);
    fftw_free(M0u);
    fftw_free(eta_M0u); 
    fftw_free(etaM0u_hat);
    fftw_free(eta_u);
    fftw_free(u2_M0u2); 
    fftw_free(damp_eta);
    fftw_free(damp_u); 
}

void define_delH1(const double* z, complex* deluH1, complex* deletaH1) {
    if (strcmp(bath,"flat") == 0){
        for (int i = 0; i < Nx; i++){		  
            deluH1[i].Re    = (Csqr[i] / g) * k[i] * z[i + 3*Nx];
            deluH1[i].Im    = -(Csqr[i] / g) * k[i] * z[i + 2*Nx];
            deletaH1[i].Re  = g * k[i] * z[i + Nx];
            deletaH1[i].Im  = -g * k[i] * z[i];           
        }
    }
    else {		
        complex* u_hat 	 = (complex*) malloc(sizeof(complex) * Nx);		
        for (int i = 0; i < Nx; i++) {//spectral
            u_hat[i].Re  = z[i + 2*Nx];
            u_hat[i].Im  = z[i + 3*Nx];
        }		
        HSS(u_hat, Csqr_v);
        free(u_hat);		
        for (int i = 0; i < Nx; i++){	
            deluH1[i].Re    = (Csqr_v[i].Im / g) * k[i];
            deluH1[i].Im    = -(Csqr_v[i].Re / g) * k[i];
            deletaH1[i].Re  = g * k[i] * z[i + Nx];
            deletaH1[i].Im  = -g * k[i] * z[i];           
        }			
    }	
}

void ifft_fft_Hsf3_v(const double* z, double* L, fftw_complex* LetaM0u_hat, fftw_complex* deletaH3_hat, fftw_complex* deluH3_hat, int N) {	
    fftw_complex* etahat_input, *eta, *uhat;
    fftw_plan     plan_iftetahat;
    fftw_complex* dxu_hat, *dxu;
    fftw_plan     plan_iftdxu_hat;
    fftw_complex* M0u_hat, *M0u;
    fftw_plan     plan_iftM0u_hat;
    fftw_complex* LetaM0u;
    fftw_plan 	  plan_LetaM0u;
      
    // Allocate input & output array 
    etahat_input   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    eta            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    uhat           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxu_hat        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    dxu            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);	
    M0u_hat        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M0u            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);	
    LetaM0u        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    plan_iftetahat = fftw_plan_dft_1d(N, etahat_input, eta, FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_iftdxu_hat= fftw_plan_dft_1d(N, dxu_hat, dxu, FFTW_BACKWARD, FFTW_ESTIMATE);	
    plan_iftM0u_hat= fftw_plan_dft_1d(N, M0u_hat, M0u, FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_LetaM0u   = fftw_plan_dft_1d(N, LetaM0u_hat, LetaM0u, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(int i = 0; i < N; i++){//spectral
        etahat_input[i][0] = z[i];
        etahat_input[i][1] = z[i + N];
        uhat[i][0]         = z[i + 2*N];
        uhat[i][1]         = z[i + 3*N];
        dxu_hat[i][0]      = prod_1ik(k[i], z[i + 2*N], z[i + 3*N], REAL);
        dxu_hat[i][1]      = prod_1ik(k[i],z[i + 2*N], z[i + 3*N], IMAG);
    }
    Op_M0(uhat, M0u_hat, N);
      
    fftw_execute(plan_iftetahat);
    fftw_execute(plan_iftdxu_hat);
    fftw_execute(plan_iftM0u_hat);
    fftw_execute(plan_LetaM0u);

    fftw_destroy_plan(plan_iftetahat);
    fftw_destroy_plan(plan_iftdxu_hat);
    fftw_destroy_plan(plan_iftM0u_hat);	
    fftw_destroy_plan(plan_LetaM0u);	

    fftw_complex  *etaM0dxetaM0u      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *etaM0dxetaM0u_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    fftw_plan     plan_etaM0dxetaM0u_hat;
    fftw_complex  *eta2_dxu           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta2_dxu_hat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    fftw_plan     plan_eta2dxu_hat;
    fftw_complex  *eta2_M0u           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex  *eta2_M0u_hat       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    fftw_plan     plan_eta2M0u_hat;
    fftw_complex  *deletaH3           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan_deletaH3_hat       = fftw_plan_dft_1d(N, deletaH3, deletaH3_hat, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create plans 
    plan_etaM0dxetaM0u_hat            = fftw_plan_dft_1d(N, etaM0dxetaM0u, etaM0dxetaM0u_hat, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_eta2dxu_hat    	      = fftw_plan_dft_1d(N, eta2_dxu, eta2_dxu_hat,  FFTW_FORWARD, FFTW_ESTIMATE);
    plan_eta2M0u_hat                  = fftw_plan_dft_1d(N, eta2_M0u, eta2_M0u_hat,  FFTW_FORWARD, FFTW_ESTIMATE);
  
    // Populate input data 
    for (int i = 0; i < N; i++){//real
        etaM0dxetaM0u[i][0]  = (eta[i][0] / N) * (LetaM0u[i][0] / N);
        etaM0dxetaM0u[i][1]  = 0;
        eta2_dxu[i][0]       = pow(eta[i][0] / N, 2) * (dxu[i][0] / N);
        eta2_dxu[i][1]       = 0;
        eta2_M0u[i][0]       = pow(eta[i][0] / N, 2) * (M0u[i][0] / N);
        eta2_M0u[i][1]       = 0;
        deletaH3[i][0]       = (M0u[i][0] / N) * ((eta[i][0] / N) * (dxu[i][0] / N) + (LetaM0u[i][0] / N)) ;
        deletaH3[i][1]       = 0 ;
    }

    // Forward DFT 
    fftw_execute(plan_etaM0dxetaM0u_hat);  
    fftw_execute(plan_eta2dxu_hat);  
    fftw_execute(plan_eta2M0u_hat); 
    fftw_execute(plan_deletaH3_hat); 

    // Free memory 
    fftw_destroy_plan(plan_etaM0dxetaM0u_hat);  
    fftw_destroy_plan(plan_eta2dxu_hat);  
    fftw_destroy_plan(plan_eta2M0u_hat); 
    fftw_destroy_plan(plan_deletaH3_hat); 

    //L is changed to be ik Op_M0
    fftw_complex  *M1etaM0dxetaM0u_hat, *M1eta2_dxu_hat;
    M1etaM0dxetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    M1eta2_dxu_hat      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    Op_M0(etaM0dxetaM0u_hat, M1etaM0dxetaM0u_hat, N);
    Op_M0(eta2_dxu_hat, M1eta2_dxu_hat, N);
	
    complex temp;
    for(int i = 0; i < N; i++){//spectral
        temp.Re = M1etaM0dxetaM0u_hat[i][0] + (1/2.) * M1eta2_dxu_hat[i][0];
        temp.Im = M1etaM0dxetaM0u_hat[i][1] + (1/2.) * M1eta2_dxu_hat[i][1];
        deluH3_hat[i][0] = prod_1ik(k[i], temp.Re, temp.Im, REAL) - (1/2.) * pow(k[i], 2) * eta2_M0u_hat[i][0];
        deluH3_hat[i][1] = prod_1ik(k[i], temp.Re, temp.Im, IMAG) - (1/2.) * pow(k[i], 2) * eta2_M0u_hat[i][1];
    }

    fftw_free(M1eta2_dxu_hat);fftw_free(M1etaM0dxetaM0u_hat);
    fftw_free(etahat_input);fftw_free(eta);fftw_free(uhat);
    fftw_free(dxu);fftw_free(dxu_hat);
    fftw_free(M0u);fftw_free(M0u_hat);	 
    fftw_free(etaM0dxetaM0u_hat);fftw_free(etaM0dxetaM0u);
    fftw_free(eta2_dxu); fftw_free(eta2_dxu_hat);
    fftw_free(eta2_M0u); fftw_free(eta2_M0u_hat);
    fftw_free(deletaH3);fftw_free(LetaM0u);
}

void read_kinematics() {
    if (strcmp(cut_temp, "default") != 0){
        cutfracwn = atoi(cut_temp);
    }

    if (strcmp(mid_temp, "default") != 0){
        mid[1] = atof(mid_temp);
    }

    if (strcmp(min_z, "default") == 0){
        zsurf1 = -2.5 * Hs_in;
        if (strcmp(initial, "user_defined") == 0){
            zsurf1 = -0.6 * Hs_in;
        }
    }
    else {
        zsurf1 = atof(min_z);
    }

    if (strcmp(max_z,"default")==0){
        zsurf2 = 2.5 * Hs_in;
        if (strcmp(initial, "user_defined") == 0){
            zsurf2 = 0.6 * Hs_in;
        }
    }
    else {
        zsurf2 = atof(max_z);
    }	

    //set partition
    NtoutIP = n / npartition;
    Ntrest  = n - (npartition) * NtoutIP;
    Nheader = NtoutIP * NxoutIP * NzoutIP;
    toutIP  = new double* [npartition];
    trest   = new double [Ntrest];
    for (int j = 0; j < npartition; j++){
        toutIP[j] = new double[NtoutIP];
        for (int tt = 0; tt < NtoutIP; tt++){
            toutIP[j][tt] = ttrace[tt + (j * NtoutIP)];
        }
    }
    for (int j = 0; j < Ntrest; j++){
        trest[j] = ttrace[npartition * NtoutIP + j];
    }

    //define max-min depth
    maxD 	= search_v(Nx, depth, MAXIMUM);
    minD 	= search_v(Nx, depth, MINIMUM);

    //only for kinematics
    if (strcmp(kinematic, "yes") == 0){
        //define xout
        indxIP1 	= closest(x, xinterv1, Nx);	
        indxIP2 	= closest(x, xinterv2, Nx);
        NxoutIP	= ((indxIP2 - indxIP1)) + 1;	
        
        xoutIP 	= new double[NxoutIP];		
        for (int i = 0; i < NxoutIP; i++){
            xoutIP[i] = x[indxIP1 + i];
        }
        
        //define zout, z at surface layer must be lower than maximum eta
        maxEta = zsurf2; //user
        minEta = zsurf1;
        
        double* zsurf  = (double*) malloc (sizeof(double) * nsurf);
        double  dzsurf = (maxEta - minEta) / (nsurf - 1);
        double  kpeak  = invers_omega(nu_peak, maxD);
        
        for (int j = 0; j < nsurf; j++){
            zsurf[j] = minEta + j * dzsurf;
        }
        
        //compute z at deeper layer using Riemann upper sum 
        double  Flux, tol;

        Flux = 1 / (kpeak * sinh(kpeak * maxD)) * (cosh(kpeak * (maxD + minEta)) - 1);
        tol  = Flux / (ndeep);	
        
        double* zdeep = (double*) malloc (sizeof(double) * ndeep);
        zdeep[0]  = -1 * (minEta - dzsurf);
        
        for (int i = 1; i < ndeep; i++) {			
            zdeep[i]  = zdeep[i - 1] + (tol / (sinh(kpeak * (maxD - zdeep[i - 1])) / sinh(kpeak * maxD)));
        }
        
        //combining surface layer and deeper layer
        NzoutIP = nsurf + ndeep;
        zoutIP  = new double[NzoutIP];	 
        for (int i = 0; i < nsurf; i++) {
            zoutIP[i] = zsurf[i];
        }
        for (int i = nsurf; i < NzoutIP; i++) {
            zoutIP[i] = -1 * zdeep[i - nsurf];
        }
        sorting(zoutIP, NzoutIP);
        
        free(zsurf);free(zdeep);

        //initiate velocity
        Velx = new double*[NxoutIP];
        Velz = new double*[NxoutIP];
        for (int j = 0; j < NxoutIP; j++){
            Velx[j] = new double[NzoutIP];
            Velz[j] = new double[NzoutIP];
        }

        for (int i = 0; i < NxoutIP; i++){
            for (int j = 0; j < NzoutIP; j++){
                Velx[i][j] = 0;
                Velz[i][j] = 0;
            }
        }
    }
}

int primeFactors(int number) {  
    // Print the number of 2s that divide n  
    while (number % 2 == 0)  
    {   
        number = number / 2;  
    }  
  
    // n must be odd at this point. So we can skip  
    // one element (Note i = i +2)  
    for (int i = 3; i <= sqrt(number); i += 2)  
    {  
        // While i divides n, print i and divide n  
        while (number % i == 0)  
        {    
            number = number / i;  
        }  
    }  
  
	if (number > 13) return 1;
	else return 0;
} 

int indexOfFirstZero(double* arr, int n) { 
    // traverse the array from left to right 
    for (int i = 0; i < n; i++) {  
        // if true, then return i 
        if (arr[i] == 0) return i; 
    }
    // 1's are not present in the array 
    return -1; 
} 

int indexOfLastZero(double* arr, int n) { 
    // traverse the array from right to left
    for (int i = n - 1; i > -1; i--) {  
        // if true, then return i 
        if (arr[i] == 0) return i;
    }  
    // 1's are not present in the array 
    return -1; 
}

void surf_velocity_z(double t, FILE* datphi_x_gauge, FILE* datphi_z_gauge, complex* phi_hatf, double* etaf, double* k, complex* eta_hatf) {
    double* dxPhi_eta = new double[Nx];
    double* dzPhi_eta = new double[Nx];

    complex* gammahat = new complex[Nx];
    double*  gamma    = new double[Nx];
    for (int j = 0; j < Nx; j++){
        gammahat[j].Re = -phi_hatf[j].Re * k[j] * tanh(k[j] * (etaf[j] + depth[j]));	
        gammahat[j].Im = -phi_hatf[j].Im * k[j] * tanh(k[j] * (etaf[j] + depth[j]));
    }
    ifft_complex_real_nx(gammahat, gamma);

    for (int j = 0; j < Nx; j++){
	dzPhi_eta[j] = -gamma[j];	
    }

    complex* part1_hat = new complex[Nx];
    complex* part2_hat = new complex[Nx];
    double* part1 = new double[Nx];
    double* part2 = new double[Nx];
    for (int j = 0; j < Nx; j++){
        part1_hat[j].Re = prod_1ik(k[j], phi_hatf[j].Re, phi_hatf[j].Im, REAL);
        part1_hat[j].Im = prod_1ik(k[j], phi_hatf[j].Re, phi_hatf[j].Im, IMAG);
        part2_hat[j].Re = prod_1ik(-k[j], eta_hatf[j].Re, eta_hatf[j].Im, REAL);
        part2_hat[j].Im = prod_1ik(-k[j], eta_hatf[j].Re, eta_hatf[j].Im, IMAG);
    }
    ifft_complex_real_nx(part1_hat, part1);
    ifft_complex_real_nx(part2_hat, part2);

    for (int j = 0; j < Nx; j++){
	dxPhi_eta[j] = part1[j] + gamma[j] * part2[j];
    }

    delete[] gammahat;delete[] gamma;	
    delete[] part1_hat;delete[] part1;
    delete[] part2_hat;delete[] part2;

    //writing data
    double dxPhig,dzPhig;
    gsl_interp_accel *acc1  = gsl_interp_accel_alloc();
    gsl_spline *spline1     = gsl_spline_alloc(gsl_interp_cspline, Nx);				
    gsl_spline *spline2     = gsl_spline_alloc(gsl_interp_cspline, Nx);
					      
    gsl_spline_init(spline1, x, dxPhi_eta, Nx);
    gsl_spline_init(spline2, x, dzPhi_eta, Nx);
		    
    fprintf(datphi_x_gauge, "%f ", t);
    fprintf(datphi_z_gauge, "%f ", t);
    for (int j = 0; j < ngauge; j++) {
        dxPhig 	= gsl_spline_eval(spline1, gauge[j], acc1);
        dzPhig 	= gsl_spline_eval(spline2, gauge[j], acc1);
        fprintf(datphi_x_gauge, "%f ", dxPhig);
        fprintf(datphi_z_gauge, "%f ", dzPhig);
    }
    fprintf(datphi_x_gauge, "\n");
    fprintf(datphi_z_gauge, "\n");	
			    
    gsl_spline_free(spline1);
    gsl_spline_free(spline2);
    gsl_interp_accel_free(acc1);
    delete[] dxPhi_eta;
    delete[] dzPhi_eta;
}

void update_damp(int index, int n, double* dampchar, double* Chi_couple) {
    for (int j = 0; j < n; j++){
        if (j < index){
            dampchar[j] = 10 * (1 - Chi_couple[j]);
        }
    }
}
