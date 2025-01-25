int rhs_HsF1(double t, const double* z, double* dz, void *params)
{
    complex *deluH1,*deletaH1;
    fftw_complex *damping_eta,*damping_u;	
    damping_eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx); 
    deluH1  	= (complex*) malloc(sizeof(complex) * Nx);
    deletaH1  	= (complex*) malloc(sizeof(complex) * Nx);
    complex* couple_etahat = new complex[Nx];
    complex* couple_uhat = new complex[Nx];

    #pragma omp parallel for
    for(int j = 0; j < Nx; j++){
        couple_etahat[j].Re = 0;
        couple_uhat[j].Re   = 0;
        couple_etahat[j].Im = 0;
        couple_uhat[j].Im   = 0;
    }

    ifft_fft_Hsf1(z, L, damping_eta, damping_u, dampchar, Nx);
    define_delH1(z, deluH1, deletaH1);//flat or varying

    //IF wind and/or friction
    if ((strcmp(friction, "yes") == 0) || (strcmp(wind, "yes") == 0) || (strcmp(wall, "yes") == 0) || (strcmp(coupling, "yes") == 0)) {
        double  *eta, *u;
        complex *etahat, *uhat;
        eta 	= (double*) malloc(sizeof(double) * Nx);
        u 	= (double*) malloc(sizeof(double) * Nx); 
        etahat 	= (complex*) malloc(sizeof(complex) * Nx);
        uhat	= (complex*) malloc(sizeof(complex) * Nx); 	
        
        //define eta and u from z
        for (int j = 0; j < Nx; j++){
          etahat[j].Re 		= z[j];
          etahat[j].Im 		= z[Nx + j];
          uhat[j].Re   		= z[2*Nx + j];
          uhat[j].Im 	 	= z[3*Nx + j]; 
        }
        ifft_complex_real_nx(etahat, eta);
        ifft_complex_real_nx(uhat, u);		
        
        //wind
        if (strcmp(wind, "yes") == 0){  
            rhs_wind(Rw_hat, x, t, Cw, uhat);
        }		
        
        //friction 
        if (strcmp(friction, "yes") == 0){
            rhs_friction(Rf_hat, depth, eta, u, Nx, Cf, k);	
        }

        //wall
        if (strcmp(wall, "yes") == 0){
            fun_wallinfluxsource(t, eta);
        }

        //coupling
        if (strcmp(coupling, "yes") == 0){
            coupling_term(eta, u, couple_etahat, couple_uhat, t, etadata, udata);
        }
        
        free(eta);
        free(u); 
        free(etahat);
        free(uhat);
    }	

    //define source function
    complex* force = new complex[Nx];
    force = define_forcing(t, Nx);

    //solve the linear equation
    #pragma omp parallel for num_threads(ncore2)
    for (int i=0;i<Nx;i++){  
        dz[i]        = (deluH1[i].Re + Swall_hat[i].Re + 2 * force[i].Re - damping_eta[i][0] + couple_etahat[i].Re) * aal[i];
        dz[i + Nx]   = (deluH1[i].Im + Swall_hat[i].Im + 2 * force[i].Im - damping_eta[i][1] + couple_etahat[i].Im) * aal[i];  
        dz[i + 2*Nx] = (deletaH1[i].Re + Rf_hat[i].Re + Rw_hat[i].Re - damping_u[i][0] + couple_uhat[i].Re) * aal[i];
        dz[i + 3*Nx] = (deletaH1[i].Im + Rf_hat[i].Im + Rw_hat[i].Im - damping_u[i][1] + couple_uhat[i].Im) * aal[i]; 
        dy[i] 	     = dz[i];
        dy[i + Nx]   = dz[i + Nx];
        dy[i + 2*Nx] = dz[i + 2*Nx];
        dy[i + 3*Nx] = dz[i + 3*Nx];
    }

    fftw_free(damping_eta);
    free(deluH1);
    fftw_free(damping_u);
    free(deletaH1);
    delete[] force;
    delete[] couple_etahat;
    delete[] couple_uhat;

    return GSL_SUCCESS;
}

int rhs_HsF2 (double t, const double* z, double* dz, void *params)
{ 
    complex *deluH1, *deletaH1;
    deluH1  	= (complex*) malloc(sizeof(complex) * Nx);
    deletaH1  	= (complex*) malloc(sizeof(complex) * Nx);

    fftw_complex  *eta_u_hat, *LetaM0u_hat, *u2_M0u2_hat, *damping_eta, *damping_u;
    eta_u_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    u2_M0u2_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    complex *deluH2, *deletaH2;	
    deluH2  	= (complex*) malloc(sizeof(complex) * Nx);	
    deletaH2  	= (complex*) malloc(sizeof(complex) * Nx);

    complex* couple_etahat = new complex[Nx];
    complex* couple_uhat   = new complex[Nx];

    #pragma omp parallel for
    for(int j = 0; j < Nx; j++){
        couple_etahat[j].Re = 0;
        couple_uhat[j].Re   = 0;
        couple_etahat[j].Im = 0;
        couple_uhat[j].Im   = 0;
    }

    #pragma omp parallel sections num_threads(ncore2)
    {
        #pragma omp section
        //for HS1		
        define_delH1(z, deluH1, deletaH1);

        #pragma omp section
        //for HS2			
        if (strcmp(bath, "flat") == 0) {		
            ifft_fft_Hsf2(z, L, eta_u_hat, LetaM0u_hat, u2_M0u2_hat, damping_eta, damping_u, dampchar, Nx);	
        }
        else {//varying bathymetry				
            ifft_fft_Hsf2_v(z, L, eta_u_hat, LetaM0u_hat, u2_M0u2_hat, damping_eta, damping_u, dampchar, Nx);		
        }
    }

    //define nonlinear term
    #pragma omp parallel for num_threads(ncore2)
    for (int i = 0; i < Nx; i++){		
        deluH2[i].Re	= -LetaM0u_hat[i][0] - prod_1ik(k[i], eta_u_hat[i][0], eta_u_hat[i][1], REAL);
        deluH2[i].Im	= -LetaM0u_hat[i][1] - prod_1ik(k[i], eta_u_hat[i][0], eta_u_hat[i][1], IMAG);
        deletaH2[i].Re  = -prod_1ik(k[i], u2_M0u2_hat[i][0], u2_M0u2_hat[i][1], REAL) / 2;
        deletaH2[i].Im  = -prod_1ik(k[i], u2_M0u2_hat[i][0], u2_M0u2_hat[i][1], IMAG) / 2;
    }		

    //nonlinear adjustment
    complex* deluH2_adj  = (complex*) malloc(sizeof(complex) * Nx);
    complex* deletaH2_adj= (complex*) malloc(sizeof(complex) * Nx);
    adjust_delH(deluH2_adj, deletaH2_adj, deluH2, deletaH2);

    //join H1 and H2
    complex *deluH12, *deletaH12;
    deluH12  	= (complex*) malloc(sizeof(complex) * Nx);
    deletaH12  	= (complex*) malloc(sizeof(complex) * Nx);
	
    #pragma omp parallel for num_threads(ncore2)
    for (int i = 0; i < Nx; i++){ 		
        deluH12[i].Re  = deluH1[i].Re + deluH2_adj[i].Re;
        deluH12[i].Im  = deluH1[i].Im + deluH2_adj[i].Im; 
        deletaH12[i].Re= deletaH1[i].Re + deletaH2_adj[i].Re;
        deletaH12[i].Im= deletaH1[i].Im + deletaH2_adj[i].Im; 	
    }

    //IF wind and/or friction and/or breaking
    if ((strcmp(friction, "yes") == 0) || (strcmp(wind, "yes") == 0) || (strcmp(wall, "yes") == 0) || (strcmp(breaking, "yes") == 0) || (strcmp(coupling, "yes") == 0)) {
        double  *eta, *u;
        complex *etahat, *uhat;
        eta 	= (double*) malloc(sizeof(double) * Nx);
        u 	= (double*) malloc(sizeof(double) * Nx); 
        etahat 	= (complex*) malloc(sizeof(complex) * Nx);
        uhat	= (complex*) malloc(sizeof(complex) * Nx); 	

        //define eta and u from z
        for (int j = 0; j < Nx; j++){
            etahat[j].Re 	= z[j];
            etahat[j].Im 	= z[Nx + j];
            uhat[j].Re   	= z[2*Nx + j];
            uhat[j].Im 	 	= z[3*Nx + j]; 
        }
        
	ifft_complex_real_nx(etahat, eta);
	ifft_complex_real_nx(uhat, u);		
	//wind
	if (strcmp(wind, "yes") == 0){  
	    rhs_wind(Rw_hat, x, t, Cw, uhat);
	}		
	//friction 
	if (strcmp(friction, "yes") == 0){
	    rhs_friction(Rf_hat, depth, eta, u, Nx, Cf, k);	
	}		
        //wall
        if (strcmp(wall, "yes") == 0){
	    fun_wallinfluxsource(t, eta);
        }
        //breaking
        if (strcmp(breaking, "yes") == 0){
	    rhs_breaking(Rb_hat, eta, etahat, u, t, x, k, pb, deluH12);	
        }						
        free(eta);
        free(u);
        free(etahat);
        free(uhat);
        //coupling
        if (strcmp(coupling, "yes") == 0){
	    coupling_term(eta, u, couple_etahat, couple_uhat, t, etadata, udata);
        }
    }	
	
    // define source Forcing
    complex* force    = new complex[Nx];
    force = define_forcing(t, Nx);

    //solve nonlinear equations
    #pragma omp parallel for num_threads(ncore2)
    for (int i = 0; i < Nx; i++){ 		
        dz[i]        = (deluH12[i].Re + Swall_hat[i].Re + 2*force[i].Re - damping_eta[i][0] + couple_etahat[i].Re) * aal[i];
        dz[i + Nx]   = (deluH12[i].Im + Swall_hat[i].Im + 2*force[i].Im - damping_eta[i][1] + couple_etahat[i].Im) * aal[i];  
        dz[i + 2*Nx] = (deletaH12[i].Re + Rf_hat[i].Re + Rw_hat[i].Re + Rb_hat[i].Re - damping_u[i][0] + couple_uhat[i].Re) * aal[i];
        dz[i + 3*Nx] = (deletaH12[i].Im + Rf_hat[i].Im + Rw_hat[i].Im + Rb_hat[i].Im - damping_u[i][1] + couple_uhat[i].Im) * aal[i];
         
        dy[i] 	     = dz[i];
        dy[i + Nx]   = dz[i + Nx];
        dy[i + 2*Nx] = dz[i + 2*Nx];
        dy[i + 3*Nx] = dz[i + 3*Nx];
    }	
	
    free(deluH1);free(deletaH1);
    free(deluH12);free(deluH2);free(deletaH12);free(deletaH2);	
    free(deluH2_adj);free(deletaH2_adj);delete[] force;

    fftw_free(u2_M0u2_hat);
    fftw_free(damping_eta);fftw_free(damping_u);	
    fftw_free(eta_u_hat);fftw_free(LetaM0u_hat);	
    delete[] couple_etahat;
    delete[] couple_uhat;

    return GSL_SUCCESS;
}

int rhs_HsF3 (double t, const double* z, double* dz, void *params)
{ 	
    complex *deluH1, *deletaH1;
    deluH1  	= (complex*) malloc(sizeof(complex) * Nx);
    deletaH1  	= (complex*) malloc(sizeof(complex) * Nx);	

    fftw_complex *eta_u_hat, *LetaM0u_hat, *u2_M0u2_hat, *damping_eta, *damping_u;
    eta_u_hat 	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    u2_M0u2_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    fftw_complex *deletaH3_hat, *deluH3_hat;
    deluH3_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    deletaH3_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    complex* couple_etahat = new complex[Nx];
    complex* couple_uhat   = new complex[Nx];

    #pragma omp parallel for
    for(int j = 0; j < Nx; j++){
        couple_etahat[j].Re = 0;
        couple_uhat[j].Re   = 0;
        couple_etahat[j].Im = 0;
        couple_uhat[j].Im   = 0;
    }
    
    #pragma omp parallel sections
    {
        //for HS1
        #pragma omp section
        define_delH1(z, deluH1, deletaH1);
        
        //for HS2 and HS3	
        #pragma omp section
        if (strcmp(bath, "flat") == 0) {		
            ifft_fft_Hsf2(z, L, eta_u_hat, LetaM0u_hat, u2_M0u2_hat, damping_eta, damping_u, dampchar, Nx);	
            ifft_fft_Hsf3(z, L, LetaM0u_hat, deletaH3_hat, deluH3_hat, Nx);
        }
        else {//varying bathymetry				
            ifft_fft_Hsf2_v(z, L, eta_u_hat, LetaM0u_hat, u2_M0u2_hat, damping_eta, damping_u, dampchar, Nx);
            ifft_fft_Hsf3_v(z, L, LetaM0u_hat, deletaH3_hat, deluH3_hat, Nx);
        }
    }

    //nonlinear only term
    complex deluH2, deletaH2, deluH3, deletaH3;	
    complex* deluH23    = (complex*) malloc(sizeof(complex) * Nx);
    complex* deletaH23  = (complex*) malloc(sizeof(complex) * Nx);

    #pragma omp parallel for private(deluH2, deletaH2, deluH3, deletaH3)
    for (int i = 0; i < Nx; i++){				
        deluH2.Re	 = -LetaM0u_hat[i][0] - prod_1ik(k[i], eta_u_hat[i][0], eta_u_hat[i][1], REAL);
        deluH2.Im	 = -LetaM0u_hat[i][1] - prod_1ik(k[i], eta_u_hat[i][0], eta_u_hat[i][1], IMAG);
        deletaH2.Re      = -prod_1ik(k[i], u2_M0u2_hat[i][0], u2_M0u2_hat[i][1], REAL) / 2;
        deletaH2.Im      = -prod_1ik(k[i], u2_M0u2_hat[i][0], u2_M0u2_hat[i][1], IMAG) / 2;

        deluH3.Re        = deluH3_hat[i][0];
        deluH3.Im        = deluH3_hat[i][1];                            
        deletaH3.Re      = -prod_1ik(k[i], deletaH3_hat[i][0], deletaH3_hat[i][1], REAL);
        deletaH3.Im      = -prod_1ik(k[i], deletaH3_hat[i][0], deletaH3_hat[i][1], IMAG);

        deluH23[i].Re    = deluH2.Re + deluH3.Re;
        deluH23[i].Im    = deluH2.Im + deluH3.Im;
        deletaH23[i].Re  = deletaH2.Re + deletaH3.Re;
        deletaH23[i].Im  = deletaH2.Im + deletaH3.Im; 
    }
    fftw_free(eta_u_hat);
    fftw_free(u2_M0u2_hat);
    fftw_free(LetaM0u_hat);

    //nonlinear adjustment
    complex* deluH23_adj   = (complex*) malloc(sizeof(complex) * Nx);
    complex* deletaH23_adj = (complex*) malloc(sizeof(complex) * Nx);
    adjust_delH(deluH23_adj, deletaH23_adj, deluH23, deletaH23);

    //join HS1, HS2, and HS3
    complex* deluH123   = (complex*) malloc(sizeof(complex) * Nx);	
    complex* deletaH123 = (complex*) malloc(sizeof(complex) * Nx);

    #pragma omp parallel for
    for (int i = 0; i < Nx; i++){ 		
        deluH123[i].Re  = deluH1[i].Re + deluH23_adj[i].Re;
        deluH123[i].Im  = deluH1[i].Im + deluH23_adj[i].Im; 
        deletaH123[i].Re= deletaH1[i].Re + deletaH23_adj[i].Re;
        deletaH123[i].Im= deletaH1[i].Im + deletaH23_adj[i].Im; 		
    }
 
    //IF wind and/or friction and/or breaking
    if ((strcmp(friction, "yes") == 0) || (strcmp(wind, "yes") == 0) || (strcmp(breaking, "yes") == 0) || (strcmp(coupling, "yes") == 0)) {
        double  *eta, *u;
        complex *etahat, *uhat;
        eta 	= (double*) malloc(sizeof(double) * Nx);
        u 	= (double*) malloc(sizeof(double) * Nx); 
        etahat 	= (complex*) malloc(sizeof(complex) * Nx);
        uhat	= (complex*) malloc(sizeof(complex) * Nx); 	
        
        //define eta and u from z
        for (int j = 0; j < Nx; j++){
            etahat[j].Re = z[j];
            etahat[j].Im = z[Nx + j];
            uhat[j].Re   = z[2*Nx + j];
            uhat[j].Im 	 = z[3*Nx + j]; 
        }
        ifft_complex_real_nx(etahat, eta);
        ifft_complex_real_nx(uhat, u);		
        //wind
        if (strcmp(wind, "yes") == 0){  
            rhs_wind(Rw_hat, x, t, Cw, uhat);
        }				
        //friction 
        if (strcmp(friction, "yes") == 0){
            rhs_friction(Rf_hat, depth, eta, u, Nx, Cf, k);	
        }			
        //breaking
        if (strcmp(breaking, "yes") == 0){
            rhs_breaking(Rb_hat, eta, etahat, u, t, x, k, pb, deluH123);	
        }			
        free(eta);free(u);free(etahat);free(uhat);
        //coupling
        if (strcmp(coupling,"yes")==0){
	        coupling_term(eta,u,couple_etahat,couple_uhat,t,etadata,udata);
        }
    }	

    // define source Forcing
    complex* force = new complex[Nx];
    force = define_forcing(t,Nx); 

    //solving equations 
    #pragma omp parallel for
    for (int i=0;i<Nx;i++){ 	
        dz[i]        = (deluH123[i].Re + 2 * force[i].Re - damping_eta[i][0] + couple_etahat[i].Re) * aal[i];
        dz[i+Nx]     = (deluH123[i].Im + 2 * force[i].Im - damping_eta[i][1] + couple_etahat[i].Im) * aal[i];  
        dz[i + 2*Nx] = (deletaH123[i].Re + Rf_hat[i].Re + Rw_hat[i].Re + Rb_hat[i].Re - damping_u[i][0] + couple_uhat[i].Re) * aal[i];
        dz[i + 3*Nx] = (deletaH123[i].Im + Rf_hat[i].Im + Rw_hat[i].Im + Rb_hat[i].Im - damping_u[i][1] + couple_uhat[i].Im) * aal[i];
         
        dy[i] 	     = dz[i];
        dy[i + Nx]   = dz[i + Nx];
        dy[i + 2*Nx] = dz[i + 2*Nx];
        dy[i + 3*Nx] = dz[i + 3*Nx];
    }	

    fftw_free(damping_eta);fftw_free(damping_u);
    fftw_free(deletaH3_hat);fftw_free(deluH3_hat);
    free(deluH1);free(deletaH1); 
    free(deluH23);free(deletaH23); 
    free(deluH123);free(deletaH123); 
    free(deluH23_adj);free(deletaH23_adj); 
    delete[] force;
    delete[] couple_etahat;
    delete[] couple_uhat;

    return GSL_SUCCESS;
}

//hs4 
int rhs_HsF4 (double t, const double* z, double* dz, void *params)
{ 	
    complex *deluH1, *deletaH1;
    deluH1  	= (complex*) malloc(sizeof(complex) * Nx);
    deletaH1  	= (complex*) malloc(sizeof(complex) * Nx);	

    fftw_complex *eta_u_hat, *LetaM0u_hat, *u2_M0u2_hat, *damping_eta, *damping_u;
    eta_u_hat 	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    LetaM0u_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    u2_M0u2_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    damping_u   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    fftw_complex *deletaH3_hat, *deluH3_hat;
    deluH3_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    deletaH3_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    
    fftw_complex *deletaH4_hat, *deluH4_hat;
    deluH4_hat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    deletaH4_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    complex* couple_etahat = new complex[Nx];
    complex* couple_uhat   = new complex[Nx];

    #pragma omp parallel for
    for(int j = 0; j < Nx; j++){
        couple_etahat[j].Re = 0;
        couple_uhat[j].Re   = 0;
        couple_etahat[j].Im = 0;
        couple_uhat[j].Im   = 0;
    }

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            define_delH1(z, deluH1, deletaH1);
        }
        
        #pragma omp section
        {
            if (strcmp(bath, "flat") == 0) {	
                ifft_fft_Hsf2(z, L, eta_u_hat, LetaM0u_hat, u2_M0u2_hat, damping_eta, damping_u, dampchar, Nx);	
                ifft_fft_Hsf3(z, L, LetaM0u_hat, deletaH3_hat, deluH3_hat, Nx);
                ifft_fft_Hsf4(z, L, LetaM0u_hat, deletaH4_hat, deluH4_hat, Nx); //problem hs4
            }
        }
    }

    //nonlinear only term
    complex deluH2, deletaH2, deluH3, deletaH3, deluH4, deletaH4;	
    complex* deluH234    = (complex*) malloc(sizeof(complex) * Nx);
    complex* deletaH234  = (complex*) malloc(sizeof(complex) * Nx);

    #pragma omp parallel for private(deluH2, deletaH2, deluH3, deletaH3)
    for (int i = 0; i < Nx; i++){				
        deluH2.Re	     = -LetaM0u_hat[i][0] - prod_1ik(k[i], eta_u_hat[i][0], eta_u_hat[i][1], REAL);
        deluH2.Im	     = -LetaM0u_hat[i][1] - prod_1ik(k[i], eta_u_hat[i][0], eta_u_hat[i][1], IMAG);
        deletaH2.Re      = -prod_1ik(k[i], u2_M0u2_hat[i][0], u2_M0u2_hat[i][1], REAL) / 2;
        deletaH2.Im      = -prod_1ik(k[i], u2_M0u2_hat[i][0], u2_M0u2_hat[i][1], IMAG) / 2;

        deluH3.Re        = deluH3_hat[i][0];
        deluH3.Im        = deluH3_hat[i][1];                            
        deletaH3.Re      = -prod_1ik(k[i], deletaH3_hat[i][0], deletaH3_hat[i][1], REAL);
        deletaH3.Im      = -prod_1ik(k[i], deletaH3_hat[i][0], deletaH3_hat[i][1], IMAG);

        deluH4.Re        = deluH4_hat[i][0];
        deluH4.Im        = deluH4_hat[i][1];                            
        deletaH4.Re      = -prod_1ik(k[i], deletaH4_hat[i][0], deletaH4_hat[i][1], REAL);
        deletaH4.Im      = -prod_1ik(k[i], deletaH4_hat[i][0], deletaH4_hat[i][1], IMAG);

        deluH234[i].Re    = deluH2.Re + deluH3.Re + deluH4.Re;
        deluH234[i].Im    = deluH2.Im + deluH3.Im + deluH4.Im;
        deletaH234[i].Re  = deletaH2.Re + deletaH3.Re + deletaH4.Re;
        deletaH234[i].Im  = deletaH2.Im + deletaH3.Im + deletaH4.Im; 
    }
    fftw_free(eta_u_hat);
    fftw_free(u2_M0u2_hat);
    fftw_free(LetaM0u_hat);

    //nonlinear adjustment
    complex* deluH234_adj   = (complex*) malloc(sizeof(complex) * Nx);
    complex* deletaH234_adj = (complex*) malloc(sizeof(complex) * Nx);
    adjust_delH(deluH234_adj, deletaH234_adj, deluH234, deletaH234);

    //join HS1, HS2, HS3 and HS4
    complex* deluH1234   = (complex*) malloc(sizeof(complex) * Nx);	
    complex* deletaH1234 = (complex*) malloc(sizeof(complex) * Nx);

    #pragma omp parallel for
    for (int i = 0; i < Nx; i++){ 		
        deluH1234[i].Re  = deluH1[i].Re + deluH234_adj[i].Re;
        deluH1234[i].Im  = deluH1[i].Im + deluH234_adj[i].Im; 
        deletaH1234[i].Re= deletaH1[i].Re + deletaH234_adj[i].Re;
        deletaH1234[i].Im= deletaH1[i].Im + deletaH234_adj[i].Im; 		
    }
 
    //IF wind and/or friction and/or breaking
    if ((strcmp(friction, "yes") == 0) || (strcmp(wind, "yes") == 0) || (strcmp(breaking, "yes") == 0) || (strcmp(coupling, "yes") == 0)) {
        double  *eta, *u;
        complex *etahat, *uhat;
        eta 	= (double*) malloc(sizeof(double) * Nx);
        u 	= (double*) malloc(sizeof(double) * Nx); 
        etahat 	= (complex*) malloc(sizeof(complex) * Nx);
        uhat	= (complex*) malloc(sizeof(complex) * Nx); 	
        
        //define eta and u from z
        for (int j = 0; j < Nx; j++){
            etahat[j].Re = z[j];
            etahat[j].Im = z[Nx + j];
            uhat[j].Re   = z[2*Nx + j];
            uhat[j].Im 	 = z[3*Nx + j]; 
        }
        ifft_complex_real_nx(etahat, eta);
        ifft_complex_real_nx(uhat, u);		
        //wind
        if (strcmp(wind, "yes") == 0){  
            rhs_wind(Rw_hat, x, t, Cw, uhat);
        }				
        //friction 
        if (strcmp(friction, "yes") == 0){
            rhs_friction(Rf_hat, depth, eta, u, Nx, Cf, k);	
        }			
        //breaking
        if (strcmp(breaking, "yes") == 0){
	        rhs_breaking(Rb_hat, eta, etahat, u, t, x, k, pb, deluH1234);	
        }			
        free(eta);free(u);free(etahat);free(uhat);
        //coupling
        if (strcmp(coupling,"yes")==0){
	        coupling_term(eta,u,couple_etahat,couple_uhat,t,etadata,udata);
        }
    }	

    // define source Forcing
    complex* force = new complex[Nx];
    force = define_forcing(t,Nx); 

    //solving equations 
    #pragma omp parallel for
    for (int i=0;i<Nx;i++){ 	
        dz[i]        = (deluH1234[i].Re + 2 * force[i].Re - damping_eta[i][0] + couple_etahat[i].Re) * aal[i];
        dz[i+Nx]     = (deluH1234[i].Im + 2 * force[i].Im - damping_eta[i][1] + couple_etahat[i].Im) * aal[i];  
        dz[i + 2*Nx] = (deletaH1234[i].Re + Rf_hat[i].Re + Rw_hat[i].Re + Rb_hat[i].Re - damping_u[i][0] + couple_uhat[i].Re) * aal[i];
        dz[i + 3*Nx] = (deletaH1234[i].Im + Rf_hat[i].Im + Rw_hat[i].Im + Rb_hat[i].Im - damping_u[i][1] + couple_uhat[i].Im) * aal[i];
         
        dy[i] 	     = dz[i];
        dy[i + Nx]   = dz[i + Nx];
        dy[i + 2*Nx] = dz[i + 2*Nx];
        dy[i + 3*Nx] = dz[i + 3*Nx];
    }	

    fftw_free(damping_eta);fftw_free(damping_u);
    fftw_free(deletaH3_hat);fftw_free(deluH3_hat);
    fftw_free(deletaH4_hat);fftw_free(deluH4_hat);
    free(deluH1);free(deletaH1); 
    free(deluH234);free(deletaH234); 
    free(deluH1234);free(deletaH1234); 
    free(deluH234_adj);free(deletaH234_adj); 
    delete[] force;
    delete[] couple_etahat;
    delete[] couple_uhat;

    return GSL_SUCCESS;
} 

void rhs_friction(complex* Rf_hat, double* depth, double* eta, double* u, int Nx, double* Cf, double* k)
//Cf is array, zero for no friction
{
    double*  Rf   = (double*) malloc(sizeof(double) * Nx);	
    if (runup == 1){
        for (int i = 0; i < Nx; i++){
	    Rf[i] = -g * pow(Cf[i], 2) * u[i] * fabs(u[i]) / pow((depth[i] + eta[i]), 4/3);
        }
    }
    else {
        double Hhere;
        for (int j = 0; j < Nx; j++){
            Hhere 	= ChiAdj[j] * eta[j] + depth[j];
            if (Hhere > Hmin_interp){
	        Rf[j]   = -g * pow(Cf[j], 2) * u[j] * fabs(u[j]) / pow(Hhere, 4/3);
            }
            else {
	        Rf[j] = 0;
            }
        }
    }	

    fft_real_complex(Rf, Rf_hat, Nx);	
    free(Rf);			
}

void rhs_wind(complex* Rw_hat, double* x, double tnow, double* Cw, complex* u_hat)
{
    int      indt;
    double*  Rwtemp, *Rw_temp;
    complex* Rwtemp_hat;

    indt 	= closest(ttrace, tnow, n);	
    Rwtemp  	= new double[Nx];
    Rw_temp 	= new double[Nx];
    Rwtemp_hat 	= new complex[Nx];

    for (int j = 0; j < Nx; j++){
        Rwtemp_hat[j].Re = omega_i(k[j], depthinf) * u_hat[j].Re;
        Rwtemp_hat[j].Im = omega_i(k[j], depthinf) * u_hat[j].Im;
    }
    ifft_complex_real_nx(Rwtemp_hat, Rwtemp);

    for (int j = 0; j < Nx; j++){
	Rw_temp[j] = Cw[indt*Nx + j] * Rwtemp[j];		
    }
    fft_real_complex(Rw_temp, Rw_hat, Nx);	
    delete[] Rw_temp; delete[] Rwtemp; delete[] Rwtemp_hat;
}

void rhs_breaking(complex* Rb_hat, double* eta, complex* etahat, double* u, double t, double* x, double* k, parbreak pb, complex* dteta_hat)
{
    double*  dteta  = new double[Nx];
    double*  Flux   = new double[Nx];
    double*  dxFlux = new double[Nx];
    double*  Rb     = new double[Nx];

    breaking_process(eta, etahat, u, t, x, pb);		
    ifft_complex_real_nx(dteta_hat, dteta);
    
    for (int j = 0; j < Nx; j++){
        Flux[j] = -B[j] * pow(delb,2) * (depth[j] + eta[j]) * dteta[j] * dteta[j];
    }

    gradient(Flux, dx, Nx, dxFlux);

    if (runup==1){
        for (int j = 0; j < Nx; j++){
	    Rb[j] = 1 / (depth[j] + eta[j]) * dxFlux[j];
        }
    }
    else {
        double Hhere;
        for (int j = 0; j < Nx; j++){
            Hhere = ChiAdj[j] * eta[j] + depth[j];
            if (Hhere > Hmin_interp){
	        Rb[j] = 1 / (Hhere) * dxFlux[j];
            }
            else {
	        Rb[j] = 0;
            }
        }
    }		

    fft_real_complex_nx(Rb, Rb_hat);	

    delete[] dteta;
    delete[] Flux;
    delete[] Rb;
    delete[] dxFlux;
}

complex* define_forcing(double t, int Nx)
{
    complex* force = new complex[Nx];
    double ppval     = gsl_spline_eval(spline, t, acc);
    double ppvalskew = gsl_spline_eval(splineskew, t, accskew); 
    
    if (strcmp(wavename, "zero") != 0) {
        if (strcmp(propagation, "Uni+") == 0) {		
            for (int i = 0; i < Nx; i++){
                force[i].Re  = 0.5 * (ppval * (gamX[i].Re) - ppvalskew * (gamX_skew[i].Re));
                force[i].Im  = 0.5 * (ppval * (gamX[i].Im) - ppvalskew * (gamX_skew[i].Im));
            }
        }
        else if (strcmp(propagation, "Uni-") == 0) {
            for (int i = 0; i < Nx; i++){
                force[i].Re  = 0.5 * (ppval * (gamX[i].Re) + ppvalskew * (gamX_skew[i].Re));
                force[i].Im  = 0.5 * (ppval * (gamX[i].Im) + ppvalskew * (gamX_skew[i].Im));
            }
        }
        else if (strcmp(propagation, "Bi") == 0) {	
            for (int i = 0; i < Nx; i++){
                force[i].Re  = ppval * (gamX[i].Re);
                force[i].Im  = ppval * (gamX[i].Im);
            }
        }
    }
    else {//IVP
        for (int i = 0; i < Nx; i++){
            force[i].Re     = 0;
            force[i].Im     = 0;
        }
    }
	    
    return force;	
}

void adjust_delH(complex* deluH_adj, complex* deletaH_adj, complex* deluH, complex* deletaH)
{
    //nonlinear adjustment for delH2
    double* iftdeluH       = (double*) malloc(sizeof(double) * Nx);
    double* iftdeluH_adj   = (double*) malloc(sizeof(double) * Nx);	
    double* iftdeletaH     = (double*) malloc(sizeof(double) * Nx);
    double* iftdeletaH_adj = (double*) malloc(sizeof(double) * Nx);	

    ifft_complex_real_nx(deluH, iftdeluH);	
    ifft_complex_real_nx(deletaH, iftdeletaH);		

    for (int i = 0; i < Nx; i++){
        iftdeluH_adj[i]   = iftdeluH[i] * ChiAdj[i];
        iftdeletaH_adj[i] = iftdeletaH[i] * ChiAdj[i];
    }
    fft_real_complex_nx(iftdeluH_adj, deluH_adj);
    fft_real_complex_nx(iftdeletaH_adj, deletaH_adj);		

    free(iftdeluH);free(iftdeluH_adj);
    free(iftdeletaH);free(iftdeletaH_adj);
}

int rhs_shore(double t, const double* z, double* dz, void *params)
{
    double  *eta, *u, *Upeta, *Upu;
    complex *etahat, *uhat, *Upetahat, *Upuhat;
    eta 	= (double*)  malloc(sizeof(double) * Nx);
    etahat 	= (complex*) malloc(sizeof(complex) * Nx);	
    u    	= (double*)  malloc(sizeof(double) * Nx);
    uhat 	= (complex*) malloc(sizeof(complex) * Nx);
    Upeta 	= (double*)  malloc(sizeof(double) * Nx);
    Upetahat    = (complex*) malloc(sizeof(complex) * Nx);	
    Upu    	= (double*)  malloc(sizeof(double) * Nx);
    Upuhat 	= (complex*) malloc(sizeof(complex) * Nx);	
	    
    //define eta from z
    #pragma omp parallel for num_threads(ncore2)
    for (int j = 0; j < Nx; j++){
        etahat[j].Re 	= z[j];
        etahat[j].Im 	= z[Nx + j]; 
        uhat[j].Re 	= z[2*Nx + j];
        uhat[j].Im 	= z[3*Nx + j];
        Upetahat[j].Re 	= Cpeak * z[j];
        Upetahat[j].Im 	= Cpeak * z[Nx + j]; 
        Upuhat[j].Re 	= Cpeak * z[2*Nx + j];
        Upuhat[j].Im 	= Cpeak * z[3*Nx + j]; 
    }
    ifft_complex_real_nx(etahat, eta);
    ifft_complex_real_nx(uhat, u);
    ifft_complex_real_nx(Upetahat, Upeta);
    ifft_complex_real_nx(Upuhat, Upu);
    free(Upetahat);
    free(Upuhat);

    //preparation
    int indDI, indDF, FlagOvertop;
    double* Hhere 	= new double[Nx];
    double* HeavH 	= new double[Nx];
    double* dampcharH 	= new double[Nx];	
    double  Hmin_shore 	= Hmin_interp / 5;

    for (int j = 0; j < Nx; j++){
        Hhere[j] = ChiAdj[j] * eta[j] + depth[j];
        if ((Hhere[j] - Hmin_shore) >= 0){
	    HeavH[j] = 1;//
        }
        else {
	    HeavH[j] = 0;//
        }
        dampcharH[j] = dampchar[j];
    }
    indDI = indexOfFirstZero(HeavH, Nx);
    indDF = indexOfLastZero(HeavH, Nx);

    //FlagOvertop
    if (indDI == -1){//no dry area
        FlagOvertop	= 0; 
        signProp	= 0;
    }
    else if(indDI == 0) {//to the left     
        FlagOvertop	= 0;
        signProp	= -1;
    }
    else if(indDF == (Nx - 1)) {//to the right
        FlagOvertop	= 0;
        signProp	= 1;
    }
    else if ((indDI > 0) && (indDF < (Nx - 1))){
        FlagOvertop	= 1;
        signProp	= 0;
    }

    //for direct model
    double *gm1 = new double[Nx]; 
    double *gp1 = new double[Nx]; 
    double *gc1 = new double[Nx];
    double *gm2 = new double[Nx]; 
    double *gp2 = new double[Nx]; 
    double *gc2 = new double[Nx];

    set_init_double(gm1, Nx);
    set_init_double(gp1, Nx);
    set_init_double(gc1, Nx);
    set_init_double(gm2, Nx);
    set_init_double(gp2, Nx);
    set_init_double(gc2, Nx);	

    if (FlagOvertop == 0){
        for (int j = 0; j < Nx; j++){
            if ((Hhere[j] >= Hmin_interp) && (Hhere[j] <= Hplus_interp)) {	
                gm1[j] 	= gsl_spline_eval(sg_m1, Hhere[j], acc_m1);
                gp1[j] 	= gsl_spline_eval(sg_p1, Hhere[j], acc_p1);
                gc1[j] 	= gsl_spline_eval(sg_c1, Hhere[j], acc_c1);
                gm2[j] 	= gsl_spline_eval(sg_m2, Hhere[j], acc_m2);
                gp2[j] 	= gsl_spline_eval(sg_p2, Hhere[j], acc_p2);
                gc2[j] 	= gsl_spline_eval(sg_c2, Hhere[j], acc_c2);				
            }
        }
    }

    //update dampchar here
    if (indDI != -1){
        for (int j = indDI; j < (indDF + 1); j++){
            dampcharH[j]  = 1;
            Upeta[j]      = Hhere[j];
            Upu[j] 	  = u[j];	
        }
    }

    complex* CCu_hat1 	       = (complex*) malloc(sizeof(complex) * Nx);
    complex* CCu_hat2 	       = (complex*) malloc(sizeof(complex) * Nx);
    complex* CCu_hat 	       = (complex*) malloc(sizeof(complex) * Nx);
    double*  CCu_ShoreChar     = (double*)  malloc(sizeof(double) * Nx);
    complex* CCu_ShoreChar_hat = (complex*) malloc(sizeof(complex) * Nx);
    double*  kC2_u  	       = (double*)  malloc(sizeof(double) * Nx);
    complex* kC2_u_hat         = (complex*) malloc(sizeof(complex) * Nx);
    HSS_runup(uhat, CCu_hat1, C2m1, C2p1, C2c1, gm1, gp1, gc1);
    HSS_runup(uhat, CCu_hat2, C2m2, C2p2, C2c2, gm2, gp2, gc2);
    for (int j = 0; j < Nx; j++){
        CCu_hat[j].Re = CCu_hat1[j].Re + CCu_hat2[j].Re;
        CCu_hat[j].Im = CCu_hat1[j].Im + CCu_hat2[j].Im;
    }
    free(CCu_hat1);
    free(CCu_hat2);
    ifft_complex_real_nx(CCu_hat, CCu_ShoreChar);	
    for (int j = 0; j < Nx; j++){
        if(Hhere[j] <= Hmin_interp){
	    CCu_ShoreChar[j] = g * Hhere[j] * u[j];
        }
    }	
    fft_real_complex_nx(CCu_ShoreChar, CCu_ShoreChar_hat);
    for (int j = 0; j < Nx;j++){
        kC2_u_hat[j].Re = k[j] * CCu_ShoreChar_hat[j].Re;
        kC2_u_hat[j].Im = k[j] * CCu_ShoreChar_hat[j].Im;		
    }
    
    ifft_complex_real_nx(kC2_u_hat, kC2_u);
    free(CCu_hat);
    free(CCu_ShoreChar);
    free(kC2_u_hat);

    //dynamic model
    complex *deluH_hat, *deletaH_hat, *term2_hat;
    double*  term2;
    deluH_hat  	 = (complex*) malloc(sizeof(complex) * Nx);
    deletaH_hat	 = (complex*) malloc(sizeof(complex) * Nx);
    term2_hat	 = (complex*) malloc(sizeof(complex) * Nx);
    term2	 = (double*) malloc(sizeof(double) * Nx);
    for (int j = 0; j < Nx; j++){
        term2[j] = (pow(u[j], 2) / 2 - pow(kC2_u[j]/g, 2) / 2) * ChiAdj[j];	
    }
    
    fft_real_complex_nx(term2, term2_hat);
    free(term2);
    free(kC2_u);
    
    //define delH	
    for (int j = 0; j < Nx; j++){
        deluH_hat[j].Re 	= -prod_1ik(k[j], 1 / g * CCu_ShoreChar_hat[j].Re, 1 /g * CCu_ShoreChar_hat[j].Im, REAL);
        deluH_hat[j].Im 	= -prod_1ik(k[j], 1 / g * CCu_ShoreChar_hat[j].Re, 1 / g *CCu_ShoreChar_hat[j].Im, IMAG);
        deletaH_hat[j].Re 	= -prod_1ik(k[j], g * z[j] + term2_hat[j].Re, g * z[j + Nx] + term2_hat[j].Im, REAL);
        deletaH_hat[j].Im 	= -prod_1ik(k[j], g * z[j] + term2_hat[j].Re, g * z[j + Nx] + term2_hat[j].Im, IMAG);
    }
    
    free(CCu_ShoreChar_hat);
    free(term2_hat);	

    //init
    complex* force          = new complex[Nx];
    complex* damping_etahat = (complex*) malloc(sizeof(complex) * Nx);
    complex* damping_uhat   = (complex*) malloc(sizeof(complex) * Nx);
    double*  damping_eta    = (double*) malloc(sizeof(double) * Nx);
    double*  damping_u      = (double*) malloc(sizeof(double) * Nx);
    complex* coef           = new complex[Nx];

    #pragma omp parallel sections num_threads(ncore2)
    {
        #pragma omp section
        {
            // define source Forcing
            force = define_forcing(t, Nx);
        }

        #pragma omp section
        {
            //define damping_eta and damping_u	
            for (int i = 0; i < Nx; i++){
                damping_eta[i] 	= Upeta[i] * dampcharH[i];
                damping_u[i]    = Upu[i] * dampcharH[i];
            }
            fft_real_complex_nx(damping_eta, damping_etahat);
            fft_real_complex_nx(damping_u, damping_uhat);
        }

        #pragma omp section
        {
            //coefficient
            for (int j = 0; j < Nx; j++){
                coef[j].Re = 1;
                coef[j].Im = signProp * k[j];
            }
        }
    }


    //IF wind and/or friction and/or breaking
    if ((strcmp(friction, "yes") == 0) || (strcmp(wind, "yes") == 0) || (strcmp(breaking, "yes") == 0)) {
        //wind
        if (strcmp(wind, "yes") == 0){  
	    rhs_wind(Rw_hat, x, t, Cw, uhat);
        }		
        //friction 
        if (strcmp(friction, "yes") == 0){
	    rhs_friction(Rf_hat, depth, eta, u, Nx, Cf, k);	
        }		
        //breaking
        if (strcmp(breaking, "yes") == 0){
	    rhs_breaking(Rb_hat, eta, etahat, u, t, x, k, pb, deluH_hat);	
        }						
    }

    //solve nonlinear equations
    #pragma omp parallel for num_threads(ncore2)
    for (int i = 0; i < Nx; i++){ 
        dz[i]        = (deluH_hat[i].Re + 2 * force[i].Re - prod_fftw(coef[i], damping_etahat[i].Re, damping_etahat[i].Im, REAL)) * aal[i];
        dz[i + Nx]   = (deluH_hat[i].Im + 2 * force[i].Im - prod_fftw(coef[i], damping_etahat[i].Re, damping_etahat[i].Im, IMAG)) * aal[i];  
        dz[i + 2*Nx] = (deletaH_hat[i].Re + Rf_hat[i].Re + Rw_hat[i].Re + Rb_hat[i].Re - prod_fftw(coef[i], damping_uhat[i].Re, damping_uhat[i].Im, REAL)) * aal[i];
        dz[i + 3*Nx] = (deletaH_hat[i].Im + Rf_hat[i].Im + Rw_hat[i].Im + Rb_hat[i].Im - prod_fftw(coef[i], damping_uhat[i].Re, damping_uhat[i].Im, IMAG)) * aal[i];
         
        dy[i] 	     = dz[i];
        dy[i + Nx]   = dz[i + Nx];
        dy[i + 2*Nx] = dz[i + 2*Nx];
        dy[i + 3*Nx] = dz[i + 3*Nx];
    }	

    delete[] coef;
    delete[] Hhere;delete[] HeavH;
    free(eta);free(u);free(etahat);free(uhat);
    free(deluH_hat);free(deletaH_hat);
    free(damping_eta);free(damping_etahat);
    free(damping_u);free(damping_uhat);	
    delete[] force;free(Upeta);free(Upu);
    delete[] dampcharH;

    return GSL_SUCCESS;
}


int rhs_shore_old(double t, const double* z, double* dz, void *params)
{
    double  *eta, *u, *Upeta, *Upu;
    complex *etahat, *uhat, *Upetahat, *Upuhat;
    eta 	= (double*)  malloc(sizeof(double) * Nx);
    etahat 	= (complex*) malloc(sizeof(complex) * Nx);	
    u    	= (double*)  malloc(sizeof(double) * Nx);
    uhat 	= (complex*) malloc(sizeof(complex) * Nx);
    Upeta 	= (double*)  malloc(sizeof(double) * Nx);
    Upetahat    = (complex*) malloc(sizeof(complex) * Nx);	
    Upu    	= (double*)  malloc(sizeof(double) * Nx);
    Upuhat 	= (complex*) malloc(sizeof(complex) * Nx);	
	    
    //define eta from z
    #pragma omp parallel for num_threads(ncore2)
    for (int j = 0; j < Nx; j++){
        etahat[j].Re 	= z[j];
        etahat[j].Im 	= z[Nx + j]; 
        uhat[j].Re      = z[2*Nx + j];
        uhat[j].Im 	= z[3*Nx + j];
        Upetahat[j].Re 	= Cpeak * z[j];
        Upetahat[j].Im 	= Cpeak * z[Nx + j]; 
        Upuhat[j].Re 	= Cpeak * z[2*Nx + j];
        Upuhat[j].Im 	= Cpeak * z[3*Nx + j]; 
    }
    ifft_complex_real_nx(etahat, eta);
    ifft_complex_real_nx(uhat, u);
    ifft_complex_real_nx(Upetahat, Upeta);
    ifft_complex_real_nx(Upuhat, Upu);
    free(Upetahat);
    free(Upuhat);

    //preparation
    double* Hhere 	= new double[Nx];
    double* HeavH 	= new double[Nx];
    double* ShoreChar 	= new double[Nx];
    double* dampcharH 	= new double[Nx];
    int indDI, indDF, FlagOvertop;

    for (int j = 0; j < Nx; j++){
        Hhere[j] = ChiAdj[j] * eta[j] + depth[j];
        if ((Hhere[j] - Hmin_interp) > 0){
	    HeavH[j] = 1;//
        }
        else {
	    HeavH[j] = 0;//
        }
        ShoreChar[j] = HeavH[j];
        dampcharH[j] = dampchar[j];
    }
    indDI = indexOfFirstZero(HeavH, Nx);
    indDF = indexOfLastZero(HeavH, Nx);

    //FlagOvertop
    if (indDI == -1){//no dry area
        FlagOvertop	 = 0; 
        signProp	 = 0;
    }
    else if(indDI == 0) {//to the left
        for (int j = indDI; j < (indDF + 1); j++){
	    ShoreChar[j] = 0;
        }   
        
        FlagOvertop	 = 0;
        signProp	 = -1;
    }
    else if(indDF == (Nx - 1)) {//to the right
        for (int j = indDI; j < (indDF + 1); j++){
	    ShoreChar[j] = 0;
        }
        FlagOvertop	= 0;
        signProp	= 1;
    }
    else if ((indDI > 0) && (indDF < (Nx - 1))){
        FlagOvertop	= 1;
        signProp	= 0;
    }

    //for direct model
    double *gm 		= new double[Nx]; 
    double *gp 		= new double[Nx]; 
    double *gc 		= new double[Nx];

    set_init_double(gm, Nx);
    set_init_double(gp, Nx);
    set_init_double(gc, Nx);
    for (int j = 0; j < Nx; j++){		
        if (ShoreChar[j] == 1){
            gm[j] 	= gsl_spline_eval(spline_sgm, Hhere[j], acc_sgm);
            gp[j] 	= gsl_spline_eval(spline_sgp, Hhere[j], acc_sgp);
            gc[j] 	= gsl_spline_eval(spline_sgc, Hhere[j], acc_sgc);
        }
    }	

    //update dampchar here
    if (indDI != -1){
        for (int j = indDI; j < (indDF + 1); j++){
            if (ShoreChar[j] == 0){
                dampcharH[j] = 1;
                Upeta[j]     = Hhere[j];
                Upu[j] 	     = u[j];	
            }
        }
    }

    complex* Cu_hat 	            = (complex*) malloc(sizeof(complex) * Nx);
    double*  Cu     	            = (double*)  malloc(sizeof(double) * Nx);
    complex* Cder_u_hat             = (complex*) malloc(sizeof(complex) * Nx);
    double*  Cder_u                 = (double*)  malloc(sizeof(double) * Nx);
    complex* Cu_ShoreChar_hat       = (complex*) malloc(sizeof(complex) * Nx);
    double*  Cu_ShoreChar           = (double*)  malloc(sizeof(double) * Nx);
    complex* ConjC_Cu_ShoreChar_hat = (complex*) malloc(sizeof(complex) * Nx);	
    
    HSS_runup(uhat, Cu_hat, Cm, Cp, Cc, gm, gp, gc);
    ifft_complex_real_nx(Cu_hat, Cu);
    
    for (int j = 0; j < Nx; j++){
        Cu_ShoreChar[j] = Cu[j] * ShoreChar[j];
    }
    
    fft_real_complex_nx(Cu_ShoreChar, Cu_ShoreChar_hat);
    HSS_runup(Cu_ShoreChar_hat, ConjC_Cu_ShoreChar_hat, Cm, Cp, Cc, gm, gp, gc);
    HSS_runup(uhat, Cder_u_hat, Cderm, Cderp, Cderc, gm, gp, gc);
    ifft_complex_real_nx(Cder_u_hat, Cder_u);
    
    delete[] gm; delete[] gp;delete[] gc;
    
    free(Cu_ShoreChar);free(Cder_u_hat);free(Cu_hat);free(Cu_ShoreChar_hat);

    //dynamic model
    complex *deluH_hat, *deletaH_hat, *Cderu_Cu_Chi_hat;
    double  *Cderu_Cu_Chi;
    Cderu_Cu_Chi  		= (double*)  malloc(sizeof(double) * Nx);
    Cderu_Cu_Chi_hat 	        = (complex*) malloc(sizeof(complex) * Nx);
    deluH_hat  			= (complex*) malloc(sizeof(complex) * Nx);
    deletaH_hat			= (complex*) malloc(sizeof(complex) * Nx);

    for (int j = 0;j < Nx; j++){
        Cderu_Cu_Chi[j] 	= Cder_u[j] * Cu[j] * ChiAdj[j];	
    }
    fft_real_complex_nx(Cderu_Cu_Chi, Cderu_Cu_Chi_hat);
    free(Cu);free(Cder_u);free(Cderu_Cu_Chi);

    //define delH
    complex* deleta_temp_h   = new complex[Nx];
    complex* deletaH_hat_tmp = new complex[Nx];
    double*  deleta_temp     = new double[Nx];
    double*  deletaH         = new double[Nx];

    #pragma omp parallel for num_threads(ncore2)
    for (int j = 0;j < Nx; j++){
        deluH_hat[j].Re     = -prod_1ik(k[j], 1 / g * ConjC_Cu_ShoreChar_hat[j].Re, 1 / g * ConjC_Cu_ShoreChar_hat[j].Im, REAL);
        deluH_hat[j].Im     = -prod_1ik(k[j], 1 / g * ConjC_Cu_ShoreChar_hat[j].Re, 1 / g * ConjC_Cu_ShoreChar_hat[j].Im, IMAG);
        deleta_temp_h[j].Re = g * z[j] + 1 / g * Cderu_Cu_Chi_hat[j].Re;
        deleta_temp_h[j].Im = g * z[j + Nx] + 1 / g * Cderu_Cu_Chi_hat[j].Im;
    }
    ifft_complex_real_nx(deleta_temp_h, deleta_temp);
    delete[] deleta_temp_h;
    free(ConjC_Cu_ShoreChar_hat);
    free(Cderu_Cu_Chi_hat);

    for (int j = 0; j < Nx; j++){
        deletaH[j] = deleta_temp[j] * ShoreChar[j];
    }	
    fft_real_complex_nx(deletaH, deletaH_hat_tmp);
    delete[] deleta_temp;
    delete[] deletaH;

    for (int j = 0;j < Nx; j++){
        deletaH_hat[j].Re 	= -prod_1ik(k[j], deletaH_hat_tmp[j].Re, deletaH_hat_tmp[j].Im, REAL);
	deletaH_hat[j].Im 	= -prod_1ik(k[j], deletaH_hat_tmp[j].Re, deletaH_hat_tmp[j].Im, IMAG);
    }
    delete[] deletaH_hat_tmp;	

    //init
    complex* force          = new complex[Nx];
    complex* damping_etahat = (complex*) malloc(sizeof(complex) * Nx);
    complex* damping_uhat   = (complex*) malloc(sizeof(complex) * Nx);
    double*  damping_eta    = (double*) malloc(sizeof(double) * Nx);
    double*  damping_u      = (double*) malloc(sizeof(double) * Nx);
    complex* coef = new complex[Nx];

    #pragma omp parallel sections num_threads(ncore2)
    {
        #pragma omp section
        {
            // define source Forcing
            force = define_forcing(t, Nx);
        }

        #pragma omp section
        {
            //define damping_eta and damping_u	
            for (int i=0; i < Nx; i++){
                damping_eta[i] 	= Upeta[i] * dampcharH[i];
                damping_u[i]    = Upu[i] * dampcharH[i];
            }
            fft_real_complex_nx(damping_eta, damping_etahat);
            fft_real_complex_nx(damping_u, damping_uhat);
        }

        #pragma omp section
        {
            //coefficient
            for (int j = 0; j < Nx; j++){
                coef[j].Re = 1;
                coef[j].Im = signProp * k[j];
            }
        }
    }

    //IF wind and/or friction and/or breaking
    if ((strcmp(friction, "yes") == 0) || (strcmp(wind, "yes") == 0) || (strcmp(breaking, "yes") == 0)) {
        //wind
        if (strcmp(wind, "yes") == 0){  
	    rhs_wind(Rw_hat, x, t, Cw, uhat);
        }		
        //friction 
        if (strcmp(friction, "yes") == 0){
	    rhs_friction(Rf_hat, depth, eta, u, Nx, Cf, k);	
        }		
        //breaking
        if (strcmp(breaking, "yes") == 0){
	    rhs_breaking(Rb_hat, eta, etahat, u, t, x, k, pb, deluH_hat);	
        }						
    }	

    //solve nonlinear equations
    #pragma omp parallel for num_threads(ncore2)
    for (int i = 0; i < Nx; i++){ 
        dz[i]        = (deluH_hat[i].Re + 2 * force[i].Re - prod_fftw(coef[i], damping_etahat[i].Re, damping_etahat[i].Im, REAL)) * aal[i];
        dz[i + Nx]   = (deluH_hat[i].Im + 2 * force[i].Im - prod_fftw(coef[i], damping_etahat[i].Re, damping_etahat[i].Im, IMAG)) * aal[i];  
        dz[i + 2*Nx] = (deletaH_hat[i].Re + Rf_hat[i].Re + Rw_hat[i].Re + Rb_hat[i].Re - prod_fftw(coef[i], damping_uhat[i].Re, damping_uhat[i].Im, REAL)) * aal[i];
        dz[i + 3*Nx] = (deletaH_hat[i].Im + Rf_hat[i].Im + Rw_hat[i].Im + Rb_hat[i].Im - prod_fftw(coef[i], damping_uhat[i].Re, damping_uhat[i].Im, IMAG)) * aal[i];
         
        dy[i] 	     = dz[i];
        dy[i + Nx]   = dz[i + Nx];
        dy[i + 2*Nx] = dz[i + 2*Nx];
        dy[i + 3*Nx] = dz[i + 3*Nx];
    }

    delete[] coef;
    delete[] Hhere;delete[] HeavH;delete[] ShoreChar;
    free(eta);free(u);free(etahat);free(uhat);
    free(deluH_hat);free(deletaH_hat);
    free(damping_eta);free(damping_etahat);
    free(damping_u);free(damping_uhat);	
    delete[] force;free(Upeta);free(Upu);
    delete[] dampcharH;

    return GSL_SUCCESS;
}

void fun_wallinfluxsource(double time, double* eta)
{
    if (time > tprevW){
        Swall_tS_eta = refl_coef * eta[idx_wall];
        if (ITERdtW < 1){
	    Swt_skew = 0;
        }
        else {
	    Swt_skew = Swt_skew + (Swall_tS_eta + Swall_tS_eta_prev) * (time - tprevW) / 2;
        }
    }	

    for (int j = 0; j < Nx; j++){
        Swall_hat[j].Re = 2 * Swt_skew * gamWall_skew_hat[j].Re;
        Swall_hat[j].Im = 2 * Swt_skew * gamWall_skew_hat[j].Im;
    }

    //update iter
    if (time > tprevW){
        ITERdtW 	    = ITERdtW + 1;
        tprevW		    = time; 
        Swall_tS_eta_prev   = Swall_tS_eta;
    }
}
