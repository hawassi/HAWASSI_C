void kinematic_modul(double t, double* eta, complex* u_hat, double* dteta, complex* dtphi_hat) 
{
    //local header
    double   nupeak = nu_peak;
    int      indx1 = indxIP1, indx2 = indxIP2;
    double   Hplus, Hmin, dhtot, DZmin, DZplus;

    /* set constant and variables */
    int const NH = 100;

    gsl_interp_accel *accHm0;
    gsl_spline *splineHm0;

      // initialize and read input data at t		
    double  gam_Hmin0[Nx], gam_Hplus0[Nx], gam_DZmin0[Nx], gam_DZplus0[Nx];	

    //define H
    double Htot[NH + 1];
    if (isvalueinarray(0, bathy, Nx)){//this is for runup
        int indxshore, Nxnew;
        double *subeta, *subD; 
        if (depth[0] > 0) {
            indxshore 	= findfirstzero(depth, Nx);
            Nxnew 	= indxshore + 1;
            subeta      = new double[Nxnew];
            subD        = (double*) malloc(sizeof(double) * Nxnew);
            
            subarrayx(Nx, eta, 0, Nxnew, subeta);
            subvector(depth, 0, Nxnew, subD);
            
            Hplus 	= search_v(Nxnew, subeta, MAXIMUM) + search_v(Nxnew, subD, MAXIMUM);
            delete[](subeta);free(subD);			
        }
        else {
            indxshore 	= findfirstnegative(depth, Nx);
            Nxnew 	= Nx - indxshore;	
            subeta 	= new double[Nxnew];
            subD 	= (double*) malloc(sizeof(double) * Nxnew);
            
            subarrayx(Nx, eta, indxshore, Nxnew, subeta);
            subvector(depth, indxshore, Nxnew, subD);
            
            Hplus 	= search_v(Nxnew, subeta, MAXIMUM) + search_v(Nxnew, subD, MAXIMUM);
            delete[](subeta);free(subD);	
        }
    }
    else {
	Hplus 	= maxEta + maxD;
    }

    Hmin 	= minEta + minD;	
    dhtot 	= (Hplus - Hmin) / NH;	

    for (int i = 0; i < (NH + 1); i++) { 
        Htot[i] = Hmin + i * dhtot;
    }		

    sinhcosh2IP(NH + 1, Htot, Hmin, Hplus, Nx, k, nupeak, gam_Hmin0, gam_Hplus0);		

    accHm0  		= gsl_interp_accel_alloc();
    splineHm0 	  	= gsl_spline_alloc(gsl_interp_cspline, NH + 1);
    gsl_spline_init(splineHm0, Htot, gam_Hmin0, NH + 1);

    gsl_interp_accel *accHp0  = gsl_interp_accel_alloc();
    gsl_spline *splineHp0     = gsl_spline_alloc(gsl_interp_cspline, NH + 1);
    gsl_spline_init(splineHp0, Htot, gam_Hplus0, NH + 1); 
      
      //define for DZ
    double 	DZtot[NH + 1];
    DZmin 	= 0;
    DZplus 	= Hplus;
    for (int i = 0; i < (NH + 1); i++) { 
        DZtot[i] = DZmin + i * (DZplus - DZmin) / NH; 
    }

    sinhcosh2IP(NH + 1, DZtot, DZmin, DZplus, Nx, k, nupeak, gam_DZmin0, gam_DZplus0);

    gsl_interp_accel *accDZm0  = gsl_interp_accel_alloc();
    gsl_spline *splineDZm0     = gsl_spline_alloc(gsl_interp_cspline, NH + 1);
    gsl_spline_init(splineDZm0, DZtot, gam_DZmin0, NH + 1);

    gsl_interp_accel *accDZp0  = gsl_interp_accel_alloc();
    gsl_spline *splineDZp0     = gsl_spline_alloc(gsl_interp_cspline, NH + 1);
    gsl_spline_init(splineDZp0, DZtot, gam_DZplus0, NH + 1);

    //allocate memory for output variable 3D    
    double (dtPhi)[NxoutIP][NzoutIP];
    double (dxPhi)[NxoutIP][NzoutIP];
    double (dzPhi)[NxoutIP][NzoutIP];
    double (Ptot)[NxoutIP][NzoutIP];
    double (dt_dxPHI)[NxoutIP][NzoutIP];
    double (dt_dzPHI)[NxoutIP][NzoutIP];            

    //calculate the interior properties at time i
    complex phihat_i[Nx];
    double  H[Nx], grad_H[Nx], grad_D[Nx];
    double  cosh_kH[Nx], tanh_kH[Nx];
    double  minEta_i, gam_Hmin, gam_Hplus;
    for (int j = 0; j < Nx; j++) {		
        H[j]			= (eta[j] + depth[j]);
        gam_Hmin 		= gsl_spline_eval(splineHm0, H[j], accHm0);
        gam_Hplus		= gsl_spline_eval(splineHp0, H[j], accHp0);
        cosh_kH[j]		= gam_Hmin * cosh(k[j] * Hmin) + gam_Hplus * cosh(k[j] * Hplus);
        tanh_kH[j]		= gam_Hmin * tanh(k[j] * Hmin) + gam_Hplus * tanh(k[j] * Hplus);
        if (k[j] == 0){
            phihat_i[j].Re 	= 0;
            phihat_i[j].Im 	= 0;
        }
        else {
            phihat_i[j].Re 	= (u_hat[j].Im)/k[j];
            phihat_i[j].Im 	= -(u_hat[j].Re)/k[j];
        }
    }//end x(j) loop		
    minEta_i 	= search_v(NxoutIP, eta, MINIMUM);	
    gradient(H, dx, Nx, grad_H);
    gradient(depth, dx, Nx, grad_D);

    #pragma omp parallel
    {
        //define new variable for computing the interior properties at each zz
        double  DZ[Nx];	
        double  DZplus, DZmin;
        double  gam_DZmin, gam_DZplus, ZZ, ratio;
        double  cosh_kDZ[Nx], sinh_kDZ[Nx];	
        complex	factor_nonifft[Nx], tempdtPhi_nonifft[Nx], tempdxPhi_nonifft[Nx], tempdzPhi_nonifft[Nx];
        double 	factor[Nx], tempdtPhi_ifft[Nx], tempdxPhi_ifft[Nx], tempdzPhi[Nx];	
        double 	tempdtPhi[Nx], tempdxPhi[Nx];	
        
        #pragma omp for
        for (int zz = 0; zz < NzoutIP; zz++) {			
            if (minD == maxD) {//flat bottom
                for (int j = 0; j < Nx; j++){
                    cosh_kDZ[j]   = cosh(k[j] * (zoutIP[zz] + depth[0]));
                    sinh_kDZ[j]   = sinh(k[j] * (zoutIP[zz] + depth[0]));
                }//end x(j)
            }//endif compute cosh sinh
            else {//nonflat bottom
                for (int j = 0; j < Nx; j++){					
                    DZ[j] = depth[j] + zoutIP[zz];
                    if (DZ[j] < 0) {DZ[j] = 0;}//DZ cannot be negative
                }//end x(j)
                
                DZplus 	= search_v(Nx, DZ, MAXIMUM);           
                DZmin 	= search_v(Nx, DZ, MINIMUM);
                if (DZmin < 0) {DZmin = 0;}
                
                if (DZmin == DZplus) {
                    for (int j = 0; j < Nx; j++){
                        cosh_kDZ[j] = cosh(k[j] * (zoutIP[zz] + depth[0]));
                        sinh_kDZ[j] = sinh(k[j] * (zoutIP[zz] + depth[0]));
                    }//endfor x(j)
                } //endif
                else {
                    for (int j = 0; j < Nx; j++){				
                        gam_DZmin   = gsl_spline_eval(splineDZm0, DZ[j], accDZm0);
                        gam_DZplus  = gsl_spline_eval(splineDZp0, DZ[j], accDZp0);	
                        if (DZ[j] == DZmin) {gam_DZmin = 1; gam_DZplus = 0;}
                        if (DZ[j] == DZplus) {gam_DZmin = 0; gam_DZplus = 1;}						
                        cosh_kDZ[j] = gam_DZmin * cosh(k[j] * DZmin) + gam_DZplus * cosh(k[j] * DZplus);
                        sinh_kDZ[j] = gam_DZmin * sinh(k[j] * DZmin) + gam_DZplus * sinh(k[j] * DZplus); 
                    }//endfor		
                }//endelse
            }//endelseif compute cosh sinh kDZ
			            
            for (int j = 0; j < Nx; j++){				
                ratio   			= cosh_kDZ[j] / cosh_kH[j];				
                factor_nonifft[j].Re 		= phihat_i[j].Re * ratio * (-k[j] * tanh_kH[j]);
                factor_nonifft[j].Im 		= phihat_i[j].Im * ratio * (-k[j] * tanh_kH[j]);
                tempdtPhi_nonifft[j].Re 	= dtphi_hat[j].Re * ratio;	
                tempdtPhi_nonifft[j].Im		= dtphi_hat[j].Im * ratio;
                tempdzPhi_nonifft[j].Re   	= phihat_i[j].Re * k[j] * sinh_kDZ[j] / cosh_kH[j];			
                tempdzPhi_nonifft[j].Im         = phihat_i[j].Im * k[j] * sinh_kDZ[j] / cosh_kH[j];
                tempdxPhi_nonifft[j].Re	  	= -k[j] * phihat_i[j].Im * ratio;	
                tempdxPhi_nonifft[j].Im	  	= k[j] * phihat_i[j].Re * ratio;							
            }//endfor x(j)
            
            ifft_complex_real_nx(factor_nonifft, factor);
            ifft_complex_real_nx(tempdtPhi_nonifft, tempdtPhi_ifft);
            ifft_complex_real_nx(tempdzPhi_nonifft, tempdzPhi);
            ifft_complex_real_nx(tempdxPhi_nonifft, tempdxPhi_ifft);	
			            
            for (int j = 0; j < Nx; j++) {
                if ((zoutIP[zz] > eta[j]) || (zoutIP[zz] < bathy[j])) {
                    tempdtPhi[j] = NAN;
                    tempdzPhi[j] = NAN;
                    tempdxPhi[j] = NAN;			
                } 
                else { 
                    tempdtPhi[j] = tempdtPhi_ifft[j] + dteta[j] * factor[j];
                    tempdxPhi[j] = tempdxPhi_ifft[j] + grad_H[j] * factor[j] + grad_D[j] * tempdzPhi[j];	
                }											
            }//end for j
            int xo = 0;
            while (x[indx1 + xo] <= xoutIP[NxoutIP - 1]) {
                dtPhi[xo][zz] = tempdtPhi[indxIP1 + xo];
                dzPhi[xo][zz] = tempdzPhi[indxIP1 + xo];
                dxPhi[xo][zz] = tempdxPhi[indxIP1 + xo];
                Ptot[xo][zz]  = -dtPhi[xo][zz] - g * zoutIP[zz] - 0.5 * (pow(dxPhi[xo][zz], 2) + 2 * dzPhi[xo][zz]);
                xo 	      = xo + 1;
            }//end xo	
        }//end zz
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
    for (int i = 0; i < NxoutIP; i++) {
        for (int j = 0; j < NzoutIP; j++) {
            dt_dxPHI[i][j] = (dxPhi[i][j] - Velx[i][j]) / (t - tprev);
            dt_dzPHI[i][j] = (dzPhi[i][j] - Velz[i][j]) / (t - tprev);		
        }		
    }	

    // save elevation data at overlay zone
    double* eta_cut = new double[NxoutIP];
    for (int ix = 0; ix < NxoutIP; ix++){
        eta_cut[ix] = eta[indxIP1 + ix];
    }

    #pragma omp parallel sections num_threads(4)
    {
        #pragma omp section
        fwrite(eta_cut, sizeof(double), NxoutIP, dateta);//save elevation at overlay	
        
        #pragma omp section				
        fwrite(Ptot, sizeof(double), NxoutIP * NzoutIP, datP);//save pressure			
        
        #pragma omp section 		
        {
            fwrite(dxPhi, sizeof(double), NxoutIP * NzoutIP, datvel);//save velocity
            fwrite(dzPhi, sizeof(double), NxoutIP * NzoutIP, datvel);
        }
        
        #pragma omp section 
        {
            fwrite(dt_dxPHI, sizeof(double), NxoutIP * NzoutIP, datacc);//save acceleration
            fwrite(dt_dzPHI, sizeof(double), NxoutIP * NzoutIP, datacc);
        }
    }
    delete[] eta_cut;
}
