double* cf(double* xi, double Xi, double del, int Nx)
{
    double* cfun = new double[Nx];
    for (int j = 0; j < Nx; j++) {
        if (xi[j] < Xi) {
            cfun[j] = 0;
        }
        else if (xi[j] == Xi) {
            cfun[j] = 0.5 * fmax((xi[j] - (Xi + del)) / fabs(xi[j] - (Xi + del)),(1 - cos((xi[j] - Xi) * Pi / del)) / 2);
        }
        else if (xi[j] > Xi) {
            cfun[j] = fmax((xi[j] - (Xi + del)) / fabs(xi[j] - (Xi + del)),(1 - cos((xi[j] - Xi) * Pi / del)) / 2);
        }
    }	
    return cfun;
}

void x_space(double* x, double* k, double* dampchar, complex* gamX, double* spatgamX, double* Ug, complex* gamX_skew) {
    double* spgam = new double[Nx];
    dampzone(dampchar, x, Xstart + dampL, Xend - dampR, dampL, dampR);

    //define dampcoef
    if (strcmp(wavename, "zero") != 0) {//IF influxing
        Cpeak  = omega_i(kp, depthinf) / kp;
        CpLeft = omega_i(kp, depth[0]) / kp;
        CpRight= omega_i(kp, depth[Nx - 1]) / kp;
    }
    else {
        Cpeak  = sqrt(g * depthinf);   
        CpLeft = sqrt(g * depth[0]); 
        CpRight= sqrt(g * depth[Nx - 1]);   
    }
    
    dampcoef = new double[Nx];
    for (int j = 0; j < (Nx / 2); j++) {
        dampcoef[j]        = 7 * CpLeft / (0.9 * dampL);
        dampcoef[j + Nx/2] = 7 * CpRight / (0.9 * dampR); 
    }	
	
    //frequency space (physical wavenumbers)
    freqspace(k, length, Nx);

    //Influxing: generation method setup;
    double fact = 1 / dx;
    double cp = sqrt(g * depthinf);
    for (int i = 0; i < Nx; i++){
        Ug[i] = group_velocity_i(k[i], depthinf);
    }

    if (strcmp(wavename, "zero") != 0) {
        if (strcmp(influxing, "Point") == 0) {
            double dista = Xinflux - Xstart;
            for (int i = 0; i < Nx; i++) {
                spatgamX[i] = 0;
                gamX[i].Re  = cos(k[i] * (Xinflux - Xstart)) * fact;
                gamX[i].Im  = -sin(k[i] * (Xinflux - Xstart)) * fact;
            }
            spatgamX[idx_inf] = fact;
        }
	else {
	    ifft_real_real(Ug, spgam, Nx);
	    circshift(spgam,int((Xinflux - Xstart) / dx), Nx);
	    
	    if (strcmp(influxing, "AreaShort") == 0) {
                double  Lb = 2 * pow(0.4 * lambda_p, 2);
                double* gaussInf = (double*) malloc(sizeof(double) * Nx);
                int 	ind1 = closest(x, Xinflux - 2 * lambda_p, Nx);
                int 	ind2 = closest(x, Xinflux + 2 * lambda_p, Nx);
                int 	i    = ind1;
                set_init_double(gaussInf, Nx);
                while (i < ind2 + 1) {
                    gaussInf[i] = exp(-pow((x[i] - Xinflux), 2) / Lb);
                    i++;
                }			
                for (int i = 0; i < Nx; i++){
	            spgam[i] = spgam[i] * gaussInf[i];
                }
                free(gaussInf);	
	    }
	    
	    fft_real_complex(spgam, gamX, Nx);
	    scaling_influx(gamX, fact, dx, cp);
	    ifft_complex_real(gamX, spatgamX, Nx);
	    delete[] spgam;		
		
	    if (strcmp(influxing, "AreaShort") == 0) {
	        double*  abs_gamX 	= new double[Nx];	
	        complex* temp 		= new complex[n];	
	        double*  gamXom 	= new double[n];
	        double*  ktemp  	= new double[Nx];
	        
	        for (int i = 0; i < Nx / 2; i++){
	            abs_gamX[i]         = sqrt(pow(gamX[Nx/2 + i].Re, 2) + pow(gamX[Nx/2 + i].Im, 2));
	            abs_gamX[Nx/2 + i]  = sqrt(pow(gamX[i].Re, 2) + pow(gamX[i].Im, 2));
	            ktemp[i]	  	= k[Nx/2 + i];
	            ktemp[Nx/2 + i] 	= k[i];
	        }						
	        
	        interp_all(ktemp, abs_gamX, Nx, &acc_gamXom, &spline_gamXom); 
	        
	        for (int i = 0; i < n; i++){				
	            if (fabs(kom[i]) > fabs(ktemp[Nx - 1])) {
		        gamXom[i] = abs_gamX[Nx - 1];
	            }
	            else {
		        gamXom[i] = gsl_spline_eval(spline_gamXom, kom[i], acc_gamXom);
	            }				
	            temp[i].Re   = insig_hat[i].Re * GVom[i] / gamXom[i] / dx;
	            temp[i].Im   = insig_hat[i].Im * GVom[i] / gamXom[i] / dx;
	        }
	        
	        ifft_complex_real(temp, insig, n);
	        cumtrapz(insig, n, dt, insig_skew);
	        delete[] temp; 
	        delete[] gamXom;
	        delete[] abs_gamX;
	        delete[] ktemp;	
	    }		
	}
	
	for (int i = 0; i < Nx; i++) {//no smoothing in gamX_skew
	    gamX_skew[i].Re = -omega_i(k[i], depthinf) * gamX[i].Im;
	    gamX_skew[i].Im = omega_i(k[i], depthinf) * gamX[i].Re;
	}			        
	
	//prepare source for ODE
	interp_all(ttrace, insig, n, &acc, &spline);//give spline and acc for ODE ppval
	interp_all(ttrace, insig_skew, n, &accskew, &splineskew);//give splineskew and accskew for ODE ppvalskew	
    }	
	
    //nonlinear adjustment for nonlinear model HsF2 and HsF3
    if ((strcmp(dynmodel, "HS1") == 0) || (strcmp(wavename, "zero") == 0)){
        ChiAdj  = new double[Nx];
        for (int j = 0; j < Nx; j++)	{
            ChiAdj[j] = 1;
        }
    }
    else
    {	
        if (strcmp(propagation, "Uni+") == 0) {
            ChiAdj    = cf(x, Xinflux, adjcoef * lambda_p, Nx);
        }
        else if (strcmp(propagation, "Uni-") == 0) {
            double* Chitemp;
            ChiAdj  = new double[Nx];
            Chitemp = cf(x, (Xinflux - adjcoef * lambda_p), adjcoef * lambda_p, Nx);
            for (int j = 0; j < Nx; j++)	{
	        ChiAdj[j] = 1 - Chitemp[j];
            }
            delete[] Chitemp;
        }
        else if (strcmp(propagation, "Bi") == 0) {
            double *ChiAdjL, *ChiAdjR;
            ChiAdj  = new double[Nx];
            ChiAdjL = cf(x, Xinflux, adjcoef * lambda_p, Nx);
            ChiAdjR = cf(x, Xinflux - adjcoef * lambda_p, adjcoef * lambda_p, Nx);
            for (int j = 0; j < Nx; j++)	{
                ChiAdj[j] = 1 - ChiAdjR[j] + ChiAdjL[j];
            }
            delete[] ChiAdjL; 
            delete[] ChiAdjR;
        }
    }
}
