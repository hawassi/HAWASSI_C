void plot_particle(marin* wave, double* ww, FILE* gp, int n, double const* spgam, double nu_peak, double* ChiAdj)
{
    double* insighat  = new double[n];
    double* insighat2 = new double[n];
    abs_fft_signal(insighat, wave, n);
    
    for (int j = 0; j < n; j++){
        insighat2[j] = pow(insighat[j], 2);
    }

    double varsig = variance(wave, n, avrg_wave);
    double intsp  = trapz(ww, insighat2, n / 2);
    double factor = varsig / intsp;

    fprintf(gp, "set multiplot\n");
    fprintf(gp, "set origin 0,0.5\n");
    fprintf(gp, "set size 0.5,0.5\n");
    fprintf(gp, "set autoscale fix\n");
    fprintf(gp, "set title 'Influx signal'\n");
    fprintf(gp, "set xlabel 'time[s]'\n");
    fprintf(gp, "set xtics font ', 8' \n");
    fprintf(gp, "set ytics font ', 8' \n");
    fprintf(gp, "plot \"-\" w l lt 3 lw 2\n");

    for(int i = 0; i < n; i++){
        fprintf(gp, "%g %g\n", wave[i].t, wave[i].eta);
    }
    fprintf(gp, "e\n"); 
    fprintf(gp, "set origin 0,0\n");
    fprintf(gp, "set size 0.5,0.5\n");
    fprintf(gp, "set xrange [%.2f:%.2f]\n", 0.0, 4 * nu_peak);
    fprintf(gp, "set title 'Variance density spectrum'\n");
    fprintf(gp, "set xlabel 'omega[rad/s]'\n");
    fprintf(gp, "plot \"-\" w l lt 3 lw 2\n");

    for(int i = 0; i < n / 2; i++){
        fprintf(gp, "%g %g\n", ww[i], factor * insighat2[i]);
    }
    fprintf(gp, "e\n");
    fprintf(gp, "set origin 0.5,0.5\n");
    fprintf(gp, "set size 0.5,0.5\n");	
    fprintf(gp, "set xrange [%lf:%lf]\n", x[0], x[Nx - 1]);
    
    if (strcmp(bath, "user_defined") == 0){
        fprintf(gp, "set yrange [-1:1]\n");
    }
    
    fprintf(gp, "set title 'Xspace'\n");
    fprintf(gp, "set xlabel 'x[m]'\n");
    fprintf(gp, "plot \"-\" w l lc 2 lw 2,\"-\" w l lc 3 lw 2,\"-\" w l lc 4 lw 2\n");
    
    for(int i = 0; i < Nx; i++){
        fprintf(gp, "%g %g\n", x[i], spgam[i]);
    }
    fprintf(gp, "e\n");
    
    for(int i = 0; i < Nx; i++){
        fprintf(gp, "%g %g\n", x[i], 1 - dampchar[i]);
    }
    fprintf(gp, "e\n");
    for(int i = 0; i < Nx; i++){
	    fprintf(gp, "%g %g\n", x[i], ChiAdj[i]);
    }
    fprintf(gp, "e\n");
    if (strcmp(bath, "user_defined") == 0){
        fprintf(gp, "plot \"-\" w l lc 5 lw 2\n");
        for(int i = 0; i < Nx; i++){
	    fprintf(gp, "%g %g\n", x[i], -depth[i] / maxD);
        }
        fprintf(gp, "e\n");
    }	
    fprintf(gp, "set origin 0.5,0.0\n");
    fprintf(gp, "set size 0.5,0.5\n");
    fprintf(gp, "set xrange [%lf:%lf]\n", 0.0, invers_omega(4 * nu_peak, depthinf));
    fprintf(gp, "set yrange [0:%lf]\n", 4 * nu_peak);	
    fprintf(gp, "set title 'Dispersion relation'\n");
    fprintf(gp, "set xlabel 'k[1/m]'\n");
    fprintf(gp, "plot \"-\" w l lc 3 lw 2\n");
    for(int i = 0; i < (int(Nx / 2) - 1); i++){	
	fprintf(gp, "%g %g\n", k[i], omega_i(k[i], depthinf));
    }
    fprintf(gp, "e\n");	
    fprintf(gp, "plot \"-\" with points pointtype 10\n");
    fprintf(gp, "%g %g\n", kp, nu_peak);	
    if(strcmp(dynmodel, "HS2") == 0){ 
	    fprintf(gp, "%g %g\n", invers_omega(2 * nu_peak, depthinf), 2 * nu_peak);			
    }	
    if(strcmp(dynmodel, "HS3") == 0){ 
        fprintf(gp, "%g %g\n", invers_omega(2 * nu_peak, depthinf), 2 * nu_peak);	
        fprintf(gp, "%g %g\n", invers_omega(3 * nu_peak, depthinf), 3 * nu_peak);		
    }
    if(strcmp(dynmodel, "HS4") == 0){ 
        fprintf(gp, "%g %g\n", invers_omega(2 * nu_peak, depthinf), 2 * nu_peak);	
        fprintf(gp, "%g %g\n", invers_omega(3 * nu_peak, depthinf), 3 * nu_peak);	
        fprintf(gp, "%g %g\n", invers_omega(4 * nu_peak, depthinf), 4 * nu_peak);		
    }
    fprintf(gp, "e\n");
    fprintf(gp, "unset multiplot\n");
    fflush(gp);

    delete[] insighat;
    delete[] insighat2;
}

void set_gnuplot_condition(FILE* gp, char* dir)
{
    fprintf(gp, "unset key\n");
    fprintf(gp, "unset mouse\n");
    fprintf(gp, "set term pngcairo size 1000,600\n");
    fprintf(gp, "set out '%s/setup.png'\n", dir);
}

void plot_initial_config(marin* wave, double* ww, int n, double* x, double* spatgamX, double* dampchar, double* k, double nu_peak, char* dir, double* ChiAdj)
{
    FILE* gp  = popen("gnuplot > /dev/null 2>&1", "w");
    set_gnuplot_condition(gp, dir);
    plot_particle(wave, ww, gp, n, spatgamX, nu_peak, ChiAdj);//plotting the influx signal 
    printf("Initialization has done.\n");
    printf("You can check the set-up in %ssetup.png\n", dir);
}

