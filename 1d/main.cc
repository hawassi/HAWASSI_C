#include "header.h"
#include "rhs.c"
#include "functions.c"
#include "influx_lib.c"
#include "x_space.c"
#include "plotting.c"
#include "breaking.c"
#include "kinematic_modul.c"
#include "coupling.c"
#include "read_input.c"

int main(int argc, char* argv[])
{
    //Checking argument
    if (argc != 4) {
        fprintf(stderr, "Usage (4 argv): Hawassi working_dir outname (integer)\n");
        return 1;
    }

    //THIS IS FOR CHECKING LICENSE DATE and MACHINES
    printf("HAWASSI v150819 with coupling module\n");
    
    //automatic default
    fftw_init_threads();
    
    // variable initialization
    argv1 = argv[1];
    strcpy(argv2, argv[1]);
    strcat(argv2, argv[2]);
    strcat(argv2, "/");
    strcpy(dispersion, "exact");
    strcpy(wind, "no");
    
    runup = 1;
    mid[0] = 3;
    mid[1] = 0.5;
    
    ncore = atoi(argv[3]);
    ncore2= ncore;
    omp_set_num_threads(ncore);

    //initialization reading data
    printf("Reading input parameters and files...\n");
    read_input_file(argv1, argv2);

    //reading input and get Nt
    if (strcmp(rampbool, "no") == 0){
        nTramp = 0;
    }
    wave = influxlib(&n, argv[1]);

    //checking prime factor of Nx
    primeNx 	= primeFactors(Nx);
    if (primeNx == 1){
        printf("Number of points in x axis has a largest prime factor>13.\nPlease adjust dx!\n");
        return 0;
    }
    
    ttrace     = new double[n];
    for (int j = 0; j < n; j++) {
        ttrace[j] = wave[j].t;
    }
    dt = ttrace[1] - ttrace[0];

    //load and define the bathymetry, friction coefficient, and gauges
    depth   = new double[Nx];
    Cf 	    = new double[Nx];
    Cw 	    = new double[n * Nx];
    Rf_hat  = new complex[Nx]; set_init_complex(Rf_hat, Nx);
    Rw_hat  = new complex[Nx]; set_init_complex(Rw_hat, Nx);	
    
    define_bath(argv[1]);	
    define_friction(argv[1]);	
    define_gauges(argv[1]);

    //get influx and spectra and other needed parameters
    set_memory();	
    if (strcmp(wavename, "zero") != 0){
        set_influx(wave, n, &nu_peak, ww, spinflux); 
        kp 		= invers_omega(nu_peak, depthinf);
        Tp_peak		= 2 * Pi / nu_peak;
        if (strcmp(wavename, "jonswap") == 0){
            nu_peak = 2 * Pi / Tp;
            kp 	= invers_omega(nu_peak, depthinf);
        }
        lambda_p  = 2 * Pi / kp;
        avrg_wave = mean(wave, n);
        Hs_in	  = 4 * sqrt(variance(wave, n, avrg_wave));
        MaxEtaInit= (Hs_in / 2) * 0.8; 	
        modified_influx(wave, n, &nu_peak, ww, spinflux, insig, insig_skew, kom, GVom);	
    }
    else {
        Hs_in      = 0;
        MaxEtaInit = 0;
    }
    
    //generation set-up
    x_space(x, k, dampchar, gamX, spatgamX, Ug, gamX_skew);
    
    gamWall_hat 	= new complex[Nx];
    gamWall_skew_hat 	= new complex[Nx];
    define_wall(gamWall_hat, gamWall_skew_hat);
    
    if (strcmp(initial, "user_defined") == 0){
        temp_wave = new double[Nx];
        temp_u    = new double[Nx];
        char  init[128];
        read_initial_file(argv1, init);
        
        complex* sp_tempwave   = new complex[Nx];
        double*  absp_tempwave = new double[Nx];
        double   spx_max;
        int      idx_spxmax;
        
        fft_real_complex(temp_wave, sp_tempwave, Nx);
        absolute_complex(absp_tempwave, sp_tempwave, Nx);	
        search_max(absp_tempwave, Nx, &spx_max, &idx_spxmax);
        
        kp 	       = fabs(k[idx_spxmax]);
        if (idx_spxmax == 0){
            kp 	       = fabs(k[1]);
        }
        
        double avrgD     = mean_array(depth, Nx);
        nu_peak 	 = omega_i(kp, avrgD);
        double max_twave = search_v(Nx, temp_wave, MAXIMUM);
        double Hs_inhere = 2 * max_twave;
        
        if (Hs_inhere > Hs_in){
            Hs_in = Hs_inhere;
        }
        MaxEtaInit  = (Hs_in / 2) * 0.8; 
	
        delete[] sp_tempwave;
        delete[] absp_tempwave;
    }   

    //IF breaking: initiate crestbreak
    Rb_hat 	= new complex[Nx]; set_init_complex(Rb_hat, Nx);	
    if (strcmp(breaking, "yes") == 0){	
        parambreak(x, Xinflux, parKBC, parTC, parTchar);

        B	        = (double*) malloc(sizeof(double) * Nx);
        CrestBreak      = (CB*) malloc(sizeof(CB) * (Nx / 4));
        CrestBreakNew 	= (CB*) malloc(sizeof(CB) * (Nx / 4));
        CrestBreakPrev 	= (CB*) malloc(sizeof(CB) * (Nx / 4));
        Break_nodes 	= (int*) malloc(sizeof(int) * (Nx / 4));
        ID_Node_end 	= (int*) malloc(sizeof(int) * (Nx / 4));
        dataBreak_nodes = new double[Nx];
        
        set_init_CB(CrestBreak, (Nx / 4));
        set_init_CB(CrestBreakNew, (Nx / 4));
        set_init_CB(CrestBreakPrev, (Nx / 4));
        
        char data_br[128];
        read_breaking_file(argv2, data_br);
    }   

    //coupling module
    if (strcmp(coupling, "yes") == 0){
        //estimate wave input in OZ
        nu_peak   = 2 * Pi / Tp;
        kp 	  = invers_omega(nu_peak, depthflat);
        lambda_p  = 2 * Pi / kp;

        //opening forcing data and interpolating
        int     Nxdata, Ntdata;
        double* xdata;
        double* tdata;
        read_constant(&Nxdata, &Ntdata, xdata, tdata);

        etadata = (double*) malloc(Nxdata * Ntdata * sizeof(double));
        udata   = (double*) malloc(Nxdata * Ntdata * sizeof(double));
        read_eta(Nxdata, Ntdata, etadata);

        alpha_couple = alpha_corr * lambda_p / (x_cpl[1] - x_cpl[0]) * sqrt(9.81 * kp * tanh(kp * depthflat)) / kp;
        idx_cpl[0] = int((x_cpl[0] - x[0]) / dx);
        idx_cpl[1] = int((x_cpl[1] - x[0]) / dx);
        
        int n_idx = idx_cpl[1] - idx_cpl[0] + 1;
        printf("Coupling zone (%g,%g) at index (%d,%d)\n", x_cpl[0], x_cpl[1], idx_cpl[0], idx_cpl[1]);
        printf("Alpha coupling %g\n", alpha_couple);

        Chi_couple = new double[Nx];
        getChi(idx_cpl, Nx);
        update_damp(idx_cpl[1], Nx, dampchar, Chi_couple); 
        
        delete[] xdata;
        delete[] tdata;

    }

    //additional nonlinear adjustment if wall is exist
    if (strcmp(wall, "yes") == 0){
        double *ChiAdj_wall;
        double *ChiAdjL, *ChiAdjR;
        wall_setup(ChiAdj_wall, ChiAdjL, ChiAdjR); 
    }	

    //setting operator
    size_t Ndim = 4 * Nx;
    
    read_kinematics();  
    define_operator(Ndim, nu_peak);   
    alias(aal, k, cutfracwn); //problem found

    if (strcmp(wavename, "zero") != 0 ){
        //information of influx signal
        printf("Peak wave length = %.2f\n", lambda_p);
        printf("Phase speed = %.2f\n", Cpeak);
        printf("Significant wave height = %.2f\n", Hs_in);
        
        plot_initial_config(wave, ww, n, x, spatgamX, dampchar, k, nu_peak, argv2, ChiAdj);
        
        char ff[256];
        read_influx_file(argv2, ff); 
    }
    
    //make global plan for FFT
    fftw_complex* inp_f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    fftw_complex* out_f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    fftw_complex* inp_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    fftw_complex* out_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    plan_ft_nx  = fftw_plan_dft_1d(Nx, inp_f, out_f, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    plan_ift_nx = fftw_plan_dft_1d(Nx, inp_b, out_b, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);   

    // time evolution
    printf("Calculating time evolution...\n"); 
    gettimeofday(&start, NULL);
    check_runup_dynmodel(runup, dynmodel, Ndim); 
    gettimeofday(&end, NULL);	
    delete[] L;
    
    rel = (((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6)/(tim[n-1] - tim[0]);  

    // saving all the output 
    printf("Results are in the directory %s\n", argv2); 
    save_all_output(argv2, argv1);

    //destroy plan of global FFT
    fftw_destroy_plan(plan_ft_nx);
    fftw_destroy_plan(plan_ift_nx);
    fftw_destroy_plan(plan_ift_c2r);
    fftw_destroy_plan(plan_ft_r2c);
    fftw_free(inp_f);fftw_free(out_f);
    fftw_free(inp_b);fftw_free(out_b);
    fftw_cleanup_threads();

    free_memory();
    delete[] insig; 
    delete[] insig_hat; 
    delete[] insig_skew; 
    delete[] Chi_couple;

    return 0;
}
