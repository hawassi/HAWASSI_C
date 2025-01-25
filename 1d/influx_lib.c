marin* influxlib(int *n, char* dir)
{		
    char params[128];
    sprintf(params, "%s/input.txt", argv1);
    FILE* fparams 	= fopen(params, "r");
    char line[256];
    if(fparams == NULL){
        printf("No such an input file!!\n");
        exit(0);
    }
    if ((strcmp(wavename, "harmonic") == 0) || (strcmp(wavename, "jonswap") == 0)){				
        for (int i = 0; i < 5; i++){
	    fgets(line, sizeof(line), fparams);
        }
        fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf\n", &T0, &Tend, &dt);
        fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf\n", &Tp, &Hs, &gamJS);
        rewind(fparams);		
        lenT 	= Tend - T0; 
        Nt 	= (lenT) / dt + 1;
        omgp 	= 2 * Pi / Tp;     //omega peak
    }	
    else if (strcmp(wavename, "zero") == 0) {
        for (int i = 0; i < 5; i++){
	    fgets(line, sizeof(line), fparams);
        }
        fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf\n", &T0, &Tend, &dt);
       	fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf\n", &Tp, &Hs, &gamJS);
        rewind(fparams);
        lenT 	= Tend - T0; 
        Nt 	= (lenT / dt) + 1;
        (*n) 	= Nt;    
        orgwave = new marin[(*n)];
        for(int i = 0; i < (*n); i++){
            orgwave[i].eta = 0;
            orgwave[i].t   = T0 + dt * i;
        }
    }
    //define orgwave
    if (strcmp(wavename, "user_defined") == 0){ 
        char influx[128];
        sprintf(influx, "%s/influx.dat", dir);	
        fpo     = fopen(influx, "r");
        if(fpo == NULL){
            printf("Can not open influx file!!\n");
            exit(0);
        }
        orgwave	= read_data(fpo, n); fclose(fpo);
    }
    else if (strcmp(wavename, "harmonic") == 0){// with ramp
        harmonic_influx(T0, Tend, Tp, dt, 0.5 * Hs, n);
    }
    else if (strcmp(wavename, "jonswap") == 0){//with ramp 		
        double* omsig = (double*) malloc(Nt * sizeof(double));
        double* halfomsig = (double*) malloc(Nt / 2 * sizeof(double)); 
        freqspace(omsig, lenT, Nt);
        for (int i = 1; i  < Nt / 2; i++){
	    halfomsig[i]   = omsig[i];
        }
        jonswap_influx(T0, Tend, Tp, dt, Nt, omsig, omgp, halfomsig, gamJS, Hs, n);
        free(omsig);
        free(halfomsig);	
    }

    if(tcoarse>1) { (*n) = int((*n) / tcoarse);}
    marin* wave = new marin[(*n)];
    for (int j = 0; j < (*n); j++) {
        wave[j].t  = orgwave[j].t;
        wave[j].eta= orgwave[j].eta;
    }
    delete[] orgwave;	

    //get others parameters
    for (int i = 0; i < 13; i++){
        fgets(line, sizeof(line), fparams);
    }
    fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf \n %*s %*s %lf\n", &Xstart, &Xend, &Xinflux);
    fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf %lf\n ", &dx, &dampL, &dampR);
    fclose(fparams);
    
    if (strcmp(dynmodel, "HS1") == 0){
	    cutfracwn = 2;
    }
    else if (strcmp(dynmodel, "HS2") == 0){
	    cutfracwn = 4;
    }
    else if (strcmp(dynmodel, "HS3") == 0){
	    cutfracwn = 6;
    }
    else if (strcmp(dynmodel, "HS4") == 0){ // by pattern
	    cutfracwn = 8;
    }

    //Domain definition and spatial discretization
    length  	= (Xend - Xstart);
    Nx		= (length / dx) + 1;
    if (Nx % 2 != 0) { 
        Nx = Nx - 1;
    }
    
    x	      = new double[Nx];		
    dx 	      = length / (Nx - 1);	
    for(int i = 0; i < Nx; i++){
	x[i]  = Xstart + dx * i;
    }
    
    idx_inf 	= closest(x, Xinflux, Nx);
    Xinflux 	= x[idx_inf];
    return wave;
}
 
void harmonic_influx(double T0, double Tend, double Tp, double dt, double ampli, int* n)
{
    Vmarin 	p_tmp; 
    double 	ramp, ramp1, ramp2;
    double 	t = T0;
    int 	i = 0;
    while(t <= Tend){
        marin 	        wave_temp;
        t		= T0 + i  * dt;
        wave_temp.t	= t;
        ramp1		= ((t - 4 * Tp) < 0? -1 : 1);
        ramp2		= (1 - cos(t * Pi / 4 / Tp)) / 2;
        ramp		= ramp1 > ramp2? ramp1 : ramp2;
        wave_temp.eta   = ampli * ramp * (sin(t * 2 * Pi / Tp));
        p_tmp.push_back(wave_temp);//ganti ke allocated memory biar lebih efisien
        i++;
    }
    
    (*n)    = p_tmp.size();    
    orgwave = new marin[(*n)];
    for(int i = 0; i < (*n); i++){
	orgwave[i] = p_tmp[i];
    }
    p_tmp.clear();
 }

double* fun_jonswap(double* omg, double omgp, double* halfomsig, double Tp, double gamJS, double Hs, int Nt)
{
    int     Nomg 	= Nt/2;                   
    double* sigma 	= (double*) malloc (Nomg * sizeof(double));

    for (int i = 0; i < Nomg; i++){
        if (omg[i] <= omgp){
	    sigma[i] 	= 0.07;
        }
        else {
	    sigma[i]	= 0.09;
        }
    }
    
    double  f, gampangkat;
    double* JS	    = (double*) malloc(Nomg * sizeof(double));
    double* JPM     = (double*) malloc(Nomg * sizeof(double));
    double* JS_NOR  = (double*) malloc(Nomg * sizeof(double));

    //double beta     = -5/4;
    double alpha    = 0.0081;
    double f_p      = omgp/(2 * Pi);

    for (int i = 1; i < Nomg; i++){
        f 	    = omg[i] / (2 * Pi);
        gampangkat  = exp(-0.5 * pow(((f / f_p - 1) / sigma[i]), 2));
        JPM[i]      = alpha * pow(g, 2) * pow((2 * Pi), (-4)) * pow(f, (-5)) * exp(-(5 / 4) * pow((f / f_p),(-4))); //Pierson - Moskowitz
        JS[i]       = alpha * pow(g, 2) * pow((2 * Pi), (-4)) * pow(f, (-5)) * exp(-(5 / 4) * pow((f / f_p),(-4))) * pow(gamJS, (gampangkat)); //JONSWAP
    }

    // calculating Hs from non-normalized JONSWAP spectrum 
    double m0;
    m0   = trapz(omg, JS, Nomg);  

    // Now normalize the spectrum based on definition of zeroth-order moment 
    // and its relation with Hs : mo = int(omg,JS(omega)) and Hs = 4*sqrt(mo); so:
    for (int i=0; i < Nomg; i++){
        JS_NOR[i]  = (pow(Hs, 2) / 16) * JS[i] / (m0);
    }
    free(sigma);
    free(JS);
    free(JPM);			

    return JS_NOR;
}

void jonswap_influx(double T0, double Tend, double Tp, double dt, int Nt, double* omsig, double omgp, double* halfomsig, double gamJS, double Hs, int* n)
{
    Vmarin  p_tmp; 
    double  tin	       = T0;
    double* JS_sp      = (double*) malloc(Nt / 2 * sizeof(double));
    double* phase      = (double*) malloc(Nt / 2 * sizeof(double));
    double* rawsig     = (double*) malloc(Nt * sizeof(double));
    double* insigh     = (double*) malloc(Nt * sizeof(double));
    complex* rawsighat = (complex*) fftw_malloc(sizeof(complex) * Nt);

    JS_sp = fun_jonswap(omsig, omgp, halfomsig, Tp, gamJS, Hs, Nt);	

    srand(time(0)); 	
    for (int i = 0; i < (Nt / 2); i++){
	    phase[i + 1]             = 2 * Pi * ((double) rand() / (RAND_MAX));
	    rawsighat[i + 1].Re      = (sqrt(JS_sp[i]) * cos(phase[i])) / dt * 2 * pow(Pi, 2);
	    rawsighat[i + 1].Im      = -(sqrt(JS_sp[i]) * sin(phase[i])) / dt * 2 * pow(Pi, 2);
	    rawsighat[Nt - 1 - i].Re = rawsighat[i + 1].Re;		
	    rawsighat[Nt - 1 - i].Im = -rawsighat[i + 1].Im;
    }
    rawsighat[0].Re = 0;		
    rawsighat[0].Im = 0;
    ifft_complex_real(rawsighat, rawsig, Nt);	

    double* ramp    = new double[Nt];
    double* endramp = new double[Nt];
    double  A, B1, B2;

    for (int j = 0; j < Nt; j++){
        tin = T0 + j * dt;
        A   = ((tin - T0) - nTramp * Tp);
        B1  = A / fabs(A);
        B2  = (1 - cos((tin - T0) * Pi / nTramp / Tp)) / 2;
        if (B1 > B2) {
          ramp[j] = B1;
        } 
        else {
          ramp[j] = B2;
        }			
    }

    double* insigt 	= new double[Nt];
    for (int j = 0; j < Nt; j++) {
        endramp[j] = ramp[Nt - 1 - j];
        insigt[j]  = ramp[j] * rawsig[j] * endramp[j];
    }		
    delete[] ramp; delete[] endramp;			
    free(phase);	

    double var_sig = var_array(insigt, Nt, mean_array(insigt, Nt));	
    for (int j = 0; j < Nt; j++) {
	insigh[j]  = insigt[j] * Hs / (4 * sqrt(var_sig));
    }
    delete[] insigt;

    for (int i = 0; i < Nt; i++){
        marin wave_temp;
        tin		= T0 + i * dt;
        wave_temp.t	= tin;
        wave_temp.eta 	= insigh[i];
        p_tmp.push_back(wave_temp);
    }
    (*n) 	= p_tmp.size();    
    orgwave = new marin[(*n)];
    for(int i = 0; i < (*n); i++){
        orgwave[i] = p_tmp[i];
    }
    p_tmp.clear();	

    free(rawsig);
    free(insigh);
    free(rawsighat);	
}
