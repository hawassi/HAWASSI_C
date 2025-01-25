void getChi(int* idx_cpl, int Nx)
{
    int nd = idx_cpl[1] - idx_cpl[0] + 1;
    int nc = nd / 3;
    int ic = idx_cpl[0] + nc;

    for (int j = 0; j < Nx; j++){
        if (j < idx_cpl[0]){
	    Chi_couple[j] = 0;
        } else if (j > (idx_cpl[1])){
	    Chi_couple[j] = 0;
        } else if ((j >= idx_cpl[0]) && (j <= ic)){
	    Chi_couple[j] = 1 - (1.0 + cos(((double) j - (double) idx_cpl[0]) / (double) (nc) * Pi)) / 2.0;
        }
        else {
	    Chi_couple[j] = 1;
        }
    }
    
    FILE* gp1 = popen("gnuplot -persistent", "w");
    fprintf(gp1, "set term x11\n");
    fprintf(gp1, "plot \"-\" w points pointtype 10\n");
    
    for (int i = idx_cpl[0] - 10; i < (idx_cpl[1] + 10); i++) {
        fprintf(gp1, "%lf %lf\n", x[i], Chi_couple[i]);
    }
    
    fprintf(gp1, "e\n");
    
    fflush(gp1);
}

void coupling_term(double* eta, double* u, complex* couple_etahat, complex* couple_uhat, double tnow, double* etadata, double* udata)
{	
    int idt = tnow / dt;

    double*  temp_eta   = new double[Nx];
    double*  temp_u     = new double[Nx];
    double*  couple_eta = new double[Nx];
    double*  couple_u   = new double[Nx];

    for(int j = 0; j < Nx; j++){
        couple_eta[j]   = 0;
        couple_u[j]     = 0;
    }
          
    for(int j = 0; j < Nx; j++){
        couple_eta[j] = etadata[idt * Nx + j];
    }

    //coupling term in real space
    for(int j = 0; j < Nx; j++){
        temp_eta[j] 	= (alpha_couple * Chi_couple[j]) * (couple_eta[j] - eta[j]);
        temp_u[j] 	= (0 * Chi_couple[j]) * (couple_u[j] - u[j]);
    }
    
    complex* temp_etahat = new complex[Nx];
    complex* temp_uhat   = new complex[Nx];
    
    fft_real_complex(temp_eta, temp_etahat, Nx);
    fft_real_complex(temp_u, temp_uhat, Nx);
	
    #pragma omp parallel for
    for(int j = 0; j < Nx; j++){
        couple_etahat[j].Re = temp_etahat[j].Re;
        couple_uhat[j].Re   = temp_uhat[j].Re;
        couple_etahat[j].Im = temp_etahat[j].Im;
        couple_uhat[j].Im   = temp_uhat[j].Im;
    }

    delete[] temp_eta;
    delete[] temp_u;
    delete[] couple_eta;
    delete[] couple_u;
    delete[] temp_etahat;
    delete[] temp_uhat;
}

void read_constant(int *Nxdata, int *Ntdata, double* xdata, double* tdata)
{
    double N, M;

    //read parameter
    char* name_constant = new char[256];
    sprintf(name_constant, "%s/CFD_01.dat", dir_force);
    FILE* fc = fopen(name_constant, "rb");
    if (fc == NULL){
        printf("File constant unable to open...\n");
    }
    else {
        fc = fopen(name_constant, "rb");
        fread(&N, sizeof(double), 1, fc);
        
        *Nxdata = (int) N;
        xdata = new double[*Nxdata];
        fread(xdata, sizeof(double), *Nxdata, fc);
        fread(&M, sizeof(double), 1, fc);
        
        *Ntdata = (int) M;
        tdata = new double[*Ntdata];
        fread(tdata, sizeof(double), *Ntdata, fc);
        
        fclose(fc);
    }
    
    delete[] name_constant;
    
    if (*Nxdata != Nx){
        printf("Nx = %d, Ninput = %d\n", Nx, *Nxdata);
        printf("Number of NX should be equal!!\n");
        exit(0);
    }
}

void read_eta(int Nxdata, int Ntdata, double* etadata)
{
    //open file elevation
    char* filename = new char[256];
    sprintf(filename, "%s/CFD_eta_01.dat", dir_force);
    FILE* feta = fopen(filename, "rb");
    if (feta == NULL){
        printf("File elevation unable to open...\n"); 
    }
    else{
        printf("File elevation opened...\n");
        fread(etadata, sizeof(double), Ntdata * Nxdata, feta);
    }
    
    fclose(feta);
    
    delete[] filename;
}

void read_u(int Nxdata, int Ntdata, double* udata)
{
    //open file velocity
    char* filename = new char[256];
    sprintf(filename, "%s/Hawassi_u_01.dat", dir_force);
    FILE* fu = fopen(filename, "rb");
    if (fu == NULL){
        printf("File velocity unable to open...\n");  
    }
    else{
        printf("File velocity opened...\n");
        fread(udata, sizeof(double), Ntdata * Nxdata, fu);
    }
    
    fclose(fu);
    
    delete[] filename;
}

