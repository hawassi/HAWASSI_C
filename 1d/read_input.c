void read_input_file(const char* argv1, const char* argv2) {
    struct stat st = {0};
    if (stat(argv2, &st) == -1) {
        mkdir(argv2, 0700);
    }
  
    char params[128], line[128];
    sprintf(params, "%s/input.txt", argv1);
    FILE* fparams = fopen(params, "r");

    if (fparams == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    
    fgets(line, sizeof(line), fparams);	
    fscanf(fparams, "%*s %*s %s \n", wavename);
    fscanf(fparams, "%*s %*s %s \n", initial);
    rewind(fparams);
    
    for (int i = 0; i < 18; i++){
        fgets(line, sizeof(line), fparams);
    }
    fscanf(fparams, "%*s %*s %s \n %*s %*s %lf \n %*s %*s %lf \n", bath, &depthflat, &adjcoef); 
    rewind(fparams);

    for (int i = 0; i < 23; i++){
	    fgets(line, sizeof(line), fparams);
    }
    fscanf(fparams, "%*s %*s %s \n", dynmodel);
    fscanf(fparams, "%*s %*s %s \n", breaking);
    fscanf(fparams, "%*s %*s %s \n", friction);
    fscanf(fparams, "%*s %*s %s \n", influxing);
    fscanf(fparams, "%*s %*s %s \n", propagation);	
    fscanf(fparams, "%*s %*s %s \n", kinematic);
    fscanf(fparams, "%*s %*s %d \n", &ngauge);
    fscanf(fparams, "%*s %*s %d \n", &npartition);
    
    fgets(line, sizeof(line), fparams);
    fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf %lf \n", &parKBC, &parTC, &parTchar);
    
    fgets(line, sizeof(line), fparams);
    fscanf(fparams, "%*s %*s %lf \n %*s %*s %lf \n", &xinterv1, &xinterv2);
    fscanf(fparams, "%*s %*s %d \n %*s %*s %d \n", &nsurf, &ndeep);
    
    fgets(line, sizeof(line), fparams);
    fscanf(fparams, "%*s %*s %s \n", cut_temp);
    fscanf(fparams, "%*s %*s %s \n", mid_temp);
    fscanf(fparams, "%*s %*s %s \n", min_z);
    fscanf(fparams, "%*s %*s %s \n", max_z);
    fscanf(fparams, "%*s %*s %s %lf", rampbool, &nTramp);
    fscanf(fparams, "%*s %*s %s \n", wall);
    fscanf(fparams, "%*s %*s %d %lf %lf\n", &walltype, &Xwall, &refl_coef);
    
    fgets(line, sizeof(line), fparams);
    fscanf(fparams, "%*s %*s %s %s\n %*s %*s %lf %lf\n", coupling, dir_force, &x_cpl[0], &x_cpl[1]);
    fscanf(fparams, "%*s %*s %lf\n", &alpha_corr);
    
    fclose(fparams);
}

void read_initial_file(const char* argv1, char* init) {
    FILE* finit;
    sprintf(init, "%s/initial.dat", argv1);	
    finit 		= fopen(init, "r");
    if(finit == NULL){
        printf("Can not open initial.dat file!!\n");
        exit(0);
    }
    for (int j = 0; j < Nx; j++){
        fscanf(finit, "%lf %lf\n", &temp_wave[j], &temp_u[j]);
    }
    fclose(finit);
}

void read_breaking_file(const char* argv2, char* data_br) {
    sprintf(data_br, "%s/breaking.dat", argv2);
    fbreak = fopen(data_br,"w");
}

void read_influx_file(const char* argv2, char* ff) {
    strcpy(ff, argv2);
    strcat(ff, "/influx.dat");
    FILE* finf = fopen(ff, "w");
    for (int j = 0; j < n; j++){
        fprintf(finf, "%g %g\n", ttrace[j], insig[j]);
    }
    fclose(finf); 
}

void read_bath_file(const char* argv1, char* input_bath) {
    FILE* fbath;
    sprintf(input_bath, "%s/bath.dat", argv1);	
    fbath	= fopen(input_bath, "r");
    
    if(fbath == NULL){
        printf("Can not open bath.dat file!!\n");
        exit(0);
    }

    for (int j = 0; j < Nx; j++){
        fscanf(fbath, "%lf \n", &depth[j]);
    }
    
    fclose(fbath);
}
