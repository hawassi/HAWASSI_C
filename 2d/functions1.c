int countlines(FILE *fp)
{
	int ch,lines=0;
	while(!feof(fp)){
		ch = fgetc(fp);
		if(ch == '\n') {
			lines++;
		}
	}
    rewind(fp);
    return lines+1; 
}

void CircularShift_left(double* arr,int SIZE)
{
	double temp=arr[0];
	for(int i=0;i<SIZE-1;i++){
		arr[i]=arr[i+1];
	}
	arr[SIZE-1]=temp;
}

void CircularShift_right(double* arr,int SIZE)
{
	double temp=arr[SIZE-1];//put the last element in the temp variable
    for(int i=SIZE-1;i>0;i--){//shift each value to the next value so that a hole can be created in the beginning
		arr[i]=arr[i-1];
	}//shift the values,index 0 is not changed because faida nahi hai
    arr[0]=temp;
}

void circshift(double* array,int shift,int N)
{
	int INDX=fabs(shift);
	if (shift>0){
		for(int i=0;i<INDX;i++) CircularShift_right(array,N); 
	}
	else {
		for(int i=0;i<INDX;i++) CircularShift_left(array,N);
	}
}

double* fun_const_substract_array(double constanta,double* array,int N)
{
	double* arrayN=new double[N];
	for(int i=0;i<N;i++){
		arrayN[i]=constanta-array[i];
	}
	return arrayN;
}

void fun_plot_2d_array(FILE* gp,double* y ,double* x ,double** var, int Ny,int Nx)
{
	fprintf(gp, "set pm3d map\n" );
    fprintf(gp, "set palette defined (0 '#000090',\
									  1 '#000fff',\
									  2 '#0090ff',\
									  3 '#0fffee',\
									  4 '#90ff70',\
									  5 '#ffee00',\
									  6 '#ff7000',\
									  7 '#ee0000',\
									  8 '#7f0000')\n");
	fprintf(gp, "splot \"-\" with pm3d\n" );
    for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
		fprintf(gp,"%g %g %g\n",y[i],x[j],var[i][j]);
		}
		fprintf(gp,"\n");
	}
	fprintf(gp,"e\n");
    fflush(gp);
}
		
void  fun_plot_line(double*x ,double* var, int Nx)//splot_opt set)
{
	FILE* gp = popen("gnuplot ","w");
    //fprintf(gp,"set title %s\n",set.title);
	//fprintf(gp,"set xlabel %s\n",set.xlabel);
	fprintf(gp, "plot \"-\" with line\n" );
 		for(int j=0;j<Nx;j++){
		fprintf(gp,"%g %g\n",x[j],var[j]);
		}
	fprintf(gp,"e\n");
    fflush(gp);
	printf("press enter to continue");	
	getchar();
	pclose(gp);
}

void  fun_plot_point(double*x ,double* var, int Nx)//splot_opt set)
{
	FILE* gp = popen("gnuplot ","w");
    //fprintf(gp,"set title %s\n",set.title);
	//fprintf(gp,"set xlabel %s\n",set.xlabel);
	fprintf(gp, "plot \"-\" with point\n" );
 		for(int j=0;j<Nx;j++){
		fprintf(gp,"%g %g\n",x[j],var[j]);
		}
	fprintf(gp,"e\n");
    fflush(gp);
	printf("press enter to continue");	
	getchar();
	pclose(gp);
}

void  fun_plot_bline(double* x,double* var, int Nx)
{
	FILE* gp = popen("gnuplot ","w");
	fprintf(gp, "plot \"-\" with points\n" );
 		for(int j=0;j<Nx;j++){
		fprintf(gp,"%g %g\n",x[j],var[j]);
		}
	fprintf(gp,"e\n");
    fflush(gp);
	pclose(gp);
}

double fun_min2darray(double** matrix,int Nx,int Ny)
{
	double min=matrix[0][0];
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			min=(min<matrix[i][j]?min:matrix[i][j]);
		}
	}
	return min;	
}

double fun_max2darray(double** matrix,int Nx,int Ny)
{
	double max=matrix[0][0];
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			max=(max>matrix[i][j]?max:matrix[i][j]);
		}
	}
	return max;	
}

void fun_plot_2d_trarray(FILE* gp,double* y ,double* x ,double** var, int Ny,int Nx)
{
	double mi=fun_min2darray(var,Nx,Ny);
	double ma=fun_max2darray(var,Nx,Ny);
	
	fprintf(gp,"set pm3d map\n" );
	if (ma-mi>1E-6){
		fprintf(gp,"set cbrange [%f:%f]\n",mi,ma);
	}
    fprintf(gp, "set palette defined (0 '#000090',\
									  1 '#000fff',\
									  2 '#0090ff',\
									  3 '#0fffee',\
									  4 '#90ff70',\
									  5 '#ffee00',\
									  6 '#ff7000',\
									  7 '#ee0000',\
									  8 '#7f0000')\n");
	fprintf(gp, "splot \"-\" with pm3d\n" );
    for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			fprintf(gp,"%g %g %g\n",x[i],y[j],var[j][i]);
		}
		fprintf(gp,"\n");
	}
	fprintf(gp,"e\n");
    fflush(gp);
}

int fun_sign(double var){
	return (var<0?-1:1);
}

int fun_heav(double var){
	return (var<0?0:1);
}

int calNpoint(double xi, double xf,double dx)
{
	return round((xf-xi)*1000/(dx*1000)+1);
}

double* fun_heaviside(double* x,int Nx, double x0)
{
	double* heav=new double[Nx];
	for(int i=0;i<Nx;i++){
		if(x[i]-x0<0){
			heav[i]=0;
		}
		else {
			heav[i]=1;
		}
	}
	return heav;
}

void set_array_val(double* array, int N, double val)
{
	for (int j=0;j<N;j++){
		array[j] = val;
	}
}

void fun_intervaltovector(double* vectvar,double varI,double varF, int NN)
{
	double dxh = (varF-varI)/(NN-1);
	for(int i=0; i<NN;i++){
		vectvar[i]=varI+(i)*dxh;	
	}
}

void fun_intervaltovector_ds(double* vectvar,double varI,double ds, int NN)
{
	for(int i=0; i<NN;i++){
		vectvar[i]=varI+(i)*ds;	
	}
}

void set_fftw_zero(fftw_complex* array,int Nx,int Ny)
{
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			array[i*Nx+j][0] = 0;
			array[i*Nx+j][1] = 0;
		}
	}
}

void set_init_double(double* B, int n)
{
	for (int j=0;j<n;j++){
		B[j]   = 0;
	}
}

double trapz(double* x, double* array, int N)
{
	double sum = 0.0;
	for ( int i=0; i < N; i++){
        if ( i == 0 || i == N-1 ) // for the first and last elements
            sum += array[i]/2;
        else
            sum += array[i]; // the rest of data
    }   
    return sum*(x[1]-x[0]); // the result
}

void idx_max_matrix(int* idx,double** array, int n, int m)
{
	double temp = array[0][0];

	idx[0] = -1;
	idx[1] = -1;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			if(array[i][j]>temp){
				temp=array[i][j];
				idx[0]=j;
				idx[1]=i;
			}
		}
	}
}

int** declare_2darray_int(int Nx, int Ny){
	int** Ar2d=new int*[Ny];
	for (int i=0;i<Ny;i++){
		Ar2d[i]=new int[Nx];
	}
	return Ar2d;
}

double** declare_2darray(int Nx, int Ny){
	double** Ar2d=new double*[Ny];
	for (int i=0;i<Ny;i++){
		Ar2d[i]=new double[Nx];
	}
	return Ar2d;
}

void free_2darray_int(int** Ar2d,int Nx, int Ny){
	for (int i=0;i<Ny;i++){
		delete[] (Ar2d[i]);
	}
	delete[] Ar2d;
}

void free_2darray(double** Ar2d,int Nx, int Ny){
	for (int i=0;i<Ny;i++){
		delete[] (Ar2d[i]);
	}
	delete[] Ar2d;
}

float** declare_2darray_float(int Nx, int Ny){
	float** Ar2d=new float*[Ny];
	for (int i=0;i<Ny;i++){
		Ar2d[i]=new float[Nx];
	}
	return Ar2d;
}

void free_2darray_float(float** Ar2d,int Nx, int Ny){
	for (int i=0;i<Ny;i++){
		delete[] (Ar2d[i]);
	}
	delete[] Ar2d;
}

complex** declare_2darray_complex(int Nx,int Ny){
	complex** Ar2dC=new complex*[Ny];
	for (int i=0;i<Ny;i++){
		Ar2dC[i]=new complex[Nx];
	}
	return Ar2dC;
}

void free_2darray_complex(complex** Ar2dC,int Nx, int Ny){
	for (int i=0;i<Ny;i++){
		delete[] (Ar2dC[i]);
	}
	delete[] (Ar2dC);
}

double** transpose(double** array,int N,int M)
{
	double** tr = declare_2darray(M,N);
	for (int i=0;i<M;i++){
		for (int j=0;j<N;j++){
			tr[j][i] = array[i][j];
		}
	}
	return tr;
}

double** fun_Heav2d(double* Y,int Ny,double bottom,double top,int Nx)
{
	double** heav = declare_2darray(Ny,Nx);
	double*  heavb;
	double*  heavt;
	
	if (bottom==0) {
		heavb = new double[Ny];
		set_array_val(heavb,Ny,1.0);
	}
	else {
		heavb = fun_heaviside(Y,Ny,(Y[0]+bottom));
	}

	if (top==0) {
		heavt = new double[Ny];
		set_array_val(heavt,Ny,0.0);
	}
	else {
		heavt = fun_heaviside(Y,Ny,(Y[Ny-1]-top));
	}

	for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			heav[i][j] = heavb[j]*(1-heavt[j]); 
		}
	}
	delete[] heavb;
	delete[] heavt;
	
	return heav;
}


void make_1d_sym(double* input_array,int Nx)
{
	int id=1;
	for (int j=Nx-1;j>Nx/2;j--){
		input_array[j] = input_array[id];
		id++;
	}
}

void make_1dcom_sym(complex* input_array,int Nx)
{
	int id=1;
	for (int j=Nx-1;j>Nx/2;j--){
		input_array[j].Re = input_array[id].Re;
		input_array[j].Im = -input_array[id].Im;
		id++;
	}
}

void fft_1d(complex* result,double* input_array,int Nx)
{
	fftw_complex *in	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_complex *out 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_plan plan_ft_h  = fftw_plan_dft_1d(Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	
	for (int j = 0; j < Nx; j++){ 
		in[j][0] = input_array[j]; //real part of input
		in[j][1] = 0; //imajiner part
	}	
	
    fftw_execute(plan_ft_h); 
    
    for (int j = 0; j < Nx; j++){
		result[j].Re = out[j][0];
		result[j].Im = out[j][1];
	}
    
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan_ft_h);	
	
}

void ifft_1d_real_real(double* result,double* input_array,int Nx)
{
	//make it symmetric first
	make_1d_sym(input_array,Nx);
	
	fftw_complex *in	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_complex *out 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_plan plan_ift_h = fftw_plan_dft_1d(Nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	
	for (int j = 0; j < Nx; j++){ 
		in[j][0] = input_array[j]; //real part of input
		in[j][1] = 0; //imajiner part
	}	
    fftw_execute(plan_ift_h); 
    
    for (int j = 0; j < Nx; j++){
		result[j] = out[j][0]/Nx;
	}
    
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan_ift_h);	
}

void ifft_1d_com_re(complex* input_array,double* result,int Nx)
{
	//make it symmetric first
	make_1dcom_sym(input_array,Nx);
	
	fftw_complex *in	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_complex *out 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_plan plan_ift_h = fftw_plan_dft_1d(Nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	
	for (int j = 0; j < Nx; j++){ 
		in[j][0] = input_array[j].Re; //real part of input
		in[j][1] = input_array[j].Im; //imajiner part
	}  
    fftw_execute(plan_ift_h); 
    
    for (int j = 0; j < Nx; j++){
		result[j] = out[j][0]/Nx;
	}
    
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan_ift_h);
}

void complex_2d_sym(complex** input_array,int Nxy,int Nt)
{
	input_array[0][0].Im = 0;
	input_array[0][Nxy/2].Im = 0;
	input_array[Nt/2][0].Im  = 0;
	input_array[Nt/2][Nxy/2].Im = 0;
	
	for (int j=1;j<Nxy/2;j++){
		input_array[0][Nxy-j].Re = input_array[0][j].Re;
		input_array[0][Nxy-j].Im = -input_array[0][j].Im;
		input_array[Nt/2][Nxy-j].Re = input_array[Nt/2][j].Re;
		input_array[Nt/2][Nxy-j].Im = -input_array[Nt/2][j].Im;
	}
	for (int i=1;i<Nt/2;i++){
		input_array[Nt-i][0].Re = input_array[i][0].Re;
		input_array[Nt-i][0].Im = -input_array[i][0].Im;
	}
	
	int iy = 1;
    int it = 1;
    for (int i=1;i<Nt;i++){
		iy=1;
		for (int j=Nxy/2+1;j<Nxy;j++){
			input_array[i][j].Re = input_array[Nt-it][Nxy/2-iy].Re;
			input_array[i][j].Im = -input_array[Nt-it][Nxy/2-iy].Im;
			iy++;
		}
		it++;
	}	
}

void fftwcomplex_2d_sym(fftw_complex* input_array,int Nxy,int Nt)
{
	input_array[0][1] = 0;
	input_array[0*Nxy+(Nxy/2)][1] = 0;
	input_array[Nt/2*Nxy][1]  = 0;
	input_array[Nt/2*Nxy+(Nxy/2)][1] = 0;
	
	for (int j=1;j<Nxy/2;j++){
		input_array[Nxy-j][0] = input_array[j][0];
		input_array[Nxy-j][1] = -input_array[j][1];
		input_array[Nt/2*Nxy+Nxy-j][0] = input_array[Nt/2*Nxy+j][0];
		input_array[Nt/2*Nxy+Nxy-j][1] = -input_array[Nt/2*Nxy+j][1];
	}
	for (int i=1;i<Nt/2;i++){
		input_array[(Nt-i)*Nxy][0] = input_array[i*Nxy][0];
		input_array[(Nt-i)*Nxy][1] = -input_array[i*Nxy][1];
	}
	
	int iy = 1;
    int it = 1;
    for (int i=1;i<Nt;i++){
		iy=1;
		for (int j=Nxy/2+1;j<Nxy;j++){
			input_array[i*Nxy+j][0] = input_array[(Nt-it)*Nxy+Nxy/2-iy][0];
			input_array[i*Nxy+j][1] = -input_array[(Nt-it)*Nxy+Nxy/2-iy][1];
			iy++;
		}
		it++;
	}	
}

void ifftw_2d_c2r(double** ifft_result,complex** input_array,int Nx,int Ny)
{	
	//make it symmetric first
	complex_2d_sym(input_array,Nx,Ny);
		
	fftw_complex *in, *out;
	fftw_plan     plan_b;

    /* Allocate input & output array */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

    /* Create plans */
    plan_b = fftw_plan_dft_2d(Ny,Nx, in, out, FFTW_BACKWARD,  FFTW_ESTIMATE);
    
    /* Populate input data */
    for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){ 
			in[i*Nx+j][0] = input_array[i][j].Re; //real part of input
			in[i*Nx+j][1] = input_array[i][j].Im; //imajiner part
		}
    }

    /* Compute inverse DFT */
    fftw_execute(plan_b);  

	for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){
			ifft_result[i][j]=out[i*Nx+j][0]/(Nx*Ny);
		}
	}
    
    /* Free memory */
    fftw_destroy_plan(plan_b);
    fftw_free(in);
    fftw_free(out);	
}

void real_2d_sym(double** input_array,int Nx,int Ny)
{	
	for (int j=1;j<Nx/2;j++){
		input_array[0][Nx-j]    = input_array[0][j];
		input_array[Ny/2][Nx-j] = input_array[Ny/2][j];
	}
	for (int i=1;i<Ny/2;i++){
		input_array[Ny-i][0] = input_array[i][0];
	}
	
	int iy = 1;
    int it = 1;
    for (int i=1;i<Ny;i++){
		iy=1;
		for (int j=Nx/2+1;j<Nx;j++){
			input_array[i][j] = input_array[Ny-it][Nx/2-iy];
			iy++;
		}
		it++;
	}	
}

complex** ifftw_2d_r2c(double** input_array,int Nx, int Ny)
{	
	//make it symmetric first
	real_2d_sym(input_array,Nx,Ny);
	
	complex** ifft_result=declare_2darray_complex(Nx,Ny);	
	fftw_complex *in, *out;
	fftw_plan     plan_b;

    /* Allocate input & output array */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

    /* Create plans */
    plan_b = fftw_plan_dft_2d(Ny,Nx, in, out, FFTW_BACKWARD,  FFTW_ESTIMATE);
  
    /* Populate input data */
    for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){ 
			in[i*Nx+j][0] = input_array[i][j]; //real part of input
			in[i*Nx+j][1] = 0; //imajiner part
		}
    }

    /* Compute inverse DFT */
    fftw_execute(plan_b);  

	for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){
			ifft_result[i][j].Re=out[i*Nx+j][0]/(Nx*Ny);
			ifft_result[i][j].Im=out[i*Nx+j][1]/(Nx*Ny);
		}
	}
    
    /* Free memory */
    fftw_destroy_plan(plan_b);
    fftw_free(in);
    fftw_free(out);	
    
    return ifft_result;
}

double** ifftw_2d_r2r(double** input_array,int Nx, int Ny)
{	
	//make it symmetric first
	real_2d_sym(input_array,Nx,Ny);
	
	double** ifft_result=declare_2darray(Nx,Ny);	
	fftw_complex *in, *out;
	fftw_plan     plan_b;

    /* Allocate input & output array */
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

    /* Create plans */
    
    plan_b = fftw_plan_dft_2d(Ny,Nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
    
    /* Populate input data */
    for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){ 
			in[i*Nx+j][0] = input_array[i][j]; //real part of input
			in[i*Nx+j][1] = 0; //imajiner part
		}
    }

    /* Compute inverse DFT */
    fftw_execute(plan_b);  

	for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){
			ifft_result[i][j]=out[i*Nx+j][0]/(Nx*Ny);
		}
	}
    
    /* Free memory */
    fftw_destroy_plan(plan_b);
    fftw_free(in);
    fftw_free(out);	
    
    return ifft_result;
}

double* fft_1d_real_real(double* input_array,int Nx)
{
	double* result = new double[Nx];
	fftw_complex *in	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_complex *out 	 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx));
	fftw_plan plan_ft_h  = fftw_plan_dft_1d(Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	
	for (int j = 0; j < Nx; j++){ 
		in[j][0] = input_array[j]; //real part of input
		in[j][1] = 0; //imajiner part
	}	
	
    fftw_execute(plan_ft_h); 
    
    for (int j = 0; j < Nx; j++){
		result[j] = out[j][0];
	}
    
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan_ft_h);	
	
	return result;
}

void fftw_2d_r2c(complex** fft_result,double** input_array,int Nx, int Ny)
{
	fftw_complex *in, *out;
	fftw_plan     plan_f;

	  /* Allocate input & output array */
	  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	  /* Create plans */
	  plan_f = fftw_plan_dft_2d(Ny,Nx,in,out, FFTW_FORWARD,  FFTW_ESTIMATE);
	  
	  /* Populate input data */
	  for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){ 
		in[i*Nx+j][0] = input_array[i][j]; //real part of input
		in[i*Nx+j][1] = 0; //imajiner part
		}
	  }

	  /* Forward DFT */
	  fftw_execute(plan_f);  

  for (int i = 0; i < Ny; i++){
	for (int j = 0; j < Nx; j++){
     fft_result[i][j].Re=out[i*Nx+j][0];
     fft_result[i][j].Im=out[i*Nx+j][1];
	}
  }
    
  /* Free memory */
  fftw_destroy_plan(plan_f);
  fftw_free(in);
  fftw_free(out);	
}

complex** fftw_2d_c2c(complex** input_array,int Nx, int Ny)
{
	complex** fft_result=declare_2darray_complex(Nx,Ny);
	fftw_complex *in, *out;
	fftw_plan     plan_f;

	  /* Allocate input & output array */
	  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));
	  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nx*Ny));

	  /* Create plans */
	  plan_f = fftw_plan_dft_2d(Ny,Nx,in,out, FFTW_FORWARD,  FFTW_ESTIMATE);
	  
	  /* Populate input data */
	  for (int i = 0; i < Ny; i++){
		for (int j = 0; j < Nx; j++){ 
		in[i*Nx+j][0] = input_array[i][j].Re; //real part of input
		in[i*Nx+j][1] = input_array[i][j].Im; //imajiner part
		}
	  }

	  /* Forward DFT */
	  fftw_execute(plan_f);  

  for (int i = 0; i < Ny; i++){
	for (int j = 0; j < Nx; j++){
     fft_result[i][j].Re=out[i*Nx+j][0];
     fft_result[i][j].Im=out[i*Nx+j][1];
	}
  }
    
  /* Free memory */
  fftw_destroy_plan(plan_f);
  fftw_free(in);
  fftw_free(out);	
  return fft_result;
}

double** fun_abs_complex(complex** array,int N,int M)
{
	double** a = declare_2darray(N,M);
	for (int i=0;i<M;i++){
		for (int j=0;j<N;j++){
			a[i][j]=sqrt(pow(array[i][j].Re,2)+pow(array[i][j].Im,2));
		}
	}
	return a;
}

double* fun_array_substract_const(double* array,double constanta, int N)
{
	double* arrayN=new double[N];
	for(int i=0;i<N;i++){
		arrayN[i]=array[i]-constanta;
	}
	return arrayN;
}

int fun_closest(double* A, int N,double Aref)
{
	int idx=0;
	double Del=fabs(Aref-A[0]);
	double DelN;
	for(int i=0;i<N;i++){
		DelN=fabs(Aref-A[i]);
		if (DelN<Del){
			Del=DelN;
			idx=i;
		}
	}
	return idx;
}

complex fun_complex_multiplication(complex A, complex B){
	complex result;
	result.Re=A.Re*B.Re-A.Im*B.Im;
	result.Im=A.Re*B.Im+A.Im*B.Re;
	return result;
}

complex fun_complex_division(complex A, complex B){
	complex result;
	result.Re=(A.Re*B.Re+A.Im*B.Im)/(pow(B.Re,2)+pow(B.Im,2));
	result.Im=(-A.Re*B.Im+A.Im*B.Re)/(pow(B.Re,2)+pow(B.Im,2));

	return result;
}

double fun_min(double* array,int N)
{
	double min;
	min=array[0];	
	for(int j=0;j<N;j++){
		min=(min<array[j]?min:array[j]);
	}
	return min;
}

double fun_max(double* array,int N)
{
	double max;
	max=array[0];	
	for(int j=0;j<N;j++){
		max=(max>array[j]?max:array[j]);
	}
	return max;
}

double fun_mean2darray(double** arrayIn,int Nx, int Ny)
{
double meanI=0;int iter=0;
for (int i=0;i<Ny;i++){
	for (int j=0;j<Nx;j++){
	meanI+=arrayIn[i][j];
    iter++;
	}
}
return (meanI/(iter+1));
}

double fun_var2darray(double** arrayIn,int Nx, int Ny)
{
	double varI=0;int iter=0;
	double meanI=fun_mean2darray(arrayIn,Nx,Ny);
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Nx;j++){
			varI+=pow(arrayIn[i][j]-meanI,2);
			iter++;
		}
	}
	return (varI/(iter+1));	
}

void gradient(double* array, double d, int n, double* grad)
{
	grad[0] 	= (array[1]-array[0])/d;
	grad[n-1] 	= (array[n-1]-array[n-2])/d;
	for (int i=1;i<n-1;i++) {grad[i]=(array[i+1]-array[i-1])/(2*d);}
}

void gradient2D(double** array, double dx, double dy, int nx, int ny, double** gradx,double** grady)
{
	for (int i=0;i<ny;i++){
		gradx[i][0] 	= (array[i][1]-array[i][0])/dx;
		gradx[i][nx-1] 	= (array[i][nx-1]-array[i][nx-2])/dx;
		for (int j=1;j<nx-1;j++) {
			gradx[i][j]=(array[i][j+1]-array[i][j-1])/(2*dx);
		}
	}

	for (int j=0;j<nx;j++){
		grady[0][j] 	= (array[1][j]-array[0][j])/dy;
		grady[ny-1][j] 	= (array[ny-1][j]-array[ny-2][j])/dy;
		for (int i=1;i<ny-1;i++) {
			grady[i][j]=(array[i+1][j]-array[i-1][j])/(2*dy);
		}
	}
}

void inner_product_array(double* shat_square,double* array1,double* array2,int n)
{
	for(int i=0;i<n;i++)shat_square[i]=array1[i]*array2[i];
}

void sorting(double* arr,int size)
{
	double temp;
	for (int i=0; i<size; i++) {
        for (int j=i+1; j<size; j++) {
            //If there is a smaller element found on right of the array then swap it.
            if(arr[j] < arr[i]) {
                temp 	= arr[i];
                arr[i] 	= arr[j];
                arr[j] 	= temp;
            }
        }
    }    
}

double mean_array(double* wave,int n)
{
	double signal=0;
	for (int i=0;i<n;i++){
		signal+=wave[i];
	}
	return (signal/n);
}

double var_array(double* wave,int n,double avrg)
{
	double var=0;
	for (int i=0;i<n;i++){
		var+=pow(wave[i]-avrg,2);
	}
	return (var/n);
}

void set_matrix_val(double** array, int N, int M, double val)
{
	#pragma parallel for
	for (int j=0;j<N;j++){
		for (int i=0;i<M;i++){
			array[i][j] = val;
		}
	}
}

void set_matrix_complex_val(complex** array, int N, int M, double val)
{
	for (int j=0;j<N;j++){
		for (int i=0;i<M;i++){
			array[i][j].Re = val;
			array[i][j].Im = val;
		}
	}
}

double omega_i(double k_i, double d_i)
{
	double om;
	if (k_i==0){
		om = 0;
	}
	else {
		om = (k_i<0?-1:1)*sqrt(grav*k_i*tanh(d_i*k_i));	
	}
		
	return om;
}

double group_velocity_i(double k_i,double d_i)
{
	double Ug_i;
	
	if(k_i==0) {
		Ug_i=sqrt(grav*d_i);
	}
	else {
		Ug_i=(k_i<0?-1:1)*sqrt(grav)/2/sqrt(k_i*tanh(d_i*k_i))*(tanh(d_i*k_i)+k_i*(1-pow(tanh(d_i*k_i),2))*d_i );
	}
		
	return Ug_i;
}

double my_f (double x, void * params)
{
	fsolve_p *p=(fsolve_p *) params;
	double om=p->w;
	double depth1=p->d;     
	return om-omega_i(x,depth1);
}
     
double my_df (double x, void * params)
{
	fsolve_p *p=(fsolve_p *) params;
	double depth1=p->d;             
	return -group_velocity_i(x,depth1);
}
     
void my_fdf (double x, void * params, double * f, double * df)
{
	fsolve_p *p=(fsolve_p *) params;
	double om=p->w;
	double depth1=p->d;
		
	*f = om-omega_i(x,depth1);
	*df = -group_velocity_i(x,depth1);
}    
 
double invers_omega(double w1,double depth1)
{
	int status;double value=0.;
    int iter = 0, max_iter = 1000;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0, x = 0.; //r_expected = sqrt (0.);
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
        status = gsl_root_fdfsolver_iterate (s);
        x0 = x;
        x = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta (x, x0, 0, 1e-9);
     
        if (status == GSL_SUCCESS){
             value=x;
        }
	}
    while (status == GSL_CONTINUE && iter < max_iter);
		gsl_root_fdfsolver_free (s);
    return value;	
}

void fun_interp1(int Nxin, double* xin, double* yin, int Nxout, double* xout, double* yout)
{
	gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
    gsl_spline *spline 	   = gsl_spline_alloc (gsl_interp_cspline,Nxin);
    gsl_spline_init (spline,xin,yin,Nxin);
	
    for (int j=0;j<Nxout;j++){
		if ((xout[j]>=xin[0])&&(xout[j]<=xin[Nxin-1])){
			yout[j] = (gsl_spline_eval(spline,xout[j],acc));
		}
		else {
			yout[j] = NAN;
		}
	}

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}

void fun_interp1_float(int Nxin, double* xin, double* yin, int Nxout, double* xout, float* yout)
{
	gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
    gsl_spline *spline 	   = gsl_spline_alloc (gsl_interp_cspline,Nxin);
    gsl_spline_init (spline,xin,yin,Nxin);
	
    for (int j=0;j<Nxout;j++){
		if ((xout[j]>=xin[0])&&(xout[j]<=xin[Nxin-1])){
			yout[j] = float(gsl_spline_eval(spline,xout[j],acc));
		}
		else{
			yout[j] = NAN;
		}
	}

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}

void fun_interp2_float(int Nxdata,int Nydata,double* xdata,double* ydata,double** data,int Nx,int Ny,double* x,double* y,float** results)
{
	double* longdata=new double[Nxdata*Nydata];
	for (int i=0;i<Nydata;i++){
	  for(int j=0;j<Nxdata;j++){
		  longdata[i*Nxdata+j]=data[i][j];
	  }
    }
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, Nxdata, Nydata);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();		
	gsl_spline2d_init(spline, xdata, ydata, longdata, Nxdata, Nydata);
	
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  if ((x[j]>=xdata[0])&&(x[j]<=xdata[Nxdata-1])&&(y[i]>=ydata[0])&&(y[i]<=ydata[Nydata-1])){
			  results[i][j] = float (gsl_spline2d_eval(spline,x[j],y[i],xacc,yacc));
		  }
		  else{
			  results[i][j] = NAN;
		  }
	  }
    }
    
    delete[] longdata;
    gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc); 
}

void fun_interp2(int Nxdata,int Nydata,double* xdata,double* ydata,double** data,int Nx,int Ny,double* x,double* y,double** results)
{
	double* longdata=new double[Nxdata*Nydata];
	for (int i=0;i<Nydata;i++){
	  for(int j=0;j<Nxdata;j++){
		  longdata[i*Nxdata+j]=data[i][j];
	  }
    }
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, Nxdata, Nydata);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();		
	gsl_spline2d_init(spline, xdata, ydata, longdata, Nxdata, Nydata);
	
	set_matrix_val(results,Nx,Ny,-9999);
	for (int i=0;i<Ny;i++){
	  for(int j=0;j<Nx;j++){
		  if ((x[j]>=xdata[0])&&(x[j]<=xdata[Nxdata-1])&&(y[i]>=ydata[0])&&(y[i]<=ydata[Nydata-1])){
			  results[i][j] = gsl_spline2d_eval(spline,x[j],y[i],xacc,yacc);
		  }
	  }
    }
    delete[] longdata;
    gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc); 
}
