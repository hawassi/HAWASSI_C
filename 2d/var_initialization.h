FILE* fpin;
FILE* fpout;

typedef struct{
	double start;
	double end;	
	double delta;
}interv;

typedef struct{
	double xstart;
	double ystart;
	double xend;
	double yend;
}lineflux;

typedef struct{
	double x;
	double y;
}coord;

typedef struct{
	double x;
	double y;
}Vec2d;

typedef struct{
	int x;
	int y;
}Vec2dint;

typedef struct{
	int hr;
	int min;
	int sec;
	int msec;
}hrms;

typedef struct{
	double Hs;
	double Tp; //peak period
	double gamJS;
	char   s[100];
	double wp; //peak freq
	double kp;
	double lambda;
	Vec2d  current;
	double mdir;
}seavar;

typedef struct{
	double Lleft;	
	double Lright;
	double Lbottom;
	double Ltop;
	double** charac;
}fbdyvar;

typedef struct{
	char Id[100];
	char Filename[100];
	double** data;
	double** pos;
	double** plus;
	double** min;
	double** Char;
	double depth;
	char type[1];
}bathyvar;

typedef struct{
	double* x;
	double* y;
	double* t;
	int Nx;
	int Ny;
	int Nt;
	double dx;
	double dy;
	double dt;
	int tcoarse;
	double Nadj;
	double* kx;
	double* ky;
	double** KK;
	double** aal;
	bathyvar bathy;
	fbdyvar  fbdy;
	double 	 U;
	double** wallchar;
}domvar;

typedef struct{
	char model[20];
	int  cutfrac;	
}evolvar;

typedef struct{
	double** Gam_min;
	double** Gam_mid;
	double** Gam_plus;	
	double  D_min;
	double  D_mid;
	double  D_plus;
	double** L2d_min;
	double** L2d_mid;
	double** L2d_plus;
	double** Om2d_min;
	double** Om2d_mid;
	double** Om2d_plus;
	double** Om2dSq_min;
	double** Om2dSq_mid;
	double** Om2dSq_plus;
	double** C2d_min;
	double** C2d_mid;
	double** C2d_plus;
	double** Csqr2d_min;
	double** Csqr2d_mid;
	double** Csqr2d_plus;
}Interp2Dvar;

typedef struct{
	double** Om2d;
	double** L2d;
	double** C2d;
	double** Csqr2d;
	double** Cpeak2d;
	double   Cpeak;
	double 	 HminDisp;
	Interp2Dvar InterpD;
}Oprtvar;

typedef struct{
	double** profile;
	double** phi;
	hrms     time_hms;
	double   time;
}outvar;

typedef struct{
	double** profile;
	double** phi;
	double** u;
	double** v;
	hrms     time_hms;	
}wavestruc;

typedef struct{
	double Re;
	double Im;
}complex;

typedef struct{
	double** eta;
	double*  x;
	double*  y;
	double*	 time;
	double*  ww;
	double*  kw;
	double 	 wmax;
	double 	 kmax;
	int 	 Nxy;
	int 	 nlines;
	int 	 idx0;
	int 	 idy0;
}waveinf;

typedef struct{
	double KBC;
	double TC;
	double Tchar;
	double delb;
	double MaxEtaInit;
	double Tstar;
	double Rsearch;
	double** Char;
}parbreak;

typedef struct{
	gsl_interp_accel* acc;
	gsl_spline* spline;
	int init;
}interp_sig;

typedef struct{
	int xindex;
	int yindex;
	int xNow;
	int yNow;
	double Ucrest;
	double theta;
	int    kwadran; 
	double time;
	int    NB;
	Vec2dint* nodes;
	Vec2dint* nodes_prev;
	int    nnode;
	int    nnode_prev;
}CBvar;

double**  newGam;
double**  SkewGam;
complex** Insig_hat;
complex** Source_hat;
complex **Sf_xhat,**Sf_yhat;
complex **Sb_xhat,**Sb_yhat;
complex **Sw_hat;
interp_sig* Spline_Signal;
interp_sig* Spline_SkewSignal;

parbreak pb; 
evolvar  evolinfo; lineflux lineinflux;
seavar   seainfo; domvar dominfo; 
interv   rad_range,rad_theta,dom_x,dom_y,dom_t;
Oprtvar  Oprt;	

//for breaking global
CBvar *CB,*CBNew,*CBPrev;
int itCB,iterdtcheck,ITERbspdt=0,ITERbdt=0;
int ITERbn=0,iterbn=0,iterNprev=0,iterNprev_prev=0;
int nbreak_nodes,nbreak_prev;
double** B;
Vec2dint*  Break_nodes;
Vec2dint*  Break_nodes_prev;
double  tbreak;

//for wall global
double tprevW;
int    iterW,n_idxwall;
Vec2d* idxWall;
int*   dirWall;
int    nwall;
double*  Xwall[2];
double*  Ywall[2];
double** dxGam1;
double** dxGam2;
double** dyGam1;
double** dyGam2;
double** Swt_skew;
double** grad_wallchar_x;
double** grad_wallchar_y;
double*  Swall_ts;
double*  Swall_prev;
double*  Swt_skewline;

double const grav =9.80665;

//move to header later
char *arg1,arg2[256];
int  ncpu,ngauge,npartition;
char 		wavename[20],initial[20],bath[20],breaking[20],wallmethod[20],
			friction[20],influxing[20],propagation[20],kinematic[20],kinematic_type[20],
			wind[20],wall[20],orientation[20],coupling[20],dir_force[256];
double parKBC,parTC,parTchar,refl_coef;
int    nsurf,ndeep,walltype;
char   cut_temp[20],mid_temp[20],min_z[20],max_z[20];
wavestruc waveinit;
waveinf   wave_influx;
double 	  meandepthinf;
struct   timeval start, end;
int*	 index_gauge[2];
double*  xgauge;
double*  ygauge;
double** ChiAdj;
double** ChiAdjPhi;
double** cfSA; 
double** cfSAwall;
double** ramp2d;
double** ramp2dori;
double** fric_coef;
double** fric_depth;
double** wall_gam;
double** wall_skewgam;
double   ramp;
double   taper;
int 	 ninterp;
int 	 flagbr = 0;
double   tprev;
FILE* fbreak;
FILE* ftes;

int Nx, Ny;
double dx, dy;
int runup_id;
double H_minShore,Dmin_input;

//runup
double 		**Cp1,**Cm1,**Cc1,**Cp2,**Cm2,**Cc2; 		
double 		**C2p1,**C2m1,**C2c1,**C2p2,**C2m2,**C2c2; 		
double 		**Om2p1,**Om2m1,**Om2c1,**Om2p2,**Om2m2,**Om2c2; 
double 		**gam_p1,**gam_m1,**gam_c1,**gam_p2,**gam_m2,**gam_c2; 				

gsl_interp_accel *acc_p1;
gsl_spline *sg_p1;

gsl_interp_accel *acc_c1;
gsl_spline *sg_c1;

gsl_interp_accel *acc_m1;
gsl_spline *sg_m1;

gsl_interp_accel *acc_p2;
gsl_spline *sg_p2;

gsl_interp_accel *acc_c2;
gsl_spline *sg_c2;

gsl_interp_accel *acc_m2;
gsl_spline *sg_m2;

typedef struct{ 
	double w;
	double d;
}fsolve_p;

void 			read_input(FILE*);
int  			calNpoint(double,double,double);
void	 		fun_intervaltovector(double*,double,double,int);
void	 		fun_intervaltovector_ds(double*,double,double,int);
void    		fun_calcseacharacteristic(seavar*,waveinf);
double**    	declare_2darray(int,int);
complex**   	declare_2darray_complex(int,int);
void    		free_2darray(double**,int,int);
void	   		free_2darray_complex(complex**,int,int);
void        	fun_plot_line(double* ,double* , int);
double*     	fun_heaviside(double*,int,double);
double			fun_min(double*,int);
double			fun_max(double*,int);
double			fun_min2darray(double**,int,int);
double			fun_max2darray(double**,int,int);
double      	fun_mean2darraywithcond(double** ,int , int ,double**, int ,double);
double 			fun_var2darray(double** , int , int);
double 			fun_var2darraywithcond(double** , int , int ,double**, int ,double);
double 			fun_mean2darray(double**,int, int );
double      	fun_exact_disp_val(double, double, double, double, double,double);
double*      	fun_exact_disp1d(double* , double,int, double);
double      	fun_exact_disp1d_val(double , double,double);
double      	fun_exact_Cp1d_val(double , double,double);
void 			fun_exact_Ug(double* ,double* ,double ,int ,double );

complex**		fft_2d(double** ,int , int ); 
double**    	ifft_2d_Re(complex**,int , int );

void			fftw_2d_r2c(complex**,double** ,int , int ); 
complex**		fftw_2d_c2c(complex** ,int , int ); 
void	    	ifftw_2d_c2r(double**,complex**,int , int );
double**    	ifftw_2d_r2r(double**,int , int );
complex**    	ifftw_2d_r2c(double**,int , int );

void	 		freqspace(double*,double ,int );
int         	fun_sign(double);
complex     	fun_exp_ix(double);
complex     	fun_complex_multiplication(complex,complex);
complex     	fun_complex_division(complex, complex);
void 			fun_gnuplot_set(FILE*  );
void 			fun_plot_2d_array(FILE* ,double* ,double* ,double** , int ,int);
void        	fun_param_init();
Oprtvar			fun_operator_setup(domvar);
void	        fun_aal2D(double**,int,double*,double*,double**,int,int);
void	        fun_fbdy2D(double**,double*,double*,int,int,fbdyvar);
void			fun_cfSA2D(double**,double*,double*,int,int,fbdyvar);
int             fun_closest(double*,int,double);
void  			ode_solver(gsl_odeiv2_system ,double* ,double* ,int ,int,int ,int ,double ,double ,int,double,outvar*);
int             rhs_AB1_2DF(double , const double* , double* ,void *);
void            fun_termshat_AB1_2DF(fftw_complex* ,fftw_complex* ,domvar ,double ,const double* );
int             rhs_AB2_2DF(double , const double* , double* ,void *);
void            fun_termshat_AB3_2DF(fftw_complex* ,fftw_complex* ,fftw_complex* ,fftw_complex* ,fftw_complex* ,
						   fftw_complex* ,fftw_complex* ,fftw_complex* ,fftw_complex* , 
						   fftw_complex*, fftw_complex* ,fftw_complex* ,domvar ,Oprtvar ,const double*);
void	        fun_bathy2d(double**,int,int,bathyvar);		
void 			funOprt_DispInterpolation_2p(Oprtvar,double**,double**,double*,double*,int,int,double,Vec2d,double);
void 			funOprt_DispInterpolation_3p(Oprtvar,double**,double**,double*,double*,int,int,double,Vec2d,double);
double 			fun_invOm(int,double*,double*,double);		
void 			fun_termshat_AB1_2DB(fftw_complex* ,fftw_complex* ,fftw_complex* ,
                         fftw_complex* ,fftw_complex* ,fftw_complex*,
						 domvar ,Oprtvar,const double*);	
void 			fun_termshat_AB2_2DF_u(double ,fftw_complex* ,fftw_complex* ,fftw_complex* ,domvar ,double ,const double* );   
