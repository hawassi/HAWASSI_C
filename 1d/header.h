#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <complex.h>
#include <vector>
#include <fftw3.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <omp.h>
#include <netdb.h> 
#include <sys/socket.h> 
#include <netinet/in.h> 
#include <arpa/inet.h> 

double const g      = 9.80665;
double const Pi     = 4.0 * atan(1.0);
double const delb   = 1.2;

int cutfracwn;

typedef struct{
    double Re;
    double Im;
}complex;

typedef struct{
    double eta;
    double t;
}marin;

typedef struct{ 
    double w;
    double d;
}fsolve_p;

typedef struct{
    int     Bindex;
    int     BindexNow;
    double  Ucrest;
    int     DirProp;
    double  Bposition;
    double  Btime; 
    int     Bcount;
}CB;

typedef struct{
    double KBC;
    double TC;
    double Tstar;
    double delb;
    int	   xbstart;
    int    xbend;
    int    xbendL;
    int    xbstartR; 
}parbreak;

typedef enum {
    REAL,
    IMAG
}num_type;

typedef enum {
    MAXIMUM,
    MINIMUM
}search_type;

parbreak pb;
CB* 	 CrestBreak;
double*  B;
double*  mid = new double[2];
int 	 Npks, Npks1, Npks2;
double   rel;
clock_t  crel;
double 	 Cpeak, CpLeft, CpRight;

int const spaceinterpnum = 1;
char 	  wavename[20], initial[20], bath[20], dynmodel[20], dispersion[20], breaking[20], friction[20], influxing[20], propagation[20], kinematic[20], wind[20], wall[20];
int       walltype;
double    lambda_p, dampL, dampR, MaxEtaInit;
double    T0, Tend, dt, Tp, Hs, gamJS, Amp, xc, sigma;
double    lenT, omgp;
int       Nt, ngauge;
double 	  adjcoef;
double 	  nu_peak;  
double    *ChiAdj, *gauge;
double    parKBC, parTC, parTchar;
int 	  lenIF, ILL, IFF;
double    Hs_in;
double    Tp_peak;
double 	  avrg_wave;
char 	  cut_temp[20], mid_temp[20], min_z[20], max_z[20];
int 	  runup;

int  	  flag = 0, iterNprev, add, indxMaxNow, NCBprev, ITERbn = 0, ITERbdt = 0, Node_end;  
CB   	  *CrestBreakPrev, *CrestBreakNew;
int  	  *ID_Node_end, *Break_nodes; 
double    *dataBreak_nodes;
char 	  argv2[256], *argv1;
struct    timeval start, end;
int 	  ncore, ncore2;
double    *temp_wave, *temp_u;
int 	  primeNx;
double 	  nTramp;
char 	  rampbool[20];

double    x_cpl[2];
int 	  idx_cpl[2];
char 	  coupling[20];
char 	  dir_force[256];
double*   Chi_couple;
double*   etadata;
double*   udata;
double    alpha_couple;
double    alpha_corr;

int    	  ITERdtW,idx_wall;
double	  Swall_tS_eta, Swall_tS_eta_prev;
double 	  tprevW, Xwall, depthXwall;
complex   *gamWall_hat, *gamWall_skew_hat, *Swall_hat;
double    refl_coef;
double    Swt_skew;

double    depthflat, Xstart, Xinflux, Xend, depthinf;
double    *Rho_k, *Rho_Om, *Refl_coef_Om, *Refl_coef_k;
double    *refl_coef_v, *refl_coef_k;
int       idx_inf, SignProp;

int       tcoarse = 1;
int       p, Nx;       
double    dx, kp, length;
double    tprev;

marin* 	  orgwave;
marin* 	  wave;
double*   insig, *insig_skew;
double 	  *depth;
double 	  *x, *k, *dampchar, *spatgamX, *Ug, *aal, *GVom, *kom, *dampcoef;
double 	  *Csqr;
complex	  *gamX, *insig_hat;
complex	  *gamX_skew;
complex	  *Rf_hat, *Rw_hat, *Rb_hat;
double 	  *ww, *spinflux, *Cf, *Cw, *ttrace;
int 	  n;

double*   L, *L0;
complex   *M0, *M1;
double	  *Csq_min, *Csq_mid, *Csq_plus, *gam_min, *gam_mid, *gam_plus;
complex*  Csqr_v;
double*   intpnu;
double*   tim;
double*   dy;
FILE* 	  fpo;
double 	  left, right;

double 	  *Cderm, *Cderp, *Cderc, Hmin_interp, Hplus_interp;
double 	  *Cm, *Cp, *Cc;
int 	  signProp, Ninterp = 101;

FILE* 	  datconst, *dateta_full, *datu_full, *datphi_full, *dateta_gauge, *datphi_x_gauge, *datphi_z_gauge, *dateta_t_gauge;
FILE* 	  fbreak;
double 	  *C2p1, *C2m1, *C2c1, *C2p2, *C2m2, *C2c2; 		

double 	  zsurf1, zsurf2;
double*   xoutIP;
double*   zoutIP;
double**  toutIP;
double*   trest;
int	  NxoutIP, NzoutIP, indxIP1, indxIP2, NtoutIP, Ntrest;
double 	  maxD, minD, maxEta, minEta, dtout;
double*   bathy;
double    **Velx, **Velz;
int 	  npartition, nsurf, ndeep, ip = 0, itloc, Nheader;
double 	  xinterv1, xinterv2, tfirst, tlast;
FILE* 	  dateta, *datvel, *datacc, *datP;

fftw_plan plan_ft_nx;
fftw_plan plan_ift_nx;
fftw_plan plan_ift_c2r;
fftw_plan plan_ft_r2c;

typedef std::vector<marin> Vmarin;
typedef std::vector<double> outCB;

gsl_interp_accel *acc;
gsl_spline *spline;

gsl_interp_accel *acc_gam;
gsl_spline *spline_gam;

gsl_interp_accel *accskew;
gsl_spline *splineskew;

gsl_interp_accel *acc_gamXom;
gsl_spline *spline_gamXom;

gsl_interp_accel *acc_sgp = gsl_interp_accel_alloc();
gsl_spline *spline_sgp    = gsl_spline_alloc(gsl_interp_cspline, Ninterp);
		
gsl_interp_accel *acc_sgc = gsl_interp_accel_alloc();
gsl_spline *spline_sgc    = gsl_spline_alloc(gsl_interp_cspline, Ninterp);  
	
gsl_interp_accel *acc_sgm = gsl_interp_accel_alloc();
gsl_spline *spline_sgm    = gsl_spline_alloc(gsl_interp_cspline, Ninterp);

gsl_interp_accel *acc_p1  = gsl_interp_accel_alloc();
gsl_spline *sg_p1 = gsl_spline_alloc(gsl_interp_cspline, 199);

gsl_interp_accel *acc_c1  = gsl_interp_accel_alloc();
gsl_spline *sg_c1 = gsl_spline_alloc(gsl_interp_cspline, 199);

gsl_interp_accel *acc_m1 = gsl_interp_accel_alloc();
gsl_spline *sg_m1 = gsl_spline_alloc(gsl_interp_cspline, 199);

gsl_interp_accel *acc_p2 = gsl_interp_accel_alloc();
gsl_spline *sg_p2 = gsl_spline_alloc(gsl_interp_cspline, 199);

gsl_interp_accel *acc_c2 = gsl_interp_accel_alloc();
gsl_spline *sg_c2 = gsl_spline_alloc(gsl_interp_cspline, 199);;

gsl_interp_accel *acc_m2 = gsl_interp_accel_alloc();
gsl_spline *sg_m2 = gsl_spline_alloc(gsl_interp_cspline, 199);

marin* influxlib(int*, char*);
int  countlines(FILE*);
void mod_influx(marin*, int, double*, double*, double*, marin*, marin*);
void set_influx(marin*, int, double*, double*, double*);
void x_space(double*, double*, double*, complex*, double*, double*, complex*);
void plot_initial_config(marin*, double*, int, double*, double*, double*, double*, double, char*, double*);
void alias(double*, double*, double);
void operator_L(double*, double*, double);
void operator_L0(double*, double*);
void operator_M0(complex*, double*);
void operator_M1(complex*, double*);

void ode_solver(gsl_odeiv_system);
void ode_solver_Ham(gsl_odeiv_system);
int rhs_HsF1 (double, const double*, double*, void *);
int rhs_HsF2 (double, const double*, double*, void *);
int rhs_HsF3 (double, const double*, double*, void *);
int rhs_HsF4 (double, const double*, double*, void *);
int rhs_shore(double, const double*, double*, void *);

void fft_real_complex(double*, complex*, int);
void fft_real_complex_nx(double*, complex*);
void ifft_complex_real(complex*, double*, int);
void ifft_complex_real_nx(complex*, double*);
void ifft_fft_Hsf1(const double*, double*, fftw_complex*, fftw_complex*, double*, int);				
void ifft_fft_Hsf2(const double*, double*, fftw_complex*, fftw_complex*, fftw_complex*, fftw_complex*, fftw_complex*, double*, int);				
void ifft_fft_Hsf3(const double*, double*, fftw_complex*, fftw_complex*, fftw_complex*, int);
void ifft_fft_Hsf4(const double*, double*, fftw_complex*, fftw_complex*, fftw_complex*, int); //hs4 flat
void ifft_fft_Hsf2_v(const double*, double*, fftw_complex*, fftw_complex*, fftw_complex*, fftw_complex*, fftw_complex*, double*, int);
void ifft_fft_Hsf3_v(const double*, double*, fftw_complex*, fftw_complex*, fftw_complex*, int);
         
void set_gnuplot_condition2(FILE*, char*);
void plot_particle2(double const*, double const*, double, FILE*, int);
void print_data(double*, double, FILE*, int);
void print_influx(FILE*, marin*, int);
void print_const_data(FILE*, int);
void set_memory();
void free_memory();

void filtering_Sig(double, double, double);
void harmonic_influx(double, double, double, double, double, int*);
void jonswap_influx(double, double, double, double, int, double*, double, double*, double, double, int*);
double* fun_jonswap(double*, double*, double, double*, double, double, double, int);
void rhs_friction(complex*, double*, double*, double*, int, double*, double*);
void rhs_wind(complex*, double*, double, double*, complex*);

double invers_omega(double, double);
double omega_i(double, double);
double group_velocity_i(double, double);
double prod_normal(complex, complex, num_type);
double prod_fftw(complex, double, double, num_type);
double prod_1ik(double, double, double, num_type);

void omega(double*, double*, double, int);
void gradient(double*, double, int, double*);
void dampzone(double*, double*, double, double, double, double);
void inner_product_real_complex(complex*, double*, complex*, int);
void kinematic_modul(double, double*, complex*, double*, complex*);
void data_saving_hawassi(FILE*, FILE*, FILE*, FILE*);
void modified_influx(marin*, int, double*, double*, double*, double*, double*, double*, double*);
int closest(double*, double, int);
void parambreak(double*, double, double, double, double); 
void breaking_process(double*, complex*, double*, double, double*, parbreak);
double mean(marin*, int);
double variance(marin*, int, double);
void set_init_CB(CB*, int);
void rhs_breaking(complex*, double*, complex*, double*, double, double*, double*, parbreak, complex*);
void logfile(char*, char*, int);
void set_init_double(double*, int);
void set_init_break(int**, int, int);
void set_init_int(int*, int);
int allvalueless(double, double*, int);
int idx_max(double*, int n);
void data_saving_break(FILE*, double*, int);
void print_xt(FILE*, int, int);
void define_bath(char*);
void define_friction(char*);
void define_wind(char*);
void define_gauges(char*);
void set_init_complex(complex*, int);
void save_all_output(char*, char*);
void ram_signal();
void smooth(double*, double*,int, int);
complex* define_forcing(double, int);
void adjust_deluH2(complex*, const double*, double*, fftw_complex*, fftw_complex*);
void adjust_delphiH2(complex*, const double*, double*, fftw_complex*, fftw_complex*);
double phase_velocity_i(double, double);
double* fun_csqr(double*, double, int);
void set_operator(int);
void adjust_delH(complex*, complex*, complex*, complex*);
void HSS(complex*, complex*);
void HSS_runup(complex*, complex*, double*, double*, double*, double*, double*, double*);
void define_operator(int, double);
void define_delH1(const double*, complex*, complex*);
void Up2IP(double);
void Up3IP(double);
void circshift(double*, int, int);
void ifft_real_real(double*, double*, int);
void read_kinematics();
//double search_vmax(int, double*);
//double search_vmin(int, double*);
double search_v(int, double*, search_type);
//double search_mmax(int, int, double**);
//double search_mmin(int, int, double**);
double search_m(int, int, double**, search_type);
double mean_array(double*, int);
void search_max(double const*, int, double*, int*);
void absolute_complex(double*, complex*, int);
double var_array(double*, int, double);
int primeFactors(int);
void Up3IP_runup(double);
int indexOfFirstZero(double*, int);
int indexOfLastZero(double* , int);
int IsExpired(char*);
void surf_velocity_z(double, FILE*, FILE*, complex*, double*, double*, complex*);
double trapz(double*, double*, int);
void fun_wallinfluxsource(double, double*);
void define_wall(complex*, complex*);
void scaling_influx(complex*, double, double, double);
double trapz_non(double*, double*, int);
double* cf(double*, double, double, int);

void getChi(int*, int);
void coupling_term(double*, double*, complex*, complex*, double, double*, double*);

void read_constant(int*, int*, double*, double*);
void read_eta(int, int, double*);
void read_u(int, int, double*);
void read_input_file(const char*, const char*);
void read_initial_file(const char*, char*);
void read_breaking_file(const char*, char*);
void read_influx_file(const char*, char*);
void read_bath_file(const char*, char*);

void check_runup_dynmodel(int, char*, size_t);
void wall_setup(double*, double*, double*);
void update_damp(int, int, double*, double*);
