void parambreak(int Nx, int Ny)
{
	if (strcmp(wavename,"zero")!=0){
		pb.MaxEtaInit= seainfo.Hs/2.0*0.9;              
		pb.Tstar 	 = pb.Tchar*seainfo.Tp;  
		pb.Rsearch   = 2*M_PI/seainfo.kp/4.0;
	}
	pb.Char = declare_2darray(Nx,Ny); 
	
	#pragma omp parallel for
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			pb.Char[i][j] = cfSA[i][j];
			if (cfSA[i][j]<1){
				pb.Char[i][j] = 0;
			}
			if (ChiAdj[i][j]<1){
				pb.Char[i][j] = 0;
			}
		}
	}

	//initiation
	ITERbdt = 0;
	CB 		= new CBvar[Nx*Ny];
	CBNew   = new CBvar[Nx*Ny];
	CBPrev  = new CBvar[Nx*Ny];
	for (int jj=0;jj<(Nx*Ny);jj++){
		CBPrev[jj].nodes = new Vec2dint[100];
		CBPrev[jj].nodes_prev = new Vec2dint[100];
	}
	B 		= declare_2darray(Nx,Ny);
	set_matrix_val(B,Nx,Ny,0.0);
	Break_nodes = new Vec2dint[Nx*Ny];
	Break_nodes_prev = new Vec2dint[Nx*Ny];
	tbreak = dominfo.t[0];
	
	if (strcmp(breaking,"yes")==0){
		char str_br[256];
		sprintf(str_br, "%sbreaking.dat",arg2);
		fbreak = fopen(str_br,"w");
	}
}

int counts(double t,double** array, int row, int col )
{
    // invariant row in [1...ROWS-2], col in [1...COLS-2]
    int cnt_high = 0 ;
    for( int i = -1 ; i < 2 ; ++i ) for( int j = -1 ; j < 2 ; ++j ){
        cnt_high += ((array[col+i][row+j]) - (array[col][row]) <= 0.001) ? 1 : 0 ;
    }
    return cnt_high;
}

int fun_peaks(double t,double** peak,int* IndPeak[2], double** EtaChecked,int ROWS,int COLS)
{
	set_matrix_val(peak,ROWS,COLS,0.0);
	int cnts;
	int idd=0;
    for( int row = 1 ; row < (ROWS-1) ; ++row ){
        for( int col = 1 ; col < (COLS-1) ; ++col ){		
			cnts = counts( t,EtaChecked, row, col ) ;
			if ((cnts==9)&&(EtaChecked[col][row]>pb.MaxEtaInit)){
				peak[col][row] = 1;
				IndPeak[0][idd] = row;
				IndPeak[1][idd] = col;
				idd++;
			}
        }
	}

	//top part
	int row = 0; int cnt_high;
    for( int col = 1 ; col < (COLS-1) ; ++col ){
        cnt_high = 0 ;
		for( int i = -1 ; i < 2 ; ++i ) for( int j = 0 ; j < 2 ; ++j ){
			cnt_high += (EtaChecked[col+i][row+j] - (EtaChecked[col][row]) <= 0.001) ? 1 : 0 ;
		}
        if ((cnt_high==6)&&(EtaChecked[col][row]>pb.MaxEtaInit)){
			peak[col][row] = 1;
			IndPeak[0][idd] = row;
			IndPeak[1][idd] = col;
			idd++;
        }
	}

	//bottom part
	row = ROWS-1; 
    for( int col = 1 ; col < (COLS-1) ; ++col ){
        cnt_high = 0 ;
		for( int i = -1 ; i < 2 ; ++i ) for( int j = -1 ; j < 1 ; ++j ){
			cnt_high += (EtaChecked[col+i][row+j] - (EtaChecked[col][row]) <= 0.001) ? 1 : 0 ;
		}
        if ((cnt_high==6)&&(EtaChecked[col][row]>pb.MaxEtaInit)){
			peak[col][row] = 1;
			IndPeak[0][idd] = row;
			IndPeak[1][idd] = col;
			idd++;
        }
	}

	//left part
	int col = 0; 
    for( int row = 1 ; row < (ROWS-1) ; ++row ){
        cnt_high = 0 ;
		for( int i = 0 ; i < 2 ; ++i ) for( int j = -1 ; j < 2 ; ++j ){
			cnt_high += (EtaChecked[col+i][row+j] - (EtaChecked[col][row]) <= 0.001) ? 1 : 0 ;
		}
        if ((cnt_high==6)&&(EtaChecked[col][row]>pb.MaxEtaInit)){
			peak[col][row] = 1;
			IndPeak[0][idd] = row;
			IndPeak[1][idd] = col;
			idd++;
        }
	}

	//right part
	col = COLS-1; 
    for( int row = 1 ; row < (ROWS-1) ; ++row ){
        cnt_high = 0 ;
		for( int i = -1 ; i < 1 ; ++i ) for( int j = -1 ; j < 2 ; ++j ){
			cnt_high += (EtaChecked[col+i][row+j] - (EtaChecked[col][row]) <= 0.001) ? 1 : 0 ;
		}
        if ((cnt_high==6)&&(EtaChecked[col][row]>pb.MaxEtaInit)){
			peak[col][row] = 1;
			IndPeak[0][idd] = row;
			IndPeak[1][idd] = col;
			idd++;
        }
	}
	
	return idd; 
}

double** fftw_to_double(fftw_complex* in,int nx,int ny)
{
	double** out = declare_2darray(nx,ny);
	for (int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			out[i][j] = in[i*nx+j][0];
		}
	}
	return out;
}

Vec2d** fun_grad(double** array, double dx, double dy, int nx, int ny)
{	
	Vec2d** grad = new Vec2d*[ny];
	for (int i=0;i<ny;i++){
		grad[i] = new Vec2d[nx];
	}
	
	#pragma omp parallel for
	for (int i=0;i<ny;i++){//for x=0 and x=end
		grad[i][0].x    = (array[i][1]-array[i][0])/(dx);
		grad[i][nx-1].x = (array[i][nx-1]-array[i][nx-2])/(dx);
		if (i==0){
			grad[i][0].y= (array[i+1][0]-array[i][0])/(dy);
			grad[i][nx-1].y= (array[i+1][nx-1]-array[i][nx-1])/(dy);
		}
		else if (i==ny-1){
			grad[i][0].y= (array[i][0]-array[i-1][0])/(dy);
			grad[i][nx-1].y= (array[i][nx-1]-array[i-1][nx-1])/(dy);
		}
		else {
			grad[i][0].y    = (array[i+1][0]-array[i-1][0])/(2.0*dy);
			grad[i][nx-1].y = (array[i+1][nx-1]-array[i-1][nx-1])/(2.0*dy);
		}
	}
	
	#pragma omp parallel for
	for (int j=0;j<nx;j++){//for y=0 and y=end
		grad[0][j].y    = (array[1][j]-array[0][j])/(dy);
		grad[ny-1][j].y = (array[ny-1][j]-array[ny-2][j])/(dy);
		if (j==0){
			grad[0][j].x= (array[0][j+1]-array[0][j])/(dx);
			grad[ny-1][j].x= (array[ny-1][j+1]-array[ny-1][j])/(dx);
		}
		else if (j==nx-1){
			grad[0][j].x= (array[0][j]-array[0][j-1])/(dx);
			grad[ny-1][j].x= (array[ny-1][j]-array[ny-1][j-1])/(dx);
		}
		else {
			grad[0][j].x    = (array[0][j+1]-array[0][j-1])/(2.0*dx);
			grad[ny-1][j].x = (array[ny-1][j+1]-array[ny-1][j-1])/(2.0*dx);
		}
	}	
	
	#pragma omp parallel for
	for (int i=1;i<ny-1;i++) {
		for (int j=1;j<nx-1;j++) {
			grad[i][j].x=(array[i][j+1]-array[i][j-1])/(2.0*dx);
			grad[i][j].y=(array[i+1][j]-array[i-1][j])/(2.0*dy);
		}
	}
	
	return grad;
}

int fun_kwadran(double ux,double uy)
{
	if ((ux>=0) && (uy>=0)){
		return 1;
	}
	else if ((ux<0) && (uy>=0)){
		return 2;
	}
	else if ((ux<0) && (uy<0)){
		return 3;
	}
	else {
		return 4;
	}
}

void breaking_criterion(double t,double* U,const double* Zhat,double** etah,fftw_complex* gradphi_x,fftw_complex* gradphi_y,int Nx,int Ny,double dx, double dy)
{
	double** PeaksChar  = declare_2darray(Nx,Ny);
	int* IndPeak[2];
	int  Npks;
	int  flag;
	IndPeak[0] = new int[Nx*Ny];
	IndPeak[1] = new int[Nx*Ny];
	
	double** EtaChecked = declare_2darray(Nx,Ny);
	#pragma omp parallel for
	for(int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			EtaChecked[i][j] = etah[i][j]*pb.Char[i][j];
		}
	}
	
	double max_eta = fun_max2darray(EtaChecked,Nx,Ny);
	
	if (max_eta>pb.MaxEtaInit) {
		Npks=fun_peaks(t,PeaksChar,IndPeak,EtaChecked,Nx,Ny);
		
		if (Npks==0){
			flag=0;
		}
		else {
			flag=1;
		}
	}
	else {
		flag = 0;
	}
	
	if (flag==1) {
		complex** HilbEta_hat = declare_2darray_complex(Nx,Ny);
		double**  HilbEta = declare_2darray(Nx,Ny);
		
		#pragma omp parallel for
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				HilbEta_hat[i][j].Re = Zhat[Nx*Ny+i*Nx+j]*fun_sign(dominfo.kx[j]+dominfo.ky[i]);
				HilbEta_hat[i][j].Im = -Zhat[i*Nx+j]*fun_sign(dominfo.kx[j]+dominfo.ky[i]);
			}
		}
		
		ifftw_2d_c2r(HilbEta,HilbEta_hat,Nx,Ny);
		
		double** gradxHilbEta = declare_2darray(Nx,Ny);
		double** gradyHilbEta = declare_2darray(Nx,Ny);
		double** gradxEta	= declare_2darray(Nx,Ny);
		double** gradyEta	= declare_2darray(Nx,Ny);
		
		gradient2D(HilbEta,dx,dy,Nx,Ny,gradxHilbEta,gradyHilbEta);		
		gradient2D(etah,dx,dy,Nx,Ny,gradxEta,gradyEta);
		
		double Kx,Ky;
		double** Kk = declare_2darray(Nx,Ny);
		
		#pragma omp parallel for
		for(int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				Kx=(etah[i][j]*gradxHilbEta[i][j]-HilbEta[i][j]*gradxEta[i][j])/(pow(etah[i][j],2)+pow(HilbEta[i][j],2));
				Ky=(etah[i][j]*gradyHilbEta[i][j]-HilbEta[i][j]*gradyEta[i][j])/(pow(etah[i][j],2)+pow(HilbEta[i][j],2));
				Kk[i][j] = sqrt(pow(Kx,2)+pow(Ky,2));
			}
		}
		
		itCB=0;
		
		int indx;
		int indy;
		double Ccrest;
		for (int j=0;j<Npks;j++){
			indx = IndPeak[0][j];
			indy = IndPeak[1][j];
			
			Ccrest = fun_exact_Cp1d_val(Kk[indy][indx],dominfo.bathy.pos[indy][indx]+etah[indy][indx],dominfo.U);
				
			if (fabs(U[indy*Nx+indx]) >= (pb.KBC*fabs(Ccrest))){	
				CB[itCB].xindex = indx;
				CB[itCB].yindex = indy;
				CB[itCB].xNow	= indx;
				CB[itCB].yNow	= indy;
				CB[itCB].Ucrest = fabs(U[indy*Nx+indx]);
				CB[itCB].time 	= t;  
				CB[itCB].theta 	= fabs(atan(gradphi_y[indy*Nx+indx][0]/gradphi_x[indy*Nx+indx][0]));
				CB[itCB].kwadran= fun_kwadran(gradphi_x[indy*Nx+indx][0],gradphi_y[indy*Nx+indx][0]);
				itCB++;
			}
		}
		
		if (itCB>0){
			CB[0].NB = itCB;
		}
		else{
			CB[0].NB = 0;
		}
		
		//if (t>ITERbdt*dominfo.dt) {
			//ITERbspdt=ITERbspdt+1;
		//}   
		
		free_2darray(Kk,Nx,Ny);	
		free_2darray(HilbEta,Nx,Ny);
		free_2darray_complex(HilbEta_hat,Nx,Ny);		
	
		free_2darray(gradxHilbEta,Nx,Ny);
		free_2darray(gradyHilbEta,Nx,Ny);
		free_2darray(gradxEta,Nx,Ny);
		free_2darray(gradyEta,Nx,Ny);
		
	}
	
	//freeing memory	
	free_2darray(EtaChecked,Nx,Ny);
	free_2darray(PeaksChar,Nx,Ny);
	delete[] (IndPeak[0]);
	delete[] (IndPeak[1]);
}

int fun_check_crest(int kk,int NCBprev)
{
	int idx,idy;
	idx = CB[kk].xindex;
	idy = CB[kk].yindex;
	
	for (int ll=0;ll<NCBprev;ll++){
		if ((CBPrev[ll].xindex==idx) && (CBPrev[ll].yindex==idy)) {
			return 0;
		}
	}
			
	for (int ll=0;ll<NCBprev;ll++){
		for (int nn=0;nn<CBPrev[ll].nnode;nn++){
			if ((CBPrev[ll].nodes[nn].x==idx) && (CBPrev[ll].nodes[nn].y==idy)) {
				return 0;
			}
		}
	}	
	
	return 1;//means new crest
}

void funBr_Char(double** Charh,double t,CBvar* CrestBreakprop,int I,int Nx,int Ny)
{
	int indxMaxPrev = CrestBreakprop[I].xNow;
	int indyMaxPrev = CrestBreakprop[I].yNow;
	double theta = CrestBreakprop[I].theta;
	int kwadran  = CrestBreakprop[I].kwadran;
	
	double xp = dominfo.x[indxMaxPrev];
	double yp = dominfo.y[indyMaxPrev];
	double xc,yc;
	
	set_matrix_val(Charh,Nx,Ny,0.0);//zeroes
	
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			xc = dominfo.x[j]-xp;
			yc = dominfo.y[i]-yp;
			if ((pow(xc,2)+pow(yc,2))<=pow(pb.Rsearch,2)){
				Charh[i][j] = 1;
			}
			else{
				Charh[i][j] = 0;
			}
		}
	}
	
	double spread = 80.0;
	double alpha,beta;
	double XXr,YYr;
	if (kwadran==1){
		alpha = theta-spread*M_PI/180.0;
		beta  = (spread*M_PI/180.0+theta)-M_PI/2.0;
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				XXr = (dominfo.x[j]-xp)*cos(alpha)+(dominfo.y[i]-yp)*sin(alpha);
				YYr = (dominfo.y[i]-yp)*cos(beta)-(dominfo.x[j]-xp)*sin(beta);
				if ((XXr<0) || (YYr<0)){
					Charh[i][j] = 0;
				}
			}
		}
	}
	else if (kwadran==2){
		alpha = spread*M_PI/180.0-theta;
		beta  = (M_PI/2.0-theta)-(spread*M_PI/180.0);
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				XXr = (dominfo.x[j]-xp)*cos(alpha)+(dominfo.y[i]-yp)*sin(alpha);
				YYr = (dominfo.y[i]-yp)*cos(beta)-(dominfo.x[j]-xp)*sin(beta);
				if ((XXr>0) || (YYr<0)){
					Charh[i][j] = 0;
				}
			}
		}
	}
    else if (kwadran==3){
		alpha = theta-spread*M_PI/180.0;
		beta  = (spread*M_PI/180.0+theta)-(M_PI/2.0);
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				XXr = (dominfo.x[j]-xp)*cos(alpha)+(dominfo.y[i]-yp)*sin(alpha);
				YYr = (dominfo.y[i]-yp)*cos(beta)-(dominfo.x[j]-xp)*sin(beta);
				if ((XXr>0) || (YYr>0)){
					Charh[i][j] = 0;
				}
			}
		}
	}
    else if (kwadran==4){
		alpha = spread*M_PI/180.0-theta;
		beta  = (M_PI/2.0-theta)-(spread*M_PI/180.0);
		for (int i=0;i<Ny;i++){
			for(int j=0;j<Nx;j++){
				XXr = (dominfo.x[j]-xp)*cos(alpha)+(dominfo.y[i]-yp)*sin(alpha);
				YYr = (dominfo.y[i]-yp)*cos(beta)-(dominfo.x[j]-xp)*sin(beta);
				if ((XXr<0) || (YYr>0)){
					Charh[i][j] = 0;
				}
			}
		}
	}	
}

void funBr_find_peak_loc(int* indexx,double t,double** eta,double** Char,int Nx,int Ny)
{
	double** EtaCheck = declare_2darray(Nx,Ny);
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			EtaCheck[i][j] = (eta[i][j])*(Char[i][j])*pb.Char[i][j];
		}
	}
	
	idx_max_matrix(indexx,EtaCheck,Nx,Ny);
	
	free_2darray(EtaCheck,Nx,Ny);
}

void funBr_BreakNodes_char(double t,double* U,double U_star,double U_F,double** eta,CBvar* CBNew,int lx,int Nx,int Ny)
{
	double** Char = declare_2darray(Nx,Ny);
	funBr_Char(Char,t,CBNew,lx,Nx,Ny);
	
	int id=0;
	int idprev=0;
	int idx,idy;
	
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			if ((U[i*Nx+j]<U_F) || (eta[i][j]<0)){
				Char[i][j] = 0;
			}
			Char[i][j] = Char[i][j]*pb.Char[i][j];
			
			if ((Char[i][j] == 1)){
				Break_nodes_prev[idprev].x = j;
				Break_nodes_prev[idprev].y = i;	
				idprev++;
				if (B[i][j] == 0){				
					//filling B
					if (U[i*Nx+j]>=(2*U_star)){
						B[i][j] = 1;
					}
					else if (U[i*Nx+j]<=U_star){
						B[i][j] = 0;
					}
					else {
						B[i][j] = ((U[i*Nx+j]/U_star)-1);
					}
					
					if (B[i][j]>0){
						//get Break_nodes when B>0
						Break_nodes[id].x = j;
						Break_nodes[id].y = i;									
						id++;
					}
				}
			}
		}
	}
	
	nbreak_nodes = id;
	nbreak_prev  = idprev;
	
	free_2darray(Char,Nx,Ny);
}

void breaking_process(double t,const double* Zhat,fftw_complex* eta,fftw_complex* gradphi_x,fftw_complex* gradphi_y,int fact)
{//eta u an v no scaling factIFFT
	Nx=dominfo.Nx;
	Ny=dominfo.Ny;
	dx=dominfo.dx;
	dy=dominfo.dy;
	double*  U=new double[Nx*Ny];
	double** etah       = declare_2darray(Nx,Ny);
	
	for (int i=0;i<Ny;i++){
		for(int j=0;j<Nx;j++){
			etah[i][j]= eta[i*Nx+j][0]/fact;
			U[i*Nx+j] = sqrt(pow(gradphi_x[i*Nx+j][0]/fact,2)+pow(gradphi_y[i*Nx+j][0]/fact,2));
		}
	}
	
	CB[0].NB = 0;//set clear CB
	set_matrix_val(B,Nx,Ny,0.0);
	breaking_criterion(t,U,Zhat,etah,gradphi_x,gradphi_y,Nx,Ny,dx,dy);//find candidat for CB
	int NCB = CB[0].NB;
	
	if ((NCB > 0) && (flagbr == 0)) {//new first crest break	
		flagbr = 1;
		for (int j=0;j<NCB;j++){
			CBNew[j].xindex = CB[j].xindex;
			CBNew[j].yindex = CB[j].yindex;
			CBNew[j].xNow	= CB[j].xNow;
			CBNew[j].yNow	= CB[j].yNow;
			CBNew[j].Ucrest = CB[j].Ucrest;
			CBNew[j].theta	= CB[j].theta;
			CBNew[j].time   = CB[j].time;
			CBNew[j].NB   	= CB[j].NB;
			CBNew[j].kwadran= CB[j].kwadran;
		}
	}
	else if ((NCB == 0) && (flagbr == 2)) {//only from previous crest break	
		int NCBprev = CBPrev[0].NB;
		for (int j=0;j<NCBprev;j++){
			CBNew[j].xindex = CBPrev[j].xindex;
			CBNew[j].yindex = CBPrev[j].yindex;
			CBNew[j].xNow	= CBPrev[j].xNow;
			CBNew[j].yNow	= CBPrev[j].yNow;
			CBNew[j].Ucrest = CBPrev[j].Ucrest;
			CBNew[j].theta	= CBPrev[j].theta;
			CBNew[j].time   = CBPrev[j].time;
			CBNew[j].NB   	= CBPrev[j].NB;
			CBNew[j].kwadran= CBPrev[j].kwadran;
		}	
		CBNew[0].NB = NCBprev;
	}
	else if ((NCB > 0) && (flagbr == 2)) {//previous and new crest break ( new)
		int NCBprev	= CBPrev[0].NB;
		
		int     dleft;
		int 	add,IdCheckNewOld;				
		int 	NCBadd = NCBprev+NCB;

		for (int j=0;j<NCBprev;j++){
			CBNew[j].xindex = CBPrev[j].xindex;
			CBNew[j].yindex = CBPrev[j].yindex;
			CBNew[j].xNow	= CBPrev[j].xNow;
			CBNew[j].yNow	= CBPrev[j].yNow;
			CBNew[j].Ucrest = CBPrev[j].Ucrest;
			CBNew[j].theta	= CBPrev[j].theta;
			CBNew[j].time   = CBPrev[j].time;
			CBNew[j].NB   	= CBPrev[j].NB;
			CBNew[j].kwadran= CBPrev[j].kwadran;
		}	
		
		add = 0;
		for (int k=0;k<NCB;k++) {
			IdCheckNewOld = fun_check_crest(k,NCBprev);//new break event or not
			
			if (IdCheckNewOld == 1) { //this is for a new break event
				CBNew[NCBprev+add].xindex  =CB[k].xindex;
				CBNew[NCBprev+add].yindex  =CB[k].yindex;
				CBNew[NCBprev+add].xNow    =CB[k].xNow;
				CBNew[NCBprev+add].yNow    =CB[k].yNow;
				CBNew[NCBprev+add].Ucrest  =CB[k].Ucrest;
				CBNew[NCBprev+add].time    =CB[k].time;
				CBNew[NCBprev+add].theta   =CB[k].theta;
				CBNew[NCBprev+add].kwadran =CB[k].kwadran;
				add=add+1;
			}
		}
		
		//update number of crest break
		CBNew[0].NB = NCBprev+add;
	}
	
	//IF NCB exist
	if ((flagbr==1) || (flagbr==2)) {
		
		CBPrev[0].NB = 0;
		int 	NCBnow 	    = CBNew[0].NB;	
		int* 	IDTerminate = (int*) malloc(sizeof(int)*NCBnow);
		int iterNprev_prev=0;
		iterNprev=0;
				
		for (int lx=0;lx<NCBnow;lx++){
			IDTerminate[lx]=0;
			if (flagbr==1){
				CBNew[lx].xNow = CBNew[lx].xindex;
				CBNew[lx].yNow = CBNew[lx].yindex;
			}
			else {					
				double** Charh  = declare_2darray(Nx,Ny);
				int* indxMaxNow = new int[2];
					
				funBr_Char(Charh,t,CBNew,lx,Nx,Ny);								
				funBr_find_peak_loc(indxMaxNow,t,etah,Charh,Nx,Ny);
				
				if (indxMaxNow[0]<0){
					IDTerminate[lx]=1;//breaking finish 
					free_2darray(Charh,Nx,Ny);
					delete[] indxMaxNow;
					continue; 
				}				
				CBNew[lx].xNow = indxMaxNow[0];
				CBNew[lx].yNow = indxMaxNow[1];											
				CBNew[lx].theta= fabs(atan(gradphi_y[indxMaxNow[1]*Nx+indxMaxNow[0]][0]/gradphi_x[indxMaxNow[1]*Nx+indxMaxNow[0]][0]));
				CBNew[lx].kwadran= fun_kwadran((gradphi_x[indxMaxNow[1]*Nx+indxMaxNow[0]][0]),(gradphi_y[indxMaxNow[1]*Nx+indxMaxNow[0]][0]));

				free_2darray(Charh,Nx,Ny);
				delete[] indxMaxNow;
			}			
			
			int 	it_indexb = 0;	
			double 	U_I,U_F,U_star;
			int 	break_index,indEND;	
			int 	lenb;
			
			tbreak 	= CBNew[lx].time;
			U_I		= CBNew[lx].Ucrest;
			U_F 	= pb.TC*U_I;
			
			if  ((t-tbreak) >= pb.Tstar) {
				U_star = U_F;
			}
			else if (((t-tbreak) >= 0) && (t-tbreak) < pb.Tstar) {
				U_star = U_I+((t-tbreak)/pb.Tstar)*(U_F-U_I);
			}
			else {
				U_star = U_I;
			}
			
			funBr_BreakNodes_char(t,U,U_star,U_F,etah,CBNew,lx,Nx,Ny);//filling matrix B, etah udah scaled	
			
			/*int idx,idy;
			for (int i=0;i<nbreak_nodes;i++){//nbreak nodes per peak
				idx = Break_nodes[i].x;
				idy = Break_nodes[i].y;
				
				if ((U[idy*Nx+idx]>=U_star)){
					IDTerminate[lx]=0;//continue breaking
					break;
				}
				else {
					IDTerminate[lx]=1;//breaking end for lx
				}
			}*/
			
			//all break nodes
			if (nbreak_prev>0){//even B=0
				for (int jj=0;jj<nbreak_prev;jj++){
					CBPrev[iterNprev_prev].nodes_prev[jj].x=Break_nodes_prev[jj].x;
					CBPrev[iterNprev_prev].nodes_prev[jj].y=Break_nodes_prev[jj].y;
				}
				CBPrev[iterNprev_prev].nnode_prev  =nbreak_prev;
				iterNprev_prev++;				
			}
			
			//define break nodes that continue breaking
			if (nbreak_nodes==0){
				IDTerminate[lx]=1;//breaking end for lx
			}
			else{
				IDTerminate[lx]=0;//continue breaking
				for (int jj=0;jj<nbreak_nodes;jj++){
					CBPrev[iterNprev].nodes[jj].x=Break_nodes[jj].x;
					CBPrev[iterNprev].nodes[jj].y=Break_nodes[jj].y;
				}
				CBPrev[iterNprev].nnode  =nbreak_nodes;
				CBPrev[iterNprev].xindex =CBNew[lx].xindex;
				CBPrev[iterNprev].yindex =CBNew[lx].yindex;
				CBPrev[iterNprev].xNow	 =CBNew[lx].xNow;
				CBPrev[iterNprev].yNow	 =CBNew[lx].yNow;
				CBPrev[iterNprev].Ucrest =CBNew[lx].Ucrest;
				CBPrev[iterNprev].time	 =CBNew[lx].time;
				CBPrev[iterNprev].theta	 =CBNew[lx].theta;
				CBPrev[iterNprev].kwadran=CBNew[lx].kwadran;
				iterNprev=iterNprev+1;
			}
			CBPrev[0].NB=iterNprev;//update NBprev CB can be more or less 
		}		

		for (int lx=0;lx<NCBnow;lx++){
			if (IDTerminate[lx] == 0){//still any break not end yet
				flagbr = 2;
				break;
			}
			else{
				flagbr = 0;
			}
		}
		
		//writing break nodes
		if ((t>dominfo.t[0]+ITERbdt*dominfo.dt) && (iterNprev_prev>0)){
			fprintf(fbreak,"%lf ",t);
			for (int lx=0;lx<iterNprev_prev;lx++){
				for (int jb=0;jb<CBPrev[lx].nnode_prev;jb++){
					fprintf(fbreak," %d %d ",CBPrev[lx].nodes_prev[jb].x,CBPrev[lx].nodes_prev[jb].y);
				}
				fflush(fbreak);
			}
			fprintf(fbreak,"\n");
		}
			
		free(IDTerminate);			
	}
	
    if (t>dominfo.t[0]+ITERbdt*dominfo.dt){
		ITERbdt=ITERbdt+1;//so no double output in the same time
	}
	
	//freeing memory
	delete[] U;
	free_2darray(etah,Nx,Ny);
}
