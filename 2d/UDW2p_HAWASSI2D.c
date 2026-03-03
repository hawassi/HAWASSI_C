/*#include <stdio.h>
#include <string.h>
extern char filename[256];
extern double seainfo.Hs;
extern double seainfo.wp;
extern char bath[20];*/

void Write_HOS_Ocean_Comment(FILE* fd)
{
	double Hs_in   = seainfo.Hs;
	double Tp_peak = 2*M_PI/seainfo.wp;
	fprintf(fd,"######################################\n");    
    fprintf(fd,"#                                     \n");
    fprintf(fd,"#  HAWASSI 2D output (no calibration)\n");
    fprintf(fd,"#  Case Hs=%f, Tp= %.4f   \n",Hs_in,Tp_peak);
    fprintf(fd,"#                                     \n");
    fprintf(fd,"######################################\n");
}

void UDW2p_writeHeader_hw2D(int headerSize,int Nx, int Ny,int Np, int Nt, double dx, double dy, double dt,
double xScale, double tScale, double** depthout)
{
	FILE* fd = fopen(filename,"w");
	Write_HOS_Ocean_Comment(fd);
	fclose(fd);

	fd = fopen(filename,"a");
	fprintf(fd,"PRECISION = single\n");
	fprintf(fd,"Software = HAWASSI\n");
	fprintf(fd,"Nheader = %d\n", headerSize);
	fprintf(fd,"Nx = %d\n", Nx);
	fprintf(fd,"Ny = %d\n", Ny);
	fprintf(fd,"Np = %d\n", Np);
	fprintf(fd,"Nt = %d\n", Nt);
	fprintf(fd,"dx = %-16.12E\n", dx/xScale);
	fprintf(fd,"dy = %-16.12E\n", dy/xScale);
	fprintf(fd,"dt = %-16.12E\n", dt/tScale);
	//fprintf(fd,"xCalib = %-16.12f\n", xCalibration);
	fprintf(fd,"xScale = %-16.12f\n", xScale);
	fprintf(fd,"tScale = %-16.12f\n", tScale);
	fprintf(fd,"WD = %-16.12f\n", xScale);
	if (strcmp(bath,"flat")!=0){
		fprintf(fd,"Bathymetry = variable\n");
		for (int i=0;i<Ny;i++){
			for (int j=0;j<Nx;j++){
				fprintf(fd,"%-10.8f ", depthout[i][j]/xScale);
			}
		}
	}
	fclose(fd);
}

void UDW2p_setUp_hw2D(int Nx, int Ny, int Np, int Nt, double dx, double dy, double dt, double xScale, double tScale, double** depthout)//add fileName argument
{
	double potScale = xScale*xScale/tScale;
	double velScale = xScale/tScale;
	int nheader  = 1;
	int recSize = (4*(Nx*2+Np*2)*Ny) + (4*(Ny*2+Np*2)*Nx);
	int sz,residu;
	
	UDW2p_writeHeader_hw2D(nheader,Nx,Ny,Np,Nt,dx,dy,dt,xScale,tScale,depthout);
	
	FILE* fd = fopen(filename,"a");
	fseek(fd, 0L, SEEK_END);
	sz = ftell(fd); 
	nheader = (sz/recSize)+1;
	
	if (nheader>1){
		UDW2p_writeHeader_hw2D(nheader,Nx,Ny,Np,Nt,dx,dy,dt,xScale,tScale,depthout);
	}
	
	while (ftell(fd)<(recSize*nheader)){
		fprintf(fd," ");
	}
	fclose(fd);	
}

//writing output every time step dt
void UDW2p_write_hw2D(char* filename,int Nx,int Ny,int Np,double **elev,double **fphi,
double **leftPhi, double **rightPhi, double **frontPhi,double **backPhi)
{	
	FILE* fd = fopen(filename,"a");
	float var;
	
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Nx;j++){
			var = float (elev[i][j]);
			fwrite(&var,sizeof(float),1,fd);//save elevation
		}
	}
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Nx;j++){
			var = float (fphi[i][j]);
			fwrite(&var,sizeof(float),1,fd);//save surface potential
		}
	}
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Np;j++){
			var = float (leftPhi[i][j]);
			fwrite(&var,sizeof(float),1,fd);
		}
	}
	for (int i=0;i<Ny;i++){
		for (int j=0;j<Np;j++){
			var = float (rightPhi[i][j]);
			fwrite(&var,sizeof(float),1,fd);
		}
	}
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Np;j++){
			var = float (frontPhi[i][j]);
			fwrite(&var,sizeof(float),1,fd);
		}
	}
	for (int i=0;i<Nx;i++){
		for (int j=0;j<Np;j++){
			var = float (backPhi[i][j]);
			fwrite(&var,sizeof(float),1,fd);
		}
	}
    
	fclose(fd);	
}

void getGQPoints(double* GQ20) 
//SPG: 20-point Gauss quadrature: sampling points and weight in the normalized interval (-1, 1) 
//Because of symmetry, only sampling points in the half internal (positive) are tabulated.
{
	double* spg = new double[10];
	spg[0] = 0.076526521133497333755;
	spg[1] = 0.227785851141645078080;
	spg[2] = 0.373706088715419560673;
	spg[3] = 0.510867001950827098004;
	spg[4] = 0.636053680726515025453;
	spg[5] = 0.746331906460150792614;
	spg[6] = 0.839116971822218823395;
	spg[7] = 0.912234428251325905868;
	spg[8] = 0.963971927277913791268;
	spg[9] = 0.993128599185094924786;

	for (int isp=0;isp<10;isp++){
		GQ20[isp] = 0.5-0.5*spg[9-isp];
		GQ20[10+isp] = 0.5+0.5*spg[isp];
	}	  
}
