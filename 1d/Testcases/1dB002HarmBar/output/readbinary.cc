#include <iostream>
#include <cassert>
#include <cstdlib>
#include <filesystem>



int main()
{	int Nxdata;
	int Ntdata;
	double* xdata;
	double* tdata;
	double Nx,Nt;
	char* name_constant = new char[256];
	sprintf(name_constant,"Constants_01.dat");
	FILE* fc = fopen(name_constant,"rb");
	if (fc==NULL){printf("File constant unable to open...\n");}
	else {fc = fopen(name_constant,"rb");
		fread(&Nx,sizeof(double),1,fc);Nxdata=(int) Nx;
		xdata = new double[Nxdata];
		fread(xdata,sizeof(double),Nxdata,fc);
		for(int j = 0; j<=Nxdata-1 ;j++){
    		printf("%f\n",xdata[j]);}
		fread(&Nt,sizeof(double),1,fc);Ntdata=(int) Nt;
		tdata = new double[Ntdata];
		fread(tdata,sizeof(double),Ntdata,fc);
		fclose(fc);}
	delete[] name_constant;
	if (Nxdata != Nx){printf("Nx=%f,Ninput=%d\n",Nx,Nxdata);
		printf("Number of NX should be equal!!\n");exit(0);}
	printf("Nx=%d\n",Nxdata);
	printf("Nt=%d\n",Ntdata);
	double* etadata=new double[Nxdata*Ntdata];
	//open file elevation
	char* filename = new char[256];
	sprintf(filename, "Hawassi_eta_01.dat");
	FILE* feta = fopen(filename,"rb");
        if (feta == NULL){printf("File elevation unable to open...\n");}
        else{printf("File elevation opened...\n");
		fread(etadata,sizeof(double),Ntdata*Nxdata,feta);}
	FILE *fptr = fopen("output.txt", "w");
    for (int i = 0; i<Nxdata ;i++){fprintf(fptr,"%f	%f\n",xdata[i],etadata[3815*Nxdata+i]);} 
    fclose(fptr);
	fclose(feta);
	delete[] filename;
	return 0;}