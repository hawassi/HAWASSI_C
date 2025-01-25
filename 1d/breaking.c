void parambreak(double* x, double Xinflux, double parKBC, double parTC, double parTchar) 
{	
    pb.KBC  = parKBC;
    pb.TC   = parTC;
    pb.delb = delb;

    if (strcmp(initial, "zero") != 0){
        pb.Tstar = 5 * sqrt(depthinf * g);//based on Kennedy 
    }
    else {
        pb.Tstar = parTchar * Tp_peak;
    }

    if (strcmp(propagation, "Uni+") == 0) {
        pb.xbstart 	= closest(x, Xinflux + 0.1 * adjcoef * lambda_p, Nx);
        pb.xbend	= Nx - 1;
        lenIF 		= pb.xbend - pb.xbstart;	
    } else if (strcmp(propagation, "Uni-") == 0) {
        pb.xbend	= closest(x, Xinflux - 0.1 * adjcoef * lambda_p, Nx);
        pb.xbstart	= 0;
        lenIF 		= pb.xbend - pb.xbstart;	
    }else if (strcmp(propagation, "Bi") == 0) {
        pb.xbstart 	= 0;
        pb.xbendL       = closest(x, Xinflux - 0.1 * adjcoef * lambda_p, Nx);
        pb.xbstartR     = closest(x, Xinflux + 0.1 * adjcoef * lambda_p, Nx);
        pb.xbend  	= Nx - 1;
        ILL 		= pb.xbendL - pb.xbstart;
        IFF 		= pb.xbend - pb.xbstartR;
    }else if (strcmp(propagation, "None") == 0) {//Init Value Problem only
        pb.xbstart	= 0;
        pb.xbend	= Nx - 1;
        lenIF 		= pb.xbend - pb.xbstart;
    }
}

int* findpeaks(double* array, int lenIF) 
{
    int   j = 1;
    int*  locs_temp = (int*) malloc(sizeof(int) * Nx);
    int   Nhere = 0;

    while (j < lenIF - 1) {
        if ((array[j]) < MaxEtaInit) {
            j++;
        }
        else {
            if ((array[j] > array[j + 1]) && (array[j] > array[j - 1])) {				
                locs_temp[Nhere + 1] = j;
                Nhere = Nhere + 1;
                j++;
            }
            else {				
	        j++;
            }
        }
    }
    
    locs_temp[0] = Nhere;
    
    return locs_temp;
}

void breaking_criterion(CB* CrestBreak, double time, double* eta, complex* etahat, double* k, double *u, parbreak pb)
{
    int     flaghere;	
    int     *indj     = new int[Nx];	
    double  *eta_temp = new double[lenIF];
    for (int i = 0; i < lenIF; i++) {
        eta_temp[i] = eta[pb.xbstart + i];
    }		

    if (valueinarraylarger(MaxEtaInit, eta_temp, lenIF) == 1) {
        if ((pb.xbendL == 0) && (pb.xbstartR == 0)) {	
            int* locsplus = (int*) malloc(sizeof(int) * Nx);				
            locsplus      = findpeaks(eta_temp, lenIF);			
            Npks 	  = locsplus[0];
            if (Npks > 0){				
                for (int i = 0; i < Npks; i++) {
	            indj[i] = locsplus[i + 1] + pb.xbstart;
                }
            }
            free(locsplus);
        }
	else { //bi-directional			
            int     *locs1, *locs2;            
            locs1   = (int*) malloc(sizeof(int) * Nx);
            locs2   = (int*) malloc(sizeof(int) * Nx);	
            
            double  *eta_temp1 = new double[ILL];
            double  *eta_temp2 = new double[IFF];
            
            for (int i = 0; i < ILL; i++) {
	        eta_temp1[i] = eta[pb.xbstart + i];
            }
            for (int i = 0; i < IFF; i++) {
	        eta_temp2[i] = eta[pb.xbstartR + i];
            }			
            
            locs1 = findpeaks(eta_temp1, ILL);
            locs2 = findpeaks(eta_temp2, IFF);
            Npks  = locs1[0] + locs2[0];	
            
            if (Npks > 0){					
                for (int i = 0; i < (locs1[0]); i++) {
	            indj[i] = locs1[i + 1] + pb.xbstart;
                }
                for (int i = locs1[0]; i < (Npks); i++) {
	            indj[i] = locs2[i - locs1[0] + 1] + pb.xbstartR;
                }				
            }
            
            free(locs1);
            free(locs2);
            delete[] eta_temp1;
            delete[] eta_temp2;
	}
	
        if (Npks == 0) {
	    flaghere = 0;//no breaking
        }
        else {
	    flaghere = 1;
        }
    }
    else {
        flaghere = 0;//no breaking
    }	

    if (flaghere == 1) {	
        double*  HxRe  	 = (double*) malloc(sizeof(double) * Nx);
        complex* Hx_hat  = (complex*) malloc(sizeof(complex) * Nx);
        
        for (int i = 0; i < Nx; i++){
            if (k[i] == 0){
                Hx_hat[i].Re = 0;
                Hx_hat[i].Im = 0;
            }
            else {
                Hx_hat[i].Re = (k[i] / fabs(k[i])) * etahat[i].Im;
                Hx_hat[i].Im = -(k[i] / fabs(k[i])) * etahat[i].Re;
            }
        }
        
        ifft_complex_real_nx(Hx_hat, HxRe);
        free(Hx_hat); 
        
        double* grad_Hx  = (double*) malloc(sizeof(double) * Nx);
        double* grad_eta = (double*) malloc(sizeof(double) * Nx);
        gradient(HxRe, dx, Nx, grad_Hx);
        gradient(eta, dx, Nx, grad_eta);		

        double  KxJ, Ccrest, Upart;
        int 	indJ;
        int 	itCB = 0; 
        
        CB* CrestBreaktemp   = new CB[Npks];
        for (int i = 0; i < Npks; i++){
            indJ 	= indj[i];
            KxJ  	= (eta[indJ] * grad_Hx[indJ] - HxRe[indJ] * grad_eta[indJ]) / (pow(eta[indJ], 2) + pow(HxRe[indJ], 2));
            Ccrest      = phase_velocity_i(KxJ, depth[indJ] + eta[indJ]);
            Upart 	= u[indJ];
          
            if (fabs(Upart) > (pb.KBC * fabs(Ccrest))) {
                CrestBreaktemp[itCB].Bindex   = indJ;
                CrestBreaktemp[itCB].BindexNow= indJ;
                CrestBreaktemp[itCB].Ucrest   = fabs(Upart);
                CrestBreaktemp[itCB].DirProp  = Upart / fabs(Upart);
                CrestBreaktemp[itCB].Bposition= x[indJ];
                CrestBreaktemp[itCB].Btime    = time;				
                itCB = itCB+1;
            }
        }
        
        free(grad_Hx);
        free(grad_eta);
        free(HxRe);
        
        if (itCB > 0){
            for (int j = 0; j < itCB; j++){
                CrestBreak[j].Bindex   = CrestBreaktemp[j].Bindex;
                CrestBreak[j].BindexNow= CrestBreaktemp[j].BindexNow;
                CrestBreak[j].Ucrest   = CrestBreaktemp[j].Ucrest;
                CrestBreak[j].DirProp  = CrestBreaktemp[j].DirProp;
                CrestBreak[j].Bposition= CrestBreaktemp[j].Bposition;
                CrestBreak[j].Btime    = CrestBreaktemp[j].Btime;
                CrestBreak[j].Bcount   = itCB;
            }
        }
        
        delete[] CrestBreaktemp;
    }//end flaghere==1	
    else {//flaghere==0
        set_init_CB(CrestBreak, (Nx / 4));	
    }
    
    delete[] eta_temp;
    delete[] indj;	
}

void breaking_process(double* eta, complex* etahat, double* u, double t, double* x, parbreak pb)
{	
    double* abs_u  = (double*) malloc (sizeof(double) * Nx);
    
    for (int j = 0; j < Nx; j++) {
        abs_u[j] = fabs(u[j]);
    }
    
    set_init_double(B, Nx);	

    breaking_criterion(CrestBreak, t, eta, etahat, k, u, pb); //detect peak KBC CrestBreak data	

    int NCB = CrestBreak[0].Bcount;
	    
    if ((NCB > 0) && (flag == 0)) {//new first crest break	
        flag = 1;
        for (int j = 0; j < NCB; j++){
            CrestBreakNew[j].Bindex   = CrestBreak[j].Bindex;
            CrestBreakNew[j].BindexNow= CrestBreak[j].BindexNow;
            CrestBreakNew[j].Ucrest   = CrestBreak[j].Ucrest;
            CrestBreakNew[j].DirProp  = CrestBreak[j].DirProp;
            CrestBreakNew[j].Bposition= CrestBreak[j].Bposition;
            CrestBreakNew[j].Btime    = CrestBreak[j].Btime;
            CrestBreakNew[j].Bcount   = CrestBreak[j].Bcount;
        }
    } else if ((NCB == 0) && (flag == 2)) {//only from previous crest break
        NCBprev = CrestBreakPrev[0].Bcount;
        for (int j = 0; j < NCBprev; j++){
            CrestBreakNew[j].Bindex   = CrestBreakPrev[j].Bindex;
            CrestBreakNew[j].BindexNow= CrestBreakPrev[j].BindexNow;
            CrestBreakNew[j].Ucrest   = CrestBreakPrev[j].Ucrest;
            CrestBreakNew[j].DirProp  = CrestBreakPrev[j].DirProp;
            CrestBreakNew[j].Bposition= CrestBreakPrev[j].Bposition;
            CrestBreakNew[j].Btime    = CrestBreakPrev[j].Btime;
            CrestBreakNew[j].Bcount   = CrestBreakPrev[j].Bcount;
        }	
    } else if ((NCB > 0) && (flag == 2)) {//previous and still any crest break
        NCB	  = CrestBreak[0].Bcount;
        NCBprev	  = CrestBreakPrev[0].Bcount;
        
        int dleft, add, IdCheckNewOld;				
        int NCBadd = NCBprev + NCB;
	    
        for (int j = 0; j < NCBprev; j++){
            CrestBreakNew[j].Bindex   = CrestBreakPrev[j].Bindex;
            CrestBreakNew[j].BindexNow= CrestBreakPrev[j].BindexNow;
            CrestBreakNew[j].Ucrest   = CrestBreakPrev[j].Ucrest;
            CrestBreakNew[j].DirProp  = CrestBreakPrev[j].DirProp;
            CrestBreakNew[j].Bposition= CrestBreakPrev[j].Bposition;
            CrestBreakNew[j].Btime    = CrestBreakPrev[j].Btime;
            CrestBreakNew[j].Bcount   = CrestBreakPrev[j].Bcount;
        }					
	    
        add = 0;		
        for (int k = 0; k < NCB; k++) { 
            IdCheckNewOld = 1;
            for (int l = 0; l < NCBprev; l++) {
                dleft = abs(ID_Node_end[l] - CrestBreakPrev[l].BindexNow);           
                if  (CrestBreakPrev[l].DirProp == 1) {
                    if ((CrestBreak[k].Bindex >= (CrestBreakPrev[l].BindexNow - 3*dleft)) && (CrestBreak[k].Bindex <= ID_Node_end[l])){
                        IdCheckNewOld = 0; //not new
                        break;
                    }
                }
                else {
                    if ((CrestBreak[k].Bindex >= CrestBreakPrev[l].BindexNow) && (CrestBreak[k].Bindex <= (ID_Node_end[l] + 3 * dleft))) {
                        IdCheckNewOld = 0; //not new
                        break;
                    }                
                }
            }
            if (IdCheckNewOld == 1) { //this is for a new break
                CrestBreakNew[NCBprev+add].Bindex   = CrestBreak[k].Bindex;
                CrestBreakNew[NCBprev+add].BindexNow= CrestBreak[k].Bindex;
                ID_Node_end[NCBprev+add]     	    = CrestBreak[k].Bindex; 
                CrestBreakNew[NCBprev+add].Ucrest   = CrestBreak[k].Ucrest;
                CrestBreakNew[NCBprev+add].Bposition= CrestBreak[k].Bposition;
                CrestBreakNew[NCBprev+add].Btime    = CrestBreak[k].Btime;
                CrestBreakNew[NCBprev+add].DirProp  = CrestBreak[k].DirProp;
                CrestBreakNew[NCBprev+add].Bcount   = CrestBreak[k].Bcount;
                add = add+1;
            }
        }//end of looping k
	    
	    //update number of crest break
        int j = 0;
        while (CrestBreakNew[j].Bindex > 0) {
            j++;
        }
        for (int i = 0; i < j; i++){
	    CrestBreakNew[i].Bcount = j;
        }
    }

    //IF NCB exist
    if ((flag == 1) || (flag == 2)) {
        int 	NCBnow 	    = CrestBreakNew[0].Bcount;	
        int 	it_indexb   = 0;		
        int* 	IDTerminate = (int*) malloc(sizeof(int) * NCBnow);
        double 	tbreak, U_I, U_F, U_star, U_FF;
        int 	break_index, indEND;	
        int 	lenb;

        iterNprev = 0;
        set_init_int(IDTerminate, NCBnow);	
        set_init_CB(CrestBreakPrev, (Nx / 4));
	    
        for (int lx = 0; lx < NCBnow; lx++) {	
            if (flag == 1) {
	        indxMaxNow = CrestBreakNew[lx].Bindex;
            }
            else {
                int  	indxMaxStart, indxMaxLoc, lenloc;
                double 	*etaloc, maxetaloc;
                
                indxMaxStart  = CrestBreakNew[lx].BindexNow; //peak from previous
                lenloc 	      = abs(ID_Node_end[lx] - indxMaxStart); 
                etaloc 	      = (double*) malloc(sizeof(double) * lenloc); 
                
                for (int j = 0; j < lenloc; j++){
	            etaloc[j] = eta[indxMaxStart + j * CrestBreakNew[lx].DirProp];
                }	
				                
                indxMaxLoc   = idx_max(etaloc, lenloc);			
                indxMaxNow   = indxMaxStart + (indxMaxLoc) * CrestBreakNew[lx].DirProp;	
                
                free(etaloc);
            }
            
	    CrestBreakNew[lx].BindexNow = indxMaxNow;
	        
            if (lx != 0) {
                if (indxMaxNow == CrestBreakNew[lx-1].BindexNow) {
                    IDTerminate[lx] = 0;   //terminate CrestBreak from previous one
                    continue;
                }
            }
            
            tbreak 	= CrestBreakNew[lx].Btime;
            U_I		= CrestBreakNew[lx].Ucrest;
            U_F 	= pb.TC * U_I;
            
            if ((t - tbreak) >= pb.Tstar) {
	        U_star = U_F;
            } else if (((t - tbreak) >= 0) && (t - tbreak) < pb.Tstar) {
	        U_star = U_I + ((t - tbreak) / pb.Tstar) * (U_F - U_I);
            }
            else {
	        U_star = U_I;
            }			
            
            if (CrestBreakNew[lx].DirProp == 1) {
	        indEND = pb.xbend;
            }
            else {
	        indEND = pb.xbstart;
            }
            
            Node_end 	  = indxMaxNow;
            lenb 	  = (indEND - indxMaxNow) / CrestBreakNew[lx].DirProp;							
	        
	    int flaghere = 0;			
            for (int ly = indxMaxNow; ly < indEND; ly += CrestBreakNew[lx].DirProp) {
                U_FF 	 = U_F;
                
                if (abs_u[ly] < U_FF ) {
                    Node_end = ly;
                    break;
                }
                
                break_index = ly;
                
                if (abs_u[break_index] > (2 * U_star)) {
                    B[break_index]  =  1;
                } else if (abs_u[break_index] <= U_star) {
                    B[break_index]  =  0;
                }
                else {
                    B[break_index]  = (abs_u[break_index] / U_star ) - 1;
                }
                
                if ((abs_u[break_index] >= U_star)){
                    flaghere = 1;
                }
                
                Break_nodes[it_indexb] = break_index;	
                it_indexb = it_indexb + 1;				
            }

            //define the IDTerminate 
            if (flaghere == 0) {
	        IDTerminate[lx] = 0;//breaking finish
            }
            else { 
                IDTerminate[lx]                         = 1;
                ID_Node_end[iterNprev]			= Node_end;
                CrestBreakPrev[iterNprev].Bindex	= CrestBreakNew[lx].Bindex;
                CrestBreakPrev[iterNprev].BindexNow	= CrestBreakNew[lx].BindexNow;
                CrestBreakPrev[iterNprev].Ucrest	= CrestBreakNew[lx].Ucrest;
                CrestBreakPrev[iterNprev].Bposition	= CrestBreakNew[lx].Bposition;
                CrestBreakPrev[iterNprev].Btime		= CrestBreakNew[lx].Btime;
                CrestBreakPrev[iterNprev].DirProp	= CrestBreakNew[lx].DirProp;
                CrestBreakPrev[iterNprev].Bcount	= CrestBreakNew[lx].Bcount;
                iterNprev                               = iterNprev + 1;				
            }				
        }//end for lx
	    
        //update CrestBreakPrev.Bcount
        int j = 0;
        while (CrestBreakPrev[j].Bindex > 0) {
            j++;
        }
        for (int i = 0; i < j; i++){
	    CrestBreakPrev[i].Bcount = j;
        }
        
        for (int i = 0; i < NCBnow; i++){
            if (IDTerminate[i] == 1){
                flag = 2;
                break;
            } 
            else {
	        flag = 0;
            }
        }				

        //allocation for data break
        if (t > (ITERbdt * dt)) {  			
            if (it_indexb > 0) {
                dataBreak_nodes[0] = t;
                dataBreak_nodes[1] = it_indexb;
                for (int j = 2; j < (it_indexb + 2); j++) {
	            dataBreak_nodes[j] = Break_nodes[j-2];
                }
                data_saving_break(fbreak, dataBreak_nodes, Nx);
                ITERbn = ITERbn + 1;
            }			
        }
        
        free(IDTerminate);
    }	

    if (t > (ITERbdt * dt)) {
        ITERbdt = ITERbdt + 1;
    }
    
    set_init_CB(CrestBreakNew, (Nx / 4));
    set_init_CB(CrestBreak, (Nx / 4));
    
    free(abs_u);
}
