#include"pressure.h"

void pressure_correction_east()
{
	double d;	/* rho*dx */
  	double dpxel, uel, apue;	/* For the interpolation work */
  	double  dpxe;	/* Rhie chow interpolated cell face velocities */ 
  	int ie;
	double rhoe;
	
	double axBFe;
	
	double xGu,xMj,xChoi;
	
	double volE;
  
  
  	for (i=1;i<=nxm;i++)
	  {
		ieast = i+1;
		iwest = i-1;
		ieasteast = i+2;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;
		
		if(i==nxm&&RightBoundType!=4)	continue;
		if(i==1&&LeftBoundType!=4)	continue;
	
	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	{
	      			/* Since, volume of the east cell equal to that of the central face, */
	      			/* interpolated cel face quantites can be written as follows */
	      
		      		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
	      			vole = dxpe * s;
	      			
	      			vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
	      			volE = (xf[ieast]-xf[ieast-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
	      	
				/* Arithmetic Interpolation of density*/
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];
	      			d = rhoe*s;
	      
	      			
	      			//Interpolation for Original Rhie Chow
	      			if(i==1)
	      			{
	      				//Extrapolation of Rhie Chow from interior cells  
	      				dpxel = 1.5*dpx[ieast][j][k]*apu[ieast][j][k]   - 0.5*dpx[ieasteast][j][k]*apu[ieasteast][j][k] ;
	      				uel   = 1.5*u[ieast][j][k][l] - 0.5*u[ieasteast][j][k][l];
					apue  = 1.5*apu[ieast][j][k]  - 0.5*apu[ieasteast][j][k];     				
	      			}
	      			
	      			else if(i==nxm) 
	      			{
	      				dpxel = 1.5*dpx[i][j][k]*apu[i][j][k]   - 0.5*dpx[iwest][j][k]*apu[iwest][j][k] ;
	      				uel   = 1.5*u[i][j][k][l] - 0.5*u[iwest][j][k][l];
					apue  = 1.5*apu[i][j][k]  - 0.5*apu[iwest][j][k];
	      			}
	      			
	      			else
	      			{
	      				dpxel = fxp * dpx[i][j][k] *apu[i][j][k]  + fxe * dpx[ieast][j][k] * apu[ieast][j][k];
	      				uel   = fxp * u[i][j][k][l] + fxe * u[ieast][j][k][l];
	      				apue  = fxp * apu[i][j][k]  + fxe * apu[ieast][j][k];
	      			}
	      			
	      			//Body force contribution
	      			axBFe = fxp * xBF[i][j][k] * apu[i][j][k]  + fxe * xBF[ieast][j][k] * apu[ieast][j][k];
	      			
	      			//xBFe[i][j][k] = fxp * xBF[i][j][k] + fxe * xBF[ieast][j][k];   
	      			xGu =  (apue*xBFe[i][j][k][l] - axBFe);
	      
	      			/* Evaluating cell face gradient and velocity */
	      			dpxe  = (  p[ieast][j][k][l] -  p[i][j][k][l])/dxpe;
	      			
	      			
	      			
	      			//Majumdhar's Correction to make solution independent of Relaxation Coefficients
	      			xMj = (1.0 - uUnderRelaxCoeff)*(ue[i][j][k][l] - fxe*u[ieast][j][k][l] - fxp*u[i][j][k][l]);
	      			
	      			//Choi's correction to make solution independent of time stepSize
	      			xChoi = (ue[i][j][k][l-1]*apue*vole*rhoe - fxe*apu[ieast][j][k]*u[ieast][j][k][l-1]*volE*rho[ieast][j][k] - fxp*apu[i][j][k]*u[i][j][k][l-1]*vol*rho[i][j][k]) * uUnderRelaxCoeff / dt; 
	      			
	      			
	      			 
	      			
	      			ue[i][j][k][l]	= uel - vole*(apue*dpxe - dpxel);
	      			ue[i][j][k][l]	= ue[i][j][k][l] + xGu;  //Adding Gu's correction
	      			ue[i][j][k][l]	= ue[i][j][k][l] + xMj;  //Adding Majumdhar correction
	      			ue[i][j][k][l]	= ue[i][j][k][l] + xChoi;  //Adding Choi's correction
	      			
	      			//printf("%e \n",xChoi);	
	      			
	      			
	      			
	      			Fe[i][j][k] = d*ue[i][j][k][l];	/* Correcting the fluxes with interpolated velocities */
	      
	      			/* Constructing the coefficients */
/*	      			if(i==nxm)*/
/*	      			{*/
/*	      				if(RightBoundType==4)*/
/*	      				{*/
	      					ae[i][j][k]  = -d*apue*s;
	      					aw[ieast][j][k]= ae[i][j][k];
/*	      				}*/
/*	      			}*/
/*	      			if(i==1)*/
/*	      			{*/
/*	      				if(LeftBoundType==4)*/
/*	      				{*/
/*	      					ae[i][j][k]  = -d*apue*s;*/
/*	      					aw[ieast][j][k]= ae[i][j][k];*/
/*	      				}*/
/*	      			}*/
	      				
/*	      			ap[i][j][k] = 0;*/
              		}
	      	  }
	     }
}

void pressure_correction_north()
{
	double d;	/* rho*dx */
  	double dpynl, vnl, apvn;	/* For the interpolation work */
  	double dpyn;	/* Rhie chow interpolated cell face velocities and gradients */ 
  	int jn;
	double rhon;
	
	double ayBFn;
	
	double yGu,yMj,yChoi;
	
	double volN;
  
  
  for (j=1;j<=nym;j++)
  {
	jnorth = j+1;
	jnorthnorth = j+2;
	jsouth = j-1;
	dypn = yc[jnorth] - yc[j];
	fyn = (yf[j] - yc[j])/dypn;	
    	fyp = 1.0 - fyn;
    	
    	if(j==nym&&TopBoundType!=4)	continue;
	if(j==1&&BottomBoundType!=4)	continue;

  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  		
	  		s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
			voln = s * dypn;
			
			vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
			volN = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[jnorth] - yf[jnorth-1]);

			/* Arithmetic Interpolation of density*/
			rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];
			d = rhon * s;
		
      			/* Since, volume of the east cell equal to that of the central face, */
      			/* interpolated cel face quantites can be written as follows */
      			if(j==1)
      			{
      				//Extrapolation of Rhie Chow from interior cells  
      				dpynl = 1.5 * dpy[i][jnorth][k]*apv[i][jnorth][k]  - 0.5 * dpy[i][jnorthnorth][k]*apv[i][jnorthnorth][k];
      				vnl   = 1.5 * v[i][jnorth][k][l] - 0.5 * v[i][jnorthnorth][k][l];
      				apvn  = 1.5 * apv[i][jnorth][k]  - 0.5 * apv[i][jnorthnorth][k];
      			}
      			
      			else if(j==nym)
      			{
      				dpynl = 1.5 * dpy[i][j][k]*apv[i][j][k]  - 0.5 * dpy[i][jsouth][k]*apv[i][jsouth][k];
      				vnl   = 1.5 * v[i][j][k][l] - 0.5 * v[i][jsouth][k][l];
      				apvn  = 1.5 * apv[i][j][k]  - 0.5 * apv[i][jsouth][k];
      			}
      			
      			else
      			{
      				dpynl = fyp * dpy[i][j][k] *apv[i][j][k]  + fyn * dpy[i][jnorth][k]*apv[i][jnorth][k];
      				vnl   = fyp * v[i][j][k][l] + fyn * v[i][jnorth][k][l];
      				apvn  = fyp * apv[i][j][k]  + fyn * apv[i][jnorth][k];
      			}
      			
      			//Body force contribution
	      		ayBFn = fyp * yBF[i][j][k] * apv[i][j][k]  + fyn * yBF[i][jnorth][k] * apv[i][jnorth][k];
	      		//yBFn[i][j][k] = fyp * yBF[i][j][k] + fyn * yBF[i][jnorth][k];
	      		yGu = (apvn*yBFn[i][j][k][l] - ayBFn);  
      	
      	
      			/* Evaluating cell face gradient and velocity */
      			dpyn  = ( p[i][jnorth][k][l] -  p[i][j][k][l])/dypn;
      			
      			//Majumdhar's Correction to make solution independent of Relaxation Coefficients
	      		yMj = (1.0 - vUnderRelaxCoeff)*(vn[i][j][k][l] - fyn*v[i][jnorth][k][l] - fyp*v[i][j][k][l]);
	      		
	      		//Choi's correction to make solution independent of time stepSize
	      		yChoi = (vn[i][j][k][l-1]*apvn*voln*rhon - fyn*apv[i][jnorth][k]*v[i][jnorth][k][l-1]*volN*rho[i][jnorth][k] - fyp*apv[i][j][k]*v[i][j][k][l-1]*vol*rho[i][j][k]) * vUnderRelaxCoeff / dt;
      			
      			
      			
      			vn[i][j][k][l]	= vnl - voln*(apvn*dpyn - dpynl);
      			vn[i][j][k][l]	= vn[i][j][k][l] + yGu;//Adding Gu's correction
      			vn[i][j][k][l]	= vn[i][j][k][l] + yMj;//Adding Majumdhar correction
      			vn[i][j][k][l]	= vn[i][j][k][l] + yChoi;	//Adding Choi's correction
      			
      			Fn[i][j][k] = d*vn[i][j][k][l];	/* Correcting the fluxes with interpolated velocities */
      	
      			/* Constructing the coefficients */
/*      			if(j==nym)*/
/*	      		{*/
/*	      			if(TopBoundType==4)*/
/*	      			{*/
/*	      				an[i][j][k]  = -d*apvn*s;*/
/*      					as[i][jnorth][k] = an[i][j][k];*/
/*	      			}*/
/*	      		}*/
/*	      		*/
/*	      		if(j==1)*/
/*	      		{*/
/*	      			if(BottomBoundType==4)*/
/*	      			{*/
	      				an[i][j][k]  = -d*apvn*s;
      					as[i][jnorth][k] = an[i][j][k];
/*	      			}*/
/*	      		}*/
      	
      			
     		}
    	}
  }
}

void pressure_correction_top()
{
	double d;	/* rho*dx */
  	double dpztl, wtl, apwt;	/* For the interpolation work */
  	double dpzt;	/* Rhie chow interpolated cell face velocities and gradients */ 
  	int kt;
	double rhot;
	
	double azBFt;
	
	double zGu,zMj,zChoi;
	
	double volT;
  
  
  for (k=1;k<=nzm;k++)
  {
	ktop = k+1;
	ktoptop = k+2;
	kbottom = k-1;
	dzpt = zc[ktop]-zc[k];
	fzt = (zf[k] - zc[k])/dzpt;
	fzp = 1.0 - fzt; 
	
	if(k==nzm&&FrontBoundType!=4)	continue;
	if(k==1&&BackBoundType!=4)	continue;
	
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  		
	  		s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
			volt = s * dzpt;
			
			vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
			volT = (xf[i]-xf[i-1])*(zf[ktop]-zf[ktop-1])*(yf[j] - yf[j-1]);

			/* Arithmetic Interpolation of density*/
			rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];
			d = rhot * s;
		
	      		/* Since, volume of the east cell equal to that of the central face, */
	      		/* interpolated cel face quantites can be written as follows */
	      		if(k==1)
	      		{
	      			//Extrapolation of Rhie Chow from interior cells  
	      			dpztl = 1.5 * dpz[i][j][ktop] * apw[i][j][ktop]  - 0.5 * dpz[i][j][ktoptop]* apw[i][j][ktoptop] ;
	      			wtl   = 1.5 * w[i][j][ktop][l] - 0.5 * w[i][j][ktoptop][l];
	      			apwt  = 1.5 * apw[i][j][ktop]  - 0.5 * apw[i][j][ktoptop];
	      		}

			else if(k==nzm)
			{
				dpztl = 1.5 * dpz[i][j][k]* apw[i][j][k]   - 0.5 * dpz[i][j][kbottom]* apw[i][j][kbottom] ;
	      			wtl   = 1.5 * w[i][j][k][l] - 0.5 * w[i][j][kbottom][l];
	      			apwt  = 1.5 * apw[i][j][k]  - 0.5 * apw[i][j][kbottom];
			}
			
			else
			{
				dpztl = fzp * dpz[i][j][k] * apw[i][j][k]  + fzt * dpz[i][j][ktop] * apw[i][j][ktop];
	      			wtl   = fzp * w[i][j][k][l] + fzt * w[i][j][ktop][l];
	      			apwt  = fzp * apw[i][j][k]  + fzt * apw[i][j][ktop];
			}
	      		
	      		//Body force contribution
	      		azBFt = fzp * zBF[i][j][k] * apw[i][j][k] + fzt * zBF[i][j][ktop] * apw[i][j][ktop];
	      		//zBFt[i][j][k] = fzp * zBF[i][j][k]  + fzt * zBF[i][j][ktop]; 
	      		zGu = (apwt*zBFt[i][j][k][l] - azBFt);  
      
	      		/* Evaluating cell face gradient and velocity */
	      		dpzt  = (p[i][j][ktop][l] -  p[i][j][k][l])/dzpt;
	      		
	      		//Majumdhar's Correction to make solution independent of Relaxation Coefficients
	      		zMj = (1 - wUnderRelaxCoeff)*(wt[i][j][k][l] - fzt*w[i][j][ktop][l] - fzp*w[i][j][k][l]);
	      		
	      		//Choi's correction to make solution independent of time stepSize
	      		zChoi = (wt[i][j][k][l-1]*apwt*volt*rhot - fzt*apw[i][j][ktop]*w[i][j][ktop][l-1]*volT*rho[i][j][ktop] - fzp*apw[i][j][k]*w[i][j][k][l-1]*vol*rho[i][j][k]) * wUnderRelaxCoeff / dt;
	      		
	      		
	      		
	      		wt[i][j][k][l]	= wtl - volt*(apwt*dpzt - dpztl);
	      		wt[i][j][k][l]	= wt[i][j][k][l] + zGu;//Adding Gu's correction
	      		wt[i][j][k][l] 	= wt[i][j][k][l] + zMj;//Adding Majumdhar correction
	      		wt[i][j][k][l] 	= wt[i][j][k][l] + zChoi;	//Adding Choi correction
	      		
	      		Ft[i][j][k] = d*wt[i][j][k][l];	/* Correcting the fluxes with interpolated velocities */
      
	      		/* Constructing the coefficients */
/*	      		if(k==nzm)*/
/*	      		{*/
/*	      			if(FrontBoundType==4)*/
/*	      			{*/
/*	      				at[i][j][k]  = -d*apwt*s;*/
/*	      				ab[i][j][ktop] = at[i][j][k];*/
/*	      			}*/
/*	      		}*/
/*	      		*/
/*	      		if(k==1)*/
/*	      		{*/
/*	      			if(BackBoundType==4)*/
/*	      			{*/
	      				at[i][j][k]  = -d*apwt*s;
	      				ab[i][j][ktop] = at[i][j][k];
/*	      			}*/
/*	      		}*/
	      		
     		}
    	}
  }
}

void pressure_correction_source()
{
	maxMassResidual=0; /* To check for continuity */
  
  	for (i=2;i<=nxm;i++)
	  {
	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	 {
		  	 	ieast 	= i + 1;
      				iwest 	= i - 1;
      				jnorth	= j + 1;
      				jsouth	= j - 1;
      				ktop  	= k + 1;
      				kbottom	= k - 1;
            
      				su[i][j][k] = Fe[iwest][j][k] - Fe[i][j][k] + Fn[i][jsouth][k] - Fn[i][j][k] + Ft[i][j][kbottom] - Ft[i][j][k] ;
      				ap[i][j][k] = - (ae[i][j][k] + aw[i][j][k] + an[i][j][k] + as[i][j][k] + at[i][j][k] + ab[i][j][k]);
      				maxMassResidual	= maxMassResidual + fabs(su[i][j][k]);	/* Checking continuity */
      				pp[i][j][k][l] = 0;
    		 }
    	}
  	 }
}

