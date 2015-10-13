#include"update.h"

//Updates the values of the pressure
void update_Pressure_Velocity()
{
	double ppe,ppw,ppn,pps,ppt,ppb; /* Linear interpolation of pprime */
	
	
  
  for(i=2;i<=nxm;i++)
  {
	dx = xf[i] - xf[i-1];
	ieast 	= i + 1;
      	iwest 	= i - 1;
	dxpe = xc[ieast] - xc[i];
	fxe = (xf[i]-xc[i])/dxpe;		
	fxp = 1.0 - fxe;
	dxpw = xc[i] - xc[iwest];
	fxw = (xf[iwest]-xc[iwest])/dxpw;

    	for(j=2;j<=nym;j++)
	  {
		dy = yf[j] - yf[j-1];
		jnorth	= j + 1;
      		jsouth	= j - 1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;	
    		fyp = 1.0 - fyn;
		dyps = yc[j] - yc[jsouth];
		fys = (yf[jsouth] - yc[jsouth])/dyps;
	
	  	for(k=2;k<=nzm;k++)
		  {
			dz = zf[k] - zf[k-1];
			ktop  	= k + 1;
		  	kbottom	= k - 1;
			dzpt = zc[ktop]-zc[k];
			fzt = (zf[k] - zc[k])/dzpt;
			fzp = 1.0 - fzt;
			dzpb = zc[k]-zc[kbottom];
			fzb = (zf[kbottom] - zc[kbottom])/dzpb;
			
      		
			p[i][j][k][l] = p[i][j][k][l] + pUnderRelaxCoeff*(pp[i][j][k][l] - ppo);
      
		  } 
      	}
  }
  
  
  for(i=2;i<=nxm;i++)
  {
	dx = xf[i] - xf[i-1];
	ieast 	= i + 1;
      	iwest 	= i - 1;
	dxpe = xc[ieast] - xc[i];
	fxe = (xf[i]-xc[i])/dxpe;		
	fxp = 1.0 - fxe;
	dxpw = xc[i] - xc[iwest];
	fxw = (xf[iwest]-xc[iwest])/dxpw;

    	for(j=2;j<=nym;j++)
	  {
		dy = yf[j] - yf[j-1];
		jnorth	= j + 1;
      		jsouth	= j - 1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;	
    		fyp = 1.0 - fyn;
		dyps = yc[j] - yc[jsouth];
		fys = (yf[jsouth] - yc[jsouth])/dyps;
	
	  	for(k=2;k<=nzm;k++)
		  {
			dz = zf[k] - zf[k-1];
			ktop  	= k + 1;
		  	kbottom	= k - 1;
			dzpt = zc[ktop]-zc[k];
			fzt = (zf[k] - zc[k])/dzpt;
			fzp = 1.0 - fzt;
			dzpb = zc[k]-zc[kbottom];
			fzb = (zf[kbottom] - zc[kbottom])/dzpb;
			
			
			vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
			
		  	ppe = fxe*pp[ieast][j][k][l]  		+ fxp*pp[i][j][k][l];	/* Interpolating pprime */
      			ppw = (1.0 - fxw)*pp[iwest][j][k][l]	+ fxw*pp[i][j][k][l];
      			ppn = fyn*pp[i][jnorth][k][l] 		+ fyp*pp[i][j][k][l];
      			pps = (1.0 - fys)*pp[i][jsouth][k][l] 	+ fys*pp[i][j][k][l];
      			ppt = fzt*pp[i][j][ktop][l] 		+ fzp*pp[i][j][k][l];
      			ppb = (1.0 - fzb)*pp[i][j][kbottom][l]	+ fzb*pp[i][j][k][l];	
      			
			
			if(i==nxm)	u[i][j][k][l] = u[i][j][k][l] - (ppe - pp[i][j][k][l])*dy*dz*apu[i][j][k];
			else if(i==2)	u[i][j][k][l] = u[i][j][k][l] - (pp[i][j][k][l] - ppw)*dy*dz*apu[i][j][k];
			else		u[i][j][k][l] = u[i][j][k][l] - (ppe - ppw)*dy*dz*apu[i][j][k];

			
			if(j==nym)	v[i][j][k][l] = v[i][j][k][l] - (ppn - pp[i][j][k][l])*dx*dz*apv[i][j][k];
			else if (j==2)	v[i][j][k][l] = v[i][j][k][l] - (pp[i][j][k][l] - pps)*dx*dz*apv[i][j][k];
			else		v[i][j][k][l] = v[i][j][k][l] - (ppn - pps)*dx*dz*apv[i][j][k];

			
			if(k==nzm)	w[i][j][k][l] = w[i][j][k][l] - (ppt - pp[i][j][k][l])*dx*dy*apw[i][j][k];
			else if(k==2)	w[i][j][k][l] = w[i][j][k][l] - (pp[i][j][k][l] - ppb)*dx*dy*apw[i][j][k];
			else    	w[i][j][k][l] = w[i][j][k][l] - (ppt - ppb)*dx*dy*apw[i][j][k];
			
			/* Arithmetic Interpolation of density*/
/*				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];*/
/*				rhow = (1.0-fxw) * rho[i][j][k] +   fxw * rho[iwest][j][k];*/

				/* Arithmetic Interpolation of density*/
/*				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];*/
/*				rhos = (1.0-fys) * rho[i][j][k] +   fys * rho[i][jsouth][k];*/
/*				*/
				/* Arithmetic Interpolation of density*/
/*				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];*/
/*				rhob = (1.0-fzb) * rho[i][j][k] +   fzb * rho[i][j][kbottom];*/
				


			u[i][j][k][l] = u[i][j][k][l] + (xBFe[i][j][k][l] + xBFe[iwest][j][k][l])*vol*apu[i][j][k]/2.0;			
      			u[i][j][k][l] = u[i][j][k][l] - (xBFe[i][j][k][l-1] + xBFe[iwest][j][k][l-1])*vol*apu[i][j][k]/2.0;
      			
      			v[i][j][k][l] = v[i][j][k][l] + (yBFn[i][j][k][l] + yBFn[i][jsouth][k][l])*vol*apv[i][j][k]/2.0;			
      			v[i][j][k][l] = v[i][j][k][l] - (yBFn[i][j][k][l-1] + yBFn[i][jsouth][k][l-1])*vol*apv[i][j][k]/2.0;
      			
      			w[i][j][k][l] = w[i][j][k][l] + (zBFt[i][j][k][l] + zBFt[i][j][kbottom][l])*vol*apw[i][j][k]/2.0;			
      			w[i][j][k][l] = w[i][j][k][l] - (zBFt[i][j][k][l-1] + zBFt[i][j][kbottom][l-1])*vol*apw[i][j][k]/2.0;
      
		  } 
      	}
  }
  
  
}


