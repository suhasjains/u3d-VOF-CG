#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void uEqnCoeff()
{
	
	
	double rhoe,rhow,rhoc;
	
	double Te,Tw;
	
	double Betae,Betaw;
	
	double SB;
	
	
	//Constructing convective fluxes along east face 
	for  (i=wcID[e][f][g]-1; i<=ecID[e][f][g]; i++)
 	{
		ieast = i+1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;
		fxp = 1.0 - fxe;
		
		if(i>=nxm)	continue;
		if(i<2)		continue;

	 	for (j=scID[e][f][g]; j<=ncID[e][f][g]; j++)
   		{
	  		for(k=bcID[e][f][g]; k<=tcID[e][f][g]; k++)
			{
		
		
				/* Evaluation of cell face area */
				s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
		
				/* Evaluation of diffusive term */
				/*Interpolating viscosity field*/	
				mue = fxp*mu[i][j][k] + fxe*mu[ieast][j][k];				
				diff = mue*s/dxpe;
		
				ieasteast 	= i+2;
				iwest  		= i-1;
		
				upwind_implicit(Fe);

      				if(scheme==1)
				{
					/* Constructing co-efficients for upwind*/
/*			      		if(i!=(wcID[e][f][g]-1))	ae[i][j][k]     =  convf - diff;*/
/*			      		if(i!=ecID[e][f][g])		aw[ieast][j][k] = -convp - diff;*/
			      		
			      		ae[i][j][k]     =  convf - diff;
			      		aw[ieast][j][k] = -convp - diff;
				}
				else if(scheme==2)
				{
					/* Constructing co-efficients for hybrid*/
			      		ae[i][j][k]     =  convf - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0);
			      		aw[ieast][j][k] = -convp - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0);
      				}
				else if(scheme==3)
				{
					/* Constructing co-efficients for power law*/
			      		ae[i][j][k]     =  convf - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0);
      					aw[ieast][j][k] = -convp - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0);
				}

/*//				Deferred Correction Upwind Flux*/
/*				if(i==ecID[e][f][g])*/
/*				{*/
/*					fuuds = convp*u[i][j][k][l] + convf*ueBA[1][j][k][b];*/
/*      					fvuds = convp*v[i][j][k][l] + convf*veBA[1][j][k][b];*/
/*      					fwuds = convp*w[i][j][k][l] + convf*weBA[1][j][k][b];*/
/*      					*/
/*      					fuuds = convp*u[i][j][k][l] + convf*u[ieast][j][k][l];*/
/*      					fvuds = convp*v[i][j][k][l] + convf*v[ieast][j][k][l];*/
/*      					fwuds = convp*w[i][j][k][l] + convf*w[ieast][j][k][l];*/
/*				}*/
/*				else*/
/*				{*/
					fuuds = convp*u[i][j][k][l] + convf*u[ieast][j][k][l];
      					fvuds = convp*v[i][j][k][l] + convf*v[ieast][j][k][l];
      					fwuds = convp*w[i][j][k][l] + convf*w[ieast][j][k][l];
/*      				}*/
      				
      				
//				Deferred Correction Central Differencing Flux
				if(i==ecID[e][f][g])
				{
					fucds = Fe[i][j][k]*(ueBA[1][j][k][b]*fxe + u[i][j][k][l]*fxp);
      					fvcds = Fe[i][j][k]*(veBA[1][j][k][b]*fxe + v[i][j][k][l]*fxp);
      					fwcds = Fe[i][j][k]*(weBA[1][j][k][b]*fxe + w[i][j][k][l]*fxp);
      				}
      				else
      				{
      					fucds = Fe[i][j][k]*(u[ieast][j][k][l]*fxe + u[i][j][k][l]*fxp);
      					fvcds = Fe[i][j][k]*(v[ieast][j][k][l]*fxe + v[i][j][k][l]*fxp);
      					fwcds = Fe[i][j][k]*(w[ieast][j][k][l]*fxe + w[i][j][k][l]*fxp);
      				}
      				
      				
      				
      				//			Deferred Correction Hybrid Flux
/*      				if(i==ecID[e][f][g])*/
/*      				{*/
/*      					fuH = convp*u[i][j][k][l] + convf*ueBA[1][j][k][b] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(ueBA[1][j][k][b] - u[i][j][k][l]);*/
/*      					fvH = convp*v[i][j][k][l] + convf*veBA[1][j][k][b] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(veBA[1][j][k][b] - v[i][j][k][l]);*/
/*      					fwH = convp*w[i][j][k][l] + convf*weBA[1][j][k][b] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(weBA[1][j][k][b] - w[i][j][k][l]);*/
/*      				}*/
/*      				else*/
/*      				{*/
      					fuH = convp*u[i][j][k][l] + convf*u[ieast][j][k][l] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(u[ieast][j][k][l] - u[i][j][k][l]);
      					fvH = convp*v[i][j][k][l] + convf*v[ieast][j][k][l]	- diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(v[ieast][j][k][l] - v[i][j][k][l]);
      					fwH = convp*w[i][j][k][l] + convf*w[ieast][j][k][l] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(w[ieast][j][k][l] - w[i][j][k][l]);
/*      				}*/
      				
      				
/*		//		Deferred Correction Power Law Flux  */
/*				fuP = convp*u[i][j][k][l] + convf*u[ieast][j][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0)*(u[ieast][j][k][l] - u[i][j][k][l]);*/
/*      				fvP = convp*v[i][j][k][l] + convf*v[ieast][j][k][l]	- diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0)*(v[ieast][j][k][l] - v[i][j][k][l]);*/
/*      				fwP = convp*w[i][j][k][l] + convf*w[ieast][j][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0)*(w[ieast][j][k][l] - w[i][j][k][l]);  */
/*      	*/
      	//			Deferred Correction QUICK Flux
				fuQ = DC_QUICK_flux_east(fuQ,u);
				fvQ = DC_QUICK_flux_east(fvQ,v);
				fwQ = DC_QUICK_flux_east(fwQ,w);

		

		//		Deferred Correction QUICK Flux when hybrid and Power Law is used to predict
				fuQHP = DC_QUICK_flux_east(fuQHP,u);
				fuQHP = fuQHP - diff*(u[ieast][j][k][l] - u[i][j][k][l]);			
		
				fvQHP = DC_QUICK_flux_east(fvQHP,v);
				fvQHP = fvQHP - diff*(v[ieast][j][k][l] - v[i][j][k][l]);
				
				fwQHP = DC_QUICK_flux_east(fwQHP,w);
				fwQHP = fwQHP - diff*(w[ieast][j][k][l] - w[i][j][k][l]);     
      	
      
      				//Adding contribution to source term
      				if(scheme==1)	v_east_source_dc(fuuds,fvuds,fwuds,fuQ,fvQ,fwQ);
      				if(scheme==2)	v_east_source_dc(fuH,fvH,fwH,fuQHP,fvQHP,fwQHP);
				if(scheme==3)	v_east_source_dc(fuP,fvP,fwP,fuQHP,fvQHP,fwQHP);
			}
  		}
	}
	
	

	
	
	//Constructing Source Terms
	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		dx = xf[i] - xf[i-1];
		ieast 	= i + 1;
      		iwest 	= i - 1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;
		dxpw = xc[i] - xc[iwest];
		fxw = (xf[iwest]-xc[iwest])/dxpw;


		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			dy = yf[j] - yf[j-1];
			jnorth	= j + 1;
      			jsouth	= j - 1;
			dypn = yc[jnorth] - yc[j];
			fyn = (yf[j] - yc[j])/dypn;	
    			fyp = 1.0 - fyn;
			dyps = yc[j] - yc[jsouth];
			fys = (yf[jsouth] - yc[jsouth])/dyps;

			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				dz = zf[k] - zf[k-1];
				ktop  	= k + 1;
		  		kbottom	= k - 1;
				dzpt = zc[ktop]-zc[k];
				fzt = (zf[k] - zc[k])/dzpt;
				fzp = 1.0 - fzt;
				dzpb = zc[k]-zc[kbottom];
				fzb = (zf[kbottom] - zc[kbottom])/dzpb;
				
				vol	= dx*dy*dz;
				
				
				
/*				if(rhoc!=rho[i][j][k])	printf("rho = %e %e \n",rhoc,rho[i][j][k]);*/
				  				
				
				/* Constructing the pressure gradient terms */
      				pe = fxe	*p[ieast][j][k][l]  	+ fxp*p[i][j][k][l];
      				pw = (1.0-fxw)	*p[iwest][j][k][l]	+ fxw*p[i][j][k][l];
      				pn = fyn	*p[i][jnorth][k][l] 	+ fyp*p[i][j][k][l];
      				ps = (1.0-fys)	*p[i][jsouth][k][l] 	+ fys*p[i][j][k][l];
      				pt = fzt	*p[i][j][ktop][l]	+ fzp*p[i][j][k][l];
      				pb = (1.0-fzb)	*p[i][j][kbottom][l]	+ fzb*p[i][j][k][l];
      
      				dpx[i][j][k] = (pe - pw)/dx;
      				dpy[i][j][k] = (pn - ps)/dy;
      				dpz[i][j][k] = (pt - pb)/dz;
      			
      			
      				/* Updating source terms with pressure gradient */
      				su[i][j][k] = su[i][j][k] - dpx[i][j][k]*vol;//*rho[i][j][k]/rhoc;
      				
			}
		}
	}
	
	
//	Adding unsteady contribution to source terms
	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				apt	= rho[i][j][k]*vol/dt;
				if(simulationTime==dt)
				{
					su[i][j][k]	= su[i][j][k] + apt*u[i][j][k][l-1];
					apu[i][j][k]	= apu[i][j][k] + apt;  
				}
				else
				{
					su[i][j][k]	= su[i][j][k] + (1.0 + Gamma)*apt*u[i][j][k][l-1] - 0.5*Gamma*apt*u[i][j][k][l-2];
					apu[i][j][k]	= apu[i][j][k] + (1.0 + 0.5*Gamma)*apt;  
				}			
			}
		}
	}
	
/*	//Adding Buoyancy Term*/
/*	for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)*/
/*	{*/

/*		for(j=scID[e][f][g];j<=ncID[e][f][g];j++)*/
/*		{*/
/*			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)*/
/*			{*/
/*				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);*/
/*				*/
/*				SB = Beta * rho[i][j][k] * vol * (T[i][j][k][l] - Tref);*/
/*				*/
/*				su[i][j][k] = su[i][j][k] + xTempGrav*SB;*/
/*				*/
/*			}*/
/*		}*/
/*	} */
/*	*/
	

	//Adding gravity contribution to source terms
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
			for(k=2;k<=nzm;k++)
			{
			
				/* Arithmetic Interpolation of density*/
				rhoe = fxp * rho[i][j][k] +   fxe * rho[ieast][j][k];
				rhow = (1.0-fxw) * rho[i][j][k] +   fxw * rho[iwest][j][k];
				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				/* Arithmetic Interpolation of Temperature*/
				Te = fxp * T[i][j][k][l] +   fxe * T[ieast][j][k][l];
				Tw = (1.0-fxw) * T[i][j][k][l] +   fxw * T[iwest][j][k][l];
				
				/* Arithmetic Interpolation of Temperature*/
				Betae = fxp * Beta[i][j][k] +   fxe * Beta[ieast][j][k];
				Betaw = (1.0-fxw) * Beta[i][j][k] +   fxw * Beta[iwest][j][k];
				
				
				xBFe[i][j][k][l] = rhoe*vol*xGrav  - rhoe*vol*xTempGrav*Betae*(Te-Tref);				//Volume effect not yet included
				xBFw[i][j][k] = rhow*vol*xGrav - rhow*vol*xTempGrav*Betaw*(Tw-Tref);
				
				//Mencinger's discontinuous Body force interpolation
				xBF[i][j][k] = fxe*xBFe[i][j][k][l] + fxw*xBFw[i][j][k];
				
				//Adding Surface Tension Body force
				xBFe[i][j][k][l] = xBFe[i][j][k][l] + xSTe[i][j][k]* vol;  		//volume effect not yet included 
				xBF[i][j][k] = xBF[i][j][k] + xST[i][j][k]*vol;
				
				//if(xST[i][j][k]!=0)	printf("%e \n",xST[i][j][k]);
				
				su[i][j][k]	= su[i][j][k] + xBF[i][j][k];
			}
		}
	}


}

