#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void wEqnCoeff()
{
	


	double rhot,rhob,rhoc;
	
	double Tt,Tb;
	
	double Betat,Betab;
	
	double SB;

	//Constructing convective fluxes along top face
	for(k=bcID[e][f][g]-1;k<=tcID[e][f][g];k++)
	{
		ktop = k+1;
		dzpt = zc[ktop]-zc[k];
		fzt = (zf[k] - zc[k])/dzpt;
		fzp = 1.0 - fzt;
		
		if(k>=nzm)	continue; 
		if(k<2)		continue;
		
		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
      
	      			ktoptop	= k+2;
	      			kbottom	= k-1;
      
				/*Interpolating viscosity field*/				
				mut = fzp*mu[i][j][k] + fzt*mu[i][j][ktop];
			      	diff = mut*s/dzpt;
      
      				upwind_implicit(Ft);
      			
				if(scheme==1)
				{
					/* Constructing co-efficients for upwind */
/*      					if(k!=(bcID[e][f][g]-1))	at[i][j][k]   	=  convf - diff;*/
/*      					if(k!=tcID[e][f][g])		ab[i][j][ktop]	= -convp - diff;*/
      					
      					at[i][j][k]   	=  convf - diff;
      					ab[i][j][ktop]	= -convp - diff;
				}
				else if(scheme==2)
				{
					/* Constructing co-efficients for hybrid */
      				at[i][j][k]	   	=  convf - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0);
      				ab[i][j][ktop]	= -convp - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0);
      				}
				else if(scheme==3)
				{
					/* Constructing co-efficients for power law */
      					at[i][j][k]	   	=  convf - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0);
      					ab[i][j][ktop]	= -convp - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0);
				}

				
				
				//Deferred Correction Upwind Flux
/*				if(k==tcID[e][f][g])*/
/*				{*/
/*					fuuds = convp*u[i][j][k][l] + convf*utBA[i][j][1][b];*/
/*      					fvuds = convp*v[i][j][k][l] + convf*vtBA[i][j][1][b];*/
/*      					fwuds = convp*w[i][j][k][l] + convf*wtBA[i][j][1][b];*/
/*      					*/
/*      					fuuds = convp*u[i][j][k][l] + convf*u[i][j][ktop][l];*/
/*      					fvuds = convp*v[i][j][k][l] + convf*v[i][j][ktop][l];*/
/*      					fwuds = convp*w[i][j][k][l] + convf*w[i][j][ktop][l];*/
/*      				}*/
/*      				else*/
/*      				{*/
      					fuuds = convp*u[i][j][k][l] + convf*u[i][j][ktop][l];
      					fvuds = convp*v[i][j][k][l] + convf*v[i][j][ktop][l];
      					fwuds = convp*w[i][j][k][l] + convf*w[i][j][ktop][l];
/*      				}*/
      				
      				
      			//		Deferred Correction Central Differencing Flux
				fucds = Ft[i][j][k]*(u[i][j][ktop][l]*fzt + u[i][j][k][l]*fzp);
      				fvcds = Ft[i][j][k]*(v[i][j][ktop][l]*fzt + v[i][j][k][l]*fzp);
      				fwcds = Ft[i][j][k]*(w[i][j][ktop][l]*fzt + w[i][j][k][l]*fzp);
      			
      			//		Deferred Correction Hybrid Flux
      				fuH = convp*u[i][j][k][l] + convf*u[i][j][ktop][l] - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0)*(u[i][j][ktop][l] - u[i][j][k][l]);
      				fvH = convp*v[i][j][k][l] + convf*v[i][j][ktop][l] - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0)*(v[i][j][ktop][l] - v[i][j][k][l]);
      				fwH = convp*w[i][j][k][l] + convf*w[i][j][ktop][l] - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0)*(w[i][j][ktop][l] - w[i][j][k][l]);
      			
/*      			//		Deferred Correction Power Law Flux*/
/*      				fuP = convp*u[i][j][k][l] + convf*u[i][j][ktop][l] - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0)*(u[i][j][ktop][l] - u[i][j][k][l]);*/
/*      				fvP = convp*v[i][j][k][l] + convf*v[i][j][ktop][l] - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0)*(v[i][j][ktop][l] - v[i][j][k][l]);*/
/*      				fwP = convp*w[i][j][k][l] + convf*w[i][j][ktop][l] - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0)*(w[i][j][ktop][l] - w[i][j][k][l]);*/
/*				*/
				
				
				
				//		Deferred Correction QUICK Flux
				fuQ = DC_QUICK_flux_top(fuQ,u);
				fvQ = DC_QUICK_flux_top(fvQ,v);
				fwQ = DC_QUICK_flux_top(fwQ,w);
				
				//		Deferred Correction QUICK Flux when hybrid and Power Law is used to predict
				fuQHP = DC_QUICK_flux_top(fuQHP,u);	
				fuQHP = fuQHP - diff*(u[i][j][ktop][l] - u[i][j][k][l]);
				
				fvQHP = DC_QUICK_flux_top(fvQHP,v);					
				fvQHP = fvQHP - diff*(v[i][j][ktop][l] - v[i][j][k][l]);
				  
				fwQHP = DC_QUICK_flux_top(fwQHP,w);					
				fwQHP = fwQHP - diff*(w[i][j][ktop][l] - w[i][j][k][l]);    
      			
      			
       
      				if(scheme==1)	v_top_source_dc(fuuds,fvuds,fwuds,fuQ,fvQ,fwQ);
      				if(scheme==2)	v_top_source_dc(fuH,fvH,fwH,fuQHP,fvQHP,fwQHP);
				if(scheme==3)	v_top_source_dc(fuP,fvP,fwP,fuQHP,fvQHP,fwQHP);
      							
				
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
				sw[i][j][k] = sw[i][j][k] - dpz[i][j][k]*vol;//*rho[i][j][k]/rhoc;
								
				
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
					sw[i][j][k]	= sw[i][j][k] + apt*w[i][j][k][l-1];
					apw[i][j][k]	= apw[i][j][k] + apt;
				}
				else
				{
					sw[i][j][k]	= sw[i][j][k] + (1.0 + Gamma)*apt*w[i][j][k][l-1] - 0.5*Gamma*apt*w[i][j][k][l-2];
					apw[i][j][k]	= apw[i][j][k] + (1.0 + 0.5*Gamma)*apt;	
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
/*				su[i][j][k] = su[i][j][k] + zTempGrav*SB;*/
/*				*/
/*			}*/
/*		}*/
/*	}*/

	//Adding gravity contribution to source terms
	for(i=2;i<=nxm;i++)
	{
		for(j=2;j<=nym;j++)
		{
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
			
			
				/* Arithmetic Interpolation of density*/
				rhot = fzp * rho[i][j][k] +   fzt * rho[i][j][ktop];
				rhob = (1.0-fzb) * rho[i][j][k] +   fzb * rho[i][j][kbottom];
				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				/* Arithmetic Interpolation of Temperature*/
				Tt = fzp * T[i][j][k][l] +   fzt * T[i][j][ktop][l];
				Tb = (1.0-fzb) * T[i][j][k][l] +   fzb * T[i][j][kbottom][l];
				
				/* Arithmetic Interpolation of Temperature*/
				Betat = fzp * Beta[i][j][k] +   fzt * Beta[i][j][ktop];
				Betab = (1.0-fzb) * Beta[i][j][k] +   fzb * Beta[i][j][kbottom];
				
				zBFt[i][j][k][l] = rhot*vol*zGrav - rhot*vol*zTempGrav*Betat*(Tt-Tref);				//Volume effect not yet included
				zBFb[i][j][k] = rhob*vol*zGrav - rhob*vol*zTempGrav*Betab*(Tb-Tref);
			
				//Mencinger's discontinuous Body force interpolation
				zBF[i][j][k] = fzt*zBFt[i][j][k][l] + fzb*zBFb[i][j][k];
				
				//Adding Surface Tension Body force
				zBFt[i][j][k][l] = zBFt[i][j][k][l] + zSTt[i][j][k]*vol;   
				zBF[i][j][k] = zBF[i][j][k] + zST[i][j][k]*vol;
				
				//if(zST[i][j][k]!=0)	printf("%e \n",zST[i][j][k]);
				
				sw[i][j][k]	= sw[i][j][k] + zBF[i][j][k];
			}
		}
	}


}

