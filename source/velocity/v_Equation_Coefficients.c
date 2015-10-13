#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void vEqnCoeff()
{
	
	
	double rhon,rhos,rhoc;
	
	double Tn,Ts;
	
	double Betan,Betas;

	double SB;

	//Constructing convective fluxes along north face
	for(j=scID[e][f][g]-1;j<=ncID[e][f][g];j++)
	{
		jnorth = j+1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;
	    	fyp = 1.0 - fyn;
	    	
	    	if(j>=nym)	continue;
	    	if(j<2)		continue;

		for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				
				s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
      
	      			jnorthnorth	= j+2;
	      			jsouth		= j-1;
      				
				/*Interpolating viscosity field*/
				mun = fyp*mu[i][j][k] + fyn*mu[i][jnorth][k];
		      		diff = mun*s/dypn;
      	
      				upwind_implicit(Fn);
      			
				if(scheme==1)
				{
					/* Constructing co-efficients for upwind */
/*      					if(j!=(scID[e][f][g]-1))	an[i][j][k]	   	=  convf - diff;*/
/*      					if(j!=ncID[e][f][g])		as[i][jnorth][k]= -convp - diff;*/
      					
      					an[i][j][k]	   	=  convf - diff;
      					as[i][jnorth][k]= -convp - diff;
				}
				else if(scheme==2)
				{
					/* Constructing co-efficients for hybrid */
      					an[i][j][k]	   	=  convf - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0);
      					as[i][jnorth][k]= -convp - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0);
      				}
				else if(scheme==3)
				{
					/* Constructing co-efficients for power law */
      					an[i][j][k]	   	=  convf - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0);
      					as[i][jnorth][k]= -convp - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0);
				}
				
				//		Deferred Correction Upwind Flux
/*				if(j==ncID[e][f][g])*/
/*				{*/
/*					fuuds = convp*u[i][j][k][l] + convf*unBA[i][1][k][b];*/
/*      					fvuds = convp*v[i][j][k][l] + convf*vnBA[i][1][k][b];*/
/*      					fwuds = convp*w[i][j][k][l] + convf*wnBA[i][1][k][b];*/
/*      					*/
/*      					fuuds = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l];*/
/*      					fvuds = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l];*/
/*      					fwuds = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l];*/
/*      				}*/
/*      				else*/
/*      				{*/
      					fuuds = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l];
      					fvuds = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l];
      					fwuds = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l];
/*      				}*/
      				
      				
      			//		Deferred Correction Central Differencing Flux
				fucds = Fn[i][j][k]*(u[i][jnorth][k][l]*fyn + u[i][j][k][l]*fyp);
      				fvcds = Fn[i][j][k]*(v[i][jnorth][k][l]*fyn + v[i][j][k][l]*fyp);
      				fwcds = Fn[i][j][k]*(w[i][jnorth][k][l]*fyn + w[i][j][k][l]*fyp);
      			
      			//		Deferred Correction hybrid Flux
      				fuH = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l] - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0)*(u[i][jnorth][k][l]-u[i][j][k][l]);
      				fvH = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l] - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0)*(v[i][jnorth][k][l]-v[i][j][k][l]);
      				fwH = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l] - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0)*(w[i][jnorth][k][l]-w[i][j][k][l]);
      			
/*      			//		Deferred Correction Power Law Flux*/
/*      				fuP = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0)*(u[i][jnorth][k][l]-u[i][j][k][l]);*/
/*      				fvP = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0)*(v[i][jnorth][k][l]-v[i][j][k][l]);*/
/*      				fwP = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0)*(w[i][jnorth][k][l]-w[i][j][k][l]);*/
/*				*/


				//		Deferred Correction QUICK Flux
				fuQ = DC_QUICK_flux_north(fuQ,u);
				fvQ = DC_QUICK_flux_north(fvQ,v);
				fwQ = DC_QUICK_flux_north(fwQ,w);

				//		Deferred Correction QUICK Flux when hybrid and power Law is used to predict
				fuQHP = DC_QUICK_flux_north(fuQHP,u);				
				fuQHP = fuQHP - diff*(u[i][jnorth][k][l]-u[i][j][k][l]);
				
				fvQHP = DC_QUICK_flux_north(fvQHP,v);		
				fvQHP = fvQHP - diff*(v[i][jnorth][k][l]-v[i][j][k][l]);
				  
				fwQHP = DC_QUICK_flux_north(fwQHP,w);
				fwQHP = fwQHP - diff*(w[i][jnorth][k][l]-w[i][j][k][l]);    
      			
       			
      				if(scheme==1)	v_north_source_dc(fuuds,fvuds,fwuds,fuQ,fvQ,fwQ);
      				if(scheme==2)	v_north_source_dc(fuH,fvH,fwH,fuQHP,fvQHP,fwQHP);
				if(scheme==3)	v_north_source_dc(fuP,fvP,fwP,fuQHP,fvQHP,fwQHP); 
      
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
      				sv[i][j][k] = sv[i][j][k] - dpy[i][j][k]*vol;//*rho[i][j][k]/rhoc;
								
				
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
					sv[i][j][k]	= sv[i][j][k] + apt*v[i][j][k][l-1]; 
					apv[i][j][k]	= apv[i][j][k] + apt;
				}
				else
				{
					sv[i][j][k]	= sv[i][j][k] + (1.0 + Gamma)*apt*v[i][j][k][l-1] - 0.5*Gamma*apt*v[i][j][k][l-2]; 
					apv[i][j][k]	= apv[i][j][k] + (1.0 + 0.5*Gamma)*apt;
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
/*				su[i][j][k] = su[i][j][k] + yTempGrav*SB;*/
/*				*/
/*			}*/
/*		}*/
/*	}*/

	//Adding gravity contribution to source terms
	for(i=2;i<=nxm;i++)
	{
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
			
				/* Arithmetic Interpolation of density*/
				rhon = fyp * rho[i][j][k] +   fyn * rho[i][jnorth][k];
				rhos = (1.0-fys) * rho[i][j][k] +   fys * rho[i][jsouth][k];
				
				vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
				
				/* Arithmetic Interpolation of Temperature*/
				Tn = fyp * T[i][j][k][l] +   fyn * T[i][jnorth][k][l];
				Ts = (1.0-fys) * T[i][j][k][l] +   fys * T[i][jsouth][k][l];
				
				/* Arithmetic Interpolation of Temperature*/
				Betan = fyp * Beta[i][j][k] +   fyn * Beta[i][jnorth][k];
				Betas = (1.0-fys) * Beta[i][j][k] +   fys * Beta[i][jsouth][k];
				
				
				yBFn[i][j][k][l] = rhon*vol*yGrav - rhon*vol*yTempGrav*Betan*(Tn-Tref);			//Volume effect not yet included
				yBFs[i][j][k] = rhos*vol*yGrav - rhos*vol*yTempGrav*Betas*(Ts-Tref);
			
				//Mencinger's discontinuous Body force interpolation
				yBF[i][j][k] = fyn*yBFn[i][j][k][l] + fys*yBFs[i][j][k];
				
				//Adding Surface Tension Body force
				yBFn[i][j][k][l] = yBFn[i][j][k][l] + ySTn[i][j][k]*vol;   
				yBF[i][j][k] = yBF[i][j][k] + yST[i][j][k]*vol;
				
				//if(yST[i][j][k]!=0)	printf("%e \n",yST[i][j][k]);
				
				rhoc = (rhon + rhos)/2.0;
				
				sv[i][j][k]	= sv[i][j][k] + yBF[i][j][k];//*rho[i][j][k]/rhoc; 
			}
		}
	}


}

