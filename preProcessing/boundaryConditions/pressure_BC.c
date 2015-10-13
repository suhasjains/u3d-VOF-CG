#include"boundaryConditions.h"

void boundary_pressurep()
{
  int js = 1, jn = ny;
  int iw = 1, ie = nx;
  int kb = 1, kt = nz;
  
  double N,D;
  
  /* Top and Bottom pressure boundary conditions */
  for (i=2;i<=nxm;i++)
  {
	for (k=2;k<=nzm;k++)
  	{
  	
		//South  	
  		j = 2;
  		
		jsouth  = j - 1;
		iwest 	= i - 1;
		kbottom	= k - 1;

		s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
		vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);


		N = Fn[i][jsouth][k] - rho[i][j][k]*v[i][j][k][l]*s + 0.5*(Ft[i][j][kbottom] - Ft[i][j][k] + Fe[iwest][j][k] - Fe[i][j][k]);
		D = rho[i][j][k]*apv[i][j][k]*vol;
		
		//printf("D = %e\n",N/D);

		if(BottomBoundType==4)	pp[i][js][k][l] = 0;
		else 			pp[i][js][k][l] = pp[i][js+1][k][l] + N/D;
		
		
		
		
		//North
		j = nym;
  		
		jsouth  = j - 1;
		iwest 	= i - 1;
		kbottom	= k - 1;

		s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
		vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
		
		
		
		N = rho[i][j][k]*v[i][j][k][l]*s - Fn[i][j][k] + 0.5*(Ft[i][j][kbottom] - Ft[i][j][k] + Fe[iwest][j][k] - Fe[i][j][k]);
		D = rho[i][j][k]*apv[i][j][k]*vol;

		
   		if(TopBoundType==4)	pp[i][jn][k][l] = 0;
   		else 			pp[i][jn][k][l] = pp[i][jn-1][k][l] + N/D;
   		
   	
   	}
  }
  
  /* East adn West boundary conditions */
  for (j=2;j<=nym;j++)
  {
  	for (k=2;k<=nzm;k++)
  	{
  		//West
  		i = 2;
  	
  		jsouth  = j - 1;
		iwest 	= i - 1;
		kbottom	= k - 1;
		
		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
		vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
		
		N = Fe[iwest][j][k] - rho[i][j][k]*u[i][j][k][l]*s + 0.5*(Fn[i][jsouth][k] -Fn[i][j][k] + Ft[i][j][kbottom] - Ft[i][j][k]); 
  		D = rho[i][j][k]*apu[i][j][k]*vol;
  	
  	
  		if(LeftBoundType==4)	pp[iw][j][k][l] = 0;
  		else 			pp[iw][j][k][l] = pp[iw+1][j][k][l] + N/D;
  	
  		
  		//East
  		i = nxm;
  		
  		jsouth  = j - 1;
		iwest 	= i - 1;
		kbottom	= k - 1;
		
		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
		vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
  		
  		
  		N = rho[i][j][k]*u[i][j][k][l]*s - Fe[i][j][k] + 0.5*( Fn[i][jsouth][k] -Fn[i][j][k] + Ft[i][j][kbottom] - Ft[i][j][k]); 
  		D = rho[i][j][k]*apu[i][j][k]*vol;
  		
  		if(RightBoundType==4)	pp[ie][j][k][l] = 0;
     		else			pp[ie][j][k][l] = pp[ie-1][j][k][l] + N/D;

     								
     	}
  }
  
  //Front and Back boundary Conditions
  for (j=2;j<=nym;j++)
  {
  	for (i=2;i<=nxm;i++)
  	{
  		//Back
  		k =2;
  		
  		jsouth  = j - 1;
		iwest 	= i - 1;
		kbottom	= k - 1;
		
		s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
		vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
  	
  		
  		N = 0.5*(Fe[iwest][j][k] - Fe[i][j][k] + Fn[i][jsouth][k] - Fn[i][j][k]) + Ft[i][j][kbottom] - rho[i][j][k]*w[i][j][k][l]*s; 
  		D = rho[i][j][k]*apw[i][j][k]*vol;
  				
  		
  		if(BackBoundType==4)	pp[i][j][kb][l] = 0;
     		else 			pp[i][j][kb][l] = pp[i][j][kb+1][l] + N/D;
     		
     		
     		
     		//Front
     		k = nzm;
     		
     		jsouth  = j - 1;
		iwest 	= i - 1;
		kbottom	= k - 1;
		
		s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
		vol = (xf[i]-xf[i-1])*(zf[k]-zf[k-1])*(yf[j] - yf[j-1]);
     		
     		
     		N = 0.5*(Fe[iwest][j][k] - Fe[i][j][k] + Fn[i][jsouth][k] - Fn[i][j][k]) + rho[i][j][k]*w[i][j][k][l]*s - Ft[i][j][k] ;
     		D = rho[i][j][k]*apw[i][j][k]*vol;
     		
     		if(FrontBoundType==4)	pp[i][j][kt][l] = 0;
     		else			pp[i][j][kt][l] = pp[i][j][kt-1][l] + N/D;
     		
  	
  	}
  }
}


void ppBoundary_Conditions()
{

	/* Bottom boundary conditions */
  	j = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
		
			as[i][j][k] = 0;
	  	}
  	  }
  	
	    
  	  /* Top boundary conditions */
  	j = nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  	
	  		an[i][j][k] = 0;

	  	}
	  }
  
  	/* West boundary conditions */
	i = 2;
  	for (j=2;j<=nym;j++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
			aw[i][j][k] = 0;
		}
    }
          
  /* East boundary conditions */
  i = nxm;
  for (j=2;j<=nym;j++)
  {
  	for(k=2;k<=nzm;k++)
  	{
			ae[i][j][k] = 0;
  	}
  }  
  
  /* Back bounday conditions */
	k = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
				
			ab[i][j][k] = 0;
		}
    }
          
  /* Front boundary conditions */
  k = nzm;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
			at[i][j][k] = 0;
		}
    }
	

}

void boundary_pressure()
{

  int js = 1, jn = ny;
  int iw = 1, ie = nx;
  int kb = 1, kt = nz;

	
  for (i=2;i<=nxm;i++)
  for (k=2;k<=nzm;k++)
  {
/*	p[i][js][k][l] = p[i][js+1][k][l];*/
/*   	p[i][jn][k][l] = p[i][jn-1][k][l];*/
   	if(BottomBoundType!=4)	p[i][js][k][l] = p[i][js][k][l] + pUnderRelaxCoeff*(pp[i][js][k][l]-ppo);
   	if(TopBoundType!=4)	p[i][jn][k][l] = p[i][jn][k][l] + pUnderRelaxCoeff*(pp[i][jn][k][l]-ppo); 
   	
  }

	
	/* East adn West boundary conditions */
  for (j=2;j<=nym;j++)
  for (k=2;k<=nzm;k++)
  {
/*     	p[iw][j][k][l] = p[iw+1][j][k][l];*/
/*     	p[ie][j][k][l] = p[ie-1][j][k][l];*/
     
     	if(LeftBoundType!=4)	p[iw][j][k][l] = p[iw][j][k][l] + pUnderRelaxCoeff*(pp[iw][j][k][l]-ppo);
     	if(RightBoundType!=4)	p[ie][j][k][l] = p[ie][j][k][l] + pUnderRelaxCoeff*(pp[ie][j][k][l]-ppo);
  }
  
  //Front and Back boundary Conditions
  for (j=2;j<=nym;j++)
  for (i=2;i<=nxm;i++)
  {
/*     	p[i][j][kb][l] = p[i][j][kb+1][l];*/
/*     	p[i][j][kt][l] = p[i][j][kt-1][l];*/
     	if(BackBoundType!=4)	p[i][j][kb][l] = p[i][j][kb][l] + pUnderRelaxCoeff*(pp[i][j][kb][l]-ppo);
     	if(FrontBoundType!=4)	p[i][j][kt][l] = p[i][j][kt][l] + pUnderRelaxCoeff*(pp[i][j][kt][l]-ppo);
  }


}

