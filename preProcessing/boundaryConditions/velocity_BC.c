#include"boundaryConditions.h"

//Velocity boundary Conditions
void velocityBoundaryConditions()	
{
	/* Bottom boundary conditions */
  	j = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	
			diff = mu[i][j-1][k]*(xf[i]-xf[i-1])*(zf[k]-zf[k-1])/(yc[j]-yc[j-1]);
		
		
			if(BottomBoundType==1)
			{
				//Wall 
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
	    			su[i][j][k]  = su[i][j][k]  + diff*u[i][j-1][k][l];
	    			sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j-1][k][l];
	    		}
	    		
	    		else if(BottomBoundType==2)
	    		{    	
	    		
	    			//symmetry
	    			apv[i][j][k] = apv[i][j][k] + diff;
	    			u[i][j-1][k][l] = u[i][j][k][l];
	    			w[i][j-1][k][l] = w[i][j][k][l];
			}		
			
			else if(BottomBoundType==3)
			{
				//zero Gradient
				u[i][j-1][k][l] = u[i][j][k][l];
	    			w[i][j-1][k][l] = w[i][j][k][l];
				v[i][j-1][k][l] = v[i][j][k][l];
			}
			
			else if(BottomBoundType==4)
			{
				//Specified Pressure
/*				u[i][j-1][k][l] = 0;//u[i][j][k][l];*/
/*	    			w[i][j-1][k][l] = 0;//w[i][j][k][l];*/
/*				v[i][j-1][k][l] = vn[i][j-1][k][l];*/

				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*vn[i][j-1][k][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j][k][l];				
			}		
	  	}
  	  }
  	
	    
  	  /* Top boundary conditions */
  	j = nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

			diff = mu[i][j+1][k]*(xf[i]-xf[i-1])*(zf[k]-zf[k-1])/(yc[j+1]-yc[j]);
		
			if(TopBoundType==1)
			{
			
				//Wall
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
	    			su[i][j][k]  = su[i][j][k]  + diff*u[i][j+1][k][l];
	    			sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j+1][k][l];
	    		}
	    		
	    		else if(TopBoundType==2)
	    		{
	    		    	//symmetry
	    			apv[i][j][k] = apv[i][j][k] + diff;
	    			u[i][j+1][k][l] = u[i][j][k][l];
	    			w[i][j+1][k][l] = w[i][j][k][l];
			}
			
			else if(TopBoundType==3)
			{
				//zero Gradient
				u[i][j+1][k][l] = u[i][j][k][l];
	    			w[i][j+1][k][l] = w[i][j][k][l];
 				v[i][j+1][k][l] = v[i][j][k][l];
			}
			
			else if(TopBoundType==4)
			{
				//Specified Pressure
/*				u[i][j+1][k][l] = 0;//u[i][j][k][l];*/
/*	    			w[i][j+1][k][l] = 0;//w[i][j][k][l];*/
/* 				v[i][j+1][k][l] = vn[i][j][k][l];*/
				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*vn[i][j][k][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j][k][l];

 			}
			
	  	}
	  }
  
  	/* West boundary conditions */
	i = 2;
  	for (j=2;j<=nym;j++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{

			diff = mu[i-1][j][k]*(yf[j]-yf[j-1])*(zf[k]-zf[k-1])/(xc[i]-xc[i-1]);
		
			if(LeftBoundType==1)
			{
				//Wall
    				apv[i][j][k] = apv[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i-1][j][k][l];
	  			sw[i][j][k]  = sw[i][j][k]  + diff*w[i-1][j][k][l];
	  		}	
	  		
	  		else if(LeftBoundType==2)
	  		{
			  	//symmetry
	  			apu[i][j][k] = apu[i][j][k] + diff;
	  			v[i-1][j][k][l] = v[i][j][k][l];
				w[i-1][j][k][l] = w[i][j][k][l];
			}
			
			else if(LeftBoundType==3)
			{
				//zero Gradient
				v[i-1][j][k][l] = v[i][j][k][l];
				w[i-1][j][k][l] = w[i][j][k][l];
				u[i-1][j][k][l] = u[i][j][k][l];
			}
			
			else if(LeftBoundType==4)
			{
				//specified Pressure
/*				v[i-1][j][k][l] = 0;//v[i][j][k][l];*/
/*				w[i-1][j][k][l] = 0;//w[i][j][k][l];*/
/*				u[i-1][j][k][l] = ue[i-1][j][k][l];*/
				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*ue[i-1][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j][k][l];
			}
		
		//printf("%e	%e\n",yf[j]-yf[j-1],dy);
		}
    }
          
  	/* East boundary conditions */
 	i = nxm;
  	for (j=2;j<=nym;j++)
  	{
  		for(k=2;k<=nzm;k++)
  		{
			diff = mu[i+1][j][k]*(yf[j]-yf[j-1])*(zf[k]-zf[k-1])/(xc[i+1]-xc[i]);
		
			if(RightBoundType==1)
			{	
				//Wall
    				apv[i][j][k] = apv[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i+1][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*w[i+1][j][k][l];
			}
	
			else if(RightBoundType==2)
			{
				//symmetry
		  		apu[i][j][k] = apu[i][j][k] + diff;
		  		v[i+1][j][k][l] = v[i][j][k][l];
				w[i+1][j][k][l] = w[i][j][k][l];
			}
	
			else if(RightBoundType==3)
			{		
				//zero Gradient
				v[i+1][j][k][l] = v[i][j][k][l];
				w[i+1][j][k][l] = w[i][j][k][l];
				u[i+1][j][k][l] = u[i][j][k][l];
			}
	
			else if(RightBoundType==4)
			{
				//specified Pressure
/*				v[i+1][j][k][l] = 0;//v[i][j][k][l];*/
/*				w[i+1][j][k][l] = 0;//w[i][j][k][l];*/
/*				u[i+1][j][k][l] = ue[i][j][k][l];*/
				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*ue[i][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j][k][l];
			}	
	
  		}
  	}  
  
 	/* Back bounday conditions */
	k = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
			diff = mu[i][j][k-1]*(yf[j]-yf[j-1])*(xf[i]-xf[i-1])/(zc[k]-zc[k-1]);
		
			if(BackBoundType==1)
			{
				//Wall
    				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k-1][l];
	  			su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k-1][l];
			}
			
			else if(BackBoundType==2)
			{
				//symmetry
				apw[i][j][k] = apw[i][j][k] + diff;
				u[i][j][k-1][l] = u[i][j][k][l];
				v[i][j][k-1][l] = v[i][j][k][l];
			}
			
			else if(BackBoundType==3)
			{
				//zero Gradient
				u[i][j][k-1][l] = u[i][j][k][l];
				v[i][j][k-1][l] = v[i][j][k][l];
				w[i][j][k-1][l] = w[i][j][k][l];
			}
			
			else if(BackBoundType==4)
			{
				//Specified Pressure
/*				u[i][j][k-1][l] = 0;//u[i][j][k][l];*/
/*				v[i][j][k-1][l] = 0;//v[i][j][k][l];*/
/*				w[i][j][k-1][l] = wt[i][j][k-1][l];*/
				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*wt[i][j][k-1][l];
			}		

		}
    	}
          
  	/* Front boundary conditions */
  	k = nzm;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{

			diff = mu[i][j][k+1]*(yf[j]-yf[j-1])*(xf[i]-xf[i-1])/(zc[k+1]-zc[k]);
			
			if(FrontBoundType==1)
			{
				//Wall
	    			apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k+1][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k+1][l];
			}	
		
			else if(FrontBoundType==2)
			{
				//symmetry
				apw[i][j][k] = apw[i][j][k] + diff;
				u[i][j][k+1][l] = u[i][j][k][l];
				v[i][j][k+1][l] = v[i][j][k][l];
			}
			
			else if(FrontBoundType==3)
			{
				//Zero Gradient
				u[i][j][k+1][l] = u[i][j][k][l];
				v[i][j][k+1][l] = v[i][j][k][l];
				w[i][j][k+1][l] = w[i][j][k][l]; 
			}
			
			else if(FrontBoundType==4)
			{
				//Specified Pressure
/*				u[i][j][k+1][l] = u[i][j][k][l];*/
/*				v[i][j][k+1][l] = v[i][j][k][l];*/
/*				w[i][j][k+1][l] = wt[i][j][k][l]; */
				apv[i][j][k] = apv[i][j][k] + diff;
    				apu[i][j][k] = apu[i][j][k] + diff;
    				apw[i][j][k] = apw[i][j][k] + diff;
    				sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k][l];
		  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k][l];
				sw[i][j][k]  = sw[i][j][k]  + diff*wt[i][j][k][l];
		  		
			}
		}
    	}
}

void setting_velocity()
{
  j=ny;						// Top boundary 
  for (i=2;i<=nxm;i++)
   for(k=2;k<=nzm;k++)
   {
    u[i][j][k][2] = uTopBoundVel;
    w[i][j][k][2] = wTopBoundVel;
   }
   
   j=1;						// Bottom boundary 
  for (i=2;i<=nxm;i++)
   for(k=2;k<=nzm;k++)
   {
    u[i][j][k][2] = uBottomBoundVel;
    w[i][j][k][2] = wBottomBoundVel;
   }
   
   k=nz;					// Front boundary 
  for (i=2;i<=nxm;i++)
   for(j=2;j<=nym;j++)
   {
    u[i][j][k][2] = uFrontBoundVel;
    v[i][j][k][2] = vFrontBoundVel;
   }
   
   k=1;						// Back boundary 
  for (i=2;i<=nxm;i++)
   for(j=2;j<=nym;j++)
   {
    u[i][j][k][2] = uBackBoundVel;
    v[i][j][k][2] = vBackBoundVel;
   }
   
   
   i=nx;					// Right boundary 
  for (j=2;j<=nym;j++)
   for(k=2;k<=nzm;k++)
   {
    v[i][j][k][2] = vRightBoundVel;
    w[i][j][k][2] = wRightBoundVel;
   }
   
   i=1;						// Left boundary 
  for (j=2;j<=nym;j++)
   for(k=2;k<=nzm;k++)
   {
    v[i][j][k][2] = vLeftBoundVel;
    w[i][j][k][2] = wLeftBoundVel;
   }
   
   
   
   
   
}

