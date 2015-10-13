#include"velocity.h"

//calculates the value of initial v velocity
void vEqn()
{
	 for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
		 		ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apv[i][j][k])*urecrv;
      			sv[i][j][k]  = sv[i][j][k] + (1 - vUnderRelaxCoeff) * ap[i][j][k] * v[i][j][k][l] ;
/*      			sv[i][j][k]	= sv[i][j][k] + rho[i][j][k]*vol*yGrav;*/
      			apv[i][j][k] = 1 / ap[i][j][k];	/* For implementation in the pressure correction approach */
		 	
		 	//printf("apv = %e\n",apv[i][j][k] );	
		 	}
 	     }
  	 }
  	 
  for(nonLinear=0;nonLinear<=vVelocityLoops;nonLinear++)
  {
  	
//  Creating matrices for ADI-TDMA algorithm
	for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
	{
		for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[j] 	= sv[i][j][k];
				
/*				if(k==tcID[e][f][g])	c[j] = c[j] - at[i][j][k]*vtBA[i][j][1][b];*/
/*				else			c[j] = c[j] - at[i][j][k]*v[i][j][ktop][l];*/
/*				if(k==bcID[e][f][g])	c[j] = c[j] - ab[i][j][k]*vbBA[i][j][1][b];*/
/*				else			c[j] = c[j] - ab[i][j][k]*v[i][j][kbottom][l];*/
/*				if(i==wcID[e][f][g])	c[j] = c[j] - aw[i][j][k]*vwBA[1][j][k][b];*/
/*				else			c[j] = c[j] - aw[i][j][k]*v[iwest][j][k][l];*/
/*				if(i==ecID[e][f][g])	c[j] = c[j] - ae[i][j][k]*veBA[1][j][k][b];*/
/*				else			c[j] = c[j] - ae[i][j][k]*v[ieast][j][k][l];*/
				
				c[j] = c[j] - at[i][j][k]*v[i][j][ktop][l];
				c[j] = c[j] - ab[i][j][k]*v[i][j][kbottom][l];
				c[j] = c[j] - aw[i][j][k]*v[iwest][j][k][l];
				c[j] = c[j] - ae[i][j][k]*v[ieast][j][k][l];
				
				c[j] = c[j]*apv[i][j][k]; 
				
				
							
				lower[j]= as[i][j][k]*apv[i][j][k];
				upper[j]= an[i][j][k]*apv[i][j][k];
				
				
			}
				mThomas(scID[e][f][g],ncID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=scID[e][f][g];thomasI<=ncID[e][f][g];thomasI++)
	  				v[i][thomasI][k][l]=x[thomasI];				
			
    	}
   }
  
  
	vMaxRes = 0;
	vTotalRes = 0;  
  	for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
	{
		for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
		{
			for(j=scID[e][f][g];j<=ncID[e][f][g];j++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
	  			vRes[i][j][k] 	= sv[i][j][k] - (an[i][j][k]*v[i][jnorth][k][l] + as[i][j][k]*v[i][jsouth][k][l] + at[i][j][k]*v[i][j][ktop][l] + ab[i][j][k]*v[i][j][kbottom][l] + ae[i][j][k]*v[ieast][j][k][l] + aw[i][j][k]*v[iwest][j][k][l] + ap[i][j][k]*v[i][j][k][l]);
	  			
	  			if(j==ncID[e][f][g])	vRes[i][j][k] = vRes[i][j][k] - an[i][j][k]*vnBA[i][1][k][b];
				else			vRes[i][j][k] = vRes[i][j][k] - an[i][j][k]*v[i][jnorth][k][l];
				if(j==scID[e][f][g])	vRes[i][j][k] = vRes[i][j][k] - as[i][j][k]*vsBA[i][1][k][b];
				else			vRes[i][j][k] = vRes[i][j][k] - as[i][j][k]*v[i][jsouth][k][l];
				if(k==tcID[e][f][g])	vRes[i][j][k] = vRes[i][j][k] - at[i][j][k]*vtBA[i][j][1][b];
				else			vRes[i][j][k] = vRes[i][j][k] - at[i][j][k]*v[i][j][ktop][l];
				if(k==bcID[e][f][g])	vRes[i][j][k] = vRes[i][j][k] - ab[i][j][k]*vbBA[i][j][1][b];
				else			vRes[i][j][k] = vRes[i][j][k] - ab[i][j][k]*v[i][j][kbottom][l];
				if(i==wcID[e][f][g])	vRes[i][j][k] = vRes[i][j][k] - aw[i][j][k]*vwBA[1][j][k][b];
				else			vRes[i][j][k] = vRes[i][j][k] - aw[i][j][k]*v[iwest][j][k][l];
				if(i==ecID[e][f][g])	vRes[i][j][k] = vRes[i][j][k] - ae[i][j][k]*veBA[1][j][k][b];
				else			vRes[i][j][k] = vRes[i][j][k] - ae[i][j][k]*v[ieast][j][k][l];
	  			
	  			
				vTotalRes += fabs(ap[i][j][k]*v[i][j][k][l]);
				if(fabs(vRes[i][j][k])>vMaxRes)	vMaxRes	= fabs(vRes[i][j][k]);
	  			
			}
    		}
   	}
	
	vMeanRes = vTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
	if(vMeanRes==0)		vNormRes = 0;
    	else			vNormRes = vMaxRes/vMeanRes;
   
	   
   	//sync();
	//printf(" v Res = %e\n",vNormRes);
   	//if(fabs(vNormRes)>vAccuracy)	vVelocityLoops = vVelocityLoops + 1;
   
 }
   
}
