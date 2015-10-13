#include"velocity.h"

//calculates the value of initial w velocity
void wEqn()
{
	for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
		 		ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apw[i][j][k])*urecrw;
      				sw[i][j][k]  = sw[i][j][k] + (1 - wUnderRelaxCoeff) * ap[i][j][k] * w[i][j][k][l] ;
/*      				sw[i][j][k]	= sw[i][j][k] + rho[i][j][k]*vol*zGrav;*/
      				apw[i][j][k] = 1 / ap[i][j][k];	/* For implementation in the pressure correction approach */
		 		
		 	}
 	     }
  	 }
  
   for(nonLinear=0;nonLinear<=wVelocityLoops;nonLinear++)
  {
  
//  Creating matrices for ADI-TDMA algorithm	 
  	 for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[k] 	= sw[i][j][k];
				
/*				if(i==wcID[e][f][g])	c[k] = c[k] - aw[i][j][k]*wwBA[1][j][k][b];*/
/*				else			c[k] = c[k] - aw[i][j][k]*w[iwest][j][k][l];*/
/*				if(i==ecID[e][f][g])	c[k] = c[k] - ae[i][j][k]*weBA[1][j][k][b];*/
/*				else			c[k] = c[k] - ae[i][j][k]*w[ieast][j][k][l];*/
/*				if(j==ncID[e][f][g])	c[k] = c[k] - an[i][j][k]*wnBA[i][1][k][b];*/
/*				else			c[k] = c[k] - an[i][j][k]*w[i][jnorth][k][l];*/
/*				if(j==scID[e][f][g])	c[k] = c[k] - as[i][j][k]*wsBA[i][1][k][b];*/
/*				else			c[k] = c[k] - as[i][j][k]*w[i][jsouth][k][l];*/
/*				*/
				
				c[k] = c[k] - aw[i][j][k]*w[iwest][j][k][l];
				c[k] = c[k] - ae[i][j][k]*w[ieast][j][k][l];
				c[k] = c[k] - an[i][j][k]*w[i][jnorth][k][l];
				c[k] = c[k] - as[i][j][k]*w[i][jsouth][k][l];
				
				
				c[k] = c[k]*apw[i][j][k];
							
				
				lower[k]= ab[i][j][k]*apw[i][j][k];
				upper[k]= at[i][j][k]*apw[i][j][k];
				
				
			}
				mThomas(bcID[e][f][g],tcID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=bcID[e][f][g];thomasI<=tcID[e][f][g];thomasI++)
	  				w[i][j][thomasI][l]=x[thomasI];				
			
    	}
  }
  
  wMaxRes = 0;
  wTotalRes = 0;
  for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
	  			wRes[i][j][k] 	= sw[i][j][k] - ap[i][j][k]*w[i][j][k][l];
	  			
	  			if(j==ncID[e][f][g])	wRes[i][j][k] = wRes[i][j][k] - an[i][j][k]*wnBA[i][1][k][b];
				else			wRes[i][j][k] = wRes[i][j][k] - an[i][j][k]*w[i][jnorth][k][l];
				if(j==scID[e][f][g])	wRes[i][j][k] = wRes[i][j][k] - as[i][j][k]*wsBA[i][1][k][b];
				else			wRes[i][j][k] = wRes[i][j][k] - as[i][j][k]*w[i][jsouth][k][l];
				if(k==tcID[e][f][g])	wRes[i][j][k] = wRes[i][j][k] - at[i][j][k]*wtBA[i][j][1][b];
				else			wRes[i][j][k] = wRes[i][j][k] - at[i][j][k]*w[i][j][ktop][l];
				if(k==bcID[e][f][g])	wRes[i][j][k] = wRes[i][j][k] - ab[i][j][k]*wbBA[i][j][1][b];
				else			wRes[i][j][k] = wRes[i][j][k] - ab[i][j][k]*w[i][j][kbottom][l];
				if(i==wcID[e][f][g])	wRes[i][j][k] = wRes[i][j][k] - aw[i][j][k]*wwBA[1][j][k][b];
				else			wRes[i][j][k] = wRes[i][j][k] - aw[i][j][k]*w[iwest][j][k][l];
				if(i==ecID[e][f][g])	wRes[i][j][k] = wRes[i][j][k] - ae[i][j][k]*weBA[1][j][k][b];
				else			wRes[i][j][k] = wRes[i][j][k] - ae[i][j][k]*w[ieast][j][k][l];
	  			
	  			
				wTotalRes += fabs(ap[i][j][k]*w[i][j][k][l]);
				if(fabs(wRes[i][j][k])>wMaxRes)	wMaxRes	= fabs(wRes[i][j][k]);
	  			
			}
    	}
   }
   
   	wMeanRes = wTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
    	if(wMeanRes==0)		wNormRes = 0;
    	else			wNormRes = wMaxRes/wMeanRes;
   
   
   	//sync();
   
	//printf(" w Res = %e\n",wNormRes);
   	//if(fabs(wNormRes)>vAccuracy)		wVelocityLoops = wVelocityLoops + 1;
 }
  

}

