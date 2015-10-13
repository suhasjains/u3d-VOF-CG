//calculates the value of initial u velocity
void uEqn()
{
	
	for (i=wcID[e][f][g];i<=ecID[e][f][g];i++)
	{
		for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
		{
			for(k=bcID[e][f][g];k<=tcID[e][f][g];k++)
			{
				ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apu[i][j][k])*urecru;
      				su[i][j][k]  = su[i][j][k] + (1 - uUnderRelaxCoeff) * ap[i][j][k] * u[i][j][k][l];
/*      				su[i][j][k]  = su[i][j][k] + rho[i][j][k]*vol*xGrav;*/
      				apu[i][j][k] = 1 / ap[i][j][k];	/* For implementation in the pressure correction approach */
				
			}
    		}
  	}


for(nonLinear=0;nonLinear<=uVelocityLoops;nonLinear++)
{

// Creating matrices for ADI-TDMA algorithm  
  for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
	{
		for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
		{
			for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[i] 	= su[i][j][k];
				 
/*				if(j==ncID[e][f][g])	c[i] = c[i] - an[i][j][k]*unBA[i][1][k][b];*/
/*				else			c[i] = c[i] - an[i][j][k]*u[i][jnorth][k][l];*/
/*				if(j==scID[e][f][g])	c[i] = c[i] - as[i][j][k]*usBA[i][1][k][b];*/
/*				else			c[i] = c[i] - as[i][j][k]*u[i][jsouth][k][l];*/
/*				if(k==tcID[e][f][g])	c[i] = c[i] - at[i][j][k]*utBA[i][j][1][b];*/
/*				else			c[i] = c[i] - at[i][j][k]*u[i][j][ktop][l];*/
/*				if(k==bcID[e][f][g])	c[i] = c[i] - ab[i][j][k]*ubBA[i][j][1][b];*/
/*				else			c[i] = c[i] - ab[i][j][k]*u[i][j][kbottom][l];*/
/*				*/


				c[i] = c[i] - an[i][j][k]*u[i][jnorth][k][l];
				c[i] = c[i] - as[i][j][k]*u[i][jsouth][k][l];
				c[i] = c[i] - at[i][j][k]*u[i][j][ktop][l];
				c[i] = c[i] - ab[i][j][k]*u[i][j][kbottom][l];
					
				
				c[i] = c[i]*apu[i][j][k];
				
				
				lower[i]= aw[i][j][k]*apu[i][j][k];
				upper[i]= ae[i][j][k]*apu[i][j][k];
				
			}
				
				mThomas(wcID[e][f][g],ecID[e][f][g]);
				
				//Equating solved value back to velocity
				for(thomasI=wcID[e][f][g];thomasI<=ecID[e][f][g];thomasI++)
	  				u[thomasI][j][k][l]=x[thomasI];				
			
    	}
  }
  
  
  
  uMaxRes = 0;
  uTotalRes = 0;
  for (j=scID[e][f][g];j<=ncID[e][f][g];j++)
	{
		for (k=bcID[e][f][g];k<=tcID[e][f][g];k++)
		{
			for(i=wcID[e][f][g];i<=ecID[e][f][g];i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				uRes[i][j][k] 	= su[i][j][k] - (ap[i][j][k]*u[i][j][k][l]);
				
				if(j==ncID[e][f][g])	uRes[i][j][k] = uRes[i][j][k] - an[i][j][k]*unBA[i][1][k][b];
				else			uRes[i][j][k] = uRes[i][j][k] - an[i][j][k]*u[i][jnorth][k][l];
				if(j==scID[e][f][g])	uRes[i][j][k] = uRes[i][j][k] - as[i][j][k]*usBA[i][1][k][b];
				else			uRes[i][j][k] = uRes[i][j][k] - as[i][j][k]*u[i][jsouth][k][l];
				if(k==tcID[e][f][g])	uRes[i][j][k] = uRes[i][j][k] - at[i][j][k]*utBA[i][j][1][b];
				else			uRes[i][j][k] = uRes[i][j][k] - at[i][j][k]*u[i][j][ktop][l];
				if(k==bcID[e][f][g])	uRes[i][j][k] = uRes[i][j][k] - ab[i][j][k]*ubBA[i][j][1][b];
				else			uRes[i][j][k] = uRes[i][j][k] - ab[i][j][k]*u[i][j][kbottom][l];
				if(i==wcID[e][f][g])	uRes[i][j][k] = uRes[i][j][k] - aw[i][j][k]*uwBA[1][j][k][b];
				else			uRes[i][j][k] = uRes[i][j][k] - aw[i][j][k]*u[iwest][j][k][l];
				if(i==ecID[e][f][g])	uRes[i][j][k] = uRes[i][j][k] - ae[i][j][k]*ueBA[1][j][k][b];
				else			uRes[i][j][k] = uRes[i][j][k] - ae[i][j][k]*u[ieast][j][k][l];
				
				
				//uTotalRes += fabs(uRes[i][j][k]);
				uTotalRes += fabs(ap[i][j][k]*u[i][j][k][l]);
				if(fabs(uRes[i][j][k])>uMaxRes)	uMaxRes	= fabs(uRes[i][j][k]);
			}
    		}
   	}
	
	uMeanRes = uTotalRes/((double)((ecID[e][f][g] - wcID[e][f][g] + 1)*(ncID[e][f][g] - scID[e][f][g] + 1)*(tcID[e][f][g] - bcID[e][f][g] + 1)));
	if(uMeanRes==0)		uNormRes = 0;
    	else			uNormRes = uMaxRes/uMeanRes;
   
	
	//sync();
	
	//printf(" u Res = %e\n",uNormRes);
	//printf(" u Max Res = %e\n",uTotalRes);
   	//if(fabs(uNormRes)>vAccuracy)		uVelocityLoops = uVelocityLoops + 1;
   	
   	
}
  
}

