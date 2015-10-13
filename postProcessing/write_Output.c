#include"write_Output.h"

//Outputs values to files
void writeOutput()
{


	//output file names
	char str[30];
	sprintf(str,"postProcessing/output/output%d.plt",t);
	
	fp=fopen(str,"w");  //opening the output file

	fprintf(fp,"title = \"Lid driven cavity\" \n");
	fprintf(fp,"variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"p\",\"T\",\"Fe\",\"Fn\",\"Ft\",\"rho\",\"mu\",\"alpha\",\"alphaS\",\"sdf\",\"xST\",\"yST\",\"zST\"\n");
	fprintf(fp,"zone i=%d, j=%d, k=%d, f=point\n",nxm+1,nym+1,nzm+1);

	for  (k=1; k<=nz; k++)
	 for (j=1; j<=ny; j++)
	  for(i=1; i<=nx; i++)
		fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",xc[i],yc[j],zc[k],u[i][j][k][l],v[i][j][k][l],w[i][j][k][l],p[i][j][k][l],T[i][j][k][l],Fe[i][j][k],Fn[i][j][k],Ft[i][j][k],rho[i][j][k],mu[i][j][k],alpha[i][j][k][l],alphaS[i][j][k][l],sdf[i][j][k],xST[i][j][k],yST[i][j][k],zST[i][j][k]);

	fclose(fp);


	sprintf(str,"postProcessing/output/output%d.csv",t);
	
	fp=fopen(str,"w");  //opening the output file

	fprintf(fp,"\"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"p\",\"T\",\"Fe\",\"Fn\",\"Ft\",\"rho\",\"mu\",\"alpha\",\"alphaS\",\"sdf\",\"xST\",\"yST\",\"zST\"\n");

	for  (k=1; k<=nz; k++)
	 for (j=1; j<=ny; j++)
	  for(i=1; i<=nx; i++)
		fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",xc[i],yc[j],zc[k],u[i][j][k][l],v[i][j][k][l],w[i][j][k][l],p[i][j][k][l],T[i][j][k][l],Fe[i][j][k],Fn[i][j][k],Ft[i][j][k],rho[i][j][k],mu[i][j][k],alpha[i][j][k][l],alphaS[i][j][k][l],sdf[i][j][k],xST[i][j][k],yST[i][j][k],zST[i][j][k]);	    


	fclose(fp);

}


