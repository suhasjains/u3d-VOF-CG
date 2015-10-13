#include"preProcessing.h"

//Inputs the default value that is required for simulation
void defaultInput()
{

		xDomainLength=0.04329626;		//Length of domain along x
		yDomainLength=0.02164813;
		zDomainLength=0.000270601625;

		xNumberOfCells=100;		//Number of cells along x
		yNumberOfCells=50;
		zNumberOfCells=1;

		startTime=0;

		endTime=1000;

		nTimeSteps=1000;

		writeOutputSteps=1;
		
		adjustableCourantNumber=0;	//1=yes
		allowableCourantNumber=2.6;
		leastCourantNumber=2.4;
		
		dtau = 10e-4;

		rho1 = 6500;			//density
		rho2 = 997.1;		
		
		mu1 = 0.001003;			//viscosity	
		mu2 = 0.001003;
		
		Kcond1 = 18;			//conductivity
		Kcond2 = 0.613;
		
		Cp1 = 540;			//Specific heat
		Cp2 = 4179;
		
		Beta1 = 0.85e-5;		//volumetric expansion coefficient
		Beta2 = 2.1e-4;
		
		fixedAdvectionVelocity = 0;		//1=yes

		uTopBoundVel=0;			//at highest value of y
		vTopBoundVel=0;
		wTopBoundVel=0;
		TopBoundTemp=0;
		TopBoundType=1;			//1=wall,2=Symmetry,3=zero Gradient,4=Specified Pressure

		uBottomBoundVel=0;		//at lowest value of y
		vBottomBoundVel=0;
		wBottomBoundVel=0;
		BottomBoundTemp=1;
		BottomBoundType=1;

		uLeftBoundVel=0;		//at Lowest value of x
		vLeftBoundVel=0;
		wLeftBoundVel=0;
		LeftBoundTemp=0;
		LeftBoundType=4;

		uRightBoundVel=0;		//at Highest value of x 
		vRightBoundVel=0;
		wRightBoundVel=0;
		RightBoundTemp=0;
		RightBoundType=4;

		uFrontBoundVel=0;		//at Highest value of z
		vFrontBoundVel=0;
		wFrontBoundVel=0;
		FrontBoundTemp=0;
		FrontBoundType=2;		
		
		uBackBoundVel=0;		//at Lowest value of z
		vBackBoundVel=0;
		wBackBoundVel=0;
		BackBoundTemp=0;
		BackBoundType=2;


		pUnderRelaxCoeff=0.7;		//Pressure under relaxation factor
		wUnderRelaxCoeff=0.3;
		vUnderRelaxCoeff=0.3;
		uUnderRelaxCoeff=0.3;
		tUnderRelaxCoeff=1;

		xGrav=0;			//Body force value in terms of acceleration
		yGrav=0;
		zGrav=0;
		
		
		xTempGrav=0;			//Buoyancy due to temperature
		yTempGrav=-9.8;
		zTempGrav=0;
		
		
		DCvalue=1;			//Deferred Correction value

		vAccuracy=10e-5;
		pAccuracy=10e-6;
		alphaAccuracy=10e-12;
		massAccuracy=10e-5;
		TAccuracy=10e-5;
		pseudoTimeAccuracy=10e-5;
		
		SOR = 1.8;

		scheme=2;			//Convective scheme 1=UQ,2=HQ,3=PQ
		

		Gamma = 1;			//1-Second order time accurate, 0-first order time accurate(Euler)
		
		nFluids = 1;  			//No of Fluids, if nFluids=1 => alpha is by default 0
		sigma = 0.0;			//Surface Tension of liquid with gas at 20deg C
/*		VOFScheme=1;			//1=CICSAM , 2=STACS*/
		curvature=3;			//1=seven point smoothing , 2 = Kernel Convolution, 3 = Height Function
		capillaryTimeStep=1.0;
		
		xNB = 1;			//No of blocks along x
		yNB = 1;
		zNB = 1;
		
		Tref = 0;			//reference temperature
		ppo  = 0;			//reference pressure



}

