#include<iostream>
#include<fstream>
#include<math.h>
#include<complex>
#include<cmath>
#include<stdlib.h>
#include<iomanip>
#include<map>
#include<vector>
#include<algorithm>
#include <list>

using namespace std;


//Preliminary
double GeV = 1.;
double TeV = 1000. * GeV;
double mV;
double f;
double gX;
double sigma0;
double Teq;
double T;
int Ncft;
double mplbar = 2.43536*1.e18;
double vhQCDref;
double TcQCDref;
double deltabeta;
int nT;
int ngX;

double v;
double m;
double kinIntNoDecay;
double kinIntWithDecay;




double iterateLog(int i, int n, double ini, double fin)
{
    return ini*pow(10,double(i)/n*log10(fin/ini));
}

double iterateLin(int i, int n, double ini, double fin)
{
    return ini+(fin-ini)*double(i)/n;
}



//We implement the rk5 iteration for the integration of ode
//We also implement a rk4 iteration embedded in the rk5 one. Then, we can compute the error and adapt the step size
//For more details: http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c16-2.pdf
void rk5(void (&system)(double*,double,double*), double* x, double &t, double &dt, double &error, double desiredErrorMin, double desiredErrorMax)
{
    
    //The iteration of yn is
    //yn+1=yn+c0*k0+c1*k1+c2*k2+c3*k3+c4*k4+c5*k5
    //where
    //k0=dt*f(xn,yn)
    //k1=dt*f(xn+a1*dt,yn+b10*k0)
    //k2=dt*f(xn+a2*dt,yn+b20*k0+b21*k1)
    //k3=dt*f(xn+a3*dt,yn+b30*k0+b31*k1+b32*k2)
    //k4=dt*f(xn+a4*dt,yn+b40*k0+b41*k1+b42*k2+b43*k3)
    //k5=dt*f(xn+a5*dt,yn+b50*k0+b51*k1+b52*k2+b53*k3+b54*k4)
    //an, bn and cn are called Butcher's coefficients
    
    //Temporary variables
    int i,j,p;

    double k[6][2];

    double xm[7][2];
    
    double dtTmp=dt;

    //First derivative of the ODE
    double dxdt[6][2];

    //Butcher's coefficients

    double a[6]={0.,1/5.,3/10.,3/5.,1.,7/8.};
    double b[6][5]={{0.,0.,0.,0.,0.},{1/5.,0.,0.,0.,0.},{3/40.,9/40.,0.,0.,0.},{3/10.,-9/10.,6/5.,0.,0.},{-11/54.,5/2.,-70/27.,35/27.,0.},{1631/55296.,175/512.,575/13824.,44275/110592.,253/4096.}};

    //Fifth-order Runge-Kuta method
    double c5[6]={37/378.,0.,250/621.,125/594.,0.,512/1771.};
    //Embedded four-order Runge kutta method
    double c4[6]={2825/27648.,0.,18575/48384.,13525/55296.,277/14336.,1/4.};


    //RK4 solution variable
    double x4[2];
    //RK5 solution variable
    double x5[2];

    int preventInstableError=0;

    do
    {
        dt = dtTmp;

        for(i=0; i<2; i++)
        {
            x4[i]=x[i];
            x5[i]=x[i];
        }

        for(i=0; i<2; i++)
        {
            xm[0][i]=x[i];
        }

        for(j=1 ; j<=6; j++)
        {
            system(xm[j-1],t+dt*a[j-1],dxdt[j-1]);

            for(i=0; i<2; i++)
            {
                xm[j][i]=x[i];
                k[j-1][i]=dt*dxdt[j-1][i];

                for(p=0; p<j && j<6 ; p++)
                {
                    xm[j][i]=xm[j][i]+b[j][p]*k[p][i];
                }
            }
        }

        for(i=0; i<2; i++)
        {
            for(j=0; j<6; j++)
            {
                x4[i]=x4[i]+c4[j]*k[j][i];
            }
        }

        for(i=0; i<2; i++)
        {
            for(j=0; j<6; j++)
            {
                x5[i]=x5[i]+c5[j]*k[j][i];
            }
        }

        //We set the "error" as the result of the difference between rk5 solution and rk4 solution
        error=fabs(pow(x5[1]*x5[1],1./2.)-pow(x4[1]*x4[1],1./2.))/pow(x5[1]*x5[1],1./2.);

        //Since error terms in rk4 and rk5 are respectly of order 5 and 6 according the step dt, the "error" we just computed scales as the step dt to the power of 5
        //Then, the new step as to be estimated as the following: newStep=oldStep*pow(desiredError/oldError,1/5.)

        if(error>desiredErrorMax)
            dtTmp=dt*pow(desiredErrorMax/error,0.2);

        if(error<desiredErrorMin && error!=0.)
            dtTmp=dt*pow(desiredErrorMin/error,0.2);

        //NOTE: DUE TO THE FACT THAT SOMETIMES ERROR WAS INSTABLE, THE PREVIOUS OPERATIONS MIGHT BE WRONG
        //That's why, we do the following: we consider that the error is correct after 3 passes through the loop

        preventInstableError++;
        if(preventInstableError>=500)
        {
            error=desiredErrorMax;
        }

    }


    //We do not exit from the loop while error is not the one expected
    while(error>desiredErrorMax || (error<desiredErrorMin && error!=0.));


    //We keep the rk5 solution
    for(i=0; i<2; i++)
    {
        x[i]=x5[i];
    }


  
//    double xp[2];
//    xp[0]=x[0];
//    xp[1]=x[1];
//    system(x,t,xp);
//    x[0]=x[0]+xp[0]*dt;
//    x[1]=x[1]+xp[1]*dt;
    
    
    //We increment time
    t=t+dt;
    
}


//beta function of scale-invariant SU(2)_D model
//ArXiv: 1306.2329, 1805.01473, 1809.11129, 2210.07075
double betaFunSU2(double sigma, double gX, double f)
{
    //return 9*pow(gX,4.)/16./(8*M_PI*M_PI);

    double Lhsigma = 2*0.13*pow(246/f,2);
    double fac = sqrt(3)*sqrt(9*pow(gX,4)+64*Lhsigma*Lhsigma);
    double logmu = log(sigma / exp(1./4.) / f);
    return (9*gX*gX + fac*tan(fac/(32*M_PI*M_PI)*logmu-atan(9*gX*gX/fac)))/(48.)/logmu;
}

//beta function of scale-invariant U(1)_D model
//ArXiv: 2311.13640
double betaFun(double sigma, double gX, double f)
{
    //return 9*pow(gX,4.)/16./(8*M_PI*M_PI);

    double Lhsigma = 2*0.13*pow(246/f,2);
    double fac = sqrt(21*pow(gX,4)+10*Lhsigma*Lhsigma);
    double logmu = log(sigma / exp(1./4.) / f);
    return (3*gX*gX + fac*tan(fac/(8*M_PI*M_PI)*logmu-atan(3*gX*gX/fac)))/(10.)/logmu;
}


//dilaton potential at zero temperature
double dilPotential(double sigma, double gX, double Ncft, double f)
{
    double beta = betaFun(sigma, gX, f);
    
    sigma = sigma + 1e-20;
    
    return beta / 4. * pow(sigma,4) * log(sigma / exp(1./4.) / f);
    
}

//J function for the finite-temperature correction
double JbHT(double x)
{
    return -pow(M_PI,4)/45. + pow(M_PI,2)/12. * x - M_PI/6. * pow(x, 3./2.) - x * x/32. * (log(x) - 5.40762);
}

double JLT(double x)
{
    return -sqrt(M_PI/2)*pow(x, 3./4.)*exp(-sqrt(x))*(1 + 15./8./sqrt(x) + 105./128./x);
}

double Jbfun(double x)
{
    return (exp(-pow(x/6.3, 4)) * JbHT(x) + (1 - exp(-pow(x/6.3, 4)))*JLT(x));
}


//Dilaton potential at finite-temperarure in the scale-invariant SU(2)_D model
double dilPotentialHTSU2(double sigma, double T, double gX, double Ncft, double f)
{
    sigma = sigma + 1e-20;
    
    double vH = 246.*GeV;
    double lambdaH = 0.13;
    double lambdaS =  betaFunSU2(sigma, gX, f);
    double QCDswitch = (1-tanh(20*(T-TcQCDref)))/2.;
    double mV = gX * sigma / 2.;
    double PiV = 5./6.*pow(gX*T,2);
    double lambdaHS  = 2*lambdaH*pow(vH/f,2);
    double xV = pow(mV / T, 2);
    double xG = lambdaS * pow( sigma / T, 2);
    double xS = 3*lambdaS * pow( sigma / T, 2);
    double d1 = lambdaS*sigma*sigma;
    double d2 = (lambdaS/2. + 3.*gX*gX/16.)*T*T;
    if (d1 < 0) {
        d1 = 0;
    }
    if (d2 < 0) {
        d2 = 0;
    }
    if(sigma>0)
    {
        double VHT = (9)*pow(T,4)/M_PI/M_PI/2.*Jbfun(xV)+ (3)*pow(T,4)/M_PI/M_PI/2.*Jbfun(xG)
        + (1)*pow(T,4)/M_PI/M_PI/2.*Jbfun(xS);
        double daisy = (3) * T / (12.*M_PI)*(pow(mV*mV,3/2.) - pow(mV*mV + PiV,3/2.));
        double VHT0 = (9+3+1)*pow(T,4)/M_PI/M_PI/2.*(-pow(M_PI,4)/45.)
        + (3) * T / (12.*M_PI)*(-pow(PiV,3/2.))
        + (4) * T / (12.*M_PI)*(-pow(3.*gX*gX/16.*T*T,3/2.));
        return VHT + daisy - VHT0 + dilPotential(sigma, gX, Ncft, f) - lambdaHS*pow(vhQCDref*sigma,2)/4.*QCDswitch;
    }
    else
    {
        return 0;
        
    }
    
}

//Dilaton potential at finite-temperarure in the scale-invariant U(1)_D model
double dilPotentialHT(double sigma, double T, double gX, double Ncft, double f)
{
    sigma = sigma + 1e-20;
    
    double vH = 246.*GeV;
    double lambdaH = 0.13;
    double lambdaS =  betaFun(sigma, gX, f);
    double QCDswitch = (1-tanh(20*(T-TcQCDref)))/2.;
    double mV = gX * sigma ;
    double PiV = 1./6.*pow(gX*T,2);
    double lambdaHS  = 2*lambdaH*pow(vH/f,2);
    double xV = pow(mV / T, 2);
    double xG = lambdaS * pow( sigma / T, 2);
    double xS = 3*lambdaS * pow( sigma / T, 2);
    double d1 = lambdaS*sigma*sigma;
    double d2 = (lambdaS/2. + 3.*gX*gX/16.)*T*T;
    if (d1 < 0) {
        d1 = 0;
    }
    if (d2 < 0) {
        d2 = 0;
    }
    if(sigma>0)
    {
        double VHT = (1)*pow(T,4)/M_PI/M_PI/2.*Jbfun(xV);
        double daisy = (1) * T / (12.*M_PI)*(pow(mV*mV,3/2.) - pow(mV*mV + PiV,3/2.));
        double VHT0 = (1)*pow(T,4)/M_PI/M_PI/2.*(-pow(M_PI,4)/45.) + (1) * T / (12.*M_PI)*(-pow(PiV,3/2.));
        return VHT + daisy - VHT0 + dilPotential(sigma, gX, Ncft, f) - lambdaHS*pow(vhQCDref*sigma,2)/4.*QCDswitch;
    }
    else
    {
        return 0;
        
    }
    
}

//QCD catalysis temperature in the scale-invariant SU(2)_D model
double TQCDfunSU2(double gX, double f)
{
    double vH = 246.*GeV;
    double lambdaH = 0.13;
    double lambdaHS  = 2*lambdaH*pow(vH/f,2);
    return 1.6*sqrt(lambdaHS)/gX*vhQCDref;
}

//QCD catalysis temperature in the scale-invariant U(1)_D model
double TQCDfun(double gX, double f)
{
    double vH = 246.*GeV;
    double lambdaH = 0.13;
    double lambdaHS  = 2*lambdaH*pow(vH/f,2);
    return sqrt(6*lambdaHS)/gX*vhQCDref;
}


//Find minimal field value of exit tunneling
double zeroIntercept(double T, double gX, double Ncft, double f)
{
    
    double min=-15.;
    double max=0.;
    double val=0;
    double prec=1.;
    for(int k=0; k<=1000 && prec>1e-10; k++)         //Newton-Raphson method
    {
        val = dilPotentialHT(pow(10, (min+max)/2.)*f, T, gX, Ncft, f);
        if(val < 0)
        {
            max = (min+max)/2.;
        }
        else
        {
            min = (min+max)/2.;
        }
        prec = max - min;
    }
    
    return pow(10,(min+max)/2.);
}

//Dilaton decay width
double sigmaDecayWidth(double gX, double f)
{
    return 3*pow(gX, 3)/(32*M_PI*pow(f, 2));
}


//Ordinary differential equations of motion of the scalar field inside the light-cone (bubble wall field trajectory during bubble propagation)
int isDecay = 1;
int isMin = 0;
void dilDiffTimeLike(double* q, double t, double* qp)
{
    double secondDerivative = ((dilPotential(q[0]*f*(1 + 1e-5 + 1e-4), gX, Ncft, f) - dilPotential(q[0]*f, gX, Ncft, f))/(q[0]*f*(1e-5 + 1e-4)) - (dilPotential(q[0]*f*(1 + 1e-5), gX, Ncft, f) - dilPotential(q[0]*f, gX, Ncft, f))/(q[0]*f*(1e-5)))/(q[0]*f*(1e-5 + 1e-4));
    if (secondDerivative < 0)
    {
        isMin = 0;
        secondDerivative = gX;
    }
    else
    {
        isMin = 1;
        secondDerivative = sqrt(secondDerivative);
    }

    qp[0]=q[1];
    qp[1]= -3/t*q[1] - isMin * isDecay * sigmaDecayWidth(secondDerivative, f) * q[1] / (T) - (dilPotentialHT(q[0]*f*(1 + 1e-5),T, gX, Ncft, f) - dilPotentialHT(q[0]*f,T, gX, Ncft, f))/(f*T*T*(q[0]*f*(1e-5)));

}

//3 for O3- and 4 for O4-symmetric bounce
double O3orO4;
//Ordinary differential equations of motion of the scalar field outside the light-cone (bounce trajectory during tunneling)
void dilDiffSpaceLike(double* q, double t, double* qp)
{
    qp[0]=q[1];
    qp[1]=-(O3orO4-1)/t*q[1] + (dilPotentialHT(q[0]*f*(1 + 1e-5 ), T, gX, Ncft, f) - dilPotentialHT(q[0]*f, T, gX, Ncft, f))/(f*T*T*(q[0]*f*(1e-5)));
    
}


//Find the radius of the bubble at nucleation
//This is the radius once the minimal value of the scalar trajecory is reached
double getMinScalar(std::map<double ,double[2]> scalarFun)
{
    std::map<double, double[2]>::iterator ptr;
    double min=1;
    double rmax;
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end(); ptr++)
    {
        if(ptr->second[0]<min)
        {
            min = ptr->second[0];
            rmax = ptr->first;
        }
    }
    return rmax;
}

//O4-symmetric bounce action
double actionO4(std::map<double ,double[2]> scalarFun, double T, double gX, double Ncft, double f)
{
    std::map<double, double[2]>::iterator ptr;
    double area;
    double rmax = getMinScalar(scalarFun);
    double sigmaMax;
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        sigmaMax = ptr-> second[0];
    }
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        area +=  (next(ptr,1)-> first - ptr-> first) * 2*M_PI*M_PI*pow(ptr-> first,3)*( pow(f/T,2) * 0.5 * pow(ptr-> second[1], 2) +  pow(T,-4) * (dilPotentialHT( ptr-> second[0] * f, T, gX, Ncft, f) - dilPotentialHT( sigmaMax * f, T, gX, Ncft, f) )  );
    }
    return area;
}


//O4-symmetric gradient energy
std::map<double, double> GradEnergyO4(std::map<double ,double[2]> scalarFun, double T, double gX, double Ncft, double f)
{
    std::map<double, double[2]>::iterator ptr;
    std::map<double, double> kinFun;
    kinFun.clear();
    double kin=0;
    double rmax = getMinScalar(scalarFun);
    
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        kinFun[ptr-> first] = kin + (next(ptr,1)-> first - ptr-> first) * 2*M_PI*M_PI*pow(ptr-> first,3)*( pow(f/T,2) * 0.5 * pow(ptr-> second[1], 2) );
        kin =  (next(ptr,1)-> first - ptr-> first) * 2*M_PI*M_PI*pow(ptr-> first,3)*( pow(f/T,2) * 0.5 * pow(ptr-> second[1], 2) );
    }
    return kinFun;
}

//Check if the validity of the O4-symmetric bounce action after it has been calculated
//See Section 6.1.2 in textbook [ArXiv] 2207.01633
double checkActionO4(std::map<double ,double[2]> scalarFun, double T, double gX, double Ncft, double f)
{
    std::map<double, double[2]>::iterator ptr;
    double kin;
    double pot;
    double rmax = getMinScalar(scalarFun);
    double sigmaMax;
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        sigmaMax = ptr-> second[0];
    }
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
        kin +=  (next(ptr,1)-> first - ptr-> first) * 2*M_PI*M_PI*pow(ptr-> first,3)*(pow(f/T,2) * 0.5 * pow(ptr-> second[1], 2) );
    for (ptr = scalarFun.begin(); ptr != scalarFun.end() && ptr-> first < rmax; ptr++)
        pot +=  (next(ptr,1)-> first - ptr-> first) * 2*M_PI*M_PI*pow(ptr-> first,3)*( pow(T,-4) * (dilPotentialHT( ptr-> second[0] * f, T, gX, Ncft, f) - dilPotentialHT( sigmaMax * f, T, gX, Ncft, f) ) );
    return (kin + 2*pot)/(kin + pot);
}

//O3-symmetric bounce action
double actionO3(std::map<double ,double[2]> scalarFun, double T, double gX, double Ncft, double f)
{
    std::map<double, double[2]>::iterator ptr;
    double area;
    double rmax = getMinScalar(scalarFun);
    double sigmaMax;
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        sigmaMax = ptr-> second[0];
    }
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        //        area +=  (next(ptr,1)-> first - ptr-> first) * 2*M_PI*M_PI*pow(ptr-> first,3)*(pow(f/T,2) * 0.5 * pow(ptr-> second[1], 2) + pow(T,-4) * dilPotentialHT( (ptr-> second[0]) * f, T, gX, Ncft, f) );
        area +=  (next(ptr,1)-> first - ptr-> first) * 4*M_PI*pow(ptr-> first,2)*( pow(f/T,2) * 0.5 * pow(ptr-> second[1], 2) +  pow(T,-4) * (dilPotentialHT( ptr-> second[0] * f, T, gX, Ncft, f) - dilPotentialHT( sigmaMax * f, T, gX, Ncft, f) )  );
    }
    return area;
}

//Check if the validity of the O3-symmetric bounce action after it has been calculated
//See Section 6.1.2 in textbook [ArXiv] 2207.01633
double checkActionO3(std::map<double ,double[2]> scalarFun, double T, double gX, double Ncft, double f)
{
    std::map<double, double[2]>::iterator ptr;
    double kin;
    double pot;
    double rmax = getMinScalar(scalarFun);
    double sigmaMax;
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
    {
        sigmaMax = ptr-> second[0];
    }
    for (ptr = scalarFun.begin(); next(ptr,1) != scalarFun.end() && ptr-> first < rmax; ptr++)
        kin +=  (next(ptr,1)-> first - ptr-> first) * 4*M_PI*pow(ptr-> first,2)*(pow(f/T,2)  * 0.5 * pow(ptr-> second[1], 2) );
    for (ptr = scalarFun.begin(); ptr != scalarFun.end() && ptr-> first < rmax; ptr++)
        pot +=  (next(ptr,1)-> first - ptr-> first) * 4*M_PI*pow(ptr-> first,2)*( pow(T,-4)  * (dilPotentialHT( ptr-> second[0] * f, T, gX, Ncft, f) - dilPotentialHT( sigmaMax * f, T, gX, Ncft, f) ) );
    return (kin + 3*pot)/(kin + pot);
}

//Analytical field value at the tunneling exit in the scale-invariant SU(2)_D model
//Sec 6.4.1 in textbook [ArXiv] 2207.01633
double chiExitAnalyticalSU2(double gX, double Ncft, double f, double Tnuc)
{
    double beta = betaFun(Tnuc, gX, f);
    double meff= sqrt(3./16.)*gX*Tnuc;
    double lambda = beta * log(exp(1./4)*f/Tnuc);
    double chiExit = 2*meff/sqrt(lambda);
    beta = betaFun(chiExit, gX, f);
    lambda = beta * log(exp(1./4)*f/chiExit);
    chiExit = 2*meff/sqrt(lambda);
    return chiExit;
}

//Analytical field value at the tunneling exit in the scale-invariant U(1)_D model
//ArXiv: 2311.13640
double chiExitAnalytical(double gX, double Ncft, double f, double Tnuc)
{
    double beta = betaFun(Tnuc, gX, f);
    double meff= sqrt(1./12.)*gX*Tnuc;
    double lambda = beta * log(exp(1./4)*f/Tnuc);
    double chiExit = 2*meff/sqrt(lambda);
    beta = betaFun(chiExit, gX, f);
    lambda = beta * log(exp(1./4)*f/chiExit);
    chiExit = 2*meff/sqrt(lambda);
    return chiExit;
}

//Critical bounce action below which nucleation becomes efficient (temporary function)
double findCriticalActionTmp(double T, double gX, double Ncft, double f)
{
    double Tc=5*gX*f/10.85/pow(100/108.75,0.25);
    double chiExit = chiExitAnalytical(gX, Ncft, f, T);
    double beta = betaFun(chiExit, gX, f);
    double rho = beta*pow(f/2.,4) + M_PI*M_PI/30.*(100.)*pow(T, 4);
    double hubble=sqrt(rho/3./mplbar/mplbar);
    return 4*log(Tc/hubble)+2*log(100./2./M_PI);
}

//Critical O4-symmetric bounce action below which nucleation becomes efficient
double findCriticalActionO4(double T, double gX, double Ncft, double f)
{
    double Tc=5*gX*f/10.85/pow(100/108.75,0.25);
    double chiExit = chiExitAnalytical(gX, Ncft, f,T);
    double beta = betaFun(chiExit, gX, f);
    double rho = beta*pow(f/2.,4) + M_PI*M_PI/30.*(100.)*pow(T, 4);
    double hubble=sqrt(rho/3./mplbar/mplbar);
    double action = findCriticalActionTmp(T, gX, Ncft, f);
    return 4*log(T/hubble)+2*log(action/2./M_PI);
}

//Critical O3-symmetric bounce action below which nucleation becomes efficient
double findCriticalActionO3(double T, double gX, double Ncft, double f)
{
    double Tc=5*gX*f/10.85/pow(100/108.75,0.25);
    double chiExit = chiExitAnalytical(gX, Ncft, f, T);
    double beta = betaFun(chiExit, gX, f);
    double rho = beta*pow(f/2.,4) + M_PI*M_PI/30.*(100.)*pow(T, 4);
    double hubble=sqrt(rho/3./mplbar/mplbar);
    double action = findCriticalActionTmp(T, gX, Ncft, f);
    return 4*log(T/hubble)+3./2.*log(action/2./M_PI);;
}

//Analytical nucleation temperature in the scale-invariant SU(2)_D model
//Sec 6.4.1 in textbook [ArXiv] 2207.01633
double TnucAnalyticalSU2(double gX, double Ncft, double f)
{
    double Tc=5*0.19*gX*f/pow(100,0.25);
    double chiExit = chiExitAnalytical(gX, Ncft, f, Tc/100.);
    double beta = betaFun(chiExit, gX, f);
    double rho = beta*pow(f/2.,4) + M_PI*M_PI/30.*(100.)*pow(T, 4);
    double Ttmp = Tc/100.;
    double lambdaHS = 2*0.13*pow(246./f,2);
    double TQCD  = 1.6*sqrt(lambdaHS)/gX*vhQCDref;
    for(int i=0; i<1; i++)
    {
        Ttmp = Tc*exp(-1018/findCriticalActionO3(Ttmp, gX, Ncft, f) / pow(gX, 3));
    }
    //return max(Ttmp,TQCD);
    return Ttmp;
}

//Analytical nucleation temperature in the scale-invariant U(1)_D model
//ArXiv: 2311.13640
double TnucAnalytical(double gX, double Ncft, double f)
{
    double Tc=5*gX*f/10.85/pow(100/108.75,0.25);
    double chiExit = chiExitAnalytical(gX, Ncft, f, Tc/1000.);
    double beta = betaFun(chiExit, gX, f);
    double rho = beta*pow(f/2.,4) + M_PI*M_PI/30.*(100.)*pow(T, 4);
    double Ttmp = Tc/100.;
    double lambdaHS = 2*0.13*pow(246./f,2);
    double TQCD  = sqrt(6*lambdaHS)/gX*vhQCDref;
    for(int i=0; i<10; i++)
    {
        Ttmp = Tc*exp(-127.3/findCriticalActionO3(Ttmp, gX, Ncft, f) / pow(gX, 3));
        if(Ttmp< 1E-100|| Ttmp > Tc)
        {
            cout<<"Stuck in false Vacuum ! (Increase gX)"<<endl;
        }
    }
    
    return Ttmp;
}

//Calculate the phase transition completion rate
double betaOHfun(double action, double action2, double T, double T2, double gX, double f, double sigma0)
{
    double delta = deltabeta;
    double der = (action2-action)/(T2 - T);
    double beta = T*der - 4 - 1.5*T*der/action;
    //cout<<"action= "<<action<<" action2 = "<<action2<<" der = "<<der<<" beta = "<<beta<<endl;
    return beta;
}

//Calculate the bounce trajectory for a given temperature
double dilEOMSpaceLike1TrajT(double sigma0, double T, double duration, double maxSteps, fstream& dat)
{
    cout<<"gX "<<gX<<" Ncft "<<Ncft<<" f "<<f<<endl;
    int n, i,j;
    
    //We fix the step size, the anticipated error
    double dt= 0.0001;
    double desiredErrorMin = 1e-8;
    double desiredErrorMax = 1e-6;
    double error;
    int N=int(duration/dt);
    time_t start=time(NULL);
    
    //We create the array in which we will store the trajectory
    double* q=new double[2];
    std::map<double ,double[2]> scalarFun;
    
    //We create the time variable
    double t;
    
    //We create variables which store the bounce action
    double action;
    double actionCrit;
    double checkAction;
    
    //We fix the initial field value in the reversed potential which corresponds to the exit tunneling point in the normal potential (it corresponds to the first try of the overshoot/undershoot method)
    q[0]=(10*zeroIntercept(T, gX, Ncft, f)+1)/5.;
    q[1]=0.;
    t=1e-50;
    //We define the variable which take the initial field value (tunneling exit)
    double val;
    //We fix an estimate of the tunneling temperature
    double TnucA = TnucAnalytical(gX, Ncft, f);
    //Check that the error parameter (d-2)K+dU is close to zero
    double DeltaSigma = 1e-4;
    double itErrorNbr = 0;
    double derivativeEnd;
    //loop for solving the nucleation temperature
    do  {
        //Check that the duration is long enough such that the field has the time to reach 0
        do {
            double zeroInter = zeroIntercept(T, gX, Ncft, f);
            cout<<"zeroInter = "<<zeroInter<<" duration = "<<duration<<endl;
            double min=zeroInter;
            double max=1;
            int itNbr = 0;
            derivativeEnd = abs(q[1]/q[0]);
            //loop for solving the initial field value (the euclidean trajectory with correct boundary condition)
            do {
                val = (min+max)/2.;
                val = sqrt(min*max);
                q[0]=val;
                q[1]=0.;
                t=1e-50;
                scalarFun.clear();
                for(n=1; t < duration && q[0] > -zeroInter; n++)
                {
                    rk5(dilDiffSpaceLike,q,t,dt,error,desiredErrorMin,desiredErrorMax);
                    scalarFun[t][0]=q[0];
                    scalarFun[t][1]=q[1];
                }
                if(q[0]<0)
                {
                    max = val;
                }
                if(q[0]>0)
                {
                    min = val;
                }
                itNbr++;
            
            //We check if the search on the initial field value (tunneling exit) has converged
            } while(abs(max-min)/(max+min)> DeltaSigma && itNbr<1e6);
            cout<<"  "<<endl;
            
        //Increase the duration in case it is not long enough
            duration = duration * 2;
            cout<<"while "<<abs(q[1]/q[0])<<" "<<abs(derivativeEnd - abs(q[1]/q[0]))/abs(q[1]/q[0])<<" error = "<< error<<endl;
        
        //We check if the duration of trajectory is long enough for the field to reach the origin
        } while( abs(q[1]/q[0]) > 0.1 && abs(derivativeEnd - abs(q[1]/q[0]))/(derivativeEnd+abs(q[1]/q[0])) > 0.01);
        
        //Assigne the exit tunneling point to the value we just calculated, and calculate the bounce action
        sigma0 = val;
        if(O3orO4 == 3)
        {
            action = actionO3(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO3(T, gX, Ncft, f);
            checkAction = checkActionO3(scalarFun, T, gX, Ncft, f);
        }
        else
        {
            action = actionO4(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO4(T, gX, Ncft, f);
            checkAction = checkActionO4(scalarFun, T, gX, Ncft, f);
        }
        itErrorNbr++;
        cout<<"DeltaSigma = "<<DeltaSigma<<endl;
        DeltaSigma = DeltaSigma/2.;
    
    //We check the error on the bounce action is small enough
    } while( checkAction > 0.1 && itErrorNbr < 1e2);
    
    //Save the field trajectory in a list in case we want to use it later (e.g. to plot it or calculate the bubble radius)
        for (std::map<double, double[2]>::iterator ptr = scalarFun.begin(); ptr != scalarFun.end(); ptr++)
        {
            double rmax = getMinScalar(scalarFun);
            if(ptr->first<rmax)
            {
                dat << ptr->first << " "<<ptr-> second[0]<<endl;
            }
    }
    dat<<"  "<<endl;
    dat<<"  "<<endl;
    //Print the results for the bounce action, the critical action, the error on the action, the initial field value and the bubble size at nucleation
    cout<<setprecision(3)<<" action "<<action<<" actionCrit "<<actionCrit<<" checkAction "<<checkAction<<" phiIni "<<val<<" Bubble size "<<getMinScalar(scalarFun)<<endl;
    //Print vector mass, the dark Higgs vev
    cout<<" gX f "<<gX*f<<" f "<<f<<endl;
    //Print symmetry of the bounce, nucleation temperature numerically calculated and analytical estimate
    cout<<"O"<<O3orO4<<" T/f "<<T/f<<" Tana/f "<<TnucA/f<<" T "<<T<<endl;
    cout<<"Tnuc analytical  "<<TnucAnalytical(gX, Ncft, f)<<" Tnuc numerical  "<<T<<endl;

//return the inital field value
    return val;
}

//Find the bounce trajectory at the nucleation temperature, deduce the nucleation temperature Tnuc and the PT completion rate beta/H
double dilEOMSpaceLike1Traj(double sigma0, double duration0, double maxSteps, fstream& dat)
{
    cout<<" gX "<<gX<<" f "<<f<<" mV "<<gX*f<<endl;
    int n, i,j;
    
    //We fix the step size and the anticipated error
    double dt= 0.0001;
    double desiredErrorMin = 1e-9;
    double desiredErrorMax = 1e-7;
    double error;
    //Create variable for the duration of the trajectory
    double duration;
    int N=int(duration/dt);
    time_t start=time(NULL);
    
    //We create the time variable
    double t;
    
    //Create array which stores the trajectory (to pass to rk5 routine)
    double* q=new double[2];
    //Create array which stores the trajectory (to store forever)
    std::map<double ,double[2]> scalarFun;
    
    //We create variable which stores the bounce action
    double action;
    double actionCrit;
    double checkAction;
    //Bounce action for T+\deltaT
    double action2;
    //The **LEGENDARY** \beta/H
    double betaOH;

    //We define the variable which take the initial field value (tunneling exit)
    double val;
    //We fix an estimate of the tunneling temperature and its prior interval of expected values over which the search will be carried
    double TnucA = TnucAnalytical(gX, Ncft, f);
    double minT = max(TnucA/(100), 1.01*TQCDfun(gX, f));
    double maxT = TnucA*100.;
    //loop for solving the nucleation temperature
    do {
        T = sqrt(minT*maxT);
        //Variable to check that the error parameter (d-2)K+dU is close to zero
        double itErrorNbr = 0;
        double DeltaSigma = 1e-4;
        //Variable to check that the derivative at the end of the trajectory is close to zero
        double derivativeEnd;
        //Check that the error parameter (d-2)K+dU is close to zero
        do  {
            duration = duration0;
            //Check that the duration is long enough such that the field has the time to reach 0
            do {
                double zeroInter = zeroIntercept(T, gX, Ncft, f);
                double min=zeroInter;
                double max=1;
                int itNbr = 0;
                derivativeEnd = abs(q[1]/q[0]);
                //loop for solving the initial field value (the euclidean trajectory with correct boundary condition)
                do {
                    val = (min+max)/2.;
                    val = sqrt(min*max);
                    q[0]=val;
                    q[1]=0.;
                    t=1e-50;
                    scalarFun.clear();
                    for(n=1; t < duration && q[0] > -zeroInter; n++)
                    {
                        rk5(dilDiffSpaceLike,q,t,dt,error,desiredErrorMin,desiredErrorMax);
                        scalarFun[t][0]=q[0];
                        scalarFun[t][1]=q[1];
                    }
                    if(q[0]<0)
                    {
                        max = val;
                    }
                    if(q[0]>0)
                    {
                        min = val;
                    }
                    itNbr++;
                    
                //We check if the search on the initial field value (tunneling exit) has converged
                } while(abs(max-min)/(max+min)>DeltaSigma && itNbr<1e6);
                //Increase the duration in case it is not long enough
                duration = duration * 1.1;
            
        //We check if the duration of trajectory is long enough for the field to reach the origin
        } while( abs(q[1]/q[0]) > 0.1 && abs(derivativeEnd - abs(q[1]/q[0]))/abs(q[1]/q[0]) > 0.01);

        //Assigne the exit tunneling point to the value we just calculated, and calculate the bounce action
        sigma0 = val;
        if(O3orO4 == 3)
        {
            action = actionO3(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO3(T, gX, Ncft, f);
            checkAction = checkActionO3(scalarFun, T, gX, Ncft, f);
        }
        else
        {
            action = actionO4(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO4(T, gX, Ncft, f);
            checkAction = checkActionO4(scalarFun, T, gX, Ncft, f);
        }
            itErrorNbr++;
            DeltaSigma = DeltaSigma/2.;
            
        //We check the error on the bounce action is small enough
        } while( checkAction > 0.01 && itErrorNbr < 1e2);
        
        if(action>actionCrit)
        {
            maxT=T;
        }
        else
        {
            minT=T;
        }
        
    //Print the results for the bounce action, the critical action, the error on the action, nucleation temperature, its error interval, the initial field value and the bubble size at nucleation
        cout<<setprecision(3)<<" action "<<action<<" actionCrit "<<actionCrit<<" checkAction "<<checkAction<<" T/f "<<T/f<<" maxT/f "<<maxT/f<<" minT/f "<<minT/f<<" phiIni "<<val<<" Bubble size "<<getMinScalar(scalarFun)<<endl;
    //We check the error on nucleation temperature is small enough
    } while((maxT - minT)/maxT > 0.1);
    
    //Compute the last trajectory more precisely (to compute beta/H)
    cout<<"   "<<endl;
    cout<<" Compute beta/H "<<endl;
    cout<<" Compute the last trajectory more precisely "<<endl;
    q[0]=val;
    q[1]=0.;
    t=1e-50;
    dt= 0.00001;
    desiredErrorMin = 1e-9;
    desiredErrorMax = 1e-7;
    double derivativeEnd;
    //Check that the error parameter (d-2)K+dU is close to zero
    double itErrorNbr = 0;
    double DeltaSigma = 1e-5;
    do  {
        duration = duration0;
        //Check that the duration is long enough such that the field has the time to reach 0
        do {
            double zeroInter = zeroIntercept(T, gX, Ncft, f);
            double min=zeroInter;
            double max=1;
            int itNbr = 0;
            derivativeEnd = abs(q[1]/q[0]);
            //loop for solving the initial field value (the euclidean trajectory with correct BC)
            do {
                val = (min+max)/2.;
                val = sqrt(min*max);
                q[0]=val;
                q[1]=0.;
                t=1e-50;
                scalarFun.clear();
                for(n=1; t < duration && q[0] > -zeroInter; n++)
                {
                    rk5(dilDiffSpaceLike,q,t,dt,error,desiredErrorMin,desiredErrorMax);
                    scalarFun[t][0]=q[0];
                    scalarFun[t][1]=q[1];
                }
                if(q[0]<0)
                {
                    max = val;
                }
                if(q[0]>0)
                {
                    min = val;
                }
                itNbr++;
                
            } while(abs(max-min)/(max+min)>DeltaSigma && itNbr<1e6);
            
            duration = duration * 1.1;
        } while( abs(q[1]/q[0]) > 0.01 && abs(derivativeEnd - abs(q[1]/q[0]))/abs(q[1]/q[0]) > 0.01);
        
        sigma0 = val;
        if(O3orO4 == 3)
        {
            action = actionO3(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO3(T, gX, Ncft, f);
            checkAction = checkActionO3(scalarFun, T, gX, Ncft, f);
        }
        else
        {
            action = actionO4(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO4(T, gX, Ncft, f);
            checkAction = checkActionO4(scalarFun, T, gX, Ncft, f);
        }
            itErrorNbr++;
            DeltaSigma = DeltaSigma/2.;
            
        } while( checkAction > 0.01 && itErrorNbr < 1e2);
        
        
    //Print the results for the bounce action, the critical action, the error on the action, nucleation temperature, its error interval, the initial field value and the bubble size at nucleation
    cout<<setprecision(3)<<" action "<<action<<" actionCrit "<<actionCrit<<" checkAction "<<checkAction<<" T/f "<<T/f<<" maxT/f "<<maxT/f<<" minT/f "<<minT/f<<" phiIni "<<val<<" Bubble size "<<getMinScalar(scalarFun)<<endl;
    
    
    //Compute the trajectory for a sligtly hotter universe (to compute beta/H)
    cout<<"   "<<endl;
    cout<<" Prepare slighlty hotter universe "<<endl;
    q[0]=val;
    q[1]=0.;
    t=1e-50;
    T=T*(1+deltabeta);
    derivativeEnd;
    //Check that the error parameter (d-2)K+dU is close to zero
    itErrorNbr = 0;
    DeltaSigma = 1e-5;
    do  {
        duration = duration0;
        //Check that the duration is long enough such that the field has the time to reach 0
        do {
            double zeroInter = zeroIntercept(T, gX, Ncft, f);
            double min=zeroInter;
            double max=1;
            int itNbr = 0;
            derivativeEnd = abs(q[1]/q[0]);
            //loop for solving the initial field value (the euclidean trajectory with correct BC)
            do {
                val = (min+max)/2.;
                val = sqrt(min*max);
                q[0]=val;
                q[1]=0.;
                t=1e-50;
                scalarFun.clear();
                for(n=1; t < duration && q[0] > -zeroInter; n++)
                {
                    rk5(dilDiffSpaceLike,q,t,dt,error,desiredErrorMin,desiredErrorMax);
                    scalarFun[t][0]=q[0];
                    scalarFun[t][1]=q[1];
                }
                if(q[0]<0)
                {
                    max = val;
                }
                if(q[0]>0)
                {
                    min = val;
                }
                itNbr++;
                
            } while(abs(max-min)/(max+min)>DeltaSigma && itNbr<1e6);
            
            duration = duration * 1.1;
            //cout<<"duration "<<duration<<endl;
        } while( abs(q[1]/q[0]) > 0.01 && abs(derivativeEnd - abs(q[1]/q[0]))/abs(q[1]/q[0]) > 0.01);
        
        sigma0 = val;
        if(O3orO4 == 3)
        {
            action2 = actionO3(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO3(T, gX, Ncft, f);
            checkAction = checkActionO3(scalarFun, T, gX, Ncft, f);
        }
        else
        {
            action2 = actionO4(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO4(T, gX, Ncft, f);
            checkAction = checkActionO4(scalarFun, T, gX, Ncft, f);
        }
            itErrorNbr++;
            DeltaSigma = DeltaSigma/2.;
            
        } while( checkAction > 0.01 && itErrorNbr < 1e2);
        
    //Print the results for the bounce action, the critical action, the error on the action, nucleation temperature, its error interval, the initial field value and the bubble size at nucleation
    cout<<setprecision(3)<<" action "<<action<<" actionCrit "<<actionCrit<<" checkAction "<<checkAction<<" T/f "<<T/f<<" maxT/f "<<maxT/f<<" minT/f "<<minT/f<<" phiIni "<<val<<" Bubble size "<<getMinScalar(scalarFun)<<endl;
    
    
    //compute beta/H
    betaOH = betaOHfun(action,action2, T/(1+deltabeta), T, gX, f, val);
    dat<<" "<<endl;
    dat<<" "<<endl;
        
    //Print gauge coupling, the dark Higgs vev, vector mass, the PT completion rate
    cout<<"   "<<endl;
    cout<<" gX "<<gX<<" f "<<f<<" mV "<<gX*f<<" beta/H "<<betaOH<<endl;
    cout<<"   "<<endl;
    //Print symmetry of the bounce, nucleation temperature numerically calculated and analytical estimate
    cout<<" O"<<O3orO4<<"-symmetric bounce with T/f "<<T/f<<" Tana/f "<<TnucA/f<<" T "<<T<<endl;
    cout<<" Tnuc analytical  "<<TnucAnalytical(gX, Ncft, f)<<" Tnuc numerical  "<<T<<endl;

//return the inital field value
    return val;
}



    
//Calculate the bounce action for different values of the temperature
void dilEOMSpaceLikeS(int nT, double duration0, double maxSteps, fstream& dat)
{
    cout<<" gX "<<gX<<" f "<<f<<" mV "<<mV<<endl;
    int n, i,j;
    
    //We fix the step size and the anticipated error
    double dt= 0.0001;
    double desiredErrorMin = 1e-8;
    double desiredErrorMax = 1e-6;
    //Create variable for the duration of the trajectory
    double duration;
    int N=int(duration/dt);
    time_t start=time(NULL);
    
    //We create the time variable
    double t;
    double error;
    
    //Create array which stores the trajectory (to pass to rk5 routine)
    double* q=new double[2];
    //Create array which stores the trajectory (to store forever)
    std::map<double ,double[2]> scalarFun;
    
    //We create variable which stores the bounce action
    double action;
    double actionCrit;
    double checkAction;
    //Bounce action for T+\deltaT
    double action2;
    //The **LEGENDARY** \beta/H
    double betaOH;
    
    //We define the variable which take the initial field value (tunneling exit)
    double val;
    //We fix an estimate of the tunneling temperature and its prior interval of expected values over which the search will be carried
    double TnucA = TnucAnalytical(gX, Ncft, f);
    double minT = max(5*TnucA, 1.05*TQCDfun(gX, f));
    minT=0.007;
    double maxT = minT*1000.;
    maxT=200.;
    //loop for solving the nucleation temperature
    for (int i=0; i <= nT; i++) {
        T = iterateLog(i, nT, minT, maxT);
        //Variable to check that the error parameter (d-2)K+dU is close to zero
        double itErrorNbr = 0;
        double DeltaSigma = 1e-4;
        //Variable to check that the derivative at the end of the trajectory is close to zero
        double derivativeEnd;
        //Check that the error parameter (d-2)K+dU is close to zero
        do  {
            duration = duration0;
            //Check that the duration is long enough such that the field has the time to reach 0
            do {
                double zeroInter = zeroIntercept(T, gX, Ncft, f);
                double min=zeroInter;
                double max=1;
                int itNbr = 0;
                derivativeEnd = abs(q[1]/q[0]);
                //loop for solving the initial field value (the euclidean trajectory with correct boundary condition)
                do {
                    val = (min+max)/2.;
                    val = sqrt(min*max);
                    q[0]=val;
                    q[1]=0.;
                    t=1e-50;
                    scalarFun.clear();
                    for(n=1; t < duration && q[0] > -zeroInter; n++)
                    {
                        rk5(dilDiffSpaceLike,q,t,dt,error,desiredErrorMin,desiredErrorMax);
                        scalarFun[t][0]=q[0];
                        scalarFun[t][1]=q[1];
                    }
                    if(q[0]<0)
                    {
                        max = val;
                    }
                    if(q[0]>0)
                    {
                        min = val;
                    }
                    itNbr++;
                //We check if the search on the initial field value (tunneling exit) has converged
                } while(abs(max-min)/(max+min)>DeltaSigma && itNbr<1e6);
                
                //Increase the duration in case it is not long enough
                duration = duration * 1.1;

            //We check if the duration of trajectory is long enough for the field to reach the origin
            } while( abs(q[1]/q[0]) > 0.1 && abs(derivativeEnd - abs(q[1]/q[0]))/abs(q[1]/q[0]) > 0.1);
            
        //Assigne the exit tunneling point to the value we just calculated, and calculate the bounce action
        sigma0 = val;
        if(O3orO4 == 3)
        {
            action = actionO3(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO3(T, gX, Ncft, f);
            checkAction = checkActionO3(scalarFun, T, gX, Ncft, f);
        }
        else
        {
            action = actionO4(scalarFun, T, gX, Ncft, f);
            actionCrit = findCriticalActionO4(T, gX, Ncft, f);
            checkAction = checkActionO4(scalarFun, T, gX, Ncft, f);
        }
            itErrorNbr++;
            DeltaSigma = DeltaSigma/2.;
            
        //We check the error on the bounce action is small enough
        } while( checkAction > 0.01 && itErrorNbr < 1e2);
        
        //Save the results in a file (dark Higgs vev, gauge coupling, nucleation temperature, bounce action, error action, critical action)
        dat<<f<<" "<<gX<<" "<<T/f<<" "<<action<<" "<<checkAction<<" "<<actionCrit<<endl;
        
        //Print the results for the maximal reheating temperature, gauge coupling bounce action, the critical action, the error on the action, nucleation temperature, its error interval, the initial field value and the bubble size at nucleation
        cout<<setprecision(3)<<"Teq "<<Teq<<" gX "<<gX<<" action "<<action<<" actionCrit "<<actionCrit<<" checkAction "<<checkAction<<" T/f "<<log10(T/f)<<" maxT/f "<<log10(maxT/f)<<" minT/f "<<log10(minT/f)<<" phiIni "<<val<<" Bubble size "<<getMinScalar(scalarFun)<<endl;
        
        //You stop the loop on temperature when S is much larger than S_crit
        if(nT>25)
        {
            if(action/actionCrit > 1.5)
            {
                i=nT;
            }
        }
        else
        {
            if(action/actionCrit > 2)
            {
                i=nT;
            }
        }
        
    }
}


//Plot trajectory
void PlotSpaceLike()
{
    fstream gnu("data.gnu", ios::out);
    gnu<<"set terminal epscairo font 'Gill Sans' rounded fontscale 0.6"<<endl;
    gnu<<"set output 'courbe.eps'"<<endl;
    gnu<<"set multiplot layout 1, 1"<<endl;
    gnu<<"set ylabel 'Φ/f'"<<endl;
    gnu<<"set border 6 back linestyle 100"<<endl;
    gnu<<"set grid linestyle 100"<<endl;
    gnu<<"set tmargin 4"<<endl;
    gnu<<"set bmargin 3.1"<<endl;;
    gnu<<"unset xrange"<<endl;
    gnu<<"set size 1, 1"<<endl;
    gnu<<"set xlabel 'Tn s   with  s^2 = r^2- t^2'"<<endl;
    gnu<<"set key bottom title 'Bounce trajectory'"<<endl;
    gnu<<std::setprecision(3)<<"set title 'f = "<<int(f/1000)<<" TeV,  g_X = "<< gX <<",  T/f = "<<T/f<<",  Φ_{0} = "<<sigma0<<"' offset 0, -0.5"<<endl;
    gnu<<"plot 'data.dat' i 0 u 1:2 w l lc 3 lw 7 title 'Φ/f'"<<endl;
    gnu<<"set key off"<<endl;
    gnu<<"unset key"<<endl;
    gnu<<"unset multiplot"<<endl;
    gnu.close();
    
    system("gnuplot data.gnu && sleep 2 && open -a preview courbe.eps");
}


//Plot scalar potential
void plotPotential()
{
    fstream dat("data.dat",ios::out);
    //We plot the scalar potential in case you are interested
    double u;
    int npt= 1000;
    for(int i=0; i<=npt; i++)
    {
        u=iterateLog(i,npt,1e-10,1.4)*f;
        //dat<<u<<" "<<dilPotential(u,gX,Ncft,f)<<endl;
        //dat<<u/f<<" "<<dilPotentialHT(u,T,gX,Ncft,f)/pow(gX*f/4.,2)<<endl;
        dat<<u/f<<" "<<dilPotentialHT(u,T,gX,Ncft,f)/pow(f,4)<<endl;
    }

    dat.close();
    
    fstream gnu("data.gnu", ios::out);
    gnu<<"set terminal epscairo font 'Gill Sans' rounded fontscale 0.6"<<endl;
    gnu<<"set output 'courbe.eps'"<<endl;
    gnu<<"set multiplot layout 2, 1"<<endl;
    gnu<<"set border 6 back linestyle 100"<<endl;
    gnu<<"set grid linestyle 100"<<endl;;
    gnu<<"set xlabel 'log_{10}(Φ/f)'"<<endl;
    gnu<<"set ylabel 'V(Φ)/f^4'"<<endl;
    gnu<<"set logscale x"<<endl;
    gnu<<"set xrange [1e-10:1e-5]"<<endl;
    gnu<<"set key bottom left box"<<endl;
    gnu<<"plot 'data.dat' i 0 u 1:2 w l lc 3 lw 7  title 'potential (zoom)'"<<endl;
    gnu<<"unset logscale x"<<endl;
    gnu<<"set xrange [1e-10:1.4]"<<endl;
    gnu<<"set xlabel 'Φ/f'"<<endl;
    gnu<<"set key top left box"<<endl;
    gnu<<"plot 'data.dat' i 0 u 1:2 w l  lc 3 lw 7  title 'potential'"<<endl;
    gnu<<"unset key"<<endl;
    gnu<<"unset multiplot"<<endl;
    gnu.close();
    
    system("gnuplot data.gnu  && sleep 4 && open -a preview courbe.eps");
}

//Plot S3/T tables
void PlotGeneral()
{
    fstream gnu("data.gnu", ios::out);
    gnu<<"set terminal epscairo font 'Gill Sans' rounded fontscale 0.6"<<endl;
    gnu<<"set output 'courbe.eps'"<<endl;
    gnu<<"set multiplot layout 1, 1"<<endl;
    gnu<<"set border 6 back linestyle 100"<<endl;
    gnu<<"set grid linestyle 100"<<endl;;
    gnu<<"set title 'Bounce action'"<<endl;
    gnu<<"set size 1, 1"<<endl;
    gnu<<"set xlabel 'T_{nuc}/f'"<<endl;
    gnu<<"set ylabel 'S3/T'"<<endl;
    gnu<<"set logscale x"<<endl;
    gnu<<"set format x '10^{%T}'"<<endl;
    gnu<<"set key top left box"<<endl;
    gnu<<"plot 'data.dat' i 0 u 3:4  w l  lc 3 lw 7   title 'bounce action'";
    
    for (int igX=1; igX<=ngX; igX++) {
        gnu<<", 'data.dat' i "<<igX<<" u 3:4 w l  title 'g_X = "<<iterateLin(igX, ngX, 0.46, 0.6)<<"'";
    }
    gnu<<""<<endl;
    gnu<<"unset key"<<endl;
    gnu<<"unset multiplot"<<endl;
    gnu.close();
    
    system("sleep 2 && gnuplot data.gnu  && sleep 2 && open -a preview courbe.eps");
}


static inline void loadbar(time_t start)
{
    double Z0 =difftime(time(NULL),start);
    int Z1=(int)(Z0/3600);
    int Z2=(int)(Z0/60)-Z1*60;
    int Z3=(int)(Z0)-Z1*3600-Z2*60;
//    cout<<"\r"<<"Elapsed time= "<<Z1<<" h "<<Z2<<" min "<<Z3<<" s "<<endl;
    cout<<"Elapsed time = "<<Z1<<" h "<<Z2<<" min "<<Z3<<" s "<<endl;
}

int main()
{
    time_t start=time(NULL);
    
    //Dummy variables (their use is discontinued and come from a old code)
    Ncft = 10;
    isDecay = 0;
    
    //O3 or O4 tunneling ??
    O3orO4=3; //3 for O3 and 4 for O4
    
    //QCD catalysis parameters
    //Higgs vev induced by QCD
    vhQCDref = 0.1*GeV;
    //QCD confinement temperature
    TcQCDref = 0.1*GeV;

    //step in temperature for computing finite difference of the action (for beta/H)
    deltabeta = 0.2;
    
    //Duration of the trajectory in unit of Tn^{-1}
    double duration = 20;
    //Max number of steps
    double maxSteps = 1e7;
    
//CHOOSE WHAT YOU WANT TO DO (change the #define MODE(1,2,3,4) below accordingy)
//MODE1: Plot the scalar potential
//MODE2: Calculate and plot the bouncing trajectory at a given temperature T
//MODE3: Compute tables of S3/T for list of values of gauge coupling constant g_X and temperature T
//MODE4: Calculate the nucleation temperature and the phase transition completion rate \beta/H.
#define MODE3
    

#ifdef MODE1
//MODE1: Plot the scalar potential
    gX=0.45;
    mV = 8*TeV;
    Teq = mV/10.85*pow(100/108.75,0.25);
    f = mV/gX;
    T=pow(10,-5.2) * f;
    plotPotential();
#endif
    
#ifdef MODE2
//MODE2: Calculate and plot the bouncing trajectory at a given temperature T
    fstream dat("data.dat",ios::out);
    gX=0.45;
    mV = 8*TeV;
    Teq = mV/10.85*pow(100/108.75,0.25);
    f = mV/gX;
    T=pow(10,-5.2) * f;
    sigma0=dilEOMSpaceLike1TrajT(0,T, duration, maxSteps, dat);
    PlotSpaceLike();
#endif
    
#ifdef MODE3
//MODE4: Compute tables of S3/T for list of values of gauge coupling constant g_X and temperature T
//Number of temperature points
    nT = 100;
//Number of gauge coupling points (ngX = 1 if one single value of gX)
    ngX = 25;
    ngX = 1;
    fstream dat("data.dat",ios::out);
    for (int igX=0; igX < ngX; igX++) {
        //for (int iOmega=0; iOmega <= nOmega; iOmega++) {
        gX = iterateLin(igX, ngX, 0.46, 0.6);
        gX=0.45;
        mV = 8*TeV;
        Teq = mV/10.85*pow(100/108.75,0.25);
        f = mV/gX;
//Calculate the bounce action as a function of temperature
        dilEOMSpaceLikeS(nT, duration, maxSteps, dat);
        dat<<" "<<endl;
        dat<<" "<<endl;
    }
//Plot with GNUPLOT
    PlotGeneral();
#endif
    
//MODE4: Calculate the nucleation temperature and the phase transition completion rate \beta/H
#ifdef MODE4
    fstream dat("data.dat",ios::out);
    gX=0.45;
    mV = 8*TeV;
    Teq = mV/10.85*pow(100/108.75,0.25);
    f = mV/gX;
    sigma0=dilEOMSpaceLike1Traj(0, duration, maxSteps, dat);
#endif


    
    return 0;
}

