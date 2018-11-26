//
//  RoFSuP.h
//  RoFSuP
//
//  Created by MacBook on 8/28/17.
//  Copyright © 2017 Casto. All rights reserved.
//

#ifndef RoFSuP_h
#include <ctime>
#include <math.h>
#include <iostream>
#include <chrono>
#include <complex>
#include <vector>
using std::vector;

double f1(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)  //Rescaled LSE Symmetric Function
{//Rescaled LSE Symmetric Dispersion Equation
    double kx = M_PI*n/w;
    double kch = sqrt((1.0-B2)*(kmd/(c-b)*kmd/(c-b)+kx*kx)/(Ep*Mu*B2-1.0)+kx*kx);
    double kchdkmd= sqrt((1.0-B2)*(1/(c-b)*1/(c-b)+kx/kmd*kx/kmd)/(Ep*Mu*B2-1.0)+kx/kmd*kx/kmd);
    return exp(-b*kch)*(Ep*kchdkmd*sinh(kch*b)*cos(kmd)-1/(c-b)*sin(kmd)*cosh(kch*b));
}
double f2(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)
{//Rescaled LSE Asymmetric Dispersion Equation
    double kx = M_PI*n/w;
    double kch = sqrt((1-B2)*(kmd/(c-b)*kmd/(c-b)+kx*kx)/(Ep*Mu*B2-1)+kx*kx);
    double kchdkmd= sqrt((1.0-B2)*(1/(c-b)*1/(c-b)+kx/kmd*kx/kmd)/(Ep*Mu*B2-1.0)+kx/kmd*kx/kmd);
    return exp(-b*kch)*(Ep*kchdkmd*cosh(kch*b)*cos(kmd)-1/(c-b)*sin(kmd)*sinh(kch*b));
}

double f3(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)
{//Rescaled LSM Symmetric Dispersion Equation
    double kx = M_PI*n/w;
    double kch = sqrt((1-B2)*(kmd/(c-b)*kmd/(c-b)+kx*kx)/(Ep*Mu*B2-1)+kx*kx);
    double kchdkmd= sqrt((1.0-B2)*(1/(c-b)*1/(c-b)+kx/kmd*kx/kmd)/(Ep*Mu*B2-1.0)+kx/kmd*kx/kmd);

    return exp(-b*kch)*(Mu*kchdkmd*sinh(kch*b)*sin(kmd)+1/(c-b)*cos(kmd)*cosh(kch*b));
}
double f4(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)
{//Rescaled LSM Asymmetric Dispersion Equation
    double kx = M_PI*n/w;
    double kch = sqrt((1-B2)*(kmd/(c-b)*kmd/(c-b)+kx*kx)/(Ep*Mu*B2-1)+kx*kx);
    double kchdkmd= sqrt((1.0-B2)*(1/(c-b)*1/(c-b)+kx/kmd*kx/kmd)/(Ep*Mu*B2-1.0)+kx/kmd*kx/kmd);
    return exp(-b*kch)*(1/Mu/kchdkmd/(c-b)*sinh(kch*b)*cos(kmd)+sin(kmd)*cosh(kch*b));
}

double YE_S(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//LSE Symmetric Eigenfunction (y-component)
    double YEs = 0.0;
    if(y >= -c && y < -b)
    {
        YEs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YEs = 1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YEs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YEs;
}
double YE_A(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//LSE Asymmetric Eigenfunction (y-component)
    double YEa = 0.0;
    if(y >= -c && y < -b)
    {
        YEa = -1.0/2.0*(1.0-exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YEa = 1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YEa = 1.0/2.0*(1.0-exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return  YEa;
}

double dyYEc_S(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{  //Derivative of LSE Symmetric Eigenfunction (y-component), placeholder for complex conjugate
    double YEs = 0.0;
    if(y >= -c && y < -b)
    {
        YEs = -kmd*1.0/2.0*(1.0+exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YEs = kch*1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YEs = kmd*1.0/2.0*(1.0+exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YEs;
}
double dyYEc_A(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Derivative of LSE Asymmetric Eigenfunction (y-component), placeholder for complex conjugate
    double YEa = 0.0;
    if(y >= -c && y < -b)
    {
        YEa = kmd*1.0/2.0*(1.0-exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YEa = kch*1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YEa = kmd*1.0/2.0*(1.0-exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YEa;
}

double IYE_S(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{ //Antiderivative of LSE Symmetric Eigenfunction (y-component)
    double YEs = 0.0;
    if(y >= -c && y < -b)
    {
        YEs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YEs = 1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YEs = -1.0/2.0*(1.0+exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c-y));
        
    }
    return YEs;
}
double IYE_A(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{ //Antiderivative of LSE Asymmetric Eigenfunction (y-component)
    double YEa = 0.0;
    if(y >= -c && y < -b)
    {
        YEa = 1.0/2.0*(1.0-exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YEa = 1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YEa = -1.0/2.0*(1.0-exp(-2.0*kch*b))/(Ep*cos(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YEa;
}

double XE(double x,int n, double w)
{//x-dependence of LSE eigenfunctions
    return sin(M_PI*n/w*x);
}
double dxXE(double x,int n, double w)
{//derivative of x-dependence of LSE eigenfunctions
    return M_PI*n/w*cos(M_PI*n/w*x);
}
//END OF E


//Begin H

double YH_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//LSM Symmetric Eigenfunction (y-component)
    double YHs = 0.0;
    if(y >= -c && y < -b)
    {
        YHs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHs = 1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YHs;
}
double YH_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//LSE Asymmetric Eigenfunction (y-component)
    double YHa = 0.0;
    if(y >= -c && y < -b)
    {
        YHa = -1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHa = 1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHa = 1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YHa;
}


double YHc_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Same as YH_S currently, placeholder for complex conjugate
    double YHs = 0.0;
    if(y >= -c && y < -b)
    {
        YHs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHs = 1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHs = 1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YHs;
}
double YHc_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Same as YH_A currently, placeholder for complex conjugate
    double YHa = 0.0;
    if(y >= -c && y < -b)
    {
        YHa = -1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHa = 1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHa = 1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*sin(kmd*(c-y));
    }
    return YHa;
}


double dyYHc_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Derivative of LSM Symmetric Eigenfunction (y-component), placeholder for complex conjugate
    double YHs = 0.0;
    if(y >= -c && y < -b)
    {
        YHs = kmd*1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHs = kch*1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHs = kmd*1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YHs;
}
double dyYHc_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{ //Derivative of LSM Asymmetric Eigenfunction (y-component), placeholder for complex conjugate
    double YHa = 0.0;
    
    if(y >= -c && y < -b)
    {
        YHa = -kmd*1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHa = kch*1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHa = kmd*1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YHa;
}
double dyYH_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Derivative of LSM Symmetric Eigenfunction (y-component)
    double YHs = 0.0;
    if(y >= -c && y < -b)
    {
        YHs = kmd*1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHs = kch*1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHs = kmd*1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YHs;
}
double dyYH_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Derivative of LSM Asymmetric Eigenfunction (y-component)
    double YHa = 0.0;
    
    if(y >= -c && y < -b)
    {
        YHa = -kmd*1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHa = kch*1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHa = kmd*1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YHa;
}


double IYH_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Antiderivative of LSM Symmetric Eigenfunction (y-component)
    double YHs = 0.0;
    if(y >= -c && y < -b)
    {
        YHs = -1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHs = 1.0/2.0*(1.0-exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHs = -1.0/2.0*(1.0+exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YHs;
}
double IYH_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2,double kch, double kmd)
{//Antiderivative of LSM Asymmetric Eigenfunction (y-component)
    double YHa = 0.0;
    if(y >= -c && y < -b)
    {
        YHa = 1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c+y));
    }
    if(y >= -b && y <= b)
    {
        YHa = 1.0/2.0*(1.0+exp(-2.0*kch*y));
    }
    if(y > b && y <= c)
    {
        YHa = -1.0/2.0*(1.0-exp(-2.0*kch*b))/(Mu*sin(kmd*(c-b)))*cos(kmd*(c-y));
    }
    return YHa;
}


double XH(double x,int n, double w)
{//x-dependence of LSM eigenfunctions
    return cos(M_PI*n/w*x);
}

//END OF H

double ky(double y, double lamb, double b, double c, double w, double Ep, double Mu, double B2, int n)
{//calculation of frequency inside and outside vaccum
    double Ky = 0;
    double kx2 = M_PI*n/w*M_PI*n/w;
    if((y < -b && y > -c) || (y > b && y < c))
        Ky = sqrt((Ep*Mu*B2-1)*lamb-kx2);
    if(y>=-b && y <= b)
        Ky = sqrt((1-B2)*lamb+kx2);
    return Ky;
}
void FindRoots( vector<vector<double> >& zeroes, double (*f)(double,double,double,double,double,double,double,double), int sN, int sI, double Ep, double Mu, double B2, double b, double c, double w, double acc)
{
    //resize eigenvalue matrix
    zeroes.resize(sI);
    for (int i = 0; i < sI; ++i)
    {
        zeroes[i].resize(sN);
    }
    double n = 1;                                        //Iteration over parameter n
    double x1 = 1.0;                                     //Lower x guess
    double x2 = 2.0;                                     //Larger x guess
    double x3 = 0.0;                                     //Current x guess
    double xlast = 0.0;                                  //Previous x guess
    double PN = 0.0;                                     //Test for (P)ositive or (N)egative value of function
    double PNlast = 0.0;                                 //Hold previous value of function
    double step = 0.0;                                   //Step size for guessing next root
    int newRoot = 0;                                     //Flag to reset root finder
    for(int nA=1;nA< sN + 1;nA++)
    {
        n = double(nA);
        x1 = step;
        x2 = x1;
        x3 = x1;
        PN = f(x1,n,Ep,Mu,B2,b,c,w);
        PNlast = f(x1,n,Ep,Mu,B2,b,c,w);
        step= .1;
        xlast=x1;
        int j = 0;
        newRoot = 1;
        //Begin rootfinding search
        while(j < sI)
        {
            if((PN < 0 && PNlast > 0)  || (PN > 0 && PNlast < 0))
            {//If a sign change was detected we have found an apprpriate bracket for root finding
                if(PN < 0)
                {
                    x1 = xlast;
                    x2 = x3;
                }
                else
                {
                    x1 = x3;
                    x2 = xlast;
                }
                while(abs(PN) > acc || PN != PN)
                {//Use bisection method to find root within requested accuracy
                    
                    if(PN > 0)
                    {
                        x1 =x3;
                    }
                    else
                    {
                        x2 = x3;
                    }
                    x3 = (x1+x2)/2.0;
                    PN = f(x3,n,Ep,Mu,B2,b,c,w);
                }
                    //Store root and reset root finder  near next root
                    zeroes[j][nA-1]= x3;
                    j++;
                    x3 += M_PI;
                    newRoot = 0;
                
            }
            
            xlast = x3;
            if(newRoot == 0)
            {//Decide which direction to search for new root
                PN= f(x3,n,Ep,Mu,B2,b,c,w);
                if((PN < 0 && f(x3+.1,n,Ep,Mu,B2,b,c,w) > PN) ||  (PN > 0 && f(x3+.1,n,Ep,Mu,B2,b,c,w) < PN))
                    step = .1;
                else
                    step = -.1;
                newRoot = 1;
            }
            x3 += step;
            PNlast = PN;
            PN = f(x3,n,Ep,Mu,B2,b,c,w);
        }
    }
}


void CalcForceMesh(vector <vector<vector<double> > >& Fz, vector <vector<vector<double> > >& Fx, vector <vector<vector<double> > >& Fy,int xdiv,int ydiv,int zetadiv, double xi, double x0, double yi, double y0, double zi, double zetamin, double zetamax, vector<double>& x,vector<double>& y,vector<double>& zeta, vector<vector<double> >& zeroesES, vector<vector<double> >& zeroesEA, vector<vector<double> >& zeroesHS, vector<vector<double> >& zeroesHA, int sN, int sI, double Ep, double Mu, double B2, double b, double c, double w, double q)
{// Calculate force at each desired spot in space
    
    //Resize vectors
    x.resize(xdiv);
    y.resize(ydiv);
    zeta.resize(zetadiv);
    zeroesES.resize(sI);
    zeroesEA.resize(sI);
    zeroesHS.resize(sI);
    zeroesHA.resize(sI);
    for (int i = 0; i < sI; ++i)
    {
        zeroesES[i].resize(sN);
        zeroesEA[i].resize(sN);
        zeroesHS[i].resize(sN);
        zeroesHA[i].resize(sN);
    }
    //Defining additional vectors to help split up computation
    // (a,b,c) component depends on value of (x,y,z)
    vector<double> aFz(xdiv);
    vector<double> aFx(xdiv);
    vector<double> aFy(xdiv);
    
    vector<double> bFzES(zetadiv);
    vector<double> bFzEA(zetadiv);
    vector<double> bFzHS(zetadiv);
    vector<double> bFzHA(zetadiv);
    vector<double> bFxES(zetadiv);
    vector<double> bFxEA(zetadiv);
    vector<double> bFxHS(zetadiv);
    vector<double> bFxHA(zetadiv);
    vector<double> bFyES(zetadiv);
    vector<double> bFyEA(zetadiv);
    vector<double> bFyHS(zetadiv);
    vector<double> bFyHA(zetadiv);
    
    vector<double> cFzES(zetadiv);
    vector<double> cFzEA(zetadiv);
    vector<double> cFzHS(zetadiv);
    vector<double> cFzHA(zetadiv);
    vector<double> cFxES(zetadiv);
    vector<double> cFxEA(zetadiv);
    vector<double> cFxHS(zetadiv);
    vector<double> cFxHA(zetadiv);
    vector<double> cFyES(zetadiv);
    vector<double> cFyEA(zetadiv);
    vector<double> cFyHS(zetadiv);
    vector<double> cFyHA(zetadiv);
    
    
    //Create Mesh in x
    if(xdiv == 1)
        x[0] = xi;
    else
        x[0] = 0;
    for(int f = 1; f < xdiv; f++)
    {
        x[f] = f*w/(xdiv-1.0);
    }
    //Create Mesh in y
    if(ydiv == 1)
        y[0] = yi;
    else
        y[0] = -c;
    for(int f = 1; f < ydiv; f++)
    {
        y[f] = f*2*c/(ydiv-1.0)-c;
    }
    //Create Mesh in zeta
    if(zetadiv == 1)
        zeta[0] = zi;
    else
        zeta[0] = zetamin;
    for(int f = 1; f < zetadiv; f++)
    {
        zeta[f] = f*(zetamax-zetamin)/(zetadiv-1.0)+zetamin;
    }
    
    
    
    //Loop over x modes
    for(int z = 1; z < sN+1; z++)
    {
        double kx2 = M_PI*z/w*M_PI*z/w;
        for(int f = 0; f < xdiv; f++)
        {
            //Calculate force components only dependent on x
            aFz[f] = XE(x0,z,w)*XE(x[f],z,w);
            aFx[f] = XE(x0,z,w)*dxXE(x[f],z,w);
            aFy[f] = aFz[f];
        }
        
        //Loop over y modes for each x mode
        for(int a = 0;a<sI;a++)
        {
            //Calculate frequency for given eigenvalue
            double kchES = sqrt((1.0-B2)*zeroesES[a][z-1]+kx2);
            double kchEA = sqrt((1.0-B2)*zeroesEA[a][z-1]+kx2);
            double kchHS = sqrt((1.0-B2)*zeroesHS[a][z-1]+kx2);
            double kchHA = sqrt((1.0-B2)*zeroesHA[a][z-1]+kx2);
            double kmdES = sqrt((Ep*Mu*B2-1.0)*zeroesES[a][z-1]-kx2);
            double kmdEA = sqrt((Ep*Mu*B2-1.0)*zeroesEA[a][z-1]-kx2);
            double kmdHS = sqrt((Ep*Mu*B2-1.0)*zeroesHS[a][z-1]-kx2);
            double kmdHA = sqrt((Ep*Mu*B2-1.0)*zeroesHA[a][z-1]-kx2);
            
            //Calculate eigenfunction normalization factors
            double x = kchES*b;
            double As = (2.0*x/(-w/2.0*(Ep*(1.0-Ep*Mu*B2)*1.0/4.0*(1.0+2.0*exp(-2.0*x)+exp(-4.0*x))*(1.0/(Ep*cos(kmdES*(c-b))))*(1.0/(Ep*cos(kmdES*(c-b))))*(c-b+sin(2.0*kmdES*(c-b))/(2.0*kmdES))*2.0*x+(1.0-B2)*b*(1.0/2.0*(1.0-exp(-4.0*x))+2.0*x*exp(-2.0*x)))));
            
            double Aas = (1.0/(-w/2.0*(Ep*(1.0-Ep*Mu*B2)*1.0/4.0*(1.0-2.0*exp(-2.0*kchEA*b)+exp(-4.0*b*kchEA))/(Ep*cos(kmdEA*(c-b)))/(Ep*cos(kmdEA*(c-b)))*(c-b+sin(2.0*kmdEA*(c-b))/(2.0*kmdEA))+(1.0-B2)*(1.0/2.0*(1.0-exp(-4.0*kchEA*b))/(2.0*kchEA)-b*exp(-2.0*b*kchEA)))));
            
            double Bs = (1.0/(-w/2.0*(Mu*(1.0-Ep*Mu*B2)*1.0/4.0*(1.0+2.0*exp(-2.0*kchHS*b)+exp(-4.0*b*kchHS))/(Mu*sin(kmdHS*(c-b)))/(Mu*sin(kmdHS*(c-b)))*(c-b-sin(2.0*kmdHS*(c-b))/(2.0*kmdHS))+(1.0-B2)*(1.0/2.0*(1.0-exp(-4.0*kchHS*b))/(2.0*kchHS)+b*exp(-2.0*b*kchHS)))));
            
            double Bas = (1.0/(-w/2.0*(Mu*(1.0-Ep*Mu*B2)*1.0/4.0*(1.0-2.0*exp(-2.0*kchHA*b)+exp(-4.0*b*kchHA))/(Mu*sin(kmdHA*(c-b)))/(Mu*sin(kmdHA*(c-b)))*(c-b-sin(2.0*kmdHA*(c-b))/(2.0*kmdHA))+(1.0-B2)*(1.0/2.0*(1.0-exp(-4.0*kchHA*b))/(2.0*kchHA)-b*exp(-2.0*b*kchHA)))));
            
            //FIX THIS
            //Calculate force component only dependent on zeta and frequency
            for(int f = 0; f < zetadiv; f++)
            {
                bFzES[f] = cos(sqrt(zeroesES[a][z-1])*abs(zeta[f]));
                bFxES[f] = sin(sqrt(zeroesES[a][z-1])*abs(zeta[f]));
                bFyES[f] = bFxES[f];
                bFzEA[f] = cos(sqrt(zeroesEA[a][z-1])*abs(zeta[f]));
                bFxEA[f] = sin(sqrt(zeroesEA[a][z-1])*abs(zeta[f]));
                bFyEA[f] = bFxEA[f];
                bFzHS[f] = cos(sqrt(zeroesHS[a][z-1])*abs(zeta[f]));
                bFxHS[f] = sin(sqrt(zeroesHS[a][z-1])*abs(zeta[f]));
                bFyHS[f] = bFxHS[f];
                bFzHA[f] = cos(sqrt(zeroesHA[a][z-1])*abs(zeta[f]));
                bFxHA[f] = sin(sqrt(zeroesHA[a][z-1])*abs(zeta[f]));
                bFyHA[f] = bFxHA[f];
            }
            
            
            //Calculate force component dependent on y
            for(int f = 0; f < ydiv; f++)
            {
                //INSIDE VACCUM
                if(abs(y[f])< b && abs(y0)<b)
                {
                    //ELECTRIC COMPONENTS
                    //SYMETRIC
                    cFzES[f] = exp(kchES*(y0+y[f]-2.0*b))*As*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_S(y0,zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES)*IYE_S(y[f],zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES))/(zeroesES[a][z-1]+kx2);
                    cFxES[f] = exp(kchES*(y0+y[f]-2.0*b))*As*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_S(y0,zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES)*IYE_S(y[f],zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES))/(zeroesES[a][z-1]+kx2)/sqrt(zeroesES[a][z-1]);
                    cFyES[f]= exp(kchES*(y0+y[f]-2.0*b))*As*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_S(y0,zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES)*YE_S(y[f],zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES))/(zeroesES[a][z-1]+kx2)/sqrt(zeroesES[a][z-1]);
                    
                    //ASYMETRIC
                    cFzEA[f] = exp(kchEA*(y0+y[f]-2.0*b))*Aas*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_A(y0,zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA)*IYE_A(y[f],zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA))/(zeroesEA[a][z-1]+kx2);
                    cFxEA[f] = exp(kchEA*(y0+y[f]-2.0*b))*Aas*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_A(y0,zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA)*IYE_A(y[f],zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA))/(zeroesEA[a][z-1]+kx2)/sqrt(zeroesEA[a][z-1]);
                    cFyEA[f] = exp(kchEA*(y0+y[f]-2.0*b))*Aas*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_A(y0,zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA)*YE_A(y[f],zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA))/(zeroesEA[a][z-1]+kx2)/sqrt(zeroesEA[a][z-1]);
                    
                    //MAGNETIC COMPONENTS
                    //SYMETRIC
                    cFzHS[f] = exp(kchHS*(y0+y[f]-2.0*b))*Bs*Mu*B2*kx2*YHc_S(y0,zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)*YH_S(y[f],zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)/(zeroesHS[a][z-1]+kx2);
                    
                    cFxHS[f] = exp(kchHS*(y0+y[f]-2.0*b))*Bs*Mu*B2*kx2*YHc_S(y0,zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)*YH_S(y[f],zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)/(zeroesHS[a][z-1]+kx2)/sqrt(zeroesHS[a][z-1]);
                    cFyHS[f] = exp(kchHS*(y0+y[f]-2.0*b))*Bs*Mu*B2*kx2*YHc_S(y0,zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)*dyYH_S(y[f],zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)/(zeroesHS[a][z-1]+kx2)/sqrt(zeroesHS[a][z-1]);
                    //ASYMETRIC
                    cFzHA[f] = exp(kchHA*(y0+y[f]-2.0*b))*Bas*Mu*B2*kx2*YHc_A(y0,zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)*YH_A(y[f],zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)/(zeroesHA[a][z-1]+kx2);
                    cFxHA[f] = exp(kchHA*(y0+y[f]-2.0*b))*Bas*Mu*B2*kx2*YHc_A(y0,zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)*YH_A(y[f],zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)/(zeroesHA[a][z-1]+kx2)/sqrt(zeroesHA[a][z-1]);
                    cFyHA[f] = exp(kchHA*(y0+y[f]-2.0*b))*Bas*Mu*B2*kx2*YHc_A(y0,zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)*dyYH_A(y[f],zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)/(zeroesHA[a][z-1]+kx2)/sqrt(zeroesHA[a][z-1]);
                    
                }
                else
                    if(abs(y[f])> b && abs(y0)<b)//Inside Dielectric
                    {
                        //ELECTRIC COMPONENTS
                        //SYMETRIC
                        cFzES[f] = exp(kchES*(y0-b))*As*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_S(y0,zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES)*IYE_S(y[f],zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES))/(zeroesES[a][z-1]+kx2);
                        cFxES[f] = exp(kchES*(y0-b))*As*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_S(y0,zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES)*IYE_S(y[f],zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES))/(zeroesES[a][z-1]+kx2)/sqrt(zeroesES[a][z-1]);
                        cFyES[f]= exp(kchES*(y0-b))*As*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*ky(y[f],zeroesES[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_S(y0,zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES)*YE_S(y[f],zeroesES[a][z-1],z,Ep,Mu,b,c,w,B2,kchES,kmdES))/(zeroesES[a][z-1]+kx2)/sqrt(zeroesES[a][z-1]);
                        
                        //ASYMETRIC
                        cFzEA[f] = exp(kchEA*(y0-b))*Aas*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_A(y0,zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA)*IYE_A(y[f],zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA))/(zeroesEA[a][z-1]+kx2);
                        cFxEA[f] = exp(kchEA*(y0-b))*Aas*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_A(y0,zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA)*IYE_A(y[f],zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA))/(zeroesEA[a][z-1]+kx2)/sqrt(zeroesEA[a][z-1]);
                        cFyEA[f] = exp(kchEA*(y0-b))*Aas*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*ky(y[f],zeroesEA[a][z-1],b,c,w,Ep,Mu,B2,z)*(dyYEc_A(y0,zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA)*YE_A(y[f],zeroesEA[a][z-1],z,Ep,Mu,b,c,w,B2,kchEA,kmdEA))/(zeroesEA[a][z-1]+kx2)/sqrt(zeroesEA[a][z-1]);
                        
                        //MAGNETIC COMPONENTS
                        //SYMETRIC
                        cFzHS[f] = exp(kchHS*(y0-b))*Bs*Mu*B2*kx2*YHc_S(y0,zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)*YH_S(y[f],zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)/(zeroesHS[a][z-1]+kx2);
                        
                        cFxHS[f] = exp(kchHS*(y0-b))*Bs*Mu*B2*kx2*YHc_S(y0,zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)*YH_S(y[f],zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)/(zeroesHS[a][z-1]+kx2)/sqrt(zeroesHS[a][z-1]);
                        cFyHS[f] = exp(kchHS*(y0-b))*Bs*Mu*B2*kx2*YHc_S(y0,zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)*dyYH_S(y[f],zeroesHS[a][z-1],z,Ep,Mu,b,c,w,B2,kchHS,kmdHS)/(zeroesHS[a][z-1]+kx2)/sqrt(zeroesHS[a][z-1]);
                        //ASYMETRIC
                        cFzHA[f] = exp(kchHA*(y0-b))*Bas*Mu*B2*kx2*YHc_A(y0,zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)*YH_A(y[f],zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)/(zeroesHA[a][z-1]+kx2);
                        cFxHA[f] = exp(kchHA*(y0-b))*Bas*Mu*B2*kx2*YHc_A(y0,zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)*YH_A(y[f],zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)/(zeroesHA[a][z-1]+kx2)/sqrt(zeroesHA[a][z-1]);
                        cFyHA[f] = exp(kchHA*(y0-b))*Bas*Mu*B2*kx2*YHc_A(y0,zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)*dyYH_A(y[f],zeroesHA[a][z-1],z,Ep,Mu,b,c,w,B2,kchHA,kmdHA)/(zeroesHA[a][z-1]+kx2)/sqrt(zeroesHA[a][z-1]);
                        
                    }
            }
            //Multiply to get total force at each point
            for(int f = 0; f < xdiv; f++)
            {
                for(int g = 0; g < ydiv; g++)
                {
                    for(int h = 0; h < zetadiv; h++)
                    {
                        Fx[f][g][h]+= 4.0*M_PI*q*aFx[f]*(bFxES[h]*cFxES[g]+bFxEA[h]*cFxEA[g]+bFxHS[h]*cFxHS[g]+bFxHA[h]*cFxHA[g]);
                        Fy[f][g][h]+= 4.0*M_PI*q*aFy[f]*(bFyES[h]*cFyES[g]+bFyEA[h]*cFyEA[g]+bFyHS[h]*cFyHS[g]+bFyHA[h]*cFyHA[g]);
                        Fz[f][g][h]+= 4.0*M_PI*q*aFz[f]*(bFzES[h]*cFzES[g]+bFzEA[h]*cFzEA[g]+bFzHS[h]*cFzHS[g]+bFzHA[h]*cFzHA[g]);
                    }
                }
            }
        }
        
        
    }
}



#define RoFSuP_h


#endif /* RoFSuP_h */

