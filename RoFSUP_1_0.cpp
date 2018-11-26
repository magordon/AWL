//
//  main.cpp
//  RoFSuP
//
//  Created by MacBook on 8/28/17.
//  Copyright © 2017 Casto. All rights reserved.
//

//
//  main.cpp
//  RoF_1.5
//
//  Created by MacBook on 8/25/17.
//  Copyright © 2017 Casto. All rights reserved.
//

#include <ctime>
#include <math.h>
#include <iostream>
#include <chrono>
#include "RoFSUP.h"
#include <vector>
#include <fstream>
using std::vector;
using std::ifstream;
int main()
{
    //Timing info for root finder
    auto start = std::chrono::high_resolution_clock::now();
    
    //Read input File
    vector<double> input;
    double inputParameter;
    ifstream inputFile ("InputFile.txt");
    if (inputFile.is_open())
    {
        while ( inputFile >> inputParameter )
        {
            input.push_back(inputParameter);
        }
        inputFile.close();
    }
    else std::cout << "Unable to open file";
    double xi =input[0];    //Position of force calculation in x if xdiv = 1
    double yi = input[1];   //Position of force calculation in y if ydiv = 1
    double zi = input[2];   //Position of force calculation in zeta if zetadiv = 1
    double x0 =input[3];    //Beam position x
    double y0 =input[4];    //Beam position y
    int xdiv = input[5];    //Force calculation mesh divisions in x
    int ydiv = input[6];    //Force calculation mesh divisions in y
    int zetadiv = input[7]; //Force calculation mesh divisions in longitudinal direction
    double Ep = input[8];   //Relative Permittivity of dielectric
    double Mu = input[9];   // Relative Permeability of dielectric
    double B = input[10];   // velocity of particle (c=1)
    double B2 = B*B;
    double b = input[11];   //Distance from cavity center to dielectric
    double c = input[12];   //Distance from cavity center to conducting wall (-c < y < c)
    double w = input[13];   //Width of Cavity (0 < x < w)
    int sN = input[14];     //Number of frequcy modes in x
    int sI = input[15];     //Number of frequency modes in y for each mode in x
    double acc = input[16]; //Percision of root finder
    double zetamin = input[17]; //Minimum value of longitudinal coordinate for force calculation
    double zetamax = input[18]; //Maximum value of longitudinal coordinate for force calculation
    double q = -2.998;      //Charge of Particle
    
    
    //Create and size vectors to hold forces.
    vector <vector<vector<double> > > Fz(xdiv);
    vector <vector<vector<double> > > Fx(xdiv);
    vector <vector<vector<double> > > Fy(xdiv);
    for (int i = 0; i < xdiv; ++i)
    {
        Fz[i].resize(ydiv);
        Fx[i].resize(ydiv);
        Fy[i].resize(ydiv);
        for (int j = 0; j < ydiv; ++j)
        {
            Fz[i][j].resize(zetadiv);
            std::fill(Fz[i][j].begin(), Fz[i][j].end(), 0);
            Fx[i][j].resize(zetadiv);
            std::fill(Fx[i][j].begin(), Fx[i][j].end(), 0);
            Fy[i][j].resize(zetadiv);
            std::fill(Fy[i][j].begin(), Fy[i][j].end(), 0);
        }
    }
    
    //Create vectors to hold force calculation coordinates
    vector<double> x(xdiv);
    vector<double> y(ydiv);
    vector<double> zeta(zetadiv);
    //Create and resize vectors to hold eigenvalues of allowed frequency modes
    vector<vector<double> > zeroesES(sI);
    vector<vector<double> > zeroesEA(sI);
    vector<vector<double> > zeroesHS(sI);
    vector<vector<double> > zeroesHA(sI);
    for (int i = 0; i < sI; ++i)
    {
        zeroesES[i].resize(sN);
        zeroesEA[i].resize(sN);
        zeroesHS[i].resize(sN);
        zeroesHA[i].resize(sN);
    }
    //Use dispersion relation to find scaled eigenvalues for each mode
    FindRoots( zeroesES, &f1, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSE Symmetric Modes
    FindRoots( zeroesEA, &f2, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSE Asymmetric Modes
    FindRoots( zeroesHS, &f3, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSM Symmetric Modes
    FindRoots( zeroesHA, &f4, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSM Asymmetric Modes
    //Rescale roots to eigenvalues of each mode
    for(int N = 0; N<sN; N++)
        {
            for(int I = 0; I< sI; I++)
            {
                zeroesES[I][N] = ((zeroesES[I][N]/(c-b)*zeroesES[I][N]/(c-b)+(N+1)*M_PI/w*(N+1)*M_PI/w)/(Ep*Mu*B2-1.0));
                zeroesEA[I][N] = ((zeroesEA[I][N]/(c-b)*zeroesEA[I][N]/(c-b)+(N+1)*M_PI/w*(N+1)*M_PI/w)/(Ep*Mu*B2-1.0));
                zeroesHS[I][N] = ((zeroesHS[I][N]/(c-b)*zeroesHS[I][N]/(c-b)+(N+1)*M_PI/w*(N+1)*M_PI/w)/(Ep*Mu*B2-1.0));
                zeroesHA[I][N] = ((zeroesHA[I][N]/(c-b)*zeroesHA[I][N]/(c-b)+(N+1)*M_PI/w*(N+1)*M_PI/w)/(Ep*Mu*B2-1.0));
            }
        }
    //Timing info for root finder
    auto finish = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time for Root Finder: " << elapsed.count() << " s\n";
    
    //Timing info for force calculator
     auto start2 = std::chrono::high_resolution_clock::now();
    
    //Calculate forces at each desired point in mesh
     CalcForceMesh(Fz,Fx,Fy,xdiv,ydiv,zetadiv,xi,x0,yi,y0,zi,zetamin,zetamax,x,y,zeta,zeroesES,zeroesEA,zeroesHS,zeroesHA,sN,sI,Ep,Mu,B2,b,c,w,q);
    
    //Output to file
     std::ofstream myfile ("FieldOutputFile.txt");
     if (myfile.is_open())
     {
     std::cout << "Writing to file.\n";
     for(int i = 0; i< zetadiv; i++)
     {
     myfile << zeta[i] << "\t" << Fx[0][0][i]<< "\t" << Fy[0][0][i]<< "\t" << Fz[0][0][i] <<"\n";
     }
     
     myfile.close();
     }
     else std::cout << "Unable to open file";
     
     //Timing info for force calculator
     auto finish2 = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> elapsed2 = finish2 - start2;
     std::cout << "Elapsed time for SuP 0.6: " << elapsed2.count() << " s\n";
    
    
    return 0;
}


