# Quast-et-al-C++_code
// The code has been changed due to an error in the original version of the code. The changes do not affect the main results or conclusions. Specifically, in the updated code, we decrease the PV population counts by 1 in lines 835 and 866. In order to compensate for resultant changes in the PV excitability, PV gap junction conductance (gelec) is increased to from 0.01 to 0.05 (line 760), PV applied current is increased from 0 to 0.1 (lines 837 and 868) and GABAa conductance between PV groups (gfsX) in increased from 0.25 to 0.5.
//  QuastEtAlCodeFinal.cpp
//  
//
//  Created by McCarthy, Michelle M on 11/11/22.
//
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream> //this is the library needed for file streaming
#include <stdlib.h>
#include <time.h> // needed to seed srand


/*-- constants --*/

#define  SVCOUNT     37 // number of state variables
#define  SVCOUNT2     41 // number of each state variable
#define  START          0
#define  END       200000 //
#define  TIMESTEP   0.05
#define  PYRCOUNT   40
#define  FSCOUNT    20
#define  FSCOUNT2   20
#define  TCCOUNT    40

/*-- function prototypes --*/

float box_muller(float m, float s);    /* normal random variate generator */
double randgauss( double min, double max, double sigma, double centre);
int IntRK4( int svCount, int svCount2, double sv[SVCOUNT][SVCOUNT2], double time, double dt, double LFPfs[21], double gTC_i[TCCOUNT][FSCOUNT], double gTC_fs2[TCCOUNT][FSCOUNT2], double gi_e[FSCOUNT][PYRCOUNT], double gi2_e[FSCOUNT2][PYRCOUNT], double gTC_e[TCCOUNT][PYRCOUNT]);
void Derive( int svCount, int svCount2, double y[SVCOUNT][SVCOUNT2], double k[SVCOUNT][SVCOUNT2], double time, double LFPfs[21], double gTC_i[TCCOUNT][FSCOUNT], double gTC_fs2[TCCOUNT][FSCOUNT2], double gi_e[FSCOUNT][PYRCOUNT], double gi2_e[FSCOUNT2][PYRCOUNT], double gTC_e[TCCOUNT][PYRCOUNT]);
double rnd();


/*-- main program --*/

using namespace std;



int main( void )

{
    
    /* Output to file */
    ofstream outputFile("results.txt", ios::out);
    if (!outputFile)
    {
        cerr << "File could not be opened" << endl;
        exit(1);
    }
    
    ofstream outputFile2("results2.txt", ios::out);
    ofstream outputFile3("results3.txt", ios::out);
    ofstream outputFile4("results4.txt", ios::out);
    ofstream outputFile5("results5.txt", ios::out);
    ofstream outputFile6("results6.txt", ios::out);
    
    
    ofstream outputFile7("results7.txt", ios::out);
    ofstream outputFile9("results9.txt", ios::out);
    ofstream outputFile12("results12.txt", ios::out);
    ofstream outputFile13("results13.txt", ios::out);
    ofstream outputFile14("results14.txt", ios::out);
    ofstream outputFile15("results15.txt", ios::out);
    
    /* Seed random number generator */
    srand ( time(NULL) );
    
    /* define variables */
    double sv[SVCOUNT][SVCOUNT2];
    
    double gTC_i[TCCOUNT][FSCOUNT]; double gTC_fs2[TCCOUNT][FSCOUNT2];
    double gi_e[FSCOUNT][PYRCOUNT]; double gi2_e[FSCOUNT2][PYRCOUNT];
    double gTC_e[TCCOUNT][PYRCOUNT]; double gTC_ep[TCCOUNT][PYRCOUNT];
    double gTC_ip[TCCOUNT][FSCOUNT]; double gTC_i2p[TCCOUNT][FSCOUNT2];
    
    double svp[PYRCOUNT]; double ispk[PYRCOUNT]; double ispktime[PYRCOUNT];
    double svp10[FSCOUNT]; double ispk10[FSCOUNT]; double ispktime10[FSCOUNT];
    double svp10_2[FSCOUNT2]; double ispk10_2[FSCOUNT2]; double ispktime10_2[FSCOUNT2]; double ispktime10_2p[FSCOUNT2][2];
    double svp30[TCCOUNT]; double ispk30[TCCOUNT]; double ispktime30[TCCOUNT];
    
    double LFPfs[21];
    
   // double wFS2m[TCCOUNT][FSCOUNT2]; double wFSm[TCCOUNT][FSCOUNT];
    
    int z; int i; int n;
    
    double time;
    int itime; int itime2;
    int ntime;
    
    double zd = z;
    double testMax;
    double testMaxI; double testMaxI2;
    
    //double r;
    //r = 0;
    
    int intTime = time; int numTime;
    numTime = intTime % 2;
    
    /**** Initial conditions for Pyramidal cells ****/
    for (z=0; z<PYRCOUNT; z++)
    {
        sv[0][z] = -63;
        sv[1][z] = 0.32*(sv[0][z]+54)/(1-exp(-(sv[0][z]+54)/4)) / (0.32*(sv[0][z]+54)/(1-exp(-(sv[0][z]+54)/4))  + 0.28*(sv[0][z]+27)/(exp((sv[0][z]+27)/5)-1));
        sv[2][z] = 0.128*exp(-(sv[0][z]+50)/18) / (0.128*exp(-(sv[0][z]+50)/18)  +  4/(1+exp(-(sv[0][z]+27)/5)));
        sv[3][z] = 0.032*(sv[0][z]+52)/(1-exp(-(sv[0][z]+52)/5)) / (0.032*(sv[0][z]+52)/(1-exp(-(sv[0][z]+52)/5))   +  0.5*exp(-(sv[0][z]+57)/40));
        sv[4][z] = 3.209*0.0001*(sv[0][z]+30)/(1-exp(-(sv[0][z]+30)/9)) / (3.209*0.0001*(sv[0][z]+30)/(1-exp(-(sv[0][z]+30)/9))   -  3.209*0.0001*(sv[0][z]+30)/(1-exp((sv[0][z]+30)/9)));
        sv[5][z] = 0.0001;
        svp[z] = 0; ispk[z] = 0; ispktime[z] = 0;
        
        for (n=0; n<TCCOUNT; n++){
            if (z == n)
                gTC_e[n][z] = 0.12;
            else
                gTC_e[n][z] = .0001;
        }
    }
    
    /**** Initial conditions for FS cells ****/
    for (z=0; z<FSCOUNT; z++)
    {
        sv[10][z] = -63;// - 20*sqrt(0.05)*box_muller(0,1);
        sv[11][z] = 0.32*(sv[10][z]+54)/(1-exp(-(sv[10][z]+54)/4)) / (0.32*(sv[10][z]+54)/(1-exp(-(sv[10][z]+54)/4))  + 0.28*(sv[10][z]+27)/(exp((sv[10][z]+27)/5)-1));
        sv[12][z] = 0.128*exp(-(sv[10][z]+50)/18) / (0.128*exp(-(sv[10][z]+50)/18)  +  4/(1+exp(-(sv[10][z]+27)/5)));
        sv[13][z] = 0.032*(sv[10][z]+52)/(1-exp(-(sv[10][z]+52)/5)) / (0.032*(sv[10][z]+52)/(1-exp(-(sv[10][z]+52)/5))   +  0.5*exp(-(sv[10][z]+57)/40));
        sv[14][z] = 0.0001;
        svp10[z] = 0; ispk10[z] = 0; ispktime10[z] = 0;
        
        for (n=0; n<PYRCOUNT; n++){
            gi_e[z][n] = 0.01; // 0.1
        }
    }
    
    
    /**** Initial conditions for FS2 cells ****/
    for (z=0; z<FSCOUNT2; z++)
    {
        sv[15][z] = -63; // - 20*sqrt(0.05)*box_muller(0,1);
        sv[16][z] = 0.32*(sv[15][z]+54)/(1-exp(-(sv[15][z]+54)/4)) / (0.32*(sv[15][z]+54)/(1-exp(-(sv[15][z]+54)/4))  + 0.28*(sv[15][z]+27)/(exp((sv[15][z]+27)/5)-1));
        sv[17][z] = 0.128*exp(-(sv[15][z]+50)/18) / (0.128*exp(-(sv[15][z]+50)/18)  +  4/(1+exp(-(sv[15][z]+27)/5)));
        sv[18][z] = 0.032*(sv[15][z]+52)/(1-exp(-(sv[15][z]+52)/5)) / (0.032*(sv[15][z]+52)/(1-exp(-(sv[15][z]+52)/5))   +  0.5*exp(-(sv[15][z]+57)/40));
        sv[19][z] = 0.0001;
        svp10_2[z] = 0; ispk10_2[z] = 0; ispktime10_2[z] = 0;
        
        for (n=0; n<PYRCOUNT; n++){
            gi2_e[z][n] = 0.01;
        }
    }
    
    for (z=0; z<21; z++)
        LFPfs[z] = 0;

    
    /**** Initial conditions for TC input cell ****/
    for (z=0; z<TCCOUNT; z++)
    {
        sv[30][z] = -63;
        sv[31][z] = 0.32*(sv[30][z]+54)/(1-exp(-(sv[30][z]+54)/4)) / (0.32*(sv[30][z]+54)/(1-exp(-(sv[30][z]+54)/4))  + 0.28*(sv[30][z]+27)/(exp((sv[30][z]+27)/5)-1));
        sv[32][z] = 0.128*exp(-(sv[30][z]+50)/18) / (0.128*exp(-(sv[30][z]+50)/18)  +  4/(1+exp(-(sv[30][z]+27)/5)));
        sv[33][z] = 0.032*(sv[30][z]+52)/(1-exp(-(sv[30][z]+52)/5)) / (0.032*(sv[30][z]+52)/(1-exp(-(sv[30][z]+52)/5))   +  0.5*exp(-(sv[30][z]+57)/40));
        sv[34][z] = 0.0001;
        sv[35][z] = 0.0001;
        
        svp30[z] = 0; ispk30[z] = 0; ispktime30[z] = 0;
        
        if (z < FSCOUNT){
            for (n=0; n<FSCOUNT; n++){
                if (z == n)
                    gTC_i[z][n] = 0.25; // 0.25
                else
                    gTC_i[z][n] = 0.001;
                
              //  wFSm[z][n] = 0;
            }
            
        }
        
        if (z >= FSCOUNT){
            for (n=0; n<FSCOUNT2; n++){
                if (z-FSCOUNT == n){
              //      gTC_fs2[z][n] = 0.15 - 0.002*n; //0.25
                    gTC_fs2[z][n] = 0.25;
                }
                else
                    gTC_fs2[z][n] = 0.001;
                
             //   wFS2m[z][n] = 0;
            }
        }
       
    }

    //*** Sum of maximal conductances of TC to E synapses ****//
    double sumgTCe[PYRCOUNT];
    for (z=0; z<PYRCOUNT; z++){
        sumgTCe[z] = 0;
        for (n=0; n<TCCOUNT; n++)
            sumgTCe[z] = sumgTCe[z] + gTC_e[n][z];
    }
    double sumgTCi[FSCOUNT];
    for (z=0; z<FSCOUNT; z++){
        sumgTCi[z] = 0;
        for (n=0; n<TCCOUNT/2; n++)
            sumgTCi[z] = sumgTCi[z] + gTC_i[n][z];
    }
    double sumgTCi2[FSCOUNT2];
    for (z=0; z<FSCOUNT2; z++){
        sumgTCi2[z] = 0;
        for (n=0; n<TCCOUNT/2; n++)
            sumgTCi2[z] = sumgTCi2[z] + gTC_fs2[n+FSCOUNT][z];
    }
    cout<<sumgTCe[0]<<"  ,  "<<sumgTCi[0]<<"  ,  "<<sumgTCi2[1]<<endl;
    
    time = START;
    
    /* start main time loop - loop until stopping time reached */
    while ( time <  END )
    {
        ntime = time; // makes time an integer
        itime = ntime % 5;  // itime is an integer mod 5 --> (fs = 200), mod2 --> fs = 500
        
        
        IntRK4(  SVCOUNT, SVCOUNT2, sv, time, TIMESTEP, LFPfs, gTC_i, gTC_fs2, gi_e, gi2_e, gTC_e);
        
        
        outputFile<<time<<","<<sv[10][18]<<"  ,  "<<sv[10][19]<<","<<sv[15][18]<<","<<sv[15][19]<<endl;
        
        
        double maxGtc_e = 0.17;
        double maxGtc_i = 0.5;
        double maxEI = 0.5;
        double maxIE = 0.5;
        double minIE = 0.001;
        double minEI = 0.001;
        double minGtc_i = 0.001;
        double minGtc_e = 0.001;
        
        double dt;
        double w;
        
        double sumMaxI = 0.5;
        
        /********** Raster Plot for TC Cells **********/
        for (z=0; z<TCCOUNT; z++)
        {
            if (time >= 1)
            {
                if (sv[30][z] > 0 && svp30[z] <= 0)
                {
                    ispk30[z] = z;
                    ispktime30[z] = time;
                    outputFile6<<ispk30[z]<<","<<ispktime30[z]<<endl;
                    
                    //*** Update TC to FS2 synaptic strengths ****//
                    
                    if (z >= FSCOUNT){
                        for (n=0; n<FSCOUNT2; n++){
                            gTC_i2p[z][n] = gTC_fs2[z][n];
                            dt = ispktime10_2[n] - ispktime30[z]; // dt will be -
                            w = -2.5*0.001*exp(dt/90);
                            gTC_fs2[z][n] = gTC_fs2[z][n] + w*gTC_fs2[z][n];
                            if (gTC_fs2[z][n] < minGtc_i) gTC_fs2[z][n] = minGtc_i;
                            //   if (gTC_fs2[z][n] > maxGtc_i) gTC_fs2[z][n] = maxGtc_i;
                            outputFile12<<z<<","<<n<<","<<ispktime30[z]<<","<<gTC_fs2[z][n]<<","<<dt<<","<<w<<","<<sumgTCi2[n]<<endl;
                            sumgTCi2[n] = sumgTCi2[n] - gTC_i2p[z][n] + gTC_fs2[z][n];
                        }
                    }
                    
                    
                    //*** Update TC to FS1 synaptic strengths ****//
                    
                    if (z < FSCOUNT){
                        for (n=0; n<FSCOUNT; n++){
                            gTC_ip[z][n] = gTC_i[z][n];
                            dt = ispktime10[z] - ispktime30[z]; // dt will be -
                            w = -2.5*0.001*exp(dt/90);
                            gTC_i[z][n] = gTC_i[z][n] + w*gTC_i[z][n];
                            if (gTC_i[z][n] < minGtc_i) gTC_i[z][n] = minGtc_i;
                            //      if (gTC_i[z][n] > maxGtc_i) gTC_i[z][n] = maxGtc_i;
                            outputFile7<<z<<","<<n<<","<<ispktime30[z]<<","<<gTC_i[z][n]<<","<<dt<<","<<w<<","<<sumgTCi[n]<<endl;
                            sumgTCi[n] = sumgTCi[n] - gTC_ip[z][n] + gTC_i[z][n];
                        }
                    }
                    
                    
                    //*** TC to E Plasticity ****//
                    for (n=0; n<PYRCOUNT; n++){
                        gTC_ep[z][n] = gTC_e[z][n];
                        dt = ispktime[n] - ispktime30[z]; // dt is -
                        w = -2.5*0.001*exp(dt/34);
                        gTC_e[z][n] = gTC_e[z][n] + w*gTC_e[z][n];
                        if (gTC_e[z][n] < minGtc_e) gTC_e[z][n] = minGtc_e;
                        if (gTC_e[z][n] > maxGtc_e) gTC_e[z][n] = maxGtc_e;
                        outputFile15<<z<<","<<n<<","<<ispktime30[z]<<","<<gTC_e[z][n]<<","<<dt<<","<<sumgTCe[8]<<endl;
                        sumgTCe[n] = sumgTCe[n] - gTC_ep[z][n] + gTC_e[z][n];
                    }
                    
                    
                    
                }
            }
        }
        
        
        /********** Raster Plot for Pyramidal Cells **********/
        for (z=0; z<PYRCOUNT; z++)
        {
            if (time >= 1)
            {
                if (sv[0][z] > 0 && svp[z] <= 0)
                {
                    ispk[z] = z;
                    ispktime[z] = time;
                    outputFile4<<ispk[z]<<","<<ispktime[z]<<endl;
                    
                    //*** FS1 to PYR plasticity ****//
                    for (n=0; n<FSCOUNT; n++){
                        dt = ispktime[z] - ispktime10[n]; // dt will be +
                        w = -5*0.001*exp(-dt/14);
                        gi_e[n][z] = gi_e[n][z] + w*gi_e[n][z];
                        if (gi_e[n][z] < minIE) gi_e[n][z] = minIE;
                        if (gi_e[n][z] > maxIE) gi_e[n][z] = maxIE;
                        outputFile9<<n<<","<<z<<","<<ispktime10[n]<<","<<gi_e[n][z]<<","<<dt<<","<<w<<endl;
                    }
                    
                    //*** FS2 to PYR plasticity ****//
                    for (n=0; n<FSCOUNT2; n++){
                        dt = ispktime[z] - ispktime10_2[n]; // dt will be +
                        w = -5*0.001*exp(-dt/14);
                        gi2_e[n][z] = gi2_e[n][z] + w*gi2_e[n][z];
                        if (gi2_e[n][z] < minIE) gi2_e[n][z] = minIE;
                        if (gi2_e[n][z] > maxIE) gi2_e[n][z] = maxIE;
                        outputFile14<<n<<","<<z<<","<<ispktime10_2[n]<<","<<gi2_e[n][z]<<","<<dt<<","<<w<<endl;
                    }
                    
                    
                    //*** TC to E plasticity ****//
                    if (sumgTCe[z] < maxGtc_e){
                        for (n=0; n<TCCOUNT; n++){
                            gTC_ep[n][z] = gTC_e[n][z];
                            dt = ispktime[z] - ispktime30[n]; // dt will be +
                            w = 5*0.001*exp(-dt/14);
                            gTC_e[n][z] = gTC_e[n][z] + w*gTC_e[n][z];
                            testMax = sumgTCe[z] - gTC_ep[n][z] + gTC_e[n][z];
                            if (testMax > maxGtc_e)
                                gTC_e[n][z] = gTC_ep[n][z];
                            if (gTC_e[n][z] < minGtc_e) gTC_e[n][z] = minGtc_e;
                            //     if (gTC_e[n][z] > maxGtc_e) gTC_e[n][z] = maxGtc_e;
                            outputFile15<<n<<","<<z<<","<<ispktime[z]<<","<<gTC_e[n][z]<<","<<dt<<","<<sumgTCe[8]<<endl;
                            sumgTCe[z] = sumgTCe[z] - gTC_ep[n][z] + gTC_e[n][z];
                        }
                    }
                    
                    
                }
            }
        }
        
        /********** Raster Plot for FS2 Cells **********/
        for (z=0; z<FSCOUNT2; z++)
        {
            if (time >= 1)
            {
                if (sv[15][z] > 0 && svp10_2[z] <= 0)
                {
                    ispk10_2[z] = z;
                    ispktime10_2p[z][1] = ispktime10_2[z];
                    ispktime10_2[z] = time;
                    ispktime10_2p[z][0] = time;
                    outputFile5<<ispk10_2[z]<<","<<ispktime10_2[z]<<endl;
                    
                    
                    //*** TC-to-FS2 plasticity ****//
                    
                    for (n=TCCOUNT/2; n<TCCOUNT; n++){
                        if (sumgTCi2[z] < sumMaxI){
                            gTC_i2p[n][z] = gTC_fs2[n][z];
                            dt = ispktime10_2[z] - ispktime30[n]; //dt will be +
                            w = 5*0.001*exp(-dt/14);
                            gTC_fs2[n][z] = gTC_fs2[n][z] + w*gTC_fs2[n][z];
                            testMaxI2 = sumgTCi2[z] - gTC_i2p[n][z] + gTC_fs2[n][z];
                            if (testMaxI2 > sumMaxI)
                                gTC_fs2[n][z] = gTC_i2p[n][z];
                            if (gTC_fs2[n][z] < minGtc_i) gTC_fs2[n][z] = minGtc_i;
                            outputFile12<<n<<","<<z<<","<<ispktime10_2[z]<<","<<gTC_fs2[n][z]<<" , "<<dt<<","<<w<<","<<sumgTCi2[z]<<endl;
                            sumgTCi2[z] = sumgTCi2[z] - gTC_i2p[n][z] + gTC_fs2[n][z];
                        }
                    }
                    
                    
                    //*** FS2-to-PYR plasticity ****//
                    for (n=0; n<PYRCOUNT; n++){
                        dt = ispktime[n] - ispktime10_2[z]; // dt will be negative
                        w = 2.5*0.001*exp(dt/34); // lateral? inhibition
                        gi2_e[z][n] = gi2_e[z][n] + w*gi2_e[z][n];
                        if (gi2_e[z][n] < minIE) gi2_e[z][n] = minIE;
                        if (gi2_e[z][n] > maxIE) gi2_e[z][n] = maxIE;
                        outputFile14<<z<<","<<n<<","<<ispktime10_2[z]<<","<<gi2_e[z][n]<<","<<dt<<","<<w<<endl;
                    }
                    
                    
                }
                
            }
        }
        
        
        
        /********** Raster Plot for FS Cells **********/
        for (z=0; z<FSCOUNT; z++)
        {
            if (time >= 1)
            {
                if (sv[10][z] > 0 && svp10[z] <= 0)
                {
                    ispk10[z] = z;
                    ispktime10[z] = time;
                    outputFile2<<ispk10[z]<<","<<ispktime10[z]<<endl;
                    
                    
                    //*** TC-to-FS1 plasticity Spread ****//
                        for (n=0; n<TCCOUNT/2; n++){

                             if (sumgTCi[z] < sumMaxI){
                                gTC_ip[n][z] = gTC_i[n][z];
                                dt = ispktime10[z] - ispktime30[n]; //dt will be +
                                    w = 5*0.001*exp(-dt/14);
                                gTC_i[n][z] = gTC_i[n][z] + w*gTC_i[n][z];
                                testMaxI = sumgTCi[z] - gTC_ip[n][z] + gTC_i[n][z];
                                   if (testMaxI > sumMaxI)
                                       gTC_i[n][z] = gTC_ip[n][z];
                                    if (gTC_i[n][z] < minGtc_i) gTC_i[n][z] = minGtc_i;
                                //     if (gTC_e[n][z] > maxGtc_e) gTC_e[n][z] = maxGtc_e;
                                    
                                outputFile7<<n<<","<<z<<","<<ispktime10[z]<<","<<gTC_i[n][z]<<" , "<<dt<<","<<w<<","<<sumgTCi[z]<<endl;
                                    
                                sumgTCi[z] = sumgTCi[z] - gTC_ip[n][z] + gTC_i[n][z];
                                 
                            }
                    }
                    
                    
                    //*** FS1-to-PYR plasticity ****//
                    for (n=0; n<PYRCOUNT; n++){
                        dt = ispktime[n] - ispktime10[z]; // dt will be negative
                        w = 2.5*0.001*exp(dt/34);
                        gi_e[z][n] = gi_e[z][n] + w*gi_e[z][n];
                        if (gi_e[z][n] < minIE) gi_e[z][n] = minIE;
                        if (gi_e[z][n] > maxIE) gi_e[z][n] = maxIE;
                        outputFile9<<z<<","<<n<<","<<ispktime10[z]<<","<<gi_e[z][n]<<","<<dt<<","<<w<<endl;
                    }
                    
                    
                }
            }
        }
        
        
        
      
        
        /********** Filtered LFP Output for FSIs **********/
        if(time >= 0)
        {
            if ( itime == 0 && itime2 != 0)
            {
                outputFile3<<time<<","<<LFPfs[0]<<endl;
            }
        }

        
        for (z=0; z<FSCOUNT; z++)
            svp10[z] = sv[10][z];
        for (z=0; z<FSCOUNT2; z++)
            svp10_2[z] = sv[15][z];
        for (z=0; z<TCCOUNT; z++)
            svp30[z] = sv[30][z];
        for (z=0; z<PYRCOUNT; z++)
            svp[z] = sv[0][z];
        
        itime2 = itime;
        
        time += TIMESTEP;
    }
    
    return 0;
}


/*---------------------- Integration Methods ----------------------
 
 *-- IntRK4() -------------------------------------------------------
 *
 *--- general purpose Runge-Kutta integrator (fourth order)
 *    code is modified from "Numerical Recipes"
 *
 *--INPUTS
 *   int     svCount    = state variable count
 *   float  *stateVar   = state variable value array
 *   double *derivative = derivative value array
 *   float  *time      = time variable address
 *   float   dt         = timestep
 *-------------------------------------------------------------------
 */

int IntRK4( int svCount, int svCount2, double sv[SVCOUNT][SVCOUNT2], double time, double dt, double LFPfs[21], double gTC_i[TCCOUNT][FSCOUNT], double gTC_fs2[TCCOUNT][FSCOUNT2], double gi_e[FSCOUNT][PYRCOUNT], double gi2_e[FSCOUNT2][PYRCOUNT], double gTC_e[TCCOUNT][PYRCOUNT])
{
    int  i;                      /* loop counter */
    int  j;/* loop counter */
    
    /*-- temporary arrays --*/
    static double y [ SVCOUNT ][ SVCOUNT2 ];
    static double k1[ SVCOUNT ][ SVCOUNT2 ];
    static double k2[ SVCOUNT ][ SVCOUNT2 ];
    static double k3[ SVCOUNT ][ SVCOUNT2 ];
    static double k4[ SVCOUNT ][ SVCOUNT2 ];
    
    
    for ( i=0; i < SVCOUNT; i++ )       /* remember initial values */
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j];
        }
    }
    
    /*-- pass 1 --*/
    Derive( SVCOUNT, SVCOUNT2, y, k1, time, LFPfs, gTC_i, gTC_fs2, gi_e, gi2_e, gTC_e);     /* get derivatives */
    
    for ( i=0; i < SVCOUNT; i++ )       /* update state variables */
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j] + k1[ i ][j] * dt/2;
        }
    }
    
    /*-- pass 2 --*/
    time += dt/2;
    Derive( SVCOUNT, SVCOUNT2, y, k2, time, LFPfs, gTC_i, gTC_fs2, gi_e, gi2_e, gTC_e);     /* get derivatives */
    
    for ( i=0; i < SVCOUNT; i++ )       /* update state variables */
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j] + k2[ i ][j] * dt/2;
        }
    }
    
    /*- pass 3 --*/
    Derive( SVCOUNT, SVCOUNT2, y, k3, time, LFPfs, gTC_i, gTC_fs2, gi_e, gi2_e, gTC_e);     /* get derivatives */
    
    for ( i=0; i < SVCOUNT; i++ )
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j] + k3[ i ][j] * dt;    /* update state variables */
        }
    }
    
    /*-- pass 4 --*/
    time += dt/2;
    Derive( SVCOUNT, SVCOUNT2, y, k4, time, LFPfs, gTC_i, gTC_fs2, gi_e, gi2_e, gTC_e);  /* get derivatives */
    for ( i=0; i < SVCOUNT; i++ )
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            sv[ i ][j] = sv[ i ][j]  + ( k1[ i ][j] + 2*k2[ i ][j] + 2*k3[ i ][j] +k4[ i ][j] ) * dt/6;
        }
    }
    
    return 1;
}


void Derive( int svCount, int svCount2, double y[SVCOUNT][SVCOUNT2], double k[SVCOUNT][SVCOUNT2], double time, double LFPfs[21], double gTC_i[TCCOUNT][FSCOUNT], double gTC_fs2[TCCOUNT][FSCOUNT2], double gi_e[FSCOUNT][PYRCOUNT], double gi2_e[FSCOUNT2][PYRCOUNT], double gTC_e[TCCOUNT][PYRCOUNT])
{
    
    int z;
    int i;
    
    
    /********** All-to-all FSI-to-FSI and Individual FSI synapse **********/
    double sgiFS[FSCOUNT]; double tauiFS = 10; // 10 for regular, 20 for a1 KO
    
    double Sum14 = 0;
    for (z=0; z<FSCOUNT; z++){
        sgiFS[z] = 2*(1+tanh(y[10][z]/4));
        k[ 14 ][z] = sgiFS[z]*(1-y[14][z])-y[14][z]/tauiFS;
        Sum14 = Sum14 + y[14][z];
    }
    
    /********** All-to-all FSI2-to-FSI2 and Individual FSI synapse **********/
    double sgiFS2[FSCOUNT2]; double tauiFS2 = 10;
    
    double Sum19 = 0;
    for (z=0; z<FSCOUNT2; z++){
        sgiFS2[z] = 2*(1+tanh(y[15][z]/4));
        k[ 19 ][z] = sgiFS2[z]*(1-y[19][z])-y[19][z]/tauiFS2;
        Sum19 = Sum19 + y[19][z];
    }
    
    /************** Sum of all FS voltages *****************/
    double SumFS = 0; // sum of all FS voltages
    for (z=0; z<FSCOUNT; z++)
        SumFS = SumFS + y[10][z];
    
    /************** Sum of all FS2 voltages *****************/
    double SumFS2 = 0; // sum of all FS voltages
    for (z=0; z<FSCOUNT2; z++)
        SumFS2 = SumFS2 + y[15][z];
    
    /******* Input from TC cells to E cells ********/
    double SumTC = 0; double SumTC1 = 0; double SumTC2 = 0;// sum of all TC input
    for (z=0; z<TCCOUNT; z++){
        SumTC = SumTC + y[34][z]; // input from all TC cells
        if (z >= PYRCOUNT/2)
            SumTC2 = SumTC2 + y[34][z]; // input from second half of TC cells
        else
            SumTC1 = SumTC1 + y[34][z]; // input from first half of TC cells
    }
    
    /******* Input from TC cells to I cells ********/
    double SumTCi = 0; double SumTC1i = 0; double SumTC2i = 0;// sum of all TC input
    for (z=0; z<TCCOUNT; z++){
        SumTCi = SumTCi + y[35][z]; // input from all TC cells
        if (z >= TCCOUNT/2)
            SumTC2i = SumTC2i + y[35][z]; // input from second half of TC cells
        else
            SumTC1i = SumTC1i + y[35][z]; // input from first half of TC cells
    }
    
    //*** gI-E **** //
    double gie[PYRCOUNT];
    for (z=0; z<PYRCOUNT; z++)
        gie[z] = 0;
    for (z=0; z<PYRCOUNT; z++){
        for (i=0; i<FSCOUNT; i++)
            gie[z] = gie[z] + gi_e[i][z]*y[14][i];
    }
    
    //*** gI2-E **** //
    double gi2e[PYRCOUNT];
    for (z=0; z<PYRCOUNT; z++)
        gi2e[z] = 0;
    for (z=0; z<PYRCOUNT; z++){
        for (i=0; i<FSCOUNT2; i++)
            gi2e[z] = gi2e[z] + gi2_e[i][z]*y[19][i];
    }
    
    //*** gTC-E **** //
    double gTCe[PYRCOUNT];
    for (z=0; z<PYRCOUNT; z++)
        gTCe[z] = 0;
    for (z=0; z<PYRCOUNT; z++){
        for (i=0; i<TCCOUNT; i++)
            gTCe[z] = gTCe[z] + gTC_e[i][z]*y[34][i];
    }
    
    //*** gTC-FS1 **** //
    double gTCi[FSCOUNT];
    for (z=0; z<FSCOUNT; z++)
        gTCi[z] = 0;
    for (z=0; z<FSCOUNT; z++){
        for (i=0; i<TCCOUNT/2; i++)
            gTCi[z] = gTCi[z] + gTC_i[i][z]*y[35][i];
    }
    
    
    //*** gTC-FS2 **** //
    double gTCfs2[FSCOUNT2];
    for (z=0; z<FSCOUNT2; z++)
        gTCfs2[z] = 0;
    for (z=0; z<FSCOUNT2; z++){
        for (i=TCCOUNT/2; i<TCCOUNT; i++){
            gTCfs2[z] = gTCfs2[z] + gTC_fs2[i][z]*y[35][i];
        }
    }
    
    
    
    
    
    /***************** Pyramidal Cells *******************/
    double iapp[PYRCOUNT];
    LFPfs[0] = 0;
    for (z=0; z<PYRCOUNT; z++)
    {
        
        double gna = 100; // mS/cm^2
        double gk = 80;
        double gm = 2; //
        double gl = 0.1;
        double Ena = 50;
        double Ek  = -100;
        double El  = -67;
        double Q10 = 3.209;
        double vhalf = -30;
        double ge_e = 0.1;
        double gi_e = 0.2;
        double Ei = -80;
        double Ee = 0;
        
        iapp[z] = 0;
        
        k[ 0 ][z] = -gna*y[2][z]*(y[1][z]*y[1][z]*y[1][z])*(y[0][z]-Ena) -gk*(y[3][z]*y[3][z]*y[3][z]*y[3][z])*(y[0][z]-Ek)
        -gm*y[4][z]*(y[0][z]-Ek)
        -gl*(y[0][z]-El)
        -gie[z]/FSCOUNT*(y[0][z] - Ei)
        -gi2e[z]/FSCOUNT2*(y[0][z] - Ei)
        -gTCe[z]*(y[0][z] - 0)
        +0*sqrt(0.05)*box_muller(0,1);
        //  + iapp[z];
        
        ///* // rate equations
        k[ 1 ][z] = 0.32*(y[0][z]+54)/(1-exp(-(y[0][z]+54)/4))*(1-y[1][z])-0.28*(y[0][z]+27)/(exp((y[0][z]+27)/5)-1)*y[1][z];     // derivative equation for m
        k[ 2 ][z] = 0.128*exp(-(y[0][z]+50)/18)*(1-y[2][z])-4/(1+exp(-(y[0][z]+27)/5))*y[2][z];   // derivative equation for h
        k[ 3 ][z] = 0.032*(y[0][z]+52)/(1-exp(-(y[0][z]+52)/5))*(1-y[3][z])-0.5*exp(-(y[0][z]+57)/40)*y[3][z];  // derivative equation for n
        k[ 4 ][z] = (Q10*0.0001*(y[0][z]-vhalf)/(1-exp(-(y[0][z]-vhalf)/9)))*(1-y[4][z]) + Q10*0.0001*(y[0][z]-vhalf)/(1-exp((y[0][z]-vhalf)/9))*y[4][z];  // derivative equation for w (for m-current) -- Mainen and Sejowski
        
        LFPfs[0] =  gie[z]/FSCOUNT + gi2e[z]/FSCOUNT2 + LFPfs[0];
        
    }
    
    
    /****************** FS Cells ************************/
    
    /*** FS equations ****/
    
    double FSgna = 100; // conductances in units of mS/cm^2
    double FSgk = 80;
    double FSgl = 0.1;
    double fsEna = 50; // reversal potentials in units of mV
    double fsEk  = -100;
    double fsEl  = -67;
    
    double gfs_fs = 0.1;
    double gfs0_fs = 1;
    double ge_i = 0.7;
    double ge_i2 = 0.7;
    double iEi = -80; // mV
    double gelec = 0;
    double iEe = 0;
    double gtc_fs1 = 0.9;
    double gtc_fs2 = 0.9;
    
    double gabaFS[FSCOUNT];
    double gabaFSX[FSCOUNT];
    double IappFS[FSCOUNT];
    double IappFS2[FSCOUNT2];
    double gabaFS2[FSCOUNT2];
    double gabaFSX2[FSCOUNT2];
    
    double gfsX = 0.5; // 0.2
    gelec = 0.05; // 0.01
    gfs_fs = 0.03; // 0.03
    gfsX = 0.25; // 0.25
    
  
    // To simulate bouts of "sleep"
    /*
    if (time > 80000 && time < 95000)
            FSgna = 10;
      //  fsEl = -90;
    
    if (time > 85000)
        gelec = 0.01 - 0.01*15/100;
    
    if (time > 125000 && time < 135000)
        FSgna = 10;
    
    if (time > 135000)
        gelec = 0.0072;
   */
    
    /*
    if (time > 60000 && time < 65000)
        FSgna = 10;
    //  fsEl = -90;
    
    if (time > 65000)
      //  gelec = 0.01 - 0.01*15/100;
        gfs_fs = 0.05;

    if (time > 100000 && time < 105000)
        FSgna = 10;
    /
    if (time > 100000)
      // gelec = 0.0072;
    gfs_fs = 0.1;
    */
    /*
    if (time > 80000 && time < 85000)
        FSgna = 10;

    if (time > 85000)
      //   gelec = 0.01 - 0.01*15/100;
        gelec = 0.0085;

    if (time > 120000 && time < 130000)
        FSgna = 10;

    if (time > 130000)
        //   gelec = 0.01 - 0.01*15/100;
        gelec = 0.001;
    */
    
    
    for (z=0; z<FSCOUNT; z++)
    {
        IappFS[z] = 0;
        
        if (time <= 20000){ //20000
            gelec = 0;
            gfs_fs = 0.03/10;
            gfsX = 0.5/10; // changed from gfsX = 0.25/10;
        }
        
        
        // FSI to FSI connections
        gabaFS[z] = gfs_fs/(FSCOUNT-1.0)*(Sum14-y[14][z])*(y[10][z]-iEi);
        gabaFSX[z] = gfsX/(FSCOUNT2)*(Sum19)*(y[10][z]-iEi);
        
        k[ 10 ][z] =
        -FSgna*y[12][z]*(y[11][z]*y[11][z]*y[11][z])*(y[10][z]-fsEna)
        -FSgk*(y[13][z]*y[13][z]*y[13][z]*y[13][z])*(y[10][z]-fsEk)
        -FSgl*(y[10][z]-fsEl)
        -gabaFS[z]
        -gabaFSX[z]
        -gelec*((FSCOUNT-1)*y[10][z]-(SumFS-y[10][z])) // all to all electrical connections
        -gTCi[z]*(y[10][z] - 0)
        +0.1;
 //       +IappFS[z];
        
        k[ 11 ][z] = 0.32*(y[10][z]+54)/(1-exp(-(y[10][z]+54)/4))*(1-y[11][z])-0.28*(y[10][z]+27)/(exp((y[10][z]+27)/5)-1)*y[11][z];     // derivative equation for m
        k[ 12 ][z] = 0.128*exp(-(y[10][z]+50)/18)*(1-y[12][z])-4/(1+exp(-(y[10][z]+27)/5))*y[12][z];   // derivative equation for h
        k[ 13 ][z] = 0.032*(y[10][z]+52)/(1-exp(-(y[10][z]+52)/5))*(1-y[13][z])-0.5*exp(-(y[10][z]+57)/40)*y[13][z];  // derivative equation for n
    }
    
    //************** FS2 Cells ************//
    for (z=0; z<FSCOUNT2; z++){
        
        if (time <= 20000){ //20000
            gelec = 0;
            gfs_fs = 0.03/10; // 0.03/10
            gfsX = 0.5/10; //0.25/10, changed from gfsX = 0.25/10;
        }
        
       
        IappFS2[z] = 0;
        
        gabaFS2[z] = gfs_fs/(FSCOUNT2-1.0)*(Sum19-y[19][z])*(y[15][z]-iEi);
        gabaFSX2[z] = gfsX/(FSCOUNT)*(Sum14)*(y[15][z]-iEi);
        
        k[ 15 ][z] =
        -FSgna*y[17][z]*(y[16][z]*y[16][z]*y[16][z])*(y[15][z]-fsEna)
        -FSgk*(y[18][z]*y[18][z]*y[18][z]*y[18][z])*(y[15][z]-fsEk)
        -FSgl*(y[15][z]-fsEl)
        -gabaFS2[z]
        -gabaFSX2[z]
        -gelec*((FSCOUNT2-1)*y[15][z]-(SumFS2-y[15][z])) // all to all electrical connections
        -gTCfs2[z]*(y[15][z] - 0)
         +0.1;
    //    +IappFS2[z];
        
        k[ 16 ][z] = 0.32*(y[15][z]+54)/(1-exp(-(y[15][z]+54)/4))*(1-y[16][z])-0.28*(y[15][z]+27)/(exp((y[15][z]+27)/5)-1)*y[16][z];     // derivative equation for m
        k[ 17 ][z] = 0.128*exp(-(y[15][z]+50)/18)*(1-y[17][z])-4/(1+exp(-(y[15][z]+27)/5))*y[17][z];   // derivative equation for hLee mill grinder
        k[ 18 ][z] = 0.032*(y[15][z]+52)/(1-exp(-(y[15][z]+52)/5))*(1-y[18][z])-0.5*exp(-(y[15][z]+57)/40)*y[18][z];  // derivative equation for n
    }
    
    
    
    /************** Input from Thalamus ************/
    double IappTC[TCCOUNT]; double TCnoise;
    
    //**** Delta input *****//
    k[ 36 ][0] = 0.3*sin(0.01*time); //0.1 is ~16 Hz
    
   //**** Monocular Deprivation *****//
    for (z=0; z<TCCOUNT; z++){
        if ((time >= 40000 && time <= END) && z < FSCOUNT){
    //    if (time > END){
            IappTC[z] = 0;
            TCnoise = 0;
        }
        else{
            IappTC[z] = 0.2; //0.2
            TCnoise = 50; //50
        }
        
        
  
        
        k[ 30 ][z] = -100*y[32][z]*(y[31][z]*y[31][z]*y[31][z])*(y[30][z]-50)
        -80*(y[33][z]*y[33][z]*y[33][z]*y[33][z])*(y[30][z]+100)
        -0.1*(y[30][z]+67)
        +TCnoise*sqrt(0.05)*box_muller(0,1)
    //    + delta
        + IappTC[z];
        
        ///* // rate equations
        k[ 31 ][z] = 0.32*(y[30][z]+54)/(1-exp(-(y[30][z]+54)/4))*(1-y[31][z])-0.28*(y[30][z]+27)/(exp((y[30][z]+27)/5)-1)*y[31][z];     // derivative equation for m
        k[ 32 ][z] = 0.128*exp(-(y[30][z]+50)/18)*(1-y[32][z])-4/(1+exp(-(y[30][z]+27)/5))*y[32][z];   // derivative equation for h
        k[ 33 ][z] = 0.032*(y[30][z]+52)/(1-exp(-(y[30][z]+52)/5))*(1-y[33][z])-0.5*exp(-(y[30][z]+57)/40)*y[33][z];  // derivative equation for n
        
        double sgTC[TCCOUNT]; double tauTC = 6;
        sgTC[z] = 5*(1+tanh(y[30][z]/4));
        k[ 34 ][z] = sgTC[z]*(1-y[34][z])-y[34][z]/tauTC;
        
        double sgTCi[TCCOUNT]; double tauTCi = 2;
        sgTCi[z] = 5*(1+tanh(y[30][z]/4));
        k[ 35 ][z] = sgTCi[z]*(1-y[35][z])-y[35][z]/tauTCi;
    }
    
}


double rnd()  {   // returns a uniform random number between 0 and 1
    return((double)rand()/RAND_MAX);
}


double randgauss( double min, double max, double sigma, double centre)
{
    double gauss = (min + (max-min) * (double)rand()/RAND_MAX); //create random domain between [min,max]
    
    //    double tmp = (random-centre)/sigma;
    //    double gauss= exp(-tmp*tmp/2); //gaussian formula
    
    return gauss;
}



float box_muller(float m, float s)    /* normal random variate generator */  //generates Gaussian random numbers
{                        /* mean m, standard deviation s */
    float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;
    
    if (use_last)                /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do {
            x1 = 2.0 * randgauss(0,1,1,0) - 1.0;
            x2 = 2.0 * randgauss(0,1,1,0) - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );
        
        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }
    
    return( m + y1 * s );
}































