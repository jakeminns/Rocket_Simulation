//
//  main.cpp
//  Rocket_Simulation
//
//  Created by Jake Minns on 08/07/2015.
//  Copyright (c) 2015 Jake Minns. All rights reserved.
//


#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Rocket Class.h"


using namespace std;

//Function Prototypes/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double delta_position(double velocity, double acceleration, double time);
double gravitational_acc(vector<double> position, double height, int var);
double drag_acc(double velocity_comp, double height, double rocket_mass);
vector<double> velocity_vec(double theta, double beta, double v_mag);
void plotResult(vector<double> xData, vector<double> yData, vector<double> zData, int dataSize);


int main() {
    
//Variable Definitions/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    vector<double> launch_position(3),position(3), tracking_position(3); // launch position, position with respect to centre of earth, position with repect to launch position.
    vector<double> velocity(3), velocity_angle(2), grav_acc(3);
    vector<double> xplot(0), yplot(0), zplot(0);
    
    double velocity_mag; // Velocity Vector Magnitude
    double position_mag; // Position Vector Magnitude, position with respect to centre of earth
    
    double const earth_radius=6371000; // Radius of earth's magnitude

    double time_iter = 0.05; // Time iteraration of mechanics calculation
    double time =0; // Time of flight

    double dv_dt; // Change in velocity with respect to time (accelaration)
    double rocket_mass; // total mass of rocket, varies with time
    
    enum co_ordinate{x,y,z}; // Enumarated list of co ordinates, x=0, y=1, z=2
    enum angle{theta,beta}; // Enumarated list of angle variables, theta=0, beta=1
    
    int iter = 0; // variable used to cycle through position vectors stored
    
//// MAKE INITIAL VELOCITY VECTOR SAME AS INITIAL LAUNCH UNIT VECTOR
//Rocket Definition/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Rocket rocket;
    
   rocket_mass = rocket.payload_mass + rocket.structural_mass + rocket.propellant_mass; // define rocket initial total mass

    // Test Imputs
    //Rocket rocket_name; // Rocket type Falcon 1 defined
    //Falc_1.isp = 100; // Falcon 1 specific impulse
   // Falc_1.payload_mass = 400; // Falcon 1 payload mass
   // Falc_1.structural_mass = 5000; // Falcon 1 structural mass
    //Falc_1.propellant_mass = 38555; // Falcon 1 propellent mass
    //Falc_1.thrust = 654000; // Falcon 1 thrust
    
//Imputs/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    velocity_mag = 1; //Initial Velocity magnitude
    
    cout << "Input launch angle between x & y plane:"<<endl;
    cin >> velocity_angle[beta];  // angle between x & y
    
    cout << "Input launch angle between z & x-y plane:"<<endl;
    cin >> velocity_angle[theta];  // angle between z & x-y plane
    
    velocity = velocity_vec(velocity_angle[theta], velocity_angle[beta], velocity_mag); //calculates initial velocity vector componants from angle and magnitude
    launch_position[z] = 6371000; // Launch Position z componant with respect to centre of earth

    
//Inisalising Variables/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    position[z] = launch_position[z];
    position_mag = sqrt(pow(position[x],2)+pow(position[y],2)+pow(position[z],2));
    xplot.push_back(position[x] - launch_position[x]); // Initial position added to xplot vector list of positions throughout flight
    yplot.push_back(position[y] - launch_position[y]); // Initial position added to yplot vector list of positions throughout flight
    zplot.push_back(position[z] - launch_position[z]); // Initial position added to zplot vector list of positions throughout flight
    

//Trajectory Calculation/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    while(position_mag>=earth_radius){ //check that projectile position does not pass through earth
   
        position_mag = sqrt(pow(position[x],2)+pow(position[y],2)+pow(position[z],2));
        
        
        for(int var = 0 ; var<3;var++){ //loop beteen x, y & z componants
            
            dv_dt = gravitational_acc(position,position[var],var)+drag_acc(velocity[var], tracking_position[var], rocket_mass)+rocket.delta_v(time,var, velocity,gravitational_acc(position,position[var],var), rocket_mass); // Addition of vectors contributing to change in velocity (acceleration)

            velocity[var] = velocity[var] + dv_dt*time_iter; // V=u+at change in velocity
            position[var] = position[var] + velocity[var]*time_iter; // Change in position according to velocity over time iteration
            tracking_position[var] = position[var] - launch_position[var]; // Updates trackng position variable, position with respect to launch position
            
            
        }
        
        rocket_mass = rocket_mass - rocket.delta_mass(rocket_mass,time); // Updates rocket mass according to fuel burnt over time iteration
        
        time += time_iter; //updates time since launch, adds time iteration to time variable
        
        xplot.push_back(tracking_position[x]);
        yplot.push_back(tracking_position[y]);
        zplot.push_back(tracking_position[z]);
        
        cout << position[0] << "   " << position[1] << "   " <<tracking_position[2]<<endl;

        iter++;
    }
    
//Plot Results/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    xplot[iter-1] =0; // remove last position entry
    yplot[iter-1] =0;
    zplot[iter-1] =0;

    
    plotResult(xplot, yplot, zplot, iter+1); // plot results
    

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<double> velocity_vec(double theta, double beta, double v_mag){ //Function to calculate velocity vector given angle and magnitude
    
    vector<double> velocity(3); // define velocity vector
    
    velocity[0] = (v_mag * cos(theta)) * cos(beta); //X velocity componant
    velocity[1] = (v_mag * cos(theta)) * sin(beta); //Y velocity componant
    velocity[2] =  v_mag * sin(theta); //Z velocity componant
    
    //cout << "velocity:   "<< "    x:  " << velocity[0] << "    y:  " << velocity[1] << "    z:  " << velocity[2] << endl;
    
    
    return velocity; // returns velocity vector
}


// For converting gravitational acceleration to consider round earth consider height as absolute position magnitude & 0,0 position is centre of earth
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double gravitational_acc(vector<double> position, double height, int var){
   
    vector<double> gravity(3);
    
    if(height == 0){ // check that height does not equal 0, avoid division of 0
       
        gravity[var] = 0;
   
    }else{

    
    double gravitational_const = 6.67384*pow(10,-11); //Define gravitational constant
    double earth_mass = 5.972*pow(10,24); // Define mass of earth
    
    double position_mag = (sqrt(pow(position[0],2)+pow(position[1],2)+pow(position[2],2))); // magnitude of position vector
    double unit_vec = position[var]/position_mag; // Define position unit vector
        
    gravity[var] = -(gravitational_const * earth_mass*unit_vec)/(pow(position_mag,2)); // Gravitational acceleration vector, = -G*M*r(unit_vector)/r^2

    }
    
    
    return gravity[var]; // Return vector gravity

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double drag_acc(double velocity_comp, double height, double rocket_mass){

    double area = 0.2; // Surface area cross section
    double drag_coef = 0.5; // Drag coefficent
    double mass_density; // air mass density
    
    
    mass_density=1.2-0.0001333333333*height; // Atmospheric density model
    
    if (mass_density<0){ // check that mass density does not drop below 0
        
        mass_density =0;
        
    }else{
    }
    
    double drag= -(0.5*mass_density*pow(velocity_comp,2)*drag_coef*area)/rocket_mass; // Delta V due to drag =  F/m = (-0.5 * p * v^2 * Cd * a)/m
    
 
    return drag;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void plotResult(vector<double> xData, vector<double> yData, vector<double> zData, int dataSize){
    
    FILE *gnuplotPipe, *tempDataFile;
    char const *tempDataFileName;
    
    double x,y,z; // x,y,z co ordinate variables

    tempDataFileName = "t"; // define temp file name
    gnuplotPipe = popen("gnuplot","w"); //
    
    if(gnuplotPipe){
        fprintf(gnuplotPipe, "splot \"%s\"\n", tempDataFileName);
        fflush(gnuplotPipe);
        tempDataFile = fopen(tempDataFileName,"w");
        
        for (int i=0; i<= dataSize; i++){
            x=xData[i];
            y=yData[i];
            z=zData[i];
            fprintf(tempDataFile,"%lf %lf %lf\n",x,y,z);
        }
        
        fclose(tempDataFile);
        printf("press enter to continue");
        
        getchar();
        remove(tempDataFileName);
        fprintf(gnuplotPipe,"exit \n");
    } else {
        printf("gnuplot not found...");
    }
}


