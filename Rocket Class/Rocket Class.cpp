//
//  Rocket Class.cpp
//  Rocket_Simulation
//
//  Created by Jake Minns on 12/07/2015.
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

Rocket::Rocket(){
    
    cout << "Insert Name of Rocket:"<<endl;
    cin >> rocket_name;
    
    cout << "The following defines the specifications for " << rocket_name << " rocket: "<< endl;
   
    cout << "Input Specific Impulse of Rocket (ISP):" << endl;
    cin >> isp;
    
    cout << "Input Payload Mass:" << endl;
    cin >> payload_mass;
    
    cout << "Input Structural Mass:" << endl;
    cin >> structural_mass;
    
    cout << "Input Propellant Mass:" << endl;
    cin >> propellant_mass;
    
    cout << "Input Rocket Thrust:" << endl;
    cin >> thrust;

}



double Rocket::delta_mass(double rocket_mass, double time)
{
 if(rocket_mass>=payload_mass+structural_mass){
    
     double d_m = thrust/(isp*9.81);
    
     return d_m;
 
 }else
    {
     return 0;
    }
}




double Rocket::delta_v(double time, int var, vector<double> velocity, double gravitatonal_acc, double rocket_mass)
{
    if(rocket_mass>=payload_mass+structural_mass){
        
        double unit_vector = velocity[var]/(sqrt(pow(velocity[0],2)+pow(velocity[1],2)+pow(velocity[2],2)));
        
        return unit_vector*(isp*(9.81))*log((payload_mass+structural_mass+propellant_mass)/(payload_mass+structural_mass));
        
    }else{
        return 0;
    }
    
    return 0;
    
}