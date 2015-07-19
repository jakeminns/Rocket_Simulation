//
//  Rocket Class.h
//  Rocket_Simulation
//
//  Created by Jake Minns on 12/07/2015.
//  Copyright (c) 2015 Jake Minns. All rights reserved.
//

#ifndef __Rocket_Simulation__Rocket_Class__
#define __Rocket_Simulation__Rocket_Class__
#include <vector>
#include <stdio.h>
using namespace std;
#endif /* defined(__Rocket_Simulation__Rocket_Class__) */

class Rocket{
    
public:

    Rocket();
    
    double payload_mass;
    double structural_mass;
    double propellant_mass;
    
    double delta_v(double time, int  var, vector<double> velocity, double gravitatonal_acc, double rocket_mass);
    double delta_mass(double rocket_mass,double time);
    
protected:
    
private:
    string rocket_name;
    double isp;
    double thrust;
    double burn_time;
    
};