/*
 current/I_Na_HH1952.hpp - Current flowing across neuronal membrane due to
                           Sodium (Na) conductance. (Hodgkin-Huxley, 1952)

 Copyright (C) 2014 Pranav Kulkarni, Collins Assisi Lab, IISER, Pune <pranavcode@gmail.com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 Brief:

 Current that flows through Sodium (Na) channel due to the potential difference
 caused by Sodium (Na) conductance across neuronal membrane. (Hodgkin-Huxley, 1952)
*/

#ifndef INCLUDED_I_NA_PY_HPP
#define INCLUDED_I_NA_PY_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_Na_PY {
 public:
 static void current(const state_type &variables, state_type &dxdt, const double t, unsigned index) {    
    double gna = 50, ena = 100;
                                                            // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned m_index = engine::neuron_index(index, "m_Na_PY");
    unsigned h_index = engine::neuron_index(index, "h_Na_PY");
    
    double v = variables[v_index];
    double m = variables[m_index];
    double h = variables[h_index];
    

                                                            // Calculating Time Constants
    double alpha_m = 0.32*vtrap(-1*(54.0+v),4.0);
    double beta_m  = 0.28*vtrap((v+27.0),5.0);
   
    double alpha_h = 0.128*exp(-(v+50.0)/18.0);
    double beta_h  = 4.0/(exp(-(27.0 + v)/5)+1);


    double m_inf=alpha_m/(alpha_m+beta_m);
    double tau_m=1/(alpha_m + beta_m);

    double h_inf=alpha_h/(alpha_h+beta_h);
    double tau_h=1/(alpha_h + beta_h);

    if (m<0) variables[m_index] = 0;  
    if (h<0) variables[h_index] = 0;   
                                                            // ODE set
    dxdt[m_index]=(m_inf-m)/tau_m;
    dxdt[h_index]=(h_inf-h)/tau_h;



                                                            // Current
     engine::neuron_value(index,"I_Na_PY",(-1*(gna*pow(m, 3)*h*(v - ena))));

  } // function current
}; // class I_Na_HH1952

} // insilico

#endif

