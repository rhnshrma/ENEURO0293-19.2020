 /*
 current/I_H_TC.hpp - Current flowing across neuronal membrane due to
                          I_h Current conductance. (Destexhe, McCormick & Sejnowski 1993 Biophys J 65:2473–2477)

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

 Current that flows as the Calcium gated I_h current. (Destexhe, McCormick & Sejnowski 1993 Biophys J 65:2473–2477)
*/

#ifndef INCLUDED_I_H_TC_HPP
#define INCLUDED_I_H_TC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_H_TC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t,const unsigned index, double It) {
    double gh = 0.1, eh = -43;
                                                            // Declaring required constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned o_index = engine::neuron_index(index, "o_h_TC");
    unsigned c_index = engine::neuron_index(index, "c_h_TC");
    unsigned p_index = engine::neuron_index(index, "p_h_TC");
    unsigned Ca_index = engine::neuron_index(index, "Ca_conc");

    double v = variables[v_index];
    double o = variables[o_index];
    double p = variables[p_index];
    double c = variables[c_index];
    double Ca = variables[Ca_index];

 

                                                            // Calculating required time constants
    double hc_inf = 1.0/(exp((v + 75.0)/5.5) + 1.0);
    double tau_s = (20 + 1000.0/(exp((v+71.5)/14.2) + exp(-(v +89.0)/11.6)));

    double beta_c=(1.0 - hc_inf)/tau_s;
    double alpha_c=hc_inf/tau_s;

    if (o<0) variables[o_index] = 0;  
    if (c<0) variables[c_index] = 0;  
    if (p<0) variables[p_index] = 0;      


                                                            // ODE set
    dxdt[o_index] = 0.0001*(1.0 - c -o) - 0.001*((1.0-p)/0.01);
    dxdt[c_index] = (beta_c*o - alpha_c*c);
    dxdt[p_index] = 0.0004*(1.0-p) - 0.004*pow((Ca/0.0002),2);
 
  


                                                            // Current
    engine::neuron_value(index, "I_H_TC", (-1*gh*(o + 2*(1.0 - c -o))*(v - eh)));


  } // function current
}; // class I_H-current_thalamo-cortical cells

} // insilico

#endif


