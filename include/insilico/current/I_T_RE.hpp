/*
 current/I_T_RE.hpp - Current flowing across neuronal membrane due to
                          Calcium (Ca) conductance. (McCarthy et al., 2008)

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


*/

#ifndef INCLUDED_I_T_RE_HPP
#define INCLUDED_I_T_RE_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_T_RE {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, unsigned index) {
    double gca = 2.3, eca = 120;
                                                                            // Declaring Constants
    unsigned v_index  = engine::neuron_index(index, "v");
    unsigned m_index  = engine::neuron_index(index, "m_t_RE");
    unsigned h_index  = engine::neuron_index(index, "h_t_RE");
    unsigned Ca_index= engine::neuron_index(index, "Ca_conc");
    
    double v  = variables[v_index];
    double m  = variables[m_index];
    double h  = variables[h_index];
    double Ca = variables[Ca_index];
    double calcium_shift = engine::neuron_value(index,"calcium_shift2");

    double vt = v + 0;


                                                                            // Calculating Time Constants
    double m_inf = 1.0/(1.0 + exp(-(vt+52.0)/7.4));
    double tau_m = 0.999 + 0.333/(exp((vt+27)/10.0) + exp(-(vt+102.0)/15.0));
   
    double h_inf = 1.0/(1.0 + exp((vt+80.0)/5.0));
    double tau_h = 28.307 + 0.33/(exp((vt+48)/4.0) + exp(-(vt+407.0)/50.0));


    eca = 2.303*8.314*309150*log(2/Ca)/(2*96489);

    if (h<0) variables[h_index] = 0;  
    if (m<0) variables[m_index] = 0;  
    if (h>1) variables[h_index] = 1;  
    if (m>1) variables[m_index] = 1;
                                                                             // ODE set
    dxdt[m_index] = (m_inf - m)/tau_m;
    dxdt[h_index] = (h_inf - h)/tau_h;

    if((v-eca)<0) { dxdt[Ca_index] = 10*(-1*calcium_shift*(gca*pow(m,2)*h*(v - eca)))/(2*96489) + (0.00024 - Ca)/3.0; }
    else { dxdt[Ca_index] = (0.00024 - Ca)/3.0; }

    // Current
    engine::neuron_value(index, "I_T_RE", (calcium_shift*(-1*(gca*pow(m,2)*h*(v - eca)))));

  } // function current
}; // class I_T_RE-Cells

} // insilico

#endif
