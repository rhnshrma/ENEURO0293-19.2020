/*
 current/I_T_TC.hpp - Current flowing across neuronal membrane due to
                          Calcium (Ca) conductance. (Destexhe et al. J. Neurophysiology 1996

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

 Current that flows through Potassium (K) channel due to the potential difference
 caused by Potassium (K) conductance across neuronal membrane. (Hodgkin-Huxley, 1952)
*/

#ifndef INCLUDED_I_T_TC_HPP
#define INCLUDED_I_T_TC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_T_TC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, unsigned index) {
    double gca = 2, eca = 120; 
                                                                                     // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned h_index = engine::neuron_index(index, "h_t_TC");
    unsigned Ca_index= engine::neuron_index(index, "Ca_conc");
    
    double v = variables[v_index];
    double h = variables[h_index];
    double Ca = variables[Ca_index];
    double vt = v + 2;
    double calcium_shift = engine::neuron_value(index,"calcium_shift2");
                                                                                            // Calculating time constants
    double h_inf = 1.0/(exp((vt + 81)/4.0) + 1.0);
    double tau_h = (30.8 + ((211.4 + exp((vt+113.2)/5.0))/(1 + exp((vt+84.0)/3.2))) )/3.74;
    eca = 8.314*309150*log(2/Ca)/(2*96489);   
/*   if(vt < -80) {
      tau_h = exp((vt+467)/66.6) / 3.74 ;
   } 
   else {
     tau_h = ( 28 + exp(-(vt+22)/10.5) ) / 3.74 ;
   }
*/
    eca = 2.303*8.314*309150*log10(2/Ca)/(2*96489);                                     // Calculate     
    double m_inf = 1.0/(1.0 + exp(-(vt+57)/6.2));
    double It = calcium_shift*(-1*(gca*pow(m_inf,2)*h*(v - eca)));


    if (h<0) variables[h_index] = 0;  
    double swtch = engine::neuron_value(index,"Ca_slow_switch");
    double eq_Ca_fast=engine::neuron_value(index,"eq_Ca_fast");
    double eq_Ca_slow=engine::neuron_value(index,"eq_Ca_slow");
    double tau_fast=engine::neuron_value(index,"tau_fast");
    double tau_slow=engine::neuron_value(index,"tau_slow");                                                                                          // ODE set
    if(It>0) { dxdt[Ca_index] = 10*(It+(engine::neuron_value(index,"I_THT_HTC")))/(2*96489) + (eq_Ca_fast - Ca)/tau_fast +swtch*(eq_Ca_slow - Ca)/tau_slow; }
    else { dxdt[Ca_index] = (eq_Ca_fast - Ca)/tau_fast +(eq_Ca_slow - Ca)/tau_slow; }
    dxdt[h_index] = (h_inf - h)/tau_h;





    // Current
    engine::neuron_value(index, "I_T_TC", ((-1*calcium_shift*(gca*pow(m_inf,2)*h*(v - eca)))));

  } // function current
}; // class I_T-current_thalamo-cortical cells

} // insilico

#endif


