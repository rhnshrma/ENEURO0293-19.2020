/*
 current/I_THT_HTC.hpp - Current flowing across neuronal membrane due to
                          Calcium (Ca) conductance. (Vijayan et al. PNAS 2012)

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

#ifndef INCLUDED_I_THT_HTC_HPP
#define INCLUDED_I_THT_HTC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_THT_HTC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, unsigned index) {
    double gca = engine::neuron_value(index,"g_THT"), eca = 120;
                                                                        // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned h_index = engine::neuron_index(index, "h_tht_HTC");
    unsigned Ca_index= engine::neuron_index(index, "Ca_conc");

    double v = variables[v_index];
    double h = variables[h_index];
    double calcium_shift = engine::neuron_value(index,"calcium_shift1");
    double Ca = variables[Ca_index];


                                                                            // Calculating Time Constants
    double h_inf = 1.0/(1.0 + exp((v+62.2)/5.5));
    double tau_h = 0.6*(0.1483*exp(-0.09398*v) + 5.284*exp(0.008855*v));

    double m_inf = 1/(1 + exp(-1*(v+40.1)/3.5));


    eca = 2.303*8.314*309150*log10(2/Ca)/(2*96485);

    double It = calcium_shift*(-1*(gca*pow(m_inf,2)*h*(v - eca)));
    

    if (h<0) variables[h_index] = 0;  
                                                                             // ODE set
   /* if(It>0) { dxdt[Ca_index] = 10*It/(2*96485) + (0.00024 - Ca)/5.0; }
    else { dxdt[Ca_index] = (0.00024 - Ca)/5.0; }*/
    dxdt[h_index] = (h_inf - h)/tau_h;

                                                                            // Current
    // Current
    engine::neuron_value(index, "I_THT_HTC", It);

  } // function current
}; // class I_T_RE-Cells

} // insilico

#endif
