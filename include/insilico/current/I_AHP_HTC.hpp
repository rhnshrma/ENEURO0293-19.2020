/*
 current/I_AHP_HTC.hpp - Current flowing across neuronal membrane due to
                          Calcium (Ca) activated Potassium conductance. (McCarthy et al., 2008)



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

#ifndef INCLUDED_I_AHP_HTC_HPP
#define INCLUDED_I_AHP_HTC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_AHP_HTC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, unsigned index) {
    double gahp = 15, eahp = -100;
                                                                        // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned m_index = engine::neuron_index(index, "m_ahp_HTC");
    unsigned Ca_index= engine::neuron_index(index, "Ca_conc");


    double v = variables[v_index];
    double m = variables[m_index];
    double Ca = variables[Ca_index];

                                                                            // Calculating Time Constants
    double m_inf = 48.0*pow(Ca,2)/(0.09 + 48*pow(Ca,2));
    double tau_m = 1/(0.09 + 48*pow(Ca,2));

    if (m>1) variables[m_index] = 1;          
    if (m<0) variables[m_index] = 0;  
    //ODE-Set
    dxdt[m_index]= (m_inf - m)/tau_m;

     // Current
    engine::neuron_value(index, "I_AHP_HTC", (-1*(gahp*pow(m,2)*(v - eahp))));

  } // function current
}; // class I_T_RE-Cells

} // insilico

#endif
