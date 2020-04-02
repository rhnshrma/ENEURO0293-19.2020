/*
 current/I_H_HTC.hpp - Current flowing across neuronal membrane due to
                          HCN(hyperpolarization activated cyclic nucleotide gated) channels (mixed cationic current). (Destexhe et al. J. Neurophysiology 1996)

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

#ifndef INCLUDED_I_H_HTC_HPP
#define INCLUDED_I_H_HTC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_H_HTC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t,const unsigned index) {
    double gh = engine::neuron_value(index,"gh"), eh = -40;
                                                                        // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned r_index = engine::neuron_index(index, "r_H_HTC");


    double v = variables[v_index];
    double r = variables[r_index];

                                                                            // Calculating Time Constants
    double r_inf = 1/(exp((v+60)/5.5)+1);
    double tau_r = 20 + 1000/(exp((v +56.5)/14.2) + exp(-1*(v+74)/11.6)) ;

    if (r<0) variables[r_index] = 0;
    if (r>1) variables[r_index] = 1;  
    dxdt[r_index] = (r_inf - r)/tau_r;

     // Current
    engine::neuron_value(index, "I_H_HTC", (-1*(gh*r*(v - eh))));

  } // function current
}; // class I_T_RE-Cells

} // insilico

#endif
