/*
 synapse/S_GABAA.hpp - Default non-specialized synapse with last spike logic


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

#ifndef INCLUDED_S_GABAA_HPP
#define INCLUDED_S_GABAA_HPP

#include "core/engine.hpp"

namespace insilico {

class S_GABAA: public Synapse {
 public:
  void ode_set(state_type& variables, state_type& dxdt, const double t, unsigned index) {
    
    unsigned g1_index = engine::synapse_index(index, "g1");

    double g1 = variables[g1_index];


    // synapse logic for delay for recently spiked neuron


    unsigned neuron_index = engine::synapse_value(index, "pre");
    unsigned v_index = engine::neuron_index(neuron_index,"v");
    double v_pre = variables[v_index];

 
    // constants from file
    double tau = engine::synapse_value(index, "tau");

    // ODE set
    dxdt[g1_index] = (2*(1 + tanh(v_pre/4))*(1-g1) -g1/tau);


  } // function ode_set
}; // class S_GABAA

} // insilico

#endif


