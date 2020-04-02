/*
 synapse/S_DefaultSynapse.hpp - Default non-specialized synapse with last spike logic

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

#ifndef INCLUDED_S_GABAB_HPP
#define INCLUDED_S_GABAB_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class S_GABAB: public Synapse {
 public:
  void ode_set(state_type& variables, state_type& dxdt, const double t, unsigned index) {
    
    unsigned R_index = engine::synapse_index(index, "R");
    unsigned G_index  = engine::synapse_index(index,"G");
    unsigned pre_neuron_index = engine::synapse_value(index, "pre");




    double R = variables[R_index];
    double G = variables[G_index];
    double last_spike = engine::neuron_value(pre_neuron_index,"last_spike");
    double delay = engine::synapse_value(index,"delay");
    double spike_duration = engine::synapse_value(index,"spike_duration");
    double T=0;

    // synapse logic for delay for recently spiked neuron


    if(((t - last_spike - delay) < (spike_duration))  &&  ((t-last_spike - delay) > 0)){
      T =  0.5;
        }
    engine::synapse_value(index,"T",T);
  

 
    // constants from file

    engine::synapse_value(index, "g1",(pow(G,4)/(100 + pow(G,4))));

    // ODE set
    dxdt[R_index] = (0.09*T*(1-R) - 0.0012*R);
    dxdt[G_index] = (0.18*R - 0.034*G);


  } // function ode_set
}; // class S_DefaultSynapse

} // insilico

#endif


