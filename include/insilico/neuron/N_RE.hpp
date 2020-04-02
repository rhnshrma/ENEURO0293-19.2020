/*
 neuron/N_rohan.hpp - Hodgkin-Huxley Squid Axon experiment (Hodgkin-Huxley, 1952)

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

#ifndef INCLUDED_N_RE_HPP
#define INCLUDED_N_RE_HPP

#include "insilico/core/engine.hpp"

#include "../current/I_Na_RE.hpp"
#include "../current/I_K_RE.hpp"
#include "../current/I_Leak_RE.hpp"
#include "../current/I_T_RE.hpp"

namespace insilico {

class N_RE: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {
    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
    std::vector<double> gsyn = engine::get_pre_neuron_values(index, "gsyn");
    double I_Syn = 0;
    std::vector<double> esyn = engine::get_pre_neuron_values(index, "esyn");

    unsigned v_index = engine::neuron_index(index, "v");
    double v = variables[v_index];

    // note the spike
    double thresh = 0, dead_time=1;
    double last_spike = engine::neuron_value(index,"last_spike");

    // associated delay for next spikes
    if((v > thresh) && (t - last_spike) > dead_time){
      engine::neuron_value(index, "last_spike", t);
    }

    // incoming synaptic currents



    for(unsigned iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn = I_Syn + variables[g1_indices.at(iterator)]*gsyn.at(iterator)*(v - esyn.at(iterator));
    }

     I_Na_RE::current(variables, dxdt, t, index);
     I_K_RE::current(variables, dxdt, t, index);
     I_T_RE::current(variables, dxdt, t, index);
     I_Leak_RE::current(variables, dxdt, t, index);


    double I_Na = engine::neuron_value(index, "I_Na_RE");
    double I_K = engine::neuron_value(index, "I_K_RE");
    double I_RE_T = engine::neuron_value(index, "I_T_RE");
    double I_RE_Leak = engine::neuron_value(index, "I_Leak_RE");
    double I_Ext = engine::neuron_value(index, "iext");

    engine::neuron_value(index,"I_Syn",I_Syn);

    
    // ODE
    dxdt[v_index] = I_Na + I_K + I_RE_Leak + I_RE_T + I_Ext - I_Syn;
  }
};

} // insilico

#endif
