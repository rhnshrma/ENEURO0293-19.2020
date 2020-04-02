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

#ifndef INCLUDED_N_HTCA_4cell_HPP
#define INCLUDED_N_HTCA_4cell_HPP

#include "insilico/core/engine.hpp"

#include "insilico/current/I_Na_TC.hpp"
#include "insilico/current/I_K_TC.hpp"
#include "insilico/current/I_Leak_HTC.hpp"
#include "insilico/current/I_T_TC.hpp"
#include "insilico/current/I_THT_HTC.hpp"
#include "insilico/current/I_AHP_HTC.hpp"
#include "insilico/current/I_H_HTC.hpp"

namespace insilico {

    class random5 {
 public:
  template<class T>
  static T rand() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(0,0.01);
    return dist(gen);
  }
};

class N_HTC_mAch_4cell: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {
    unsigned v_index = engine::neuron_index(index, "v");
    //unsigned v_pre_index = engine::neuron_index(1,"v");
    
    double v = variables[v_index];
    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
    std::vector<double> gsyn = engine::get_pre_neuron_values(index, "gsyn");
    double I_Syn = 0;
    double I_Syn_GJ = 0;
    std::vector<double> esyn = engine::get_pre_neuron_values(index, "esyn");
    std::vector<double> v_pre = engine::get_pre_neuron_values(index, "v_pre");

    for(std::vector<unsigned>::size_type iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn = I_Syn + variables[g1_indices.at(iterator)]*gsyn.at(iterator)*(v - esyn.at(iterator)); 
    }

     I_Syn_GJ = engine::neuron_value(index,"g_gj")*(v-variables[v_pre_index]); 
    

    //if (t>10){ engine::neuron_value(index,"iext",0);}

    I_Na_TC::current(variables, dxdt, t, index);
    I_K_TC::current(variables, dxdt, t, index);
    I_THT_HTC::current(variables, dxdt, t, index);
    I_AHP_HTC::current(variables, dxdt, t, index);
    I_H_HTC::current(variables, dxdt, t, index);
    I_T_TC::current(variables, dxdt, t, index);
    I_Leak_HTC::current(variables, dxdt, t, index);
    
    double thresh = 0, dead_time=1;
    double last_spike = engine::neuron_value(index,"last_spike");
     
    // associated delay for next spikes
    if((v > thresh) && (t - last_spike) > dead_time){
      engine::neuron_value(index, "last_spike", t);
    }


    double I_Na = engine::neuron_value(index,"I_Na_TC");
    double I_K = engine::neuron_value(index,"I_K_TC");
    double I_HTC_Leak = engine::neuron_value(index,"I_Leak_HTC");
    double I_TLT = engine::neuron_value(index,"I_T_TC");
    double I_H = 0;//engine::neuron_value(index,"I_H_HTC"); 
    double I_THT =  engine::neuron_value(index,"I_THT_HTC");
    double I_AHP = engine::neuron_value(index,"I_AHP_HTC");
    double I_Ext = engine::neuron_value(index, "iext");
    engine::neuron_value(index,"I_Syn_GJ",I_Syn_GJ);
    //engine::neuron_value(index,"v_p",variables[v_pre_index]);

    I_Ext=0;
    // ODE
    dxdt[v_index] = I_TLT + I_Na + I_K + I_HTC_Leak + I_H + I_THT + I_AHP + I_Ext - I_Syn - I_Syn_GJ;
  }
};

} // insilico

#endif
