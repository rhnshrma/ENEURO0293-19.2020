/*
 neuron/N_TC.hpp - Hodgkin-Huxley Squid Axon experiment (Hodgkin-Huxley, 1952)

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

#ifndef INCLUDED_N_TC_HPP
#define INCLUDED_N_TC_HPP

#include "insilico/core/engine.hpp"

#include "../current/I_Na_TC.hpp"
#include "../current/I_K_TC.hpp"
#include "../current/I_Leak_TC.hpp"
#include "../current/I_H_TC.hpp"
#include "../current/I_T_TC.hpp"
#include <random>

namespace insilico {

    class random2 {
 public:
  template<class T>
  static T rand() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(0,1);
    return dist(gen);
  }
};

class N_TC : public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {
    
    unsigned v_index = engine::neuron_index(index, "v");

        double v = variables[v_index];
    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
    std::vector<double> gsyn = engine::get_pre_neuron_values(index, "gsyn");
    double I_Syn = 0;
    std::vector<double> esyn = engine::get_pre_neuron_values(index, "esyn");

    for(std::vector<unsigned>::size_type iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn = I_Syn + variables[g1_indices.at(iterator)]*gsyn.at(iterator)*(v - esyn.at(iterator)); 
    }

    double thresh = 0, dead_time=1;
     double last_spike = engine::neuron_value(index,"last_spike");

    // associated delay for next spikes
    if((v > thresh) && (t - last_spike) > dead_time){
      engine::neuron_value(index, "last_spike", t);
    }

     I_Na_TC::current(variables, dxdt, t, index);
     I_K_TC::current(variables, dxdt, t, index);
     I_T_TC::current(variables, dxdt, t, index);
     double I_TC_T = engine::neuron_value(index, "I_T_TC");

     I_H_TC::current(variables, dxdt, t, index, I_TC_T);
     I_Leak_TC::current(variables, dxdt, t, index);

    double I_Na = engine::neuron_value(index, "I_Na_TC");
    double I_K = engine::neuron_value(index, "I_K_TC");
    double I_Leak = engine::neuron_value(index, "I_Leak_TC");
    double I_H = engine::neuron_value(index, "I_H_TC");

    double I_Ext = engine::neuron_value(index, "iext");
    double T_last = engine::neuron_value(index,"T_last");
    double T_next = engine::neuron_value(index,"T_next");
    double gs = engine::neuron_value(index,"gs");


    if((t-T_next) > 0)
    {
        double z = random2::rand<double>();
        double next_T = -1*(log(1-z)/10);
        engine::neuron_value(index,"T_last",T_next);
        engine::neuron_value(index,"T_next",(T_next+next_T));

     }

    double I_epsp = ((t < (T_last + 2)) ? (gs*exp(T_last - t)*v):0 );

    //std::cout<<"I_Na "<< I_Na << " I_K " << I_K << " I_H " << I_H << " I_TC_T " <<  I_TC_T<< " I_Leak " <<  I_Leak << " T_last " <<  T_last << " exp(T_last - t) " << exp(T_last - t) << " I_epsp " <<  I_epsp  << std::endl;

    // ODE
    dxdt[v_index] =  I_Na + I_K + I_H + I_TC_T + I_Leak + I_Ext - I_Syn - I_epsp;
  }
};

} // insilico

#endif
