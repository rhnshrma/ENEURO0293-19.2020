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

#ifndef INCLUDED_N_HTCA1_HPP
#define INCLUDED_N_HTCA1_HPP

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

class N_HTC_mAch_1: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned v_pre_index = engine::neuron_index(1,"v");
    
    double v = variables[v_index];


    double I_Syn_GJ = 0;
    I_Syn_GJ = engine::neuron_value(index,"g_gj")*(v-variables[v_pre_index]); 
    
    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
    std::vector<double> gsyn = engine::get_pre_neuron_values(index, "gsyn");
    std::vector<double> T = engine::get_pre_neuron_values(index, "T");
    double I_Syn = 0;
    double I_Syn_I = 0;
    double I_Syn_E = 0;
    std::vector<double> esyn = engine::get_pre_neuron_values(index, "esyn");
    //std::vector<double> v_pre = engine::get_pre_neuron_values(index, "v_pre");

    for(std::vector<unsigned>::size_type iterator = 0; iterator < g1_indices.size(); ++iterator) {
      double I = variables[g1_indices.at(iterator)]*gsyn.at(iterator)*(v - esyn.at(iterator));
      I_Syn = I_Syn + I;
      if (I>0){I_Syn_I = I_Syn_I + I; }
      else {I_Syn_E = I_Syn_E + I; }  
    }

    
    

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


    double T_sum = 0;
    for(std::vector<unsigned>::size_type iterator = 0; iterator < T.size(); ++iterator) {
       T_sum = T_sum + T[iterator];
     }
    double T_avg = T_sum/T.size();
    engine::neuron_value(index,"T",T_avg);

/*    double T_last_epsp = engine::neuron_value(index,"T_last_epsp");
    double T_next_epsp = engine::neuron_value(index,"T_next_epsp");
    double T_last_ipsp = engine::neuron_value(index,"T_last_ipsp");
    double T_next_ipsp = engine::neuron_value(index,"T_next_ipsp");
    double gs = engine::neuron_value(index,"gs");
    
    if((t-T_next_epsp) > 0)
    {
        double z = random2::rand<double>();
        double next_T_epsp = -1*(log(1-z)*engine::neuron_value(index,"t_epsp")); //Converting Uniform Ditribution to Poisson(http://stackoverflow.com/questions/2106503/pseudorandom-number-generator-exponential-distribution)
        engine::neuron_value(index,"T_last_epsp",T_next_epsp);
        engine::neuron_value(index,"T_next_epsp",(T_next_epsp+next_T_epsp));

     }

    double I_epsp = ((t < (T_last_epsp + 2)) ? (gs*exp(T_last_epsp - t)*v):0 );

        if((t-T_next_ipsp) > 0)
    {
        double z = random2::rand<double>();
        double next_T_ipsp = -1*(log(1-z)*engine::neuron_value(index,"t_ipsp")); //Converting Uniform Ditribution to Poisson(http://stackoverflow.com/questions/2106503/pseudorandom-number-generator-exponential-distribution)
        engine::neuron_value(index,"T_last_ipsp",T_next_ipsp);
        engine::neuron_value(index,"T_next_ipsp",(T_next_epsp+next_T_ipsp));

     }

    double I_ipsp = ((t < (T_last_ipsp + 2)) ? (1.0*gs*exp(T_last_ipsp - t)*(v+85) ):0 );*/

    engine::neuron_value(index,"I_Syn_I",I_Syn_I);
    engine::neuron_value(index,"I_Syn_E",I_Syn_E);

    double I_Na = engine::neuron_value(index,"I_Na_TC");
    double I_K = engine::neuron_value(index,"I_K_TC");
    double I_HTC_Leak = engine::neuron_value(index,"I_Leak_HTC");
    double I_TLT = engine::neuron_value(index,"I_T_TC");
    double I_H = engine::neuron_value(index,"I_H_HTC"); 
    double I_THT =  engine::neuron_value(index,"I_THT_HTC");
    double I_AHP = engine::neuron_value(index,"I_AHP_HTC");
    double I_Ext = engine::neuron_value(index, "iext");
    engine::neuron_value(index,"I_Syn_I",I_Syn_I);
    engine::neuron_value(index,"I_Syn_E",I_Syn_E);


    // ODE
    dxdt[v_index] = I_TLT + I_Na + I_K + I_HTC_Leak + I_H + I_THT + I_AHP + I_Ext - I_Syn - I_Syn_GJ /*- I_epsp - I_ipsp*/;
  }
};

} // insilico

#endif

