/*
 neuron/N_PY.hpp - Hodgkin-Huxley Squid Axon experiment (Hodgkin-Huxley, 1952)

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

#ifndef INCLUDED_N_PY_HPP
#define INCLUDED_N_PY_HPP

#include "insilico/core/engine.hpp"

#include "I_Na_PY.hpp"
#include "I_K_PY.hpp"
#include "I_Leak_PY.hpp"



namespace insilico {

/*class random {
 public:
  template<class T>
  static T rand() {
    random_device rd;
    mt19937_64 gen(rd());
    normal_distribution<double> dist(0,0.1);
    return dist(gen);
  } // function rand
};

class random2 {
 public:
  template<class T>
  static T rand() {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> dist(0,1);
    return dist(gen);
  } // function rand
 
};*/ // class random

class N_PY : public Neuron {



 public:
    void ode_set(state_type &variables, state_type &dxdt, const double t, unsigned index) {
   

    

    unsigned v_index = engine::neuron_index(index, "v");
        
    double v = variables[v_index];
    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
    std::vector<double> gsyn = engine::get_pre_neuron_values(index, "gsyn");
    double I_Syn = 0;
    std::vector<double> esyn = engine::get_pre_neuron_values(index, "esyn");

    //double T_last = engine::neuron_value(index,"T_last");
    

    for(std::vector<unsigned>::size_type iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn = I_Syn + variables[g1_indices.at(iterator)]*gsyn.at(iterator)*(v - esyn.at(iterator)); 
    }
    
    //ODE_Set
    I_Na_PY::current(variables, dxdt, t, index);
    I_K_PY::current(variables, dxdt, t, index);
    I_Leak_PY::current(variables, dxdt, t, index);

    double I_Na_PY = engine::neuron_value(index, "I_Na_PY");
    double I_K_PY = engine::neuron_value(index, "I_K_PY");
    double I_Leak_PY = engine::neuron_value(index, "I_Leak_PY");
    double I_Ext = engine::neuron_value(index, "iext");

/*    if(T_next<t)
    {
        T_last = T_next;
        T_next = T_last  - 0.1*log(random2::rand<double>());
        engine::neuron_value(index,"T_next",T_next);
        engine::neuron_value(index,"T_last",T_last);
   
    }


    I_Ext = I_Ext + ran;
    double I_poission = 0.02*exp(-(t-T_last)/2);
    
    I_Ext = I_Ext + I_poission;*/
 //std::cout << " I_Na_PY " << I_Na_PY << " I_K_PY " << I_K_PY << " I_ext " << I_Ext;
    // ODE
    dxdt[v_index] = I_Na_PY + I_K_PY + I_Leak_PY + I_Ext - I_Syn;
  }
};

} // insilico

#endif
