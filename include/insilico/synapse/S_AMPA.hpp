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

#ifndef INCLUDED_S_AMPA_HPP
#define INCLUDED_S_AMPA_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"
#include <random>


namespace insilico {



class random4 {
 public:
  template<class T>
  static T rand() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(0,1);
    return dist(gen);
  }
};

class S_AMPA: public Synapse {
 public:
  void ode_set(state_type& variables, state_type& dxdt, const double t, unsigned index) {
    
    unsigned g1_index = engine::synapse_index(index, "g1");
    unsigned pre_neuron_index = engine::synapse_value(index, "pre");





    double g1             =    variables[g1_index];
    double last_spike     =    engine::neuron_value(pre_neuron_index,"last_spike");
    double delay          =    engine::synapse_value(index,"delay");
    double spike_duration =    engine::synapse_value(index,"spike_duration");
    double T              =    0;
    double prob           =    engine::synapse_value(index,"prob");
    double to_spike       =    engine::synapse_value(index,"to_spike");
    double spiking        =    engine::synapse_value(index,"spiking");

    // synapse logic for delay for recently spiked neuron


    if(((t - last_spike - delay) < (spike_duration))  &&  ((t-last_spike - delay) > 0))
    {

        if (spiking<0)
        {
            double randVar = random4::rand<double>();
            if (randVar < prob)
            {
                engine::synapse_value(index,"to_spike",1);
                T = 0.5;
            }
            else
            {
                engine::synapse_value(index,"to_spike",-1);
                T = 0; 
            }
            engine::synapse_value(index,"spiking",1);
        }
        else
        {
            if (to_spike>0)
            {
                T = 0.5;
            }
            else
            {
                T = 0;
            }

        }
    }
    else
    {
        engine::synapse_value(index,"spiking",-1);
        T = 0;
    }


 
    // constants from file


    // ODE set
    dxdt[g1_index] = (0.98*T*(1-g1) - 0.180*g1);


  } // function ode_set
}; // class S_AMPA

} // insilico

#endif


