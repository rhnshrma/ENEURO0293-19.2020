/*
 current/I_Na_HH1952.hpp - Current flowing across neuronal membrane due to
                           Sodium (Na) conductance. (Hodgkin-Huxley, 1952)

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

/*
 Brief:

 Current that flows through Sodium (Na) channel due to the potential difference
 caused by Sodium (Na) conductance across neuronal membrane. (Hodgkin-Huxley, 1952)
*/

#ifndef INCLUDED_I_Leak_HTC_HPP
#define INCLUDED_I_Leak_HTC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_Leak_HTC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, unsigned index) {    

  	double gleak=0.01, gkleak=engine::neuron_value(index,"gkleak"), eleak=-70, ekleak=-100;

    unsigned v_index = engine::neuron_index(index, "v");

    
    double v = variables[v_index];



    										// Current
    // Current
    engine::neuron_value(index, "I_Leak_HTC", (-1*(gleak*(v-eleak) + gkleak*(v-ekleak))));

  } // function current
}; // class I_Leak_HTC

} // insilico

#endif

