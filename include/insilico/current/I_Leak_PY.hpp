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

#ifndef INCLUDED_I_Leak_PY_HPP
#define INCLUDED_I_Leak_PY_HPP

#include "insilico/core/engine.hpp"

namespace insilico {

class I_Leak_PY {
 public:
  static void current(const state_type &variables, state_type &dxdt, const double t, unsigned index) {    

    unsigned v_index = engine::neuron_index(index, "v");

    
    double v = variables[v_index];



    										// Current
     engine::neuron_value(index,"I_Leak_PY",(-1*(0.1*(v+61.0))));

  } // function current
}; // class I_Leak_PY

} // insilico

#endif

