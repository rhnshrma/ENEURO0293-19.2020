/*
 current/I_K_PY.hpp - Current flowing across neuronal membrane due to
                          Potassium (K) conductance. (McCarthy et al., 2008)

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

 Current that flows through Potassium (K) channel due to the potential difference
 caused by Potassium (K) conductance across neuronal membrane. (Hodgkin-Huxley, 1952)
*/

#ifndef INCLUDED_I_K_TC_HPP
#define INCLUDED_I_K_TC_HPP

#include "insilico/core/engine.hpp"
#include "vtrap.cpp"

namespace insilico {

class I_K_TC {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, const unsigned index) {
    double gk = 10, ek = -100;
                                                            // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned n_index = engine::neuron_index(index, "n_K_TC");

    double v = variables[v_index];
    double n = variables[n_index];

    double vt = v + 25;


                                                            // Calculating time constants
    double alpha_n = 0.032*vtrap((15.0-vt),5.0);
    double beta_n = 0.5*(exp((10-vt)/40.0));

    double n_inf = alpha_n/(alpha_n+beta_n);
    double tau_n = 1/(alpha_n + beta_n);


    if (n<0) variables[n_index] = 0;
    if (n>1) variables[n_index] = 1;
                                                            // ODE set
    dxdt[n_index] = (n_inf-n)/tau_n;


                                                            // Current
    engine::neuron_value(index, "I_K_TC", (-1*(gk*pow(n,4)*(v - ek))));

  } // function current
}; // class I_potassium_thalamo-cortical cells

} // insilico

#endif
