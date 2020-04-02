/*
 current/I_M_LTS.hpp - Current flowing through slow potassium channel. (Mcarthy-Brown-Kopell, 2008)

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

 Current that flows through the Muscaranic K-Channel (Mcarthy-Brown-Kopell, 2008)
*/

#ifndef INCLUDED_I_M_LTS_HPP
#define INCLUDED_I_M_LTS_HPP

#include "core/engine.hpp"
#include "script/vtrap.cpp"

namespace insilico {

class I_M_LTS {
 public:
  static double current(const state_type &variables, state_type &dxdt, const double t, unsigned index) {    
    double gm = 2, em = -100, qs=3.209;
                                                          // Declaring Constants
    unsigned v_index = engine::neuron_index(index, "v");
    unsigned m_index = engine::neuron_index(index, "m_m_LTS");
    
    double v = variables[v_index];
    double m = variables[m_index];
    

                                                          // Calculating time constants
    double alpha_m = 0.0001*qs*vtrap(-1*(30.0+v),9.0);
    double beta_m  = 0.0001*qs*vtrap((v+30.0),9.0);

    double m_inf=alpha_m/(alpha_m+beta_m);
    double tau_m=1/(alpha_m + beta_m);



    if (m<0) variables[m_index] = 0;  
                                                          // ODE set
    dxdt[m_index]=(m_inf-m)/tau_m;

    
                                                          // Current
    return  (-1*(gm*m*(v - em)));


  } // function current
}; // class I_M_LTS

} // insilico

#endif

