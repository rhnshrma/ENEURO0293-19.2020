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

#ifndef INCLUDED_N_HTC_HPP
#define INCLUDED_N_HTC_HPP

#include "insilico/core/engine.hpp"

#include "insilico/current/I_Na_TC.hpp"
#include "insilico/current/I_K_TC.hpp"
#include "insilico/current/I_Leak_HTC.hpp"
#include "insilico/current/I_T_TC.hpp"
#include "insilico/current/I_THT_HTC.hpp"
#include "insilico/current/I_AHP_HTC.hpp"
#include "insilico/current/I_H_HTC.hpp"

namespace insilico {

class N_HTC_mGluR: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {

/*    std::ostringstream strm; 
    strm << index;*/
    unsigned v_index = engine::neuron_index(index, "v");
/*    std::cerr << "N1 "  << "\n" ;     
*/    double v = variables[v_index];
/*     std::cerr << strm;
*//*     std::cerr << "\n before \n";
*/    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
/*     std::cerr << "after \n";
*//*    strm.flush();
*//*    std::cerr << "N2 "  << g1_indices.size() <<"\n" ; 
*/    std::vector<double> gsyn = engine::get_pre_neuron_values(index, "gsyn");
    double I_Syn = 0;
/*    strm.flush();
*//*    strm <<   gsyn.size();      
*//*    std::cerr << "N3 "  << "\n" ; 
*/    std::vector<double> esyn = engine::get_pre_neuron_values(index, "esyn");
/*    strm.flush();
*//*    strm <<  esyn.size();   
*//*    std::cerr << "N4 "  << "\n" ; 
*/    std::vector<double> v_pre = engine::get_pre_neuron_values(index, "v_pre_index");
/*    strm.flush();
*//*    strm << v_pre.size();
*//*    std::cerr << "N5"  << "\n" ; 
*/    double I_Syn_gj = 0;



    /*for(std::vector<unsigned>::size_type iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn = I_Syn + variables[g1_indices.at(iterator)]*gsyn.at(iterator)*(v - esyn.at(iterator)); 
           
    }*/

       std::cerr << "N11"  << "\n"   ; 

    for(std::vector<unsigned>::size_type iterator = 0; iterator < v_pre.size(); ++iterator) {

      I_Syn_gj = I_Syn_gj +   engine::neuron_value(index,"g_gj")*(v - variables[v_pre.at(iterator)]);
 
    }

     engine::neuron_value(index,"v_p",g1_indices.size());
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

 std::cerr << "N7"  << "\n";

     double I_Na = engine::neuron_value(index,"I_Na_TC");
     double I_K = engine::neuron_value(index,"I_K_TC");
     double I_HTC_Leak = engine::neuron_value(index,"I_Leak_HTC");
     double I_TLT = engine::neuron_value(index,"I_T_TC");
     double I_H = engine::neuron_value(index,"I_H_HTC");
     double I_THT =  engine::neuron_value(index,"I_THT_HTC");
     double I_AHP = engine::neuron_value(index,"I_AHP_HTC");
     double I_Ext = engine::neuron_value(index, "iext");
 std::cerr << "N8"  << "\n";
     engine::neuron_value(index,"I_Syn",I_Syn);
     engine::neuron_value(index,"I_Syn_gj",I_Syn_gj);
   std::cerr << "N9"  << "\n"  ;
    // ODE
    dxdt[v_index] = I_TLT + I_Na + I_K + I_HTC_Leak + I_H + I_THT + I_AHP + I_Ext - I_Syn;
  }
};

} // insilico

#endif
