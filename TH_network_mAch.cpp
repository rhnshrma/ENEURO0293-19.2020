  /*
  main.cpp - insilico's example using neuron and synapse for illustrations

  Copyright (C) 2014 Collins Assisi, Collins Assisi Lab, IISER, Pune
  Copyright (C) 2014 Arun Neru, Collins Assisi Lab, IISER, Pune <areinsdel@gmail.com>
  Copyright (C) 2014-2015 Pranav Kulkarni, Collins Assisi Lab, IISER, Pune <pranavcode@gmail.com>

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


#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include "insilico/core.hpp"

#include "insilico/neuron/N_RE_mAch.hpp"
#include "insilico/neuron/N_TC_mAch.hpp"
#include "insilico/neuron/N_HTC_mAch_1.hpp"
#include "insilico/neuron/N_HTC_mAch_2.hpp"
#include "insilico/synapse/S_GJ.hpp"
#include "insilico/synapse/S_AMPA.hpp"
#include "insilico/synapse/S_GABAB.hpp"
#include "insilico/synapse/S_2GABAA.hpp"
#include <boost/numeric/odeint.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <iostream>
#include <fstream>
#include <map>
#include <random>

using namespace boost::numeric::odeint;
using namespace insilico;
using namespace std;

class random {
 public:
  template<class T>
  static T rand(T mean, T sd) {
    random_device rd;
    mt19937_64 gen(rd());
    normal_distribution<> dist(mean, sd);
    return dist(gen);
  } // function rand
};

class random2 {
 public:
  template<class T>
  static T rand() {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> dist(0,0.1);
    return dist(gen);
  }
};

class stoch_driver {
 public:
  void operator()(state_type &variables, state_type &dxdt, const double time) {
	  //loop through the indices of the state vector;
	  //find the corresponding variable string;
	  //get the variableStr_rev and variableStr_pow;
	  
	  std::string varStr;
	  std::string varStrCpy;
	  double rev_pot, mag, power, randVal; 
	  unsigned nid, v_i;
	  bool nerror = false;
	  for(unsigned index = 0 ; index < engine::neuron_end_list_ids.back() ; ++index ) {
	    varStr = engine::variable_name_from_index(index);

	    varStrCpy = varStr;
 	    varStrCpy.append("_stoch");
 	    nid = engine::neuron_id_from_index(index, nerror);
        
	     if(!nerror) {
		    if(engine::neuron_value(nid, varStrCpy) > 0) {
			    varStrCpy = varStr;
		    	//varStrCpy.append("_pow");
		    	//power = engine::neuron_value(nid, varStrCpy);
		     	//varStrCpy = varStr;
		    	//varStrCpy.append("_rev");
		    	//rev_pot = engine::neuron_value(nid, varStrCpy);
		    	varStrCpy = varStr;
		    	varStrCpy.append("_mag");
		    	mag = engine::neuron_value(nid, varStrCpy);
		    	randVal = random::rand<double>(0,0.1);
		    	//v_i = engine::neuron_index(nid, "v");
			    //dxdt[index] = mag*randVal*pow((variables[v_i]-rev_pot),power);
			    dxdt[index] = mag*randVal;
		    }  
		}
	  }
  }
};

typedef double time_type;
class stochastic_euler
{
public:
	typedef std::vector< double > state_type;
	typedef std::vector< double > deriv_type;	
	typedef double value_type;
	typedef double time_type;
	typedef unsigned short order_type;
	typedef boost::numeric::odeint::stepper_tag stepper_category;
	
	static order_type order( void ) { return 1; }
		
	template < class System >
	void do_step( System system , state_type &variables, time_type t, time_type dt ) const
	{
	  unsigned res = variables.size();
		deriv_type det(res), stoch(res) ;
		system.first( variables , det, t ) ;
		system.second( variables, stoch, t );
    /*#pragma omp parallel for schedule(runtime)*/

		for(unsigned i=0 ; i < variables.size() ; ++i) {
			variables[i] += dt * det[i] + sqrt( dt ) * stoch[i];
		}
	}
};

int main(int argc, char **argv) {
  configuration::initialize(argc, argv);
  
  configuration::observe("v");
  configuration::observe("r_H_HTC");
  configuration::observe("h_tht_HTC");
  configuration::observe("m_ahp_HTC");
  configuration::observe("Ca_conc");
  configuration::observe_skipiters(40);
  double time_simul=500.0;

  if(insilico::engine::time_specified){
  	time_simul = insilico::engine::simulation_time;
  }
  
/*  for (unsigned i = 0; i < cmds.size(); i=i+2){
  	std::string cmd=cmds[i];
  	if(cmd.at(1)=='t')
  	{
  		time_simul=insilico::string_to_double(cmds[i+1]);
  	}
  }*/
  
  //engine::generate_neuron<N_Stellate_HR2005>(10);
  //engine::generate_neuron<N_InterNeuron_Wang96>(10);
  //engine::generate_synapse<S_DefaultSynapse>(40);
	
  engine::generate_neuron<N_HTC_mAch_1>(1);
  engine::generate_neuron<N_HTC_mAch_2>(1);
  engine::generate_neuron<N_TC_mAch>(8);
  engine::generate_neuron<N_RE_mAch>(10);
  engine::generate_synapse<S_AMPA>(100);
  engine::generate_synapse<S_2GABAA>(206);
  engine::generate_synapse<S_GABAB>(100);
  

  //  engine::generate_neuron<N_InterNeuron_Wang96>(1);

  state_type variables = engine::get_variables();
  integrate_const(stochastic_euler(), make_pair(engine::driver(), stoch_driver()) , variables,
                  0.0,time_simul, 0.01, configuration::observer());

  configuration::finalize();
  return 0;
}

