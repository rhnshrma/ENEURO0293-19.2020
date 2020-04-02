// Reference : Pankaj Sah(1999)

#ifndef I_SAHP
#define I_SAHP

#include "insilico/core/engine.hpp"

namespace insilico {

class I_sahp {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, int index)
 {
    double G_sahp = 0.266, E_sahp  = -72.0;
    unsigned v_index = insilico::engine::neuron_index(index, "v");
    unsigned Ca_conc_index = insilico::engine::neuron_index(index, "Ca_conc");
    unsigned R_index = insilico::engine::neuron_index(index, "I_sahp_R");
    unsigned CA1R_index = insilico::engine::neuron_index(index, "I_sahp_CA1R");
    unsigned CA2R_index = insilico::engine::neuron_index(index, "I_sahp_CA2R");
    unsigned CA3R_index = insilico::engine::neuron_index(index, "I_sahp_CA3R");
    unsigned CA4R_index = insilico::engine::neuron_index(index, "I_sahp_CA4R");
    unsigned O_index = insilico::engine::neuron_index(index, "I_sahp_O");  
  
    double Ca_conc = variables[Ca_conc_index];
    double v = variables[v_index];
    
    double R = variables[R_index] = std::max<double>(variables[R_index],0.000001);
    double CA1R = variables[CA1R_index] = std::max<double>(variables[CA1R_index],0.000001);
    double CA2R = variables[CA2R_index] = std::max<double>(variables[CA2R_index],0.000001);
    double CA3R = variables[CA3R_index] = std::max<double>(variables[CA3R_index],0.000001);
    double CA4R = variables[CA4R_index] = std::max<double>(variables[CA4R_index],0.000001);
    double O = variables[O_index] = std::max<double>(variables[O_index],0.000001);

    double ru = 0.5E-3;
    double rb = 4.0;
    double ro = 0.600;
    double rc = 0.400;
    

    // ODE set
    dxdt[R_index] = (ru*CA1R) - (4.0*rb*R*Ca_conc); 
    dxdt[CA1R_index] = (4.0*rb*R*Ca_conc) - (ru*CA1R) + (2.0*ru*CA2R) - (3.0*rb*Ca_conc*CA1R);
    dxdt[CA2R_index] = (3.0*rb*CA1R*Ca_conc) - (2.0*ru*CA2R) + (3.0*ru*CA3R) - (2.0*rb*CA2R*Ca_conc);
    dxdt[CA3R_index] = (2.0*rb*CA2R*Ca_conc) - (3.0*ru*CA3R) + (4.0*ru*CA4R) - (rb*CA3R*Ca_conc);
    dxdt[CA4R_index] = (rb*CA3R*Ca_conc) - (4.0*ru*CA4R) + (rc*O) - (ro*CA4R);
    dxdt[O_index] = (ro*CA4R) - (rc*O);
  
    double fraction_of_open_channels = (O/(R + CA1R + CA2R + CA3R + CA4R + O));

   // std::cout<<O<<"\t"<<CA1R<<"\t"<<CA2R<<"\t"<<CA3R<<"\t"<<CA4R<<"\t"<<R<<"\t"<<Ca_conc<<"\n \n"; 

    double I_sahp =   G_sahp * fraction_of_open_channels  * (v - E_sahp);
    
    insilico::engine::neuron_value(index,"fraction_of_open_channels",fraction_of_open_channels);
    insilico::engine::neuron_value(index, "I_sahp", I_sahp);
    insilico::engine::neuron_value(index, "R",R);
    insilico::engine::neuron_value(index, "CA1R",CA1R);
    insilico::engine::neuron_value(index, "CA2R",CA2R);
    insilico::engine::neuron_value(index, "CA3R",CA3R);
    insilico::engine::neuron_value(index, "CA4R",CA4R);
    insilico::engine::neuron_value(index, "O",O);
    
   
  } // function current
}; // class I_SAHP

} // insilico

#endif
