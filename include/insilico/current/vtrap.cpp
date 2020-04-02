
#ifndef INCLUDED_vtrap_HPP
#define INCLUDED_vtrap_HPP

namespace insilico 
{
	
	
  double vtrap(double x,double y)
  {
        double val=0;
        if(std::abs(x/y) < 0.000001)
        { 
            val = y*(1 - (x/(y))/2);
        }
        else 
        { 
            val = (x)/(exp((x)/y) - 1.0); 
        }
        return val;
  } 


}

#endif