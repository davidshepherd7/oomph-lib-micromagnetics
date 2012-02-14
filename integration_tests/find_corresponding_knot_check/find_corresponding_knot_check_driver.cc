#include "../../variable_quadrature.h"

using namespace oomph;

unsigned next(const unsigned &order)
{
  if(order ==2) return 4;

  return 2*order -1;
}

int main()
{

  unsigned max_order = 49;

  VariableClenshawCurtis quad_scheme;
  quad_scheme.set_dim(1);

  for(unsigned order=2; order <=max_order; order=next(order))
    {
      quad_scheme.set_order(order);
      std::cout << "Low order = " << order
		<< ". High order = " << max_order << std::endl;

      for(unsigned kn=0; kn<quad_scheme.nweight(); kn++)
	{
	  std::cout << kn << " = "
		    << quad_scheme.find_corresponding_knot(kn,max_order)
		    << std::endl;
	}
    }

}
