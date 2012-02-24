
# include <fstream>
# include <vector>
# include <functional>

# include "./quadrule.h"

/*
  Compile using:
  g++ generate_quadrature_rules_driver.cc quadrule.cc -Wall -g -Wconversion --std=c++0x

  To get other rules search replace legendre_dr_compute with
  clenshaw_curtis_compute or fejer2_compute (or other rules) as appropriate.

  Details of the function used (from quadrule.cpp):

  Purpose:

  LEGENDRE_DR_COMPUTE: Gauss-Legendre quadrature by Davis-Rabinowitz method.

  Discussion:

  The integral:

  Integral ( -1 <= X <= 1 ) F(X) dX

  The quadrature rule:

  Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )

  ............

  Parameters:

  Input, int ORDER, the order.
  ORDER must be greater than 0.

  Output, double XTAB[ORDER], the abscissas.

  Output, double WEIGHT[ORDER], the weights.
*/

void output_weight(const std::vector<double> &weight, const unsigned &length,
		   std::ofstream &output_stream)
{
  output_stream << "{";
  for(unsigned i=0; i<length; i++)
    output_stream << weight[i] << ",";
  output_stream << "}," << std::endl;
}

void output_knot(const std::vector<std::vector<double> > &knot, const unsigned &length,
		 const unsigned &dim, std::ofstream &output_stream)
{
  output_stream << "{";
  for(unsigned i=0; i<length; i++)
    {
      output_stream << "{";
      for(unsigned j=0; j<dim; j++)
	output_stream << knot[i][j] << ",";
      output_stream << "}, ";
    }
  output_stream << "}," << std::endl;
}

void quad_compute_weights(const std::function<void(unsigned,double*,double*)>
			 &one_d_weights_compute,
			  std::ofstream &rules_stream, const unsigned &max_order)
{

  //============================================================
  /// Pre-calculate weights in 1D
  //============================================================

  std::vector<std::vector<double> > weight_1d(max_order+1);
  for(unsigned order=1; order<=max_order; order++)
    {
      double temp_weight[order], temp_knot[order];
      one_d_weights_compute(order,temp_knot,temp_weight);

      weight_1d[order].resize(order);
      for(unsigned i=0; i<order; i++)
	weight_1d[order][i] = temp_weight[i];
    }


  //=================================================================
  // Calculate + output weights in 1 to 3D
  //=================================================================

  // Write some array structure stuff and dim 0
  rules_stream << "{" << std::endl
	       << "// Dim = 0" << std::endl
	       << "{{}}," << std::endl
	       << std::endl;


  // Calculate and output the weights for dimension 1
  rules_stream << "{" << std::endl
	       << "//Dim = 1" << std::endl
	       << "{}, // order 0" << std::endl;
  for(unsigned order=1; order<=max_order; order++)
    {
      // Dump weights
      output_weight(weight_1d[order],order,rules_stream);
    }
  rules_stream << "}" << std::endl
	       << std::endl;


  // Calculate and output weights for dimension 2
  rules_stream << "{" << std::endl
	       << "//Dim = 2" << std::endl
	       << "{}, // order 0" << std::endl;
  for(unsigned order=1; order<=max_order; order++)
    {
      // Calculate the tensor product to get dimension 2
      std::vector<double> weight_2d(order*order,0.0);
      for(unsigned i=0; i<order; i++)
	for(unsigned j=0; j<order; j++)
	  weight_2d[order*i + j] = weight_1d[order][i]*weight_1d[order][j];

      // Dump weights
      output_weight(weight_2d,order*order,rules_stream);
    }
  rules_stream << "}" << std::endl
	       << std::endl;


  // Calculate and output weights for dimension 3
  rules_stream << "{" << std::endl
	       << "//Dim = 3" << std::endl
	       << "{}, // order 0" << std::endl;
  for(unsigned order=1; order<=max_order; order++)
    {
      // Calculate the tensor product to get dimension 3
      std::vector<double> weight_3d(order*order*order,0.0);
      for(unsigned i=0; i<order; i++)
	for(unsigned j=0; j<order; j++)
	  for(unsigned k=0; k<order; k++)
	    {
	      weight_3d[order*order*i + order*j + k]
		= weight_1d[order][i]*weight_1d[order][j]*weight_1d[order][k];
	    }

      // Dump weights
      output_weight(weight_3d,order*order*order,rules_stream);
    }
  rules_stream << "}" << std::endl
	       << std::endl;

  rules_stream << "}" << std::endl
	       << std::endl;
}

void quad_compute_knots(const std::function<void(unsigned,double*,double*)>
			&one_d_knots_compute,
			std::ofstream &rules_stream, const unsigned &max_order)
{

  //============================================================
  /// Pre-calculate knots in 1D
  //============================================================

  std::vector<std::vector<double> > knot_1d(max_order+1);
  for(unsigned order=1; order<=max_order; order++)
    {
      double temp_weight[order], temp_knot[order];
      one_d_knots_compute(order,temp_knot,temp_weight);

      knot_1d[order].resize(order);
      for(unsigned i=0; i<order; i++)
	knot_1d[order][i] = temp_knot[i];
    }


  //=================================================================
  // Calculate + output knots in 1 to 3D
  //=================================================================

  // Write some array structure stuff and dim 0
  rules_stream << "{" << std::endl
	       << "// Dim = 0" << std::endl
	       << "{{{}}}," << std::endl
	       << std::endl;


  // Calculate and output the knots (dimension 1)
  rules_stream << "{" <<std::endl
  	       << "// Dim = 1" << std::endl
  	       << "{{}} // order 0" << std::endl;
  for(unsigned order=1; order<=max_order; order++)
    {
      // Convert to a vector of vectors (to use output function)
      std::vector<std::vector<double> > knot_1d_out;
      for(unsigned i=0; i<order; i++)
  	knot_1d_out.push_back(std::vector<double>(1,knot_1d[order][i]));

      // Dump knots
      output_knot(knot_1d_out, order, 1, rules_stream);
    }
  rules_stream << "}" << std::endl
	       << std::endl;


  // Calculate and output the knots (dimension 2)
  rules_stream << "{" <<std::endl
	       << "// Dim = 2" << std::endl
	       << "{{}} // order 0" << std::endl;
  for(unsigned order=1; order<=max_order; order++)
    {
      // Get tensor product for 2d
      std::vector<std::vector<double> > knot_2d(order*order);
      for(unsigned i=0; i<order; i++)
	for(unsigned j=0; j<order; j++)
	  {
	    knot_2d[i*order +j] = std::vector<double>(2,0.0);
	    knot_2d[i*order +j][0] = knot_1d[order][i];
	    knot_2d[i*order +j][1] = knot_1d[order][j];
	  }

      // Output
      output_knot(knot_2d, order*order, 2, rules_stream);
    }
  rules_stream << "}" << std::endl
	       << std::endl;


  // Calculate and output the knots (dimension 3)
  rules_stream << "{" <<std::endl
	       << "// Dim = 3" << std::endl
	       << "{{}} // order 0" << std::endl;
  for(unsigned order=1; order<=max_order; order++)
    {
      // Get tensor product for 3d
      std::vector<std::vector<double> > knot_3d(order*order*order);
      for(unsigned i=0; i<order; i++)
	for(unsigned j=0; j<order; j++)
	  for(unsigned k=0; k<order; k++)
	    {
	      knot_3d[i*order*order + j*order + k] = std::vector<double>(3,0.0);
	      knot_3d[i*order*order + j*order + k][0] = knot_1d[order][i];
	      knot_3d[i*order*order + j*order + k][1] = knot_1d[order][j];
	      knot_3d[i*order*order + j*order + k][1] = knot_1d[order][k];
	    }

      // Output
      output_knot(knot_3d, order*order*order, 3, rules_stream);
    }
  rules_stream << "}" << std::endl
	       << std::endl;

  rules_stream << "}" << std::endl
	       << std::endl;


  // Done so close output file
  rules_stream.close();
}

int main()
{
  const unsigned max_order = 5;

  // Set up filename
  char filename[100];
  sprintf(filename,"VariableOrderQuadrature.cc");
  std::ofstream rules_stream;
  rules_stream.precision(16);
  rules_stream.open(filename);

  rules_stream << "// Automatically generated by generate_quadrature_rules_driver,"
	       << " based on QUADRULE." << std::endl
	       << "//See https://github.com/davidshepherd7/oomph-lib-additions/tree/master/generate_quadrature_rules ."
	       << std::endl;


  // List of functions to compute weights and knots With
  std::vector<std::string> class_name;
  std::vector<std::function<void(unsigned,double*,double*)> > function;

  class_name.push_back("VariableOrderGaussLegendre");
  function.push_back(legendre_dr_compute);

  class_name.push_back("VariableOrderClenshawCurtis");
  function.push_back(clenshaw_curtis_compute);

  class_name.push_back("VariableOrderFejerSecond");
  function.push_back(fejer2_compute);


  // Create and output to file the quadrature rules for all QElement schemes:
  for(unsigned i=0; i<function.size(); i++)
    {
      // Weights:
      rules_stream << "const " << class_name[i] << "::weights_data_structure" << std::endl
		   << class_name[i] << "::Weights =" << std::endl;
      quad_compute_weights(function[i], rules_stream, max_order);

      // Knots
      rules_stream << "const " << class_name[i] << "::knots_data_structure" << std::endl
		   << class_name[i] << "::Knots =" << std::endl;
      quad_compute_knots(function[i], rules_stream, max_order);
    }

}
