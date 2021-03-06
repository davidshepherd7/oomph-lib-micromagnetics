
# include <fstream>
# include <vector>
# include <functional>

# include "./quadrule.h"

/*
  Compile using:
  g++ generate_quadrature_rules_driver.cc quadrule.cc -Wall -g -Wconversion --std=c++0x

  Sorry this code is such a mess, it was written to work for all dimensions but
  the resulting file was far too large so it is rather overcomplicating things
  for just 1D calculations.

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


namespace oomph_compute_weights
{

  // Quadrature rule function pointer typedef
  typedef void (*quadrature_compute_fn)(int,double*,double*);

  void quad_compute_weights_knots(const quadrature_compute_fn &one_d_weights_compute,
                                  const std::vector<unsigned> order_list,
                                  std::ofstream &rules_stream)
  {

    //============================================================
    /// Calculations
    //============================================================

    std::vector<std::vector<double> > weight_1d(order_list.size());
    std::vector<std::vector<double> > knot_1d(order_list.size());
    for(unsigned i=0; i<order_list.size(); i++)
      {
        unsigned order = order_list[i];

        // Compute weights and knots for this order
        double temp_weight[order], temp_knot[order];
        one_d_weights_compute(order,temp_knot,temp_weight);

        // Move weights to a vector
        weight_1d[i].resize(order);
        for(unsigned j=0; j<order; j++)
          weight_1d[i][j] = temp_weight[j];

        // Move knots to a vector
        knot_1d[i].resize(order);
        for(unsigned j=0; j<order; j++)
          knot_1d[i][j] = temp_knot[j];
      }


    //=================================================================
    // Output weights
    //=================================================================
    rules_stream << "double weights_data_array[] =" << std::endl << "{" << std::endl;
    for(unsigned i=0; i<order_list.size(); i++)
      {
        unsigned order = order_list[i];

        // Dump weights
        for(unsigned j=0; j<order; j++)
          rules_stream << weight_1d[i][j] << ",";
        rules_stream << std::endl;
      }
    rules_stream << "};" << std::endl << std::endl;


    //=================================================================
    // Output knots
    //=================================================================
    rules_stream << "double knots_data_array[] =" << std::endl << "{" << std::endl;
    for(unsigned i=0; i<order_list.size(); i++)
      {
        unsigned order = order_list[i];

        // Dump knots
        for(unsigned j=0; j<order; j++)
          rules_stream << knot_1d[i][j] << ",";
        rules_stream << std::endl;
      }
    rules_stream << "};" << std::endl << std::endl;


    //============================================================
    /// Output list of orders used and number of orders used
    //============================================================
    rules_stream << "unsigned order_array[] = {";
    for(unsigned i=0; i<order_list.size(); i++)
      rules_stream << order_list[i] << ",";
    rules_stream << "};" << std::endl;

    rules_stream << "unsigned order_array_length = " <<  order_list.size() << ";" << std::endl;

  }

}

using namespace oomph_compute_weights;

int main()
{

  // Create the list of order_list to compute rules for (1-50 + some powers of 2)
  std::vector<unsigned> order_list;
  for(unsigned i=1; i<51; i++)
    order_list.push_back(i);
  order_list.push_back(64); order_list.push_back(128); order_list.push_back(256);

  // // Code for testing (few orders)
  // std::vector<unsigned> order_list;
  // order_list.push_back(1);
  // order_list.push_back(2);
  // order_list.push_back(3);

  // for(unsigned i=0; i<order_list.size(); i++)
  //   std::cout << order_list[i] << ", ";



  // Set up filename
  char filename[100];
  sprintf(filename,"VariableOrderQuadrature.cc");
  std::ofstream rules_stream;
  rules_stream.precision(16);
  rules_stream.open(filename);

  // Create list of quadrature types to compute rules for
  std::vector<std::string> class_name;
  std::vector<quadrature_compute_fn> function;

  class_name.push_back("VariableOrderGaussLegendre");
  function.push_back(legendre_dr_compute);

  class_name.push_back("VariableOrderClenshawCurtis");
  function.push_back(clenshaw_curtis_compute);

  class_name.push_back("VariableOrderFejerSecond");
  function.push_back(fejer2_compute);


  // Start of file info:
  rules_stream << "// Automatically generated by generate_quadrature_rules_driver,"
               << " based on QUADRULE." << std::endl
               << "//See https://github.com/davidshepherd7/oomph-lib-additions/tree/master/generate_quadrature_rules ."
               << std::endl;

  rules_stream << "weights_data_structure::weights_data_structure(const unsigned scheme, const bool weights)"
               << std::endl
               << "{\nswitch(scheme)\n{" << std::endl;

  // Create and output to file the quadrature rules for all QElement schemes:
  for(unsigned i=0; i<function.size(); i++)
    {
      rules_stream << "case " << i << " : \n{\n" << std::endl;

      // Compute and output as arrays
      quad_compute_weights_knots(function[i], order_list, rules_stream);

      // Write the arrays into the object
      rules_stream << "\n// Construct the weights or knots as requested\nif(weights)\nconstruction_helper(weights_data_array, order_array, order_array_length);\nelse\nconstruction_helper(knots_data_array, order_array, order_array_length);\n"
                   << std::endl;

      rules_stream << "break;\n}\n\n" << std::endl;
    }

  rules_stream << "default : " <<std::endl
               << "throw OomphLibError(\"The given quadrature scheme does not exist.\",
OOMPH_CURRENT_FUNCTION,\n"
               << "\"OOMPHEXCEPTIONLOCATION\");" << std::endl;

  rules_stream << "}\n}\n\n" << std::endl;

  // Create the function calls to actually create the static objects
  rules_stream << "// Call the constructors." << std::endl;
  for(unsigned i=0; i<function.size(); i++)
    {
      rules_stream << "const weights_data_structure " << class_name[i] << "::"
                   << "Weights(" << i << ",1);" <<std::endl;

      rules_stream << "const weights_data_structure " << class_name[i] << "::"
                   << "Knots(" << i << ",0);" <<std::endl;
    }

  rules_stream.close();
}
