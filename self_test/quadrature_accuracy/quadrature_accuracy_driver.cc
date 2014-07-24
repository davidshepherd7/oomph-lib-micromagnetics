
#include "generic.h"
#include "micromag.h"

using namespace oomph;
using namespace Factories;
using namespace MManipulation;

inline double x_func(const double& x, const double& y)
{
  return x + 0.3333*y;
}

class LinearFunction : public ElementalFunction
{
  double call(const GeneralisedElement* ele_pt, CachingMMInterpolator* intp_pt) const override
    {
      double x = x_func(intp_pt->x()[0], intp_pt->x()[1]);
      return 0.7*x;
    }
};

class QuadraticFunction : public ElementalFunction
{
  double call(const GeneralisedElement* ele_pt, CachingMMInterpolator* intp_pt) const override
  {
    double x = x_func(intp_pt->x()[0], intp_pt->x()[1]);
    return std::sqrt(2)*x*x - x;
  }
};

class CubicFunction : public ElementalFunction
{
  double call(const GeneralisedElement* ele_pt, CachingMMInterpolator* intp_pt) const override
  {
    double x = x_func(intp_pt->x()[0], intp_pt->x()[1]);
    return std::sqrt(3)*x*x*x + std::sqrt(2)*x*x - x;
  }
};

class QuarticFunction : public ElementalFunction
{
  double call(const GeneralisedElement* ele_pt, CachingMMInterpolator* intp_pt) const override
  {
    double x = x_func(intp_pt->x()[0], intp_pt->x()[1]);
    return std::sqrt(5)*x*x*x*x + std::sqrt(3)*x*x*x + std::sqrt(2)*x*x - x;
  }
};

inline void set_quadratures(Mesh* mesh_pt,
                            NodalQuadratureFactoryFctPt quad_factory_pt,
                            Vector<Integral*>& quad_pts)
{
  const double mean_elemental_volume = 0; // don't use the hacky rescaling
                                          // for now, don't need this.

  const unsigned n_ele = mesh_pt->nelement();
  for(unsigned ele=0; ele<n_ele; ele++)
    {
      FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);

      // Create an integration scheme as specified by the factory function.
      // If the scheme is null then use gaussian
      Integral* quad_pt = quad_factory_pt(ele_pt, mean_elemental_volume);
      if(quad_pt == 0)
        {
          quad_pt = gauss_integration_factory(ele_pt->dim(),
                                              ele_pt->nnode_1d(),
                                              ele_pt->element_geometry());
        }
      ele_pt->set_integration_scheme(quad_pt);
      quad_pts.push_back(quad_pt);
    }
}

inline Integral* high_order_gauss_factory(const FiniteElement* ele_pt,
                                          const double& mean_size)
{
  return Factories::gauss_integration_factory(ele_pt->dim(), 4,
                                              ele_pt->element_geometry());
}

int main()
{
  // build
  // ============================================================


  const double scaling = 1.0;
  const unsigned refinement = 5;

  TimeStepper* ts_pt = time_stepper_factory("steady");

  // construct meshes
  Vector<std::pair<Mesh*, std::string> > mesh_pts;
  mesh_pts.push_back(std::make_pair(llg_mesh_factory("st_square", refinement, ts_pt, scaling, 2),
                                    "st_square"));
  mesh_pts.push_back(std::make_pair(llg_mesh_factory("sq_square", refinement, ts_pt, scaling, 2),
                                    "sq_square"));
  mesh_pts.push_back(std::make_pair(llg_mesh_factory("ut_square", refinement, ts_pt, scaling, 2),
                                    "ut_square"));
  mesh_pts.push_back(std::make_pair(llg_mesh_factory("sq_annular", refinement, ts_pt, scaling, 2),
                                    "sq_annular"));

  Vector<std::pair<NodalQuadratureFactoryFctPt, std::string> > quad_fact_pts;
  quad_fact_pts.push_back(std::make_pair(&high_order_gauss_factory,
                                         "exact "));
  quad_fact_pts.push_back(std::make_pair(nodal_quadrature_factory_factory("gauss"),
                                         "gauss "));
  quad_fact_pts.push_back(std::make_pair(nodal_quadrature_factory_factory("lnodal"),
                                         "lnodal"));
  // quad_fact_pts.push_back(std::make_pair(nodal_quadrature_factory_factory("nodal"),
  //                                        "nodal"));

  // Vector to store pointers to quadratures so we can delete them later
  Vector<Integral*> quad_pts;

  Vector<std::pair<ElementalFunction*, std::string> > functions;
  functions.push_back(std::make_pair(new LinearFunction, "LinearFunction"));
  functions.push_back(std::make_pair(new QuadraticFunction, "QuadraticFunction"));
  functions.push_back(std::make_pair(new CubicFunction, "CubicFunction"));
  functions.push_back(std::make_pair(new QuarticFunction, "QuarticFunction"));


  // run
  // ============================================================

  // ??ds the testing and output code here is awful but c++ makes it so
  // hard to do things in a nicely functional way. So we check the and
  // output the results in the same code as the calculations, all mixed
  // up...

  //??ds use map?
  Vector<double> results;

  // high precision
  std::cout << std::setprecision(12);

  // loop parameters
  const unsigned n_msh = mesh_pts.size();
  const unsigned n_func = functions.size();

  // Do the calculations!
  for(unsigned i_func=0; i_func<n_func; i_func++)
    {
      std::cout << functions[i_func].second << ":" << std::endl;
      for(unsigned msh=0; msh<n_msh; msh++)
        {
          results.clear();

          for(unsigned j=0; j<quad_fact_pts.size(); j++)
            {

              // Set the quadratures
              set_quadratures(mesh_pts[msh].first, quad_fact_pts[j].first,
                              quad_pts);

              // Integrate
              double result = integrate_over_mesh(functions[i_func].first,
                                                  mesh_pts[msh].first,
                                                  0);
              std::cout << mesh_pts[msh].second << " "
                        << quad_fact_pts[j].second << " "
                        << result << std::endl;

              results.push_back(result);
            }

          // Check nodal value is close to high order gauss
          double relerr = std::abs((results[0] - results[2])/results[0]);
          if(relerr > 1e-3)
            {
              std::cerr << "nodal quadrature relerr of " << relerr
                        << "larger than tol when integrating "
                        <<  functions[i_func].second
                        << " on mesh " << mesh_pts[msh].second << std::endl;

              // fail
              return i_func*n_func + msh + 1;
            }
          else
            {
              std::cout << "relerr of " << relerr << " ok." << std::endl << std::endl;
            }

        }
      std::cout << std::endl;
    }


  // Clean up
  for(unsigned msh=0; msh<n_msh; msh++)
    {
      delete mesh_pts[msh].first; mesh_pts[msh].first = 0;
    }
  for(unsigned i_func=0; i_func<n_func; i_func++)
    {
      delete functions[i_func].first; functions[i_func].first = 0;
    }
  for(unsigned j=0; j<quad_pts.size(); j++)
    {
      delete quad_pts[j]; quad_pts[j] = 0;
    }

  return 0;
}
