
#include "micromag_factories.h"

#include "llg_problem.h"
#include "boundary_element_handler.h"

namespace oomph
{
  namespace factories
  {

    void bem_handler_factory(const BemBoundaryData& bem_boundaries,
                             const unsigned& phi_index,
                             const unsigned& phi_1_index,
                             const CornerDataInput& input_corner_data,
                             bool use_hlib,
                             BoundaryElementHandler& new_bem_handler)
    {

      if(use_hlib)
        {
          std::string err = "not yet implemented";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      // Get the first finite element we can find in the a bulk mesh on
      // which we are going to construct our bem mesh. Assuimg that all
      // elements in all the meshes are the same type...
      FiniteElement* sample_fele_pt =
        bem_boundaries[0].second->finite_element_pt(0);

      // Figure out which element type we should use in the bem mesh
      // (based on the element type used in the bulk mesh) and store the
      // function needed to create them.
      new_bem_handler.Bem_element_factory_fpt = LLGFactories::
        bem_element_factory_factory(sample_fele_pt);

      // Create an integration scheme ??ds move this outside somewhere...
      new_bem_handler.integration_scheme_pt() = LLGFactories::
        variable_order_integrator_factory(sample_fele_pt);

      // Set indexes to look in for phi/phi1 variables
      new_bem_handler.set_input_index(phi_1_index);
      new_bem_handler.set_output_index(phi_index);

      // Copy in the list of boundaries to operate on
      new_bem_handler.Bem_boundaries = bem_boundaries;

      // Now build it
      new_bem_handler.build(input_corner_data);
    }


    BoundaryElementHandler* bem_handler_factory
    (const BemBoundaryData& bem_boundaries,
     const unsigned& phi_index,
     const unsigned& phi_1_index,
     const CornerDataInput& input_corner_data,
     bool use_hlib)
    {
      // Create with new, fill in with factory
      BoundaryElementHandler* bem_handler_pt = new BoundaryElementHandler;
      bem_handler_factory(bem_boundaries, phi_index, phi_1_index,
                          input_corner_data, use_hlib, *bem_handler_pt);
      return bem_handler_pt;
    }


  }
}
