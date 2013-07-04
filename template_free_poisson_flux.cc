
#include "template_free_poisson_flux.h"

namespace oomph
{


  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  TFPoissonFluxEquations::
  TFPoissonFluxEquations(FiniteElement* const &bulk_el_pt,
                         const int &face_index) :
    FaceElement()
  {
# ifdef PARANOID
    // Check that the element is not a refineable 3d element
    if((bulk_el_pt->dim()==3)
       && (dynamic_cast<RefineableElement*>(bulk_el_pt) != 0))
      {
        //Issue a warning
        OomphLibWarning("This flux element will not work correctly if nodes are hanging\n",
                        OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index,this);

    // Initialise the prescribed-flux function pointer to zero
    Flux_fct_pt = 0;

    // Extract the dimension of the problem from the dimension of
    // the first node
    Dim = this->node_pt(0)->ndim();

    // Set up U_index_poisson. Initialise to zero, which probably won't change
    // in most cases, oh well, the price we pay for generality.
    U_index_poisson = dynamic_cast<TFPoissonEquations*>(bulk_el_pt)
      ->u_index_poisson();
  }


  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  void TFPoissonFluxEquations::
  fill_in_generic_residual_contribution_poisson_flux(
                                                     Vector<double> &residuals, DenseMatrix<double> &jacobian,
                                                     const unsigned& flag)
  {
    //Find out how many nodes there are
    const unsigned n_node = nnode();

    //Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);

    //Set the value of Nintpt
    const unsigned n_intpt = integral_pt()->nweight();

    //Set the Vector to hold local coordinates
    Vector<double> s(Dim-1);

    //Integers to hold the local equation and unknown numbers
    int local_eqn=0;

    // Locally cache the index at which the variable is stored
    const unsigned u_index_poisson = U_index_poisson;

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {

        //Assign values of s
        for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

        //Get the integral weight
        double w = integral_pt()->weight(ipt);

        //Find the shape and test functions and return the Jacobian
        //of the mapping
        double J = shape_and_test(s,psif,testf);

        //Premultiply the weights and the Jacobian
        double W = w*J;

        //Need to find position to feed into flux function, initialise to zero
        Vector<double> interpolated_x(Dim,0.0);

        //Calculate velocities and derivatives
        for(unsigned l=0;l<n_node;l++)
          {
            //Loop over velocity components
            for(unsigned i=0;i<Dim;i++)
              {
                interpolated_x[i] += nodal_position(l,i)*psif[l];
              }
          }

        //Get the imposed flux
        double flux = 0, elemental_flux = 0;
        get_flux(interpolated_x, flux);
        get_elemental_flux(s, elemental_flux);
        flux += elemental_flux;

        //Now add to the appropriate equations

        //Loop over the test functions
        for(unsigned l=0;l<n_node;l++)
          {
            local_eqn = nodal_local_eqn(l,u_index_poisson);
            /*IF it's not a boundary condition*/
            if(local_eqn >= 0)
              {
                //Add the prescribed flux terms
                residuals[local_eqn] -= flux*testf[l]*W;

                // Imposed traction doesn't depend upon the solution,
                // --> the Jacobian is always zero, so no Jacobian
                // terms are required
              }
          }
      }
  }
}
