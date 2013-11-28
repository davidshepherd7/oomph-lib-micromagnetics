#include "./micromagnetics_boundary_element.h"


namespace oomph
{
  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  MicromagBEMElementEquations::
  MicromagBEMElementEquations(FiniteElement* const bulk_el_pt, const int &face_index)
    : FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index,this);
  }

  /// Get boundary element matrix contributions for this element using
  /// an adaptive scheme.
  void MicromagBEMElementEquations::
  fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix)
    const
  {
    // Find out dimension of element
    const unsigned el_dim = dim(), node_dim = nodal_dimension();

    //Find out how many nodes there are
    const unsigned n_element_node = nnode();

    // Cast pointer to the a variable order integration scheme
    BaseVariableOrderQuadrature* v_int_pt =
      checked_dynamic_cast<BaseVariableOrderQuadrature* >(integral_pt());


    // Start of adaptive integration scheme
    //====================================

    // Setup error parameters for adaptive integration
    double reltol = 1e-8, reldiff;

    // Set up storage to keep track of orders used:
    //std::vector<unsigned> order_list = {2,4,8,16,32,64,128,256,};
    unsigned array_list[] = {2,4,8,16,32,64,128,256,};
    // unsigned array_list[] = {128,256};
    std::vector<unsigned> order_list(array_list,array_list+7);
    unsigned n_order = order_list.size();
    Vector<unsigned> order_is_calculated(n_order,0);

    // Set up storage for data at each knot for each order
    Vector< Vector< Vector<double> > > x_kn(n_order), normal(n_order);
    Vector< Vector<double> > jw(n_order);
    Vector< Vector<ShapeWithDeepCopy> > test(n_order), psi(n_order);

    // Loop over all source nodes on the boundary
    unsigned n_boundary_node = boundary_mesh_pt()->nnode();
    for(unsigned i_sn=0; i_sn<n_boundary_node; i_sn++)
      {
        // Get coordinates of source node
        Vector<double> source_node_x(node_dim,0.0);
        boundary_mesh_pt()->node_pt(i_sn)->position(source_node_x);

        // Reset the quadrature order
        unsigned i_order = 0;

        // Temporary storage for the values computed
        Vector<double> temp_bm_prev(n_element_node,0.0), temp_bm(n_element_node,5000);
        do
          {
            // Copy the previous result into another vector (for comparison later)
            // and initialise current result vector.
            temp_bm_prev = temp_bm;
            temp_bm.assign(n_element_node,0.0);

            ////////////////////////////////////////
            // Precalculations (if needed)

            // Check if this order of integration has been used yet, if not calculate things
            if( order_is_calculated[i_order] != 1)
              {
                // Get new order
                unsigned new_order = order_list[i_order];

                // Get the number of knots at this order
                unsigned n_knot = v_int_pt->nweight(new_order);

                // Calculate and push back values at each knot
                for(unsigned kn=0; kn<n_knot; kn++)
                  {
                    // Get the local coordinate
                    Vector<double> s(el_dim,0.0);
                    for(unsigned j=0; j<el_dim; j++)
                      s[j] = v_int_pt->knot(kn,j,new_order);

                    // Get the unit normal
                    Vector<double> local_normal(node_dim,0.0);
                    outer_unit_normal(s,local_normal);
                    normal[i_order].push_back(local_normal);

                    // Get the Jacobian*weight
                    double J = J_eulerian(s);
                    jw[i_order].push_back(J * v_int_pt->weight(kn,new_order));
#ifdef PARANOID
                    // Check that the Jacobian of the transformation isn't
                    // inverted or singular.
                    if(std::abs(J) < 1e-6)
                      {
                        std::ostringstream error_msg;
                        error_msg << "Singular Jacobian!";
                        throw OomphLibError(error_msg.str(),
                                            OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
                      }
                    if(J < 0.0)
                      {
                        std::ostringstream error_msg;
                        error_msg << "Inverted Jacobian!";
                        throw OomphLibError(error_msg.str(),
                                            OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
                      }
#endif


                    // Get shape and test(=shape) functions
                    ShapeWithDeepCopy psi_local(n_element_node);
                    shape(s,psi_local);
                    psi[i_order].push_back(psi_local);
                    test[i_order].push_back(psi_local);

                    // Get the global coordinate
                    Vector<double> x_local(node_dim,0.0);
                    for(unsigned l=0; l<n_element_node; l++)
                      {
                        for(unsigned j=0; j<node_dim; j++)
                          x_local[j] += nodal_position(l,j)*psi[i_order][kn][l];
                      }
                    x_kn[i_order].push_back(x_local);

                  }

                // Record that the new order has been calculated
                order_is_calculated[i_order] = 1;
              }

            ////////////////////////////////////////////////////////////
            // The calculation itself

            unsigned order = order_list[i_order];

            double n_knot = v_int_pt->nweight(order);
            for(unsigned kn=0;kn<n_knot;kn++)
              {
                double dgreendn =
                  green_normal_derivative_opt(source_node_x, x_kn[i_order][kn],
                                              normal[i_order][kn]);

                // Loop over test functions, i.e. local/target nodes, adding contributions
                for(unsigned l=0; l<n_element_node; l++)
                  {
                    temp_bm[l] += dgreendn *test[i_order][kn][l] *jw[i_order][kn];
                  }
              }
            ////////////////////////////////////////////////////////////

            // Compare with previous result
            reldiff = max_rel_error(temp_bm, temp_bm_prev);

            // Go to the next order
            i_order++;

            // If we've hit the order limit without being accurate enough give an error
            if((i_order >= n_order) && (reldiff>reltol))
              {
                std::ostringstream error_msg;
                error_msg << "Quadrature order not high enough, returning with relative error "
                          << reldiff << " (reltol is set to " << reltol << ")." << std::endl;
                throw OomphLibWarning(error_msg.str(),
                                      "MicromagBEMElementEquations::fill_in_be_contribution_adaptive",
                                      OOMPH_EXCEPTION_LOCATION);
              }
          }
        // Repeat while the difference is too large (or at least run twice)
        while( (reldiff>reltol) || (i_order == 1) );

        // Add the values in the temp vector to the real boundary matrix
        for(unsigned l=0; l<n_element_node; l++)
          {
            boundary_matrix(l,i_sn) += temp_bm[l];
          }

      } // End of loop over source nodes

  } // End of function

  // Given all the values for a certain order, dump them out
  void MicromagBEMElementEquations::
  dump_values(const Vector< Vector<double> > &x_kn,
              const Vector<double> &source_node_x,
              const Vector<Vector<double> > &normal) const
  {
    // Get the values
    Vector<double> gnd(x_kn.size(),0.0);
    Vector< Vector<double> > s_kn(x_kn.size());
    for(unsigned kn=0; kn<x_kn.size(); kn++)
      {
        gnd[kn] = green_normal_derivative_opt(source_node_x, x_kn[kn], normal[kn]);
      }

    // // Output
    // std::cout << "For source node at " << source_node_x << std::endl;

    // for(unsigned kn=0; kn<x_kn.size(); kn++)
    //   std::cout << "s = " << s_kn[kn] << ", green normal deriv = " << gnd[kn] << std::endl;

  }


  //============================================================
  /// Get the maximum elementwise difference between two vectors relative to the
  /// appropriate element of v1.
  //============================================================
  double MicromagBEMElementEquations::
  max_rel_error(const Vector<double> &v1, const Vector<double> &v2) const
  {
    // Get the element-wise relative difference between the two vectors
    Vector<double> diff;
    VectorOps::relative_abs_vector_diff(v1, v2, diff);

    // Return the maximum
    return *(std::max_element(diff.begin(),diff.end()));
  }

  //======================================================================
  ///
  //======================================================================
  double MicromagBEMElementEquations::
  green_normal_derivative(const Vector<double>& x,
                          const Vector<double>& y,
                          const Vector<double>& n) const
  {
    // Get dimensions
    const unsigned node_dim = nodal_dimension();

#ifdef PARANOID
    // Check that n is a unit vector
    if(std::abs(two_norm(n) - 1.0) > 1e-10)
      {
        throw OomphLibError("n is not a unit vector",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Calculate the distance r from x to y.
    Vector<double> r; VectorOps::vector_diff(y,x,r);
    double norm_r = VectorOps::two_norm(r);

    //??ds if x is a potentially singular point CHEAT
    //??ds could go horribly wrong, probably fine so long as not actually singular
    //??ds only works in 2d anyway

    // This is approximately true because the r and n unit vectors
    // become perpendicular as x approaches y. Hence n.r = 0 and
    // the whole function is zero.
    if(norm_r < 1e-6) return 0;

    // Calculate the unit vector in direction of r
    Vector<double> unit_r = r; VectorOps::normalise(unit_r);

    // Calculate dot product of n with the unit vector in direction r.
    double ndotr = VectorOps::dot(n,unit_r);

    // See write up for details of calculation.
    int mdm1 = - (int(node_dim) - 1);
    return (-1.0/Pi) * ndotr * std::pow((2.0*norm_r), mdm1); //??ds had *2 here for some reason...
  }

  double MicromagBEMElementEquations::
  green_normal_derivative_opt(const Vector<double>& x,
                              const Vector<double>& y,
                              const Vector<double>& n) const
  {
    // Get dimensions
    const unsigned node_dim = nodal_dimension();

#ifdef PARANOID
    // Check that n is a unit vector
    if(std::abs(two_norm(n) - 1.0) > 1e-10)
      {
        throw OomphLibError("n is not a unit vector",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Calculate the distance r from x to y.
    double norm_r = 0;
    double r[3]; // max possible DIM is 3 right?
    for(unsigned j=0; j<node_dim; j++)
      {
        r[j] = y[j] - x[j];
        norm_r += std::pow(r[j], 2);
      }
    norm_r = std::sqrt(norm_r);

    //??ds if x is a potentially singular point CHEAT
    //??ds could go horribly wrong, probably fine so long as not actually singular
    //??ds only works in 2d anyway

    // This is approximately true because the r and n unit vectors
    // become perpendicular as x approaches y. Hence n.r = 0 and
    // the whole function is zero.
    if(norm_r < 1e-6) return 0;

    double ndotr = 0;
    for(unsigned j=0; j<node_dim; j++)
      {
        // Calculate dot product of n with the unit vector in direction r.
        double rj_unit = r[j] / norm_r;
        ndotr += (n[j] * rj_unit);
      }

    // See write up for details of calculation.
    int mdm1 = - (int(node_dim) - 1);
    double result = (-1.0/Pi) * ndotr * std::pow((2.0*norm_r), mdm1);


    // ??ds if we had c++11 we could do a really nice optimisation here
    // which according to cachegrind would ~double speed of building the
    // BEM. Set node_dim to a const_expr and the calls to std::pow(...)
    // optimise out to either x*x or just x. Alternatively we could
    // template these elements, but then we lose some nice high level
    // stuff.

#ifdef PARANOID
    if(green_normal_derivative(x,y,n) != result)
      {
        std::string error_msg = "Optimised version gave different value to normal version!";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    return result;
  }

  /*
    hybrid FEM/BEM method from:
    T. R. Koehler, D. R. Fredkin, "Finite Element Methods for
    Micromagnetics", IEEE Trans. Magn. 28 (1992) 1239-1244

    integrations for matrix elements calculated according to:
    D. A. Lindholm, "Three-Dimensional Magnetostatic Fields from
    Point-Matched Integral Equations with Linearly Varying Scalar
    Sources", IEEE Trans. Magn. MAG-20 (1984) 2025-2032.

    Code is slightly transformed from magpar.
  */
  int Bele(const std::vector<double>& bvert, const std::vector<double>& facv1,
           const std::vector<double>& facv2, const std::vector<double>& facv3,
           std::vector<double>& matele)
  {
    using namespace magpar;
    using namespace MathematicalConstants;

    // Definitions from magpar headers
    const unsigned ND = 3; /**< space dimensions (no of cartesian coordinates) */
    const unsigned D_EPS = 1e-15*100; /**< threshold for equality of two real numbers */
    // const unsigned MU0 = 12.566371e-7; /**< mu0 = 4*M_PI*1e-7 = 12.566371e-7 Tm/A (=Vs/Am) */
    // const unsigned GAMMA = 2.210173e5; /**< gamma = mu0*g*|e|/(2*me) [m/As] = 1.758799e+11 [1/(Ts)] (cf. Diplomarbeit Scholz S. 14) */

    std::vector<double> rr(ND,0.0), zeta(ND,0.0);
    std::vector<double> rho1(ND,0.0),rho2(ND,0.0),rho3(ND,0.0);
    std::vector<double> eta1(ND,0.0),eta2(ND,0.0),eta3(ND,0.0);
    std::vector<double> xi1(ND,0.0),xi2(ND,0.0),xi3(ND,0.0);
    std::vector<double> gamma1(ND,0.0),gamma2(ND,0.0),gamma3(ND,0.0);
    std::vector<double> p(ND,0.0);
    double a(0.0), omega(0.0), eta1l,eta2l,eta3l, s1,s2,s3, rho1l,rho2l,rho3l, zetal;

    matele[0]=matele[1]=matele[2]=0.0;

    /* get coordinates of face's vertices */
    my_dcopy(ND,facv1,1,rho1,1);
    my_dcopy(ND,facv2,1,rho2,1);
    my_dcopy(ND,facv3,1,rho3,1);

    /* calculate edge vectors and store them in xi_j */
    my_dcopy(ND,rho2,1,xi1,1);
    my_daxpy(ND,-1.0,rho1,1,xi1,1);
    my_dcopy(ND,rho3,1,xi2,1);
    my_daxpy(ND,-1.0,rho2,1,xi2,1);
    my_dcopy(ND,rho1,1,xi3,1);
    my_daxpy(ND,-1.0,rho3,1,xi3,1);

    /* calculate zeta direction */
    douter(ND,xi1,xi2,zeta);

    /* calculate area of the triangle */
    zetal=my_dnrm2(ND,zeta,1);
    a=0.5*zetal;

    /* renorm zeta */
    my_dscal(ND,1.0/zetal,zeta,1);

    /* calculate s_j and normalize xi_j */
    s1=my_dnrm2(ND,xi1,1);
    my_dscal(ND,1.0/s1,xi1,1);
    s2=my_dnrm2(ND,xi2,1);
    my_dscal(ND,1.0/s2,xi2,1);
    s3=my_dnrm2(ND,xi3,1);
    my_dscal(ND,1.0/s3,xi3,1);

    douter(ND,zeta,xi1,eta1);
    douter(ND,zeta,xi2,eta2);
    douter(ND,zeta,xi3,eta3);

    gamma1[0]=gamma3[1]=my_ddot(ND,xi2,1,xi1,1);
    gamma1[1]=my_ddot(ND,xi2,1,xi2,1);
    gamma1[2]=gamma2[1]=my_ddot(ND,xi2,1,xi3,1);

    gamma2[0]=gamma3[2]=my_ddot(ND,xi3,1,xi1,1);
    gamma2[2]=my_ddot(ND,xi3,1,xi3,1);

    gamma3[0]=my_ddot(ND,xi1,1,xi1,1);

    /* get R=rr, the location of the source vertex (and the singularity) */
    rr=bvert;

    // If very close to the current surface element (triangle) then result is zero, done.
    double d = PointFromPlane(rr,rho1,rho2,rho3);
    if (fabs(d)<D_EPS) return(0);

    /* calculate rho_j */
    my_daxpy(ND,-1.0,rr,1,rho1,1);
    my_daxpy(ND,-1.0,rr,1,rho2,1);
    my_daxpy(ND,-1.0,rr,1,rho3,1);

    /* zetal gives ("normal") distance of R from the plane of the triangle */
    zetal=my_ddot(ND,zeta,1,rho1,1);

    /* skip the rest if zetal==0 (R in plane of the triangle)
       -> omega==0 and zetal==0 -> matrix entry=0
    */
    if (std::abs(zetal)<=D_EPS) {
      return(0);
    }

    rho1l=my_dnrm2(ND,rho1,1);
    rho2l=my_dnrm2(ND,rho2,1);
    rho3l=my_dnrm2(ND,rho3,1);

    double t_nom,t_denom;
    t_nom=
      rho1l*rho2l*rho3l+
      rho1l*my_ddot(ND,rho2,1,rho3,1)+
      rho2l*my_ddot(ND,rho3,1,rho1,1)+
      rho3l*my_ddot(ND,rho1,1,rho2,1);
    t_denom=
      std::sqrt(2.0*
                (rho2l*rho3l+my_ddot(ND,rho2,1,rho3,1))*
                (rho3l*rho1l+my_ddot(ND,rho3,1,rho1,1))*
                (rho1l*rho2l+my_ddot(ND,rho1,1,rho2,1))
                );

    /* catch special cases where the argument of acos
       is close to -1.0 or 1.0 or even outside this interval

       use 0.0 instead of D_EPS?
       fixes problems with demag field calculation
       suggested by Hiroki Kobayashi, Fujitsu
    */
    if (t_nom/t_denom<-1.0) {
      omega=(zetal >= 0.0 ? 1.0 : -1.0)*2.0*Pi;
    }
    /* this case should not occur, but does - e.g. examples1/headfield */
    else if (t_nom/t_denom>1.0) {
      return(0);
    }
    else {
      omega=(zetal >= 0.0 ? 1.0 : -1.0)*2.0*std::acos(t_nom/t_denom);
    }

    eta1l=my_ddot(ND,eta1,1,rho1,1);
    eta2l=my_ddot(ND,eta2,1,rho2,1);
    eta3l=my_ddot(ND,eta3,1,rho3,1);

    p[0]=log((rho1l+rho2l+s1)/(rho1l+rho2l-s1));
    p[1]=log((rho2l+rho3l+s2)/(rho2l+rho3l-s2));
    p[2]=log((rho3l+rho1l+s3)/(rho3l+rho1l-s3));

    matele[0]=(eta2l*omega-zetal*my_ddot(ND,gamma1,1,p,1))*s2/(8.0*Pi*a);
    matele[1]=(eta3l*omega-zetal*my_ddot(ND,gamma2,1,p,1))*s3/(8.0*Pi*a);
    matele[2]=(eta1l*omega-zetal*my_ddot(ND,gamma3,1,p,1))*s1/(8.0*Pi*a);

    return(0);
  }


  bool MicromagBEMElementEquations::
  normals_match(const Vector<unsigned> &node_list) const
  {

#warning using a long, expensive test to check if nodes need to be swapped
    // ??Ds fix this!

    // Get global nodal positions
    Vector<Vector<double> > x_tn(3);
    for(unsigned l=0, nl=node_list.size(); l<nl; l++)
      {
        x_tn[l].assign(3,0.0);
        node_pt(node_list[l])->position(x_tn[l]);
      }


    // Only works in 3D and for triangles (3 nodes)
    const unsigned node_dim = 3;
    const unsigned n_node = 3;

    /* First some pre-calculations to get everything ready. */

    // Calculate the length of the triangle sides and the unit vectors along the
    // triangle sides.
    Vector<Vector<double> > side_direction(n_node);
    for(unsigned i=0; i < n_node; i++)
      {
        // Get the next node around the triangle (i.e. 0 -> 1 -> 2 -> 0 etc.).
        unsigned next_node = (i+1)%n_node;

        // Get the vector along this side (xi in Lindholm1984)
        vector_diff(x_tn[next_node],x_tn[i],side_direction[i]);
      }

    // Calculate the normal to triangle. Assuming flat element => same
    // everywhere so we can just take the cross product of any two vectors in
    // the plane. Use the first two edge vectors.
    Vector<double> normal(node_dim,0.0);
    cross(side_direction[0],side_direction[1],normal);

    //??ds
    Vector<double> oomph_normal(node_dim, 0.0);
    Vector<double> s(dim(), 0.3);
    outer_unit_normal(s, oomph_normal);

    return VectorOps::dot(oomph_normal, normal) > 0;
  }


  //======================================================================
  ///
  //======================================================================
  void MicromagBEMElementEquations::
  fill_in_be_contribution_analytic(DenseMatrix<double> &boundary_matrix) const
  {
#ifdef PARANOID
    // Check 3d
    if(nodal_dimension() !=3)
      {
        std::ostringstream err;
        err << "Analytic calculation of boundary matrix only works in 3D.";
        throw OomphLibError(err.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(nnode_1d() != 2)
      {
        std::ostringstream err;
        err << "Analytic calculation of boundary matrix only "
            << "works for linear (i.e. nnode_1d = 2) elements.";
        throw OomphLibError(err.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    // ??ds check no hanging nodes - not sure what to do there yet
#endif

    // List of the node numbers for the three nodes that we are taking to
    // be the three corners of the triangle.
    Vector<unsigned> node_list(3,0);
    node_list[0] = 0; node_list[1] = 1; node_list[2] = 2;

    //??ds should be possible to get rid of this!
    if(!(normals_match(node_list)))
      {
        std::swap(node_list[0], node_list[1]);
        // std::cout << "swapping" << std::endl;
      }

    if(!(normals_match(node_list)))
      {
        std::cout <<  "!!! still not the same normal!" << std::endl;
      }

    // Get global nodal positions
    Vector<Vector<double> > x_nds(3);
    for(unsigned l=0, nl=node_list.size(); l<nl; l++)
      {
        x_nds[l].assign(3,0.0);
        node_pt(node_list[l])->position(x_nds[l]);
      }

    // Evaluate the integrals and store in boundary_matrix
    analytic_integral_dgreendn_triangle(x_nds, node_list, boundary_matrix);

  }

  //======================================================================
  ///
  //??ds could pass in unit normal since same for all triangle sub-elements
  //======================================================================
  void MicromagBEMElementEquations::
  analytic_integral_dgreendn_triangle(const Vector<Vector<double> >& x_tn,
                                      const Vector<unsigned>& l,
                                      DenseMatrix<double>& boundary_matrix) const
  {
    // Only works in 3D and for triangles (3 nodes)
    const unsigned node_dim = 3;
    const unsigned n_node = 3;

    /* First some pre-calculations to get everything ready. */

    // Calculate the length of the triangle sides and the unit vectors along the
    // triangle sides.
    Vector<double> side_length(n_node,0.0);
    Vector<Vector<double> > side_direction(n_node);
    for(unsigned i=0; i < n_node; i++)
      {
        // Get the next node around the triangle (i.e. 0 -> 1 -> 2 -> 0 etc.).
        unsigned next_node = (i+1)%n_node;

        // Get the vector along this side (xi in Lindholm1984)
        vector_diff(x_tn[next_node],x_tn[i],side_direction[i]);

        // Get the length of this side (s in Lindholm1984)
        side_length[i] = VectorOps::two_norm(side_direction[i]);
      }

    // Calculate the normal to triangle. Assuming flat element => same
    // everywhere so we can just take the cross product of any two vectors in
    // the plane. Use the first two edge vectors.
    Vector<double> normal(node_dim,0.0);
    cross(side_direction[0],side_direction[1],normal);

    // Calculate area of triangle using the cross product of the (unnormalised)
    // side_direction vectors, already calculated since it is the normal.
    double area = VectorOps::two_norm(normal)/2.0;

    // Get unit vectors of the above now that we have the area.
    Vector<double> unit_normal = normal;
    normalise(unit_normal);
    for(unsigned i=0; i<n_node; i++) normalise(side_direction[i]);

    // Calculate gamma
    Vector< Vector<double> > gamma(n_node);
    for(unsigned i=0; i<3; i++)
      {
        gamma[i].assign(n_node,0.0); // Initialise gamma[i]
        unsigned next_node = (i+1)%n_node; // Get next triangle vertex
        for(unsigned j=0; j<n_node; j++)
          {
            gamma[i][j] = dot(side_direction[next_node], side_direction[j]);
          }
      }

    /* Now evaluate the integral for every boundary node and add to
       boundary matrix. */
    for(unsigned i_sn=0; i_sn < boundary_mesh_pt()->nnode(); i_sn++)
      {
        // Get position of this node
        Vector<double> x_sn(node_dim,0.0);
        boundary_mesh_pt()->node_pt(i_sn)->position(x_sn);

        // Calculate rho (vector from each node in the triangle to source node)
        // and it's length.
        Vector<Vector<double> > rho(n_node); Vector<double> rhol(n_node);
        for(unsigned i_nd=0; i_nd<n_node; i_nd++)
          {
            vector_diff(x_tn[i_nd],x_sn,rho[i_nd]);
            rhol[i_nd] = VectorOps::two_norm(rho[i_nd]);
          }

        // Calculate zeta: the distance between the element and the source node
        // in the normal direction to the element.
        double zeta = dot(unit_normal,rho[0]);

        // If source node is in the plane of the element (i.e. zeta ~ 0) then
        // n.r = 0, nothing to calculate or add so we can move on to the next
        // source node.
        double tol = 1e-8;
        if( std::abs(zeta) < tol ) continue;

        // Calculate "P" (see paper Lindholm1984) for each node in the triangle
        Vector<double> P(n_node,0.0);
        for(unsigned i=0; i<n_node; i++)
          {
            unsigned next_node = (i+1)%n_node;
            P[i] = std::log( (rhol[i] + rhol[next_node] + side_length[i])
                             /(rhol[i] + rhol[next_node] - side_length[i]) );
          }


        // Calculate the solid angle (see Lindholm 1984)
        // I think it's impossible for the denominator to be zero since that
        // would imply either a zero rho vector or two anti-parallel rho vectors
        // which only happens for a source within the triangle. Both of these
        // cases would give zeta < 1e-8.
        double ratio = (rhol[0]*rhol[1]*rhol[2]
                        + rhol[0] * dot(rho[1],rho[2])
                        + rhol[1] * dot(rho[0],rho[2])
                        + rhol[2] * dot(rho[1],rho[0]))
          /
          std::sqrt(2
                    *(rhol[1]*rhol[2] + dot(rho[1],rho[2]) )
                    *(rhol[1]*rhol[0] + dot(rho[1],rho[0]) )
                    *(rhol[0]*rhol[2] + dot(rho[0],rho[2]) )
                    );

        // Round-off errors can cause ratio to be out of range for acos so
        // we need to check it.
        if(ratio > 1.0) ratio = 1.0;
        else if(ratio < -1.0) ratio = -1.0;

        int sign = (zeta > 0.0 ? +1.0 : -1.0);
        double omega = sign * 2 * std::acos(ratio);

        // Calculate eta: the unit vector normal to the side. Use it to
        // calculate etal: the distance in the diretion of eta to the source
        // node.
        Vector<Vector<double> > eta(3); Vector<double> etal(3,0.0);
        for(unsigned i=0; i<n_node; i++)
          {
            eta[i].assign(3,0.0);
            cross(unit_normal,side_direction[i],eta[i]);
            etal[i] = dot(eta[i],rho[i]);
          }

        // //??temp calculate with magpar code for testing
        // Vector<double> magpar_contribution(3,0.0);
        // Bele(x_sn, x_tn[0], x_tn[1], x_tn[2], magpar_contribution);

        /* Now put it all together and add the contribution to the boundary
           element matrix */
        for(unsigned i_tn=0; i_tn < n_node; i_tn++)
          {
            unsigned next_node = (i_tn+1)%n_node;

            // Add contribution to the appropriate value in the boundary matrix
            double oomph_contribution = (side_length[next_node]/(8*Pi*area))
              *((etal[next_node] * omega) - (zeta * dot(gamma[i_tn],P)));

            // if(std::abs(oomph_contribution - magpar_contribution[i_tn]) > 1e-7)
            //   {
            //     std::cout << std::abs(oomph_contribution - magpar_contribution[i_tn])
            //               << std::endl;
            //   }

            boundary_matrix(l[i_tn],i_sn) += oomph_contribution;

            // Note that Lindholm does include the factor of 1/(4*pi) in his
            // operator. Also note that he does NOT include the solid angle
            // factor on surfaces, the solid angle factor above is unrelated.
          }

      }
  }


  //====================================================================
  // Force building of templates
  //====================================================================
  template class QMicromagBEMElement<2,2>;
  template class QMicromagBEMElement<2,3>;
  template class QMicromagBEMElement<3,2>;
  template class QMicromagBEMElement<3,3>;

  template class TMicromagBEMElement<2,2>;
  template class TMicromagBEMElement<2,3>;
  template class TMicromagBEMElement<3,2>;
  template class TMicromagBEMElement<3,3>;

}
