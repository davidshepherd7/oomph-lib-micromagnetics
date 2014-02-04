
#include "../../src/generic/double_vector.h"
#include "../../src/generic/Vector.h"
#include "../../src/generic/matrices.h"
#include "../../src/generic/oomph_utilities.h"

#include "../../src/generic/explicit_timesteppers.h"

#include "cayley_integrators.h"
#include "vector_helpers.h"


namespace oomph
{
  // safe to have th below functioins in global oomph namespace because
  // they are not in the header so not callable outside this file.

  void cayley(const DoubleVector& v, DenseDoubleMatrix& out)
  {
#ifdef PARANOID
    if(v.nrow() != 3)
      {
        std::string err = "Only for 3D vectors";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    out.resize(3, 3);

    // From wikipedia! http://en.wikipedia.org/wiki/Cayley_transform factor
    // of 1/2 comes in from application to llg, see e.g. Efficiency of the
    // Geometric Integration of Landau–Lifshitz–Gilbert Equation Based on
    // Cayley Transform by Oriano Bottauscio and Alessandra Manzin
    const double w = 1, x=v[0]/2, y=v[1]/2, z=v[2]/2;
    const double K = w*w + x*x + y+y +z*z;

    out(0, 0) = (w*w + x*x - y*y - z*z)/K;
    out(0, 1) = (2*(x*y - w*z))/K;
    out(0, 2) = (2*(w*y + x*z))/K;
    out(1, 0) = (2*(x*y + w*z))/K;
    out(1, 1) = (w*w - x*x + y*y - z*z)/K;
    out(1, 2) = (2*(y*z - w*x))/K;
    out(2, 0) = (2*(x*z - w*y))/K;
    out(2, 1) = (2*(w*x + y*z))/K;
    out(2, 2) = (w*w - x*x - y*y + z*z)/K;
  }


  DoubleVector derivs_to_omega(const DoubleVector& global_m,
                               const DoubleVector& global_dmdt)
  {
    DoubleVector global_omega(global_dmdt.distribution_pt());

    LinearAlgebraDistribution local_dist(0, 3, false);

    for(unsigned k = 0, nk = global_m.nrow(); k < nk; k += 3)
      {
        // pull out m and dmdt
        DoubleVector local_dmdt(local_dist), local_m(local_dist);
        for(unsigned j=0; j<3; j++) {local_m[j] = global_m[k+j];}
        for(unsigned j=0; j<3; j++) {local_dmdt[j] = global_dmdt[k+j];}


        // This is very dodgy! We need h + damp (mxh) but all we have is
        // mdot = m x ( h + damp (mxh))... so we use the
        // sort-of-inverse-cross-product which should be ok because we only
        // use it in other cross products
        for(unsigned j=0; j<3; j++)
          {
            global_omega[k+j] = -VectorOps::opt_cross(j, local_dmdt.values_pt(),
                                                      local_m.values_pt());
          }
      }

    return global_omega;
  }

  DoubleVector cayley_step(const DoubleVector& global_m_n,
                           const DoubleVector& global_omega,
                           const double& dt)
  {
    DoubleVector global_m_np1(global_m_n.distribution_pt());

    LinearAlgebraDistribution local_dist(0, 3, false);
    for(unsigned k = 0, nk = global_m_n.nrow(); k < nk; k += 3)
      {
        // pull out m and omega
        DoubleVector local_m_n(local_dist), local_omega(local_dist);
        for(unsigned j=0; j<3; j++) {local_m_n[j] = global_m_n[k+j];}
        for(unsigned j=0; j<3; j++) {local_omega[j] = dt * global_omega[k+j];}

        // Take the time step
        DenseDoubleMatrix cayley_operator;
        DoubleVector local_m_np1(local_dist);
        cayley(local_omega, cayley_operator);
        cayley_operator.multiply(local_m_n, local_m_np1);

        // put back m
        for(unsigned j=0; j<3; j++) {global_m_np1[k+j] = local_m_np1[j];}

#ifdef PARANOID
        Vector<double> mvec(3);
        mvec[0] = local_m_n[0];
        mvec[1] = local_m_n[1];
        mvec[2] = local_m_n[2];
        if(std::abs(VectorOps::two_norm(mvec) - 1) > 1e-3)
          {
            std::ostringstream err;
            err << "Magnetisation not near unit length, possibly not actually the magnetisation?";
            err << "here it is: " << mvec << std::endl;
            err << "length is: " << VectorOps::two_norm(mvec) << std::endl;
            throw OomphLibError(err.str(), OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
#endif
      }

    return global_m_np1;
  }


  void CayleyEuler::timestep(ExplicitTimeSteppableObject* const &object_pt,
                           const double &dt)
  {
    object_pt->actions_before_explicit_timestep();
    object_pt->actions_before_explicit_stage();

    // Get values
    DoubleVector global_dmdt, global_m_n;
    object_pt->get_inverse_mass_matrix_times_residuals(global_dmdt);
    object_pt->get_dofs(global_m_n);
    DoubleVector global_omega_n = derivs_to_omega(global_m_n, global_dmdt);

    // Try to check that all dofs are magnetisation. Should always be 3
    // dofs per node so should be a multiple of 3.
#ifdef PARANOID
    if(global_m_n.nrow() % 3 != 0)
      {
        std::string err = "Only valid for purely magnetisation dofs!";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Take a step
    DoubleVector global_m_np1 = cayley_step(global_m_n, global_omega_n, dt);

    // Update object
    object_pt->set_dofs(global_m_np1);
    object_pt->time() += dt;

    object_pt->actions_after_explicit_stage();
    object_pt->actions_after_explicit_timestep();
  }


  void CayleyRK2::timestep(ExplicitTimeSteppableObject* const &object_pt,
                         const double &dt)
  {
    object_pt->actions_before_explicit_timestep();
    object_pt->actions_before_explicit_stage();

    // Get values
    DoubleVector global_dmdt_n, global_m_n;
    object_pt->get_inverse_mass_matrix_times_residuals(global_dmdt_n);
    object_pt->get_dofs(global_m_n);
    DoubleVector global_omega_n = derivs_to_omega(global_m_n, global_dmdt_n);

    // take intermediate step and update dofs
    DoubleVector global_m_star = cayley_step(global_m_n, global_omega_n, dt);
    object_pt->set_dofs(global_m_star);
    object_pt->time() += dt;

    // Get intermediate values
    DoubleVector global_dmdt_star;
    object_pt->get_inverse_mass_matrix_times_residuals(global_dmdt_star);
    DoubleVector global_omega_star = derivs_to_omega(global_m_star,
                                                     global_dmdt_star);

    // Get omega to use in real step
    DoubleVector global_omega_squiggle = global_omega_star;
    global_omega_squiggle += global_omega_n;
    global_omega_squiggle /= 2;


    // Take real step and update dofs
    DoubleVector global_m_np1 = cayley_step(global_m_n, global_omega_squiggle, dt);
    object_pt->set_dofs(global_m_np1);

    object_pt->actions_after_explicit_stage();
    object_pt->actions_after_explicit_timestep();
  }
}
