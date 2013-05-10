#ifndef OOMPH_MAGNETOSTATIC_FIELD_FLUX_ELEMENT_H
#define OOMPH_MAGNETOSTATIC_FIELD_FLUX_ELEMENT_H

/*
  description of file goes here
*/

#include "generic.h"
#include "./template_free_poisson_flux.h"

using namespace oomph;

namespace oomph
{

  // ============================================================
  /// Poisson flux element with source function overridden to get m.n from
  /// magnetic face elements (hacky but oh well). Note that Poisson elements
  /// actually invert the sign of the residuals (and Jacobian) as compared to
  /// my derivation but this shouldn't matter.
  // ============================================================
  class MagnetostaticFieldFluxEquations :
    public TFPoissonFluxEquations
  {
  public:

    /// Real constructor
    MagnetostaticFieldFluxEquations(FiniteElement* const &bulk_el_pt,
                                    const int &face_index)
      : TFPoissonFluxEquations(bulk_el_pt, face_index)
    {}

    /// Destructor
    ~MagnetostaticFieldFluxEquations() {}

    /// Get flux from both mdotn and supplied flux function.
    void get_elemental_flux(const Vector<double> &face_s, double &flux) const
    {
      // Get contribution from normal component of magnetisation.
      flux = mdotn(face_s);
    }

    /// Calculate the dot product of m with unit normal
    double mdotn(const Vector<double> &face_s) const
    {
      Vector<double> normal(this->nodal_dimension(),0.0);
      this->outer_unit_normal(face_s,normal);

      // Cast bulk element to MagnetostaticFieldEquations
      MagnetostaticFieldEquations* field_ele_pt =
        checked_dynamic_cast<MagnetostaticFieldEquations*>(this->bulk_element_pt());

      // Get magnetisation from bulk magnetics element.
      unsigned dim = this->nodal_dimension();
      Vector<double> m, bulk_s(dim,0.0);
      this->get_local_coordinate_in_bulk(face_s, bulk_s);
      field_ele_pt->micromag_element_pt()->interpolated_m_micromag(bulk_s,m);
      // ??ds we have assumed (again) that the poisson and magnetic
      //elements are in the same places. This is pretty stupid :(

      // Take dot product up to the dimension of the unit normal (m could
      // have a higher dimension since we always have 3 components of m. If
      // so then the dot product contribution of the extra components is
      // zero because the component of the normal in that direction is
      // obviously zero.).
      double mdotn = 0;
      for(unsigned j=0, nj=normal.size(); j<nj; j++)
        {
          mdotn += m[j] * normal[j];
        }

      return mdotn;
    }


  private:

    /// Inacessible default constructor
    MagnetostaticFieldFluxEquations() {}

    /// Inaccessible copy constructor
    MagnetostaticFieldFluxEquations(const MagnetostaticFieldFluxEquations &dummy)
    {BrokenCopy::broken_copy("MagnetostaticFieldFluxEquations");}

    /// Inaccessible assignment operator
    void operator=(const MagnetostaticFieldFluxEquations &dummy)
    {BrokenCopy::broken_assign("MagnetostaticFieldFluxEquations");}
  };


  /// \short Class to pull in the geometrical element code. Bit strange to
  /// use the weird templated FaceGeometry class but it seems to work.
  template<unsigned DIM, unsigned NNODE_1D>
  class QMagnetostaticFieldFluxElement :
    public virtual FaceGeometry<class QElement<DIM, NNODE_1D> >,
    public virtual MagnetostaticFieldFluxEquations
  {
    public:
    /// Real constructor
    QMagnetostaticFieldFluxElement(FiniteElement* const &bulk_el_pt,
                                   const int &face_index)
      : FaceGeometry<class QElement<DIM, NNODE_1D> >(),
        MagnetostaticFieldFluxEquations(bulk_el_pt, face_index)
    {}
  };

  /// \short Class to pull in the geometrical element code. Bit strange to
  /// use the weird templated FaceGeometry class but it seems to work.
  template<unsigned DIM, unsigned NNODE_1D>
  class TMagnetostaticFieldFluxElement :
    public virtual FaceGeometry<class TElement<DIM, NNODE_1D> >,
    public virtual MagnetostaticFieldFluxEquations
  {
    public:
    /// Real constructor
    TMagnetostaticFieldFluxElement(FiniteElement* const &bulk_el_pt,
                                   const int &face_index)
      : FaceGeometry<class TElement<DIM, NNODE_1D> >(),
        MagnetostaticFieldFluxEquations(bulk_el_pt, face_index)
    {}
  };


} // End of oomph namespace

#endif
