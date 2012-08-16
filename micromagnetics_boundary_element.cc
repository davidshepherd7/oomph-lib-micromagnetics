/*
  description of file goes here
*/


using namespace oomph;
using namespace MathematicalConstants;


namespace oomph
{

  /// Assign the static boundary mesh pointer for the bem elements. We have to do
  /// this seperately to the class because it is a static variable.
  template<class ELEMENT,unsigned DIM>
  Mesh* MicromagFaceElement<ELEMENT,DIM>::Boundary_mesh_pt=0;

} // End of oomph namespace
