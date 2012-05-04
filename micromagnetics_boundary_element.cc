
//??ds this file should have more in it eventually...

/// Assign the static boundary mesh pointer for the bem elements. We have to do
/// this seperately to the class because it is a static variable.
template<class ELEMENT,unsigned DIM>
Mesh* MicromagFaceElement<ELEMENT,DIM>::Boundary_mesh_pt=0;
