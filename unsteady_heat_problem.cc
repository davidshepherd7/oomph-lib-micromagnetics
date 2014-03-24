// Meshes for mesh factory
#include "../../src/meshes/simple_rectangular_quadmesh.h"
#include "../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../src/meshes/simple_cubic_tet_mesh.h"
#include "../../src/meshes/simple_cubic_mesh.h"
#include "../../src/meshes/tetgen_mesh.h"
#include "../../src/meshes/triangle_mesh.h"

#include "./multi_mesh.h"
#include "./single_element_mesh.h"

#include "unsteady_heat_problem.h"

void UnsteadyHeatProblem::build(Vector<Mesh*>& bulk_mesh_pts)
{

  // Call the underlying build to deal with adding meshes and time stepper
  MyProblem::build(bulk_mesh_pts);

  // Set the boundary conditions for this problem:
  for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
    {
      unsigned n_bound = mesh_pt(msh)->nboundary();
      for(unsigned b=0; b<n_bound; b++)
        {
          unsigned n_node = mesh_pt(msh)->nboundary_node(b);
          for (unsigned n=0;n<n_node;n++)
            {
              // Pinned to initial condition which is exact solution
              mesh_pt(msh)->boundary_node_pt(b,n)->pin(0);
            }
        }
    }


  // Complete the build of all elements so they are fully functional
  //----------------------------------------------------------------

  // Find number of elements in mesh
  // Set the boundary conditions for this problem:
  for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
    {
      unsigned n_element = mesh_pt(msh)->nelement();

      // Loop over the elements to set up element-specific
      // things that cannot be handled by constructor
      for(unsigned i=0;i<n_element;i++)
        {
          // Upcast from FiniteElement to the present element
          UnsteadyHeatEquationsBase *el_pt
            = dynamic_cast<UnsteadyHeatEquationsBase*>(mesh_pt(msh)->element_pt(i));

          //Set the source function pointer
          el_pt->source_fct_pt() = Source_fct_pt;
        }
    }


  // Finish building
  // ============================================================

  // Build the global mesh
  this->build_global_mesh();

  // Number the equations
  this->assign_eqn_numbers();


  // Write out some stuff
  oomph_info << "Number of equations: " << ndof() << std::endl;
  oomph_info << "Number of sub meshes: " << this->nsub_mesh() << std::endl;

}

namespace UnsteadyHeatFactories
{
  /// \short Make a mesh as specified by an input argument. Refined
  /// according to the given refinement level (in some way appropriate
  /// for that mesh type). Assumption: this will be passed into a
  /// problem, which will delete the pointer when it's done.
  Mesh* mesh_factory(const std::string& _mesh_name,
                     int refinement_level,
                     TimeStepper* time_stepper_pt,
                     double scaling_factor,
                     unsigned nnode1d)
  {
    // Ignore case in mesh names
    const std::string mesh_name = StringConversion::to_lower(_mesh_name);

    // Refinement always roughly the same for structured meshes
    unsigned nx = 5 * std::pow(2, refinement_level-1);

    // Make the mesh and store a pointer to it
    Mesh* mesh_pt = 0;
    if(mesh_name == "sq_square" && nnode1d == 2)
      {
        double lx = 1.0;
        mesh_pt = new SimpleRectangularQuadMesh<QUnsteadyHeatElement<2,2> >
          (nx, nx, lx, lx, time_stepper_pt);
      }
    else if(mesh_name == "st_square" && nnode1d == 2)
      {
        double lx = 1.0;
        mesh_pt = new SimpleRectangularTriMesh<TUnsteadyHeatElement<2,2> >
          (nx, nx, lx, lx, time_stepper_pt);

        mesh_pt->setup_boundary_element_info();

        // Turn off triangle refinement dump stuff (breaks Micromag
        // elements).
        checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
          disable_triangulateio_restart();
      }
    else if(mesh_name == "single-element" && nnode1d == 2)
      {
        mesh_pt = new SingleElementMesh<QUnsteadyHeatElement<2,2> >(time_stepper_pt);
      }
    else if(mesh_name == "ut_square" && nnode1d == 2)
      {
        mesh_pt = new TriangleMesh<TUnsteadyHeatElement<2, 2> >
          ("./meshes/square." + to_string(refinement_level) + ".node",
           "./meshes/square." + to_string(refinement_level) + ".ele",
           "./meshes/square." + to_string(refinement_level) + ".poly",
           time_stepper_pt);

        // Turn off triangle refinement dump stuff (breaks Micromag
        // elements).
        checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
          disable_triangulateio_restart();
      }
    else if(mesh_name == "st_cubeoid" && nnode1d == 2)
      {
        // nmag cubeoid
        double lx = 1, ly = lx, lz = 3*lx;
        unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
        mesh_pt = new SimpleCubicTetMesh<TUnsteadyHeatElement<3, 2> >
          (nx, ny, nz, lx, ly, lz, time_stepper_pt);

        mesh_pt->setup_boundary_element_info();
      }
    else if(mesh_name == "sqt_cubeoid" && nnode1d == 2)
      {
        Mesh* qmesh_pt = mesh_factory("sq_cubeoid", refinement_level,
                                      time_stepper_pt,
                                      1, nnode1d);
        //??ds memory leak, fix? Can't delete this mesh or nodes will
        //go...

        TetMeshBase* tmesh_pt = new TetMeshBase;
        ElementFactoryFctPt factory_fpt =
          MeshCreationHelpers::new_element<TUnsteadyHeatElement<3, 2> >;

        MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

        mesh_pt = tmesh_pt;
      }
    else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
      {
        mesh_pt = new TetgenMesh<TUnsteadyHeatElement<3, 2> >
          ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
           "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
           "./meshes/cubeoid." + to_string(refinement_level) + ".face",
           time_stepper_pt);
      }
    else if(mesh_name == "ut_mumag4" && nnode1d == 2)
      {
        mesh_pt = new TetgenMesh<TUnsteadyHeatElement<3, 2> >
          ("./meshes/mumag4." + to_string(refinement_level) + ".node",
           "./meshes/mumag4." + to_string(refinement_level) + ".ele",
           "./meshes/mumag4." + to_string(refinement_level) + ".face",
           time_stepper_pt);
      }
    else if(mesh_name == "st_mumag4" && nnode1d == 2)
      {
        mesh_pt = new SimpleCubicTetMesh<TUnsteadyHeatElement<3, 2> >
          (5*nx, std::ceil(1.25*nx), 1, 500, 125, 3, time_stepper_pt);

        mesh_pt->setup_boundary_element_info();
      }
    else if(mesh_name == "sq_mumag4" && nnode1d == 2)
      {
        unsigned this_nx = refinement_level;

        mesh_pt = new SimpleCubicMesh<QUnsteadyHeatElement<3, 2> >
          (5*this_nx, std::ceil(1.25*this_nx), 2, 500, 125, 3, time_stepper_pt);

        mesh_pt->setup_boundary_element_info();
      }
    else if(mesh_name == "sqt_mumag4" && nnode1d == 2)
      {
        Mesh* qmesh_pt = mesh_factory("sq_mumag4", refinement_level,
                                      time_stepper_pt,
                                      1, nnode1d);
        //??ds memory leak, fix? Can't delete this mesh or nodes will
        //go...

        // Convert to tet mesh
        TetMeshBase* tmesh_pt = new TetMeshBase;
        ElementFactoryFctPt factory_fpt =
          MeshCreationHelpers::new_element<TUnsteadyHeatElement<3, 2> >;
        MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

        mesh_pt = tmesh_pt;
      }
    else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
      {
        double lx = 1, ly = lx, lz = 3*lx;
        mesh_pt = new SimpleCubicMesh<QUnsteadyHeatElement<3, 2> >
          (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
      }
    else if(mesh_name == "ut_sphere" && nnode1d == 2)
      {
        mesh_pt = new TetgenMesh<TUnsteadyHeatElement<3, 2> >
          ("./meshes/sphere." + to_string(refinement_level) + ".node",
           "./meshes/sphere." + to_string(refinement_level) + ".ele",
           "./meshes/sphere." + to_string(refinement_level) + ".face",
           time_stepper_pt);
      }
    else
      {
        throw OomphLibError("Unrecognised mesh name " + mesh_name,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    // Scale the mesh as requested
    scale_mesh(scaling_factor, mesh_pt);

    // Done: pass out the mesh pointer
    return mesh_pt;
  }

}
