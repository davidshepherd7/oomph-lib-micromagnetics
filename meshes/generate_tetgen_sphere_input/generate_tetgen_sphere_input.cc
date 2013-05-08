//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.90. August 3, 2009.
//LIC//
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <vector>
#include <map>
#include <set>


using namespace std;

//=====================================================================
/// Class for the spherical inclusion
//=====================================================================
class SphereMeshInput
{
  //Storage for the radius
  double Radius;

  //Storage for the centre
  vector<double> Centre;

  //Storage for the refinement level
  unsigned N_level;

  //Storage for the points
  vector<vector<double> > Vertex;

  //Storage for the Faces
  vector<vector<unsigned> > Facet;

public:

  ///Constructor take centre, radius and number of levels of refinement
  SphereMeshInput(const double &radius,const double &x_cen=0.0,
                     const double &y_cen=0.0, const double &z_cen=0.0,
                     const unsigned &n_level=0)
  {
    //Set the radius
    Radius = radius;
    //Set the centre
    Centre.resize(3);
    Centre[0] = x_cen; Centre[1] = y_cen; Centre[2] = z_cen;
    //Set the number of refinement levels
    N_level = n_level;

    //Golden ratio
    const double phi = 0.5*(1.0 + sqrt(5.0));

    //Set up the basic data structure for an icosahedron
    Vertex.resize(12);

    Vertex[0].resize(3);
    Vertex[0][0] = 0.0;
    Vertex[0][1] = 1.0;
    Vertex[0][2] = phi;
    Vertex[1].resize(3);
    Vertex[1][0] = 0.0;
    Vertex[1][1] = -1.0;
    Vertex[1][2] = phi;
    Vertex[2].resize(3);
    Vertex[2][0] = 0.0;
    Vertex[2][1] = 1.0;
    Vertex[2][2] = -phi;
    Vertex[3].resize(3);
    Vertex[3][0] = 0.0;
    Vertex[3][1] = -1.0;
    Vertex[3][2] = -phi;
    Vertex[4].resize(3);
    Vertex[4][0] = 1.0;
    Vertex[4][1] = phi;
    Vertex[4][2] = 0.0;
    Vertex[5].resize(3);
    Vertex[5][0] = -1.0;
    Vertex[5][1] = phi;
    Vertex[5][2] = 0.0;
    Vertex[6].resize(3);
    Vertex[6][0] = 1.0;
    Vertex[6][1] = -phi;
    Vertex[6][2] = 0.0;
    Vertex[7].resize(3);
    Vertex[7][0] = -1.0;
    Vertex[7][1] = -phi;
    Vertex[7][2] = 0.0;
    Vertex[8].resize(3);
    Vertex[8][0] = phi;
    Vertex[8][1] = 0.0;
    Vertex[8][2] = 1.0;
    Vertex[9].resize(3);
    Vertex[9][0] = phi;
    Vertex[9][1] = 0.0;
    Vertex[9][2] = -1.0;
    Vertex[10].resize(3);
    Vertex[10][0] = -phi;
    Vertex[10][1] = 0.0;
    Vertex[10][2] = 1.0;
    Vertex[11].resize(3);
    Vertex[11][0] = -phi;
    Vertex[11][1] = 0.0;
    Vertex[11][2] = -1.0;

    //Set up the connectivity
    Facet.resize(20);
    Facet[0].resize(3);
    Facet[0][0] = 0;
    Facet[0][1] = 1;
    Facet[0][2] = 8;

    Facet[1].resize(3);
    Facet[1][0] = 0;
    Facet[1][1] = 10;
    Facet[1][2] = 1;

    Facet[2].resize(3);
    Facet[2][0] = 0;
    Facet[2][1] = 5;
    Facet[2][2] = 10;

    Facet[3].resize(3);
    Facet[3][0] = 0;
    Facet[3][1] = 4;
    Facet[3][2] = 5;

    Facet[4].resize(3);
    Facet[4][0] = 0;
    Facet[4][1] = 8;
    Facet[4][2] = 4;

    Facet[5].resize(3);
    Facet[5][0] = 5;
    Facet[5][1] = 11;
    Facet[5][2] = 10;

    Facet[6].resize(3);
    Facet[6][0] = 5;
    Facet[6][1] = 2;
    Facet[6][2] = 11;

    Facet[7].resize(3);
    Facet[7][0] = 4;
    Facet[7][1] = 2;
    Facet[7][2] = 5;

    Facet[8].resize(3);
    Facet[8][0] = 4;
    Facet[8][1] = 9;
    Facet[8][2] = 2;

    Facet[9].resize(3);
    Facet[9][0] = 8;
    Facet[9][1] = 9;
    Facet[9][2] = 4;

    Facet[10].resize(3);
    Facet[10][0] = 6;
    Facet[10][1] = 9;
    Facet[10][2] = 8;

    Facet[11].resize(3);
    Facet[11][0] = 1;
    Facet[11][1] = 6;
    Facet[11][2] = 8;

    Facet[12].resize(3);
    Facet[12][0] = 1;
    Facet[12][1] = 7;
    Facet[12][2] = 6;

    Facet[13].resize(3);
    Facet[13][0] = 10;
    Facet[13][1] = 7;
    Facet[13][2] = 1;

    Facet[14].resize(3);
    Facet[14][0] = 10;
    Facet[14][1] = 11;
    Facet[14][2] = 7;

    Facet[15].resize(3);
    Facet[15][0] = 11;
    Facet[15][1] = 3;
    Facet[15][2] = 7;

    Facet[16].resize(3);
    Facet[16][0] = 7;
    Facet[16][1] = 3;
    Facet[16][2] = 6;

    Facet[17].resize(3);
    Facet[17][0] = 6;
    Facet[17][1] = 3;
    Facet[17][2] = 9;

    Facet[18].resize(3);
    Facet[18][0] = 9;
    Facet[18][1] = 3;
    Facet[18][2] = 2;

    Facet[19].resize(3);
    Facet[19][0] = 2;
    Facet[19][1] = 3;
    Facet[19][2] = 11;

    //Snap all vertices onto the sphere of the given
    //radius.
    {
      unsigned n_vertex = Vertex.size();
      for(unsigned n=0;n<n_vertex;n++)
        {
          double r2 = 0.0;
          for(unsigned i=0;i<3;i++)
            {
              r2 += (Vertex[n][i]*Vertex[n][i]);
            }
          for(unsigned i=0;i<3;i++)
            {
              //Scale to new radius
              Vertex[n][i] *= Radius/sqrt(r2);
            }
        }
    }


    //If we are refining then we need to split all the edges and
    //construct centre nodes and then the new facets
    for(unsigned l=0;l<N_level;++l)
      {
        //Find out the existing number of vertices and facets
        unsigned vertex_counter=Vertex.size();
        unsigned n_facet = Facet.size();
        unsigned facet_counter=0;

        //Storage for the new facets
        vector<vector<unsigned> > new_facet;
        //Local storage for the edge nodes index by the set of vertices
        //bounding the edge. The set will be sorted, so this gives a unique
        //key for each edge
        std::map<std::set<unsigned>, unsigned> edge_split;
        //Local storage for the key to the mape
        std::set<unsigned> edge_index;
        //Storage for a new vertex
        vector<double> new_vertex(3);
        //Storage for the indices of  edge vertices
        vector<unsigned> edge_vertex_index(3);

        //Loop over all the facets
        for(unsigned f=0;f<n_facet;f++)
          {
            //Loop over the edges
            for(unsigned e=0;e<3;e++)
              {
                //Set up the indices of the vertices bounding each edge
                unsigned vertex_index[2]={Facet[f][e],Facet[f][(e+1)%3]};
                //Clear the old set and then insert the indices
                edge_index.clear();
                edge_index.insert(vertex_index[0]);
                edge_index.insert(vertex_index[1]);

                //Get the edge vertex index from the map
                edge_vertex_index[e] = edge_split[edge_index];
                //If the edge has not been split, split it and add the vertex
                if(edge_vertex_index[e]==0)
                  {
                    for(unsigned i=0;i<3;i++)
                      {
                        new_vertex[i] = 0.5*(Vertex[vertex_index[0]][i] +
                                             Vertex[vertex_index[1]][i]);
                      }
                    //Add the vertex to the vector
                    Vertex.push_back(new_vertex);
                    //Store its index for the edge
                    edge_split[edge_index] = vertex_counter;
                    //and as the edge index
                    edge_vertex_index[e] = vertex_counter;
                    ++vertex_counter;
                  }
              }

            //We have now split all the edges, so output the facet information
            new_facet.resize(facet_counter+4);
            {
              new_facet[facet_counter].resize(3);
              new_facet[facet_counter][0] = Facet[f][0];
              new_facet[facet_counter][1] = edge_vertex_index[0];
              new_facet[facet_counter][2] = edge_vertex_index[2];
              ++facet_counter;

              new_facet[facet_counter].resize(3);
              new_facet[facet_counter][0] = edge_vertex_index[0];
              new_facet[facet_counter][1] = Facet[f][1];
              new_facet[facet_counter][2] = edge_vertex_index[1];
              ++facet_counter;


              new_facet[facet_counter].resize(3);
              new_facet[facet_counter][0] = edge_vertex_index[1];
              new_facet[facet_counter][1] = Facet[f][2];
              new_facet[facet_counter][2] = edge_vertex_index[2];
              ++facet_counter;

              new_facet[facet_counter].resize(3);
              new_facet[facet_counter][0] = edge_vertex_index[0];
              new_facet[facet_counter][1] = edge_vertex_index[1];
              new_facet[facet_counter][2] = edge_vertex_index[2];
              ++facet_counter;
            }
          } //End of split of the facets

        //Now copy over the new facets into the old facets
        Facet.resize(facet_counter);
        for(unsigned f=0;f<facet_counter;++f)
          {
            Facet[f] = new_facet[f];
          }

        //Snap all vertices onto the sphere of given
        //radius after the level has been
        //completed
        unsigned n_vertex = Vertex.size();
        for(unsigned n=0;n<n_vertex;n++)
          {
            double r2 = 0.0;
            for(unsigned i=0;i<3;i++)
              {
                //and calculation the contribution to r squared
                r2 += (Vertex[n][i]*Vertex[n][i]);
              }
            //Then move it back
            for(unsigned i=0;i<3;i++)
              {
                //Scale to new radius
                Vertex[n][i] *= Radius/sqrt(r2);
              }
          }
      }

    //Finally, translate to the appropriate centre position
    {
      unsigned n_vertex = Vertex.size();
      for(unsigned n=0;n<n_vertex;n++)
        {
          for(unsigned i=0;i<3;i++) {Vertex[n][i] += Centre[i];}
        }
    }


  }

  //Write out the data structure
  void write_facets(vector<vector<double> > &vertex,
                    vector<vector<unsigned> > &facet)
  {
    vertex = Vertex;
    facet = Facet;
  }

};


int main(int argc, char* argv[])
{

  double radius;
  unsigned refinement_level;

  std::string usage("Usage: ./generate_tetgen_sphere_input radius refinement_level");
    usage += " > your_output_file";

  if (argc == 3)
    {
      radius = std::atof(argv[1]);
      refinement_level = std::atoi(argv[2]);
    }
  else
    {
      std::cerr << "Wrong number of inputs, should be two"
                << " but got " << argc - 1 << "."
                << std::endl
                << usage
                << std::endl;
      return 1;
    }

  // Make a sphere centred at the origin
  SphereMeshInput sphere(radius, 0.0, 0.0, 0.0, refinement_level);

  vector<vector<double> > points;
  vector<vector<unsigned> > facets;

  sphere.write_facets(points, facets);

  // Output the points to file
  {
    // First line: [number of nodes] [dimension = 3] [number of attributes = 0]
    // [number of boundary markers = 1]
    std::cout << points.size() << " " << 3 << " 0 1" <<std::endl;


    // Remaining lines: [node index] [x] [y] [z] [[attributes]] [[boundary marker]]
    for(unsigned i_pt=0; i_pt < points.size(); i_pt++)
      {
        std::cout << i_pt + 1 << " " // index of node
                << points[i_pt][0] << " "
                << points[i_pt][1] << " "
                << points[i_pt][2] << " "
                << "1" // all nodes are on the same boundary (boundary number 1)
                << std::endl;
      }


    // Output the facets to file

    // First line: <# of facets> <boundary markers (0 or 1)>
    std::cout << facets.size() << " " << "1" <<std::endl;

    for(unsigned i_pt=0; i_pt < facets.size(); i_pt++)
      {
        // Remaining lines:  <# of polygons> [# of holes] [boundary marker]
        std::cout << "1 0 1" << std::endl;

        //each polygon: <# of corners> <corner 1> <corner 2> ... <corner #>
        //              ...
        std::cout << "3 "
                << facets[i_pt][0] +1 << " "
                << facets[i_pt][1] +1<< " "
                << facets[i_pt][2] +1<< " "
                << std::endl;
      }

    std::cout << "0" <<std::endl;
    // // Output list of holes: don't mesh the external region
    // std::cout << "1" <<std::endl;
    // std::cout << "1 " << 3*r << " " << 3*r << " " << 3*r << std::endl;
  }

  return 0;
}
