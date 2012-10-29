/*
    This file is part of magpar.

    Copyright (C) 2002-2010 Werner Scholz

    www:   http://www.magpar.net/
    email: magpar(at)magpar.net

    magpar is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    magpar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with magpar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* #include "field.h" */
/* #include "util/util.h" */

#include <cmath>
#include <vector>

#include <iostream>

// Include the appropriate version of the pretty print header depending on if we
// are using c++0x or not
#ifdef __GXX_EXPERIMENTAL_CXX0X__
  #include "../prettyprint.hpp"
#else
  #include "../prettyprint98.hpp"
#endif



/* Macros taken from magpar */

/** \name Defines of numbers
*/
/*@{*/
#define ND 3       /**< space dimensions (no of cartesian coordinates) */
#define NV 4       /**< number of vertices(=degrees of freedom) per element */
#define NF 4       /**< number of faces per element */
#define NN 3       /**< number of vertices per face */
#define NP 19      /**< number of property data */
/*@}*/


/** \name Defines of indicators
*/
/*@{*/
#define C_BND -1   /**< indicator for boundary node/face */
#define C_INT -2   /**< indicator for interior node/face */
#define C_UNK -4   /**< indicator for unknown state */
/*@}*/

#define D_EPS 1e-15*100 /**< threshold for equality of two real numbers */

/** \name Defines for physical constants
*/
/*@{*/
#define MU0 (12.566371e-7) /**< mu0 = 4*M_PI*1e-7 = 12.566371e-7 Tm/A (=Vs/Am) */
#define GAMMA (2.210173e5) /**< gamma = mu0*g*|e|/(2*me) [m/As] = 1.758799e+11 [1/(Ts)] (cf. Diplomarbeit Scholz S. 14) */
/*@}*/
#define PI (3.14159)

/* value of damping constant for which M is locked/rigid */
#define RIGID_M_ALPHA 999

/* replace all calls to cblas with simple macros here
   because we are just coping with ND=3 vectors
   thus, big overhead for calling real cblas functions!
*/
#define my_daxpy(a,b,c,d,e,f) {(e)[0]+=b*(c)[0];(e)[1]+=b*(c)[1];(e)[2]+=b*(c)[2];}
#define my_dcopy(a,b,c,d,e)   {(d)[0]=(b)[0];(d)[1]=(b)[1];(d)[2]=(b)[2];}
#define my_dnrm2(a,b,c)       std::sqrt((b)[0]*(b)[0]+(b)[1]*(b)[1]+(b)[2]*(b)[2])
#define my_dscal(a,b,c,d)     {(c)[0]*=b;(c)[1]*=b;(c)[2]*=b;}
#define my_ddot(a,b,c,d,e)    ((b)[0]*(d)[0]+(b)[1]*(d)[1]+(b)[2]*(d)[2])
#define douter(a,b,c,d)       {(d)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1];(d)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2];(d)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0];}



double PointFromPlane( const std::vector<double>& x, const std::vector<double>& v1,
		       const std::vector<double>& v2, const std::vector<double>& v3)
{
  std::vector<double>  ab(ND,0.0), ac(ND,0.0);
  std::vector<double>  n(ND,0.0);

  /*calculate edge vectors */
  my_dcopy(ND,v1,1,ab,1);
  my_daxpy(ND,-1.0,v2,1,ab,1);
  my_dcopy(ND,v1,1,ac,1);
  my_daxpy(ND,-1.0,v3,1,ac,1);

  /* calculate normal vector */
  douter(ND,ab,ac,n);

  /* calculate distance */
  return my_ddot(ND,x,1,n,1)-my_ddot(ND,v1,1,n,1);

  return(0);
}

/*
  hybrid FEM/BEM method from:
  T. R. Koehler, D. R. Fredkin, "Finite Element Methods for
  Micromagnetics", IEEE Trans. Magn. 28 (1992) 1239-1244

  integrations for matrix elements calculated according to:
  D. A. Lindholm, "Three-Dimensional Magnetostatic Fields from
  Point-Matched Integral Equations with Linearly Varying Scalar
  Sources", IEEE Trans. Magn. MAG-20 (1984) 2025-2032.
*/
int Bele(const std::vector<double>& bvert, const std::vector<double>& facv1,
	 const std::vector<double>& facv2, const std::vector<double>& facv3,
	 std::vector<double>& matele)
{
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

  std::cout << t_nom/t_denom << std::endl;
  /* catch special cases where the argument of acos
     is close to -1.0 or 1.0 or even outside this interval

     use 0.0 instead of D_EPS?
     fixes problems with demag field calculation
     suggested by Hiroki Kobayashi, Fujitsu
  */
  if (t_nom/t_denom<-1.0) {
    omega=(zetal >= 0.0 ? 1.0 : -1.0)*2.0*PI;
  }
  /* this case should not occur, but does - e.g. examples1/headfield */
  else if (t_nom/t_denom>1.0) {
    return(0);
  }
  else {
    omega=(zetal >= 0.0 ? 1.0 : -1.0)*2.0*std::acos(t_nom/t_denom);
  }

  // //ok
  // std::cout << omega << std::endl;

  eta1l=my_ddot(ND,eta1,1,rho1,1);
  eta2l=my_ddot(ND,eta2,1,rho2,1);
  eta3l=my_ddot(ND,eta3,1,rho3,1);

  // //ok
  // std::cout << eta1 << std::endl;
  // std::cout << eta2 << std::endl;
  // std::cout << eta3 << std::endl;

  p[0]=log((rho1l+rho2l+s1)/(rho1l+rho2l-s1));
  p[1]=log((rho2l+rho3l+s2)/(rho2l+rho3l-s2));
  p[2]=log((rho3l+rho1l+s3)/(rho3l+rho1l-s3));

  //ok std::cout << p << std::endl;

  // ok
  //   std::cout << zetal << std::endl;

  // //ok
  // std::cout << gamma1 << " " << gamma2 << " " << gamma3 << std::endl;

  // //ok
  // std::cout << a << std::endl;

  matele[0]=(eta2l*omega-zetal*my_ddot(ND,gamma1,1,p,1))*s2/(8.0*PI*a);
  matele[1]=(eta3l*omega-zetal*my_ddot(ND,gamma2,1,p,1))*s3/(8.0*PI*a);
  matele[2]=(eta1l*omega-zetal*my_ddot(ND,gamma3,1,p,1))*s1/(8.0*PI*a);

  return(0);
}


//======================================================================
/// Supporting functions for oomph-lib calculation
//======================================================================


inline double dot(const std::vector<double>& a, const std::vector<double>& b)
{
#ifdef PARANOID
  if (a.size() != b.size())
    throw OomphLibError("std::vectors must be the same length", "mod_diff",
			OOMPH_EXCEPTION_LOCATION);
#endif
  double temp = 0;
  for(unsigned i=0; i<a.size(); i++)
    temp += a[i] * b[i];
  return temp;
}


/// Calculate the cross product of vectors A and B, store the result in
/// vector output. NOTE: the cross product is only valid for 3-dimensional
/// vectors
void cross(std::vector<double>& A, std::vector<double>& B, std::vector<double>& output)
{
  output[0] = A[1]*B[2] - A[2]*B[1];
  output[1] = A[2]*B[0] - A[0]*B[2];
  output[2] = A[0]*B[1] - A[1]*B[0];
}


inline double mod(const std::vector<double>& a)
{
  return std::sqrt(dot(a,a));
}

inline void vector_diff(const std::vector<double>& a, const std::vector<double>& b,
			std::vector<double>& diff)
{
#ifdef PARANOID
  if (a.size() != b.size())
    throw OomphLibError("std::vectors must be the same length", "mod_diff",
			OOMPH_EXCEPTION_LOCATION);
#endif
  diff.assign(a.size(),0.0);
  for(unsigned i=0; i<a.size(); i++)
    diff[i] = a[i] - b[i];
}

inline void abs_vector_diff(const std::vector<double>& a, const std::vector<double>& b,
			    std::vector<double>& diff)
{
  vector_diff(a,b,diff);
  for(unsigned i=0; i<a.size(); i++)
    diff[i] = std::abs(diff[i]);
}

inline double mod_diff(const std::vector<double>& a, const std::vector<double>& b)
{
  std::vector<double> diff(a.size(),0.0);
  vector_diff(a,b,diff);
  return mod(diff);
}

inline void normalise(std::vector<double>& a)
{
  double length = mod(a);
  for(unsigned i=0; i<a.size(); i++)
    a[i] /= length;
}




//======================================================================
/// oomph-lib calculation of dgreendn:
//======================================================================
void analytic_integral_dgreendn_triangle(const std::vector<double>& x_sn,
					 const std::vector<std::vector<double> >& x_nds,
					 std::vector<double>& boundary_matrix)
{
  // Only works in 3D and for triangles (3 nodes)
  const unsigned node_dim = 3;
  const unsigned n_node = 3;

  /* First some pre-calculations to get everything ready. */

  // Calculate the length of the triangle sides and the unit vectors along the
  // triangle sides.
  std::vector<double> side_length(n_node,0.0);
  std::vector<double> side_length2(n_node,0.0);
  std::vector<std::vector<double> > side_direction(n_node);
  for(unsigned i=0; i < n_node; i++)
    {
      // Get the next node around the triangle (i.e. 0 -> 1 -> 2 -> 0 etc.).
      unsigned next_node = (i+1)%n_node;

      // Get the vector along this side (xi in Lindholm1984)
      vector_diff(x_nds[next_node],x_nds[i],side_direction[i]);

      // Get the length of this side (s in Lindholm1984)
      side_length[i] = mod(side_direction[i]);
    }

  // Calculate the non-unit normal to triangle. Assuming flat element => same
  // everywhere so we can just take the cross product of any two vectors in
  // the plane. Use the first two edge vectors.
  std::vector<double> unit_normal(node_dim,0.0);
  cross(side_direction[0],side_direction[1],unit_normal);

  // std::vector<double> unit_normal2(node_dim,0.0);
  // cross(side_direction[1],side_direction[2],unit_normal2);
  // normalise(unit_normal2);

  // std::vector<double> unit_normal3(node_dim,0.0);
  // cross(side_direction[2],side_direction[0],unit_normal3);
  // normalise(unit_normal3);

  // if ((mod_diff(unit_normal,unit_normal2) > 1e-8) ||
  // 	(mod_diff(unit_normal,unit_normal3) > 1e-8))
  //   {
  // 	std::cout << unit_normal << std::endl;
  // 	std::cout << unit_normal2 << std::endl;
  // 	std::cout << unit_normal3 << std::endl;
  // 	std::cout << std::endl;
  //   }

  // Calculate area of triangle using the cross product of the (unnormalised)
  // side_direction vectors, already calculated since it is the normal.
  double area = mod(unit_normal)/2;

  // //ok
  // std::cout << area<< std::endl;

  // Normalise the unit normal and side direction vectors now that we have the
  // area.
  normalise(unit_normal);
  for(unsigned i=0; i<n_node; i++) normalise(side_direction[i]);

  // Calculate gamma
  std::vector< std::vector<double> > gamma(n_node);
  for(unsigned i=0; i<3; i++)
    {
      gamma[i].assign(n_node,0.0); // Initialise gamma[i]
      unsigned next_node = (i+1)%n_node; // Get next triangle vertex
      for(unsigned j=0; j<n_node; j++)
	gamma[i][j] = dot(side_direction[next_node],side_direction[j]);
    }

  // // ok
  // std::cout << gamma[0] << gamma[1] << gamma[2] << std::endl;

  // /* Now evaluate the integral for every boundary node and add to
  //    boundary matrix. */
  // for(unsigned i_sn=0; i_sn < boundary_mesh_pt()->nnode(); i_sn++)
  //   {
      // // Get position of this node
      // std::vector<double> x_sn(node_dim,0.0);
      // boundary_mesh_pt()->node_pt(i_sn)->position(x_sn);

      // Calculate rho (vector from each node in the triangle to source node)
      // and it's length.
      std::vector<std::vector<double> > rho(n_node); std::vector<double> rhol(n_node);
      for(unsigned i_nd=0; i_nd<n_node; i_nd++)
	{
	  vector_diff(x_nds[i_nd],x_sn,rho[i_nd]);
	  rhol[i_nd] = mod(rho[i_nd]);
	}

      // Calculate zeta: the distance between the element and the source node
      // in the normal direction to the element.
      double zeta = dot(unit_normal,rho[0]);

      // If source node is in the plane of the element (i.e. zeta ~ 0) then
      // n.r = 0, nothing to calculate or add so we can move on to the next
      // source node.
      double tol = 1e-10;
      if( (std::abs(zeta) < tol) )
	{
	  // std::cout << zeta <<  " " << std::abs(zeta) <<  std::endl;
	  return;
	}

      // if( (rhol[0] < tol) || (rhol[1] < tol) || (rhol[2] < tol) )
      //   {
      //     std::cout << "rho0 = " << rho[0] << std::endl;
      //     std::cout << "rho1 = " << rho[1] << std::endl;
      //     std::cout << "rho2 = " << rho[2] << std::endl;
      //     std::cout << "unit_normal = " << unit_normal << std::endl;
      //     std::cout << "x_nds0 = " << x_nds[0] << std::endl;
      //     std::cout << "x_nds1 = " << x_nds[1] << std::endl;
      //     std::cout << "x_nds2 = " << x_nds[2] << std::endl;
      //     std::cout << "x_sn = " << x_sn << std::endl;
      //   }

      // Calculate "P" (see paper Lindholm1984) for each node in the triangle
      std::vector<double> P(n_node,0.0);
      for(unsigned i=0; i<n_node; i++)
	{
	  unsigned next_node = (i+1)%n_node;
	  P[i] = std::log( (rhol[i] + rhol[next_node] + side_length[i])
			   /(rhol[i] + rhol[next_node] - side_length[i]) );
	}

      // ok!      std::cout << P << std::endl;

      // Calculate the solid angle (see Lindholm 1984)
      double ratio = (rhol[0]*rhol[1]*rhol[2]
		      + rhol[0] * dot(rho[1],rho[2])
		      + rhol[1] * dot(rho[2],rho[0])
		      + rhol[2] * dot(rho[0],rho[1]))
	/
	std::sqrt(2.0
		  *(rhol[1]*rhol[2] + dot(rho[1],rho[2]) )
		  *(rhol[2]*rhol[0] + dot(rho[2],rho[0]) )
		  *(rhol[0]*rhol[1] + dot(rho[0],rho[1]) )
		  );
      // ok
      // std::cout << "ratio: " << ratio << std::endl;

      // double t_nom = rhol[0]*rhol[1]*rhol[2]
      // 	+ rhol[0] * dot(rho[1],rho[2])
      // 	+ rhol[1] * dot(rho[2],rho[0])
      // 	+ rhol[2] * dot(rho[0],rho[1]);

      // double t_denom = 	std::sqrt(2.0
      // 		  *(rhol[1]*rhol[2] + dot(rho[1],rho[2]) )
      // 		  *(rhol[2]*rhol[0] + dot(rho[2],rho[0]) )
      // 		  *(rhol[0]*rhol[1] + dot(rho[0],rho[1]) )
      // 		  );

      // std::cout << "t/t: " << t_nom/t_denom << std::endl;

      // t_nom=
      // 	rho1l*rho2l*rho3l+
      // 	rho1l*my_ddot(ND,rho2,1,rho3,1)+
      // 	rho2l*my_ddot(ND,rho3,1,rho1,1)+
      // 	rho3l*my_ddot(ND,rho1,1,rho2,1);
      // t_denom=
      // 	std::sqrt(2.0*
      // 		  (rho2l*rho3l+my_ddot(ND,rho2,1,rho3,1))*
      // 		  (rho3l*rho1l+my_ddot(ND,rho3,1,rho1,1))*
      // 		  (rho1l*rho2l+my_ddot(ND,rho1,1,rho2,1))
      // 		  );
      int sign = (zeta > 0.0 ? +1.0 : -1.0);

      //ok
      // std::cout << zeta << std::endl;

      // Round-off errors can cause ratio to be out of range for inverse cos
      // so we need to check it.
      if(ratio > 1.0) ratio = 1.0;
      else if(ratio < -1.0) ratio = -1.0;
      double omega = sign * 2 * std::acos(ratio);

      // // ok
      // std::cout << omega << std::endl;

      // Calculate eta: the unit vector normal to the side. Use it to
      // calculate etal: the distance in the diretion of eta to the source
      // node.
      std::vector<std::vector<double> > eta(3); std::vector<double> etal(3,0.0);
      for(unsigned i=0; i<n_node; i++)
	{
	  eta[i].assign(3,0.0);
	  cross(unit_normal,side_direction[i],eta[i]);
	  etal[i] = dot(eta[i],rho[i]);
	}

      // ok
      // std::cout << eta[0] << std::endl;
      // std::cout << eta[1] << std::endl;
      // std::cout << eta[2] << std::endl;


      /* Now put it all together and add the contribution to the boundary
	 element matrix */
      for(unsigned i_tn=0; i_tn < n_node; i_tn++)
	{
	  unsigned next_node = (i_tn+1)%n_node;

	  // Add contribution to the appropriate value in the boundary matrix
	  boundary_matrix[i_tn] = (side_length[next_node]/(8*PI*area))
	    *( //??temp: not sure wether to add this here or in driver
	      (etal[next_node] * omega)
	      //??temp change to + because I think definition of green's fn is opposite sign
	      - (zeta * dot(gamma[i_tn],P))
	       );


	  //std::cout << boundary_matrix(l[i_tn],i_sn) << std::endl;

	  //??ds do we HAVE to add solid angle contrib here? or can we add it in
	  //driver code like now?
	}
      //}
}

int main()
{

  // set up triangle corner loacations and source location
  std::vector<double> t1(ND,0.0);
  std::vector<double> t2(ND,0.0); t2[1] = 1.0;
  std::vector<double> t3(ND,0.0); t3[0] = 1.0;
  std::vector<double> source_node(ND,0.5);

  // Get magpar calculation result
  std::vector<double> magpar_result(NN,0.0);
  Bele(source_node,t1,t2,t3,magpar_result);

  std::cout << "magpar gives:" << magpar_result << std::endl;

  // Get my (oomph-lib) result
  std::vector<double> oomph_result(NN,0.0);
  std::vector<std::vector<double> > t123(NN);
  t123[0] = t1; t123[1] = t2; t123[2] = t3;
  analytic_integral_dgreendn_triangle(source_node,t123,oomph_result);

  std::cout << "oomph-lib gives:" << oomph_result << std::endl;

  return 0;
}
