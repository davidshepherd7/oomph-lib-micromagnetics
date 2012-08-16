#ifndef OOMPH_VECTOR_HELPERS_H
#define OOMPH_VECTOR_HELPERS_H

#include "generic.h"

using namespace oomph;

namespace VectorOps
{

  // Probably not always best/fastest because not optimised for dimension but
  // useful...
  inline double dot(const Vector<double>& a, const Vector<double>& b)
  {
#ifdef PARANOID
    if (a.size() != b.size())
      throw OomphLibError("Vectors must be the same length", "mod_diff",
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
    void cross(Vector<double>& A, Vector<double>& B, Vector<double>& output)
    {
#ifdef PARANOID
      if((A.size() != 3) || (B.size() != 3))
	{
	  std::ostringstream error_msg;
	  error_msg << "Cross product only defined for vectors of length 3.";
	  throw OomphLibError(error_msg.str(),
			      "VectorOps::cross",
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif
      output.assign(3,0.0);
      output[0] = A[1]*B[2] - A[2]*B[1];
      output[1] = A[2]*B[0] - A[0]*B[2];
      output[2] = A[0]*B[1] - A[1]*B[0];
    }


  inline double mod(const Vector<double>& a)
  {
    return std::sqrt(dot(a,a));
  }

  inline void vector_diff(const Vector<double>& a, const Vector<double>& b,
			      Vector<double>& diff)
  {
#ifdef PARANOID
    if (a.size() != b.size())
      throw OomphLibError("Vectors must be the same length", "mod_diff",
    			  OOMPH_EXCEPTION_LOCATION);
#endif
    diff.assign(a.size(),0.0);
    for(unsigned i=0; i<a.size(); i++)
      diff[i] = a[i] - b[i];
  }

  inline void abs_vector_diff(const Vector<double>& a, const Vector<double>& b,
			      Vector<double>& diff)
  {
    vector_diff(a,b,diff);
    for(unsigned i=0; i<a.size(); i++)
      diff[i] = std::abs(diff[i]);
  }

  inline double mod_diff(const Vector<double>& a, const Vector<double>& b)
  {
    Vector<double> diff(a.size(),0.0);
    vector_diff(a,b,diff);
    return mod(diff);
  }

  inline void normalise(Vector<double>& a)
  {
    double length = mod(a);
    for(unsigned i=0; i<a.size(); i++)
      a[i] /= length;
  }

}

#endif
