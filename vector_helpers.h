#ifndef OOMPH_VECTOR_HELPERS_H
#define OOMPH_VECTOR_HELPERS_H

#include "./prettyprint98.hpp"

#include <cmath>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/matrices.h"

// Note: all functions must be declared inline so that the compiler doesn't
// die! Stupid compiler...

namespace VectorOps
{
  using namespace oomph;
  using namespace StringConversion;

  inline void check_lengths_match(const Vector<double> &a, const Vector<double> &b)
  {
#ifdef PARANOID
    if (a.size() != b.size())
      {
        std::string err = "Vectors must be the same length. ";
        err += "len(a) = " + to_string(a.size()) + ", len(b) = " + to_string(b.size());
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
  }

  inline void check_lengths_match(const Vector<double>& a, const unsigned& len)
  {
#ifdef PARANOID
    if (a.size() != len)
      {
        std::string err = "Vectors must be the same length. ";
        err += "len(a) = " + to_string(a.size()) + ", array len = "
          + to_string(len);
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
  }

  // Probably not always best/fastest because not optimised for dimension but
  // useful...
  inline double dot(const Vector<double>& a, const Vector<double>& b)
  {
    check_lengths_match(a,b);
    double temp = 0;
    for(unsigned i=0, ni=a.size(); i<ni; i++)
      {
        temp += a[i] * b[i];
      }
    return temp;
  }
  inline double dot(const double* const a, const double* const b,
                    const unsigned& len)
  {
    double temp = 0;
    for(unsigned i=0, ni=len; i<ni; i++)
      {
        temp += a[i] * b[i];
      }
    return temp;
  }
  inline double dot(const Vector<double>& a, const double* const b,
                    const unsigned& len)
  {
    check_lengths_match(a, len);
    double temp = 0;
    for(unsigned i=0, ni=len; i<ni; i++)
      {
        temp += a[i] * b[i];
      }
    return temp;
  }
  inline double dot(const double* const a, const Vector<double>& b,
                    const unsigned& len)
  {
    check_lengths_match(b, len);
    double temp = 0;
    for(unsigned i=0, ni=len; i<ni; i++)
      {
        temp += a[i] * b[i];
      }
    return temp;
  }

  /// Cross product using "proper" output (move semantics means this is ok
  /// nowadays).
  inline Vector<double> cross(const Vector<double>& A, const Vector<double>& B)
  {
#ifdef PARANOID
    if((A.size() != 3) || (B.size() != 3))
      {
        std::string err = "Cross product only defined for vectors of length 3.";
        err += "len(a) = " + to_string(A.size()) + ", len(b) = " + to_string(B.size());
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    Vector<double> output(3,0.0);
    output[0] = A[1]*B[2] - A[2]*B[1];
    output[1] = A[2]*B[0] - A[0]*B[2];
    output[2] = A[0]*B[1] - A[1]*B[0];
    return output;
  }

  // Collection of cross product functions: get a single row of the cross
  // product. Can use double* or vectors for compatability.
  inline double opt_cross(unsigned i, const double* A, const double* B)
  {
    unsigned ia = ((i + 1) % 3), ib = ((i + 2) % 3);
    return A[ia]*B[ib] - A[ib]*B[ia];
  }
  inline double opt_cross(unsigned i, const Vector<double>& A, const double* B)
  {
    unsigned ia = ((i + 1) % 3), ib = ((i + 2) % 3);
    return A[ia]*B[ib] - A[ib]*B[ia];
  }
  inline double opt_cross(unsigned i, const double* A, const Vector<double>& B)
  {
    unsigned ia = ((i + 1) % 3), ib = ((i + 2) % 3);
    return A[ia]*B[ib] - A[ib]*B[ia];
  }
  inline double opt_cross(unsigned i, const Vector<double>& A,
                          const Vector<double>& B)
  {
    unsigned ia = ((i + 1) % 3), ib = ((i + 2) % 3);
    return A[ia]*B[ib] - A[ib]*B[ia];
  }

  inline double opt_double_cross(unsigned i,
                                 const Vector<double>& A,
                                 const Vector<double>& B,
                                 const Vector<double>& C)
  {
    if(i==0)
      {
        return A[1]*(B[0]*C[1] - B[1]*C[0]) - A[2]*(B[2]*C[0] - C[2]*B[0]);
      }
    else if(i==1)
      {
        return A[2]*(B[1]*C[2] - B[2]*C[1]) - A[0]*(B[0]*C[1] - C[0]*B[1]);
      }
    else if(i==2)
      {
        return A[0]*(B[2]*C[0] - B[0]*C[2]) - A[1]*(B[1]*C[2] - C[1]*B[2]);
      }
    else
      {
        throw OomphLibError("i index out of range", OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
  }

  inline double opt_double_cross(unsigned i,
                                 const double* A,
                                 const double* B,
                                 const Vector<double>& C)
  {
    if(i==0)
      {
        return A[1]*(B[0]*C[1] - B[1]*C[0]) - A[2]*(B[2]*C[0] - C[2]*B[0]);
      }
    else if(i==1)
      {
        return A[2]*(B[1]*C[2] - B[2]*C[1]) - A[0]*(B[0]*C[1] - C[0]*B[1]);
      }
    else if(i==2)
      {
        return A[0]*(B[2]*C[0] - B[0]*C[2]) - A[1]*(B[1]*C[2] - C[1]*B[2]);
      }
    else
      {
        throw OomphLibError("i index out of range", OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
  }





  /// Calculate the cross product of vectors A and B, store the result in
  /// vector output. NOTE: the cross product is only valid for 3-dimensional
  /// vectors
  inline void cross(const Vector<double>& A, const Vector<double>& B, Vector<double>& output)
  {output = cross(A, B);}

  inline double two_norm(const Vector<double>& a)
  {
    return std::sqrt(dot(a,a));
  }


  inline Vector<double> vector_diff(const Vector<double>& a, const Vector<double>& b)
  {
    check_lengths_match(a,b);
    unsigned ni = a.size();
    Vector<double> diff(ni, 0.0);
    for(unsigned i=0; i<ni; i++) {diff[i] = a[i] - b[i];}
    return diff;
  }
  inline void vector_diff(const Vector<double>& a, const Vector<double>& b,
                          Vector<double>& diff)
  {diff = vector_diff(a, b);}


  inline double max_vector_diff(const Vector<double> &a,
                                const Vector<double> &b)
    {
      VectorOps::check_lengths_match(a, b);
      double diff = 0;
      for(unsigned i=0, ni=a.size(); i<ni; i++)
        {
          double this_diff = std::abs(a[i] - b[i]);
          if(this_diff > diff) diff = this_diff;
        }
      return diff;
    }


  /// \short Get the (smallest) angle between two vectors.
  inline double angle_diff(const Vector<double> &a, const Vector<double> &b)
    {
      // Use the dot product formula:
      double temp = dot(a, b) / (two_norm(a) * two_norm(b));

      // Be safe for slightly wrong floating point values
      if(temp > 1.0)
        {
          if(temp < 1.0 + 1e-12)
            {
              temp = 1.0;
            }
          else
            {
              throw OomphLibError(to_string(temp) +" is out of range",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
        }

      return std::acos(temp);
    }


  inline Vector<double> abs_vector_diff(const Vector<double>& a, const Vector<double>& b)
  {
    check_lengths_match(a,b);
    unsigned ni = a.size();
    Vector<double> diff(ni, 0.0);
    for(unsigned i=0; i<ni; i++) {diff[i] = std::abs(a[i] - b[i]);}
    return diff;
  }
  inline void abs_vector_diff(const Vector<double>& a, const Vector<double>& b,
                                 Vector<double>& diff)
  {diff = abs_vector_diff(a, b);}


  inline bool numerical_zero(const double &a, const double& tol=1e-10)
  {
    return std::abs(a) < tol;
  }

  inline Vector<double> relative_abs_vector_diff(const Vector<double>& a,
                                                 const Vector<double>& b)
  {
    check_lengths_match(a,b);
    unsigned ni = a.size();
    Vector<double> diff(ni, 0.0);
    for(unsigned i=0; i<ni; i++)
      {
        // if a[i] is not zero then just do it normally
        if( !(numerical_zero(a[i])))
          {
            diff[i] = std::abs( (a[i] - b[i]) / a[i] );
          }

        // If a is zero but b isn't then relative error is large
        else if( !(numerical_zero(b[i]))) diff[i] = 1.0;

        // If both values are roughly zero then there is no error
        else diff[i] = 0.0;
      }

    return diff;
  }

  inline void relative_abs_vector_diff(const Vector<double>& a,
                                       const Vector<double>& b,
                                       Vector<double>& diff)
    {diff = relative_abs_vector_diff(a,b);}

  inline double two_norm_diff(const Vector<double>& a, const Vector<double>& b)
  {
    Vector<double> diff;
    vector_diff(a,b,diff);
    return two_norm(diff);
  }

  inline void normalise(Vector<double>& a)
  {
    double length = two_norm(a);
    for(unsigned i=0, ni=a.size(); i<ni; i++)
      {
        a[i] /= length;
      }
  }


  // Equivalent to std::find but for floating point values. Return -1 if
  // not found.
  inline int fp_find(double search_value, const Vector<double> &vec,
                     double tol=1e-12)
    {

      int found_location = -1;
      for(unsigned j=0, nj=vec.size(); j<nj; j++)
        {
          if(std::abs(vec[j] - search_value) < tol)
            {
              found_location = j;
              break;
            }
        }

      return found_location;
    }



  inline void rowstart2rowindex(const Vector<int>& row_start,
                                Vector<int>& row_index)
  {
    // Initialise
    int nrow = row_start.back();
    row_index.reserve(nrow);

    int row = 0, i_row_index = 0;
    for(int i_row_start=0; i_row_start < int(row_start.size()); i_row_start++)
      {
        int next_row_start = row_start[i_row_start + 1];
        while(i_row_index < next_row_start)
          {
            row_index.push_back(row);
            i_row_index++;
          }
        row++;
      }
  }

  /// Construct a row start vector from a row index vector. The row index
  /// vector must be sorted.
  inline void rowindex2rowstart(const Vector<int>& row_index,
                                const unsigned& nrow,
                                Vector<int>& row_start)
  {
#ifdef PARANOID
    // Check that the row index vector is sorted. Use this weird thing
    // instead of std::is_sorted because we don't have c++11. Looks for
    // adjacent elements s.t. e1 > e2, i.e. not sorted.
    if(std::adjacent_find(row_index.begin(), row_index.end(),
                          std::greater<double>()) != row_index.end())
      {
        std::string err = "Row index vector must be sorted!";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    row_start.clear();
    row_start.reserve(row_index.back()+1);

    row_start.push_back(0);

    unsigned i=0;
    for(int row = 0; row < int(nrow); row++)
      {
        int count = 0;

        // Otherwise: count the number of entries, the next row starts
        // that number of entries further along.
        while((i < row_index.size()) && (row == row_index[i]))
          {
            count++;
            i++;
          }
        row_start.push_back(row_start.back() + count);
      }

    // Final entry of row_start vector:
    row_start.push_back(row_index.size());
  }


  inline void diag_cr_matrix(CRDoubleMatrix& cr_matrix, const unsigned& n,
                             const double& value)
  {
    Vector<int> row_index(n), col_index(n), row_start;
    Vector<double> values(n,value);

    for(unsigned i=0; i<n; i++)
      {
        row_index[i] = i;
        col_index[i] = i;
      }

    rowindex2rowstart(row_index, n, row_start);
    cr_matrix.build(n, values, col_index, row_start);
  }

  inline Vector<double> random_vector(unsigned length, double max_value=100)
    {
      Vector<double> a(length);

      for(unsigned i=0; i<length; i++)
        {
          a[i] = double(rand() % int(10000*max_value)) / 10000.0;
        }
      return a;
    }


  inline void random_single_element_per_row_cr_matrix
  (CRDoubleMatrix& cr_matrix, const unsigned& n,
   const unsigned& nnz, const double& max_value=100)
  {
    Vector<int> row_index(nnz), col_index(nnz), row_start;
    Vector<double> values(nnz);

    for(unsigned j=0; j<nnz; j++)
      {
        row_index[j] = j;
        col_index[j] = (rand() % n);
        values[j] = double(rand() % int(1000*max_value)) / 1000.0;
      }


    std::sort(row_index.begin(), row_index.end());

    // std::cout << row_index << std::endl;
    // std::cout << col_index << std::endl;
    // std::cout << values << std::endl;

    rowindex2rowstart(row_index, n, row_start);
    // std::cout << row_start << std::endl;

    cr_matrix.build(n, values, col_index, row_start);
  }

  inline void get_as_indicies(const DoubleMatrixBase &matrix,
                              Vector<double> &values,
                              Vector<int> &col_index,
                              Vector<int> &row_index)
  {
    for(unsigned i=0; i< matrix.nrow(); i++)
      {
        for(unsigned j=0; j< matrix.ncol(); j++)
          {
            if(matrix(i,j) != 0)
              {
                row_index.push_back(i);
                col_index.push_back(j);
                values.push_back(matrix(i,j));
              }
          }
      }
  }

  /// Get three vectors of the values, column indicies and row indidices of
  /// entries in a CR matrix. Sorted by row index then col index. Optimised
  /// version for CR matrices
  inline void get_as_indicies(const CRDoubleMatrix &matrix,
                              Vector<double> &values,
                              Vector<int> &col_index,
                              Vector<int> &row_index)
  {
    // Reserve space so we don't re-allocate the vectors inside the loop.
    int nval = matrix.nnz();
    values.reserve(nval);
    col_index.reserve(nval);
    row_index.reserve(nval);

    // Loop over rows
    const int* row_start = matrix.row_start();
    for(int i = 0; i < int(matrix.nrow_local()); i++)
      {
        // Throw col/val pairs in this row into a map (to get them sorted)
        std::map<int, double> col_val_map;
        for(int entry=row_start[i]; entry<row_start[i+1]; entry++)
          {
            col_val_map[matrix.column_index()[entry]] = matrix.value()[entry];
          }

        // Pull them out of the map and into the output vectors
        std::map<int, double>::const_iterator it;
        for(it=col_val_map.begin(); it!=col_val_map.end(); ++it)
          {
            row_index.push_back(i);
            col_index.push_back(it->first);
            values.push_back(it->second);
          }
      }
  }

  inline double rel_dense_matrix_diff(DenseMatrix<double> &mat1,
                                      DenseMatrix<double> &mat2)
  {
    if((mat1.nrow() != mat2.nrow())
       || (mat1.ncol() != mat2.ncol()))
      {
        std::cout <<  "Different number of rows/cols" << std::endl;
        return 2.341e200;
      }

    double total_diff = 0.0;
    for(unsigned i=0; i< mat1.nrow(); i++)
      {
        for(unsigned j=0; j< mat1.ncol(); j++)
          {
            double val = std::abs(mat1(i,j));
            if(!numerical_zero(val))
              {
                total_diff += mat1(i,j) - mat2(i,j) / val;
              }
          }
      }

    return total_diff / double( mat1.nrow() * mat1.ncol());
  }


  inline bool numerically_close(const Vector<double> &x1,
                                const Vector<double> &x2,
                                const double& tol=1e-10)
  {
    return numerical_zero(two_norm_diff(x1,x2));
  }

  template <typename T>
  inline T mean(const Vector<T> &vec)
  {
    return std::accumulate(vec.begin(), vec.end(), 0.0) / double(vec.size());
  }

  template <typename T>
  inline T max(const Vector<T> &vec)
  {
    return *std::max_element(vec.begin(), vec.end());
  }

  template <typename T>
  inline double stddev(const Vector<T> &vec)
  {
    double vec_mean = mean(vec);
    double sum_square_deviations = 0.0;
    unsigned vec_size = vec.size();
    for(unsigned i=0; i<vec_size; i++)
      {
        sum_square_deviations +=
          std::pow(vec[i] - vec_mean,2);
      }

    return std::sqrt(sum_square_deviations / double(vec_size));
  }


  /// Check if a vector contains any duplicate values
  template <typename T>
  inline bool contains_duplicates(const std::vector<T> &v)
  {
    // Construct a set (which has no duplicates by definition) and compare
    // sizes.
    return std::set<T>(v.begin(), v.end()).size() != v.size();
  }


  inline void check_matrices_compatible(const CRDoubleMatrix& m1,
                                        const CRDoubleMatrix& m2)
  {
#ifdef PARANOID
    if(m1.nrow() != m2.nrow())
      {
        std::string err = "The two matrices have a different number of rows";
        err += "matrix 1 has " + to_string(m1.nrow()) + " rows ";
        err += "but matrix 2 has " + to_string(m2.nrow()) +".";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
    if(m1.ncol() != m2.ncol())
      {
        std::string err = "The two matrices have a different number of cols";
        err += "matrix 1 has " + to_string(m1.ncol()) + " cols ";
        err += "but matrix 2 has " + to_string(m2.ncol()) +".";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }

    if(*(m1.distribution_pt()) != *(m2.distribution_pt()))
      {
        std::string err = "The two matrices have different distributions.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
  }

  //??ds move to class?
  inline void check_matrix_built(const CRDoubleMatrix& m)
  {
#ifdef PARANOID
    if(!m.built())
      {
        std::string err = "This matrix is not built.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
  }

  /// \short Add two matrices together and put the result into mout. The
  /// two input matrices must have the same shape and distribution. The
  /// output matrix will have the same shape and distribution as the
  /// inputs. It is safe to use one of the inputs as the output.
  inline void cr_matrix_add(const CRDoubleMatrix& m1, const CRDoubleMatrix& m2,
                            CRDoubleMatrix& mout)
  {
    // Basic idea is:
    // 1. Get both matrices in coordinate form to make them easier to work
    // with.
    // 2. Merge the entries into one using something similar to the merge
    // from mergesort. (Two counters, one for each matrix. Compare row/col
    // of current entry in each matrx. Insert the smallest and increment
    // that counter. Loop)
    // 3. Convert back to compressed row format and build the new matrix.

    // Possible optimisations if needed:
    // 1. Use row start format directly
    // 2. Use C-arrays instead of vectors and use build without copy function.


    // Paranoid tests
    check_matrix_built(m1);
    check_matrix_built(m2);
    check_matrices_compatible(m1, m2);

    // Get the values of the matrices in sorted coordinate form
    Vector<double> v1, v2;
    Vector<int> c1, c2, r1, r2;
    get_as_indicies(m1, v1, c1, r1);
    get_as_indicies(m2, v2, c2, r2);

    // Get the sizes of the input matrix coordinate vectors
    int n1 = c1.size(), n2 = c2.size();

    // Create output matrix vectors and reserve space. Get enough space for
    // the case where m1 and m2 have no overlapping values because it's an
    // easy to calculate upper bound for the amount we actually need.
    Vector<double> vs;
    Vector<int> cs, rs;
    vs.reserve(n1 + n2);
    cs.reserve(n1 + n2);
    rs.reserve(n1 + n2);

    // Loop over all entries in the two lists of coordinates while at least
    // one of them hasn't reached the end
    int i1 = 0, i2 = 0;
    while(i1 < n1 && i2 < n2)
      {
        // If the current entry in matrix 1 is first (i.e. on an earlier
        // row or column) then insert it and increment its counter.
        if(r1[i1] < r2[i2] || (r1[i1] == r2[i2] && c1[i1] < c2[i2]))
          {
            vs.push_back(v1[i1]);
            cs.push_back(c1[i1]);
            rs.push_back(r1[i1]);
            i1++;
          }
        // If the current entry in matrix 2 is first then insert it
        // and increment its counter.
        else if(r1[i1] > r2[i2] || (r1[i1] == r2[i2] && c1[i1] > c2[i2]))
          {
            vs.push_back(v2[i2]);
            cs.push_back(c2[i2]);
            rs.push_back(r2[i2]);
            i2++;
          }
        // If they are both the same row and col then add them together and
        // insert, increment both counters.
        else if(r1[i1] == r2[i2] && c1[i1] == c2[i2])
          {
            vs.push_back(v1[i1] + v2[i2]);
            cs.push_back(c1[i1]);
            rs.push_back(r1[i1]);
            i1++;
            i2++;
          }
        // Never get here I hope
        else
          {
            std::string err = "Never get here!";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
      }

    // Copy over the rest
    while(i1 < n1)
      {
        vs.push_back(v1[i1]);
        cs.push_back(c1[i1]);
        rs.push_back(r1[i1]);
        i1++;
      }
    while(i2 < n2)
      {
        vs.push_back(v2[i2]);
        cs.push_back(c2[i2]);
        rs.push_back(r2[i2]);
        i2++;
      }


    // Convert the row index vector to row starts
    Vector<int> sum_row_start;
    rowindex2rowstart(rs, m1.nrow(), sum_row_start);

    // Copy out dist and ncol so that we can't accidentally delete them
    // before using them if the output matrix is one of the inputs. Safe
    // to use m1 for distribution, ncol because we checked that they are
    // the same before.
    LinearAlgebraDistribution dist = *m1.distribution_pt();
    int ncol = m1.ncol();

    // Finally, build the output matrix. Because we do this only at the very
    // end, after all processing has been done, it is safe to output the
    // result into one of the input matrices if needed.
    mout.build(&dist, ncol, vs, cs, sum_row_start);

  }


}

#endif
