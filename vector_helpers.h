#ifndef OOMPH_VECTOR_HELPERS_H
#define OOMPH_VECTOR_HELPERS_H

#include "../../src/generic/Vector.h"
#include "./prettyprint98.hpp"
#include <numeric>
#include <functional>

using namespace oomph;
using namespace StringConversion;

namespace VectorOps
{

  void check_lengths_match(const Vector<double> &a, const Vector<double> &b)
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

  // Probably not always best/fastest because not optimised for dimension but
  // useful...
  inline double dot(const Vector<double>& a, const Vector<double>& b)
  {
    check_lengths_match(a,b);
    double temp = 0;
    for(unsigned i=0; i<a.size(); i++)
      temp += a[i] * b[i];
    return temp;
  }

  /// Cross product using "proper" output (move semantics means this is ok
  /// nowadays).
  Vector<double> cross(const Vector<double>& A, const Vector<double>& B)
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

  /// Calculate the cross product of vectors A and B, store the result in
  /// vector output. NOTE: the cross product is only valid for 3-dimensional
  /// vectors
  void cross(const Vector<double>& A, const Vector<double>& B, Vector<double>& output)
  {output = cross(A, B);}

  inline double two_norm(const Vector<double>& a)
  {
    return std::sqrt(dot(a,a));
  }


  inline Vector<double> vector_diff(const Vector<double>& a, const Vector<double>& b)
  {
    check_lengths_match(a,b);
    Vector<double> diff(a.size(), 0.0);
    for(unsigned i=0; i<a.size(); i++) {diff[i] = a[i] - b[i];}
    return diff;
  }
  inline void vector_diff(const Vector<double>& a, const Vector<double>& b,
                          Vector<double>& diff)
  {diff = vector_diff(a, b);}


  /// \short Get the (smallest) angle between two vectors.
  double angle_diff(const Vector<double> &a, const Vector<double> &b)
    {
      // Use the dot product formula:
      return std::acos(dot(a, b) / (two_norm(a) * two_norm(b)));
    }


  inline Vector<double> abs_vector_diff(const Vector<double>& a, const Vector<double>& b)
  {
    check_lengths_match(a,b);
    Vector<double> diff(a.size(), 0.0);
    for(unsigned i=0; i<a.size(); i++) {diff[i] = std::abs(a[i] - b[i]);}
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
    Vector<double> diff(a.size(), 0.0);
    for(unsigned i=0; i<a.size(); i++)
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
    Vector<double> diff(a.size(),0.0);
    vector_diff(a,b,diff);
    return two_norm(diff);
  }

  inline void normalise(Vector<double>& a)
  {
    double length = two_norm(a);
    for(unsigned i=0; i<a.size(); i++)
      {
        a[i] /= length;
      }
  }


  // Equivalent to std::find but for floating point values. Return -1 if
  // not found.
  int fp_find(double search_value, const Vector<double> &vec, double tol=1e-12)
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



  void rowstart2rowindex(const Vector<int>& row_start, Vector<int>& row_index)
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

  // Warning: must be sorted
  void rowindex2rowstart(const Vector<int>& row_index, Vector<int>& row_start)
  {
    row_start.clear();
    row_start.reserve(row_index.back()+1);

    row_start.push_back(0);

    int i=0;
    for(int row = 0; row < row_index.back() + 1; row++)
      {
        int count = 0;

        // If our row is smaller than row_index[i] then it has no entries
        if (row < row_index[i])
          {
            row_start.push_back(row_start.back() + count);
          }
        else
          {
            // Otherwise: count the number of entries
            while((unsigned(i) < row_index.size()) && (row == row_index[i]))
              {
                count++;
                i++;
              }
            row_start.push_back(row_start.back() + count);
          }
      }


  }

  void diag_cr_matrix(CRDoubleMatrix& cr_matrix, const unsigned& n,
                      const double& value)
  {
    Vector<int> row_index(n), col_index(n), row_start;
    Vector<double> values(n,value);

    for(unsigned i=0; i<n; i++)
      {
        row_index[i] = i;
        col_index[i] = i;
      }

    rowindex2rowstart(row_index,row_start);
    cr_matrix.build(n, values, col_index, row_start);
  }

  void random_single_element_per_row_cr_matrix(CRDoubleMatrix& cr_matrix,
                                               const unsigned& n,
                                               const unsigned& nnz,
                                               const double& max_value=100)
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

    rowindex2rowstart(row_index, row_start);
    // std::cout << row_start << std::endl;

    cr_matrix.build(n, values, col_index, row_start);
  }

  void get_as_indicies(DoubleMatrixBase &matrix,
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

  double rel_dense_matrix_diff(DenseMatrix<double> &mat1, DenseMatrix<double> &mat2)
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


  bool numerically_close(const Vector<double> &x1,
                         const Vector<double> &x2,
                         const double& tol=1e-10)
  {
    return numerical_zero(two_norm_diff(x1,x2));
  }

  template <typename T>
  T mean(const Vector<T> &vec)
  {
    return std::accumulate(vec.begin(), vec.end(), 0.0) / double(vec.size());
  }

  template <typename T>
  T max(const Vector<T> &vec)
  {
    return *std::max_element(vec.begin(), vec.end());
  }

  template <typename T>
  double stddev(const Vector<T> &vec)
  {
    double vec_mean = mean(vec);
    double sum_square_deviations = 0.0;
    for(unsigned i=0; i<vec.size(); i++)
      {
        sum_square_deviations +=
          std::pow(vec[i] - vec_mean,2);
      }

    return std::sqrt(sum_square_deviations / double(vec.size()));
  }


  /// Check if a vector contains any duplicate values
  template <typename T> bool contains_duplicates(const std::vector<T> &v)
  {
    // Construct a set (which has no duplicates by definition) and compare
    // sizes.
    return std::set<T>(v.begin(), v.end()).size() != v.size();
  }

}

#endif
