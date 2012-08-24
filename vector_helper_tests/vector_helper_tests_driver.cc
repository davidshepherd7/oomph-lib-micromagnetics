#include "../vector_helpers.h"
#include "../prettyprint98.hpp"

using namespace oomph;
using namespace VectorOps;

bool equal(const Vector<int>& a, const Vector<int>& b)
{
  if(a.size() != b.size())
    return false;

  for(unsigned i=0; i< a.size(); i++)
    {
      if(a[i] != b[i])
	return false;
    }

  return true;
}

int main()
{
  Vector<int> row_index, row_start, new_row_index;
  row_index.push_back(0);
  row_index.push_back(1);
  row_index.push_back(6);
  row_index.push_back(6);
  row_index.push_back(6);
  row_index.push_back(7);
  row_index.push_back(7);


  // Convert to CR form and back again
  rowindex2rowstart(row_index,row_start);
  rowstart2rowindex(row_start,new_row_index);

  if(!(equal(row_index,new_row_index)))
    {
      std::cout << "Failed rowindex2rowstart or inverse"  << std::endl;
      std::cout << row_index << std::endl;
      std::cout << row_start << std::endl;
      std::cout << new_row_index << std::endl;
      return 1;
    }

  // Try making a matrix from it
  Vector<double> values(row_index.size(), 1.0);
  Vector<int> cols;
  cols.push_back(0);
  cols.push_back(1);
  cols.push_back(4);
  cols.push_back(5);
  cols.push_back(6);
  cols.push_back(1);
  cols.push_back(7);

  LinearAlgebraDistribution linalg(0,8,false);
  CRDoubleMatrix test(&linalg);
  test.build(8,values,cols,row_start);

  test.sparse_indexed_output(std::cout);

  return 0;
}
