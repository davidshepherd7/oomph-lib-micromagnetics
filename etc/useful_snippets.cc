
/// Output stuff from main()
// ============================================================

// dump mesh for testing
unsigned n_nd = problem.mesh_pt()->nnode();
std::ofstream mesh_plot;
mesh_plot.open("./mesh_points");
for(unsigned nd=0; nd<n_nd; nd++)
  {
    for(unsigned j=0; j<dim; j++)
      mesh_plot << problem.mesh_pt()->node_pt(nd)->x(j) << " ";
    mesh_plot << std::endl;
  }
mesh_plot.close();

// dump boundary for testing
unsigned b = 0;
std::ofstream bound_plot;
bound_plot.open("./bound_points");
for(unsigned nd=0; nd<problem.mesh_pt()->nboundary_node(b); nd++)
  {
    for(unsigned j=0; j<dim; j++)
      bound_plot << problem.mesh_pt()->boundary_node_pt(b,nd)->x(j) << " ";
    bound_plot << std::endl;
  }
bound_plot.close();

// dump initial Jacobian and residuals for tests
DenseDoubleMatrix jacobian;
DoubleVector residuals;
problem.get_jacobian(residuals,jacobian);
std::ofstream matrix_file;
matrix_file.precision(16);
char filename[100];
sprintf(filename,"results/jacobian");
matrix_file.open(filename);
jacobian.output(matrix_file);
matrix_file.close();
std::ofstream residual_file;
residual_file.precision(16);
char filename2[100];
sprintf(filename2,"results/residual");
residual_file.open(filename2);
residuals.output(residual_file);
residual_file.close();

// dump boundary matrix
std::ofstream bem_file;
bem_file.precision(16);
char bem_filename[100];
sprintf(bem_filename,"results/bem");
bem_file.open(bem_filename);
problem.boundary_matrix_pt()->output(bem_file);
bem_file.close();
