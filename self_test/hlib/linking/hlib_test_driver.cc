
#include "generic.h"


#include "bem3d.h"
#include "cluster.h"
#include "h2conversion.h"
#include "h2virtual.h"
#include "h2apriori.h"
#include "h2arithmetics.h"
#include "hca.h"
#include "laplacebem.h"
#include "surfacebem.h"
#include "supermatrix.h"
#include "factorisations.h"
#include "hcoarsening.h"
#include "graphics.h"
#include "krylov.h"
#include "norm.h"

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Solution: u(x,y,z) = x^2-z^2 */

static double
dirichlet_quadratic(const double *x, const double *n)
{
  (void) n;
  return x[0] * x[0] - x[2] * x[2];
}


int main()
{

  std::cout << "working" << std::endl;

  pmacrogrid3d mg;
  pbemgrid3d gr;
  pclustertree ct;
  pcluster root;
  pblockcluster blk;
  psurfacebemfactory bfactory;
  psupermatrix Kmod = NULL;
  psupermatrix M = NULL;
  psupermatrix M_compress = NULL;
  psupermatrix M_full = NULL;
  psupermatrix M_precond = NULL;
  psparsematrix Mass;
  peval_powerdata epd;
  ppowerdata pd;
  double *x0, *x1, *x2;
  char buf[80], filename[80], operator_char, geometry;
  double norm, normfull, eta;
  int acak, output, linear;
  int compression_method, order_interpolation;
  int reference, lower;
  double h_coarsen;
  double aca_eps;
  struct tms buffer;
  clock_t t1, t2;
  time_t t;
  long clk_tck;
  int split, nfdeg, qmin=2, rows, cols;
  int i, nmin, maxiter,j;
  size_t size_m;
  int clusters;
  int testdis=0, dist_adaptive;
  ph2recompression h2r;
  EvalMode mode;
  BEMBasis basisfunc=HLIB_CONSTANT_BASIS;


  clk_tck = sysconf(_SC_CLK_TCK);

  (void) time(&t);
  srand(t);

  split = 16;
  filename[0] = '\0';
  operator_char = 'd';
  output = 0;

  aca_eps = 1e-4;
  acak = 200;
  eta = 2.0;
  nfdeg = 3;
  nmin = 16;
  compression_method = 4;
  order_interpolation = 4;
  linear = 1;
  testdis = 1;
  h_coarsen = -1;
  dist_adaptive = 1;
  qmin = 0;
  geometry = 's';

  basisfunc=(linear ? HLIB_LINEAR_BASIS : HLIB_CONSTANT_BASIS);

  (void) printf("==============================================\n");
  (void) printf("#  Various BEM 3D approximations             #\n");
  (void) printf("==============================================\n");

  (void) printf("# Split edges into how many parts? (or filename for file input) [%d]\n",
		split);
  // (void) fgets(buf, 80, stdin);
  // remove_newline(buf);
  sprintf(buf, "8");

  if(sscanf(buf, "%d", &split) == 1) {
    assert(split >= 0);
    filename[0] = '\0';
    geometry = 's';
  }
  else if(tolower(buf[0]) == 's') {
    sscanf(buf+1, "%d", &split);
    assert(split >= 0);
    filename[0] = '\0';
    geometry = 's';
  }
  else if(tolower(buf[0]) == 'b') {
    sscanf(buf+1, "%d", &split);
    assert(split >= 0);
    filename[0] = '\0';
    geometry = 'b';
  }
  else if(tolower(buf[0]) == 'k') {
    sscanf(buf+1, "%d", &split);
    assert(split >= 0);
    (void) strcpy(filename,"kurbel25744.tri");
    geometry = 'f';
  }
  else if(buf[0] != '\0') {
    split = 0;
    (void) strcpy(filename, buf);
    geometry = 'f';
  }

  (void) printf("# D)ouble or S)ingle layer potential? [%c]\n",
		operator_char);
  // (void) fgets(buf, 80, stdin);
  // remove_newline(buf);
  sscanf(buf, "%c", &operator_char);
  operator_char = tolower(operator_char);

  (void) printf("# Initial compression by\n");
  (void) printf("# 0 ACA\n# 1 ACA+\n# 2 Interpolation\n");
  (void) printf("# 3 HCA(I)\n# 4 HCA(II)\n# 5 H2\n# 6 HCA-H2 ? [%d]\n",
		compression_method);
  // (void) fgets(buf, 80, stdin);
  sscanf(buf, "%d", &compression_method);

  switch(compression_method){
  case 0:
  case 1:
    (void) printf("# Accuracy for ACA? [%g]\n",
		  aca_eps);
    (void) fgets(buf, 80, stdin);
    sscanf(buf, "%lf", &aca_eps);
    break;
  case 2:
  case 5:
  case 6:
    (void) printf("# Order of the interpolation? [%d]\n",
		  order_interpolation);
    (void) fgets(buf, 80, stdin);
    sscanf(buf, "%d", &order_interpolation);
    (void) printf("# Accuracy for recompression? [%g]\n",
		  aca_eps);
    (void) fgets(buf, 80, stdin);
    sscanf(buf, "%lf", &aca_eps);
    break;
  default:
  case 3:
  case 4:
    (void) printf("# Accuracy for HCA? [%g]\n",
		  aca_eps);
    // (void) fgets(buf, 80, stdin);
    sscanf(buf, "%lf", &aca_eps);
    (void) printf("# Order of the interpolation for HCA? [%d]\n",
		  order_interpolation);
    // (void) fgets(buf, 80, stdin);
    sscanf(buf, "%d", &order_interpolation);
    break;
  }

  (void) printf("# Nearfield quadrature order? [%d]\n",
		nfdeg);
  // (void) fgets(buf, 80, stdin);
  sprintf(buf, "3");
  sscanf(buf, "%d", &nfdeg);

  (void) printf("# Output matrices: 0=No or 1=Rank or 2=SVD? [%d]\n",
		output);
  // (void) fgets(buf, 80, stdin);
  sscanf(buf, "%d", &output);

  switch(geometry) {
  default:
    (void) fprintf(stderr,
		   "Illegal geometry %c (%d), assuming sphere\n",
		   geometry, (int) geometry);
  case 's':
    (void) printf("# Creating surface approximation for sphere\n");
    mg = ellipsoid_macrogrid3d(1.0, 1.0, 1.0);
    gr = macro2_bemgrid3d(mg, split, 1);
    break;

  case 'b':
    (void) printf("# Creating surface approximation for box\n");
    mg = box_macrogrid3d(1.0, 1.0, 1.0);
    gr = macro2_bemgrid3d(mg, split, 1);
    break;

  case 'f':
    (void) printf("# Reading grid file \"%s\"\n", filename);
    mg = NULL;
    gr = read_bemgrid3d(filename);
    break;
  }

  if(geometry=='f' && split>0){
    (void) printf("  %d vertices, %d triangles\n",
		  gr->vertices, gr->triangles);
    for(i=0; i<split; i++)
      gr = refine_bemgrid3d(gr);
  }
  (void) printf("  %d vertices, %d triangles\n",
		gr->vertices, gr->triangles);
  prepare_bemgrid3d(gr);

  (void) printf("# Creating cluster tree\n");

  if(basisfunc == HLIB_CONSTANT_BASIS){
    ct = buildcluster_bemgrid3d(gr, HLIB_REGULAR, nmin,0);
  }else{
    ct = buildvertexcluster_bemgrid3d(gr, HLIB_REGULAR, nmin,0);
  }
  root = ct->root;
  clusters = count_cluster(root);
  rows = root->size;
  cols = root->size;
  (void) printf("  %d clusters\n",clusters);

  x0 = (double *) malloc((size_t) sizeof(double) * cols);
  assert(x0 != NULL);
  x1 = (double *) malloc((size_t) sizeof(double) * rows);
  assert(x1 != NULL);
  x2 = (double *) malloc((size_t) sizeof(double) * cols);
  assert(x2 != NULL);

  for(reference=0; reference<2; reference++){
    if(reference==1){
      (void) printf("# Creating reference supermatrix\n");
      order_interpolation++;
      nfdeg++;
      aca_eps *= 0.1;
    }else{
      (void) printf("# Creating supermatrix\n");
    }
    t1 = times(&buffer);
    if(operator_char == 'd'){
      bfactory = new_surfacebemfactory_dlp(gr, basisfunc, ct, basisfunc, ct,
					   nfdeg, nfdeg, order_interpolation, 0.5);
      lower = 0;
    }else{
      bfactory = new_surfacebemfactory_slp(gr, basisfunc, ct, basisfunc, ct,
					   nfdeg, nfdeg, order_interpolation);
      lower = 1;
    }

    bfactory->dist_adaptive = dist_adaptive;
    bfactory->qmin = qmin;
    bfactory->q2->minorder = qmin;
    bfactory->disp = (stdout_is_terminal() ? eta_progressbar : 0);
    switch(compression_method){
    case 0:
      M = onthefly_hca_coarsen_supermatrix(root, root, bfactory,aca_eps*0.9, acak,
					   aca_eps, 1, HLIB_ACA,lower,eta,0);
      break;
    case 1:
      M = onthefly_hca_coarsen_supermatrix(root, root, bfactory,aca_eps*0.01, acak,
					   aca_eps, 1, HLIB_ACAP,lower,eta,0);

      break;
    case 2:
      M = onthefly_hca_coarsen_supermatrix(root, root, bfactory,aca_eps*0.9, acak,
					   aca_eps, 1, HLIB_INTERPOL,lower,eta,0);
      break;
    case 3:
      M = onthefly_hca_coarsen_supermatrix(root, root, bfactory,
					   (operator_char=='d' ? aca_eps/30.0 : aca_eps/4.0),
					   acak, aca_eps, 1, HLIB_HCAI,lower,eta,0);
      break;
    default:
    case 4:
      M = onthefly_hca_coarsen_supermatrix(root, root, bfactory,
					   (operator_char=='d' ? aca_eps/30.0 : aca_eps/4.0),
					   acak, aca_eps, 1, HLIB_HCAII,lower,eta,0);
      /*
        blk = build_blockcluster(root, root, HLIB_MAXADMISSIBILITY,
        HLIB_BLOCK_INHOMOGENEOUS, eta, nmin);
        M = coarsen_hca_from_blockcluster(blk, bfactory,
        (operator_char=='d' ? aca_eps/30.0 : aca_eps/4.0),
        aca_eps,
        1, acak);
      */

      break;
    case 5:
      blk = build_blockcluster(root, root, HLIB_MAXADMISSIBILITY,
			       HLIB_BLOCK_INHOMOGENEOUS, eta, nmin);
      h2r = newrecompression_surfacebem(bfactory, aca_eps, 0.0, acak,
					HLIB_EUCLIDEAN_RELATIVE);
      M = virtual2_supermatrix(root, root, h2r, blk);
      break;
    case 6:
      blk = build_blockcluster(root, root, HLIB_MAXADMISSIBILITY,
			       HLIB_BLOCK_INHOMOGENEOUS, eta, nmin);
      h2r = newrecompression_surfacebem(bfactory, aca_eps, 0.0, acak,
					HLIB_EUCLIDEAN_RELATIVE);
      h2r->cstrategy = HLIB_COEFF_ACA;
      M = virtual2_supermatrix(root, root, h2r, blk);
      break;
    }
    t2 = times(&buffer);

    if(reference==1){
      order_interpolation--;
      nfdeg--;
      aca_eps /= 0.1;
    }else{
    }

    (void) printf("  Fill: %.1f sec.,",((double)t2-t1)/clk_tck);
    size_m = getsize_supermatrix(M);
    (void) printf("  Storage: %.1f MB (%.1f KB/DoF)\n",
		  size_m / 1024.0 / 1024.0, size_m / 1024.0 / rows);

    normfull = norm2_supermatrix(M);
    printf("  |M| = %g\n",normfull);

    if(operator_char=='d'){
      fill_vector(cols, x0, 1.0);
      eval_supermatrix(M, x0, x1);
      norm = norm2_lapack(rows, x1);
      (void) printf("  |(M+I/2)*1| / (|M| |1|) = %g\n",
		    norm/normfull/sqrt((double)root->size));

      if(reference==0){
	Kmod = new_supermatrix(1, 1, rows, cols, NULL,
			       new_rkmatrix(1, rows, cols), NULL);
	fill_vector(rows, Kmod->r->a, sqrt(1.0 / rows));
	fill_vector(cols, Kmod->r->b, sqrt(1.0 / cols));
	Kmod->r->kt=1;
	norm = norm2_supermatrix(M);
	scale_lapack(Kmod->rows, Kmod->r->a,
		     norm / norm2_lapack(rows, Kmod->r->a));
      }
    }else {
      Kmod = 0x0;
    }

    fill_vector(cols, x0, 1.0);

    t1 = times(&buffer);
    maxiter = 1 + 8192 / rows;
    for(i=0; i<maxiter; i++){
      if(operator_char=='d')
	mvm_supermatrix(M, HLIB_EVAL_NOTTRANSPOSED, HLIB_EVAL_DEFAULT,
			1.0, x0, 0.0, x1);
      else
	mvm_supermatrix(M, HLIB_EVAL_NOTTRANSPOSED, HLIB_EVAL_SYMM_LOWER,
			1.0, x0, 0.0, x1);
    }
    t2 = times(&buffer);
    (void) printf("  M*v: %.3f sec.\n",
		  (double) (t2-t1) / clk_tck / maxiter);
    if(output>0 && reference==0){
      (void) printf("# Print matrix structure to \"M.ps\"\n");
      if(output==2){
	outputsvd_supermatrix(M, "M.ps");
      }else{
	outputrank_supermatrix(M, "M.ps");
      }
    }
    if(reference==0) M_compress = M;
    if(reference==1) M_full = M;
  }

  epd = new_eval_powerdata(rows);
  epd->typ = 1;
  epd->e_A.A = M_full;
  epd->e_B.A = M_compress;
  if(operator_char=='s'){
    epd->e_A.mode_A = HLIB_EVAL_SYMM_LOWER;
    epd->e_B.mode_A = HLIB_EVAL_SYMM_LOWER;
  }
  pd = new_powerdata(rows,cols,5);
  pd->data = (void*) epd;
  pd->eval = eval_power_general;
  if(operator_char=='d'){
    pd->filter=1;
    epd->e_A.A_Add = Kmod;
    epd->e_B.A_Add = Kmod;
  }
  (void) printf("  apxerror = %g\n",norm2_power(pd)/normfull);

  printf("# Preparing preconditioner...\n");
  if(operator_char=='d'){
    printf("# Coarsen\n");
    M_precond = coarsenh2_supermatrix(M_compress,0.1,1,0);
    addrk2_supermatrix(Kmod->r, M_precond, 0, 0);
    printf("# LU\n");
    ludecomposition_supermatrix(M_precond);
  }else{
    printf("# Coarsen\n");
    M_precond = coarsenh2low_supermatrix(M_compress,0.05,1,0);
    printf("# Cholesky\n");
    choleskydecomposition_supermatrix(M_precond);
  }

  if(testdis==1){
    printf("# Estimate inversion error:\n");
    epd->typ = 0;
    epd->e_A.A = M_full;
    epd->e_B.A = M_compress;
    if(operator_char=='s')epd->e_B.A_Pre_Cholesky = M_precond;
    if(operator_char=='d')epd->e_B.A_Pre_LU = M_precond;
    (void) printf("  |I-MM^-1| = %g\n",norm2_power(pd));
  }

  /* test approximation for one right-hand side */
  t1 = times(&buffer);
  Mass = buildmassmatrix_surfacebem(1, bfactory);
  t2 = times(&buffer);
  fillfunctional_surfacebem(x2, 1, dirichlet_quadratic, bfactory);
  solve_conjgrad_sparsematrix(Mass, x2, x0, 1e-8, 50,
			      NULL, HLIB_PREC_JACOBI, NULL,
			      HLIB_EVAL_DEFAULT, 0);
  if(operator_char=='d') mvm_supermatrix(M_full, HLIB_EVAL_NOTTRANSPOSED,
                                         HLIB_EVAL_DEFAULT, 1.0, x0, 0.0, x1);
  else mvm_supermatrix(M_full, HLIB_EVAL_NOTTRANSPOSED,
		       HLIB_EVAL_SYMM_LOWER, 1.0, x0, 0.0, x1);
  printf("# Solving system for a quadratic rhs...\n");
  fill_vector(cols, x2, 0.0);
  if(operator_char=='d')
    j = solve_gmres_supermatrix(M_compress, x1, x2, 1e-8, 250, M_precond,
				HLIB_PREC_LU,
				Kmod, HLIB_EVAL_DEFAULT, 0);
  else j = solve_gmres_supermatrix(M_compress, x1, x2, 1e-8, 250, M_precond,
				   HLIB_PREC_CHOLESKY,
				   Kmod, HLIB_EVAL_SYMM_LOWER, 0);
  if(operator_char=='d'){
    norm = 0.0; for(i=0; i<cols; i++) norm += x2[i];
    norm /= (double)cols; for(i=0; i<cols; i++) x2[i] -= norm;
  }
  norm = l2norm_surfacebem(x0,NULL,1,bfactory);
  for(i=0; i<cols; i++) x2[i] -= x0[i];
  printf("  |x - x_H|/|x| = %g (%d steps)\n",
	 l2norm_surfacebem(x2,NULL,1,bfactory)/norm,j);


  h_coarsen = aca_eps*sqrt(10.0);
  mode = HLIB_EVAL_DEFAULT;
  if(operator_char=='s') mode = HLIB_EVAL_SYMM_LOWER;

  for(;h_coarsen > 0 && h_coarsen < 1.1; h_coarsen *= sqrt(10.0)){
    if(h_coarsen > 0.9) h_coarsen = 1.0;
    t1 = times(&buffer);
    if(operator_char=='d')
      M = coarsenh2_supermatrix(M_compress,h_coarsen,1,0);
    else
      M = coarsenh2low_supermatrix(M_compress,h_coarsen,1,0);

    t2 = times(&buffer);
    (void) printf("  Coarsen (eps=%.3f): %.1f sec.,",h_coarsen,
		  (double) (t2-t1) / clk_tck);
    (void) printf(" Storage: %.1f KB/DoF\n",
		  getsize_supermatrix(M) / 1024.0 / rows);
    t1 = times(&buffer);
    maxiter = 1 + 65536 / rows;
    for(i=0; i<maxiter; i++)
      eval_supermatrix(M, x0, x2);
    t2 = times(&buffer);
    (void) printf(" M*v: %.3f sec.\n",
		  (double) (t2-t1) / clk_tck / maxiter);

    epd->typ = 0;
    epd->e_B.A = M;
    epd->e_B.A_Cholesky = 0;
    epd->e_B.A_LU = 0;
    epd->e_B.mode_A = mode;
    if(operator_char=='s')epd->e_B.A_Pre_Cholesky = M_precond;
    if(operator_char=='d')epd->e_B.A_Pre_LU = M_precond;
    (void) printf("  |I-MM^-1| = %g\n",norm2_power(pd));

    if(operator_char=='d') addrk2_supermatrix(Kmod->r, M, 0, 0);
    rk_arithmetics_set_adaptive(0);
    t1 = times(&buffer);
    if(operator_char == 'd') {
      ludecomposition_supermatrix(M);
    } else {
      choleskydecomposition_supermatrix(M);
    }

    t2 = times(&buffer);
    (void) printf("  %.1f seconds for setup\n",(double) (t2-t1) / clk_tck);

    if(operator_char == 'd') {
      epd->e_B.A = 0;
      epd->e_B.A_Cholesky = 0;
      epd->e_B.A_Pre_Cholesky = 0;
      epd->e_B.A_LU = M;
      epd->e_B.A_Pre_LU = M;
      epd->e_A.A_Add = Kmod;
      epd->e_B.A_Add = 0;
      (void) printf("  |I-MM^~1| = %g\n",norm2_power(pd));
      epd->e_A.A_Add = 0;
    } else {
      epd->e_B.A = 0;
      epd->e_B.A_Cholesky = M;
      epd->e_B.A_Pre_Cholesky = M;
      epd->e_B.A_LU = 0;
      epd->e_B.A_Pre_LU = 0;
      (void) printf("  |I-MM^~1| = %g\n",norm2_power(pd));
    }

    fill_vector(cols, x2, 0.0);

    t1 = times(&buffer);
    if(operator_char == 'd') {
      (void) printf("# GMRES solve (LU-Preconditioner)");
      i = solve_gmres_supermatrix(M_compress, x1, x2, 1e-6, 5000, M,
				  HLIB_PREC_LU, Kmod, mode, 0);
    } else {
      (void) printf("# GMRES solve (Cholesky-Preconditioner)");
      i = solve_gmres_supermatrix(M_compress, x1, x2, 1e-6, 5000, M,
				  HLIB_PREC_CHOLESKY, Kmod, mode, 0);
    }

    t2 = times(&buffer);
    (void) printf("  %.1f seconds, %d steps\n", (double) (t2-t1) / clk_tck,i);
  }

  fill_vector(cols, x2, 0.0);
  t1 = times(&buffer);
  (void) printf("# GMRES solve (No Preconditioner)");
  i = solve_gmres_supermatrix(M_compress, x1, x2, 1e-6, 500, 0x0,
			      HLIB_PREC_DEFAULT, Kmod, mode, 0);
  t2 = times(&buffer);
  (void) printf("  %.1f seconds, %d steps\n", (double) (t2-t1) / clk_tck,i);
  fill_vector(cols, x2, 0.0);
  t1 = times(&buffer);

  return 0;
}
