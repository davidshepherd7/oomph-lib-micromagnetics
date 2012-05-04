#ifndef OOMPH_MAGPAR_REQUIREMENTS_H
#define OOMPH_MAGPAR_REQUIREMENTS_H

namespace magpar
{
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

}

#endif
