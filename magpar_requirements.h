#ifndef OOMPH_MAGPAR_REQUIREMENTS_H
#define OOMPH_MAGPAR_REQUIREMENTS_H

namespace magpar
{
  /* Macros taken from magpar */

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


  inline double PointFromPlane( const std::vector<double>& x, const std::vector<double>& v1,
			 const std::vector<double>& v2, const std::vector<double>& v3)
  {
    const unsigned ND = 3;

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
