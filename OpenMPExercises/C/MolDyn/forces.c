
/*
 *  Compute forces and accumulate the virial and the potential
 */
  extern double epot, vir;

  void
  forces(int npart, double x[], double f[], double side, double rcoff){

    vir    = 0.0;
    epot   = 0.0;

    double sideh  = 0.5*side;
    double rcoffs = rcoff*rcoff;

    for (int i=0; i<npart*3; i+=3) {
      double xi  = x[i];
      double yi  = x[i+1];
      double zi  = x[i+2];
      double fxi = 0.0;
      double fyi = 0.0;
      double fzi = 0.0;

      for (int j=i+3; j<npart*3; j+=3) {
        double xx = xi-x[j];
        double yy = yi-x[j+1];
        double zz = zi-x[j+2];
        if (xx<-sideh) xx += side;
        if (xx> sideh) xx -= side;
        if (yy<-sideh) yy += side;
        if (yy> sideh) yy -= side;
        if (zz<-sideh) zz += side;
        if (zz> sideh) zz -= side;
        double rd = xx*xx+yy*yy+zz*zz;

        if (rd<=rcoffs) {
          double rrd      = 1.0/rd;
          double rrd3     = rrd*rrd*rrd;
          double rrd6     = rrd3*rrd3;
          double r148     = rrd6*rrd-0.5*rrd3*rrd;

          epot    += rrd6-rrd3;
          vir     -= rd*r148;
          
	  fxi     += xx*r148;
          fyi     += yy*r148;
          fzi     += zz*r148;

          f[j]    -= xx*r148;
          f[j+1]  -= yy*r148;
          f[j+2]  -= zz*r148;
        }
      }
      f[i]     += fxi;
      f[i+1]   += fyi;
      f[i+2]   += fzi;
    }
  }
