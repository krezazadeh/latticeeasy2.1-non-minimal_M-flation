/*
MODEL.H file for the TWOFLDM model.
Written by Gary Felder (gfelder@physics.stanford.edu)
Last Modified 3/7/01

This file contains the model specific functions and definitions for the model:
V = 1/2 m^2 phi^2 + 1/2 g^2 phi^2 chi^2
Note that this model is simply a special case of the NFLDM model. Choosing nflds=2 for that model is equivalent to running this one.
Comments: The inflaton potential is taken to be pure m^2 phi^2. The inflaton decays (generally resonantly) to the single field chi.

We use the default rescalings (see below), which give A=1/f0, B=m, r=3/2, s=0
The potential in program units is given by V_{pr} = 1/2 phi_{pr}^2 + 1/2 gm (f0^2/a^3) phi_{pr}^2 chi_{pr}^2
*/

/*
General comments about the model.h file:
This file contains the following functions - all called externally
modelinfo(FILE *info_) outputs information about the model and model-specific parameters to a file.
modelinitialize() performs any model-specific initialization
potential_energy(int term, double *field_values) calculates the average potential energy density term by term. The variable num_potential_terms (just above this function) specifies how many separate potential terms are used in this model.
dvdf(int fld, int i, int j, int k) calculates the potential term in the equation of motion, dV/dfield, for the field fld at the lattice point (i,j,k)
effective_mass(double mass_sq[], double *field_values) calculates the square masses of the fields and puts them in the array mass_sq. The parameter beginning tells the function to use initial field values - if this parameter is zero then the field quantities will be calculated dynamically.
model_output(int flush, char *ext_) allows each model to include its own specialized output function(s). The parameter flush is set to 1 when infrequent calculations are being performed and 0 otherwise. The string ext_ gives the extension for output filenames.
*/

/* - This section should be copied into parameters.h when using this model
// ---Adjustable parameters for TWOFLDM model--- //
const double m2=1.e-12; // Square mass of inflaton
const double gm=2.5e+5; // Resonance parameter g^2/m^2
*/

// Rescaling parameters.
// The program variables ("_pr") are defined as
//   f_pr = rescale_A a^rescale_r f (f=field value)
//   x_pr = rescale_B x (x=distance)
//   dt_pr = rescale_B a^rescale_s dt (t=time)
// The constants beta, cpl, and f0 are used to set the variable rescalings rescale_A, rescale_B, rescale_r, and rescale_s.
// These rescaling constants may be reset independently; their settings in terms of beta, cpl, and f0 are suggestions. See the documentation for more details.
// These rescalings are intrinsic to the model and probably shouldn't be changed for individual runs. Adjustable parameters are stored in "parameters.h"
const double beta=2.; // Exponent of the dominant term in the potential
const double cpl=m2; // Coefficient of the dominant term in the potential (up to numerical factors - see documentation)
// const double f0=0.2220064125323017; // Initial value of phi in Planck units, typically the point at which phi'=0.

// Important:
const double f0=0.03618503456905715; // Initial value of phi in Planck units, typically the point at which phi'=0.

// const double f0=0.02;

// By default these are automatically set to A=1/f0, B=sqrt(cpl) f0^(-1+beta/2), R=6/(2+beta), S=3(2-beta)/(2+beta). They may be adjusted to different values, but the relationship S=2R-3 must be maintained for the program equations to remain correct.
const double rescale_A=1./f0;
const double rescale_B=sqrt(cpl)*pow(f0,-1.+beta/2.);
const double rescale_r=6./(2.+beta);
// The value of S in terms of R SHOULD NOT be changed.
const double rescale_s=2.*rescale_r-3.;

// Other global variables
// The array model_vars is intended to hold any model-specific, non-constant global variables.
// These variables should be initialized in modelinitialize() below
// Even if you're not using any, num_model_vars should be at least 1 to keep some compilers happy.
const int num_model_vars=1;
// Model specific variables: None are defined for this model.
// extern double model_vars[num_model_vars];

// const double g = sqrt(gm*m2);

const double MP = 1.0/sqrt(8.0*pi);

const double mu = 0.0001*MP;
const double xi = 100.0;
const double lambdaeff = 3.791890825833465e-6;

const double m = sqrt(lambdaeff/2.0)*mu;
const double kappaeff = (3.0*lambdaeff*mu)/4.0;
// const double w = -3.0;
const double w = -20.0;

const double a1 = 0.8650264449449371;
const double a2 = -6.101596802925804e-10;
const double a3 = 44675.15750901987;
const double a4 = -2.8090457802473342e-6;
const double a5 = 1.0713648950984436e8;
const double a6 = 0.00008770340690798705;
const double a7 = 2.113453867392008e10;
const double a8 = 0.013473990464323156;
const double a9 = 4.02530086096828e11;
const double a10 = -0.2268006148624882;
const double b1 = -1.0746328227949008e-9;
const double b2 = 107261.80678049699;
const double b3 = -0.000012496561706468416;
const double b4 = 4.862285357438475e8;
const double b5 = -0.000909039918686635;
const double b6 = 1.6244842526172632e11;
const double b7 = 0.1725029929745898;
const double b8 = 4.680184755624326e12;
const double b9 = -2.4778702328840896;
const double b10 = -2.35138402687282e12;

double phiJ(double PHIph);
double dphiJ(double PHIph);
double d2phiJ(double PHIph);

// Macros to make the equations more readable: The values of fld are 0=Phi,1=Chi
#define PHI FIELD(0)
#define CHI FIELD(1)
#define GMFA (gm*pw2(f0)/a/a/a) // The combination g^2 f0^2/(m^2 a^3) occurs frequently

// Model specific details about the run to be output to an information file
inline void modelinfo(FILE *info_)
{
  // Name and description of model
  fprintf(info_,"Two Field M Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2 + 1/2 g^2 phi^2 chi^2\n\n");

  // Model specific parameter values
  fprintf(info_,"gm (g^2/m^2) = %f\n",gm);
  fprintf(info_,"m^2 = %e\n",m2);
}

// Perform any model specific initialization
// This function is called twice, once before initializing the fields (which_call=1) and once after (which_call=2)
inline void modelinitialize(int which_call)
{
  if(which_call==1)
  {
    if(nflds!=2)
    {
      printf("Number of fields for TWOFLDM model must be 2. Exiting.\n");
      exit(1);
    }
  }
}

const int num_potential_terms=2; // Number of terms that are calculated separately in the potential



inline double potential_energy(int term, double *field_values)
{
  DECLARE_INDICES
  double potential=0.;

  double PHIph;
  double CHIph;

  if(field_values==NULL) // If no values are given calculate averages on the lattice
  {
    // Loop over grid to calculate potential term
    LOOP{

        PHIph = (f0*PHI)/pow(a,1.5);
        CHIph = (CHI*f0)/pow(a,1.5);

        // CHIph = 0.0;

      if(term==0){
          potential += (
              (pow(a,3)*((lambdaeff*pw2(mu - phiJ(PHIph))*
              pw2(phiJ(PHIph)))/
              (4.*pw2(1 + (xi*pw2(phiJ(PHIph)))/
              pw2(MP))) +
              (avCHIph2*(pw2(m) +
              2*kappaeff*w*phiJ(PHIph) +
              (lambdaeff*(-1 + w)*w*
              pw2(phiJ(PHIph)))/2.))/
              (2.*(1 + (xi*pw2(phiJ(PHIph)))/
              pw2(MP)))))/(pw2(f0)*m2)
          );
      }
      else if(term==1){
          potential += (
              (pow(a,3)*pw2(CHIph)*
              (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
              (lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph)))/2.))
              /(2.*pw2(f0)*m2*
              (1 + (xi*pw2(phiJ(PHIph)))/pw2(MP)))
          );
      }
    }

    // Convert sum to average
    potential /= gridsize;
  }
  else // If field values are given then use them instead
  {
      PHIph = (f0*field_values[0])/pow(a,1.5);
      CHIph = (field_values[1]*f0)/pow(a,1.5);

      // CHIph = 0.0;

    if(term==0){
      potential = (
          (pow(a,3)*((lambdaeff*pw2(mu - phiJ(PHIph))*
          pw2(phiJ(PHIph)))/
          (4.*pw2(1 + (xi*pw2(phiJ(PHIph)))/
          pw2(MP))) +
          (avCHIph2*(pw2(m) +
          2*kappaeff*w*phiJ(PHIph) +
          (lambdaeff*(-1 + w)*w*
          pw2(phiJ(PHIph)))/2.))/
          (2.*(1 + (xi*pw2(phiJ(PHIph)))/
          pw2(MP)))))/(pw2(f0)*m2)
      );
    }
    else if(term==1){
      potential = (
          (pow(a,3)*pw2(CHIph)*
          (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
          (lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph)))/2.))
          /(2.*pw2(f0)*m2*
          (1 + (xi*pw2(phiJ(PHIph)))/pw2(MP)))
      );
    }
  }

  // Include numerical coefficients
  // if(term==0) // 1/2 m^2 phi^2
  //   potential *= .5;
  // else if(term==1) // 1/2 g^2 phi^2 chi^2
  //   potential *= .5*GMFA;

  return (potential);
}




// Potential terms in the equations of motion, dV/dfield, evaluated at point (i,j,k)
// See documentation for details on the normalization of these terms
inline double dvdf(int fld, INDEXLIST)
{
  // static double gmfa = GMFA; // Since this same value gets used at every point only calculate it once per time
  // static double tprevious = t;
  // if(t!=tprevious)
  //   gmfa = GMFA;

  double PHIph;
  double CHIph;
  double dPHIph;
  double dCHIph;

  PHIph = (f0*PHI)/pow(a,1.5);
  CHIph = (CHI*f0)/pow(a,1.5);
  dPHIph = f0/pow(a,1.5);
  dCHIph = f0/pow(a,1.5);

  // CHIph = 0.0;
  // dCHIph = 0.0;

  if(fld==0) // Phi
  {
      // CHIph = 0.0;
      // dCHIph = 0.0;
    return(
        (pow(a,3)*dPHIph*pw2(MP)*dphiJ(PHIph)*
        (2*(avCHIph2 + pw2(CHIph))*kappaeff*
        pow(MP,4)*w +
        pw2(MP)*(lambdaeff*pw2(MP)*
        (pw2(mu) + avCHIph2*(-1 + w)*w) -
        2*avCHIph2*pw2(m)*xi +
        pw2(CHIph)*
        (lambdaeff*pw2(MP)*(-1 + w)*w -
        2*pw2(m)*xi))*phiJ(PHIph) -
        3*lambdaeff*pow(MP,4)*mu*
        pw2(phiJ(PHIph)) +
        (-2*avCHIph2*pw2(m)*pw2(xi) +
        pw2(CHIph)*xi*
        (lambdaeff*pw2(MP)*(-1 + w)*w -
        2*pw2(m)*xi) +
        lambdaeff*pw2(MP)*
        (2*pw2(MP) - pw2(mu)*xi +
        avCHIph2*(-1 + w)*w*xi))*
        pow(phiJ(PHIph),3) +
        xi*(lambdaeff*pw2(MP)*mu -
        2*avCHIph2*kappaeff*w*xi -
        2*pw2(CHIph)*kappaeff*w*xi)*
        pow(phiJ(PHIph),4)))/
        (2.*pw2(f0)*m2*pow(pw2(MP) +
        xi*pw2(phiJ(PHIph)),3))
     );
 }
  else // Chi
    return(
        (pow(a,3)*CHIph*dCHIph*
        (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
        (lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph)))/2.))
        /(pw2(f0)*(m2 +
        (m2*xi*pw2(phiJ(PHIph)))/pw2(MP)))
     );
}

// Calculate effective mass squared and put it into the array mass_sq[] (used for initial conditions and power spectra)
// See documentation for normalization of these terms
// When setting initial conditions field values will be supplied in the array field_values. Otherwise this function will calculate them on the lattice.
inline void effective_mass(double mass_sq[], double *field_values)
{
  DECLARE_INDICES
  int fld;
  // double fldsqrd[nflds]; // Square value of field
  double mass_sq_temp[nflds]; // Square value of field
  double correction; // Used to adjust masses by the appropriate power of the scale factor. (See documentation.)

  double PHIph;
  double CHIph;
  double dPHIph;
  double dCHIph;

  // Loop over fields to find mean-square value
  if(field_values==NULL){ // If no values are given calculate averages on the lattice
      for(fld=0;fld<nflds;fld++){
          LOOP{
                  PHIph = (f0*PHI)/pow(a,1.5);
                  CHIph = (CHI*f0)/pow(a,1.5);
                  dPHIph = f0/pow(a,1.5);
                  dCHIph = f0/pow(a,1.5);

                  // CHIph = 0.0;
                  // dCHIph = 0.0;

                  if(fld==0){
                      mass_sq_temp[fld] += (
                          (pow(a,3)*(-4*pw2(dPHIph)*lambdaeff*pow(MP,4)*
                          pw2(dphiJ(PHIph))*(mu - phiJ(PHIph))*
                          phiJ(PHIph)*(pw2(MP) -
                          xi*pw2(phiJ(PHIph)))*
                          (pw2(MP) + xi*pw2(phiJ(PHIph))) -
                          4*pw2(CHIph)*pw2(dPHIph)*pw2(MP)*w*
                          xi*pw2(dphiJ(PHIph))*phiJ(PHIph)*
                          (2*kappaeff + lambdaeff*(-1 + w)*phiJ(PHIph))*
                          pw2(pw2(MP) + xi*pw2(phiJ(PHIph)))
                          + lambdaeff*pow(MP,4)*pw2(phiJ(PHIph))*
                          pw2(pw2(MP) + xi*pw2(phiJ(PHIph)))*
                          (pw2(dPHIph)*pw2(dphiJ(PHIph)) +
                          pw2(dPHIph)*d2phiJ(PHIph)*
                          (-mu + phiJ(PHIph))) +
                          pw2(CHIph)*pw2(MP)*w*
                          pow(pw2(MP) + xi*pw2(phiJ(PHIph)),3)*
                          (pw2(dPHIph)*lambdaeff*(-1 + w)*
                          pw2(dphiJ(PHIph)) +
                          pw2(dPHIph)*d2phiJ(PHIph)*
                          (2*kappaeff +
                          lambdaeff*(-1 + w)*phiJ(PHIph))) +
                          2*pw2(CHIph)*pw2(MP)*xi*
                          (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
                          (lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph)))/
                          2.)*(pw2(MP) + xi*pw2(phiJ(PHIph)))*
                          (-(pw2(dPHIph)*pw2(dphiJ(PHIph))*
                          (pw2(MP) - 3*xi*pw2(phiJ(PHIph))))
                          - pw2(dPHIph)*d2phiJ(PHIph)*
                          phiJ(PHIph)*
                          (pw2(MP) + xi*pw2(phiJ(PHIph)))) +
                          lambdaeff*pow(MP,4)*pw2(mu - phiJ(PHIph))*
                          (pw2(dPHIph)*d2phiJ(PHIph)*phiJ(PHIph)*
                          (pow(MP,4) -
                          pw2(xi)*pow(phiJ(PHIph),4)) +
                          pw2(dPHIph)*pw2(dphiJ(PHIph))*
                          (pow(MP,4) -
                          8*pw2(MP)*xi*pw2(phiJ(PHIph)) +
                          3*pw2(xi)*pow(phiJ(PHIph),4))) -
                          avCHIph2*(pw2(MP) +
                          xi*pw2(phiJ(PHIph)))*
                          (4*pw2(dPHIph)*pw2(MP)*w*xi*
                          pw2(dphiJ(PHIph))*phiJ(PHIph)*
                          (2*kappaeff +
                          lambdaeff*(-1 + w)*phiJ(PHIph))*
                          (pw2(MP) + xi*pw2(phiJ(PHIph))) -
                          w*pw2(pow(MP,3) + MP*xi*pw2(phiJ(PHIph)))*
                          (pw2(dPHIph)*lambdaeff*(-1 + w)*
                          pw2(dphiJ(PHIph)) +
                          pw2(dPHIph)*d2phiJ(PHIph)*
                          (2*kappaeff +
                          lambdaeff*(-1 + w)*phiJ(PHIph))) -
                          2*pw2(MP)*xi*
                          (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
                          (lambdaeff*(-1 + w)*w*
                          pw2(phiJ(PHIph)))/2.)*
                          (-(pw2(dPHIph)*pw2(dphiJ(PHIph))*
                          (pw2(MP) -
                          3*xi*pw2(phiJ(PHIph)))) -
                          pw2(dPHIph)*d2phiJ(PHIph)*phiJ(PHIph)*
                          (pw2(MP) + xi*pw2(phiJ(PHIph)))))
                          ))/(2.*pw2(f0)*m2*
                          pow(pw2(MP) + xi*pw2(phiJ(PHIph)),4))
                      );
                  }
                  else if(fld==1){
                      mass_sq_temp[fld] += (
                          (pow(a,3)*pw2(dCHIph)*pw2(MP)*
                          (2*pw2(m) + 4*kappaeff*w*phiJ(PHIph) +
                          lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph))))/
                          (2.*pw2(f0)*m2*(pw2(MP) +
                          xi*pw2(phiJ(PHIph))))
                      );
                  }
            }
            mass_sq_temp[fld] /= gridsize;
      }
  }
  else{ // If field values are given then use them instead
      PHIph = (f0*field_values[0])/pow(a,1.5);
      CHIph = (field_values[1]*f0)/pow(a,1.5);
      dPHIph = f0/pow(a,1.5);
      dCHIph = f0/pow(a,1.5);

      // CHIph = 0.0;
      // dCHIph = 0.0;

      for(fld=0;fld<nflds;fld++){
          if(fld==0){
              mass_sq_temp[fld] += (
                  (pow(a,3)*(-4*pw2(dPHIph)*lambdaeff*pow(MP,4)*
                  pw2(dphiJ(PHIph))*(mu - phiJ(PHIph))*
                  phiJ(PHIph)*(pw2(MP) -
                  xi*pw2(phiJ(PHIph)))*
                  (pw2(MP) + xi*pw2(phiJ(PHIph))) -
                  4*pw2(CHIph)*pw2(dPHIph)*pw2(MP)*w*
                  xi*pw2(dphiJ(PHIph))*phiJ(PHIph)*
                  (2*kappaeff + lambdaeff*(-1 + w)*phiJ(PHIph))*
                  pw2(pw2(MP) + xi*pw2(phiJ(PHIph)))
                  + lambdaeff*pow(MP,4)*pw2(phiJ(PHIph))*
                  pw2(pw2(MP) + xi*pw2(phiJ(PHIph)))*
                  (pw2(dPHIph)*pw2(dphiJ(PHIph)) +
                  pw2(dPHIph)*d2phiJ(PHIph)*
                  (-mu + phiJ(PHIph))) +
                  pw2(CHIph)*pw2(MP)*w*
                  pow(pw2(MP) + xi*pw2(phiJ(PHIph)),3)*
                  (pw2(dPHIph)*lambdaeff*(-1 + w)*
                  pw2(dphiJ(PHIph)) +
                  pw2(dPHIph)*d2phiJ(PHIph)*
                  (2*kappaeff +
                  lambdaeff*(-1 + w)*phiJ(PHIph))) +
                  2*pw2(CHIph)*pw2(MP)*xi*
                  (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
                  (lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph)))/
                  2.)*(pw2(MP) + xi*pw2(phiJ(PHIph)))*
                  (-(pw2(dPHIph)*pw2(dphiJ(PHIph))*
                  (pw2(MP) - 3*xi*pw2(phiJ(PHIph))))
                  - pw2(dPHIph)*d2phiJ(PHIph)*
                  phiJ(PHIph)*
                  (pw2(MP) + xi*pw2(phiJ(PHIph)))) +
                  lambdaeff*pow(MP,4)*pw2(mu - phiJ(PHIph))*
                  (pw2(dPHIph)*d2phiJ(PHIph)*phiJ(PHIph)*
                  (pow(MP,4) -
                  pw2(xi)*pow(phiJ(PHIph),4)) +
                  pw2(dPHIph)*pw2(dphiJ(PHIph))*
                  (pow(MP,4) -
                  8*pw2(MP)*xi*pw2(phiJ(PHIph)) +
                  3*pw2(xi)*pow(phiJ(PHIph),4))) -
                  avCHIph2*(pw2(MP) +
                  xi*pw2(phiJ(PHIph)))*
                  (4*pw2(dPHIph)*pw2(MP)*w*xi*
                  pw2(dphiJ(PHIph))*phiJ(PHIph)*
                  (2*kappaeff +
                  lambdaeff*(-1 + w)*phiJ(PHIph))*
                  (pw2(MP) + xi*pw2(phiJ(PHIph))) -
                  w*pw2(pow(MP,3) + MP*xi*pw2(phiJ(PHIph)))*
                  (pw2(dPHIph)*lambdaeff*(-1 + w)*
                  pw2(dphiJ(PHIph)) +
                  pw2(dPHIph)*d2phiJ(PHIph)*
                  (2*kappaeff +
                  lambdaeff*(-1 + w)*phiJ(PHIph))) -
                  2*pw2(MP)*xi*
                  (pw2(m) + 2*kappaeff*w*phiJ(PHIph) +
                  (lambdaeff*(-1 + w)*w*
                  pw2(phiJ(PHIph)))/2.)*
                  (-(pw2(dPHIph)*pw2(dphiJ(PHIph))*
                  (pw2(MP) -
                  3*xi*pw2(phiJ(PHIph)))) -
                  pw2(dPHIph)*d2phiJ(PHIph)*phiJ(PHIph)*
                  (pw2(MP) + xi*pw2(phiJ(PHIph)))))
                  ))/(2.*pw2(f0)*m2*
                  pow(pw2(MP) + xi*pw2(phiJ(PHIph)),4))
              );
          }
          else if(fld==1){
              mass_sq_temp[fld] += (
                  (pow(a,3)*pw2(dCHIph)*pw2(MP)*
                  (2*pw2(m) + 4*kappaeff*w*phiJ(PHIph) +
                  lambdaeff*(-1 + w)*w*pw2(phiJ(PHIph))))/
                  (2.*pw2(f0)*m2*(pw2(MP) +
                  xi*pw2(phiJ(PHIph))))
              );
          }
      }
  }

  mass_sq[0] = mass_sq_temp[0];
  mass_sq[1] = mass_sq_temp[1];

  // Put in scale factor correction. This calculation should be the same for all models.
  if(expansion>0) // If there's no expansion don't bother with this.
  {
    correction = pow(a,2.*rescale_s+2.);
    for(fld=0;fld<nflds;fld++)
      mass_sq[fld] *= correction;
  }
}

// Model-specific output functions
// inline void model_output(int flush,char *ext_){}

inline double phiJ(double PHIph) {
    return (
        (a1*PHIph + a2*pw2(PHIph) + a3*pow(PHIph,3) +
        a4*pow(PHIph,4) + a5*pow(PHIph,5) + a6*pow(PHIph,6) +
        a7*pow(PHIph,7) + a8*pow(PHIph,8) + a9*pow(PHIph,9) +
        a10*pow(PHIph,10))/
        (1 + b1*PHIph + b2*pw2(PHIph) + b3*pow(PHIph,3) +
        b4*pow(PHIph,4) + b5*pow(PHIph,5) + b6*pow(PHIph,6) +
        b7*pow(PHIph,7) + b8*pow(PHIph,8) + b9*pow(PHIph,9) +
        b10*pow(PHIph,10))
    );
}

inline double dphiJ(double PHIph) {
    return (
        (1 + (xi*pw2(phiJ(PHIph)))/pw2(MP))/
        sqrt(1 + (xi*(1 + 6*xi)*pw2(phiJ(PHIph)))/
        pw2(MP))
    );
}

inline double d2phiJ(double PHIph) {
    return (
        (xi*phiJ(PHIph)*(pw2(MP) +
        xi*pw2(phiJ(PHIph)))*
        (pw2(MP)*(1 - 6*xi) +
        xi*(1 + 6*xi)*pw2(phiJ(PHIph))))/
        pw2(pow(MP,3) + MP*xi*(1 + 6*xi)*pw2(phiJ(PHIph)))
    );
}

///////////// Model-specific output functions///////////////////

//// Calculate gravity waves ////
// These functions must be declared here so the compiler will recognize them in the function below
void fftrn(double f[], double fnyquist[], int ndims, int size[], int forward); // Do a Fourier transform of an ndims dimensional array of real numbers

// Increments a grid location accounting for periodic wrapping
inline int incr(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}
// Decrements a grid location accounting for periodic wrapping
inline int decr(int i)
{
  return( (i==0) ? N-1 : i-1 );
}

//this takes the derivative of the field at a point of a field in a particular direction
inline double deriv(int i, int j, int k, int fld, int direction)
{
  if(direction==1)
    return((f[fld][incr(i)][j][k]-f[fld][decr(i)][j][k])/(2*dx));
  else if(direction==2)
    return((f[fld][i][incr(j)][k]-f[fld][i][decr(j)][k])/(2*dx));
  else if(direction==3)
    return((f[fld][i][j][incr(k)]-f[fld][i][j][decr(k)])/(2*dx));
  else
  {
    printf("Invalid direction! Must be 1, 2, or 3. Exiting.");
    exit(1);
  }
}

#if SMULTIDIREC
static double cosIntegrand[9][2][N/2+1][2], sinIntegrand[9][2][N/2+1][2]; //x=0, y=1, z=2, xy=3, yz=4, xz=5, -xy=6, -yz=7, -xz=8
static double q11[N][N][N], q12[N][N][N], q13[N][N][N]; //these q arrays are only necessary if all the directions are calculated
#else
static double cosIntegrand[1][2][N/2+1][2], sinIntegrand[1][2][N/2+1][2]; //x=0
#endif
static double fnyquist[N][2*N]; // Used by FFT routine to store modes with k=nyquist
static double q22[N][N][N], q23[N][N][N], q33[N][N][N]; //the arrays that stores the fourier transform of the derivatives of the field (used in calculating pi_ab)

// Calculate metric perturbations h11(p), h13(p), h22(p) and h23(p).
inline void gravity_waves(int flush ,char *ext_)
{
  int i,j,k;
  int arraysize[] = {N,N,N}; // Array of grid size in all dimensions - used by FFT routine
  int fld; // Array for looping over fields
  double px; // Components of momentum p in correct (program) units
  double dp = 2.*pi/L; // Size of grid spacing in momentum space
  double insidecoeff, outsidecoeff; // Coefficients inside and outside the integrand
  double pi_ab[9][2][2];  // The complex term Pi_ab(t,p1), re=0, im=1.
  double deriv1, deriv2, deriv3; //used to calculate q array
  double fmode, xk, omegah2; // Final calculated quantities
  double tempcoeff1,tempcoeff2; //insidecoeff*cosine (1), insidecoeff*sine (2)
  int numdirections = (SMULTIDIREC ? 9 : 1); // Number of directions in k space being calculated (1 or 9 depending on the parameter SMULTIDIREC)

  static FILE *gravtimes_, *gravitywaves_[9]; //pointers to files that record data
  char name_[500];

  static double tlast=t; // Records the last time at which the quantities here were calculated

  double astar = exp(5.);
  double w = 0.0;

  double Hphys = ad*rescale_B*pow(a,rescale_s-1.);

  if(tlast==t)
  {
    for(j=0;j<numdirections;j++) //creating and opening the gravity files
    {
      sprintf(name_,"gravity%d%s", j, ext_);
      gravitywaves_[j] = fopen(name_,"w");
    } //end for over numdirections

    sprintf(name_,"gravtimes%s",ext_); //creating and opening the gravtimes file (containes the flush times)
    gravtimes_ = fopen(name_,"w");

    for(j=0;j<numdirections;j++)
    {
      for(k=0;k<=1;k++)
      {
	for(i=0;i<=N/2;i++) //for loop to initialize the real and imaginary integrand arrays
	{
	  cosIntegrand[j][k][i][0]=0;
	  cosIntegrand[j][k][i][1]=0;
	  sinIntegrand[j][k][i][0]=0;
	  sinIntegrand[j][k][i][1]=0;
	} //end for over i
      } //end for over k
    } //end for over numdirections
    return; // The first time this function is called simply record the time and exit.
  }//end if

  insidecoeff = (t-tlast)/a; // The (t-tlast) term is delta_t, for the integration (coefficient inside the integral)

  for(fld=0;fld<nflds;fld++)
  {
    ////Calculate q array////

    LOOP //loop over the grid to initialize the q arrays
    {
      deriv1=deriv(i,j,k,fld,1);
      deriv2=deriv(i,j,k,fld,2);
      deriv3=deriv(i,j,k,fld,3);

#if SMULTIDIREC//compiler conditional
      q11[i][j][k]=pw2(deriv1);
      q12[i][j][k]=deriv1*deriv2;
      q13[i][j][k]=deriv1*deriv3;
#endif

      q22[i][j][k]=pw2(deriv2);
      q23[i][j][k]=deriv2*deriv3;
      q33[i][j][k]=pw2(deriv3);
    } //end LOOP

    //taking the fourier transform of each element in q
#if SMULTIDIREC //compiler conditional
    fftrn((double *)q11, (double *)fnyquist, NDIMS, arraysize, 1);
    fftrn((double *)q12, (double *)fnyquist, NDIMS, arraysize, 1);
    fftrn((double *)q13, (double *)fnyquist, NDIMS, arraysize, 1);
#endif
    fftrn((double *)q22, (double *)fnyquist, NDIMS, arraysize, 1);
    fftrn((double *)q23, (double *)fnyquist, NDIMS, arraysize, 1);
    fftrn((double *)q33, (double *)fnyquist, NDIMS, arraysize, 1);

    ///// Loop over the x component of p (i) for the x axis and y component of p (i in q13) for the y axis /////
    for (i=1;i<N/2;i++)
    {
      px = dp*(double)i;

     ///// Calculate components of the tensor h_ab /////

     ///// Calculate Pi_ab /////
      //x-axis
      //pi_ab[0][0][0,1] ==> pi22x
      //pi_ab[0][1][0,1] ==> pi23x
      pi_ab[0][0][0] = .5*(q22[i][0][0]-q33[i][0][0]);
      pi_ab[0][0][1] = .5*(q22[i][0][1]-q33[i][0][1]);
      pi_ab[0][1][0] = q23[i][0][0];
      pi_ab[0][1][1] = q23[i][0][1];

#if SMULTIDIREC //compiler conditional
      //y-axis
      //pi_ab[1][0][0,1] ==> pi11y
      //pi_ab[1][1][0,1] ==> pi13y
      pi_ab[1][0][0]=.5*(q11[0][i][0]-q33[0][i][0]);
      pi_ab[1][0][1]=.5*(q11[0][i][1]-q33[0][i][1]);
      pi_ab[1][1][0]=q13[0][i][0];
      pi_ab[1][1][1]=q13[0][i][1];

      //z-axis
      //pi_ab[2][0][0,1] ==> pi11z
      //pi_ab[2][1][0,1] ==> pi12z
      pi_ab[2][0][0]=.5*(q11[0][0][2*i]-q22[0][0][2*i]);
      pi_ab[2][0][1]=.5*(q11[0][0][2*i+1]-q22[0][0][2*i+1]);
      pi_ab[2][1][0]=q12[0][0][2*i];
      pi_ab[2][1][1]=q12[0][0][2*i+1];

      //xy diagonal
      //pi_ab[3][0][0,1] ==> pi13xy
      //pi_ab[3][1][0,1] ==> pi11xy
      pi_ab[3][0][0] = .5*(q13[i][i][0]-q23[i][i][0]);
      pi_ab[3][0][1] = .5*(q13[i][i][1]-q23[i][i][1]);
      pi_ab[3][1][0] = .125*(q11[i][i][0]+q22[i][i][0])-.25*(q12[i][i][0]+q33[i][i][0]);
      pi_ab[3][1][1] = .125*(q11[i][i][1]+q22[i][i][1])-.25*(q12[i][i][1]+q33[i][i][1]);

      //yz diagonal
      //pi_ab[4][0][0,1] ==> pi13yz
      //pi_ab[4][1][0,1] ==> pi22yz
      pi_ab[4][0][0]=.5*(q13[0][i][2*i]-q12[0][i][2*i]);
      pi_ab[4][0][1]=.5*(q13[0][i][2*i+1]-q12[0][i][2*i+1]);
      pi_ab[4][1][0]=.125*(q22[0][i][2*i]+q33[0][i][2*i])-.25*(q11[0][i][2*i]+q23[0][i][2*i]);//changed .5->.125 and fixed misplaced parentheses
      pi_ab[4][1][1]=.125*(q22[0][i][2*i+1]+q33[0][i][2*i+1])-.25*(q11[0][i][2*i+1]+q23[0][i][2*i+1]);//changed .5->.125 and fixed misplaced parentheses

      //xz diagonal
      //pi_ab[5][0][0,1] ==> pi12xz
      //pi_ab[5][1][0,1] ==> pi11xz
      pi_ab[5][0][0]=.5*(q12[i][0][2*i]-q23[i][0][2*i]);
      pi_ab[5][0][1]=.5*(q12[i][0][2*i+1]-q23[i][0][2*i+1]);
      pi_ab[5][1][0]=.125*(q11[i][0][2*i]+q33[i][0][2*i])-.25*(q13[i][0][2*i]+q22[i][0][2*i]);
      pi_ab[5][1][1]=.125*(q11[i][0][2*i+1]+q33[i][0][2*i+1])-.25*(q13[i][0][2*i+1]+q22[i][0][2*i+1]);

      //-xy diagonal
      //pi_ab[6][0][0,1] ==> pi13
      //pi_ab[6][1][0,1] ==> pi11
      pi_ab[6][0][0]=.5*(q13[N-i][i][0]+q23[N-i][i][0]);
      pi_ab[6][0][1]=.5*(q13[N-i][i][1]+q23[N-i][i][1]);
      pi_ab[6][1][0]=.125*(q11[N-i][i][0]+q22[N-i][i][0])+.25*(q12[N-i][i][0]-q33[N-i][i][0]);
      pi_ab[6][1][1]=.125*(q11[N-i][i][1]+q22[N-i][i][1])+.25*(q12[N-i][i][1]-q33[N-i][i][1]);

      /////had -yz and -xz switched, now fixed!
      //-yz diagonal
      //pi_ab[7][0][0,1] ==> pi13
      //pi_ab[7][1][0,1] ==> pi22
      pi_ab[7][0][0]=.5*(q12[0][N-i][2*i]+q13[0][N-i][2*i]);
      pi_ab[7][0][1]=.5*(q12[0][N-i][2*i+1]+q13[0][N-i][2*i+1]);
      pi_ab[7][1][0]=.25*(q23[0][N-i][2*i]-q11[0][N-i][2*i])+.125*(q22[0][N-i][2*i]+q33[0][N-i][2*i]);//changed .5->.125
      pi_ab[7][1][1]=.25*(q23[0][N-i][2*i+1]-q11[0][N-i][2*i+1])+.125*(q22[0][N-i][2*i+1]+q33[0][N-i][2*i+1]);//changed .5->.125

      //-xz diagonal
      //pi_ab[8][0][0,1] ==> pi12
      //pi_ab[8][1][0,1] ==> pi11
      pi_ab[8][0][0]=.5*(q12[N-i][0][2*i]+q23[N-i][0][2*i]);
      pi_ab[8][0][1]=.5*(q12[N-i][0][2*i+1]+q23[N-i][0][2*i+1]);
      pi_ab[8][1][0]=.125*(q11[N-i][0][2*i]+q33[N-i][0][2*i])+.25*(q13[N-i][0][2*i]-q22[N-i][0][2*i]);
      pi_ab[8][1][1]=.125*(q11[N-i][0][2*i+1]+q33[N-i][0][2*i+1])+.25*(q13[N-i][0][2*i+1]-q22[N-i][0][2*i+1]);

#endif//end if(SMULTIDIREC)

      //first dimension of the cosIntegrand or sinIntegrand corresponds to the direction
      //second dimension signifies which energy density array (h array) it will correspond to
      //third dimension is the k value
      //fourth dimension signifies whether that array is the real (=0) or imaginary (=1) part of the complex number
      for(j=0;j<numdirections;j++)
      {
	if(j==3) px *= sqrt(2.0);
	tempcoeff1 = insidecoeff*cos(px*t);
	tempcoeff2 = insidecoeff*sin(px*t);
	for(k=0;k<=1;k++)
	{
	  cosIntegrand[j][k][i][0] += tempcoeff1*pi_ab[j][k][0];
	  cosIntegrand[j][k][i][1] += tempcoeff1*pi_ab[j][k][1];
	  sinIntegrand[j][k][i][0] += tempcoeff2*pi_ab[j][k][0];
	  sinIntegrand[j][k][i][1] += tempcoeff2*pi_ab[j][k][1];
	} //end of for loop over numdirections
      } //end of for loop over k
    } // End of loop over i
  } // End of loop over fld

  if(flush) // output gravitywaves infrequently
  {
    fprintf(gravtimes_,"%f\n",t);
    outsidecoeff = pw2(rescale_B)/pw2(pw2(rescale_A))*L*L*L/pow((double)N,6)*16*pi*pi; //coefficients outside the integral

    for(i = 1; i < N/2; i++)
    {
      px = dp*(double)i;
      // The following formulas are only valid for lambda phi^4!!!
      // f (the physical wave number of the mode today)
      // fmode = (4.0e+10*px*rescale_B);
      // fmode = (4.e10*rescale_B*px)/(pow(a/astar,(3.*(1. + w))/4.)*astar*pow((pow(a,-2.*rescale_r + 2.*rescale_s)*pow(rescale_B,2.)*totalenergy)/pow(rescale_A,2.),0.25));
      fmode = (1.3434181619042972e11*rescale_B*px)/sqrt(Hphys);
      for(j=0;j<numdirections;j++)
      {
	if(j==3) //the value of px and fmode are sqrt(2) times larger for the diagonal than for the axes
	{
	  px *= sqrt(2.0);
	  fmode *= sqrt(2.0);
	}
	// X_k - The factors of 2 for the axes accounts for the fact that we are only calculating two components aand using the symmatries of the matricies to calculate the rest. the diagonals need a factor of 4 and 8 because different elements in their matrices go to zero. A more detailed derivation can be found in the runsMultDiag.txt file.
	if(j<3){ //the axes
	  xk = outsidecoeff*px*px*px*2*(pw2(cosIntegrand[j][0][i][0])+pw2(cosIntegrand[j][0][i][1])
				       +pw2(sinIntegrand[j][0][i][0])+pw2(sinIntegrand[j][0][i][1])
				       +pw2(cosIntegrand[j][1][i][0])+pw2(cosIntegrand[j][1][i][1])
				       +pw2(sinIntegrand[j][1][i][0])+pw2(sinIntegrand[j][1][i][1]));
	}
	else{ //the diagonals
	  xk = outsidecoeff*px*px*px*4*(pw2(cosIntegrand[j][0][i][0])+pw2(cosIntegrand[j][0][i][1])
				       +pw2(sinIntegrand[j][0][i][0])+pw2(sinIntegrand[j][0][i][1])
				    +2*(pw2(cosIntegrand[j][1][i][0])+pw2(cosIntegrand[j][1][i][1])
				       +pw2(sinIntegrand[j][1][i][0])+pw2(sinIntegrand[j][1][i][1])));
	}
	// Omega_gw h^2
	// omegah2 = 9.3e-6*xk/pow((2*pi),3);
    // omegah2 = (9.3e-6*pow(a,-4. + 2.*rescale_r - 2.*rescale_s)*pow(rescale_A,2.)*pow(a/astar,1. - 3.*w)*xk)/(pow(rescale_B,2.)*totalenergy);
    // omegah2 = (9.3e-6*xk/pow((2*pi),3))/((3./(8.*pi))*pow(a,4.)*pw2(Hphys));
    // omegah2 = 1.0e-37*(9.3e-6*xk/pow((2*pi),3))/((3./(8.*pi))*pow(a,4.)*pw2(Hphys)*pow(rescale_B,4));
    // omegah2 = 1.0e-3*(9.3e-6*xk/pow((2*pi),3))/((3./(8.*pi))*pow(a,4.)*pw2(Hphys));
    omegah2 = 1.0e-3*(9.3e-6*xk/pow((2*pi),3))/((3./(8.*pi))*pow(a,4.)*pw2(Hphys));
	fprintf(gravitywaves_[j],"%e\t%e\n",fmode,omegah2);
      }
    }
    for(j=0;j<numdirections;j++)
    {
      fprintf(gravitywaves_[j],"\n\n");
      fflush(gravitywaves_[j]);
    }
    fflush(gravtimes_);
  }

  tlast=t; // Used to keep track of time interval for integration
} //end gravitywaves function

// Outputs the rms value <sqrt[sigma_i^2]>
inline void sigmarms(int flush, char *ext_)
{
  static FILE *sigmarms_;
  char name_[500];
  DECLARE_INDICES
  int fld;
  double sigmarms,sigmasq;

  static int first=1;
  if(first) // Open output files
  {
    sprintf(name_,"sigmarms%s",ext_);
    sigmarms_=fopen(name_,mode_);
    first=0;
  }

  sigmarms=0.;
  // Calculate rms
  LOOP
  {
    sigmasq=0.;
    for(fld=1;fld<nflds;fld++)
      // sigmasq += pw2(SIGMAF);
      sigmasq += pw2(CHI);
    sigmarms += sqrt(sigmasq);
  }
  sigmarms = sigmarms/(double)gridsize; // Convert sum to average
  fprintf(sigmarms_,"%f %e\n",t,sigmarms*rescaling);
  if(flush)
    fflush(sigmarms_);
}

inline void model_output(int flush,char *ext_)
{
  if(sgravitywaves && t>tgravitywaves)
    gravity_waves(flush,ext_);
  if(ssigmarms)
    sigmarms(flush,ext_);
}

// #undef PHI
// #undef SIGMA

#undef PHI
#undef CHI
