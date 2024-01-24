#include "atmosphere.h"
#include "interp2d.h"
#include "fvect.h"


const double EARTHRADIUS = 6360e3;
const double ATMOSRADIUS = 6460e3;
const double ATMOSHEIGHT = ATMOSRADIUS - EARTHRADIUS;
const double ATMOSRADIUS_M = 6395e3;
const unsigned int TRANSMITTANCE_W = 256;
const unsigned int TRANSMITTANCE_H = 64;

// Scale heights (m)
// Rayleigh scattering
const double HR_MS = 8820;          // Midlatitude summer
const double HR_MW = 8100;          // Midlatitude winter
const double HR_SS = 8550;          // Subarctic summer
const double HR_SW = 7600;          // Subarctic winter
const double HR_T  = 9050;           // Tropical
// Mie scattering
const double HMS_CC = 2385;         // Continental clean
// Mie absorption
const double HMA_CC = 4300;         // Continental clean

const double SOLOMG = 6.7967e-5;    // .533 apex angle
const int WVLSPAN = 400;            // 380nm to 780nm

// Aerosol optical depth at 550nm for continental clean
const double AOD_CC_550 = 0.05352;

/** 380nm to 780nm at 20nm intervals **/
// Extraterrestrial solar W/m^2/nm
const SCOLOR EXTSOL = {1.11099345, 1.75363199, 1.67126921, 1.97235667, 2.01863015,
    1.93612482, 1.87310758, 1.88792588, 1.86274213, 1.84433248, 1.80370942, 1.73020456,
    1.67036654, 1.57482149, 1.53646383, 1.45515319, 1.38207435, 1.3265141 , 1.26648111,
    1.20592043};

// Rayleight scattering coefficients at sea level in m^-1

const SCOLOR BR0_MS = {4.63889967e-05, 3.76703293e-05, 3.09193474e-05, 2.56223081e-05,
    2.14198405e-05, 1.80483398e-05, 1.53177709e-05, 1.30867298e-05, 1.12487269e-05,
    9.72370830e-06, 8.44932642e-06, 7.37769845e-06, 6.47116653e-06, 5.70005327e-06,
    5.04091603e-06, 4.47430240e-06, 3.98549839e-06, 3.56178512e-06, 3.19293761e-06,
    2.87072599e-06};
const SCOLOR BR0_MW = {5.03857564e-05, 4.09159104e-05, 3.35832809e-05, 2.78298619e-05,
    2.32653203e-05, 1.96033395e-05, 1.66375116e-05, 1.42142496e-05, 1.22178890e-05,
    1.05614786e-05, 9.17729921e-06, 8.01334246e-06, 7.02870602e-06, 6.19115557e-06,
    5.47522872e-06, 4.85979708e-06, 4.32887894e-06, 3.86865960e-06, 3.46803311e-06,
    3.11806054e-06};
const SCOLOR BR0_SS = {4.73789047e-05, 3.84741872e-05, 3.15791442e-05, 2.61690698e-05,
    2.18769246e-05, 1.84334785e-05, 1.56446412e-05, 1.33659912e-05, 1.14887667e-05,
    9.93120528e-06, 8.62962901e-06, 7.53513326e-06, 6.60925660e-06, 5.82168833e-06,
    5.14848557e-06, 4.56978082e-06, 4.07054608e-06, 3.63779107e-06, 3.26107261e-06,
    2.93198523e-06};
const SCOLOR BR0_SW = {5.30623659e-05, 4.30894595e-05, 3.53673035e-05, 2.93082494e-05,
    2.45012287e-05, 2.06447149e-05, 1.75213353e-05, 1.49693438e-05, 1.28669320e-05,
    1.11225291e-05, 9.66481886e-06, 8.43902999e-06, 7.40208736e-06, 6.52004427e-06,
    5.76608570e-06, 5.11796089e-06, 4.55883914e-06, 4.07417187e-06, 3.65226316e-06,
    3.28369923e-06};
const SCOLOR BR0_T = {4.55376661e-05, 3.69790036e-05, 3.03519157e-05, 2.51520876e-05,
    2.10267438e-05, 1.77171168e-05, 1.50366593e-05, 1.28465622e-05, 1.10422904e-05,
    9.54525887e-06, 8.29426444e-06, 7.24230298e-06, 6.35240772e-06, 5.59544593e-06,
    4.94840517e-06, 4.39219003e-06, 3.91235655e-06, 3.49641927e-06, 3.13434084e-06,
    2.81804244e-06};

// OPAC Mie scattering coefficient at sea level in m^-1
// continental clean
const SCOLOR BM0_CC = {2.8278e-05, 2.6766e-05, 2.5349e-05, 2.3945e-05, 2.2771e-05,
    2.1593e-05, 2.0508e-05, 1.9494e-05, 1.8466e-05, 1.7606e-05, 1.6766e-05,
    1.5987e-05, 1.5268e-05, 1.4521e-05, 1.3935e-05, 1.3328e-05, 1.2748e-05,
    1.2181e-05, 1.1641e-05, 1.112e-05};
const SCOLOR BM0_CA = {7.8237e-05, 7.3966e-05, 6.9974e-05, 6.6022e-05,
    6.2729e-05, 5.9427e-05, 5.6391e-05, 5.3561e-05, 5.0694e-05, 4.8299e-05,
    4.596e-05, 4.3797e-05, 4.1802e-05, 3.9726e-05, 3.8103e-05, 3.6423e-05,
    3.4817e-05, 3.3253e-05, 3.1763e-05, 3.033e-05};
const SCOLOR AM0_CC = {7.12e-7, 6.83e-07, 6.56e-07, 6.3e-07, 6.09e-07, 5.87e-07,
    5.79e-07, 5.84e-07, 5.89e-07, 5.71e-07, 5.54e-07, 5.47e-07, 5.51e-07, 5.54e-07,
    5.4e-07, 5.26e-07, 5.24e-07, 5.35e-07, 5.45e-07, 5.52e-07};
const SCOLOR AM0_CA = {7.671e-06, 7.3e-06, 6.961e-06, 6.625e-06, 6.347e-06,
    6.068e-06, 5.839e-06, 5.654e-06, 5.466e-06, 5.262e-06, 5.063e-06, 4.909e-06,
    4.8e-06, 4.685e-06, 4.538e-06, 4.386e-06, 4.279e-06, 4.218e-06, 4.159e-06, 4.104e-06};


static void
swap12(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}


static double
linear_interp(double x, const double *xp, const double *fp, int size)
{
    if (x <= xp[0]) {
        return fp[0];
    } else if (x >= xp[size - 1]) {
        return fp[size - 1];
    } else {
        int j = 0;
        while (x > xp[j + 1] && j < size - 2) j++;

        double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
        return fp[j] + slope * (x - xp[j]);
    }
}


static int
hit_sphere(const double r, const double ct, const double radius, double *roots)
{
    float A = 1;
    float B = -2 * r * ct;
    float C = r * r - radius * radius;

    if (!quadratic(roots, A, B, C))
        return 0;

    if (roots[0] > roots[1])
        swap12(&roots[0], &roots[1]);

    return 1;
}

void
square2disk2(float *phi, float *r, double seedx, double seedy)
{

   double a = 2.*seedx - 1;   /* (a,b) is now on [-1,1]^2 */
   double b = 2.*seedy - 1;

   if (a > -b) {     /* region 1 or 2 */
       if (a > b) {  /* region 1, also |a| > |b| */
           *r = a;
           *phi = (M_PI/4.) * (b/a);
       }
       else       {  /* region 2, also |b| > |a| */
           *r = b;
           *phi = (M_PI/4.) * (2. - (a/b));
       }
   }
   else {        /* region 3 or 4 */
       if (a < b) {  /* region 3, also |a| >= |b|, a != 0 */
            *r = -a;
            *phi = (M_PI/4.) * (4. + (b/a));
       }
       else       {  /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
            *r = -b;
            if (b != 0.)
                *phi = (M_PI/4.) * (6. - (a/b));
            else
                *phi = 0.;
       }
   }
   *r *= 0.9999999999999;	/* prophylactic against MS sin()/cos() impl. */
}

static
int get_tau_rct(float x, float y, float* r, float* muS) {
  float xr = x / TRANSMITTANCE_W;
  float yr = y / TRANSMITTANCE_H;

  double *ds;
  // map xr yr to unit disk centered around origin
  square2disk2(muS, r, xr, yr);
  // *r = EARTHRADIUS + yr * yr * ATMOSHEIGHT;
  // *muS = -0.15 + tan(1.5 * xr) / tan(1.5) * (1.0 + 0.15);
  return 1;
}


float
compute_optical_depth(const float scale_height, const float radius, const float ct)
{
    double seg = 0;
    float result = 0;
    double roots[2] = {0};
    unsigned int i;

    if (!hit_sphere(radius, ct, ATMOSRADIUS, roots) || roots[1] < 0)
        return 0;
    const float tmax = roots[1];
    const float apprx_seglen = 200;
    const unsigned int nsamp = (tmax / apprx_seglen) ? tmax / apprx_seglen : 2;
    const float seglen = tmax / nsamp;
    double optical_depth = exp(-(radius - EARTHRADIUS) / scale_height);
    for (i = 0; i < nsamp; ++i) {
        seg += seglen;
        const float _r = sqrt(radius * radius + seg * seg + 2 * radius * seg * ct);
        const float hr = exp(-(_r - EARTHRADIUS) / scale_height);
        result += (optical_depth + hr) * seglen * 0.5;
        optical_depth = hr;
    }
    return result;
}

static
void
set_transmittance_interpolator(FILE *fp, INTERP2* trans_ip, float* transmittance)
{
    float r;
    float ct;
    for (unsigned int j = 0; j < TRANSMITTANCE_H; ++j) {
        for (unsigned int i = 0; i < TRANSMITTANCE_W; ++i) {
            float x = (float)i / TRANSMITTANCE_W;
            float y = (float)j / TRANSMITTANCE_H;
            trans_ip->spt[i+TRANSMITTANCE_W*j][0] = x;
            trans_ip->spt[i+TRANSMITTANCE_W*j][1] = y;
            get_tau_rct(x+0.5, y+0.5, &r, &ct);
            transmittance[i+TRANSMITTANCE_W*j] = compute_optical_depth(HR_MS, r, ct);
            fprintf(fp, "%f %f %f\n", x, y, transmittance[i+j*TRANSMITTANCE_W]);
        }
    }
}

static
void
set_inscatter_interpolator(INTERP2* ins_ip, double* inscatter)
{
    for (unsigned int j = 0; j < TRANSMITTANCE_H; ++j) {
        for (unsigned int i = 0; i < TRANSMITTANCE_W; ++i) {
            double x = (double)i / TRANSMITTANCE_W;
            double y = (double)j / TRANSMITTANCE_H;
            ins_ip->spt[i+j*i][0] = x;
            ins_ip->spt[i+j*i][1] = y;
            inscatter[i+j*i] = compute_optical_depth(HR_MS, ATMOSRADIUS, y);
        }
    }
}

static
int
get_tau_xy(double *x, double *y, const double r, const double ct) {
	// *x = sqrt((r - EARTHRADIUS) / (ATMOSRADIUS - EARTHRADIUS));
	// *y = atan((ct + 0.15) / (1.0 + 0.15) * tan(1.5)) / (1.5);
    double r_norm = (r - EARTHRADIUS) / ATMOSHEIGHT;
    double sin_theta = sqrt(1 - ct * ct);

    double dx = r_norm * sin_theta;
    double dy = r_norm * ct;
    double sq[2];
    disk2square(sq, dx, dy);
    *x = sq[0];
    *y = sq[1];
    return 1;
}


static
float
interpolate_transmittance(INTERP2* ip, double r, double ct, float* transmittance)
{
    float wt[ip->ns];
    double x, y;
    int i;
    get_tau_xy(&x, &y, r, ct);
    if (!interp2_weights(wt, ip, x, y))
        exit(1);
    float result = 0;
    printf("tau300: %f\n", transmittance[16383]);
    for (i=0; i < ip->ns; i++) {
        // printf("wt: %f, transmittance: %f\n", wt[i], transmittance[i]);
        result += wt[i] * transmittance[i];
    }
    return result;
}


static
int
load_transmittance(FILE *fp, INTERP2* ip, float* transmittance)
{
    for (unsigned int j = 0; j < TRANSMITTANCE_H; ++j) {
        for (unsigned int i = 0; i < TRANSMITTANCE_W; ++i) {
            float x, y;
            fscanf(fp, "%f %f %f\n", &x, &y, &transmittance[i+j*TRANSMITTANCE_W]);
            ip->spt[i+j*TRANSMITTANCE_W][0] = x;
            ip->spt[i+j*TRANSMITTANCE_W][1] = y;
        }
    }
    return 1;
}


int
precompute(int sorder)
{
    const int nsamp = TRANSMITTANCE_H * TRANSMITTANCE_W;
    float transmittance[nsamp];
    INTERP2 *tau_ip = interp2_alloc(nsamp);
    INTERP2 *insr_ip = interp2_alloc(nsamp);
    INTERP2 *insm_ip = interp2_alloc(nsamp);

    interp2_free(tau_ip);
    char *tname = "tau.dat";
    FILE *fp;
    // Load transmittance data
    if ((fp = fopen(tname, "r")) != NULL) {
        load_transmittance(fp, tau_ip, transmittance);
    } else {
        fp = fopen(tname, "w");
        set_transmittance_interpolator(fp, tau_ip, transmittance);
        // save transmittance data to file
    };
    fclose(fp);
    double r = 6361e3;
    double ct = 0.5;
    float tau = interpolate_transmittance(tau_ip, r, ct, transmittance);
    printf("tau: %f\n", tau);
    return 1;
}


int
main(int argc, char *argv[])
{
    precompute(1);
    return 1;
}
