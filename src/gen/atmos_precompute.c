#include "atmosphere.h"
#include "data.h"
#include "interp2d.h"
#include "fvect.h"
#include <stdio.h>

char *progname;

const double EARTHRADIUS = 6360e3;
const double ATMOSRADIUS = 6420e3;
const double ATMOSHEIGHT = ATMOSRADIUS - EARTHRADIUS;
const double ATMOSRADIUS_M = 6395e3;
const unsigned int TRANSMITTANCE_W = 256;
const unsigned int TRANSMITTANCE_H = 64;
const float START_WVL = 380.0;
const float END_WVL = 780.0;
const unsigned int NSAMP = 20;

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

const int TRANSMITTANCE_INTEGRAL_SAMPLES = 500;
const int IRRADIANCE_INTEGRAL_SAMPLES = 32;
const int INSCATTER_INTEGRAL_SAMPLES = 50;
const int INSCATTER_SPHERICAL_INTEGRAL_SAMPLES = 16;

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


const unsigned int RES_R = 32;
const unsigned int RES_MU = 128;
const unsigned int RES_MU_S = 32;
const unsigned int RES_NU = 8;

static void swap12(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}


static double linear_interp(double x, const double *xp, const double *fp, int size)
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


static int hit_sphere(double *roots, const double r, const double ct, const double radius)
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

 
double limit(double r, double mu) {
  double d_out = -r * mu + sqrt(r * r * (mu * mu - 1.0) + (ATMOSRADIUS+1e3) * (ATMOSRADIUS+1e3));
  double delta_sq = r * r * (mu * mu - 1.0) + EARTHRADIUS * EARTHRADIUS;
  if (delta_sq >= 0.0) {
    double d_in = -r * mu - sqrt(delta_sq);
    if (d_in >= 0.0) {
      d_out = d_out < d_in ? d_out : d_in;
    }
  }
  return d_out;
}


static double compute_optical_depth(const double scale_height, const double radius, const double ctheta)
{
    double seg = 0;
    double result = 0;
    double roots[2] = {0, 0};
    unsigned int i;

    // if (!hit_sphere(roots, radius, ctheta, ATMOSRADIUS) || roots[1] < 0)
    //     return 0;
    // const double tmax = roots[1];
    const double tmax = limit(radius, ctheta);
    FILE *fp = fopen("tmax.txt", "a");
    fprintf(fp, " %f\n", tmax);
    fclose(fp);
    const double apprx_seglen = 200.;
    const unsigned int nsamp = (tmax / apprx_seglen) ? tmax / apprx_seglen : 2;
    const double seglen = tmax / nsamp;
    double optical_depth = exp(-(radius - EARTHRADIUS) / scale_height);
    for (i = 0; i < nsamp; ++i) {
        seg += seglen;
        const double _r = sqrt(radius * radius + seg * seg + 2 * radius * seg * ctheta);
        const double hr = exp(-(_r - EARTHRADIUS) / scale_height);
        result += (optical_depth + hr) * seglen * 0.5;
        optical_depth = hr;
    }
    double ct_horizon = -sqrt(1.0 - (EARTHRADIUS / radius) * (EARTHRADIUS / radius));
    return ctheta < ct_horizon ? 1e9 : result;
}


static int compute_transmittance(double *sumr, const double r, const double ct)
{
    const double taur = compute_optical_depth(HR_MS, r, ct);
    const double taum = compute_optical_depth(HMS_CC, r, ct);
    FILE *fp = fopen("tau.txt", "a");
    fprintf(fp, " %f %f\n", taur, taum);
    fclose(fp);
    for (int i = 0; i < NSAMP; ++i) {
        sumr[i] = exp(-(taur * BR0_MS[i] + taum * BM0_CC[i]));
    }
    return 1;
}


static
void
write_transmittance_data(char *fname)
{
    double radi;
    double ctheta;
    double tau[NSAMP];
    FILE *fp = fopen(fname, "w");

    /* .dat file header */
    fprintf(fp, "3\n");
    fprintf(fp, "0 %d %d\n", TRANSMITTANCE_W, TRANSMITTANCE_W);
    fprintf(fp, "0 %d %d\n", TRANSMITTANCE_H, TRANSMITTANCE_H);
    fprintf(fp, "0 %d %d\n", NSAMP, NSAMP);

    for (unsigned int j = 0; j < TRANSMITTANCE_H; ++j) {
        for (unsigned int i = 0; i < TRANSMITTANCE_W; ++i) {

            /* x, y to radius and cos(theta) */
            float xr = (i + 0.5) / TRANSMITTANCE_W;
            float yr = (j + 0.5) / TRANSMITTANCE_H;
            radi = EARTHRADIUS + yr * yr * ATMOSHEIGHT;
            ctheta = -0.15 + tan(1.5 * xr) / tan(1.5) * 1.15;

            // FILE *fp2 = fopen("tau.txt", "a");
            // fprintf(fp2, "%f %f %f %f", i+0.5, j+0.5, radi, ctheta);
            // fclose(fp2);

            compute_transmittance(tau, radi, ctheta);
            for (int k = 0; k < NSAMP; ++k) {
                fprintf(fp, "%f\n", tau[k]);
            }
        }
    }
    fclose(fp);
}


static
int
interpolate_transmittance(double *result, DATARRAY* dp, double radius, double ctheta)
{
    double x, y;
    double pt[3];
    int i;
    pt[1] = sqrt((radius - EARTHRADIUS) / ATMOSHEIGHT) * TRANSMITTANCE_H;
    pt[0] = (atan((ctheta + 0.15) / 1.15 * tan(1.5)) / 1.5) * TRANSMITTANCE_W;
    for (i = 0; i < NSAMP; ++i) {
        pt[2] = i;
        printf("pt: %f %f %f\n", pt[0], pt[1], pt[2]);
        result[i] = datavalue(dp, pt);
    }
    return 1;
}


static void set_layer(int layer, double* r, double* dmin, double* dmax, double* dminp, double* dmaxp) {
    double u = (layer * layer) / ((RES_R - 1.0) * (RES_R - 1.0));
    *r = sqrt(EARTHRADIUS * EARTHRADIUS + u * (ATMOSRADIUS * ATMOSRADIUS - EARTHRADIUS * EARTHRADIUS)) + (layer == 0 ? 0.01 : (layer == RES_R - 1 ? -0.001 : 0.0));
    *dmin = ATMOSRADIUS - *r;
    *dmax = sqrt((*r) * (*r) - EARTHRADIUS * EARTHRADIUS) + sqrt(ATMOSRADIUS * ATMOSRADIUS - EARTHRADIUS * EARTHRADIUS);
    *dminp = *r - EARTHRADIUS;
    *dmaxp = sqrt((*r) * (*r) - EARTHRADIUS * EARTHRADIUS);
}


void compute_inscatter1_integrand(DATARRAY *tau_dp, double radius, double ct_view, double ct_sun, double ct_view_sun, double xj, double *rayleigh, double *mie) {
    double iradius = sqrt(radius * radius + xj * xj + 2 * radius * xj * ct_view);
    double isun_theta = (ct_view_sun * xj + ct_sun * radius) / iradius;
    iradius = iradius > EARTHRADIUS ? iradius : EARTHRADIUS;
    if (isun_theta >= -sqrt(1.0 - EARTHRADIUS * EARTHRADIUS / (iradius * iradius))) {
        double tau[NSAMP];
        double taui[NSAMP];
        interpolate_transmittance(tau, tau_dp, radius, ct_view);
        interpolate_transmittance(taui, tau_dp, iradius, isun_theta);
        for (int i = 0; i < NSAMP; ++i) {
            double ti = tau[i] * taui[i];
            rayleigh[i] = ti * exp(-(iradius - EARTHRADIUS) / HR_MS);
            mie[i] = ti * exp(-(iradius - EARTHRADIUS) / HR_MS);
        }
    } else {
        for (int i = 0; i < NSAMP; ++i) {
            rayleigh[i] = 0;
            mie[i] = 0;
        }
    }
}


void ComputeInscatter1(
    const TransmittanceTexture& transmittance_sampler, Length r, Length dmin,
    Length dmax, Length dminp, Length dmaxp, vec2 xy,
    IrradianceSpectrum* rayleigh, IrradianceSpectrum* mie) {
  Number mu, muS, nu;
  Number x = xy.x - 0.5;
  Number y = xy.y - 0.5;
  if (y < RES_MU / 2.0) {
    Length d = (1.0 - y / (RES_MU / 2.0 - 1.0)) * dmaxp;
    d = min(max(dminp, d), dmaxp * 0.999);
    *mu = (Rg * Rg - r * r - d * d) / (2.0 * r * d);
    *mu = min(*mu, -sqrt(1.0 - (Rg / r) * (Rg / r)) - 0.001);
  } else {
    Length d = ((y - RES_MU / 2.0) / (RES_MU / 2.0 - 1.0)) * dmax;
    d = min(max(dmin, d), dmax * 0.999);
    *mu = (Rt * Rt - r * r - d * d) / (2.0 * r * d);
  }
  *muS = mod(x, RES_MU_S) / (RES_MU_S - 1.0);
  // paper formula
  // muS = -(0.6 + log(1.0 - muS * (1.0 -  exp(-3.6)))) / 3.0;
  // better formula
  *muS = tan((2.0 * *muS - 1.0 + 0.26) * 1.1 * rad) / tan(1.26 * 1.1 * rad);
  *nu = -1.0 + floor(x / RES_MU_S) / (RES_NU - 1) * 2.0;

  DimensionlessSpectrum ray_sum = DimensionlessSpectrum(0.0);
  DimensionlessSpectrum mie_sum = DimensionlessSpectrum(0.0);
  Length dx = Limit(r, mu) / INSCATTER_INTEGRAL_SAMPLES;
  DimensionlessSpectrum rayi;
  DimensionlessSpectrum miei;
  ComputeInscatter1Integrand(
      transmittance_sampler, r, mu, muS, nu, 0.0 * m, &rayi, &miei);
  for (int i = 1; i <= INSCATTER_INTEGRAL_SAMPLES; ++i) {
    Length xj = i * dx;
    DimensionlessSpectrum rayj;
    DimensionlessSpectrum miej;
    ComputeInscatter1Integrand(
        transmittance_sampler, r, mu, muS, nu, xj, &rayj, &miej);
    ray_sum += (rayi + rayj);
    mie_sum += (miei + miej);
    rayi = rayj;
    miei = miej;
  }
  *rayleigh = ray_sum * (dx / 2) * RayleighScattering() * SolarSpectrum();
  *mie = mie_sum * (dx / 2) * MieScattering() * SolarSpectrum();
}

void compute_inscatter1(DATARRAY *tau_dp, double radius, double dmin, double dmax, double dminp, double dmaxp, double i, double j, double *rayleigh, double *mie) {
    double tau[NSAMP];

    double i = i - 0.5;
    double j = j - 0.5;
    if (j < RES_MU / 2.0) {
        double d = (1.0 - j / (RES_MU / 2.0 - 1.0)) * dmaxp;
        d = min(max(dminp, d), dmaxp * 0.999);
        double ct_view = (EARTHRADIUS * EARTHRADIUS - radius * radius - d * d) / (2.0 * radius * d);
        ct_view = min(ct_view, -sqrt(1.0 - (EARTHRADIUS / radius) * (EARTHRADIUS / radius)) - 0.001);
    } else {
        double d = ((j - RES_MU / 2.0) / (RES_MU / 2.0 - 1.0)) * dmax;
        d = min(max(dmin, d), dmax * 0.999);
        double ct_view = (ATMOSRADIUS * ATMOSRADIUS - radius * radius - d * d) / (2.0 * radius * d);
    }
    ct_sun = mod(i, RES_MU_S) / (RES_MU_S - 1.0);
    ct_sun = tan((2.0 * muS - 1.0 + 0.26) * 1.1 ) / tan(1.26 * 1.1);
    nu = -1.0 + floor(i / RES_MU_S) / (RES_NU - 1) * 2.0;

    double ray_sum = {0.0};
    double mie_sum = {0.0};
    double dx = limit(radius, ct_view) / INSCATTER_INTEGRAL_SAMPLES;
    double *rayi;
    double *miei;
    compute_inscatter1_integrand(tau_dp, radius, ct_view, ct_sun, ct_view_sun, 0, rayi, miei)
    for (int i=0; i<=INSCATTER_INTEGRAL_SAMPLES; ++i) {
        double xj = i * dx;
        double *rayj;
        double *miej;
        compute_inscatter1_integrand(tau_dp, radius, ct_view, ct_sun, ct_view_sun, xj, rayj, miej);
        for (int j=0; j<NSAMP; ++j) {
            ray_sum[j] = rayi[j] + rayj[j];
            mie_sum[j] = miei[j] + miej[j];
            rayi[j] = rayj[j];
            miei[j] = miej[j];
        }
    }
    for (int k=0; k<NSAMP; ++k) {
        rayleigh[k] = ray_sum[k] * (dx / 2) * BR0_MS[k] * EXT_SOL[k];
        mie[k] = mie_sum[k] * (dx / 2) * BM0_CC[k] * EXT_SOL[k];
    }
}


void write_inscatter1(const DATARRAY *tau_dp, char *rname, char *mname) {
    FILE *rfp = fopen(rname, "w");
    FILE *mfp = fopen(mname, "w");
    fprintf(rfp, "5\n");
    fprintf(rfp, "0 %d %d\n", RES_R, RES_R);
    fprintf(rfp, "0 %d %d\n", RES_MU, RES_MU);
    fprintf(rfp, "0 %d %d\n", RES_MU_S, RES_MU_S);
    fprintf(rfp, "0 %d %d\n", RES_NU, RES_NU);
    fprintf(rfp, "0 %d %d\n", NSAMP, NSAMP);
    fprintf(mfp, "5\n");
    fprintf(mfp, "0 %d %d\n", RES_R, RES_R);
    fprintf(mfp, "0 %d %d\n", RES_MU, RES_MU);
    fprintf(mfp, "0 %d %d\n", RES_MU_S, RES_MU_S);
    fprintf(mfp, "0 %d %d\n", RES_NU, RES_NU);
    fprintf(mfp, "0 %d %d\n", NSAMP, NSAMP);
    for (unsigned int k=0; k< RES_R; ++k) {
        double radius, dmin, dmax, dminp, dmaxp;
        set_layer(k, &radius, &dmin, &dmax, &dminp, &dmaxp);
        for (unsigned int j=0; j< RES_MU; ++j) {
            for (unsigned int i=0; i < RES_MU_S; ++i) {
                for (unsigned int l=0; l < RES_NU; ++l) {
                    double rayleigh[NSAMP];
                    double mie[NSAMP];
                    compute_inscatter1(tau_dp, radius, dmin, dmax, dminp, dmaxp, i+0.5, j+0.5, rayleigh, mie);
                    for (int m = 0; m < NSAMP; ++m) {
                        fprintf(rfp, "%f\n", rayleigh[m]);
                        fprintf(mfp, "%f\n", mie[m]);
                    }
                }
            }
        }
    }
}


int
precompute(int sorder)
{
    const int nsamp = TRANSMITTANCE_H * TRANSMITTANCE_W;
    float transmittance[nsamp];
    DATARRAY *tau_dp;
    DATARRAY *inscat1r_dp;
    DATARRAY *inscat1m_dp;


    printf("setting interpolators...\n");
    char *tname = "tau.dat";
    if (getpath(tname, getrlibpath(), R_OK) == NULL) {
        printf("writing transmittance data to file...\n");
        write_transmittance_data(tname);
    };
    tau_dp = getdata(tname);

    // double r = 6360e3;
    // double ct = 0.5;
    // double tau[NSAMP];
    // interpolate_transmittance(tau, tau_dp, r, ct);
    // for (int i = 0; i < NSAMP; ++i) {
    //     printf("tau: %f\n", tau[i]);
    // }
    //
    char *i1rname = "inscat1r.dat";
    char *i1mname = "inscat1m.dat";
    if (getpath(i1rname, getrlibpath(), R_OK) == NULL) {
        printf("writing transmittance data to file...\n");
        write_inscat1_data(i1rname);
    };
    tau_dp = getdata(tname);
    
    return 1;
}


int
main(int argc, char *argv[])
{

    progname = "atmos_precompute";
    printf("precomputing...\n");
    precompute(1);
    return 1;
}
