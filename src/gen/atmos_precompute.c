#include "atmosphere.h"
#include "data.h"
#include "interp2d.h"
#include "rtmath.h"
#include <math.h>
#include <stdio.h>

char *progname;

const double ER = 6360e3;       /* Eearth radius in meters */
const double AR = 6420e3;       /* Atmosphere radius in meters */
const double AH = AR - ER;      /* Atmosphere thickness in meters */
const double HMAX = 875671;     /* Maximum atmosphere thickness in meters: sqrt(AR^2 - ER^2) */

// const double AR_M = 6395e3;

const float START_WVL = 380.0;
const float END_WVL = 780.0;
const unsigned int NSSAMP = 20;
const double SUNRAD = 0.004625;

// The cosine of the maximum Sun zenith angle for which atmospheric scattering
// must be precomputed (for maximum precision, use the smallest Sun zenith
// angle yielding negligible sky light radiance values. For instance, for the
// Earth case, 102 degrees is a good choice - yielding mu_s_min = -0.2).
const double MU_S_MIN = -0.2;

const double GROUND_ALBEDO = 0.2;

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


const double MIE_G = 0.76;

const unsigned int RES_R = 32;
const unsigned int RES_MU = 128;
const unsigned int RES_MU_S = 32;
const unsigned int RES_NU = 8;

const unsigned int TRANSMITTANCE_W = 256;
const unsigned int TRANSMITTANCE_H = 64;

const unsigned int IRRADIANCE_W = 64;
const unsigned int IRRADIANCE_H = 16;



typedef struct {
    double width;
    double exp_term;
    double exp_scale;
    double linear_term;
    double constant_term;
} DensityProfileLayer;


typedef struct {
    DensityProfileLayer layers[2];
} DensityProfile;

DensityProfile rayleigh_density = {
    .layers = {
        {
            .width = 0.0,
            .exp_term = 1.0,
            .exp_scale = -1 / HR_MS,
            .linear_term = 0.0,
            .constant_term = 0.0
        },
    }
};

DensityProfile mie_density = {
    .layers = {
        {
            .width = 0.0,
            .exp_term = 1.0,
            .exp_scale = -1 / HMS_CC,
            .linear_term = 0.0,
            .constant_term = 0.0
        },
    }
};


static inline double dmin(double x, double y) {
    return x < y ? x : y;
}

static inline double dmax(double x, double y) {
    return x > y ? x : y;
}


double
smoothstep(double edge0, double edge1, double x)
{
    x = (x - edge0) / (edge1 - edge0);
    x = x < 0.0 ? 0.0 : x > 1.0 ? 1.0 : x;
    return x * x * (3.0 - 2.0 * x);
}


double
limit(double r, double mu)
{
    double d_out = -r * mu + sqrt(r * r * (mu * mu - 1.0) + (AR+1e3) * (AR+1e3));
    double delta_sq = r * r * (mu * mu - 1.0) + ER * ER;
    if (delta_sq >= 0.0) {
        double d_in = -r * mu - sqrt(delta_sq);

        if (d_in >= 0.0) {
            d_out = d_out < d_in ? d_out : d_in;
        }
    }
    return d_out;
}


double
distance_to_space(double r, double mu)
{
    if (r > AR) {
        return 0;
    }
    if (mu < -1.0 || mu > 1.0) {
        return 0;
    }
    double descriminant = r * r * (mu * mu - 1.0) + AR * AR;
    double distance = -r * mu + sqrt(descriminant);
    return distance < 0.0 ? 0.0 : distance;
}


double
distance_to_earth(double r, double mu)
{
    if (r < ER) {
        return 0;
    }
    if (mu < -1.0 || mu > 1.0) {
        return 0;
    }
    double descriminant = r * r * (mu * mu - 1.0) + ER * ER;
    double distance = -r * mu - sqrt(descriminant);
    return distance < 0.0 ? 0.0 : distance;
}


static
int
ray_intersects_ground(double r, double mu)
{
    return mu < 0.0 && r * r * (mu * mu - 1.0) + ER * ER >= 0.0;
}


double
get_layer_density(DensityProfileLayer *layer, double altitude)
{
    double density = layer->exp_term * exp(layer->exp_scale * altitude) + layer->linear_term * altitude + layer->constant_term;
    return density < 0.0 ? 0 : density > 1.0 ? 1 : density;
}


double
get_profile_density(DensityProfile *profile, double altitude)
{
    return altitude < profile->layers[0].width ? get_layer_density(&(profile->layers[0]), altitude) : get_layer_density(&(profile->layers[1]), altitude);
}


static
double
compute_optical_depth(const double scale_height, const double radius, const double ctheta)
{
    double seg = 0;
    double result = 0;
    double roots[2] = {0, 0};
    unsigned int i;

    // if (!hit_sphere(roots, radius, ctheta, AR) || roots[1] < 0)
    //     return 0;
    // const double tmax = roots[1];
    const double tmax = limit(radius, ctheta);
    // FILE *fp = fopen("tmax.txt", "a");
    // fprintf(fp, " %f\n", tmax);
    // fclose(fp);
    const double apprx_seglen = 200.;
    const unsigned int nsamp = (tmax / apprx_seglen) ? tmax / apprx_seglen : 2;
    const double seglen = tmax / nsamp;
    double optical_depth = exp(-(radius - ER) / scale_height);
    for (i = 0; i < nsamp; ++i) {
        seg += seglen;
        const double _r = sqrt(radius * radius + seg * seg + 2 * radius * seg * ctheta);
        const double hr = exp(-(_r - ER) / scale_height);
        result += (optical_depth + hr) * seglen * 0.5;
        optical_depth = hr;
    }
    double ct_horizon = -sqrt(1.0 - (ER / radius) * (ER / radius));
    return ctheta < ct_horizon ? 1e9 : result;
}


static
int
compute_transmittance(const double r, const double ct, double *sumr)
{
    const double taur = compute_optical_depth(HR_MS, r, ct);
    const double taum = compute_optical_depth(HMS_CC, r, ct);
    // FILE *fp = fopen("tau.txt", "a");
    // fprintf(fp, " %f %f\n", taur, taum);
    // fclose(fp);
    for (int i = 0; i < NSSAMP; ++i) {
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
    double tau[NSSAMP];
    FILE *fp = fopen(fname, "w");

    /* .dat file header */
    fprintf(fp, "3\n");
    fprintf(fp, "0 %d %d\n", TRANSMITTANCE_W, TRANSMITTANCE_W);
    fprintf(fp, "0 %d %d\n", TRANSMITTANCE_H, TRANSMITTANCE_H);
    fprintf(fp, "0 %d %d\n", NSSAMP, NSSAMP);

    for (unsigned int j = 0; j < TRANSMITTANCE_H; ++j) {
        for (unsigned int i = 0; i < TRANSMITTANCE_W; ++i) {

            /* x, y to radius and cos(theta) */
            float xr = (i + 0.5) / TRANSMITTANCE_W;
            float yr = (j + 0.5) / TRANSMITTANCE_H;
            radi = ER + yr * yr * AH;
            ctheta = -0.15 + tan(1.5 * xr) / tan(1.5) * 1.15;

            // FILE *fp2 = fopen("tau.txt", "a");
            // fprintf(fp2, "%f %f %f %f", i+0.5, j+0.5, radi, ctheta);
            // fclose(fp2);

            compute_transmittance(radi, ctheta, tau);
            for (int k = 0; k < NSSAMP; ++k) {
                fprintf(fp, "%f\n", tau[k]);
            }
        }
    }
    fclose(fp);
}


double
get_texture_coord_from_unit_range(double x, int texture_size)
{
    return 0.5 / texture_size + x * (1.0 - 1.0 / texture_size);
}


double
get_unit_range_from_texture_coord(double u, int texture_size)
{
    return (u - 0.5 / texture_size) / (1.0 - 1.0 / texture_size);
}


static
int
interpolate_transmittance(DATARRAY* dp, double radius, double ctheta, double *result)
{
    double x, y;
    double pt[3];
    int i;
    pt[1] = sqrt((radius - ER) / AH) * TRANSMITTANCE_H;
    pt[0] = (atan((ctheta + 0.15) / 1.15 * tan(1.5)) / 1.5) * TRANSMITTANCE_W;
    for (i = 0; i < NSSAMP; ++i) {
        pt[2] = i;
        // printf("pt: %f %f %f\n", pt[0], pt[1], pt[2]);
        result[i] = datavalue(dp, pt);
    }
    return 1;
}


static
void
get_transmittance(DATARRAY *tau_dp, double r, double mu, double d, int ray_r_mu_intersects_ground, double *result)
{
    double r_d = sqrt(d * d + 2.0 * r * mu * d + r * r);
    r_d = r_d > AR ? AR : r_d < ER ? ER : r_d;
    double mu_d = (r * mu + d) / r_d;
    mu_d = mu_d > 1.0 ? 1.0 : mu_d < -1.0 ? -1.0 : mu_d;
    double result1[NSSAMP];
    double result2[NSSAMP];
    double v;

    if (ray_r_mu_intersects_ground) {
        interpolate_transmittance(tau_dp, r_d, -mu_d, result1);
        interpolate_transmittance(tau_dp, r, -mu, result2);
    } else {
        interpolate_transmittance(tau_dp, r, mu, result1);
        interpolate_transmittance(tau_dp, r_d, mu_d, result2);
    }
    for (int i = 0; i < NSSAMP; ++i) {
        v = result1[i] / result2[i];
        result[i] = v > 1 ? 1.0 : v;
    }
}


static
void
get_transmittance_to_sun(DATARRAY *tau_dp, double r, double mu_s, double *result)
{
    double sin_theta_h = ER / r;
    double cos_theta_h = (1.0 - sin_theta_h * sin_theta_h);
    cos_theta_h = cos_theta_h < 0.0 ? 0.0 : -sqrt(cos_theta_h);
    double v[NSSAMP];
    interpolate_transmittance(tau_dp, r, mu_s, v);
    double st = smoothstep(-sin_theta_h * SUNRAD, sin_theta_h * SUNRAD, mu_s - cos_theta_h);
    for (int i = 0; i < NSSAMP; ++i) {
        result[i] = v[i] * st;
    }
}


void
compute_single_scattering_integrand(
        DATARRAY *tau_dp,
        double r,
        double mu,
        double mu_s,
        double nu,
        double d,
        int ray_r_mu_intersects_ground,
        double *rayleigh,
        double *mie)
{
    double r_d = sqrt(d * d + 2.0 * r * mu * d + r * r);
    r_d =  r_d > AR ? AR : r_d < ER ? ER : r_d;
    double mu_s_d = (r * mu_s + d * nu) / r_d;
    mu_s_d = mu_s_d > 1.0 ? 1.0 : mu_s_d < -1.0 ? -1.0 : mu_s_d;
    // double transmittance[NSSAMP];
    double tau_r_mu[NSSAMP];
    double tau_sun[NSSAMP];
    get_transmittance(tau_dp, r, mu, d, ray_r_mu_intersects_ground, tau_r_mu);
    get_transmittance_to_sun(tau_dp, r_d, mu_s_d, tau_sun);
    double rayleigh_profile_density = get_profile_density(&rayleigh_density, r_d - ER);
    double mie_profile_density = get_profile_density(&mie_density, r_d - ER);
    double transmittance;
    for (int i=0; i< NSSAMP; ++i) {
        transmittance = tau_r_mu[i] * tau_sun[i];
        rayleigh[i] = transmittance * rayleigh_profile_density;
        mie[i] = transmittance * mie_profile_density;
    }
}


double
distance_to_nearst_atmosphere_boundary(double r, double mu, int ray_r_mu_intersects_ground)
{
    if (ray_r_mu_intersects_ground) {
        return distance_to_earth(r, mu);
    } else {
        return distance_to_space(r, mu);
    }
}


void
compute_single_scattering(DATARRAY *tau_dp, double r, double mu, double mu_s, double nu, int ray_r_mu_intersects_ground, double *rayleigh, double *mie)
{
    const int SAMPLE_COUNT = 50;
    double dx = distance_to_nearst_atmosphere_boundary(r, mu, ray_r_mu_intersects_ground) / SAMPLE_COUNT;
    double rayleigh_sum[NSSAMP];
    double mie_sum[NSSAMP];
    for (int i=0; i <= SAMPLE_COUNT; ++i) {
        double d_i = i * dx;
        double rayleigh_i[NSSAMP];
        double mie_i[NSSAMP];
        compute_single_scattering_integrand(tau_dp, r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
        double weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        for (int j=0; j<NSSAMP; ++j) {
            rayleigh_sum[j] += rayleigh_i[j] * weight_i;
            mie_sum[j] += mie_i[j] * weight_i;
        }
    }
    for (int k=0; k<NSSAMP; ++k) {
        rayleigh[k] = rayleigh_sum[k] * dx * EXTSOL[k] * BR0_MS[k];
        mie[k] = mie_sum[k] * dx * EXTSOL[k] * BM0_CC[k];
    }
}


static
void
compute_direct_irradiance(DATARRAY *tau_dp, double r, double mu_s, double *result)
{
    // Approximate the average of the cosine factor mu_s over the visible fraction of the sun disc
    double average_cosine_factor = mu_s < -SUNRAD ? 0.0 : (mu_s > SUNRAD ? mu_s : ((mu_s + SUNRAD) * (mu_s + SUNRAD)) / (4.0 * SUNRAD));
    double transmittance[NSSAMP];
    interpolate_transmittance(tau_dp, r, mu_s, transmittance);
    for (int i = 0; i < NSSAMP; ++i) {
        result[i] = EXTSOL[i] * transmittance[i] * average_cosine_factor;
    }
}


static
double
rayleigh_phase_function(double nu)
{
    double k = 3.0 / (16.0 * PI);
    return k * (1.0 + nu * nu);
}


static
double
mie_phase_function(double g, double nu)
{
    double k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}


static
void
to_scattering_uvwz(double r, double mu, double mu_s, double nu, int ray_r_mu_intersects_ground, double *u, double *v, double *w, double *z)
{
    double rho = r * r - ER * ER;
    rho = rho < 0.0 ? 0.0 : sqrt(rho);
    double u_r = get_texture_coord_from_unit_range(rho / HMAX, RES_R);

    double r_mu = r * mu;
    double discriminant = r_mu * r_mu - r * r + ER * ER;

    double u_mu;
    if (ray_r_mu_intersects_ground) {
        double d = -r_mu - sqrt(discriminant);
        double d_min = r - ER;
        double d_max = rho;
        u_mu = 0.5 - 0.5 * get_texture_coord_from_unit_range(d_max == d_min ? 0.0 : (d - d_min) / (d_max - d_min), RES_MU / 2);
    } else {
        double d = -r_mu + sqrt(discriminant + HMAX * HMAX);
        double d_min = AR - r;
        double d_max = rho + HMAX;
        u_mu = 0.5 + 0.5 * get_texture_coord_from_unit_range((d - d_min) / (d_max - d_min), RES_MU / 2);
    }

    double d = distance_to_space(ER, mu_s);
    double d_min = AH;
    double d_max = HMAX;
    double a = (d - d_min) / (d_max - d_min);
    double D = distance_to_space(ER, MU_S_MIN);
    double A = (D - d_min) / (d_max - d_min);
    double u_mu_s = get_texture_coord_from_unit_range(dmax(1.0 - a / A, 0.0), RES_MU_S);
    double u_nu = (nu + 1.0) * 0.5;
    *u = u_nu;
    *v = u_mu_s;
    *w = u_mu;
    *z = u_r;
}


static
void
from_scattering_uvwz(double u, double v, double w, double z, double *r, double *mu, double *mu_s, double *nu, int *ray_r_mu_intersects_ground)
{
    double rho = HMAX * get_unit_range_from_texture_coord(w, RES_R);
    *r = sqrt(rho * rho + ER * ER);

    if (z < 0.5) {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
        // we can recover mu:
        double d_min = *r - ER;
        double d_max = rho;
        double d = d_min + (d_max - d_min) * get_unit_range_from_texture_coord(1.0 - 2.0 * z, RES_MU / 2);
        if (d == 0.0) {
            *mu = -1.0;
        } else {
            *mu = -(rho * rho + d * d) / (2.0 * *r * d);
            *mu = *mu < -1.0 ? -1.0 : *mu > 1.0 ? 1.0 : *mu;
            *ray_r_mu_intersects_ground = 1;
        }
    } else {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon) - from which we can recover mu:
        double d_min = AR - *r;
        double d_max = rho + HMAX;
        double d = d_min + (d_max - d_min) * get_unit_range_from_texture_coord(2.0 * z - 1.0, RES_MU / 2);
        if (d == 0.0) {
            *mu = 1.0;
        } else {
            *mu = -(HMAX * HMAX - rho * rho - d * d) / (2.0 * *r * d);
            *mu = *mu < -1.0 ? -1.0 : *mu > 1.0 ? 1.0 : *mu;
            *ray_r_mu_intersects_ground = 0;
        }
    }
    double x_mu_s = get_unit_range_from_texture_coord(v, RES_MU_S);
    double d_min = AH;
    double d_max = HMAX;
    double D = distance_to_space(ER, MU_S_MIN);
    double A = (D - d_min) / (d_max - d_min);
    double a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
    double d = d_min + dmin(a, A) * (d_max - d_min);
    if (d == 0.0) {
        *mu_s = 1.0;
    } else {
        *mu_s = -(HMAX * HMAX - d * d) / (2.0 * ER * d);
        *mu_s = *mu_s < -1.0 ? -1.0 : *mu_s > 1.0 ? 1.0 : *mu_s;
    }
    *nu = u * 2.0 - 1.0;
}


void
write_single_scattering_data(DATARRAY *tau_dp, char *rpath, char *mpath)
{
    FILE *rfp = fopen(rpath, "w");
    FILE *mfp = fopen(mpath, "w");
    fprintf(rfp, "5\n");
    fprintf(rfp, "0 %d %d\n", RES_R, RES_R);
    fprintf(rfp, "0 %d %d\n", RES_MU, RES_MU);
    fprintf(rfp, "0 %d %d\n", RES_MU_S, RES_MU_S);
    fprintf(rfp, "0 %d %d\n", RES_NU, RES_NU);
    fprintf(rfp, "0 %d %d\n", NSSAMP, NSSAMP);
    fprintf(mfp, "5\n");
    fprintf(mfp, "0 %d %d\n", RES_R, RES_R);
    fprintf(mfp, "0 %d %d\n", RES_MU, RES_MU);
    fprintf(mfp, "0 %d %d\n", RES_MU_S, RES_MU_S);
    fprintf(mfp, "0 %d %d\n", RES_NU, RES_NU);
    fprintf(mfp, "0 %d %d\n", NSSAMP, NSSAMP);
    for (unsigned int i=0; i< RES_R; ++i) {
        for (unsigned int j=0; j< RES_MU; ++j) {
            for (unsigned int k=0; k < RES_MU_S; ++k) {
                for (unsigned int l=0; l < RES_NU; ++l) {
                    double rayleigh[NSSAMP];
                    double mie[NSSAMP];
                    double r, mu, mu_s, nu;
                    int ray_r_mu_intersects_ground;
                    from_scattering_uvwz(i+0.5, j+0.5, k+0.5, l+0.5, &r, &mu, &mu_s, &nu, &ray_r_mu_intersects_ground);
                    compute_single_scattering(tau_dp, r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
                    for (int m = 0; m < NSSAMP; ++m) {
                        fprintf(rfp, "%f\n", rayleigh[m]);
                        fprintf(mfp, "%f\n", mie[m]);
                    }
                }
            }
        }
    }
}



static
void
interpolate_scattering(DATARRAY *dp, double r, double mu, double mu_s, double nu, int ray_r_mu_intersects_ground, double *result)
{
    double u, v, w, z;
    to_scattering_uvwz(r, mu, mu_s, nu, ray_r_mu_intersects_ground, &u, &v, &w, &z);
    for (int i = 0; i < NSSAMP; ++i) {
        double pt[5] = {u, v, w, z, i};
        result[i] = datavalue(dp, pt);
    }
}


static
void
to_irradiance_uv(double r, double mu_s, double *u, double *v)
{
    double x_r = (r - ER) / AH;
    double x_mu_s = mu_s * 0.5 + 0.5;
    *u = get_texture_coord_from_unit_range(x_mu_s, IRRADIANCE_W);
    *v = get_texture_coord_from_unit_range(x_r, IRRADIANCE_H);
}


static
void
from_irradiance_uv(double u, double v, double *r, double *mu_s)
{
    double x_mu_s = get_unit_range_from_texture_coord(u, IRRADIANCE_W);
    double x_r = get_unit_range_from_texture_coord(v, IRRADIANCE_H);
    *r = x_r * AH + ER;
    *mu_s = x_mu_s * 2.0 - 1.0;
    *mu_s = *mu_s < -1.0 ? -1.0 : *mu_s > 1.0 ? 1.0 : *mu_s;
}


static
void
get_irradiance(DATARRAY *dp, double r, double mu_s, double *result)
{
    double u, v;
    to_irradiance_uv(r, mu_s, &u, &v);
    double pt[2] = {u, v};
    for (int i = 0; i < NSSAMP; ++i) {
        double pt[3] = {u, v, i};
        result[i] = datavalue(dp, pt);
    }
}


static
void
get_scattering(
        DATARRAY *scat1r,
        DATARRAY *scat1m,
        DATARRAY *scat,
        double r,
        double mu,
        double mu_s,
        double nu,
        int ray_r_mu_intersects_ground,
        int scattering_order,
        double *result)
{
    if (scattering_order == 1) {
        double rayleigh[NSSAMP];
        double mie[NSSAMP];
        interpolate_scattering(scat1r, r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh);
        interpolate_scattering(scat1m, r, mu, mu_s, nu, ray_r_mu_intersects_ground, mie);
        double rayleigh_phase = rayleigh_phase_function(nu);
        double mie_phase = mie_phase_function(MIE_G, nu);
        for (int i = 0; i < NSSAMP; ++i) {
            result[i] = rayleigh[i] * rayleigh_phase + mie[i] * mie_phase;
        }
    } else {
        interpolate_scattering(scat, r, mu, mu_s, nu, ray_r_mu_intersects_ground, result);
    }
}


static
void
compute_multiple_scattering(
    DATARRAY *tau_dp,
    DATARRAY *scattering_density_dp,
    double r, double mu, double mu_s, double nu,
    int ray_r_mu_intersects_ground, double *result)
{

    int i,j;
    // Number of intervals for the numerical integration.
    const int SAMPLE_COUNT = 50;
    // The integration step, i.e. the length of each integration interval.
    double dx = distance_to_nearst_atmosphere_boundary(r, mu, ray_r_mu_intersects_ground) / SAMPLE_COUNT;
    // Integration loop.
    for (i = 0; i <= SAMPLE_COUNT; ++i) {
        double d_i = i * dx;

        // The r, mu and mu_s parameters at the current integration point (see the
        // single scattering section for a detailed explanation).
        double r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
        r_i = r_i > AR ? AR : r_i < ER ? ER : r_i;
        double mu_i = (r * mu + d_i) / r_i;
        mu_i = mu_i > 1.0 ? 1.0 : mu_i < -1.0 ? -1.0 : mu_i;
        double mu_s_i = (r * mu_s + d_i * nu) / r_i;
        mu_s_i = mu_s_i > 1.0 ? 1.0 : mu_s_i < -1.0 ? -1.0 : mu_s_i;

        // The Rayleigh and Mie multiple scattering at the current sample point.
        double rayleigh_mie_i[NSSAMP];
        double rayleigh_mie_s[NSSAMP];
        double transmittance[NSSAMP];
        interpolate_scattering(scattering_density_dp, r_i, mu_i, mu_s_i, nu, ray_r_mu_intersects_ground, rayleigh_mie_s);
        get_transmittance(tau_dp, r, mu, d_i, ray_r_mu_intersects_ground, transmittance);
        for (j = 0; j < NSSAMP; ++j)
            rayleigh_mie_i[j] = rayleigh_mie_s[j] * transmittance[j] * dx;
        // Sample weight (from the trapezoidal rule).
        double weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        for (j = 0; j < NSSAMP; ++j)
            result[j] += rayleigh_mie_i[j] * weight_i;
    }
}


static
void
compute_scattering_density(
        DATARRAY *tau_dp,
        DATARRAY *scat1r,
        DATARRAY *scat1m,
        DATARRAY *scat,
        DATARRAY *irrad_dp,
        double r,
        double mu,
        double mu_s,
        double nu,
        int scattering_order,
        double *result)
{
    // Compute unit direction vectors for the zenith, the view direction omega and
    // and the sun direction omega_s, such that the cosine of the view-zenith
    // angle is mu, the cosine of the sun-zenith angle is mu_s, and the cosine of
    // the view-sun angle is nu. The goal is to simplify computations below.
    double zenith_direction[3] = {0.0, 0.0, 1.0};
    double omega[3] = {sqrt(1.0 - mu * mu), 0.0, mu};
    double sun_dir_x = omega[0] == 0.0 ? 0.0 : (nu - mu * mu_s) / omega[0];
    double sun_dir_y = sqrt(dmax(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
    double omega_s[3] = {sun_dir_x, sun_dir_y, mu_s};
    const int SAMPLE_COUNT = 16;
    const double dphi = PI / SAMPLE_COUNT;
    const double dtheta = PI / SAMPLE_COUNT;

    // Nested loops for the integral over all the incident directions omega_i.
    for (int l = 0; l < SAMPLE_COUNT; ++l) {
        double theta = (l + 0.5) * dtheta;
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        int ray_r_theta_intersects_ground = ray_intersects_ground(r, cos_theta);

        // The distance and transmittance to the ground only depend on theta, so we
        // can compute them in the outer loop for efficiency.
        double distance_to_ground = 0.0;
        double transmittance_to_ground[NSSAMP];
        double ground_albedo;
        if (ray_r_theta_intersects_ground) {
            distance_to_ground = distance_to_earth(r, cos_theta);
            get_transmittance(tau_dp, r, cos_theta, distance_to_ground, 1, transmittance_to_ground);
            ground_albedo = GROUND_ALBEDO;
        }

        for (int m = 0; m < 2 * SAMPLE_COUNT; ++m) {
            double phi = (m + 0.5) * dphi;
            double omega_i[3] = {cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta};
            double domega_i = dtheta * dphi * sin(theta);

            // The radiance L_i arriving from direction omega_i after n-1 bounces is
            // the sum of a term given by the precomputed scattering texture for the
            // (n-1)-th order:
            double nu1 = fdot(omega_s, omega_i);
            double incident_radiance[NSSAMP];
            get_scattering(scat1r, scat1m, scat, r, omega_i[2], mu_s, nu1, ray_r_theta_intersects_ground, scattering_order - 1, incident_radiance);

            // and of the contribution from the light paths with n-1 bounces and whose
            // last bounce is on the ground. This contribution is the product of the
            // transmittance to the ground, the ground albedo, the ground BRDF, and
            // the irradiance received on the ground after n-2 bounces.
            double ground_normal[3] = {zenith_direction[0] * r + omega_i[0] * distance_to_ground, zenith_direction[1] * r + omega_i[1] * distance_to_ground, zenith_direction[2] * r + omega_i[2] * distance_to_ground};
            normalize(ground_normal);
            double ground_irradiance[NSSAMP];
            get_irradiance(irrad_dp, ER, fdot(ground_normal, omega_s), ground_irradiance);
            for (int i = 0; i < NSSAMP; ++i) {
                incident_radiance[i] += transmittance_to_ground[i] * ground_albedo * (1.0 / PI) * ground_irradiance[i];
            }

            // The radiance finally scattered from direction omega_i towards direction
            // -omega is the product of the incident radiance, the scattering
            // coefficient, and the phase function for directions omega and omega_i
            // (all this summed over all particle types, i.e. Rayleigh and Mie).
            double nu2 = fdot(omega, omega_i);
            double rayleigh_density_ = get_profile_density(&rayleigh_density, r - ER);
            double mie_density_ = get_profile_density(&mie_density, r - ER);
            for (int j = 0; j < NSSAMP; ++j) {
                result[j] += incident_radiance[j] * (BR0_MS[j] * rayleigh_density_ * rayleigh_phase_function(nu2) + mie_density_ * mie_phase_function(MIE_G, nu2) * BM0_CC[j]) * domega_i;
            }
        }
    }
}


static
void
write_direct_irradiance_data(DATARRAY *tau_dp, char *fname)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "3\n");
    fprintf(fp, "0 %d %d\n", IRRADIANCE_W, IRRADIANCE_W);
    fprintf(fp, "0 %d %d\n", IRRADIANCE_H, IRRADIANCE_H);
    fprintf(fp, "0 %d %d\n", NSSAMP, NSSAMP);
    for (unsigned int j = 0; j < IRRADIANCE_H; ++j) {
        for (unsigned int i = 0; i < IRRADIANCE_W; ++i) {
            double r, mu_s;
            double result[NSSAMP];
            from_irradiance_uv(i + 0.5, j + 0.5, &r, &mu_s);
            compute_direct_irradiance(tau_dp, r, mu_s, result);
            for (int k = 0; k < NSSAMP; ++k) {
                fprintf(fp, "%f\n", result[k]);
            }
        }
    }
}


static
void
write_initial_irradiance_data(char *fname)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "3\n");
    fprintf(fp, "0 %d %d\n", IRRADIANCE_W, IRRADIANCE_W);
    fprintf(fp, "0 %d %d\n", IRRADIANCE_H, IRRADIANCE_H);
    fprintf(fp, "0 %d %d\n", NSSAMP, NSSAMP);
    for (unsigned int j = 0; j < IRRADIANCE_H; ++j) {
        for (unsigned int i = 0; i < IRRADIANCE_W; ++i) {
            for (int k = 0; k < NSSAMP; ++k) {
                fprintf(fp, "%f\n", 0.0);
            }
        }
    }
}


static
void
compute_indirect_irradiance(DATARRAY *scat1r, DATARRAY *scat1m, DATARRAY *scat, double r, double mu_s, int scattering_order, double *result)
{
    const int SAMPLE_COUNT = 32;
    const double dphi = PI / SAMPLE_COUNT;
    const double dtheta = PI / SAMPLE_COUNT;

    double omega_s[3] = {sqrt(1.0 - mu_s * mu_s), 0.0, mu_s};
    for (int j = 0; j < SAMPLE_COUNT / 2; ++j) {
        double theta = (j + 0.5) * dtheta;
        for (int i = 0; i < 2 * SAMPLE_COUNT; ++i) {
            double phi = (i + 0.5) * dphi;
            double omega[3] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
            double domega = dtheta * dphi * sin(theta);
            double nu = fdot(omega, omega_s);
            double result1[NSSAMP];
            get_scattering(
                    scat1r,
                    scat1m,
                    scat,
                    r,
                    omega[2],
                    mu_s,
                    nu,
                    0, /* ray_r_theta_intersects_ground */
                    scattering_order,
                    result1);
            for (int k = 0; k < NSSAMP; ++k) {
                result[k] += result1[k] * omega[2] * domega;
            }
        }
    }
}


DATARRAY *
allocate_3d_datarray(char *name, int ri, int rj, int rk)
{
    DATARRAY *dp;
	if ((dp = (DATARRAY *)malloc(sizeof(DATARRAY))) == NULL)
		goto memerr;
    int asize = ri*rj*rk;
    dp->name = name;
    dp->type = DATATY;
    dp->nd = 3;
    dp->dim[0].org = 0;
    dp->dim[0].siz = ri;
    dp->dim[0].ne = ri;
    dp->dim[1].org = 0;
    dp->dim[1].siz = rj;
    dp->dim[1].ne = rj;
    dp->dim[2].org = 0;
    dp->dim[2].siz = rk;
    dp->dim[2].ne = rk;
    // dp->arr.p = malloc(sizeof(DATATYPE)*dp->dim[0].ne*dp->dim[1].ne*dp->dim[2].ne);
	if ((dp->arr.d = (DATATYPE *)malloc(asize*sizeof(DATATYPE))) == NULL)
		goto memerr;
    return(dp);

    memerr:
        fprintf(stderr, "Memory allocation error in allocate_3d_datarray\n");
        return(NULL);
}


DATARRAY *
allocate_5d_datarray(char *name, int ri, int rj, int rk, int rl, int rm)
{
    DATARRAY *dp;
    int asize = ri*rj*rk*rl*rm;
	if ((dp = (DATARRAY *)malloc(sizeof(DATARRAY))) == NULL)
		goto memerr;
    dp->name = name;
    dp->type = DATATY;
    dp->nd = 5;
    dp->dim[0].org = 0;
    dp->dim[0].siz = ri;
    dp->dim[0].ne = ri;
    dp->dim[1].org = 0;
    dp->dim[1].siz = rj;
    dp->dim[1].ne = rj;
    dp->dim[2].org = 0;
    dp->dim[2].siz = rk;
    dp->dim[2].ne = rk;
    dp->dim[3].org = 0;
    dp->dim[3].siz = rl;
    dp->dim[3].ne = rl;
    dp->dim[4].org = 0;
    dp->dim[4].siz = rm;
    dp->dim[4].ne = rm;
    if ((dp->arr.d = (DATATYPE *)malloc(asize*sizeof(DATATYPE))) == NULL)
        goto memerr;
    return(dp);
    memerr:
        fprintf(stderr, "Memory allocation error in allocate_5d_datarray\n");
        return(NULL);
}


/* set value by coordinate */
void
setvalue(DATARRAY *dp, int *pt, float val)
{
    int i;
    int idx;
    int ni = dp->nd - 1;
    for (i = 0; i < ni; i++) {
        idx += pt[i]*dp->dim[ni-i].ne;
    }
    idx += pt[i+1];
    dp->arr.d[idx] = val;
}


/* save to a .dat file */
void
savedata(DATARRAY *dp)
{
    if (dp == NULL)
        return;
    FILE *fp = fopen(dp->name, "w");
    int i, j;
    int nvals;
    for (i = 0; i < dp->nd; i++) {
        nvals *= dp->dim[i].ne;
    }
    fprintf(fp, "%d\n", dp->nd);
    for (i = 0; i < dp->nd; i++) {
        fprintf(fp, "%f %f %d\n", dp->dim[i].org, dp->dim[i].siz, dp->dim[i].ne);
    }
    for (i = 0; i < nvals; i++) {
        fprintf(fp, "%f\n", dp->arr.d[i]);
    }
    fclose(fp);
}


void
increment_dp(DATARRAY *dp1, DATARRAY *dp2)
{
    if (dp1 == NULL || dp2 == NULL)
        return;

    if (dp1->nd != dp2->nd)
        return;

    int i;
    int nvals;
    for (i = 0; i < dp1->nd; i++) {
        nvals *= dp1->dim[i].ne;
    }
    for (i = 0; i < nvals; i++) {
        dp1->arr.d[i] += dp2->arr.d[i];
    }
}


static
int
precompute(int sorder)
{
    unsigned int scattering_order;

    if (sorder < 2) {
        printf("scattering order must be at least 2\n");
        return 0;
    }

    const int nsamp = TRANSMITTANCE_H * TRANSMITTANCE_W;
    float transmittance[nsamp];
    unsigned int idx, i, j, k, l, m;


    DATARRAY *tau_dp = allocate_3d_datarray("transmittance.dat", TRANSMITTANCE_W, TRANSMITTANCE_H, NSSAMP);
    DATARRAY *delta_irradiance_dp = allocate_3d_datarray("delta_irradiance.dat", IRRADIANCE_W, IRRADIANCE_H, NSSAMP);
    DATARRAY *irradiance_dp = allocate_3d_datarray("irradiance.dat", IRRADIANCE_W, IRRADIANCE_H, NSSAMP);

    DATARRAY *delta_rayleigh_scattering_dp = allocate_5d_datarray("delta_rayleigh_scattering.dat", RES_R, RES_MU, RES_MU_S, RES_NU, NSSAMP);
    DATARRAY *delta_mie_scattering_dp = allocate_5d_datarray("delta_mie_scattering.dat", RES_R, RES_MU, RES_MU_S, RES_NU, NSSAMP);
    DATARRAY *scattering_dp = allocate_5d_datarray("scattering.dat", RES_R, RES_MU, RES_MU_S, RES_NU, NSSAMP);
    DATARRAY *delta_multiple_scattering_dp = allocate_5d_datarray("delta_multiple_scattering.dat", RES_R, RES_MU, RES_MU_S, RES_NU, NSSAMP);
    DATARRAY *delta_scattering_density_dp = allocate_5d_datarray("delta_scattering_density.dat", RES_R, RES_MU, RES_MU_S, RES_NU, NSSAMP);


    printf("computing transmittance...\n");
    // if (getpath(tau_path, getrlibpath(), R_OK) == NULL) {
    //     printf("writing transmittance data to file...\n");
    //     write_transmittance_data(tau_path);
    // };
    // tau_dp = getdata(tau_path);
    double radi, ctheta;
    double tau[NSSAMP];
    for (j = 0; j < TRANSMITTANCE_H; ++j) {
        for (i = 0; i < TRANSMITTANCE_W; ++i) {

            /* x, y to radius and cos(theta) */
            float xr = (i + 0.5) / TRANSMITTANCE_W;
            float yr = (j + 0.5) / TRANSMITTANCE_H;
            radi = ER + yr * yr * AH;
            ctheta = -0.15 + tan(1.5 * xr) / tan(1.5) * 1.15;

            // FILE *fp2 = fopen("tau.txt", "a");
            // fprintf(fp2, "%f %f %f %f", i+0.5, j+0.5, radi, ctheta);
            // fclose(fp2);

            idx = j * TRANSMITTANCE_W * NSSAMP + i * NSSAMP;
            compute_transmittance(radi, ctheta, tau);
            for (k = 0; k < NSSAMP; ++k) {
                tau_dp->arr.d[idx + k] = tau[k];
            }
        }
    }

    // double r = 6360e3;
    // double ct = 0.5;
    // double tau[NSSAMP];
    // interpolate_transmittance(tau, tau_dp, r, ct);
    // for (int i = 0; i < NSSAMP; ++i) {
    //     printf("tau: %f\n", tau[i]);
    // }
    //
    // Compute the direct irradiance, store it in delta_irradiance_texture, and
    // initialize irradiance_texture_ with zeros (we don't want the direct
    // irradiance in irradiance_texture_, but only the irradiance from the sky).
    printf("computing direct irradiance...\n");
    for (j = 0; j < IRRADIANCE_H; ++j) {
        for (i = 0; i < IRRADIANCE_W; ++i) {
            double r, mu_s;
            double result[NSSAMP];
            idx = j * IRRADIANCE_W * NSSAMP + i * NSSAMP;
            from_irradiance_uv(i + 0.5, j + 0.5, &r, &mu_s);
            compute_direct_irradiance(tau_dp, r, mu_s, result);
            for (k = 0; k < NSSAMP; ++k) {
                delta_irradiance_dp->arr.d[idx+k] = result[k];
                irradiance_dp->arr.d[idx+k] = 0.0;
            }
        }
    }


    // Compute the rayleigh and mie single scattering, and store them in
    // delta_rayleigh_scattering_texture and delta_mie_scattering_texture, as well
    // as in scattering_texture.
    // if (getpath(ins1r_path, getrlibpath(), R_OK) == NULL) {
    //     printf("writing transmittance data to file...\n");
    //     write_single_scattering_data(tau_dp, ins1r_path, ins1m_path);
    // };
    // delta_rayleigh_scattering_dp = getdata(ins1r_path);
    // delta_mie_scattering_dp = getdata(ins1m_path);
    // scattering_dp = getdata(ins1r_path);
    printf("computing single scattering...\n");
    for (i=0; i < RES_R; ++i) {
        for (j=0; j < RES_MU; ++j) {
            for (k=0; k < RES_MU_S; ++k) {
                for (l=0; l < RES_NU; ++l) {
                    double rayleigh[NSSAMP];
                    double mie[NSSAMP];
                    double r, mu, mu_s, nu;
                    int ray_r_mu_intersects_ground;
                    unsigned int index = i * RES_MU * RES_MU_S * RES_NU * NSSAMP + j * RES_MU_S*RES_NU*NSSAMP + k*RES_NU*NSSAMP + l*NSSAMP;
                    from_scattering_uvwz(i+0.5, j+0.5, k+0.5, l+0.5, &r, &mu, &mu_s, &nu, &ray_r_mu_intersects_ground);
                    compute_single_scattering(tau_dp, r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
                    for (m = 0; m < NSSAMP; ++m) {
                        delta_rayleigh_scattering_dp->arr.d[index+m] = rayleigh[m];
                        delta_mie_scattering_dp->arr.d[index+m] = mie[m];
                        scattering_dp->arr.d[index+m] = rayleigh[m];
                    }
                }
            }
        }
    }




    printf("computing multiple scattering...\n");
    // Compute the 2nd, 3rd and 4th order of scattering, in sequence.
    for (scattering_order = 2; scattering_order <= sorder; ++scattering_order) {
        // Compute the scattering density, and store it in delta_scattering_density_texture.
        for (unsigned int l = 0; l < RES_R; ++l) {
            for (unsigned int k = 0; k < RES_MU; ++k) {
                for (unsigned int j = 0; j < RES_MU_S; ++j) {
                    for (unsigned int i = 0; i < RES_NU; ++i) {
                        double scattering_density[NSSAMP];
                        compute_scattering_density(
                                tau_dp,
                                delta_rayleigh_scattering_dp,
                                delta_mie_scattering_dp,
                                delta_multiple_scattering_dp,
                                delta_irradiance_dp,
                                l, k, j, i,
                                scattering_order,
                                scattering_density);
                        idx = l * RES_MU * RES_MU_S * RES_NU * NSSAMP + k * RES_MU_S*RES_NU*NSSAMP + j*RES_NU*NSSAMP + i*NSSAMP;
                        for (int m = 0; m < NSSAMP; ++m) {
                            delta_scattering_density_dp->arr.d[idx + m] = scattering_density[m];
                        }
                    }
                }
            }
        }

        // Compute the indirect irradiance, store it in delta_irradiance_texture and accumulate it in irradiance_texture_.
        for (unsigned int j = 0; j < IRRADIANCE_H; ++j) {
            for (unsigned int i = 0; i < IRRADIANCE_W; ++i) {
                double delta_irradiance[NSSAMP];
                double r, mu_s;
                from_irradiance_uv(i + 0.5, j + 0.5, &r, &mu_s);
                compute_indirect_irradiance(
                        delta_rayleigh_scattering_dp,
                        delta_mie_scattering_dp, delta_multiple_scattering_dp,
                        r, mu_s, scattering_order-1, delta_irradiance);
                idx = j * IRRADIANCE_W * NSSAMP + i * NSSAMP;
                for (int k = 0; k < NSSAMP; ++k) {
                    delta_irradiance_dp->arr.d[idx + k] = delta_irradiance[k];
                }
            }
        }

        increment_dp(irradiance_dp, delta_irradiance_dp);

        // Compute the multiple scattering, store it in delta_multiple_scattering_texture, and accumulate it in scattering_texture_.
        for (unsigned int l = 0; l < RES_R; ++l) {
            for (unsigned int k = 0; k < RES_MU; ++k) {
                for (unsigned int j = 0; j < RES_MU_S; ++j) {
                    for (unsigned int i = 0; i < RES_NU; ++i) {
                        double delta_multiple_scattering[NSSAMP];
                        double r, mu, mu_s, nu;
                        int ray_r_mu_intersects_ground;
                        from_scattering_uvwz(l+0.5, k+0.5, j+0.5, i+0.5, &r, &mu, &mu_s, &nu, &ray_r_mu_intersects_ground);
                        compute_multiple_scattering(tau_dp, delta_scattering_density_dp, r, mu, mu_s, nu, ray_r_mu_intersects_ground, delta_multiple_scattering);
                        idx = l * RES_MU * RES_MU_S * RES_NU * NSSAMP + k * RES_MU_S*RES_NU*NSSAMP + j*RES_NU*NSSAMP + i*NSSAMP;
                        for (m = 0; m < NSSAMP; ++m) {
                            delta_multiple_scattering_dp->arr.d[idx + m] = delta_multiple_scattering[m];
                            scattering_dp->arr.d[idx + m] += delta_multiple_scattering[m] * (1.0 / rayleigh_phase_function(nu));
                        }
                    }
                }
            }
        }
    }


    // transmittance_texture_->Save(cache_directory_ + "transmittance.dat");
    savedata(tau_dp);
    savedata(scattering_dp);
    savedata(delta_mie_scattering_dp);
    savedata(irradiance_dp);

    freedata(tau_dp);
    freedata(delta_rayleigh_scattering_dp);
    freedata(delta_mie_scattering_dp);
    freedata(scattering_dp);
    freedata(delta_multiple_scattering_dp);
    freedata(delta_scattering_density_dp);
    freedata(delta_irradiance_dp);
    freedata(irradiance_dp);
    return 1;
}


int main(int argc, char *argv[])
{
    progname = "atmos_precompute";
    printf("precomputing...\n");
    if (!precompute(2)) {
        printf("precompute failed\n");
        return 1;
    }
    return 0;
}
