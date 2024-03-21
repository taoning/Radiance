#include "sun.h"
#include "fvect.h"
#include "data.h"
#include "color.h"

#define NSSAMP 20

#define ER 6360.0e3
#define AR 6420.0e3
#define AH 60.0e3

const double MU_S_MIN = -0.2;

const double GROUND_ALBEDO = 0.2;

// Scale heights (m)
// Rayleigh scattering
const double HR_MS = 8820; // Midlatitude summer
const double HR_MW = 8100; // Midlatitude winter
const double HR_SS = 8550; // Subarctic summer
const double HR_SW = 7600; // Subarctic winter
const double HR_T = 9050;  // Tropical
// Mie scattering
const double HMS_CC = 2385; // Continental clean
// Mie absorption
const double HMA_CC = 4300; // Continental clean

const double SOLOMG = 6.7967e-5; // .533 apex angle
const int WVLSPAN = 400;         // 380nm to 780nm

// Aerosol optical depth at 550nm for continental clean
const double AOD_CC_550 = 0.05352;

/** 380nm to 780nm at 20nm intervals **/
// Extraterrestrial solar W/m^2/nm
const float EXTSOL[NSSAMP] = {1.11099345, 1.75363199, 1.67126921, 1.97235667,
                              2.01863015, 1.93612482, 1.87310758, 1.88792588,
                              1.86274213, 1.84433248, 1.80370942, 1.73020456,
                              1.67036654, 1.57482149, 1.53646383, 1.45515319,
                              1.38207435, 1.3265141,  1.26648111, 1.20592043};

// Rayleight scattering coefficients at sea level in m^-1
const float BR0_MS[NSSAMP] = {
    4.63889967e-05, 3.76703293e-05, 3.09193474e-05, 2.56223081e-05,
    2.14198405e-05, 1.80483398e-05, 1.53177709e-05, 1.30867298e-05,
    1.12487269e-05, 9.72370830e-06, 8.44932642e-06, 7.37769845e-06,
    6.47116653e-06, 5.70005327e-06, 5.04091603e-06, 4.47430240e-06,
    3.98549839e-06, 3.56178512e-06, 3.19293761e-06, 2.87072599e-06};
const float BR0_MW[NSSAMP] = {
    5.03857564e-05, 4.09159104e-05, 3.35832809e-05, 2.78298619e-05,
    2.32653203e-05, 1.96033395e-05, 1.66375116e-05, 1.42142496e-05,
    1.22178890e-05, 1.05614786e-05, 9.17729921e-06, 8.01334246e-06,
    7.02870602e-06, 6.19115557e-06, 5.47522872e-06, 4.85979708e-06,
    4.32887894e-06, 3.86865960e-06, 3.46803311e-06, 3.11806054e-06};

const float BR0_SS[NSSAMP] = {
    4.73789047e-05, 3.84741872e-05, 3.15791442e-05, 2.61690698e-05,
    2.18769246e-05, 1.84334785e-05, 1.56446412e-05, 1.33659912e-05,
    1.14887667e-05, 9.93120528e-06, 8.62962901e-06, 7.53513326e-06,
    6.60925660e-06, 5.82168833e-06, 5.14848557e-06, 4.56978082e-06,
    4.07054608e-06, 3.63779107e-06, 3.26107261e-06, 2.93198523e-06};
const float BR0_SW[NSSAMP] = {
    5.30623659e-05, 4.30894595e-05, 3.53673035e-05, 2.93082494e-05,
    2.45012287e-05, 2.06447149e-05, 1.75213353e-05, 1.49693438e-05,
    1.28669320e-05, 1.11225291e-05, 9.66481886e-06, 8.43902999e-06,
    7.40208736e-06, 6.52004427e-06, 5.76608570e-06, 5.11796089e-06,
    4.55883914e-06, 4.07417187e-06, 3.65226316e-06, 3.28369923e-06};
const float BR0_T[NSSAMP] = {
    4.55376661e-05, 3.69790036e-05, 3.03519157e-05, 2.51520876e-05,
    2.10267438e-05, 1.77171168e-05, 1.50366593e-05, 1.28465622e-05,
    1.10422904e-05, 9.54525887e-06, 8.29426444e-06, 7.24230298e-06,
    6.35240772e-06, 5.59544593e-06, 4.94840517e-06, 4.39219003e-06,
    3.91235655e-06, 3.49641927e-06, 3.13434084e-06, 2.81804244e-06};

// OPAC Mie scattering coefficient at sea level in m^-1
// continental clean
const float BM0_CC[NSSAMP] = {2.8278e-05, 2.6766e-05, 2.5349e-05, 2.3945e-05,
                              2.2771e-05, 2.1593e-05, 2.0508e-05, 1.9494e-05,
                              1.8466e-05, 1.7606e-05, 1.6766e-05, 1.5987e-05,
                              1.5268e-05, 1.4521e-05, 1.3935e-05, 1.3328e-05,
                              1.2748e-05, 1.2181e-05, 1.1641e-05, 1.112e-05};

const float BM0_CA[NSSAMP] = {7.8237e-05, 7.3966e-05, 6.9974e-05, 6.6022e-05,
                              6.2729e-05, 5.9427e-05, 5.6391e-05, 5.3561e-05,
                              5.0694e-05, 4.8299e-05, 4.596e-05,  4.3797e-05,
                              4.1802e-05, 3.9726e-05, 3.8103e-05, 3.6423e-05,
                              3.4817e-05, 3.3253e-05, 3.1763e-05, 3.033e-05};

const float AM0_CC[NSSAMP] = {7.12e-7,  6.83e-07, 6.56e-07, 6.3e-07,  6.09e-07,
                              5.87e-07, 5.79e-07, 5.84e-07, 5.89e-07, 5.71e-07,
                              5.54e-07, 5.47e-07, 5.51e-07, 5.54e-07, 5.4e-07,
                              5.26e-07, 5.24e-07, 5.35e-07, 5.45e-07, 5.52e-07};

const float AM0_CA[NSSAMP] = {
    7.671e-06, 7.3e-06,   6.961e-06, 6.625e-06, 6.347e-06, 6.068e-06, 5.839e-06,
    5.654e-06, 5.466e-06, 5.262e-06, 5.063e-06, 4.909e-06, 4.8e-06,   4.685e-06,
    4.538e-06, 4.386e-06, 4.279e-06, 4.218e-06, 4.159e-06, 4.104e-06};

const double MIE_G = 0.76;

const int TRANSMITTANCE_TEXTURE_WIDTH = 256;
const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;

const int SCATTERING_TEXTURE_R_SIZE = 32;
const int SCATTERING_TEXTURE_MU_SIZE = 128;
const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
const int SCATTERING_TEXTURE_NU_SIZE = 8;

// const int SCATTERING_TEXTURE_R_SIZE = 16;
// const int SCATTERING_TEXTURE_MU_SIZE = 64;
// const int SCATTERING_TEXTURE_MU_S_SIZE = 16;
// const int SCATTERING_TEXTURE_NU_SIZE = 4;

const int IRRADIANCE_TEXTURE_WIDTH = 64;
const int IRRADIANCE_TEXTURE_HEIGHT = 16;

extern const double HR_MS;
extern const double HR_MW;
extern const double HR_SS;
extern const double HR_SW;
extern const double HR_T;
extern const double HMS_CC;
extern const double HMA_CC;

extern const float EXTSOL[NSSAMP];
extern const float BR0_MS[NSSAMP];
extern const float BR0_MW[NSSAMP];
extern const float BR0_SS[NSSAMP];
extern const float BR0_SW[NSSAMP];
extern const float BR0_T[NSSAMP];
extern const float BM0_CC[NSSAMP];
extern const float AM0_CC[NSSAMP];
extern const float BM0_CA[NSSAMP];
extern const float AM0_CA[NSSAMP];


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


typedef struct {
  double aod;
  double ccover;
} Weather;


typedef struct {
  const double scaleheight;
  const float *scattering_coefficients0;
} RayleighAtmos;

typedef struct {
  const double scaleheight_s;
  const double scaleheight_a;
  const float *scoeffs0;
  const float *acoeffs0;
  const double aod;
} MieAtmos;

typedef struct {
  const double alpha;
  const double beta;
} AngstromCoefficients;



typedef struct {
  const RayleighAtmos *rayleigh;
  const MieAtmos *mie;
} Atmosphere;

int compute_sundir(const Datetime *dts, const Location *lc, const int tsolar, double sundir[3]);

int precompute();

int get_sky_radiance();
