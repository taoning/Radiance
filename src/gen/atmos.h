#ifndef ATMOS_H
#define ATMOS_H

#include <assert.h> // to be removed
#include <math.h>
#include <stdio.h>
#include <time.h> // to be removed

#include "color.h"
#include "data.h"
#include "fvect.h"
#include "paths.h"
#include "rtio.h"
#include "rtmath.h"
#include "sun.h"

#define NSSAMP 20

typedef struct {
  double width;
  double exp_term;
  double exp_scale;
  double linear_term;
  double constant_term;
} DensityProfileLayer;

typedef struct {
  DensityProfileLayer layers[3];
} DensityProfile;

typedef struct {
  DensityProfile rayleigh_density;
  DensityProfile cloud_density;
  DensityProfile ozone_density;
  const float *beta_r0;
  const float beta_c;
  float beta_scale;
  float cloud_cover;
  DATARRAY *beta_m;
} Atmosphere;

extern const double ER;
extern const double AH;
extern const double HR_MS;
extern const double HR_MW;
extern const double HR_SS;
extern const double HR_SW;
extern const double HR_T;
extern const int WVLSPAN;
extern const float EXTSOL[NSSAMP];
extern const float BR0_MS[NSSAMP];
extern const float BR0_MW[NSSAMP];
extern const float BR0_SS[NSSAMP];
extern const float BR0_SW[NSSAMP];
extern const float BR0_T[NSSAMP];
extern const float BCLOUD;
extern const double AOD0_CA;
extern const double SOLOMG;

extern void get_sky_radiance(DATARRAY *tau_dp, DATARRAY *scat_dp,
                             DATARRAY *scat1m_dp, FVECT camera, FVECT view_ray,
                             double shadow_length, FVECT sundir,
                             float *transmittance, float *result);

extern int compute_sundir(const int year, const int month, const int day,
                          const double hour, const int tsolar,
                          double sundir[3]);

extern int precompute(const int sorder, const Atmosphere *atmos,
                      int num_threads);

#endif // ATMOS_H
