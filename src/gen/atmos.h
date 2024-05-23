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
#include "resolu.h"
#include "rtio.h"
#include "rtmath.h"
#include "sun.h"
#include "view.h"

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <pthread.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define THREAD HANDLE
#define CREATE_THREAD(thread, func, context)                                   \
  ((*(thread) = CreateThread(NULL, 0, (LPTHREAD_STRAT_ROUTINE)(func),          \
                             (context), 0, NULL)) != NULL)
#define THREAD_RETURN DWORD WINAPI
#define THREAD_JOIN(thread) WaitForSingleObject(thread, INFINITE)
#define THREAD_CLOSE(thread) CloseHandle(thread)
#else
#define THREAD pthread_t
#define CREATE_THREAD(thread, func, context)                                   \
  (pthread_create((thread), NULL, (func), (context)) == 0)
#define THREAD_RETURN void *
#define THREAD_JOIN(thread) pthread_join(thread, NULL)
#endif

#define NSSAMP 20

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
  DensityProfile rayleigh_density;
  DensityProfile mie_density;
  DensityProfile ozone_density;
  const float *beta_r0;
  float beta_scale;
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
extern const double AOD0_CA;
extern const double SOLOMG;

extern void get_sky_radiance(DATARRAY *tau_dp, DATARRAY *scat_dp,
                             DATARRAY *scat1m_dp, FVECT camera, FVECT view_ray,
                             double shadow_length, FVECT sundir,
                             float *transmittance, float *result);

extern int compute_sundir(const int year, const int month, const int day,
                          const double hour, const int tsolar,
                          double sundir[3]);

extern int precompute(const int sorder, const Atmosphere *atmos, int num_threads);

#endif // ATMOS_H
