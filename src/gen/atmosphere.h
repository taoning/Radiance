#ifndef _ATMOSPHERE_H
#define _ATMOSPHERE_H

#include "color.h"
#include "paths.h"
#include "rtio.h"
#include "rtmath.h"

#define NSSAMP 20

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
  const double scaleheight;
  const float *scoeffs0;
} RayleighAtmos;

typedef struct {
  const double scaleheight_s;
  const double scaleheight_a;
  const float *scoeffs0;
  const float *acoeffs0;
  const double aod;
} MieAtmos;

typedef struct {
  const RayleighAtmos *rayleigh;
  const MieAtmos *mie;
} Atmosphere;

int tracesky(float *scolor, const FVECT dir, const FVECT sundir,
             const Atmosphere *atmos);
int tracesun(float *suncolor, const FVECT sundir, const Atmosphere *atmos);

#endif // _ATMOSPHERE_H
