#ifndef _ATMOSPHERE_H
#define _ATMOSPHERE_H

#include "color.h"
#include "rtmath.h"

extern const double HR_MS;
extern const double HR_MW;
extern const double HR_SS;
extern const double HR_SW;
extern const double HR_T;
extern const double HMS_CC;
extern const double HMA_CC;

extern const SCOLOR EXTSOL;
extern const SCOLOR BR0_MS;
extern const SCOLOR BR0_MW;
extern const SCOLOR BR0_SS;
extern const SCOLOR BR0_SW;
extern const SCOLOR BR0_T;
extern const SCOLOR BM0_CC;
extern const SCOLOR AM0_CC;
extern const SCOLOR BM0_CA;
extern const SCOLOR AM0_CA;


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

int tracesky(SCOLOR scolor, const FVECT dir, const FVECT sundir, const Atmosphere *atmos);
int tracesun(SCOLOR suncolor, const FVECT sundir, const Atmosphere *atmos);


#endif // _ATMOSPHERE_H
