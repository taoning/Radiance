#include "atmos.h"
#include "paths.h"
#include "rtmath.h"

char *progname;

inline void vectorize(double altitude, double azimuth, FVECT v) {
  v[1] = cos(altitude);
  v[0] = (v)[1] * sin(azimuth);
  v[1] *= cos(azimuth);
  v[2] = sin(altitude);
}

int rh_init(int rhsubdiv, float *rh_palt, float *rh_pazi, float *rh_dom) {
  const int NROW = 7;
  static const int tnaz[7] = {30, 30, 24, 24, 18, 12, 6};
  const double alpha = (PI / 2.) / (NROW * rhsubdiv + .5);
  int p, i, j;
  /* allocate patch angle arrays */
  int nskypatch = 0;
  for (p = 0; p < NROW; p++)
    nskypatch += tnaz[p];
  nskypatch *= rhsubdiv * rhsubdiv;
  nskypatch += 2;
  rh_palt = (float *)malloc(sizeof(float) * nskypatch);
  rh_pazi = (float *)malloc(sizeof(float) * nskypatch);
  rh_dom = (float *)malloc(sizeof(float) * nskypatch);
  if ((rh_palt == NULL) | (rh_pazi == NULL) | (rh_dom == NULL)) {
    fprintf(stderr, "gensskymtx: out of memory in rh_init()\n");
    exit(1);
  }
  rh_palt[0] = -PI / 2.; /* ground & zenith patches */
  rh_pazi[0] = 0.;
  rh_dom[0] = 2. * PI;
  rh_palt[nskypatch - 1] = PI / 2.;
  rh_pazi[nskypatch - 1] = 0.;
  rh_dom[nskypatch - 1] = 2. * PI * (1. - cos(alpha * .5));
  p = 1; /* "normal" patches */
  for (i = 0; i < NROW * rhsubdiv; i++) {
    const float ralt = alpha * (i + .5);
    const int ninrow = tnaz[i / rhsubdiv] * rhsubdiv;
    const float dom =
        2. * PI * (sin(alpha * (i + 1)) - sin(alpha * i)) / (double)ninrow;
    for (j = 0; j < ninrow; j++) {
      rh_palt[p] = ralt;
      rh_pazi[p] = 2. * PI * j / (double)ninrow;
      rh_dom[p++] = dom;
    }
  }
  return nskypatch;
}

/* Resize daylight matrix (GW) */
float *resize_dmatrix(float *mtx_data, int nsteps, int npatch) {
  if (mtx_data == NULL)
    mtx_data = (float *)malloc(sizeof(float) * NSSAMP * nsteps * npatch);
  else
    mtx_data =
        (float *)realloc(mtx_data, sizeof(float) * NSSAMP * nsteps * npatch);
  if (mtx_data == NULL) {
    fprintf(stderr, "%s: out of memory in resize_dmatrix(%d,%d)\n", progname,
            nsteps, npatch);
    exit(1);
  }
  return (mtx_data);
}

/* Add in solar direct to nearest sky patches (GW) */
void add_direct(const double altitude, const double azimuth, const float cc,
                const int nskypatch, const float *rh_alt, const float *rh_azi,
                float *parr) {
  const int NSUNPATCH = 4;
  int nsuns = NSUNPATCH;
  FVECT svec;
  double near_dprod[NSUNPATCH];
  int near_patch[NSUNPATCH];
  double wta[NSUNPATCH], wtot;
  int i, j, p;

  if ((1 - cc) <= 1e-4)
    return;
  /* identify nsuns closest patches */
  if (nsuns > NSUNPATCH)
    nsuns = NSUNPATCH;
  else if (nsuns <= 0)
    nsuns = 1;
  for (i = nsuns; i--;)
    near_dprod[i] = -1.;
  vectorize(altitude, azimuth, svec);
  for (p = 1; p < nskypatch; p++) {
    FVECT pvec;
    double dprod;
    vectorize(rh_alt[p], rh_azi[p], pvec);
    dprod = DOT(pvec, svec);
    for (i = 0; i < nsuns; i++)
      if (dprod > near_dprod[i]) {
        for (j = nsuns; --j > i;) {
          near_dprod[j] = near_dprod[j - 1];
          near_patch[j] = near_patch[j - 1];
        }
        near_dprod[i] = dprod;
        near_patch[i] = p;
        break;
      }
  }
  wtot = 0; /* weight by proximity */
  for (i = nsuns; i--;)
    wtot += wta[i] = 1. / (1.002 - near_dprod[i]);
  /* add to nearest patch radiances */
  for (i = nsuns; i--;) {
    float *pdest = parr + 3 * near_patch[i];
    float val_add = wta[i] * dir_illum / (WHTEFFICACY * wtot);

    val_add /= (fixed_sun_sa > 0) ? fixed_sun_sa : rh_dom[near_patch[i]];
    *pdest++ += val_add * suncolor[0];
    *pdest++ += val_add * suncolor[1];
    *pdest++ += val_add * suncolor[2];
  }
}

void calc_sky_patch_radiance(const int nskypatch, const float *rh_palt,
                             const float *rh_pazi, float *radiance) {
  int i;
  for (i = 0; i < nskypatch; i++) {
    printf("rh_palt[%d] = %f\n", i, rh_palt[i]);
    printf("rh_pazi[%d] = %f\n", i, rh_pazi[i]);
    get_sky_radiance(rh_palt[i], rh_pazi[i], radiance);
  }
}

int main(int argc, char *argv[]) {
  progname = argv[0];
  float *mtx_data = NULL;
  int mtx_offset = 0;
  char buf[256];
  int rhsubdiv = 1;
  int input = 0;
  double rotation = 0.0;
  double elevation = 0;
  int ntsteps = 0;  /* number of time steps */
  int tstorage = 0; /* number of allocated time steps */
  int nstored = 0;  /* number of time steps in matrix */
  int sun_hours_only = 0;
  int mo, da;
  double hr, aod, cc;
  float *rh_palt, *rh_pazi, *rh_dom;
  for (int i = 1; i < argc && argv[i][0] == '-'; i++) {
    switch (argv[i][1]) {
    case 'm':
      rhsubdiv = atoi(argv[++i]);
      break;
    case 'r':
      rotation = atof(argv[++i]);
      break;
    case 'u': /* solar hours only */
      sun_hours_only = 1;
      break;
    }
  }
  /* read weather tape header */
  if (scanf("place %[^\r\n] ", buf) != 1)
    goto fmterr;
  if (scanf("latitude %lf\n", &s_latitude) != 1)
    goto fmterr;
  if (scanf("longitude %lf\n", &s_longitude) != 1)
    goto fmterr;
  if (scanf("time_zone %lf\n", &s_meridian) != 1)
    goto fmterr;
  if (scanf("site_elevation %lf\n", &elevation) != 1)
    goto fmterr;
  if (scanf("weather_data_file_units %d\n", &input) != 1)
    goto fmterr;
  int nskypatch = rh_init(rhsubdiv, rh_palt, rh_pazi, rh_dom);
  printf("nskypatch = %d\n", nskypatch);
  /* initial allocation */
  mtx_data = resize_dmatrix(mtx_data, tstorage = 2, nskypatch);
  while (scanf("%d %d %lf %lf %lf\n", &mo, &da, &hr, &aod, &cc) == 5) {
    FVECT sundir;
    compute_sundir(0, mo, da, hr, 0, sundir);
    if (sun_hours_only && sundir[2] <= 0.)
      continue; /* skipping nighttime points */
    mtx_offset = NSSAMP * nskypatch * nstored;
    nstored += !nstored;
    /* make space for next row */
    if (nstored > tstorage) {
      tstorage += (tstorage >> 1) + nstored + 7;
      mtx_data = resize_dmatrix(mtx_data, tstorage, nskypatch);
    }
    ntsteps++; /* keep count of time steps */
  }
  return 0;
userr:
  fprintf(stderr,
          "Usage: %s [-v][-h][-A][-d|-s|-n][-u][-D file [-M modfile]][-r "
          "deg][-m N][-g r g b][-c r g b][-o{f|d}][-O{0|1}] [tape.wea]\n",
          progname);
  exit(1);
fmterr:
  fprintf(stderr, "%s: weather tape format error in header\n", progname);
  exit(1);
writerr:
  fprintf(stderr, "%s: write error on output\n", progname);
  exit(1);
}
