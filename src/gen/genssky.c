#include <assert.h> // to be removed
#include <math.h>
#include <stdio.h>

#include "atmos.h"
#include "color.h"
#include "data.h"
#include "fvect.h"
#include "paths.h"
#include "resolu.h"
#include "rtio.h"
#include "rtmath.h"
#include "sun.h"
#include "view.h"

char *progname;

const double D65EFF = 203.; /* standard illuminant D65 */

// Relative daylight spectra where CCT = 6415K for overcast;
const double D6415[NSSAMP] = {0.63231, 1.06171, 1.00779, 1.36423, 1.34133,
                              1.27258, 1.26276, 1.26352, 1.22201, 1.13246,
                              1.0434,  1.05547, 0.98212, 0.94445, 0.9722,
                              0.82387, 0.87853, 0.82559, 0.75111, 0.78925};

inline static double wmean2(double a, double b, double x) {
  return a * (1 - x) + b * x;
}

inline static double wmean(double a, double x, double b, double y) {
  return (a * x + b * y) / (a + b);
}

// from gensky.c
static double get_overcast_brightness(const double dz, const double sundir[3]) {
  double zenithbr;
  if (sundir[2] < 0) {
    zenithbr = 0;
  } else {
    zenithbr = (8.6 * sundir[2] + .123) * 1000.0 / D65EFF;
  }
  double groundbr = zenithbr * 0.777778;
  return wmean(pow(dz + 1.01, 10), zenithbr * (1 + 2 * dz) / 3,
               pow(dz + 1.01, -10), groundbr);
}

int gen_spect_sky(DATARRAY *tau_clear, DATARRAY *scat_clear,
                  DATARRAY *scat1m_clear, double cloud_cover, FVECT sundir,
                  char *outname) {

  double rmu;
  FVECT view_point = {0, 0, ER + 1};
  double point_radius = VLEN(view_point);
  double sun_ct = fdot(view_point, sundir) / point_radius;
  char hsrfile[256];
  if (!snprintf(hsrfile, sizeof(hsrfile), "%s.hsr", outname)) {
    fprintf(stderr, "Error creating header file name\n");
    return 0;
  }
  FILE *hfp = fopen(hsrfile, "w");
  float wvsplit[4] = {380, 480, 588,
                      780}; /* RGB wavelength limits+partitions (nm) */
  newheader("RADIANCE", hfp);
  fputncomp(NSSAMP, hfp);
  fputwlsplit(wvsplit, hfp);
  fputformat(SPECFMT, hfp);
  fputc('\n', hfp);
  int yres = 32;
  int xres = 32;
  RESOLU rs = {PIXSTANDARD, xres, yres};
  RREAL loc[2];
  double d;
  FVECT rorg = {0};
  FVECT rdir = {0};
  VIEW vw = {VT_ANG, {0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}, 1.,
             180.,   180.,         0.,           0.,           0.,
             0.,     {0., 0., 0.}, {0., 0., 0.}, 0.,           0.};
  setview(&vw);
  fputsresolu(&rs, hfp);
  CNDX[3] = NSSAMP;
  for (unsigned int j = 0; j < yres; ++j) {
    for (unsigned int i = 0; i < xres; ++i) {
      double u, v, w, z;
      double r, mu, mu_s, nu;
      SCOLR sclr;
      pix2loc(loc, &rs, i, j);
      d = viewray(rorg, rdir, &vw, loc[0], loc[1]);
      rmu = fdot(view_point, rdir);
      mu = rmu / point_radius;
      nu = fdot(rdir, sundir);
      to_scattering_uvwz(point_radius, mu, sun_ct, nu, 0, &u, &v, &w, &z);
      double pt[4] = {z, w, v, u};
      SCOLOR result = {0};
      get_sky_radiance(scat_clear, scat1m_clear, nu, pt, result);
      for (int k = 0; k < NSSAMP; ++k) {
        result[k] *= WVLSPAN;
      }
      if (cloud_cover > 0) {
        double skybr = get_overcast_brightness(rdir[2], sundir);
        for (int k = 0; k < NSSAMP; ++k) {
          result[k] = wmean2(result[k], skybr * D6415[k], cloud_cover);
        }
      }
      scolor2scolr(sclr, result, 20);
      putbinary(sclr, LSCOLR, 1, hfp);
    }
  }
  fclose(hfp);

  // Get solar radiance
  double nus = fdot(sundir, sundir);
  double us, vs, ws, zs;
  to_scattering_uvwz(point_radius, sun_ct, sun_ct, nus, 0, &us, &vs, &ws, &zs);
  double pts[4] = {zs, ws, vs, us};
  float trans_sun[NSSAMP] = {0};
  float sky_radiance_sun[NSSAMP] = {0};
  get_sky_transmittance(tau_clear, point_radius, sun_ct, trans_sun);
  get_sky_radiance(scat_clear, scat1m_clear, nus, pts, sky_radiance_sun);
  double sun_radiance[NSSAMP] = {0};
  for (int i = 0; i < NSSAMP; ++i) {
    sun_radiance[i] = sky_radiance_sun[i] + trans_sun[i] * EXTSOL[i] / SOLOMG;
  }
  if (cloud_cover > 0) {
    double skybr = get_overcast_brightness(sundir[2], sundir);
    for (int i = 0; i < NSSAMP; ++i) {
      sun_radiance[i] =
          wmean2(sun_radiance[i], D6415[i] * skybr / WVLSPAN, cloud_cover);
    }
  }
  char radfile[256];
  if (!snprintf(radfile, sizeof(radfile), "%s.rad", outname)) {
    fprintf(stderr, "Error creating rad file name\n");
    return 0;
  }
  FILE *rfp = fopen(radfile, "w");
  fprintf(rfp, "void spectrum sunrad\n0\n0\n22 380 780 ");
  for (int i = 0; i < NSSAMP; ++i) {
    fprintf(rfp, "%.1f ", sun_radiance[i] * WVLSPAN);
  }
  fprintf(rfp, "\n\nsunrad light solar\n0\n0\n3 1 1 1\n\n");
  fprintf(rfp, "solar source sun\n0\n0\n4 %f %f %f 0.533\n\n", sundir[0],
          sundir[1], sundir[2]);
  fprintf(rfp,
          "void specpict skyfunc\n8 noop %s fisheye.cal fish_u fish_v -rx 90 "
          "-mx\n0\n0\n\n",
          hsrfile);
  fprintf(rfp, "skyfunc glow sky_glow\n0\n0\n4 1 1 1 0\n\n");
  fprintf(rfp, "sky_glow source sky\n0\n0\n4 0 0 1 180\n\n");
  return 1;
}

int main(int argc, char *argv[]) {
  progname = argv[0];
  const double arctic_circle_latitude = 67.;
  const double tropic_latitude = 23.;
  const int summer_start_month = 4;
  const int summer_end_month = 9;
  int month, day;
  double hour;
  FVECT sundir;
  int num_threads = 1;
  int sorder = 4;
  int year = 0;
  int tsolar = 0;
  double hsm = 0.0;
  double ccover = 0.0;
  double aod = AOD0_CA;
  char *outname = "out";
  char dirname[PATH_MAX];
  snprintf(dirname, sizeof(dirname), "ssky_%f", aod);
  char *mie_ca_path = getpath("mie_ca.dat", getrlibpath(), R_OK);

  if (argc < 5) {
    fprintf(stderr, "Usage: %s month day hour -y year -a lat -o lon -m tz\n",
            argv[0]);
    return 0;
  }
  if (!strcmp(argv[1], "-ang")) {
    float altitude = atof(argv[2]) * (M_PI / 180);
    float azimuth = atof(argv[3]) * (M_PI / 180);
    month = 0;
    sundir[0] = -sin(azimuth) * cos(altitude);
    sundir[1] = -cos(azimuth) * cos(altitude);
    sundir[2] = sin(altitude);
  } else {
    month = atoi(argv[1]);
    day = atoi(argv[2]);
    hour = atof(argv[3]);
    if (!compute_sundir(year, month, day, hour, tsolar, sundir)) {
      fprintf(stderr, "Cannot compute solar angle\n");
      exit(1);
    }
  }
  for (int i = 4; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'a':
        s_latitude = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'b':
        hsm = atof(argv[++i]);
        break;
      case 'c':
        ccover = atof(argv[++i]);
        break;
      case 'd':
        aod = atof(argv[++i]);
        break;
      case 'i':
        sorder = atoi(argv[++i]);
        break;
      case 'm':
        s_meridian = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'o':
        s_longitude = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'n':
        num_threads = atoi(argv[++i]);
        break;
      case 'y':
        year = atoi(argv[++i]);
        break;
      case 'r':
        outname = argv[++i];
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", argv[i]);
        exit(1);
      }
    }
  }

  // Determine if it's summer for the given hemisphere
  int is_northern_hemisphere = (s_latitude >= 0);
  int is_summer = (month >= summer_start_month && month <= summer_end_month);
  if (!is_northern_hemisphere) {
    is_summer = !is_summer;
  }

  Atmosphere clear_atmos = {
      .ozone_density = {.layers =
                            {
                                {.width = 25000.0,
                                 .exp_term = 0.0,
                                 .exp_scale = 0.0,
                                 .linear_term = 1.0 / 15000.0,
                                 .constant_term = -2.0 / 3.0},
                                {.width = AH,
                                 .exp_term = 0.0,
                                 .exp_scale = 0.0,
                                 .linear_term = -1.0 / 15000.0,
                                 .constant_term = 8.0 / 3.0},
                            }},
      .rayleigh_density = {.layers =
                               {
                                   {.width = AH,
                                    .exp_term = 1.0,
                                    .exp_scale = -1.0 / HR_MS,
                                    .linear_term = 0.0,
                                    .constant_term = 0.0},
                               }},
      .beta_r0 = BR0_MS,
      .beta_scale = aod / AOD0_CA,
      .beta_m = NULL,
  };

  // Set rayleigh density profile
  if (fabs(s_latitude) > arctic_circle_latitude) {
    if (is_summer) {
      clear_atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_SS;
      clear_atmos.beta_r0 = BR0_SS;
    } else {
      clear_atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_SW;
      clear_atmos.beta_r0 = BR0_SW;
    }
  } else if (fabs(s_latitude) > tropic_latitude) {
    if (is_summer) {
      clear_atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_MS;
      clear_atmos.beta_r0 = BR0_MS;
    } else {
      clear_atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_MW;
      clear_atmos.beta_r0 = BR0_MW;
    }
  } else {
    clear_atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_T;
    clear_atmos.beta_r0 = BR0_T;
  }

  // Load mie density data
  DATARRAY *mie_ca_dp = getdata(mie_ca_path);
  if (mie_ca_dp == NULL) {
    fprintf(stderr, "Error reading mie data\n");
    return 0;
  }
  clear_atmos.beta_m = mie_ca_dp;

  DpPaths clear_paths = {
      .tau = "tau_clear.dat",
      .scat = "scat_clear.dat",
      .scat1m = "scat1m_clear.dat",
      .irrad = "irrad_clear.dat",
  };

  if (getpath(clear_paths.tau, getrlibpath(), R_OK) == NULL ||
      getpath(clear_paths.scat, getrlibpath(), R_OK) == NULL ||
      getpath(clear_paths.scat1m, getrlibpath(), R_OK) == NULL ||
      getpath(clear_paths.irrad, getrlibpath(), R_OK) == NULL) {
    printf("# Precomputing...\n");
    if (!precompute(sorder, &clear_paths, &clear_atmos, num_threads)) {
      printf("Clear precompute failed\n");
      return 0;
    }
  }

  DATARRAY *tau_clear_dp = getdata(clear_paths.tau);
  DATARRAY *scat_clear_dp = getdata(clear_paths.scat);
  DATARRAY *scat1m_clear_dp = getdata(clear_paths.scat1m);

  if (!gen_spect_sky(tau_clear_dp, scat_clear_dp, scat1m_clear_dp, ccover,
                     sundir, outname)) {
    fprintf(stderr, "gen_spect_sky failed\n");
    exit(1);
  }

  freedata(mie_ca_dp);
  freedata(tau_clear_dp);
  freedata(scat_clear_dp);
  freedata(scat1m_clear_dp);

  return 1;
}
