#include <assert.h> // to be removed
#include <math.h>
#include <stdio.h>
#include <time.h> // to be removed

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

int gen_spect_sky(DATARRAY *tau_dp, DATARRAY *scat_dp, DATARRAY *scat1m_dp,
                  FVECT sundir, char *outname) {
  FVECT camera = {0, 0, ER + 1};
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
      pix2loc(loc, &rs, i, j);
      d = viewray(rorg, rdir, &vw, loc[0], loc[1]);
      float trans[NSSAMP] = {0};
      SCOLOR results = {0};
      get_sky_radiance(tau_dp, scat_dp, scat1m_dp, camera, rdir, 0.0, sundir,
                       trans, results);
      for (int k = 0; k < NSSAMP; ++k) {
        results[k] *= WVLSPAN;
      }
      SCOLR sclr;
      scolor2scolr(sclr, results, 20);
      putbinary(sclr, LSCOLR, 1, hfp);
    }
  }
  fclose(hfp);
  float trans_sun[NSSAMP] = {0};
  float sky_radiance_sun[NSSAMP] = {0};
  get_sky_radiance(tau_dp, scat_dp, scat1m_dp, camera, sundir, 0, sundir,
                   trans_sun, sky_radiance_sun);
  double sun_radiance[NSSAMP] = {0};
  for (int i = 0; i < NSSAMP; ++i) {
    sun_radiance[i] = sky_radiance_sun[i] + trans_sun[i] * EXTSOL[i] / SOLOMG;
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
  int month, day;
  int num_threads = 1;
  int sorder = 4;
  int year = 0;
  double hour;
  int tsolar = 0;
  double hsm = 0.0;
  double ccover = 0.0;
  double aod = AOD0_CA;
  char *outname = "out";
  float angs_mie[NSSAMP] = {0};
  FVECT sundir;
  char *mie_ca_path = getpath("mie_ca.dat", getrlibpath(), R_OK);
  char *tau_path = "ssky_transmittance.dat";
  char *scat_path = "ssky_scattering.dat";
  char *scat1m_path = "ssky_delta_mie_scattering.dat";
  char *irrad_path = "ssky_irradiance.dat";
  if (argc < 4) {
    fprintf(stderr, "Usage: %s month day hour -y year -a lat -o lon -m tz\n",
            argv[0]);
    return 0;
  }
  if (argc == 4) {
    return 1;
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

  const double arctic_circle_latitude = 67;
  const double tropic_latitude = 23;
  const int summer_start_month = 4;
  const int summer_end_month = 9;

  // Determine if it's summer for the given hemisphere
  int is_northern_hemisphere = (s_latitude >= 0);
  int is_summer = (month >= summer_start_month && month <= summer_end_month);
  if (!is_northern_hemisphere) {
    is_summer = !is_summer;
  }

  Atmosphere atmos = {
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
      .cloud_density =
          {
              .layers =
                  {
                      {.width = 1000.0,
                       .exp_term = 0.0,
                       .exp_scale = 0.0,
                       .linear_term = 0.0,
                       .constant_term = 0.0},
                      {.width = 1100.0,
                       .exp_term = 0.0,
                       .exp_scale = 0.0,
                       .linear_term = 0.0,
                       .constant_term = 1.0},
                      {.width = 0.0,
                       .exp_term = 0.0,
                       .exp_scale = 0.0,
                       .linear_term = 0.0,
                       .constant_term = 0.0},
                  },

          },
      .beta_r0 = BR0_MS,
      .beta_c = BCLOUD,
      .beta_scale = aod / AOD0_CA,
      .cloud_cover = ccover,
      .beta_m = NULL,
  };

  // RayleighAtmos rayleigh_atmos;
  if (fabs(s_latitude) > arctic_circle_latitude) {
    if (is_summer) {
      atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_SS;
      atmos.beta_r0 = BR0_SS;
    } else {
      atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_SW;
      atmos.beta_r0 = BR0_SW;
    }
  } else if (fabs(s_latitude) > tropic_latitude) {
    if (is_summer) {
      atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_MS;
      atmos.beta_r0 = BR0_MS;
    } else {
      atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_MW;
      atmos.beta_r0 = BR0_MW;
    }
  } else {
    atmos.rayleigh_density.layers[0].exp_scale = -1.0 / HR_T;
    atmos.beta_r0 = BR0_T;
  }

  DATARRAY *mie_ca_dp = getdata(mie_ca_path);
  if (mie_ca_dp == NULL) {
    fprintf(stderr, "Error reading mie data\n");
    return 0;
  }
  atmos.beta_m = mie_ca_dp;

  if ((getpath(scat_path, getrlibpath(), R_OK)) == NULL) {
    printf("# Precomputing...\n");
    if (!precompute(sorder, &atmos, num_threads)) {
      printf("precompute failed\n");
      return 0;
    }
  }
  DATARRAY *tau_dp = getdata(tau_path);
  DATARRAY *scat_dp = getdata(scat_path);
  DATARRAY *scat1m_dp = getdata(scat1m_path);
  if (tau_dp == NULL || scat_dp == NULL || scat1m_dp == NULL) {
    fprintf(stderr, "Error reading transmittance data\n");
    return 0;
  }

  if (!gen_spect_sky(tau_dp, scat_dp, scat1m_dp, sundir, outname)) {
    fprintf(stderr, "gen_spect_sky failed\n");
    exit(1);
  }

  freedata(mie_ca_dp);
  freedata(tau_dp);
  freedata(scat_dp);
  freedata(scat1m_dp);

  return 1;
}
