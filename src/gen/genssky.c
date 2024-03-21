
#include "atmos.h"
#include "rtio.h"
#include "rtmath.h"



int gen_spect_sky(DATARRAY *tau_dp, DATARRAY *scat_dp, DATARRAY *scat1m_dp,
                  DATARRAY *irrad_dp, FVECT sundir) {
  const int angres = 3;
  const int nthetas = 90 / angres + 1;
  const int nphis = 360 / angres + 1;
  FVECT camera = {0, 0, ER + 1};
  char *skyfile = "skyspec2.dat";
  if (remove(skyfile) == 0) {
    perror("overwrite skyspec2.dat");
  }
  FILE *ssdat = fopen(skyfile, "w");
  fprintf(ssdat, "3\n0 90 %d\n0 360 %d\n380 780 20\n", nthetas, nphis);

  const unsigned int max_phi = 361;
  const unsigned int max_theta = 90;
  // 0 is south, -90 is east
  // double azimuthd = azimuth * 180.0 / M_PI;
  for (unsigned int t = 0; t <= max_theta; ++t) {
    for (int p = 0; p < max_phi; ++p) {
      SCOLOR skycolor;
      double phi = positive_modulo(270 - p, 360) / 180.0 * M_PI;
      double theta = t / 180.0 * PI;
      FVECT dir = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
      float trans[NSSAMP] = {0};
      float results[NSSAMP] = {0};
      get_sky_radiance(tau_dp, scat_dp, scat1m_dp, camera, dir, 0, sundir,
                       trans, results);
      for (int i = 0; i < NSSAMP; ++i) {
        fprintf(ssdat, "%f\n", results[i] * WVLSPAN);
      }
      p += angres - 1;
    }
    t += angres - 1;
  }
  fclose(ssdat);
  float trans_sun[NSSAMP] = {0};
  float sky_radiance_sun[NSSAMP] = {0};
  get_sky_radiance(tau_dp, scat_dp, scat1m_dp, camera, sundir, 0, sundir,
                       trans_sun, sky_radiance_sun);
  double sun_radiance[NSSAMP] = {0};
  for (int i = 0; i < NSSAMP; ++i) {
    sun_radiance[i] = sky_radiance_sun[i] + trans_sun[i] * EXTSOL[i] / SOLOMG;
  }
  printf("void spectrum sunrad\n0\n0\n22 380 780 ");
  for (int i = 0; i < 20; ++i) {
    printf("%.1f ", sun_radiance[i] * WVLSPAN);
  }
  printf("\n\nsunrad light solar\n0\n0\n3 1 1 1\n\n");
  printf("solar source sun\n0\n0\n4 %f %f %f 0.533\n\n", sundir[0], sundir[1],
         sundir[2]);
  printf("void specdata skyfunc\n5 noop %s . \"Acos(Dz)/DEGREE\" "
         "\"mod(atan2(-Dx, -Dy)/DEGREE,360)\"\n0\n0\n\n",
         skyfile);
  printf("skyfunc glow sky_glow\n0\n0\n4 1 1 1 0\n\n");
  printf("sky_glow source sky\n0\n0\n4 0 0 1 180\n\n");
  return 1;
}

static int parse_options(int argc, char *argv[], Datetime *dt, Location *lc,
                         Weather *w) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s month day hour -y year -a lat -o lon -m tz\n",
            argv[0]);
    return 0;
  }
  int c;
  dt->month = atoi(argv[1]);
  dt->day = atoi(argv[2]);
  dt->hour = atof(argv[3]);
  if (argc == 4) {
    return 1;
  }
  for (int i = 4; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'a':
        lc->latitude = atof(argv[i + 1]);
        break;
      // angstrom coefficient
      case 'b':
        w->aod = atof(argv[i + 1]);
        break;
      case 'c':
        w->ccover = atof(argv[i + 1]);
        break;
      case 'd':
        w->aod = atof(argv[i + 1]);
        break;
      case 'm':
        lc->meridian = atof(argv[i + 1]);
        break;
      case 'o':
        lc->longitude = atof(argv[i + 1]);
        break;
      case 'y':
        dt->year = atoi(argv[i + 1]);
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", argv[i]);
        exit(1);
      }
    }
  }
  return 1;
}

RayleighAtmos build_rayleigh_atmosphere(const Datetime *dt,
                                        const Location *lc) {
  const double arctic_circle_latitude = 67;
  const double tropic_latitude = 23;
  const int summer_start_month = 4;
  const int summer_end_month = 9;

  // Determine if it's summer for the given hemisphere
  int is_northern_hemisphere = (lc->latitude >= 0);
  int is_summer =
      (dt->month >= summer_start_month && dt->month <= summer_end_month);
  if (!is_northern_hemisphere) {
    is_summer = !is_summer;
  }

  if (fabs(lc->latitude) > arctic_circle_latitude) {
    if (is_summer) {
      return (RayleighAtmos){HR_SS, BR0_SS};
    } else {
      return (RayleighAtmos){HR_SW, BR0_SW};
    }
  } else if (fabs(lc->latitude) > tropic_latitude) {
    if (is_summer) {
      return (RayleighAtmos){HR_MS, BR0_MS};
    } else {
      return (RayleighAtmos){HR_MW, BR0_MW};
    }
  } else {
    return (RayleighAtmos){HR_T, BR0_T};
  }
}

int main(int argc, char *argv[]) {
  progname = "atmos_precompute";
  char *tau_path = "transmittance.dat";
  char *scat_path = "scattering.dat";
  char *scat1m_path = "delta_mie_scattering.dat";
  char *irrad_path = "irradiance.dat";
  if ((getpath(irrad_path, getrlibpath(), R_OK)) == NULL) {
    printf("precomputing...\n");
    if (!precompute(4)) {
      printf("precompute failed\n");
      return 1;
    }
  }
  DATARRAY *tau_dp = getdata(tau_path);
  DATARRAY *scat_dp = getdata(scat_path);
  DATARRAY *scat1m_dp = getdata(scat1m_path);
  DATARRAY *irrad_dp = getdata(irrad_path);
  double sun_irradiance[NSSAMP] = {0};
  Datetime dt = {0};
  Location lc = {0};
  Weather w = {0};
  FVECT sundir;
  if (!parse_options(argc, argv, &dt, &lc, &w)) {
    exit(1);
  }

  const RayleighAtmos rayleigh_atmos = build_rayleigh_atmosphere(&dt, &lc);

  // Always use the continental clean model, fix later.
  const MieAtmos mie_atmos = {HMS_CC, HMA_CC, BM0_CC, AM0_CC, w.aod};

  const Atmosphere atmos = {&rayleigh_atmos, &mie_atmos};

  if (!compute_sundir(&dt, &lc, 0, sundir)) {
    fprintf(stderr, "Cannot compute solar angle\n");
    exit(1);
  }

  if (!gen_spect_sky(tau_dp, scat_dp, scat1m_dp, irrad_dp, sundir)) {
    fprintf(stderr, "gen_spect_sky failed\n");
    exit(1);
  }
  return 1;
}
