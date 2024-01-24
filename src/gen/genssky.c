
#include "rtio.h"
#include "atmosphere.h"
#include "sun.h"


typedef struct {
    double aod;
    double ccover;
} Weather;


static int
positive_modulo(int i, int n) {
    return (i % n + n) % n;
}

int gen_spect_sky(const FVECT sundir, const Atmosphere *atmos)
{
    const int angres = 3;
    const int nthetas = 90 / angres + 1;
    const int nphis = 360 / angres + 1;
    char *skyfile = "skyspec.dat";
    if (remove(skyfile) == 0) {
        perror("overwrite skyspec.dat");
    }
    FILE *ssdat = fopen(skyfile, "a");
    fprintf(ssdat, "3\n0 90 %d\n0 360 %d\n380 780 20\n", nthetas, nphis);

    const unsigned int max_phi = 361;
    const unsigned int max_theta = 90;
    // 0 is south, -90 is east
    // double azimuthd = azimuth * 180.0 / M_PI;
    for (unsigned int t = 0; t <= max_theta; ++t) {
        for (int p = 0; p < max_phi; ++p) {
            SCOLOR skycolor;
            double phi = positive_modulo(270 - p, 360) / 180.0 * M_PI;
            double theta = t / 180.0 * M_PI;
            FVECT dir = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
            if (!tracesky(skycolor, dir, sundir, atmos)) {
                fprintf(stderr, "tracesky failed\n");
                exit(1);
            }
            for (int i = 0; i < 20; ++i) {
                fprintf(ssdat, "%f\n", skycolor[i]);
            }
            p += angres - 1;
        }
        t += angres - 1;
    }
    SCOLOR sunrad = {0};
    tracesun(sunrad, sundir, atmos);
    printf("void spectrum sunrad\n0\n0\n22 380 780 ");
    for (int i = 0; i < 20; ++i) {
        printf("%.1f ", sunrad[i]);
    }
    printf("\n\nsunrad light solar\n0\n0\n3 1 1 1\n\n");
    printf("solar source sun\n0\n0\n4 %f %f %f 0.533\n\n", sundir[0], sundir[1], sundir[2]);
    printf("void specdata skyfunc\n5 noop skyspec.dat . \"Acos(Dz)/DEGREE\" \"mod(atan2(-Dx, -Dy)/DEGREE,360)\"\n0\n0\n\n");
    printf("skyfunc glow sky_glow\n0\n0\n4 1 1 1 0\n\n");
    printf("sky_glow source sky\n0\n0\n4 0 0 1 180\n\n");
    fclose(ssdat);
    return 1;
}


static int parse_options(int argc, char *argv[], Datetime *dt, Location *lc, Weather *w) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s month day hour -y year -a lat -o lon -m tz\n", argv[0]);
        return 0;
    }
    int c;
    dt->month = atoi(argv[1]);
    dt->day = atoi(argv[2]);
    dt->hour = atof(argv[3]);
    if (argc == 4) {
        return 1;
    }
    for (int i=4; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'y':
                    dt->year = atoi(argv[i+1]);
                    break;
                case 'a':
                    lc->latitude = atof(argv[i+1]);
                    break;
                case 'o':
                    lc->longitude = atof(argv[i+1]);
                    break;
                case 'm':
                    lc->meridian = atof(argv[i+1]);
                    break;
                case 'd':
                    w->aod = atof(argv[i+1]);
                    break;
                case 'c':
                    w->ccover = atof(argv[i+1]);
                    break;
                default:
                    fprintf(stderr, "Unknown option %s\n", argv[i]);
                    exit(1);
            }
        }
    }
    return 1;
}

RayleighAtmos build_rayleigh_atmosphere(const Datetime *dt, const Location *lc) {
    const double arctic_circle_latitude = 67;
    const double tropic_latitude = 23;
    const int summer_start_month = 4;
    const int summer_end_month = 9;

    // Determine if it's summer for the given hemisphere
    int is_northern_hemisphere = (lc->latitude >= 0);
    int is_summer = (dt->month >= summer_start_month && dt->month <= summer_end_month);
    if (!is_northern_hemisphere) {
        is_summer = !is_summer;
    }

    if (fabs(lc->latitude) > arctic_circle_latitude) {
        if (is_summer) {
            return (RayleighAtmos) {HR_SS, BR0_SS};
        } else {
            return (RayleighAtmos) {HR_SW, BR0_SW};
        }
    } else if (fabs(lc->latitude) > tropic_latitude) {
        if (is_summer) {
            return (RayleighAtmos) {HR_MS, BR0_MS};
        } else {
            return (RayleighAtmos) {HR_MW, BR0_MW};
        }
    } else {
        return (RayleighAtmos) {HR_T, BR0_T};
    }
}


int main(int argc, char *argv[])
{
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

    if (!compute_sundir(sundir, &dt, &lc, 0)) {
        fprintf(stderr, "Cannot compute solar angle\n");
        exit(1);
    }
    if (!gen_spect_sky(sundir, &atmos)) {
        fprintf(stderr, "gen_spect_sky failed\n");
        exit(1);
    }
    return 1;
}
