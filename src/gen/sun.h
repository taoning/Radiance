/* RCSid $Id: sun.h,v 2.3 2019/11/07 23:15:07 greg Exp $ */
/*
 * Header file for solar position calculations
 */

#ifndef _RAD_SUN_H_
#define _RAD_SUN_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int year;
  int month;
  int day;
  float hour;
} Datetime;

typedef struct {
  double latitude;
  double longitude;
  double meridian;
} Location;
/* sun calculation constants */
extern double s_latitude;
extern double s_longitude;
extern double s_meridian;

extern int jdate(int month, int day);
extern double stadj(const int jd, const double s_meridian,
                    const double s_longitude);
extern double sdec(const int jd);
extern double salt(const double sd, const double st, const double s_latitude);
extern double sazi(const double sd, const double st, const double s_latitude);

extern double mjdate(const int year, const int month, const int day,
                     double hour, const double s_meridian);
extern double msdec(const double mjd, double *stp, const double s_longitude);

extern int compute_sundir(const Datetime *dt, const Location *loc,
                          const int tsolar, double sundir[3]);

#ifdef __cplusplus
}
#endif

#endif /* _RAD_SUN_H_ */
