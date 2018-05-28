/*
* SDR observation to ITU data conversion v180118
* allow user to fill a fixed value for reference delay
*/
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

void write_header(FILE *, char *, int, char *);
void quadfit(int, double *, double *, double *, double *, double *, int *, int *);
void date2doy(int , int, int, int *);
void doy2date(int, int, int *, int *);
int date2mjd(int, int, int);
void mjd2date(int, int *, int *, int *);
int doy2mjd(int, int);
void mjd2doy(int, int *, int *);

int main(int argc, char *argv[])
{
  FILE *fid, *fh;
  char prefix[80], path[80], filename[80], paramfile[80], data[25][25], str[200], *pch;
  char ch[10], loc[10], rem[10], li[10], ci[10], s[10], calr[10], esd[10], esig[10];
  int i, mjd, year, month, day, hr, sec, smp, atl, nobs, offset = 0, stmin, pn, rc;
  int hh, cnt;
  double hd, fm, gd, *y, *x, *w, res, rms, fc, fr = 3.0e+3, refd[2] = {0.0};
  // check if the input arguments are correct
  if (argc != 7 && argc != 8)
  {
    printf("input format: obs2ITU rx50_folder prefix year month day hour\n(e.g. ./obs2ITU /home/time/rx50 TWTL 2018 1 18 10)\n");
    return 1;
  }
  // check if the param file exists
  sscanf(argv[1], "%s", path);
  sscanf(argv[2], "%s", prefix);
  if (argc == 8) sscanf(argv[7], "%s", paramfile);
  else sprintf(paramfile, "%s/bin/ITU.param", path);
  if ((fh = fopen(paramfile, "r")) == NULL)
  {
    printf("no parameter file: %s\n", paramfile);
    return 1;
  }
  else
    fclose(fh);

  sscanf(argv[3], "%d", &year);
  sscanf(argv[4], "%d", &month);
  sscanf(argv[5], "%d", &day);
  sscanf(argv[6], "%d", &hr);
  mjd = date2mjd(year, month, day);

  // write header
  sprintf(filename, "%s/ITU/%s%2d.%03d", path, prefix, mjd / 1000, mjd % 1000);
  if ((fid = fopen(filename, "r")) == NULL)
  {
    fid = fopen(filename, "w");
    write_header(fid, paramfile, mjd, prefix);
  }
  fclose(fid);

  // read reference delays
  fh = fopen(paramfile, "r");
  while (fgets(str, 200, fh))
  {
    if (strncmp("[REFDELAY_ns_CH_A]", str, 18) == 0) fscanf(fh, "%lf", &refd[0]);
    if (strncmp("[REFDELAY_ns_CH_B]", str, 18) == 0) fscanf(fh, "%lf", &refd[1]);
  }
  fclose(fh);

  y = (double *)malloc(sizeof(double) * 3600);
  w = (double *)malloc(sizeof(double) * 3600);
  mjd2date(mjd, &year, &month, &day);

  // read measurements
  fh = fopen(paramfile, "r");
  while (fgets(str, 200, fh))
  {
    if (strncmp("[RX]", str, 4) == 0)
    {
      fgets(str, 200, fh);
      while (strncmp("[END]", str, 5) != 0)
      {
        i = 0;
			  pch = (char *)strtok(str, " ;\r\n");
			  if (pch != NULL)
				  strcpy(data[i++], pch);
			  while ((pch = (char *)strtok(NULL, " ;\r\n")) != NULL)
				  strcpy(data[i++], pch);
				if (i == 14 && (str[0] == 'A' || str[0] == 'B'))
				{
				  strcpy(ch, data[0]); // channel A or B
				  strcpy(loc, data[1]); // LOC station
				  strcpy(rem, data[2]); // REM station
				  sscanf(data[3], "%d", &stmin); // minute of STTIME
				  if (stmin < 0 || stmin > 59)
				  {
				    printf("error STTIME: %d\n", stmin);
				    continue;
				  }
				  sscanf(data[4], "%d", &pn); // PRN no.
				  if (pn < 0 || pn > 31)
				  {
				    printf("error PN : %d\n", pn);
				    continue;
				  }
				  sscanf(data[5], "%d", &rc); // Chip Rate (kcps)
				  if (rc != 1000 && rc != 2500)
				  {
				    printf("error chip rate: %d\n", rc);
				    continue;
				  }
				  sscanf(data[6], "%lf", &fm); // IF frequency (Hz)
				  sscanf(data[7], "%d", &nobs); // nominal track length
				  nobs += 1;
				  if (stmin + nobs / 60 > 60)
				  {
				    printf("error STTIME: %d not a full track\n", stmin);
				    continue;
				  }
				  x = (double *)malloc(sizeof(double) * nobs);
	        for (i = 0; i < nobs; i++)
		        x[i] = (double)i - ((double)(nobs - 1) / 2.0);
				  strcpy(li, data[8]);
				  strcpy(ci, data[9]);
				  strcpy(s, data[10]);
				  strcpy(calr, data[11]);
				  strcpy(esd, data[12]);
				  strcpy(esig, data[13]);

				  sprintf(filename, "%s/result/%4d%02d%02d%02d.ch%s.pn%02d.%04dkcps.dat", path, year, month, day, hr, ch, pn, rc);
	        if ((fid = fopen(filename, "r")) != NULL)
	        {
		        memset(y, 0, sizeof(double) * 3600);
		        memset(w, 0, sizeof(double) * 3600);
		        cnt = 0;
		        while (fgets(str, 200, fid))
		        {
			        i = 0;
			        pch = (char *)strtok(str, " ;\r\n");
			        if (pch != NULL)
				        strcpy(data[i++], pch);
			        while ((pch = (char *)strtok(NULL, " ;\r\n")) != NULL)
				        strcpy(data[i++], pch);
			        sscanf(data[0], "%d", &hh);
			        if (hh != hr)
			         continue;
              sscanf(data[3], "%lf", &fc);
			        sscanf(data[5], "%d", &sec);
			        sscanf(data[7], "%lf", &hd);
			        sscanf(data[8], "%lf", &gd);
			        i = sec % 3600 + offset;
			        if (fc > fm - fr && fc < fm + fr)
			        {
			          y[i] = hd * 1.0e+9 + gd;
			          w[i] = 1.0;
			          cnt++;
			        }
		        }
		        fclose(fid);
		        if (cnt * 2 > nobs)
		        {
			        quadfit(nobs, x, y + stmin * 60, w + stmin * 60, &res, &rms, &smp, &atl);
		          if (rms > 0.0 && rms < 10.0 && smp * 2 > nobs)
		          {
		            sprintf(filename, "%s/ITU/%s%2d.%03d", path, prefix, mjd / 1000, mjd % 1000);
			          fid = fopen(filename, "a");
			          fprintf(fid, "%6s %6s %2s %5d %02d%02d00 %3d", loc, rem, li, mjd, hr, stmin, nobs - 1);
			          fprintf(fid, "  %14.12lf %5.3lf %3d %3d", res * 1.0e-9, rms, smp, atl);
			          if (ch[0] == 'A' && refd[0] != 0.0)
			            fprintf(fid, " %+15.12lf 99999 %3s %1s %9s %9s %5s 999 999 9999\n", refd[0] * 1.0e-9, ci, s, calr, esd, esig);
			          else if (ch[0] == 'B' && refd[1] != 0.0)
			            fprintf(fid, " %+15.12lf 99999 %3s %1s %9s %9s %5s 999 999 9999\n", refd[1] * 1.0e-9, ci, s, calr, esd, esig);
			          else
			            fprintf(fid, " 999999999999999 99999 %3s %1s %9s %9s %5s 999 999 9999\n", ci, s, calr, esd, esig);
			          fclose(fid);
		          }
		        }
	        }
	        else
	          printf("no observation file: %s\n", filename);
	        free(x);
				}
        fgets(str, 200, fh);
      }
    }
  }
  fclose(fh);
	return 0;
}

void write_header(FILE *fid, char *paramfile, int mjd, char *prefix)
{
  FILE *fh;
  char str[200];
  fprintf(fid, "* %s%2d.%03d\n", prefix, mjd / 1000, mjd % 1000);
  fh = fopen(paramfile, "r");
  while (fgets(str, 200, fh))
  {
    if (strncmp("[HE]", str, 4) == 0)
    {
      fgets(str, 200, fh);
      while (strncmp("[END]", str, 5) != 0)
      {
        fprintf(fid, "%s", str);
        fgets(str, 200, fh);
      }
    }
  }
  fclose(fh);
	return;
}

void quadfit(int nobs, double *x, double *y, double *w, double *res, double *rms, int *cnt, int *atl)
{
	double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0, i = 0.0;
	double c0, c1, c2, det, temp, *z, mid, q1, q2, sigma3;
	int k, j, nobs2 = 0;
	*res = 0.0;
	*rms = 0.0;
	*cnt = 0;
	*atl = -1;
	z = (double *)malloc(sizeof(double) * nobs);
	for (k = 0; k < nobs; k++)
	{
		if (w[k] > 0.0)
		{
			z[nobs2] = y[k];
			nobs2++;
		}
	}
	if (2.0 * nobs2 > (double)nobs)
	{
		for (k = 0; k < nobs2 - 1; k++)
		{
			for (j = 0; j < nobs2 - k - 1; j++)
			{
				if (z[j] >= z[j + 1])
				{
					temp = z[j + 1];
					z[j + 1] = z[j];
					z[j] = temp;
				}
			}
		}
		mid = z[(int)((nobs2 - 1) / 2)];
		q1  = z[(int)((nobs2 - 1) / 4)];
		q2  = z[(int)(3 * (nobs2 - 1) / 4)];
		sigma3 = 3.0 * (q2 - q1) / 1.349;
		for (k = 0; k < nobs; k++)
		{
			if (w[k] > 0.0)
			{
				while (fabs(y[k] - mid) > 2.0e+6)
					(y[k] > mid) ? (y[k] -= 4.0e+6) : (y[k] += 4.0e+6);
				if (y[k] < mid + sigma3 && y[k] > mid - sigma3)
					y[k] -= mid;
				else
					w[k] = 0.0;
			}
		}
		// |a b c|   |c2|   |d|
		// |b c f| * |c1| = |e|
		// |c f i|   |c0|   |g|
		for (k = 0; k < nobs; k++)
		{
			if (w[k] > 0.0)
			{
				if (*atl < 0)
					nobs2 = k;
				*atl = k;
				a += x[k] * x[k] * x[k] * x[k];
				b += x[k] * x[k] * x[k];
				c += x[k] * x[k];
				f += x[k];
				i += 1.0;
				d += y[k] * x[k] * x[k];
				e += y[k] * x[k];
				g += y[k];
			}
		}
		*atl -= nobs2;
		det = a * (c * i - f * f) + b * (c * f - b * i) + c * (b * f - c * c);
		c2 = ((c * i - f * f) * d + (c * f - b * i) * e + (b * f - c * c) *g) / det;
		c1 = ((c * f - b * i) * d + (i * a - c * c) * e + (c * b - a * f) *g) / det;
		c0 = ((b * f - c * c) * d + (c * b - a * f) * e + (a * c - b * b) *g) / det;
		*res = c0 + mid + c1 * (mid * 1.0e-9) + c2 * (mid * 1.0e-9) * (mid * 1.0e-9);
		for (k = 0; k < nobs; k++)
			*rms += w[k] * (y[k] - c2 * x[k] * x[k] - c1 * x[k] - c0) * (y[k] - c2 * x[k] * x[k] - c1 * x[k] - c0);
		*rms = sqrt(*rms / i);
		*cnt = (int)i;
	}
	return;
}

void date2doy(int year, int month, int day, int *doy)
{
	int i, mday[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	if ((year % 4 == 0 && year % 100 != 0 ) || year % 400 == 0)
		mday[2] = 29;
  *doy = day;
  for (i = 1; i < month; i++)
  	*doy += mday[i];
	return;
}

void doy2date(int year, int doy, int *month, int *day)
{
	int mday[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	if ((year % 4 == 0 && year % 100 != 0 ) || year % 400 == 0)
		mday[2] = 29;
	*month = 1;
	*day = doy;
	while (*day > mday[*month])
  {
    *day -= mday[*month];
    (*month)++;
  }
	return;
}

/*
  date to MJD
*/
int date2mjd(int year, int month, int day)
{
  int mday[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int i, leap_days;
  leap_days = (year - 1) / 4 - (year - 1) / 100 + (year - 1) / 400;
  if ((year % 4 == 0 && year % 100 != 0 ) || year % 400 == 0)
    mday[2] = 29;
  for (i = 1; i < month; i++)
    day = day + mday[i];
  return (year - 1) * 365 + day + leap_days - 678576;
}
/*
  MJD to date
*/
void mjd2date(int mjd, int *year, int *month, int *day)
{
  int mday [13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  // 400 years = 146097 days, 100 years = 36524 days, 4 years = 1461 days, 3 normal years = 1095
  *year = 1; *month = 1;
  *day = mjd + 678576; //total days from 0001/01/01
  *year += *day / 146097 * 400; *day = *day % 146097;
  *year += *day /  36524 * 100; *day = *day %  36524;
  *year += *day /   1461 *   4; *day = *day %   1461;

  if (*day <= 1095)
  {
  	*year += *day / 365;
  	*day = *day % 365;
  }
  else
  {
  	*year += 3;
  	*day -= 1095;
  }

  if (*day == 0)
  {
  	(*year)--;
  	*month = 12;
  	*day = 31;
  	return;
  }

  if (( *year % 4 == 0 && *year % 100 != 0 ) || *year % 400 == 0)
  	mday[2] = 29;

  while (*day > mday[*month])
  {
    *day -= mday[*month];
    (*month)++;
  }
  return;
}
/*
  DOY (day of year) to MJD
*/
int doy2mjd(int year, int doy)
{
  int leap_days;
  leap_days = (year - 1) / 4 - (year - 1) / 100 + (year - 1) / 400;
  return (year - 1) * 365 + doy + leap_days - 678576;
}
/*
  MJD to DOY (day of year)
*/
void mjd2doy(int mjd, int *year, int *doy)
{
	// 400 years = 146097 days, 100 years = 36524 days, 4 years = 1461 days
	*year = 1;
	*doy = mjd + 678576;
	*year += *doy / 146097 * 400; *doy = *doy % 146097;
	*year += *doy /  36524 * 100; *doy = *doy %  36524;
	*year += *doy /   1461 *   4; *doy = *doy %   1461;

	if (*doy <= 1095)
  {
  	*year += *doy / 365;
  	*doy = *doy % 365;
  }
  else
  {
  	*year += 3;
  	*doy -= 1095;
  }

	if (*doy == 0)
	{
		(*year)--;
		if (( *year % 4 == 0 && *year % 100 != 0 ) || *year % 400 == 0)
			*doy = 366;
		else
			*doy = 365;
	}

	return;
}
