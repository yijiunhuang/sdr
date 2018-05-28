/*
* The software-defined radio receiver V170215
* The authors would thank NICT for original source codes,
* thanks NICT, KRISS, NTSC, PTB, OP, NIST, VNIIFTRI, INRIM, AOS, and NIM for
* their technical support and feedback, and
* thank BIPM and CCTF WG on TWSTFT for coordinating the interational activities.
*/

#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_fit.h>
#include <pthread.h>
#include <cublas_v2.h>
#include "sample_v170215.hpp"

typedef struct _channel_info
{
  bool is_first; // 0x1 = first track; 0x0 = not first track
  bool is_trk;   // 0x1 = being tracked; 0x0 = no signal
  bool is_sic;   // 0x1 = successive interference cancellation (SIC); 0x0 = normal
  bool is_chA;   // 0x1 = physical ch A; 0x0 = ch B
  int prnno;     // PRN number
  int rc;        // number of chip per second (chips)
  int pt;        // coarse code head position (samples)
  int clen;      // number of chip per code length
  int *code;        // PRN code sequence {-1, 0, 1}
  int *cuda_code;   // PRN code sequence {-1, 0, 1} on GPU memory
  cufftDoubleComplex *cuda_prn;      // freq domain code samples
  cufftDoubleComplex *cuda_prn_t;    // time domain code waveform
  cufftDoubleComplex *cuda_prn_acq;  // freq domain code waveform (ACQ)
  double fc_start; // central freq. from parameter file
  double fc;   // carrier freq. (Hz)
  double df;
  double gd;   // code phase (ns)
  double dg;   // code phase rate (ns/s)
  double phi;  // carrier phase (cycle)
  double last_phi; // last carrier phase (cycle)
  double peak; // peak value
  double prn_power; // power of the reference signal (V^2)
  double fmin; // lower (negative) cut-off freq (Hz) of the low-pass filter
  double fmax; // upper (positive) cut-off freq (Hz) of the low-pass filter
  double range; // acquisition frequency range (Hz)
  double step;  // acquisition frequency step (Hz)
  double snr_min; // minimum required SNR (dB)
  int n_stop;
  int n_start;
} channel_info;

typedef struct _fx_result
{
  FILE *fout;
  int count;        // available code arrival time measurements per second
  int count_sic;    // available code arrival time measurements, SIC, per second
  int *pidx;
  double *ttag_gd;  // time reference of the code arrival time
  double *gd;       // code arrival time (ns) on ttag_gd
  double *gd_sic;   // code arrival time, SIC (ns) on ttag_gd
  double *ttag_phi; // time reference of the carrier phase
  double *phi;      // carrier phase (cycle) on ttag_phi
  double *phi_raw;  // carrier phase (cycle) on the reference time
  double *w;        // indicator for code, 1 = available; 0 = no use
  double *w_sic;    // indicator for code, SIC, 1 = available; 0 = no use
  double *amp;      // signal amplitude (V)
  double *signal;   // signal power (V^2)
  double *noise;    // noise power (V^2)
} fx_result;

typedef struct _system_info
{
  int nch;      // number of channels
  int ntrk;     // number of channels in tracking state
  int nsic;     // number of channels for SIC
  int nobs;     // number of samples per code length
  int sps;      // samples per second
  int portion;  // number of codes per second
  double fs;    // sampling frequency (Hz)
  double period;// code length (s)
} system_info;

/* function prototype */
int CAcode(channel_info *);
int NICTcode(channel_info *);
int SATREcode(channel_info *);
double average(int, double *, double *);
double kth_smallest(double *, int, int);

void current_time(int, int *, int *, int *, int *, int *, int *, int *, int *);
void date2doy(int , int, int, int *);
void doy2date(int, int, int *, int *);
int date2mjd(int, int, int);
void mjd2date(int, int *, int *, int *);
int doy2mjd(int, int);
void mjd2doy(int, int *, int *);

__global__ void SIC(double *, int, int, int, int, double, cufftDoubleComplex *, int, double, double);
__global__ void PRN_sampling(int, int *, cufftDoubleComplex *, int, double, int, double);
__global__ void down_conversion(int, double, double, double *, cufftDoubleComplex *);
__global__ void down_conversion2(int, double, double, char *, cufftDoubleComplex *, bool);
__global__ void cross_spectrum(int, cufftDoubleComplex *, cufftDoubleComplex *, double, double, double, cufftDoubleComplex *);
__global__ void binary_to_waveform(int, char *, double *, bool);
__global__ void conv(int, int, cufftDoubleComplex *, double *, double *);

/* MAIN PROCESS */
int main(int argc, char *argv[])
{
  // check input arguments
  if (argc != 1 && argc != 7)
  {
    printf("use ./rx50 or ./rx50 yyyy MM dd hh mm ss\n");
    return 1;
  }

  FILE *fparam, *flog;
  samplebuffer sb;
  pthread_t tsampling;
  timetag *ttag;
  system_info sinfo; sinfo.sps = 50000000; sinfo.nch = 0;
  channel_info *cinfo;
  fx_result *fxres;
  cufftHandle plan, plan_trk, plan_acq, plan_sic;
  cufftDoubleComplex *cuda_xcor_fx, *cuda_robs, *cuda_robs_acq;

  cudaSetDevice(0); // identify the GPU card for computation
  cublasHandle_t handle;
  cublasCreate(&handle);

  char str[200], str2[200], result_dir[80], filename[80], paramfile[80], prn_type[10], ch[2], id[2], *pch;
  double c0, c1, chisq, peak, c00, c01, c11;
  double phi, x1, x2, x3, x4, x5, rx_power, *rpow;
  double fcc, fbest, frange, fstep, flow, fhigh;
  double global_peak, init_peak, *cuda_ddbs;
  int i, k, p, peak_idx, blocks, global_peak_idx, cnt;
  char *cuda_samples;
  double *cuda_xcor_conv, *xcor_conv, *cuda_xcor_phi, *xcor_phi;
  int ii, jj, idx, pn, rc, imax;
  double *res, tmp, snr_min, fc_start, ftrmax;
  bool is_chA, is_sic;
  struct tm *timeinfo;
  timeinfo = (struct tm *)malloc(sizeof(struct tm));
  
  // adjustable variables
  int nobs = 20; // maximum delay spread (+-samps)
  int nch_max = 8; // maximum number of software channels

  srand(time(NULL)); // random pick-up for ACQ
  sprintf(prn_type, "SATRE");
  sprintf(result_dir, "../result");
  sprintf(paramfile, "satre.param");

  // assign system parameters
  if ((fparam = fopen(paramfile, "r")) == NULL)
  {
    printf("no such parameter file : %s\n", paramfile);
    return 1;
  }
  fgets(str, 200, fparam); // discard the first raw

  // assign channel parameters and function result
  ttag = (timetag *)malloc(sizeof(timetag));
  cinfo = (channel_info *)malloc(sizeof(channel_info) * nch_max);
  fxres = (fx_result *)malloc(sizeof(fx_result) * nch_max);

  // modify parameters
  if (strcmp(prn_type, "GNSS") == 0)
    sinfo.period = 0.001; // 1ms
  else if (strcmp(prn_type, "NICT") == 0)
    sinfo.period = 0.002; // 2ms
  else if (strcmp(prn_type, "SATRE") == 0)
    sinfo.period = 0.004; // 4ms
  else
  {
    printf("PRN type error: %s", prn_type);
    return 1;
  }
  sinfo.nobs = (int)round((double)sinfo.sps * sinfo.period);
  sinfo.portion = (int)round(1.0 / sinfo.period);
  sinfo.fs = (double)sinfo.sps;
  cufftPlan1d(&plan, sinfo.nobs, CUFFT_Z2Z, 1);
  cufftPlan1d(&plan_acq, sinfo.nobs * 2, CUFFT_Z2Z, 1);

  //printf("sampling frequency : %.0lf Hz\n", sinfo.fs);
  //printf("code period : %lf sec\n", sinfo.period);
  //printf("nobs : %d\n", sinfo.nobs);
  //printf("portions : %d\n", sinfo.portion);

  // assign memories
  for (i = 0; i < nch_max; i++)
  {
    cinfo[i].prnno = 0;
    cinfo[i].fc_start = 0.0;
    cinfo[i].rc = 0;
    cinfo[i].fmax = 0.0;
    cinfo[i].fmin = 0.0;
    cinfo[i].is_chA = 0x0;
    cinfo[i].is_sic = 0x0;
    cinfo[i].range = 1.0;
    cinfo[i].step = 1.0;
    cinfo[i].snr_min = 0.0;
    cinfo[i].is_trk = 0x0;
    cinfo[i].is_first = 0x0;
    cinfo[i].fc = 0.0;
    cinfo[i].pt = 0;
    cinfo[i].gd = 0.0;
    cinfo[i].dg = 0.0;
    cinfo[i].phi = 0.0;
    cinfo[i].n_stop = 0;
    cinfo[i].n_start = 0;
    cinfo[i].clen = 0;
    cinfo[i].peak = 0.0;
    cinfo[i].last_phi = 0.0;
    cudaMalloc((void **)&cinfo[i].cuda_prn, sizeof(cufftDoubleComplex) * sinfo.nobs);
    cudaMalloc((void **)&cinfo[i].cuda_prn_t, sizeof(cufftDoubleComplex) * sinfo.nobs);
    cudaMalloc((void **)&cinfo[i].cuda_prn_acq, sizeof(cufftDoubleComplex) * sinfo.nobs * 2);
    fxres[i].pidx     = (int *)malloc(sizeof(int) * sinfo.portion);
    fxres[i].gd       = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].gd_sic   = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].ttag_gd  = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].phi      = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].ttag_phi = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].signal   = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].noise    = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].w        = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].w_sic    = (double *)malloc(sizeof(double) * sinfo.portion);
    fxres[i].amp      = (double *)malloc(sizeof(double) * sinfo.portion);
  }

  // initialize the buffer of samples
  for (i = 0; i < MAXBUF; i++)
  {
    sb.ttag[i] = (timetag *)malloc(sizeof(timetag));
    memset(sb.ttag[i], 0x0, sizeof(timetag));
    sb.buf[i] = (char *)malloc(sizeof(char) * sps * 2);
    if (sb.buf[i] == NULL)
    {
      printf("sample buffer error\n");
      return 1;
    }
    memset(sb.buf[i], 0, sizeof(char) * sps * 2);
  }

  // set time stamp
  if (argc == 1) sb.is_tt = 0x0;
  else sb.is_tt = 0x1;
  if (sb.is_tt == 0x1)
  {
    sscanf(argv[1], "%d", &timeinfo->tm_year);
    timeinfo->tm_year -= 1900;
    sscanf(argv[2], "%d", &timeinfo->tm_mon);
    timeinfo->tm_mon -= 1;
    sscanf(argv[3], "%d", &timeinfo->tm_mday);
    sscanf(argv[4], "%d", &timeinfo->tm_hour);
    sscanf(argv[5], "%d", &timeinfo->tm_min);
    sscanf(argv[6], "%d", &timeinfo->tm_sec);
    sb.tt = mktime(timeinfo);
  }
  free(timeinfo);
  sb.readd = -1;
  sb.writed = -1;

  cudaMalloc((void **)&cuda_robs_acq, sizeof(cufftDoubleComplex) * sinfo.nobs * 2);
  cudaMalloc((void **)&cuda_samples, sizeof(char) * sinfo.sps * 2);

  rpow = (double *)malloc(sizeof(double) * nch_max);
  cudaMalloc((void **)&cuda_robs, sizeof(cufftDoubleComplex) * sinfo.nobs * nch_max);
  cudaMalloc((void **)&cuda_xcor_fx, sizeof(cufftDoubleComplex) * sinfo.nobs * nch_max);
  cudaMalloc((void **)&cuda_ddbs, sizeof(double) * sinfo.nobs);

  cudaMalloc((void **)&cuda_xcor_conv, sizeof(double) * (nobs * 2 + 1));
  cudaMalloc((void **)&cuda_xcor_phi, sizeof(double) * (nobs * 2 + 1));
  xcor_conv = (double *)malloc(sizeof(double) * (nobs * 2 + 1));
  xcor_phi = (double *)malloc(sizeof(double) * (nobs * 2 + 1));
  res = (double *)malloc(sizeof(double) * sinfo.portion);

  // create a thread and start sampling
  pthread_mutex_init(&sb.mymutex, NULL);
  pthread_cond_init(&sb.mycond, NULL);
  pthread_create(&tsampling, NULL, sampling, (void *)&sb);

////////////////////////////// infinite observation //////////////////////////

  while (1)
  {
    // Get samples from buffer
    pthread_mutex_lock(&sb.mymutex);
    sb.readd++;
    if (sb.readd > sb.writed || sb.writed == -1)
      pthread_cond_wait(&sb.mycond, &sb.mymutex); // wait until the first sample
    cudaMemcpy(cuda_samples, sb.buf[sb.readd % MAXBUF], sizeof(char) * sps * 2, cudaMemcpyHostToDevice);
    memcpy(ttag, sb.ttag[sb.readd % MAXBUF], sizeof(timetag));
    pthread_mutex_unlock(&sb.mymutex);

    // Start: set new parameters and refresh every 30 seconds
    if (sinfo.nch == 0 || ttag->second == 17 || ttag->second == 47)
    {
      if ((fparam = fopen(paramfile, "r")) != NULL)
      {
        i = 0;
        while (fgets(str, 200, fparam))
        {
          if (str[0] == '#') continue;
          strcpy(str2, str);
          k = 0;
			    pch = (char *)strtok(str2, " ;\r\n");
			    if (pch != NULL) k++;
			    while ((pch = (char *)strtok(NULL, " ;\r\n")) != NULL) k++;
          if ((str[0] == 'A' || str[0] == 'B') && (str[2] == 'N' || str[2] == 'S') && k == 9)
          {
            if (str[0] == 'A') is_chA = 0x1; // physical channel A
            if (str[0] == 'B') is_chA = 0x0; // physical channel B
            if (str[2] == 'S') is_sic = 0x1; // with SIC
            if (str[2] == 'N') is_sic = 0x0; // normal
            sscanf(str, "%s %s %d %lf %d %lf %lf %lf %lf", ch, id, &pn, &fc_start, &rc, &ftrmax, &frange, &fstep, &snr_min);

            if (cinfo[i].is_chA == is_chA && cinfo[i].is_sic == is_sic && cinfo[i].prnno == pn && cinfo[i].fc_start == fc_start && cinfo[i].rc == rc * 1000 && cinfo[i].fmax == ftrmax * 1.0e+3 && cinfo[i].range >= frange && cinfo[i].range < frange * 2.0 && cinfo[i].step >= fstep && cinfo[i].step < fstep * 2.0 && cinfo[i].snr_min == snr_min)
            {
              i++;
              if (i < nch_max) continue;
              else break;
            }
            else if (pn >= 0 && pn <= 31 && fc_start > 1.0e+7 && fc_start < 1.0e+8 && (rc == 2500 || rc == 1000) && ftrmax > 0.0 && ftrmax < 1.0e+5 && frange > 0.0 && frange < 1.0e+6 && frange > fstep && snr_min > -100.0)
            {
              cinfo[i].prnno = pn;
              cinfo[i].fc_start = fc_start;
              cinfo[i].rc = rc * 1000;
              cinfo[i].fmax = ftrmax * 1.0e+3;
              cinfo[i].fmin = -cinfo[i].fmax;
              cinfo[i].is_chA = is_chA;
              cinfo[i].is_sic = is_sic;
              cinfo[i].range = 1.0;
              while (cinfo[i].range < frange) cinfo[i].range *= 2.0;
              cinfo[i].step = 1.0;
              while (cinfo[i].step < fstep) cinfo[i].step *= 2.0;
              cinfo[i].snr_min = snr_min;
              cinfo[i].is_trk = 0x0;
              cinfo[i].is_first = 0x1;
              cinfo[i].fc = 0.0;
              cinfo[i].pt = 0;
              cinfo[i].gd = 0.0;
              cinfo[i].dg = 0.0;
              cinfo[i].phi = 0.0;
              cinfo[i].n_stop = (int)floor(cinfo[i].fmax * sinfo.period);
              cinfo[i].n_start = (int)floor(cinfo[i].fmin * sinfo.period);
              cinfo[i].clen = 0;
              cinfo[i].peak = 0.0;
              cinfo[i].last_phi = 0.0;
    
              // obtain code length
              if (strcmp(prn_type, "GNSS") == 0) cinfo[i].clen = 1023; // GNSS
              else if (strcmp(prn_type, "NICT") == 0) cinfo[i].clen = 4095; // NICT modem
              else if (strcmp(prn_type, "SATRE") == 0)
              {
                if (cinfo[i].rc == 2500000) cinfo[i].clen = 10000; // SATRE modem
                else if (cinfo[i].rc == 1000000) cinfo[i].clen = 4000; // SATRE modem
                else
                {
                  printf("code rate error (new parameter): %d\n", cinfo[i].rc);
                  i++;
                  if (i < nch_max) continue;
                  else break;
                }
              }

              // assign memory for PRN code
              cinfo[i].code = (int *)malloc(sizeof(int) * cinfo[i].clen);
              cudaMalloc((void **)&cinfo[i].cuda_code, sizeof(int) * cinfo[i].clen);
              if (strcmp(prn_type, "GNSS") == 0) CAcode(&cinfo[i]); // GNSS
              else if (strcmp(prn_type, "NICT") == 0) NICTcode(&cinfo[i]); // NICT modem
              else if (strcmp(prn_type, "SATRE") == 0) SATREcode(&cinfo[i]); // SATRE modem

              cudaMemcpy(cinfo[i].cuda_code, cinfo[i].code, sizeof(int) * cinfo[i].clen, cudaMemcpyHostToDevice);
              // assign GPU memory for PRN samples and spectrum
              cudaMemset(cinfo[i].cuda_prn_acq, 0x0, sizeof(cufftDoubleComplex) * sinfo.nobs * 2);
              blocks = sinfo.nobs / 1000;
              PRN_sampling<<<blocks, 1000>>>(sinfo.nobs, cinfo[i].cuda_code, cinfo[i].cuda_prn, cinfo[i].rc, sinfo.fs, cinfo[i].clen, 0.0); // cuda_prn = +- 1
              cudaMemcpy(cinfo[i].cuda_prn_t, cinfo[i].cuda_prn, sizeof(cufftDoubleComplex) * sinfo.nobs, cudaMemcpyDeviceToDevice); // cuda_prn_t = +- 1
              cudaMemcpy(cinfo[i].cuda_prn_acq, cinfo[i].cuda_prn, sizeof(cufftDoubleComplex) * sinfo.nobs, cudaMemcpyDeviceToDevice);
              cufftExecZ2Z(plan, cinfo[i].cuda_prn, cinfo[i].cuda_prn, CUFFT_FORWARD);
              cufftExecZ2Z(plan_acq, cinfo[i].cuda_prn_acq, cinfo[i].cuda_prn_acq, CUFFT_FORWARD);

              // compute reference signal power Pc (baseband filtered PRN)
              cublasDznrm2(handle, abs(cinfo[i].n_start), (cuDoubleComplex *)cinfo[i].cuda_prn + sinfo.nobs - abs(cinfo[i].n_start), 1, &peak);
              cublasDznrm2(handle, abs(cinfo[i].n_stop), (cuDoubleComplex *)cinfo[i].cuda_prn, 1, &rx_power);
              cinfo[i].prn_power = pow(peak / (double)sinfo.nobs, 2) + pow(rx_power / (double)sinfo.nobs, 2); // power in V^2

              //printf("#%02d, reference signal power: %lf dBm\n", cinfo[i].prnno, 10.0 * log10(cinfo[i].prn_power * 1000.0 / 50.0));

              // initialize memories
              memset(fxres[i].pidx    , 0x0, sizeof(int)    * sinfo.portion);
              memset(fxres[i].gd      , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].gd_sic  , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].ttag_gd , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].phi     , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].ttag_phi, 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].signal  , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].noise   , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].w       , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].w_sic   , 0x0, sizeof(double) * sinfo.portion);
              memset(fxres[i].amp     , 0x0, sizeof(double) * sinfo.portion);

              // LOG
              sprintf(filename, "%s/%4d%02d%02d%02d.log", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour);
              flog = fopen(filename, "a");
              fprintf(flog, "%02d:%02d:%02d + %11.9lf set param: #%02d %8.0lf %4d %5.0lf %5.0lf %5.0lf %3.0lf\n", ttag->hour, ttag->minute, ttag->second, ttag->fsec, cinfo[i].prnno, cinfo[i].fc_start, cinfo[i].rc / 1000, cinfo[i].fmax * 1.0e-3, cinfo[i].range, cinfo[i].step, cinfo[i].snr_min);
              fclose(flog);
              i++;
              if (i < nch_max) continue;
              else break;
            }
          }
        }
        fclose(fparam);
        sinfo.nch = i;
      }
    } // End: set new parameters and refresh every 30 seconds

    // Start: two-stage acquisition
    for (p = 0; p < 2; p++)
    {
      for (i = 0; i < sinfo.nch; i++)
      {
        if ((p == 0 && cinfo[i].is_trk == 0x0) || (p == 1 && cinfo[i].is_trk == 0x0 && cinfo[i].fc != 0.0))
        {
          if (p == 0) 
          {
            fbest = cinfo[i].fc_start;
            frange = cinfo[i].range;
            fstep = cinfo[i].step;
          }
          else
          {
            fbest = cinfo[i].fc;
            frange = cinfo[i].range / 2.0;
            fstep = cinfo[i].step / 2.0;
          }
          init_peak = 0.0;
          rx_power = 0.0;
          idx =  2 * (rand() % (sinfo.portion - 1)) * sinfo.nobs;
          while (1)
          {
            flow = fbest - frange;
            fhigh = fbest + frange;
            for (fcc = flow; fcc <= fhigh; fcc += fstep)
            {
              blocks = sinfo.nobs * 2 / 1000;
              // down convert
              down_conversion2<<<blocks, 1000>>>(sinfo.nobs * 2, fcc / sinfo.fs, 0.0, cuda_samples + idx, cuda_robs_acq, cinfo[i].is_chA);
              cufftExecZ2Z(plan_acq, cuda_robs_acq, cuda_robs_acq, CUFFT_FORWARD);

              // compute reception power (Px) in the last search period
              if (fstep >= 1.0 && fstep < 2.0 && rx_power == 0.0)
              {
                cublasDznrm2(handle, abs(cinfo[i].n_start * 2), (cuDoubleComplex *)cuda_robs_acq + 2 * sinfo.nobs - abs(cinfo[i].n_start * 2), 1, &peak);
                cublasDznrm2(handle, abs(cinfo[i].n_stop * 2) , (cuDoubleComplex *)cuda_robs_acq, 1, &rx_power);
                rx_power = pow(peak / 2.0 / (double)sinfo.nobs, 2) + pow(rx_power / 2.0 / (double)sinfo.nobs, 2); // reception power (V^2)
              }

              // perform cross correlation
              cross_spectrum<<<blocks, 1000>>>(sinfo.nobs * 2, cuda_robs_acq, cinfo[i].cuda_prn_acq, sinfo.fs / (double)sinfo.nobs / 2.0, cinfo[i].fmax, cinfo[i].fmin, cuda_robs_acq);
              cufftExecZ2Z(plan_acq, cuda_robs_acq, cuda_robs_acq, CUFFT_INVERSE);

              // find peak
              cublasIzamax(handle, sinfo.nobs, (cuDoubleComplex *)cuda_robs_acq, 1, &peak_idx);
              peak_idx -= 1; // cublasIzamax() returns [1:nfft], but we want [0:nfft-1]
              cublasDznrm2(handle, 1, (cuDoubleComplex *)cuda_robs_acq + peak_idx, 1, &peak);

              if (init_peak == 0.0)
              {
                init_peak = peak;
                global_peak = peak;
              }
              if (peak > global_peak)
              {
                global_peak = peak;
                fbest = fcc;
                global_peak_idx = peak_idx;
              }
            }
            // update to finer search range
            frange = fstep;
            fstep = fstep / 2.0;
            if (fstep < 1.0)
              break;
          }
          global_peak = 8.0 * global_peak * global_peak / cinfo[i].prn_power; // signal power (V^2)
          snr_min = pow(10.0, cinfo[i].snr_min / 10.0);
          if ((1.0 + snr_min) * global_peak > snr_min * rx_power) // say signal exists if SNR > required SNR
          {
            if (p == 0)
            {
              cinfo[i].fc = fbest;
              cinfo[i].gd = (double)global_peak_idx * 1.0e+9 / sinfo.fs;
              cinfo[i].pt = global_peak_idx;
              cinfo[i].peak = global_peak;
            }
            else
            {
              cinfo[i].fc = floor((fbest + cinfo[i].fc) / 2.0);
              cinfo[i].gd = (double)(global_peak_idx + cinfo[i].pt) * 1.0e+9 / sinfo.fs / 2.0;
              cinfo[i].pt = (global_peak_idx + cinfo[i].pt) / 2;
              cinfo[i].peak = (cinfo[i].peak + global_peak) / 2.0;
              cinfo[i].is_trk = 0x1;
            }
            // LOG
            sprintf(filename, "%s/%4d%02d%02d%02d.log", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour);
            flog = fopen(filename, "a");
            if (p == 0) fprintf(flog, "%02d:%02d:%02d + %11.9lf ACQ1: PRN#%2d: %3d %8.0lf %7.0lf %6d %8.3lf %8.3lf\n", ttag->hour, ttag->minute, ttag->second, ttag->fsec, cinfo[i].prnno, idx / 2 / sinfo.nobs, cinfo[i].fc, cinfo[i].gd, cinfo[i].pt, 10.0 * log10(cinfo[i].peak * 1000.0 / 50.0), 10.0 * log10(rx_power * 1000.0 / 50.0));
            else
            {
              if (cinfo[i].is_chA == 0x1) printf("A: #%02d %4.1lf Mcps Code Lock\n", cinfo[i].prnno, (double)cinfo[i].rc * 1.0e-6);
              else printf("B: #%02d %4.1lf Mcps Code Lock\n", cinfo[i].prnno, (double)cinfo[i].rc * 1.0e-6);
              fprintf(flog, "%02d:%02d:%02d + %11.9lf ACQ2: PRN#%2d: %3d %8.0lf %7.0lf %6d %8.3lf %8.3lf\n", ttag->hour, ttag->minute, ttag->second, ttag->fsec, cinfo[i].prnno, idx / 2 / sinfo.nobs, cinfo[i].fc, cinfo[i].gd, cinfo[i].pt, 10.0 * log10(cinfo[i].peak * 1000.0 / 50.0), 10.0 * log10(rx_power * 1000.0 / 50.0));
            }
            fclose(flog);
            // LOG
          }
          else cinfo[i].fc = 0.0;
        }
      }
    } // End: two-stage acquisition

    sinfo.ntrk = 0;
    sinfo.nsic = 0;
    for (i = 0; i < sinfo.nch; i++)
    {
      if (cinfo[i].is_trk == 0x1)
      {
      	sinfo.ntrk++;
      	if (cinfo[i].is_sic == 0x1)
      	  sinfo.nsic++;
      }
    }
    if (sinfo.ntrk == 0) // if no desired PN code here, re-do ACQ for the next samples
    {
      printf("acquisition: %d , no signal\n", sinfo.nch);
      continue;
    }

    ttag->mjd = date2mjd(ttag->year, ttag->month, ttag->day);
    cufftPlan1d(&plan_trk, sinfo.nobs, CUFFT_Z2Z, sinfo.ntrk);
    cufftPlan1d(&plan_sic, sinfo.nobs, CUFFT_Z2Z, sinfo.nsic);
    for (i = 0; i < sinfo.nch; i++)
    {
      memset(fxres[i].w, 0x0, sizeof(double) * sinfo.portion);
      memset(fxres[i].w_sic, 0x0, sizeof(double) * sinfo.portion);
      fxres[i].count = 0;
      fxres[i].count_sic = 0;
    }

    // Start: Tracking
    for (p = 0; p < sinfo.portion - 1; p++)
    {
      // Start: measure code phase, amplitude, and dphi
      cnt = 0;
      blocks = sinfo.nobs / 1000;
      for (i = 0; i < sinfo.nch; i++)
      {
        if (cinfo[i].is_trk == 0x1)
        {
          phi = fmod(cinfo[i].fc * (double)(p * sinfo.nobs + cinfo[i].pt) / sinfo.fs, 1.0);
          idx = 2 * (p * sinfo.nobs + cinfo[i].pt);
          down_conversion2<<<blocks, 1000>>>(sinfo.nobs, cinfo[i].fc / sinfo.fs, phi, cuda_samples + idx, cuda_robs + cnt * sinfo.nobs, cinfo[i].is_chA);
          cnt++;
        }
      }
      cufftExecZ2Z(plan_trk, cuda_robs, cuda_robs, CUFFT_FORWARD);

      // compute reception power (Px)
      cnt = 0;
      for (i = 0; i < sinfo.nch; i++)
      {
        if (cinfo[i].is_trk == 0x1)
        {
          cublasDznrm2(handle, abs(cinfo[i].n_start), (cuDoubleComplex *)cuda_robs + cnt * sinfo.nobs + sinfo.nobs - abs(cinfo[i].n_start), 1, &peak);
          cublasDznrm2(handle, abs(cinfo[i].n_stop), (cuDoubleComplex *)cuda_robs + cnt * sinfo.nobs, 1, &rx_power);
          rpow[i] = pow(rx_power / (double)sinfo.nobs, 2) +  pow(peak / (double)sinfo.nobs, 2);
          cnt++;
        }
      }

      // perform cross correlation
      cnt = 0;
      for (i = 0; i < sinfo.nch; i++)
      {
        if (cinfo[i].is_trk == 0x1)
        {
          cross_spectrum<<<blocks, 1000>>>(sinfo.nobs, cuda_robs + cnt * sinfo.nobs, cinfo[i].cuda_prn, sinfo.fs / (double)sinfo.nobs, cinfo[i].fmax, cinfo[i].fmin, cuda_xcor_fx + cnt * sinfo.nobs);
          cnt++;
        }
      }
      cufftExecZ2Z(plan_trk, cuda_xcor_fx, cuda_xcor_fx, CUFFT_INVERSE);

      cnt = 0;
      for (i = 0; i < sinfo.nch; i++)
      {
        if (cinfo[i].is_trk == 0x1)
        {
          fxres[i].ttag_gd[p] = (double)p * sinfo.period;
          fxres[i].ttag_phi[p] = (double)p * sinfo.period + (double)cinfo[i].pt / sinfo.fs;
          conv<<<1, 1000>>>(sinfo.nobs, nobs, cuda_xcor_fx + cnt * sinfo.nobs, cuda_xcor_conv, cuda_xcor_phi);
          cublasIdamax(handle, nobs * 2 + 1, cuda_xcor_conv, 1, &peak_idx); // peak
          peak_idx -= 1;
          fxres[i].pidx[p] = peak_idx - nobs;
          if (peak_idx - nobs >= 0)
            cublasDznrm2(handle, 1, (cuDoubleComplex *)cuda_xcor_fx + peak_idx - nobs + cnt * sinfo.nobs, 1, &peak);
          else
            cublasDznrm2(handle, 1, (cuDoubleComplex *)cuda_xcor_fx + peak_idx - nobs + (cnt + 1) * sinfo.nobs, 1, &peak);

          // compute SNR
          // peak = A * Pc
          // amplitude: A (in V) = peak / Pc
          // reception power: Px (in V^2) = A^2 * Pc + Pn = rx_power
          // reference signal power: Pc (in V^2) = cinfo[i].prn_power
          // signal power: Ps (in V^2) = A^2 * Pc = peak * peak / Pc = fxres[i].signal
          // noise power: Pn (in V^2) = Px - Ps = fxres[i].noise
          // SNR = A^2 * Pc / Pn = 1.0 / ((Px * Pc / peak) - 1.0)
          fxres[i].amp[p] = peak / cinfo[i].prn_power; // amplitude, used for SIC and code head decision
          fxres[i].signal[p] = peak * peak / cinfo[i].prn_power; // signal power (in V^2)
          if (rpow[i] > fxres[i].signal[p]) fxres[i].noise[p] = rpow[i] - fxres[i].signal[p]; // noise power (in V^2)
          else fxres[i].noise[p] = 5.0e-7; // if Pn > Px, then noise power can be ignored (assigned by a small value)

          // proceed if peak_idx is within the range of the delay spread, and SNR meets min SNR requirement
          snr_min = pow(10.0, cinfo[i].snr_min / 10.0);
          if (peak_idx - 2 >= 0 && peak_idx + 2 < nobs * 2 + 1 && (1.0 + snr_min) * fxres[i].signal[p] > snr_min * rpow[i])
          {
            cudaMemcpy(xcor_conv, cuda_xcor_conv, sizeof(double) * (nobs * 2 + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(xcor_phi, cuda_xcor_phi, sizeof(double) * (nobs * 2 + 1), cudaMemcpyDeviceToHost);

            // compute the carrier phase
            fxres[i].phi[p] = xcor_phi[peak_idx];

            // compute the code phase
            x1 = xcor_conv[peak_idx - 2]; x2 = xcor_conv[peak_idx - 1]; x3 = xcor_conv[peak_idx];
            x4 = xcor_conv[peak_idx + 1]; x5 = xcor_conv[peak_idx + 2];
            // perform Narrow Correlator
            // fxres[i].gd[p] = ((x2 - x4) / (x2 - 2.0 * x3 + x4) / 2.0 + (double)(cinfo[i].pt + peak_idx - nobs)) * 1.0e+9 / sinfo.fs;
            // perform High Resolution Correlator
            fxres[i].gd[p] = ((x2 - x4) / (x2 - 2.0 * x3 + x4) - (x1 - x5) / (x1 - 2.0 * x3 + x5) + (double)(cinfo[i].pt + peak_idx - nobs)) * 1.0e+9 / sinfo.fs;
            if (fxres[i].gd[p] > sinfo.period * 1.0e+9) fxres[i].gd[p] -= (sinfo.period * 1.0e+9);
            else if (fxres[i].gd[p] < 0.0) fxres[i].gd[p] += (sinfo.period * 1.0e+9);
            fxres[i].w[p] = 1.0;
            fxres[i].count++;
          }
          cnt++;
        }
      } // End: measure code phase, amplitude, and dphi
      // Start: Successive Inteference Cancellation (SIC) for the (p - 1)th obs
      blocks = sinfo.nobs / 500;
      if (p > 1)
      {
        // use local phase and aquired frequency to down convert, and then accumulate
        cnt = 0;
        for (i = 0; i < sinfo.nch; i++)
        {
          if (cinfo[i].is_sic == 0x1 && cinfo[i].is_trk == 0x1) // the i-th desired signal
          {
            if (cinfo[i].is_chA == 0x1) idx = 2 * ((p - 1) * sinfo.nobs + cinfo[i].pt);
            else idx = 2 * ((p - 1) * sinfo.nobs + cinfo[i].pt);
            blocks = sinfo.nobs / 1000;
            binary_to_waveform<<<blocks, 1000>>>(sinfo.nobs, cuda_samples + idx, cuda_ddbs, cinfo[i].is_chA);
            for (k = 0; k < sinfo.nch; k++)
            {
              if (i != k && cinfo[i].is_chA == cinfo[k].is_chA && cinfo[k].is_trk) // the k-th interference
              {
                if (cinfo[k].pt > cinfo[i].pt) // desired signal comes earlier than the interference, use p - 1 & p - 2
                {
                  if (fxres[k].w[p - 2] != 0.0) // cancel head, p - 2
                    SIC<<<blocks, 1000>>>(cuda_ddbs, sinfo.nobs, 0, sinfo.nobs + cinfo[i].pt - cinfo[k].pt, sinfo.nobs - 1               , fxres[k].amp[p - 2], cinfo[k].cuda_prn_t, fxres[k].pidx[p - 2], (cinfo[k].fc + cinfo[k].df) / sinfo.fs, fxres[k].phi[p - 2] + cinfo[k].fc * (double)((p - 2) * sinfo.nobs + cinfo[k].pt) / sinfo.fs - cinfo[k].df * sinfo.period);
                  if (fxres[k].w[p - 1] != 0.0) // cancel tail, p - 1
                    SIC<<<blocks, 1000>>>(cuda_ddbs, sinfo.nobs, cinfo[k].pt - cinfo[i].pt, 0, sinfo.nobs + cinfo[i].pt - cinfo[k].pt - 1, fxres[k].amp[p - 1], cinfo[k].cuda_prn_t, fxres[k].pidx[p - 1], (cinfo[k].fc + cinfo[k].df) / sinfo.fs, fxres[k].phi[p - 1] + cinfo[k].fc * (double)((p - 1) * sinfo.nobs + cinfo[k].pt) / sinfo.fs - cinfo[k].df * sinfo.period);
                }
                else // desired signal comes later than the interference, use p - 1 & p
                {
                  if (fxres[k].w[p - 1] != 0.0) // cancel head, p - 1
                    SIC<<<blocks, 1000>>>(cuda_ddbs, sinfo.nobs, 0, cinfo[i].pt - cinfo[k].pt             , sinfo.nobs - 1               , fxres[k].amp[p - 1], cinfo[k].cuda_prn_t, fxres[k].pidx[p - 1], (cinfo[k].fc + cinfo[k].df) / sinfo.fs, fxres[k].phi[p - 1] + cinfo[k].fc * (double)((p - 1) * sinfo.nobs + cinfo[k].pt) / sinfo.fs - cinfo[k].df * sinfo.period);
                  if (fxres[k].w[p] != 0.0) // cancel tail, p
                    SIC<<<blocks, 1000>>>(cuda_ddbs, sinfo.nobs, sinfo.nobs + cinfo[k].pt - cinfo[i].pt, 0, cinfo[i].pt - cinfo[k].pt - 1, fxres[k].amp[p]    , cinfo[k].cuda_prn_t, fxres[k].pidx[p]    , (cinfo[k].fc + cinfo[k].df) / sinfo.fs, fxres[k].phi[p]     + cinfo[k].fc * (double)(p       * sinfo.nobs + cinfo[k].pt) / sinfo.fs - cinfo[k].df * sinfo.period);
                }
              }
            }
            phi = fmod(cinfo[i].fc * (double)((p - 1) * sinfo.nobs + cinfo[i].pt) / sinfo.fs, 1.0);
            down_conversion<<<blocks, 1000>>>(sinfo.nobs, cinfo[i].fc / sinfo.fs, phi, cuda_ddbs, cuda_robs + cnt * sinfo.nobs);
            cnt++;
          }
        }

        // perform cross correlation
        cufftExecZ2Z(plan_sic, cuda_robs, cuda_robs, CUFFT_FORWARD);
        cnt = 0;
        for (i = 0; i < sinfo.nch; i++)
        {
          if (cinfo[i].is_sic == 0x1 && cinfo[i].is_trk == 0x1)
          {
            cross_spectrum<<<blocks, 1000>>>(sinfo.nobs, cuda_robs + cnt * sinfo.nobs, cinfo[i].cuda_prn, sinfo.fs / (double)sinfo.nobs, cinfo[i].fmax, cinfo[i].fmin, cuda_xcor_fx + cnt * sinfo.nobs);
            cnt++;
          }
        }
        cufftExecZ2Z(plan_sic, cuda_xcor_fx, cuda_xcor_fx, CUFFT_INVERSE);

        cnt = 0;
        for (i = 0; i < sinfo.nch; i++)
        {
          if(cinfo[i].is_sic == 0x1 && cinfo[i].is_trk == 0x1)
          {
            conv<<<1, 1000>>>(sinfo.nobs, nobs, cuda_xcor_fx + cnt * sinfo.nobs, cuda_xcor_conv, cuda_xcor_phi);
            cublasIdamax(handle, nobs * 2 + 1, cuda_xcor_conv, 1, &peak_idx); // peak
			      peak_idx -= 1;
            cudaMemcpy(xcor_conv, cuda_xcor_conv, sizeof(double) * (nobs * 2 + 1), cudaMemcpyDeviceToHost);
            cudaMemcpy(xcor_phi, cuda_xcor_phi, sizeof(double) * (nobs * 2 + 1), cudaMemcpyDeviceToHost);
            if (peak_idx - 2 >= 0 && peak_idx + 2 < nobs * 2 + 1) // proceed if peak_idx is in the range of the delay spread
            {
              // compute the code phase
              fxres[i].ttag_gd[p - 1] = (double)(p - 1) * sinfo.period;
              x1 = xcor_conv[peak_idx - 2]; x2 = xcor_conv[peak_idx - 1]; x3 = xcor_conv[peak_idx];
              x4 = xcor_conv[peak_idx + 1]; x5 = xcor_conv[peak_idx + 2];
              // perform Narrow Correlator
              // fxres[i].gd_sic[p - 1] = ((x2 - x4) / (x2 - 2.0 * x3 + x4) / 2.0 + (double)(cinfo[i].pt + peak_idx - nobs)) * 1.0e+9 / sinfo.fs;
              // perform High Resolution Correlator
              fxres[i].gd_sic[p - 1] = ((x2 - x4) / (x2 - 2.0 * x3 + x4) - (x1 - x5) / (x1 - 2.0 * x3 + x5) + (double)(cinfo[i].pt + peak_idx - nobs)) * 1.0e+9 / sinfo.fs;
              // remove delay code phase ambiguity
              if (fxres[i].gd_sic[p - 1] > sinfo.period * 1.0e+9) fxres[i].gd_sic[p - 1] -= (sinfo.period * 1.0e+9);
              else if (fxres[i].gd_sic[p - 1] < 0.0) fxres[i].gd_sic[p - 1] += (sinfo.period * 1.0e+9);
              fxres[i].w_sic[p - 1] = 1.0;
              fxres[i].count_sic++;
            }
            cnt++;
          }
        }
      } // End: SIC
    } // End: Tracking

    cufftDestroy(plan_trk);
    cufftDestroy(plan_sic);

    if (sinfo.ntrk != 0) printf("acquisition: %d , lock: %d \n", sinfo.nch - sinfo.ntrk, sinfo.ntrk);
    for (i = 0; i < sinfo.nch; i++)
    {
      if (cinfo[i].is_trk == 0x1)
      {
        if (fxres[i].count_sic * 2 > sinfo.portion && cinfo[i].is_sic == 0x1 && cinfo[i].is_first == 0x0) // SIC result
        {
          // Start: update and record
          if (cinfo[i].is_chA == 0x1) sprintf(filename, "%s/%4d%02d%02d%02d.chA.pn%02d.%04dkcps.sic.dat", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour, cinfo[i].prnno, cinfo[i].rc / 1000);
          else sprintf(filename, "%s/%4d%02d%02d%02d.chB.pn%02d.%04dkcps.sic.dat", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour, cinfo[i].prnno, cinfo[i].rc / 1000);
          fxres[i].fout = fopen(filename, "a");
          fprintf(fxres[i].fout, "%2d %2d %2d ", ttag->hour, ttag->minute, ttag->second);

          // filter out the obs larger than 3 sigma
          ii = 0;
          memset(res, 0x0, sizeof(double) * sinfo.portion);
          for (p = 0; p < sinfo.portion; p++)
          {
            if (fxres[i].w_sic[p] > 0.0)
            {
              res[ii] = fxres[i].gd_sic[p];
              ii++;
            }
          }
          c0 = kth_smallest(res, ii, ii / 2);
          chisq = (kth_smallest(res, ii, ii * 3 / 4) - kth_smallest(res, ii, ii / 4)) / 1.349;
          fxres[i].count_sic = 0;
          for (p = 0; p < sinfo.portion; p++)
          {
            if (fxres[i].w_sic[p] != 0.0)
              if (fabs(fxres[i].gd_sic[p] - c0) > 3.0 * chisq)
                fxres[i].w_sic[p] = 0.0;
            if (fxres[i].w_sic[p] > 0.0) fxres[i].count_sic++;
          }

          // code phase
          gsl_fit_wlinear(fxres[i].ttag_gd, 1, fxres[i].w_sic, 1, fxres[i].gd_sic, 1, sinfo.portion, &c0, &c1, &c00, &c01, &c11, &chisq);
          chisq = sqrt(chisq / (double)fxres[i].count_sic);
          fprintf(fxres[i].fout, "%5d %3d %14.6lf %14.6lf %8.3lf\n", ttag->hour * 3600 + ttag->minute * 60 + ttag->second, fxres[i].count_sic, c0 + 0.5 * c1, c1, chisq);
          fclose(fxres[i].fout); // close output file

          // End: update and record
          // reset all the measurements
          memset(fxres[i].gd_sic , 0, sizeof(double) * sinfo.portion);
        }

        // Normal result
        if (fxres[i].count * 2 > sinfo.portion)
        {
          // find head of the code sequence 
          imax = 0;
          tmp = 0.0;
          if (cinfo[i].is_first == 0x0)
          {
            memset(res, 0x0, sizeof(double) * sinfo.portion);
            if (cinfo[i].rc == 1000000) // for 1Mcps SATRE code
            {
              fxres[i].amp[sinfo.portion - 1] = (fxres[i].amp[0] + fxres[i].amp[sinfo.portion - 2]) / 2.0;
              for (ii = 0; ii < sinfo.portion; ii++)
              {
                for (jj = 0; jj < 16; jj++)
                {
                  idx = (ii + jj) % sinfo.portion;
                  res[ii] += fxres[i].amp[idx];
                }
              }
              tmp = res[0];
              for (ii = 1; ii < sinfo.portion; ii++)
              {
                if (res[ii] < tmp)
                {
                  imax = ii;
                  tmp = res[ii];
                }
              }
              for (ii = -1; ii < 17; ii++) // -1 ~ 16
              {
                jj = (ii + imax + sinfo.portion) % sinfo.portion;
                fxres[i].w[jj] = 0.0;
              }
            }
            else if (cinfo[i].rc == 2500000) // for 2.5Mcps SATRE code
            {
              tmp = kth_smallest(fxres[i].gd, sinfo.portion, sinfo.portion / 2); // find the median among the obs
              for (ii = 0; ii < sinfo.portion; ii++) // fill the median
              {
                if (fxres[i].w[ii] == 0.0)
                  fxres[i].gd[ii] = tmp;
              }
              for (ii = 0; ii < sinfo.portion - 1; ii++)
                res[ii] = fxres[i].gd[ii] - fxres[i].gd[ii + 1];
              res[sinfo.portion - 1] = fxres[i].gd[sinfo.portion - 1] - fxres[i].gd[0];
              tmp = res[0];
              for (ii = 1; ii < sinfo.portion; ii++)
              {
                if (res[ii] > tmp)
                {
                  imax = ii;
                  tmp = res[ii];
                }
              }
              fxres[i].w[imax] = 1.0;
              fxres[i].gd[imax] -= 200.0;
              fxres[i].w[(imax + 1) % sinfo.portion] = 1.0;
              fxres[i].gd[(imax + 1) % sinfo.portion] += 200.0;
            }
          }

          // filter out the obs larger than 3 sigma
          ii = 0;
          memset(res, 0x0, sizeof(double) * sinfo.portion);
          for (p = 0; p < sinfo.portion; p++)
          {
            if (fxres[i].w[p] > 0.0)
            {
              res[ii] = fxres[i].gd[p];
              ii++;
            }
          }
          c0 = kth_smallest(res, ii, ii / 2);
          chisq = (kth_smallest(res, ii, ii * 3 / 4) - kth_smallest(res, ii, ii / 4)) / 1.349;
          fxres[i].count = 0;
          for (p = 0; p < sinfo.portion; p++)
          {
            if (fxres[i].w[p] != 0.0)
            {
              if (fabs(fxres[i].gd[p] - c0) < 3.0 * chisq)
              {
                fxres[i].count++;
                while (fabs(fxres[i].phi[p] - cinfo[i].last_phi) > 0.25) // BPSK phase adjustment
                {
                  if (fxres[i].phi[p] > cinfo[i].last_phi) fxres[i].phi[p] -= 0.5;
                  else fxres[i].phi[p] += 0.5;
                }
                cinfo[i].last_phi = fxres[i].phi[p];
              }
              else fxres[i].w[p] = 0.0;
            }
          }

          // Start: update and record
          if (cinfo[i].is_first == 0x0)
          {
            if (cinfo[i].is_chA == 0x1) sprintf(filename, "%s/%4d%02d%02d%02d.chA.pn%02d.%04dkcps.dat", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour, cinfo[i].prnno, cinfo[i].rc / 1000);
            else sprintf(filename, "%s/%4d%02d%02d%02d.chB.pn%02d.%04dkcps.dat", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour, cinfo[i].prnno, cinfo[i].rc / 1000);
            fxres[i].fout = fopen(filename, "a");
            fprintf(fxres[i].fout, "%2d %2d %2d ", ttag->hour, ttag->minute, ttag->second);
          }

          // apply linear fit, and then update fc and pt
          gsl_fit_wlinear(fxres[i].ttag_phi, 1, fxres[i].w, 1, fxres[i].phi, 1, sinfo.portion, &c0, &c1, &c00, &c01, &c11, &chisq);
          cinfo[i].fc += floor(c1);
          cinfo[i].df = c1 - floor(c1);
          cinfo[i].phi = fmod(c0 + 1000.0, 1.0);
          if (cinfo[i].is_chA == 0x1) printf("A: #%02d %4.1lf Mcps %12.3lf Hz ", cinfo[i].prnno, (double)cinfo[i].rc * 1.0e-6, cinfo[i].fc + cinfo[i].df);
          else printf("B: #%02d %4.1lf Mcps %12.3lf Hz ", cinfo[i].prnno, (double)cinfo[i].rc * 1.0e-6, cinfo[i].fc + cinfo[i].df);
          if (cinfo[i].is_first == 0x0) fprintf(fxres[i].fout, "%14.6lf %11.8lf ", cinfo[i].fc + c1 - floor(c1), cinfo[i].phi);

          // code phase
          gsl_fit_wlinear(fxres[i].ttag_gd, 1, fxres[i].w, 1, fxres[i].gd, 1, sinfo.portion, &c0, &c1, &c00, &c01, &c11, &chisq);
          chisq = sqrt(chisq / (double)fxres[i].count);
          cinfo[i].gd = c0 + 0.5 * c1;
          cinfo[i].dg = c1;
          cinfo[i].pt = (int)round((c0 + c1) * sinfo.fs / 1.0e+9);
          printf("%13.3lf (%5.3lf) ns ", cinfo[i].gd + (double)imax * sinfo.period * 1.0e+9, chisq / sqrt((double)fxres[i].count));
          if (cinfo[i].is_first == 0x0) fprintf(fxres[i].fout, "%5d %3d %5.3lf %14.6lf %11.6lf %8.3lf ", ttag->hour * 3600 + ttag->minute * 60 + ttag->second, fxres[i].count, (double)imax * sinfo.period, cinfo[i].gd, cinfo[i].dg, chisq);

          // signal power in dBm
          cinfo[i].peak = average(sinfo.portion, fxres[i].signal, fxres[i].w);
          if (cinfo[i].is_first == 0x0) fprintf(fxres[i].fout, "%7.3lf ", 10.0 * log10(cinfo[i].peak * 1000.0 / 50.0));

          // noise power in dBm
          peak = average(sinfo.portion, fxres[i].noise, fxres[i].w);
          printf("S/N %6.2lf dB\n", 10.0 * log10(cinfo[i].peak / peak)); // show S/N in screen
          if (cinfo[i].is_first == 0x0) fprintf(fxres[i].fout, "%7.3lf\n", 10.0 * log10(peak * 1000.0 / 50.0));

          if (cinfo[i].is_first == 0x0) fclose(fxres[i].fout); // close output file

          // End: update and record
          // reset all the measurements
          memset(fxres[i].phi    , 0, sizeof(double) * sinfo.portion);
          memset(fxres[i].gd     , 0, sizeof(double) * sinfo.portion);
          memset(fxres[i].noise  , 0, sizeof(double) * sinfo.portion);
          memset(fxres[i].signal , 0, sizeof(double) * sinfo.portion);
          cinfo[i].is_first = 0x0;
        }
        else
        {
          // LOG
          sprintf(filename, "%s/%4d%02d%02d%02d.log", result_dir, ttag->year, ttag->month, ttag->day, ttag->hour);
          flog = fopen(filename, "a");
          if (cinfo[i].is_chA == 0x1)
          {
            //printf("Ch. A, PRN#%2d, count = %d / %d ,loss of lock\n", cinfo[i].prnno, fxres[i].count, sinfo.portion);
            fprintf(flog, "%02d:%02d:%02d + %11.9lf loss of lock: Ch. A, PRN#%2d, count = %d / %d\n", ttag->hour, ttag->minute, ttag->second, ttag->fsec, cinfo[i].prnno, fxres[i].count, sinfo.portion);
          }
          else
          {
            //printf("Ch. B, PRN#%2d, count = %d / %d ,loss of lock\n", cinfo[i].prnno, fxres[i].count, sinfo.portion);
            fprintf(flog, "%02d:%02d:%02d + %11.9lf loss of lock: Ch. B, PRN#%2d, count = %d / %d\n", ttag->hour, ttag->minute, ttag->second, ttag->fsec, cinfo[i].prnno, fxres[i].count, sinfo.portion);
          }
          fclose(flog);
          // LOG
          cinfo[i].fc = 0.0;
          cinfo[i].is_first = 0x1;
          cinfo[i].is_trk = 0x0;
          cinfo[i].dg = 0.0;
          cinfo[i].phi = 0.0;
          cinfo[i].last_phi = 0.0;
          // reset all the measurements
          memset(fxres[i].phi    , 0, sizeof(double) * sinfo.portion);
          memset(fxres[i].gd     , 0, sizeof(double) * sinfo.portion);
          memset(fxres[i].noise  , 0, sizeof(double) * sinfo.portion);
          memset(fxres[i].signal , 0, sizeof(double) * sinfo.portion);
          //continue;
        }
      }
    }
  }
  return 0;
}

double kth_smallest(double *a, int n, int k) // find the k-th smallest value in the array a[]
{
  int i, j, l, m;
  double x, *b, t;
  b = (double *)malloc(sizeof(double) * n);
  memcpy(b, a, sizeof(double) * n);
  l = 0; m = n - 1;
  while (l < m)
  {
    x = b[k]; i = l; j = m;
    do
    {
      while (b[i] < x) i++;
      while (b[j] > x) j--;
      if (i <= j)
      {
        t = b[i]; b[i] = b[j]; b[j] = t;
        i++; j--;
      }
    } while (i <= j);
    if (j < k) l = i;
    if (k < i) m = j;
  }
  x = b[k];
  free(b);
  return x;
}

int SATREcode(channel_info *c)
{
  int i, j, shift_back, reg[14];
  int tap[32][14] = {{  1,	1,  0,	0,  0,	0,  1,	0,  1,	0,  0,	0,  0,	0}, /* PN # 0 */
                     {  1,	1,  0,	0,  0,	0,  1,	0,  0,	0,  0,	1,  0,	0}, /* PN # 1 */
                     {  1,	0,  1,	0,  1,	0,  0,	0,  0,	0,  0,	0,  0,	1}, /* PN # 2 */
                     {  1,	0,  0,	0,  0,	0,  0,	0,  0,	1,  0,	1,  0,	1}, /* PN # 3 */
                     {  1,	0,  0,	0,  0,	0,  1,	0,  0,	0,  0,	1,  1,	0}, /* PN # 4 */
                     {  1,	1,  0,	0,  0,	0,  0,	0,  0,	0,  0,	1,  1,	0}, /* PN # 5 */
                     {  1,	1,  0,	0,  0,	0,  1,	0,  0,	0,  1,	0,  0,	0}, /* PN # 6 */
                     {  1,	1,  1,	0,  0,	0,  0,	0,  0,	0,  0,	0,  1,	0}, /* PN # 7 */
                     {  1,	0,  0,	1,  0,	0,  1,	1,  0,	0,  0,	0,  1,	1}, /* PN # 8 */
                     {  1,	1,  0,	0,  0,	0,  0,	0,  0,	1,  0,	1,  1,	1}, /* PN # 9 */
                     {  1,	0,  0,	1,  0,	0,  0,	0,  1,	1,  0,	0,  1,	1}, /* PN #10 */
                     {  1,	0,  0,	1,  0,	0,  1,	1,  0,	0,  0,	0,  1,	1}, /* PN #11 */
                     {  1,	0,  0,	0,  1,	0,  0,	0,  1,	1,  1,	0,  0,	1}, /* PN #12 */
                     {  1,	0,  1,	1,  0,	0,  0,	1,  0,	1,  0,	0,  0,	1}, /* PN #13 */
                     {  1,	0,  1,	0,  0,	0,  0,	0,  1,	0,  1,	1,  1,	0}, /* PN #14 */
                     {  1,	0,  0,	0,  0,	1,  0,	1,  1,	0,  0,	1,  1,	0}, /* PN #15 */
                     {  1,  0,  0,  0,  0,  1,  1,  1,  0,  0,  1,  0,  1,  0}, /* PN #16 */
                     {  1,  0,  0,  1,  1,  1,  0,  0,  0,  0,  1,  0,  1,  0}, /* PN #17 */
                     {  1,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0,  1,  0}, /* PN #18 */
                     {  1,  0,  0,  0,  1,  0,  1,  1,  1,  0,  0,  0,  1,  0}, /* PN #19 */
                     {  1,  1,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  1,  0}, /* PN #20 */
                     {  1,  0,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  1,  0}, /* PN #21 */
                     {  1,  0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  1,  0,  0}, /* PN #22 */
                     {  1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  1,  1,  0,  0}, /* PN #23 */
                     {  1,  1,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0}, /* PN #24 */
                     {  1,  0,  0,  1,  0,  1,  1,  0,  1,  0,  0,  1,  0,  0}, /* PN #25 */
                     {  1,  0,  0,  0,  1,  1,  0,  0,  1,  1,  1,  0,  0,  0}, /* PN #26 */
                     {  1,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1,  0,  0,  0}, /* PN #27 */
                     {  1,  1,  0,  0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0}, /* PN #28 */
                     {  1,  1,  1,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0}, /* PN #29 */
                     {  1,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  0}, /* PN #30 */
                     {  1,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0}};/* PN #31 */
  if (c->prnno < 0 || c->prnno > 31)
  {
    printf("no such prn number, %d\n", c->prnno);
    return 1;
  }

  for (i = 0; i < 14; i++)
    reg[i] = 1;
  for (i = 0; i < c->clen; i++)
  {
    if (c->clen == 4000 && i < 10)
      c->code[i] = 0; // 0
    else
      c->code[i] = 1 - 2 * reg[0]; // -1, 1
    shift_back = reg[0];
    for (j = 1; j < 14; j++)
      if (tap[c->prnno][j] == 1)
        shift_back ^= reg[j];
    for (j = 1; j < 14; j++)
      reg[j - 1] = reg[j];
    reg[13] = shift_back;
  }
  return 0;
}

int NICTcode(channel_info *c)
{
  int i, j, shift_back, reg[12];
  int tap[8][12] = {{ 1,	0,  0,	0,  0,	0,  1,	0,  1,	0,  0,	1},	 /* PN #0 */
                    { 1,	1,  0,	1,  1,	1,  0,	1,  0,	0,  1,	1},	 /* PN #1 */
                    { 1,	0,  0,	1,  0,	0,  0,	0,  0,	1,  1,	0},	 /* PN #2 */
                    { 1,	0,  1,	1,  1,	0,  1,	1,  1,	0,  1,	0},	 /* PN #3 */
                    { 1,	1,  0,	0,  0,	1,  0,	1,  1,	1,  0,	0},	 /* PN #4 */
                    { 1,	0,  0,	1,  1,	1,  0,	1,  0,	0,  0,	1},	 /* PN #5 */
                    { 1,	1,  1,	0,  0,	1,  0,	1,  1,	1,  0,	1},	 /* PN #6 */
                    { 1,	0,  1,	0,  1,	1,  1,	0,  1,	1,  1,	0}}; /* PN #7 */
  if (c->prnno < 0 || c->prnno > 7)
  {
    printf("no such prn number, %d\n", c->prnno);
    return 1;
  }

  for (i = 0; i < 12; i++)
    reg[i] = 1;
  for (i = 0; i < 4095; i++)
  {
    shift_back = reg[0];
    for (j = 1; j < 12; j++)
      if (tap[c->prnno][j] == 1) shift_back ^= reg[j];

    c->code[i] = shift_back; // 0, 1

    for (j = 1; j < 12; j++)
      reg[j - 1] = reg[j];
    reg[11] = shift_back;
  }
  return 0;
}

/* GNSS Coarse Aquisition code */
int CAcode(channel_info *c)
{
  int i, j, res, shift_back;
  int code_phase[37][2] = {{ 2,  6}, { 3,  7}, { 4,  8}, { 5,  9}, { 1,  9},
                           { 2, 10}, { 1,  8}, { 2,  9}, { 3, 10}, { 2,  3},
                           { 3,  4}, { 5,  6}, { 6,  7}, { 7,  8}, { 8,  9},
                           { 9, 10}, { 1,  4}, { 2,  5}, { 3,  6}, { 4,  7},
                           { 5,  8}, { 6,  9}, { 1,  3}, { 4,  6}, { 5,  7},
                           { 6,  8}, { 7,  9}, { 8, 10}, { 1,  6}, { 2,  7},
                           { 3,  8}, { 4,  9}, { 5, 10}, { 4, 10}, { 1,  7},
                           { 2,  8}, { 4, 10}};
  int *G1_reg, *G2_reg;
  if (c->prnno < 1 || c->prnno > 37)
  {
    printf("no such prn number, %d\n", c->prnno);
    return 1;
  }
  G1_reg = (int *)malloc(sizeof(int) * 10);
  G2_reg = (int *)malloc(sizeof(int) * 10);
  for (i = 0; i < 10; i++)
  {
    G1_reg[i] = 1;
    G2_reg[i] = 1;
  }
  for (i = 0; i < 1023; i++)
  {
    res = G1_reg[9] ^ G2_reg[code_phase[c->prnno - 1][0] - 1] ^ G2_reg[code_phase[c->prnno - 1][1] - 1];
    c->code[i] = res; // 0, 1

    /* G1 */
    shift_back = G1_reg[2] ^ G1_reg[9];
    for (j = 9; j >= 1; j--)
      G1_reg[j] = G1_reg[j - 1];
    G1_reg[0] = shift_back;

    /* G2 */
    shift_back = G2_reg[1] ^ G2_reg[2] ^ G2_reg[5] ^ G2_reg[7] ^ G2_reg[8] ^ G2_reg[9];
    for (j = 9; j >= 1; j--)
      G2_reg[j] = G2_reg[j - 1];
    G2_reg[0] = shift_back;
  }
  free(G1_reg);
  free(G2_reg);
  return 0;
}

double average(int nobs, double *x, double *w)
{
	int i, count = 0;
	double res = 0.0;
	if (nobs == 0)
	  return 0.0;
	for (i = 0; i < nobs; i++)
	{
		if (w[i] > 0.0)
		{
			res += x[i];
			count++;
		}
	}
	return res / (double)count;
}

__global__ void SIC(double *obs, int nobs, int start_pt, int init_pt, int end_pt, double amp, cufftDoubleComplex *prn, int pidx, double ff, double phi)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= start_pt && i <= start_pt + end_pt - init_pt)
    obs[i] -= amp * prn[i - start_pt + init_pt].x * cos(2.0 * PI * (ff * (double)(i - start_pt + init_pt) + phi));
}

__global__ void PRN_sampling(int nobs, int *code, cufftDoubleComplex *prn, int rc, double fs, int clen, double delay)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x, idx;
  if (i < nobs)
  {
    idx = (int)floor(fmod(((double)i / fs - delay * 1.0e-9) * (double)rc, (double)clen));
    if (idx < 0)
      idx += clen;
    else if (idx >= clen)
      idx -= clen;
    prn[i].x = (double)code[idx];
    prn[i].y = 0.0;
  }
}

__global__ void cross_spectrum(int nobs, cufftDoubleComplex *obs, cufftDoubleComplex *prn, double df, double fmax, double fmin, cufftDoubleComplex *robs)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x, idx;
  if (i < nobs)
  {
    idx = (i >= nobs / 2) ? i - nobs : i;
    if ((double)idx * df < fmax && (double)idx * df > fmin && idx != 0)
    {
      robs[i].x = (obs[i].x * prn[i].x + obs[i].y * prn[i].y) / (double)nobs / (double)nobs;
      robs[i].y = (obs[i].y * prn[i].x - obs[i].x * prn[i].y) / (double)nobs / (double)nobs;
    }
    else
    {
      robs[i].x = 0.0;
      robs[i].y = 0.0;
    }
  }
}

/* down convertion, cuda_obs -> cuda_dobs
 * ff : digital frequency
 * phi : initial phase
 */
__global__ void down_conversion(int nobs, double ff, double phi, double *obs, cufftDoubleComplex *dobs)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nobs)
  {
    dobs[i].x = 1.4142135624 * obs[i] * cos(-1.0 * (2.0 * PI * ff * (double)i + phi * 2.0 * PI));
    dobs[i].y = 1.4142135624 * obs[i] * sin(-1.0 * (2.0 * PI * ff * (double)i + phi * 2.0 * PI));
  }
}

__global__ void down_conversion2(int nobs, double ff, double phi, char *samp, cufftDoubleComplex *dobs, bool is_chA)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nobs)
  {
    if (is_chA == 0x1)
    {
      dobs[i].x = 1.4142135624 * (double)samp[i * 2] * cos(-1.0 * (2.0 * PI * ff * (double)i + phi * 2.0 * PI)) / 128.0;
      dobs[i].y = 1.4142135624 * (double)samp[i * 2] * sin(-1.0 * (2.0 * PI * ff * (double)i + phi * 2.0 * PI)) / 128.0;
    }
    else
    {
      dobs[i].x = 1.4142135624 * (double)samp[i * 2 + 1] * cos(-1.0 * (2.0 * PI * ff * (double)i + phi * 2.0 * PI)) / 128.0;
      dobs[i].y = 1.4142135624 * (double)samp[i * 2 + 1] * sin(-1.0 * (2.0 * PI * ff * (double)i + phi * 2.0 * PI)) / 128.0;
    }
  }
}

__global__ void binary_to_waveform(int nobs, char *samp, double *obs, bool is_chA)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nobs)
  {
    if (is_chA == 0x1) obs[i] = (double)samp[2 * i] / 128.0;
    else obs[i] = (double)samp[2 * i + 1] / 128.0;
  }
}

__global__ void conv(int snobs, int nobs, cufftDoubleComplex *in, double *out, double *phi) // Cartesian to polar coordinate conversion
{
  int i = blockIdx.x * blockDim.x + threadIdx.x, idx;
  if (i <= 2 * nobs)
  {
    idx = (i - nobs < 0) ? i - nobs + snobs : i - nobs;
    out[i] = (in[idx].x * in[idx].x + in[idx].y * in[idx].y);
    phi[i] = atan2(in[idx].y, in[idx].x) / 2.0 / PI; // phase in cycle
  }
}

/*
	current time (UTC + 0)
*/
void current_time(int offset, int *mjd, int *year, int *month, int *day, int *doy, int *hour, int *minute, int *second)
{
	time_t sec;
	sec = (int)time(NULL) + offset;
	*mjd = 40587 + sec / 86400;
	*hour = (sec % 86400) / 3600;
	*minute = (sec % 3600) / 60;
	*second = sec % 60;
	mjd2doy(*mjd, year, doy);
	mjd2date(*mjd, year, month, day);
//	return 40587.0 + (double)time(NULL) / 86400.0;
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
