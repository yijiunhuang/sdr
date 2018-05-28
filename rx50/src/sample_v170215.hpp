#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/thread/thread.hpp>
#include <iostream>
#include <fstream>
#include <csignal>
#include <complex>
#include <pthread.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#define PI 3.141592653589793
#define MAXBUF 10

const int sps = 50000000;
const double fs = 5.0e7;

typedef struct _timetag
{
  int mjd;
  int year;
  int doy;
  int month;
  int day;
  int sod;  /* second of day */
  int hour;
  int minute;
  int second;
  double fsec; /* frantional sec */
} timetag;

typedef struct _samplebuffer
{
  int readd, writed;
  char *buf[MAXBUF];
  bool is_tt;
  time_t tt;
  timetag *ttag[MAXBUF];
  pthread_mutex_t mymutex;
  pthread_cond_t mycond;
} samplebuffer;

void *sampling(void *);
void transmit_worker(short *, int, uhd::tx_streamer::sptr, uhd::tx_metadata_t);
