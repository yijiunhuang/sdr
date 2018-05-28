/*
* Input signal statistics V160824
*/
#include "sample_v170215.hpp"
#include <time.h>
#include <boost/thread/thread.hpp>

void current_time(int, int *, int *, int *, int *, int *, int *, int *, int *);
void date2doy(int , int, int, int *);
void doy2date(int, int, int *, int *);
int date2mjd(int, int, int);
void mjd2date(int, int *, int *, int *);
int doy2mjd(int, int);
void mjd2doy(int, int *, int *);

// sampling subroutine that should be called by a thread
int main()
{
	// variables to be set by "boost::program_options"
  std::string args = "addr=192.168.10.2", subdev = "A:AB";
  
	// create a usrp device
	uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
	
	// always select the subdevice first, the channel mapping affects the other settings
	usrp->set_rx_subdev_spec(subdev); // A:AB, A:BA, A:A, A:B
	std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;
	
	// get usrp parameter range informaion
	uhd::meta_range_t sr = usrp->get_tx_rates(0);
  //std::cout << boost::format("TX Sample Rate Range: %s") % sr.to_pp_string() << std::endl;
	sr = usrp->get_rx_rates(0);
	//std::cout << boost::format("RX Sample Rate Range: %s") % sr.to_pp_string() << std::endl;
	
	uhd::freq_range_t fr = usrp->get_tx_freq_range(0);
	std::cout << boost::format("TX Frequency Range: %s") % fr.to_pp_string() << std::endl;
	fr = usrp->get_rx_freq_range(0);
	std::cout << boost::format("RX Frequency Range: %s") % fr.to_pp_string() << std::endl;
	
	uhd::gain_range_t gr = usrp->get_tx_gain_range(0);
	std::cout << boost::format("TX Gain Range: %s") % gr.to_pp_string() << std::endl;
	gr = usrp->get_rx_gain_range(0);
	std::cout << boost::format("RX Gain Range: %s") % gr.to_pp_string() << std::endl;
	
	uhd::meta_range_t mr = usrp->get_tx_bandwidth_range(0);
  std::cout << boost::format("TX Bandwidth Range: %s") % mr.to_pp_string() << std::endl;
	mr = usrp->get_rx_bandwidth_range(0);
	std::cout << boost::format("RX Bandwidth Range: %s") % mr.to_pp_string() << std::endl;
	
	// Default Parameters
	std::cout << boost::format("Default TX Sample Rate: %f Msps") % (usrp->get_tx_rate(0) / 1.0e+6) << std::endl;
  std::cout << boost::format("Default TX Frequency: %f MHz") % (usrp->get_tx_freq(0) / 1.0e+6) << std::endl;
  std::cout << boost::format("Default TX Gain: %f") % usrp->get_tx_gain(0) << std::endl;
  std::cout << boost::format("Default TX Bandwidth: %f MHz") % (usrp->get_tx_bandwidth(0) / 1.0e+6) << std::endl;
	std::cout << boost::format("Default RX Sample Rate: %f Msps") % (usrp->get_rx_rate(0) / 1.0e+6) << std::endl;
	std::cout << boost::format("Default RX Frequency: %f MHz") % (usrp->get_rx_freq(0) / 1.0e+6) << std::endl;
	std::cout << boost::format("Default RX Gain: %f") % usrp->get_rx_gain(0) << std::endl;
	std::cout << boost::format("Default RX Bandwidth: %f MHz") % (usrp->get_rx_bandwidth(0) / 1.0e+6) << std::endl;
	std::cout << boost::format("Master Clock Rate: %f MHz") % (usrp->get_master_clock_rate(0) / 1.0e+6) << std::endl;
	std::cout << boost::format("PPS Source: %s") % usrp->get_time_source(0) << std::endl;
	std::cout << boost::format("Clock Source: %s") % usrp->get_clock_source(0) << std::endl;
	
	
	// Lock mboard clocks and external PPS
  usrp->set_clock_source("external");
  std::cout << boost::format("Actual Clock Source: %s") % usrp->get_clock_source(0) << std::endl;
  
	usrp->set_time_source("external"); //none, external, _external_, mimo
	std::cout << boost::format("Actual PPS Source: %s") % usrp->get_time_source(0) << std::endl;
	
	// set the sample rate, provided (in Msps) 0.2, 0.4, 0.5, 0.8, 1.0, 1.25, 2.0, 2.5, 3.125, 4.0, 5.0, 6.25, 10.0, 12.5, 20.0, 50.0(only A:A or A:B)
	usrp->set_tx_rate(1000000);
  std::cout << boost::format("Actual TX Sample Rate: %f Msps...") % (usrp->get_tx_rate(0)/1e6) << std::endl;
	usrp->set_rx_rate(50000000);
	std::cout << boost::format("Actual RX Sample Rate: %f Msps...") % (usrp->get_rx_rate(0)/1e6) << std::endl;
	
	//set the center frequency
	usrp->set_tx_freq(0.0);
  std::cout << boost::format("Actual TX Frequency: %f MHz...") % (usrp->get_tx_freq(0)/1e6) << std::endl;
	usrp->set_rx_freq(0.0);
	std::cout << boost::format("Actual RX Frequency: %f MHz...") % (usrp->get_rx_freq(0)/1e6) << std::endl;
	
	//set the rf gain
	usrp->set_tx_gain(0.0);
  std::cout << boost::format("Actual TX Gain: %f dB...") % usrp->get_tx_gain(0) << std::endl;
	usrp->set_rx_gain(0.0);
	std::cout << boost::format("Actual RX Gain: %f dB...") % usrp->get_rx_gain(0) << std::endl;
	
	// set the IF filter bandwidth
	usrp->set_tx_bandwidth(2.5e+8); // fixed for Basic TX and Basic RX: 250.0MHz
  std::cout << boost::format("Actual TX Bandwidth: %f MHz...") % (usrp->get_tx_bandwidth(0)/1.0e+6) << std::endl;
	usrp->set_rx_bandwidth(2.5e+8); // fixed for Basic TX and Basic RX: 250.0MHz
	std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % (usrp->get_rx_bandwidth(0)/1.0e+6) << std::endl;
	
	// get the time onboard
	uhd::time_spec_t lastpps;
	lastpps = usrp->get_time_last_pps();
	
	// allow for some setup time
	//boost::this_thread::sleep(boost::posix_time::seconds(1));
	
	// create a receive streamer
	// stream_args(host sample format, device wire format)
	// format: sc = short complex, fc = float complex
	// host format: fc64 = complex double, fc32 = complex float, sc16 = complex short (map to +-32768), sc8 = complex char (map to +-128)
	// over-the-wire format: sc8(up to 50Msps for host) or sc16(up to 25Msps)
	uhd::stream_args_t stream_args("sc8", "sc8");
	uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
	uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);
	
	// create a RX metadata
	uhd::rx_metadata_t md;
	
	// create a buffer
	//std::vector<std::complex<short> > buff((size_t)sps);
	FILE *fout;
  char cmd[80];
  char *buf;
  int i, res[256] = {0};
  double power = 0.0;
  size_t num_rx_samps;
  buf = (char *)malloc(sizeof(char) * sps); // previous, IQIQ..
  //ttag = (timetag *)malloc(sizeof(timetag));
	
  // setup continously streaming
  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.num_samps = sps;
	
  // last PPS + 2.0 seconds in the future to receive
  stream_cmd.stream_now = true;
  //stream_cmd.time_spec = lastpps + uhd::time_spec_t(1.0);
  
  // issue the streaming command (start sampling)
  usrp->issue_stream_cmd(stream_cmd);
  
  // control	
  while (1)
  {
    num_rx_samps = rx_stream->recv(buf, sps / 2, md);
    num_rx_samps = rx_stream->recv(buf, sps / 2, md);
    if (num_rx_samps >= sps / 2)
    {
      uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
      usrp->issue_stream_cmd(stream_cmd);
      break;
    }
  }
  // show the time tag
  //ttag->mjd = date2mjd(ttag->year, ttag->month, ttag->day);
  int mjd, yr, mon, day, doy, hr, min, sec;
  current_time(0, &mjd, &yr, &mon, &day, &doy, &hr, &min, &sec);
  printf("%5d %4d/%03d %4d/%02d/%02d %02d:%02d:%02d\n", mjd, yr, doy, yr, mon, day, hr, min, sec);
  // result output
  
  fout = fopen("plt.txt", "w");
  fprintf(fout, "set xlabel 'mV'\n");
  fprintf(fout, "set ylabel 'Histogram'\n");
  fprintf(fout, "set xtics (\"-1000\" 0, \"-800\" 25.6, \"-600\" 51.2, \"-400\" 76.8, \"-200\" 102.4, \"0\" 128, \"200\" 153.6, \"400\" 179.2, \"600\" 204.8, \"800\" 230.4, \"1000\" 256)\n");
  fprintf(fout, "set grid\n");
  fprintf(fout, "set terminal png large\n");
  fprintf(fout, "set output 'hist_%4d%02d%02d%02d%02d%02d.png'\n", yr, mon, day, hr, min, sec);
  fprintf(fout, "set multiplot\n");
  fprintf(fout, "set size 0.5, 1\n");
  
  double mean = 0.0, var;
  power = 0.0;
  for (i = 0; i < sps; i += 2) // channel A
  {
    power += (pow((double)buf[i], 2.0) / 25000000.0); // n210 range: [80-FF:00-7F] ->(char)-> [-128:127] -> -1.0V ~ +1.0V
    mean += ((double)buf[i] / 25000000.0);
    res[(int)buf[i] + 128]++;
  }
  var = sqrt(power);
  FILE *fo = fopen("hist_A.txt", "w");
  for (i = 0; i < 256; i++)
    fprintf(fo, "%3d %.12lf %.12lf\n", i, (double)res[i] * 2.0 / (double)sps, exp(-pow((double)(i - 128) - mean, 2.0) / 2.0 / power) / var / sqrt(2.0 * 3.14159));
  fclose(fo);
  
  fprintf(fout, "set origin 0, 0\n");
  fprintf(fout, "set title 'Ch. A Power %7.3lf dBm'\n", 10.0 * log10(power * 1000.0 / (128.0 * 128.0 * 50.0)));
  fprintf(fout, "plot [51:205] 'hist_A.txt' using 1:2 t'mean %7.3lf mV' with boxes,\\\n", mean * 1000.0 / 128.0);
  fprintf(fout, "'hist_A.txt' using 1:3 t'std %7.3lf mV' with lines lt 3\n", var * 1000.0 / 128.0);
  
  mean = 0.0, var = 0.0; power = 0.0;
  memset(res, 0x0, sizeof(int) * 256);
  for (i = 1; i < sps; i += 2) // channel B
  {
    power += (pow((double)buf[i], 2.0) / 25000000.0); // n210 range: [80-FF:00-7F] ->(char)-> [-128:127] -> -1.0V ~ +1.0V
    mean += ((double)buf[i] / 25000000.0);
    res[(int)buf[i] + 128]++;
  }
  var = sqrt(power);
  fo = fopen("hist_B.txt", "w");
  for (i = 0; i < 256; i++)
    fprintf(fo, "%3d %.12lf %.12lf\n", i, (double)res[i] * 2.0 / (double)sps, exp(-pow((double)(i - 128) - mean, 2.0) / 2.0 / power) / var / sqrt(2.0 * 3.14159));
  fclose(fo);
  fprintf(fout, "set ylabel ''\n");
  fprintf(fout, "set origin 0.5, 0\n");
  fprintf(fout, "set title 'Ch. B Power %7.3lf dBm'\n", 10.0 * log10(power * 1000.0 / (128.0 * 128.0 * 50.0)));
  fprintf(fout, "plot [51:205] 'hist_B.txt' using 1:2 t'mean %7.3lf mV' with boxes,\\\n", mean * 1000.0 / 128.0);
  fprintf(fout, "'hist_B.txt' using 1:3 t'std %7.3lf mV' with lines lt 3\n", var * 1000.0 / 128.0);

  fprintf(fout, "set nomultiplot\n");
  fprintf(fout, "set output\n");
  fclose(fout);
  sprintf(cmd, "gnuplot plt.txt");
  system(cmd);
  printf("hist_%4d%02d%02d%02d%02d%02d.png saved\n", yr, mon, day, hr, min, sec);
  free(buf);
  //free(ttag);
  return 0;
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
