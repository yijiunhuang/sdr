#include "sample_v170215.hpp"

// sampling subroutine that should be called by a thread
void *sampling(void *arg)
{
  samplebuffer *sb = (samplebuffer *)arg;
  
	// variables to be set by "boost::program_options"
  std::string args = "addr=192.168.10.2", subdev = "A:AB";
  
	// create a usrp device
	uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
	
	// always select the subdevice first, the channel mapping affects the other settings
	usrp->set_rx_subdev_spec(subdev); // A:AB, A:BA, A:A, A:B
  usrp->set_clock_source("external"); // 10 MHz source : internal, external
	usrp->set_time_source("external"); // PPS source: none, external, _external_, mimo
	usrp->set_rx_rate(sps);
	usrp->set_rx_freq(0.0); // set the center frequency = 0 : intermediate frequecy sampling
	usrp->set_rx_gain(0.0);
	usrp->set_rx_bandwidth(2.5e+8); // fixed for Basic RX: 250.0MHz
	
	// create a receive streamer
	// stream_args(host sample format, device wire format)
	// format: sc = short complex, fc = float complex
	// host format: fc64 = complex double, fc32 = complex float, sc16 = complex short (map to +-32768)
	// over-the-wire format: sc8(up to 50Msps for host) or sc16(up to 25Msps)
	uhd::stream_args_t stream_args("sc8", "sc8");
	uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
	uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

  uhd::rx_metadata_t md; // create a RX metadata
  char *buff; // create a buffer
  size_t num;
  double dt;
  timetag *ttag;
  struct timeval pctime;
  struct timezone tz;
  time_t current_time, lasttime;
  struct tm *timeinfo;
  uhd::time_spec_t boradtime;
  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);

  timeinfo = (struct tm *)malloc(sizeof(struct tm));
  buff = (char *)malloc(sizeof(char) * sps * 2); // previous, IQIQ..
  ttag = (timetag *)malloc(sizeof(timetag));

  if (sb->is_tt == 0x1) // set time stamp manually
    current_time = sb->tt + 1;
  else
  {
    lasttime = 0;
    while (1) // get time from PC
    {
      time(&current_time);
      if (lasttime == 0) lasttime = current_time;
      if (lasttime != current_time) break;
    }
    printf("set PC time stamp to the board ...\n");
  }
  while (1)
  {
    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS;// set continous streaming
    stream_cmd.num_samps = 0;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = uhd::time_spec_t(current_time + 2); // streaming samples upon the second PPS pulse
    usrp->set_time_next_pps(uhd::time_spec_t(current_time + 1)); // the time stamp that will be set to the board
    sleep(1); // set the time stamp to the board
    usrp->issue_stream_cmd(stream_cmd); // start sampling
    while (1)
    {
      num = rx_stream->recv(buff, sps, md);
      if (md.error_code == 0x0)
      {
        if (num == sps)
        {
          // take the time difference PC - board
          boradtime = usrp->get_time_now();// get onboard time
          gettimeofday(&pctime, &tz);
          dt = (double)(pctime.tv_sec - boradtime.get_full_secs()) + (double)pctime.tv_usec * 1.0e-6 - boradtime.get_frac_secs();
          // display time stamp onboard
          current_time = boradtime.get_full_secs();
          timeinfo = gmtime(&current_time);
          if (timeinfo->tm_min == 0 && timeinfo->tm_sec == 7) // when reset? xx:00:07
          {
            if (sb->is_tt == 0x1) // manual mode
            {
              printf("%4d/%02d/%02d %02d:%02d:%02d UTC (manual)\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
              printf("reset sampling\n");
            }
            else // correct time stamp according to NTP
            {
              printf("%4d/%02d/%02d %02d:%02d:%02d UTC\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
              if (fabs(fmod(dt, 1.0)) < 0.1 || fabs(fmod(dt, 1.0)) > 0.9)
              {
                current_time += (int)round(dt);
                if ((int)round(dt) == 0) printf("reset sample streaming\n");
                else printf("reset sampling, and implement %+d s\n", (int)round(dt));
              }
            }
            stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
            usrp->issue_stream_cmd(stream_cmd);
            break;
          }
          else
          {
            if (sb->is_tt == 0x1) // manual mode
              printf("%4d/%02d/%02d %02d:%02d:%02d UTC (manual)\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
            else
            {
              if (fabs(dt) < 0.1) // say "NTP is alive"
                printf("%4d/%02d/%02d %02d:%02d:%02d UTC\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
              else if (fabs(fmod(dt, 1.0)) < 0.1 || fabs(fmod(dt, 1.0)) > 0.9) // say "NTP is alive and time stamp will be corrected"
                printf("%4d/%02d/%02d %02d:%02d:%02d UTC, a %+d s will be implemented\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, (int)round(dt));
              else // say "NTP is failed", and display alert
                printf("%4d/%02d/%02d %02d:%02d:%02d UTC, PC - board = %9.3lf s, please check NTP\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, dt);
            }
          }
          // write samples and time stamp to the buffer
          current_time = md.time_spec.get_full_secs();
          timeinfo = gmtime(&current_time);
          ttag->year = timeinfo->tm_year + 1900;
          ttag->month = timeinfo->tm_mon + 1;
          ttag->day = timeinfo->tm_mday;
          ttag->hour = timeinfo->tm_hour;
          ttag->minute = timeinfo->tm_min;
          ttag->second = timeinfo->tm_sec;
          ttag->fsec = md.time_spec.get_frac_secs();
          pthread_mutex_lock(&sb->mymutex);
          sb->writed++;
          while (sb->writed - sb->readd > MAXBUF)
            sb->readd += MAXBUF;
          memcpy(sb->buf[sb->writed % MAXBUF], buff, sizeof(char) * sps * 2);
          memcpy(sb->ttag[sb->writed % MAXBUF], ttag, sizeof(timetag));
          if (sb->readd == sb->writed)
            pthread_cond_signal(&sb->mycond);
          pthread_mutex_unlock(&sb->mymutex);
        }
      }
      else if (md.error_code == 0x1) //  ERROR_CODE_TIMEOUT
        continue;
      else
      {
        time(&current_time);
        stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
        usrp->issue_stream_cmd(stream_cmd);
        break;
      }
    }
  }
  pthread_exit(NULL);
}
