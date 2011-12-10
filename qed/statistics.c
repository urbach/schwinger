//
//  statistics.c
//  qed
//
//  Created by Daniel Alm on 30.11.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>

#include <string.h>
#include <math.h>

#include "statistics.h"

void reset_statistics_data(statistics_data *data)
{
  memset(data, 0, sizeof(statistics_data));
}

void add_statistics_entry(statistics_data *data, double entry)
{
  data->N += 1;
  data->sum += entry;
  data->square_sum += entry * entry;
}

void calculate_statistics_data(statistics_data *data)
{
  int N = data->N;
  data->mean = data->sum / N;
  data->error = sqrt((data->square_sum / N - data->mean * data->mean) / (N - 1));
}

void calculate_statistics_array(statistics_data *data, int array_size)
{
  for (int i = 0; i < array_size; i ++)
    calculate_statistics_data(&data[i]);
}

void print_statistics_data(statistics_data *data, const char *name, double factor)
{
  calculate_statistics_data(data);
  printf("\t %-38s%2.6lf +/- %2.6lf (e = %02i%%)\n", name, data->mean * factor, data->error * factor, (int)(100.0*data->error/data->mean));
}

void print_statistics_array(statistics_data *data, const char *name, int array_size, double factor)
{
  for (int i = 0; i < array_size; i ++)
  {
    char name_buffer[1000];
    sprintf(name_buffer, "%s[%i]:", name, i);
    print_statistics_data(&data[i], name_buffer, factor);
  }
}

double autocorrelation_Gamma(double *measurements, double mean, int N, int n)
{
  double result = 0;
  // if n is negative, we have to start the loop at a later point
  int loopstart = (n > 0 ? 0 : -n);
  int loopend = (n > 0 ? N - n : N);
  for (int i = loopstart; i < loopend; i ++)
    result += (measurements[i] - mean) * (measurements[i + n] - mean);
  
  return result / (loopstart - loopend);
}

double autocorrelation_time(double *measurements, int N)
{
  double mean = 0;
  for (int i = 0; i < N; i ++)
    mean += measurements[i];
  mean /= N;
  
  double Gamma0 = 0.5 * autocorrelation_Gamma(measurements, mean, N, 0);
  double result = Gamma0;
  
  //printf("%g\n", Gamma0);
  for (int i = 1; i < N; i ++)
  {
    result += autocorrelation_Gamma(measurements, mean, N, i);
    //printf("%g\n", autocorrelation_Gamma(measurements, mean, N, i));
  }
  
  return result / Gamma0;
}
