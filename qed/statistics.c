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

void print_statistics_data(statistics_data *data, const char *name)
{
  calculate_statistics_data(data);
  printf("\t %-28s%2.6lf +/- %2.6lf (e = %02i%%)\n", name, data->mean, data->error, (int)(100.0*data->error/data->mean));
}

void print_statistics_array(statistics_data *data, const char *name, int array_size)
{
  for (int i = 0; i < array_size; i ++)
  {
    char name_buffer[1000];
    sprintf(name_buffer, "%s[%i]:", name, i);
    print_statistics_data(&data[i], name_buffer);
  }
}
