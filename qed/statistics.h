//
//  statistics.h
//  qed
//
//  Created by Daniel Alm on 30.11.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef qed_statistics_h
#define qed_statistics_h

typedef struct statistics_data {
  int N;
  double sum, square_sum, mean, error;
} statistics_data;

void reset_statistics_data(statistics_data *data);

void add_statistics_entry(statistics_data *data, double entry);
void calculate_statistics_data(statistics_data *data);
void calculate_statistics_array(statistics_data *data, int array_size);

void print_statistics_data(statistics_data *data, const char *name);
void print_statistics_array(statistics_data *data, const char *name, int array_size);

#endif
