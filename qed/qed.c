#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>

#include "statistics.h"
#include "lattice.h"
#include "hmc.h"
#include "fields.h"
#include "measurements.h"
#include "2MN_integrator.h"

/* global variables */
int g_thermalize = 500;      //Number of MC updates for thermalization
int g_measurements = 1000;  //Number of measurements (statistically independent configurations)
int g_intermediate = 0;    //Number of MC updates between the measurements
/* extern in hmc.h      */
double ham = 0, ham_old = 0;
/* The critical mass at beta=1.0 is about 0.32 */
double g_mass = -0.06;//3125;//23125;       //Fermion mass
double g_musqr = 0.7;

double beta = 1.0;       //Coupling constant for the gauge field

double thermalize_min_acc = 0.7; // minimum ratio of accepted thermalization steps

void echo_sim_params();
void save_gauge(const char *filename);
void load_gauge(const char *filename);

/* The following we need for the recursive integration scheme */
/* the integration step numbers for the different time scales */
int n_steps[3];
/* list of function pointers to the momentum update functions */
up_m up_momenta[3];
int no_timescales = 3;
double tau = 1.;

int main(int argc, char **argv) 
{
  clock_t clock_start = clock();
  
  n_steps[0] = 0;
  n_steps[1] = 0;
  n_steps[2] = 0;
  
  struct option long_options[] =
  {
    /* These options set a flag. */
    {"thermalize", required_argument, NULL, 0},
    {"measurements", required_argument, NULL, 0},
    {"intermediate", required_argument, NULL, 0},
    {"mass", required_argument, NULL, 0},
    {"musqr", required_argument, NULL, 0},
    {"beta", required_argument, NULL, 0},
    {"n_steps_0", required_argument, NULL, 0},
    {"n_steps_1", required_argument, NULL, 0},
    {"n_steps_2", required_argument, NULL, 0},
    {"tau", required_argument, NULL, 0},
    {"thermalize_min_acc", required_argument, NULL, 0},
    {"no_timescales", required_argument, NULL, 0},
    {0, 0, 0, 0}
  };
  
  while (1)
  {
    int option_index = 0;
    
    int c = getopt_long(argc, argv, "", long_options, &option_index);
    
    if (c == -1)
      break;
    
    if (c != 0
        || !optarg)
      continue;
    
    const char *optionName = long_options[option_index].name;
    double optionDoubleValue = strtod(optarg, NULL);
    printf("%s = %f\n", optionName, optionDoubleValue);
    
    if (strcmp(optionName, "thermalize") == 0)
      g_thermalize = optionDoubleValue;
    else if (strcmp(optionName, "measurements") == 0)
      g_measurements = optionDoubleValue;
    else if (strcmp(optionName, "intermediate") == 0)
      g_intermediate = optionDoubleValue;
    else if (strcmp(optionName, "mass") == 0)
      g_mass = optionDoubleValue;
    else if (strcmp(optionName, "musqr") == 0)
      g_musqr = optionDoubleValue;
    else if (strcmp(optionName, "beta") == 0)
      beta = optionDoubleValue;
    else if (strcmp(optionName, "n_steps_0") == 0)
      n_steps[0] = optionDoubleValue;
    else if (strcmp(optionName, "n_steps_1") == 0)
      n_steps[1] = optionDoubleValue;
    else if (strcmp(optionName, "n_steps_2") == 0)
      n_steps[2] = optionDoubleValue;
    else if (strcmp(optionName, "tau") == 0)
      tau = optionDoubleValue;
    else if (strcmp(optionName, "thermalize_min_acc") == 0)
      thermalize_min_acc = optionDoubleValue;
    else if (strcmp(optionName, "no_timescales") == 0)
      no_timescales = optionDoubleValue;
  }
  
  // setup integration parameters
  /* most inner steps per next outer step  */
  /* this is typically the gauge action    */
  if (!n_steps[0])
    n_steps[0] = 12;
  
  /* The gauge action is always done on the most inner timescale */
  /* these are defined in 2MN_integrator.h                       */
  up_momenta[0] = &update_momenta_gauge;
  if (no_timescales < 3)
  {
    up_momenta[1] = &update_momenta_PF1;
    g_musqr = 0.;
    if (!n_steps[1])
      n_steps[1] = 20;
  }
  else
  {
    up_momenta[1] = &update_momenta_PF1;
    if (!n_steps[1])
      n_steps[1] = 1;
    up_momenta[2] = &update_momenta_PF2;
    if (!n_steps[2])
      n_steps[2] = 6;
  }
  
  int i, l;
  int accepted = 0;        //Total number of accepted configurations
  int total_updates = 0;   //Total number of updates
  /* Containers for the mean plaquette, Polyakov loop and chiral condensate */
  /* and corresponding errors etc. */
  double total_cgiterations1 = 0, total_cgiterations2 = 0;
  /* Initialize the random number generator */
  rlxd_init(2, 123456); 
  /* Initialize the lattice geometry */
  init_lattice(X1, X2);
  /* Initialize the fields */
  coldstart();
  /* Print out the run parameters */
  echo_sim_params(); 
  
  /* thermalization */
  printf("\n Thermalization: \n\n");
  char dump_filename[1000];
  sprintf(dump_filename, "dumps/beta=%lf mass=%lf X1=%i X2=%i.dmp", beta, g_mass, X1, X2);
  
  if (access(dump_filename, F_OK) != -1)
  {
    // dump exists, so just load it
    load_gauge(dump_filename);
    calculatelinkvars();
  }
  else
  {
    for(i=0; i<g_thermalize; i++)
    {
      accepted += update();
      printf("\t Step %04i,\t mp = %2.4lf,\t pl = %2.4lf,\t cc = %2.4lf\t dh = %2.4lf,\tcg1 = %d,\tcg2 = %d\n", i, mean_plaquette(), polyakov_loop(), chiral_condensate(), ham-ham_old, g_cgiterations1, g_cgiterations2);
      g_cgiterations1 = 0;
      g_cgiterations2 = 0;
    }
    if (accepted < g_thermalize * thermalize_min_acc)
    {
      printf("Acceptance ratio %0.2f too low, aborting...\n", accepted / (double)g_thermalize);
      exit(1);
    }
    save_gauge(dump_filename);
  }
  accepted = 0;
  
  // re-initialize the random number generator,
  // because the thermalization may have been omitted
  rlxd_init(2, 123456);
  
  /* measure the iterations only during real simulation, not thermalization */
  R = 0;               //Total number of rejected configurations
  
  printf("\n Generation: \n\n");
  
  statistics_data gauge_force_statistics; reset_statistics_data(&gauge_force_statistics);
  statistics_data PF1_force_statistics; reset_statistics_data(&PF1_force_statistics);
  statistics_data PF2_force_statistics; reset_statistics_data(&PF2_force_statistics);
  statistics_data mp_statistics; reset_statistics_data(&mp_statistics);
  statistics_data pl_statistics; reset_statistics_data(&pl_statistics);
  statistics_data cc_statistics; reset_statistics_data(&cc_statistics);
  statistics_data tc_statistics; reset_statistics_data(&tc_statistics);
  statistics_data dh_statistics; reset_statistics_data(&dh_statistics);
  statistics_data expdh_statistics; reset_statistics_data(&expdh_statistics);
  statistics_data pion_correlation_statistics[X2];
  statistics_data pcac_correlation_statistics[X2];
  
  double *mp_measurements = malloc(g_measurements * sizeof(double));
  double *tc_measurements = malloc(g_measurements * sizeof(double));
  for (int t = 0; t < X2; t ++)
  {
    reset_statistics_data(&pion_correlation_statistics[t]);
    reset_statistics_data(&pcac_correlation_statistics[t]);
  }
  for(i = 0; i < g_measurements; i ++)
  {
    g_cgiterations1 = 0;
    g_cgiterations2 = 0;
    gauge_force = 0;
    PF1_force = 0;
    PF2_force = 0;
    /* do g_intermediate updates before measurement */
    int accepted_cur = 0;
    for (l=0; l<g_intermediate; l++)
      accepted_cur += update();
    accepted_cur += update();
    accepted += accepted_cur;
    /* Measurements */
    double mp  = mean_plaquette();
    double pl  = polyakov_loop();
    double cc  = chiral_condensate();
    double tc  = topological_charge();
    double dh  = ham_old - ham;
    
    if (accepted_cur >= 1)
    {
      add_statistics_entry(&gauge_force_statistics, gauge_force);
      add_statistics_entry(&PF1_force_statistics, PF1_force);
      add_statistics_entry(&PF2_force_statistics, PF2_force);
    }
    
    
    double C_X2_2 = 0;
    for (int t = 0; t < X2; t ++)
    {
      double pion_correlation = pion_correlation_function(t);
      add_statistics_entry(&pion_correlation_statistics[t], pion_correlation);
      add_statistics_entry(&pcac_correlation_statistics[t], pcac_correlation_function(t));
      if (t == X2 / 2)
        C_X2_2 = pion_correlation;
    }
    
    add_statistics_entry(&mp_statistics, mp);
    add_statistics_entry(&pl_statistics, pl);
    add_statistics_entry(&cc_statistics, cc);
    add_statistics_entry(&tc_statistics, tc);
    add_statistics_entry(&dh_statistics, dh);
    add_statistics_entry(&expdh_statistics, exp(dh));
    
    mp_measurements[i] = mp;
    tc_measurements[i] = tc;
    
    total_cgiterations1 += g_cgiterations1;
    total_cgiterations2 += g_cgiterations2;
    
    printf("\t Step %04i,\t mp = %2.8lf,\t pl = %2.4lf,\t cc = %2.4lf,\t tc = %2.1lf,\t C(%i) = %2.10lf,\t dh = %2.8lf,\t cg1 = %d,\t cg2 = %d,\t acc = %d\n", i, mp, pl, cc, tc, X2 / 2, C_X2_2, -dh, g_cgiterations1, g_cgiterations2, accepted_cur);
  }
  
  /* Some output for diagnostics */
  printf("\n\n Algorithm configuration:\n");
  printf("\t X_1 = %d, beta = %f, g_mass = %f, g_musqr = %f\n", X1, beta, g_mass, g_musqr);
  printf("\t g_measurements = %d, g_intermediate = %d, n_steps[0] = %d, n_steps[1] = %d, n_steps[2] = %d\n", g_measurements, g_intermediate, n_steps[0], n_steps[1], n_steps[2]);
  
  total_updates = g_measurements*(g_intermediate + 1);
  printf("\n\n Algorithm performance:\n");
  printf("\t Acceptance rate:                      %2.2lf\n\n", (double)accepted/(double)total_updates);
  printf("\t Inner loop CG iterations per update:  %2.2lf\n", (double)total_cgiterations1/(double)total_updates);
  printf("\t Outer loop CG iterations per update:  %2.2lf\n", (double)total_cgiterations2/(double)total_updates);
  // FIXME this is wrong for 1 and 2 scales!
  printf("\t Inner loop CG iterations per solve:   %2.2lf\n", (double)total_cgiterations1/(double)total_updates/(n_steps[2]*n_steps[1] + 1));
  printf("\t Outer loop CG iterations per solve:   %2.2lf\n\n", (double)total_cgiterations2/(double)total_updates/(n_steps[2] + 1));
  printf("\t Runtime / seconds:                    %2.2lf\n", (double)(clock() - clock_start)/(double)CLOCKS_PER_SEC);
  
  print_statistics_data(&gauge_force_statistics, "Gauge force per update:", tau);
  print_statistics_data(&PF1_force_statistics, "Inner loop force per update:", tau);
  print_statistics_data(&PF2_force_statistics, "Outer loop force per update:", tau);
  print_statistics_data(&gauge_force_statistics, "Gauge force per application:", tau/(n_steps[2]*n_steps[1]*n_steps[0]));
  print_statistics_data(&PF1_force_statistics, "Inner loop force per application:", tau/(n_steps[2]*n_steps[1]));
  print_statistics_data(&PF2_force_statistics, "Outer loop force per application:", tau/(n_steps[2]));
  
  /* Measurement results */
  printf("\n\n Mean values:\n");
  print_statistics_data(&mp_statistics, "Plaquette:", 1);
  printf("\t Plaquette autocorrelation time:       %g\n", autocorrelation_time(mp_measurements, g_measurements));
  print_statistics_data(&pl_statistics, "Polyakov loop:", 1);
  print_statistics_data(&cc_statistics, "Chiral Condensate:", 1);
  print_statistics_data(&tc_statistics, "Topological Charge:", 1);
  printf("\t Topological Charge autocorrelation time:       %g\n", autocorrelation_time(tc_measurements, g_measurements));
  print_statistics_data(&dh_statistics, "-Delta H:", 1);
  print_statistics_data(&expdh_statistics, "exp(-Delta H):", 1);
  printf("\n");
  print_statistics_array(pion_correlation_statistics, "Pion Correlation", X2, 1);
  printf("\n");
  print_statistics_array(pcac_correlation_statistics, "PCAC Correlation", X2, 1);
  
  free(mp_measurements);
  free(tc_measurements);
  
  free(left1);
  free(left2);
  free(right1);
  free(right2);
  
  return 0;
}

void echo_sim_params()
{
  printf("Hybrid Monte-Carlo for two-dimensional QED with Wilson fermions\n\n");
  printf("Run parameters (adjust in qed.c !!!):\n");
  printf("\t Beta:                            %2.4lf\n",  beta);
  printf("\t Mass:                            %2.4lf\n",  g_mass);
  printf("\t Lattice size:                    %i x %i\n", X1, X2);
  printf("\t Thermalization steps:            %i\n",      g_thermalize);
  printf("\t Number of measurements:          %i\n",      g_measurements);
  printf("\t MC updates between measurements: %i\n",      g_intermediate);
  printf("\n\n");
}

void save_gauge(const char *filename)
{
  FILE *file;
  
  file = fopen(filename, "w");
  fprintf(file, "%.16lg\n\n", s_g_old);
  for (int i = 0; i < GRIDPOINTS; i ++)
    fprintf(file, "%.16lg %.16lg\n", gauge1[i], gauge2[i]);
  fclose(file);
}
void load_gauge(const char *filename)
{
  FILE *file;
  
  file = fopen(filename, "r");
  fscanf(file, "%lg", &s_g_old);
  for (int i = 0; i < GRIDPOINTS; i ++)
  {
    double a, b;
    fscanf(file, "%lg %lg", &a, &b);
    gauge1[i] = gauge1_old[i] = a;
    gauge2[i] = gauge2_old[i] = b;
  }
  fclose(file);
}
