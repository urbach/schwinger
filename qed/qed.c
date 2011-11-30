#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>
#include <string.h>

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
double g_mass = -0.23125;       //Fermion mass
double g_musqr = 0.21;

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
int no_timescales = 2;
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
    n_steps[0] = 8;
  
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
      n_steps[1] = 2;
    up_momenta[2] = &update_momenta_PF2;
    if (!n_steps[2])
      n_steps[2] = 2;
  }
  
  int i, l;
  int accepted = 0;        //Total number of accepted configurations
  int total_updates = 0;   //Total number of updates
  /* Containers for the mean plaquette, Polyakov loop and chiral condensate */
  /* and corresponding errors etc. */
  double mp = 0.0, dmp = 0.0, amp = 0.0, pl = 0.0, dpl = 0.0, apl = 0.0, cc = 0.0, dcc = 0.0, acc = 0.0; 
  double emp, epl, ecc;
  // Containers for the mean value of Delta H and exp(-Delta H)
  double dh = 0, ddh = 0, adh = 0, edh = 0;
  double expdh = 0, dexpdh = 0, aexpdh = 0, eexpdh = 0;
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
  gauge_force_sum = 0;
  PF1_force_sum = 0;
  PF2_force_sum = 0;
  for(i=0; i<g_measurements; i++)
  {
    g_cgiterations1 = 0;
    g_cgiterations2 = 0;
    /* do g_intermediate updates before measurement */
    int accepted_cur = 0;
    for (l=0; l<g_intermediate; l++)
      accepted_cur += update();
    accepted_cur += update();
    accepted += accepted_cur;
    /* Measurements */
    mp  = mean_plaquette();
    pl  = polyakov_loop();
    cc  = chiral_condensate();
    dh  = ham_old - ham;
    amp  += mp;    //Mean plaquette
    dmp  += mp*mp; // ... and its dispersion - for error estimation
    apl  += pl;    //Polyakov loop
    dpl  += pl*pl;
    acc  += cc;    //Chiral condensate
    dcc  += cc*cc;
    adh += dh;
    ddh += dh * dh;
    expdh = exp(dh);
    aexpdh += expdh;
    dexpdh += expdh * expdh;
    
    total_cgiterations1 += g_cgiterations1;
    total_cgiterations2 += g_cgiterations2;
    
    printf("\t Step %04i,\t mp = %2.4lf,\t pl = %2.4lf,\t cc = %2.4lf,\t dh = %2.4lf,\tcg1 = %d,\tcg2 = %d,\tacc = %d\n", i, mp, pl, cc, -dh, g_cgiterations1, g_cgiterations2, accepted_cur);
  };
  
  /* Some output for diagnostics */
  printf("\n\n Algorithm configuration:\n");
  printf("\t X_1 = %d, beta = %f, g_mass = %f, g_musqr = %f\n", X1, beta, g_mass, g_musqr);
  printf("\t g_measurements = %d, g_intermediate = %d, n_steps[0] = %d, n_steps[1] = %d, n_steps[2] = %d\n", g_measurements, g_intermediate, n_steps[0], n_steps[1], n_steps[2]);
  
  total_updates = g_measurements*(g_intermediate + 1);
  printf("\n\n Algorithm performance:\n");
  printf("\t Acceptance rate:                      %2.2lf\n\n", (double)accepted/(double)total_updates);
  printf("\t Inner loop CG iterations per update:  %2.2lf\n", (double)total_cgiterations1/(double)total_updates);
  printf("\t Outer loop CG iterations per update:  %2.2lf\n", (double)total_cgiterations2/(double)total_updates);
  printf("\t Inner loop CG iterations per solve:   %2.2lf\n", (double)total_cgiterations1/(double)total_updates/(n_steps[2]*n_steps[1]));
  printf("\t Outer loop CG iterations per solve:   %2.2lf\n\n", (double)total_cgiterations2/(double)total_updates/(n_steps[2]));
  printf("\t Gauge force per update:               %2.6lf\n", (double)gauge_force_sum/(double)total_updates*tau);
  printf("\t Inner loop force per update:          %2.6lf\n", (double)PF1_force_sum/(double)total_updates*tau);
  printf("\t Outer loop force per update:          %2.6lf\n", (double)PF2_force_sum/(double)total_updates*tau);
  printf("\t Gauge force per application:          %2.6lf\n", (double)gauge_force_sum/(double)total_updates*tau/(n_steps[2]*n_steps[1]*n_steps[0]));
  printf("\t Inner loop force per application:     %2.6lf\n", (double)PF1_force_sum/(double)total_updates*tau/(n_steps[2]*n_steps[1]));
  printf("\t Outer loop force per application:     %2.6lf\n\n", (double)PF2_force_sum/(double)total_updates*tau/(n_steps[2]));
  printf("\t Runtime / seconds:                    %2.2lf\n", (double)(clock() - clock_start)/(double)CLOCKS_PER_SEC);
  
  /* Measurement results */
  amp = amp/(double)g_measurements;
  apl = apl/(double)g_measurements;
  acc = acc/(double)g_measurements;
  adh = adh/(double)g_measurements;
  aexpdh = aexpdh/(double)g_measurements;
  dmp = dmp/(double)g_measurements - amp*amp;
  dpl = dpl/(double)g_measurements - apl*apl;
  dcc = dcc/(double)g_measurements - acc*acc;
  ddh = ddh/(double)g_measurements - adh*adh;
  dexpdh = dexpdh/(double)g_measurements - aexpdh*acc;
  emp = sqrt(dmp/(double)(g_measurements - 1)); //Errors of corresponding mean values
  epl = sqrt(dpl/(double)(g_measurements - 1));
  ecc = sqrt(dcc/(double)(g_measurements - 1));
  edh = sqrt(ddh/(double)(g_measurements - 1));
  eexpdh = sqrt(dexpdh/(double)(g_measurements - 1));
  
  printf("\n\n Mean values:\n");
  printf("\t Plaquette:                  %2.6lf +/- %2.6lf (e = %02i%%)\n", amp, emp, (int)(100.0*emp/amp));
  printf("\t Polyakov loop:              %2.6lf +/- %2.6lf (e = %02i%%)\n", apl, epl, (int)(100.0*epl/apl));
  printf("\t Chiral Condensate:          %2.6lf +/- %2.6lf (e = %02i%%)\n", acc, ecc, (int)(100.0*ecc/acc));
  printf("\t -Delta H:                   %2.6lf +/- %2.6lf (e = %02i%%)\n", adh, edh, (int)(100.0*edh/adh));
  printf("\t exp(-Delta H):              %2.6lf +/- %2.6lf (e = %02i%%)\n", aexpdh, eexpdh, (int)(100.0*eexpdh/aexpdh));
  printf("\n\n");
  
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
