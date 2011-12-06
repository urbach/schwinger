CC = gcc -std=c99 -Werror -Wall -O3

all: qed

rec_lf_integrator.o: rec_lf_integrator.c rec_lf_integrator.h leapfrog.h Makefile
	$(CC) -c $< -o $@


fields.o: fields.c fields.h lattice.h linalg.h complex/complex.h rand/ranlxd.h Makefile
	$(CC) -c $< -o $@

hmc.o: hmc.c hmc.h lattice.h fields.h linalg.h complex/complex.h rand/ranlxd.h rand/gauss.h Makefile
	$(CC) -c $< -o $@

lattice.o: lattice.c lattice.h
	$(CC) -c $< -o $@

statistics.o: statistics.c statistics.h
	$(CC) -c $< -o $@

linalg.o: linalg.c linalg.h lattice.h complex/complex.h Makefile
	$(CC) -c $< -o $@

ranlxd.o: rand/ranlxd.c rand/ranlxd.h Makefile
	$(CC) -c $< -o $@ -I rand

gauss.o: rand/gauss.c rand/gauss.h Makefile
	$(CC) -c $< -o $@ -I rand

measurements.o: measurements.c measurements.h fields.h lattice.h linalg.h hmc.h complex/complex.h Makefile
	$(CC) -c $< -o $@

dirac.o: dirac.c dirac.h lattice.h linalg.h complex/complex.h rand/ranlxd.h Makefile
	$(CC) -c $< -o $@

leapfrog.o: leapfrog.c leapfrog.h hmc.h Makefile
	$(CC) -c $< -o $@

2MN_integrator.o: 2MN_integrator.c 2MN_integrator.h hmc.h Makefile
	$(CC) -c $< -o $@

leapfrog2.o: leapfrog2.c leapfrog2.h 2MN_integrator.h hmc.h Makefile
	$(CC) -c $< -o $@

qed.o: qed.c measurements.h fields.h lattice.h statistics.h linalg.h hmc.h complex/complex.h Makefile rand/ranlxd.h
	$(CC) -c $< -o $@

test1.o: test1.c measurements.h fields.h lattice.h linalg.h hmc.h complex/complex.h Makefile rand/ranlxd.h
	$(CC) -c $< -o $@

test2.o: test2.c measurements.h fields.h lattice.h linalg.h hmc.h complex/complex.h Makefile rand/ranlxd.h
	$(CC) -c $< -o $@

qed: fields.o qed.o hmc.o lattice.o statistics.o linalg.o measurements.o ranlxd.o gauss.o dirac.o leapfrog.o 2MN_integrator.o leapfrog2.o rec_lf_integrator.o Makefile
	$(CC) qed.o fields.o hmc.o lattice.o statistics.o linalg.o measurements.o ranlxd.o gauss.o dirac.o leapfrog.o 2MN_integrator.o leapfrog2.o rec_lf_integrator.o -o qed -lm

test1: fields.o test1.o hmc.o lattice.o linalg.o measurements.o ranlxd.o gauss.o dirac.o leapfrog.o 2MN_integrator.o Makefile
	$(CC) test1.o fields.o hmc.o lattice.o linalg.o measurements.o ranlxd.o gauss.o dirac.o leapfrog.o 2MN_integrator.o -o test1 -lm

test2: fields.o test2.o hmc.o lattice.o linalg.o measurements.o ranlxd.o gauss.o dirac.o leapfrog.o 2MN_integrator.o Makefile
	$(CC) test2.o fields.o hmc.o lattice.o linalg.o measurements.o ranlxd.o gauss.o dirac.o leapfrog.o 2MN_integrator.o -o test2 -lm

clean:
	rm -f *.o qed test1 test2

