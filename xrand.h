//
// Definitions for the random number generation
//

#ifndef __XRAND__
#define __XRAND__

void rnd_init (unsigned seed);
int rnd_i ();
unsigned rnd_u ();
int rnd_ri (int rng);
double rnd_01d ();
double rnd_ned (double lam);
double rnd_nedi (double lam);

#endif
