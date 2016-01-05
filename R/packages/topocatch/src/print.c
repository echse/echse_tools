#include <R.h>

/*
void F77_SUB(pdbl)(double *x, int *newline) {
  if (*newline == 0) {
    printf("%f",*x);
  } else {
    printf("%f\n",*x);
  }
}
*/

void F77_SUB(printstate)(double *x) {
  printf("[1] Done: %0.0f%%\n",*x);
}

