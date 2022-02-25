#include "nl-partsol.h"

void WriteGnuplot(Matrix X, Matrix Y, double X_0, double X_n, int Instant,
                  int Number, char *Title) {
  /* Initialize GnuPlot*/
  FILE *gnuplotPipe = popen("gnuplot -persist", "w");
  if (gnuplotPipe == NULL) {
    printf("Error opening pipe to GNU plot. Check if you have it! \n");
    exit(0);
  }
  if (gnuplotPipe) {
    fprintf(gnuplotPipe, "set style data lines\n");
    fprintf(gnuplotPipe,
            "set terminal png nocrop enhanced size 1280,720; set output "
            "'%s_%d.png'\n",
            Title, Instant);
    fprintf(gnuplotPipe, "set title '%s in time : %d'\n", Title, Instant);
    fprintf(gnuplotPipe, "set xrange [%f:%f]; set yrange [-1.2:1.2]\n", X_0,
            X_n);
    fprintf(gnuplotPipe, "set grid \n");
    fprintf(gnuplotPipe, "plot '-' with lines lw 2 lc rgb 'black' \n");
    for (int i = 0; i < Number; i++) {
      fprintf(gnuplotPipe, "%f %f\n", X.nV[i], Y.nV[i]);
    }
    fprintf(gnuplotPipe, "e\n");
    fprintf(gnuplotPipe, "\n");
    fflush(gnuplotPipe);
    fprintf(gnuplotPipe, "exit \n");
    pclose(gnuplotPipe);
  }
}
