#include <stdlib.h>
#include <stdio.h>
#define NUM_POINTS 5
#define NUM_COMMANDS 2

int main()
{
    char * commandsForGnuplot[] = {"set title \"TITLEEEEE\"", "plot 'data.temp'"};
    double xvals[NUM_POINTS] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double yvals[NUM_POINTS] = {5.0 ,3.0, 1.0, 3.0, 5.0};
    FILE * temp = fopen("data.temp", "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    for (i=0; i < NUM_POINTS; i++)
    {
    fprintf(temp, "%lf %lf \n", xvals[i], yvals[i]); //Write the data to a temporary file
    }

    for (i=0; i < NUM_COMMANDS; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    return 0;
}

fprintf(gnuplotPipe, "plot '-' \n");
int i;

for (int i = 0; i < NUM_POINTS; i++)
{
  fprintf(gnuplotPipe, "%lf %lf\n", xvals[i], yvals[i]);
}

fprintf(gnuplotPipe, "e");


int frame = 0;
event movie(i += 2){
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot '-' w l lw 5 t 'C', '-' w l lw 3 t 'CC'\n");
  foreach()
    fprintf(gnuplotPipe, "%g %g\n",x, C[]);
  fprintf(gnuplotPipe, "e\n");
  foreach()
    fprintf(gnuplotPipe, "%g %g\n",x, CC[]);
  fprintf(gnuplotPipe, "e\n");
  frame++;
}
