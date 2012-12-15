
#include <time.h>

#define DEBUGSTART(s)             \
{FILE *fd; fd=fopen("debug","w"); \
fprintf(fd,s); fclose(fd);}       


#define DEBUG(s)                  \
{FILE *fd; fd=fopen("debug","a"); \
fprintf(fd,s); fclose(fd);}       


#define CLOCK(x,s)                                              \
{clock_t start = clock();                                       \
  {x;}                                                          \
clock_t finish = clock();                                       \
double duration = (double)(finish - start) / CLOCKS_PER_SEC;    \
FILE *fd; fd=fopen("debug","a");                                \
fprintf(fd,"ExecTime [%s] : %5.3f\n",s,duration); fclose(fd);}  

long mem_free();
void mem_print(float n);
