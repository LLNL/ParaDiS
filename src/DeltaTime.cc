#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "DeltaTime.h"
#include "UTC_Time.h"

//------------------------------------------------------------------------------------------------------------
// This file implements the deltaTime debug output that can be saved by the ParaDiS application.
//------------------------------------------------------------------------------------------------------------

static void File_Save (const char *path, const char *str)
{
   if (path && str)
   {
      FILE *fd = (FILE *) ( path ? fopen(path,"w") : 0 );

      if (!fd) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, path ); }

      if (fd)
      {
         fwrite(str,strlen(str),1,fd);
         fclose(fd);
      }
   }
}

// (alternate form)


//------------------------------------------------------------------------------------------------------------
// dtime_gnu - this sequence of GNUplot commands is saved to an output file in the debug directory.
// It can be used to visualize the deltaTime data file. 
//------------------------------------------------------------------------------------------------------------

static const char *dtime_gnu = 
{
   "#set term x11 font  'helvetica,12'\n"
   "set term pngcairo size 800,480\n"
   "set title font 'helvetica,14'\n"
   "\n"
   "#-----------------------------------------------------------------------\n"
   "reset\n"
   "set out 'deltaTime_simt.png'\n"
   "unset key\n"
   "#set key box lt -1 lw 1\n"
   "#set key bottom right\n"
   "set xlabel 'Simulation Cycle'\n"
   "set ylabel 'Simulation Time (sec)'\n"
   "set logscale y\n"
   "set format   y '1e%L'\n"
   "set title 'ParaDiS Simulation Time'\n"
   "plot 'deltaTime.dat' u 1:2 w l lt 1\n"
   "\n"
   "#-----------------------------------------------------------------------\n"
   "reset\n"
   "set out 'deltaTime_dt.png'\n"
   "unset key\n"
   "#set key box lt -1 lw 1\n"
   "#set key bottom right\n"
   "set logscale y\n"
   "set format   y '1e%L'\n"
   "set xlabel 'Simulation Cycle'\n"
   "set ylabel 'Simulation dT (secs)'\n"
   "set title 'ParaDiS Timestep'\n"
   "plot 'deltaTime.dat' u 1:3 w l lt  1 notitle               , \\\n"
   "     'deltaTime.dat' u 1:3 w l lt -1 notitle smooth sbezier\n"
   "\n"
   "#-----------------------------------------------------------------------\n"
   "reset\n"
   "set out 'deltaTime_cpu.png'\n"
   "unset key\n"
   "#set key box lt -1 lw 1\n"
   "#set key bottom right\n"
   "set xlabel 'Simulation Cycle'\n"
   "set ylabel 'Cycle CPU Time (secs/cycle)'\n"
   "set title 'ParaDiS Cycle CPU Timing'\n"
   "plot 'deltaTime.dat' u 1:5 w l lt  1 notitle               , \\\n"
   "     'deltaTime.dat' u 1:5 w l lt -1 notitle smooth sbezier\n"
   "\n"
   "#-----------------------------------------------------------------------\n"
   "reset\n"
   "set out 'deltaTime_simc.png'\n"
   "unset key\n"
   "#set key box lt -1 lw 1\n"
   "#set key bottom right\n"
   "set xdata time\n"
   "set xlabel 'Simulation CPU Time (hours)'\n"
   "set ylabel 'Simulation Time (sec)'\n"
   "set logscale y\n"
   "set format   y '1e%L'\n"
   "set title 'ParaDiS Accumulated Simulation Time'\n"
   "plot 'deltaTime.dat' u ($4/60.0):2 w l lt 1\n"
   "\n"
   "quit\n"
};

static const char *dtime_desc = 
{
   "Column    Description\n"
   "------    -------------------------------------------\n"
   "1         simulation cycle number\n"
   "2         simulation time      (seconds)\n"
   "3         cycle time           (seconds)\n"
   "4         cycle deltaTT        (seconds)\n"
   "5         simulation cpu time  (seconds)\n"
   "6         cycle  time          (seconds)\n"
};

static const char *tstep_desc = 
{
      "Column    Description\n"
      "------    -------------------------------------------\n"
      "1         simulation cycle number\n"
      "2         domain index   \n"
      "3         node   index   \n"
      "4         node   position (X) \n"
      "5         node   position (Y) \n"
      "6         node   position (Z) \n"
};

//------------------------------------------------------------------------------------------------------------

typedef struct DeltaTime
{
   int       cycle    ;
   real8     timeNow  ;
   real8     deltaTT  ;
   long long cpu_dt   ;
   long long cpu_tot  ;
} DeltaTime;

static DeltaTime *dtime             = 0;     //< dtime             = points to an array of static deltaTime structs
static int        dtime_n           = 0;     //< dtime_n           = size of the deltaTime array
static long long  dtime_cycle_start = 0LL;   //< dtime_cycle_start = UTC time of application start
static long long  dtime_cycle_prev  = 0LL;   //< dtime_cycle_prev  = UTC time of previous cycle
static long long  dtime_cycle_curr  = 0LL;   //< dtime_cycle_prev  = UTC time of current cycle

//------------------------------------------------------------------------------------------------------------
// DeltaTime_Init()
// 
// Will initialize the deltaTime buffer and save the gnuplot and column description files
//------------------------------------------------------------------------------------------------------------

static void DeltaTime_Init (const Home_t *home)
{
   Param_t *param = ( home ? home->param : 0 );

   if (home && param && (home->myDomain==0) && (param->savedt==1) )
   {
      dtime_cycle_start = Current_UTC_Time();
      dtime_cycle_prev  = dtime_cycle_start;
      dtime_cycle_curr  = dtime_cycle_start;

      // save the gnuplot command and column descriptor files...

      File_Save("debug/deltaTime.gnu" , dtime_gnu );
      File_Save("debug/deltaTime.desc", dtime_desc);

      // allocate space to buffer the delta time structures...

      dtime_n = param->savedtfreq;
      dtime   = (DeltaTime *) ( (dtime_n>0) ? calloc(1,dtime_n*sizeof(DeltaTime)) : 0 );
   }
}

//------------------------------------------------------------------------------------------------------------
// DeltaTime_Append()
//
// Given an array of deltaTime structures, will append N of them to the deltaTime save file
//------------------------------------------------------------------------------------------------------------

static void DeltaTime_Append (DeltaTime *dt, const int n)
{
   if (dt && (n>0))
   {
      FILE *fd = fopen("debug/deltaTime.dat","a");
   
      if (fd)
      {
         for (int i=0; (i<n); i++)
            fprintf(fd," %8d  %12.4e %12.4e %9.2lf %6.2lf\n", 
               dt[i].cycle    ,
               dt[i].timeNow  ,
               dt[i].deltaTT  ,
               (double) (dt[i].cpu_tot)/1.0e6,
               (double) (dt[i].cpu_dt )/1.0e6 );
   
         fclose(fd);
      }
   }
}

//------------------------------------------------------------------------------------------------------------
// DeltaTime_Save()
//
// Manages the output to the deltaTime debug output. Note that this output is saved for every time step
// irrespective of what the savedt save frequency is set. Static values are used to determine if this is
// the first time called and if so - the savedt data structures are initialized.  Once initialized, individual
// cycle values are stored on the array. Output is only appended to the output file when the buffer is full.
//
// Note that you really shouldn't be saving deltaTime output for production runs - the resulting files can 
// get rather large - (BMW)
//------------------------------------------------------------------------------------------------------------

void DeltaTime_Save (const Home_t *home)
{
   Param_t *param = ( home ? home->param : 0 );

   if ( home && param && (home->myDomain==0) && (param->savedt==1) )
   {
      // If first time through, initialize the array

      static int init=1;
   
      if (init) { DeltaTime_Init(home); init=0; }

      dtime_cycle_curr = Current_UTC_Time();

      // accumulate the current values to the delta time array...

      static int indx=0;

      if (dtime) 
      {
         dtime[indx].cycle     = home->cycle     ;
         dtime[indx].timeNow   = param->timeNow  ;
         dtime[indx].deltaTT   = param->deltaTT  ;
         dtime[indx].cpu_dt    = dtime_cycle_curr - dtime_cycle_prev ;
         dtime[indx].cpu_tot   = dtime_cycle_curr - dtime_cycle_start;
         indx++;
      }

      dtime_cycle_prev = dtime_cycle_curr;

      // Only save the values out when the buffer is full...

      if (dtime && (indx==dtime_n))
      {
         DeltaTime_Append(dtime,dtime_n);
         indx=0;
      }
   }
}

//------------------------------------------------------------------------------------------------------------
// DeltaTime_Append_Timestep()
//
// Whenever the timestep has been cut, record the node ID and position where it was cut.

static void Timestep_Cut_Save_GNU (const Home_t *home, const char *path)
{
   Param_t *param = ( home ? home->param : 0 );

   FILE *fd = (FILE *) ( path ? fopen(path,"w") : 0 );

   if (param && fd)
   {
      fprintf(fd,"#set term x11 font  'helvetica,12'\n" );
      fprintf(fd,"#reset\n" );
      fprintf(fd,"#set view 60,30\n" );
      fprintf(fd,"unset key\n" );
      fprintf(fd,"set xlabel 'X' \n" );
      fprintf(fd,"set ylabel 'Y' \n" );
      fprintf(fd,"set zlabel 'Z' \n" );
      fprintf(fd,"set ticslevel 0 \n" );
      fprintf(fd,"set xrange[%4d:%4d]\n", (int) -(param->Lx/2.0), (int) (param->Lx/2.0) );
      fprintf(fd,"set yrange[%4d:%4d]\n", (int) -(param->Ly/2.0), (int) (param->Ly/2.0) );
      fprintf(fd,"set zrange[%4d:%4d]\n", (int) -(param->Lz/2.0), (int) (param->Lz/2.0) );
      fprintf(fd,"\n" );
      fprintf(fd,"stats 'timeStep.dat' nooutput\n" );
      fprintf(fd,"n = int(STATS_records) - 1\n" );
      fprintf(fd,"\n" );
      fprintf(fd,"do for [i=1:n]                                    \\\n" );
      fprintf(fd,"{                                                 \\\n" );
      fprintf(fd,"   splot 'timeStep.dat' every ::::i u 4:5:6                             w points  pt 6 ps 0.5 lc rgb 'blue' , \\\n" );
      fprintf(fd,"         'timeStep.dat' every ::::i u 4:5:6:(sprintf('(%%d,%%d)',$2,$3))  w labels  offset char 3,0           ; \\\n" );
      fprintf(fd,"   pause 0.1;  \\\n" );
      fprintf(fd,"}\n" );
      fprintf(fd,"\n" );
      fprintf(fd,"#quit\n" );

      fclose(fd);
   }
}

//------------------------------------------------------------------------------------------------------------

void Timestep_Cut_Save (const Home_t *home, Node_t *node)
{
   static int init=1;

   Param_t *param = ( home ? home->param : 0 );

   if ( home && param && (home->myDomain==0) && (param->savedt==1) )
   {
      if (init) { File_Save("debug/timeStep.desc", tstep_desc); 
                  Timestep_Cut_Save_GNU(home, "debug/timeStep.gnu");
                  init=0; }

      if (node)
      { 
         FILE *fd = fopen("debug/timeStep.dat","a");
      
         if (fd)
         {
            fprintf(fd,"%4d "                     , home->cycle );
            fprintf(fd,"%3d %-4d "                , node->myTag.domainID, node->myTag.index );
            fprintf(fd,"%15.8lf %15.8lf %15.8lf\n", node->x, node->y, node->z );
            fclose (fd);
         }
      }
   }
}
