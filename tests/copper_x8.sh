#!/bin/bash 

rm -rf copper_x8_results
rm -rf copper_x8_maqao

# leaves paradis results in    ./copper_x8_results
# leaves performance traces in ./copper_x8_maqao
#--------------------------------------------------------------------------------

#maqao=/collab/usr/global/tools/maqao/chaos_5_x86_64_ib/maqao

# if maqao is available - run the test with maqao performance analysis tool.
# else just execute a normal MPI run of the copper test.

#if [ -x $maqao ]; then
#   mpirun -n 8 $maqao lprof xp=copper_x8_maqao -- ./bin/paradis .copper_x8.ctrl 
#
#   $maqao lprof     xp=copper_x8_maqao of=html
#   $maqao lprof -df xp=copper_x8_maqao
#   $maqao lprof -dl xp=copper_x8_maqao
#else
#   mpirun -n 8 ./bin/paradis ./tests/copper_x8.ctrl 
#fi

mpirun -n 8 ../bin/paradis copper_x8.ctrl 

# Here is an example of running the test using valgrind/memcheck...
#
# Note that valgrind seriously impacts performance. You will need to adjust
# the ctrl file maxstep to 200-500 to run the memory check or you will be
# waiting for a very long time...

#srun -ppdebug -n 8 -t 30 memcheck_all ./bin/paradis copper_x8.ctrl 
