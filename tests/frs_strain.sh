#!/bin/bash 

rm -rf frs_strain*.log frs_strain*.rslts  slurm*.out

mode="serial"
mode="parallel"

exe=../bin/paradis
dat=frs_strain.data

if [ $mode = "serial" ]; then
   echo "executing serial tests..."
   ( $exe -d $dat frs_strain_trap.ctrl      | tee -a frs_strain_trap.log       )
   ( $exe -d $dat frs_strain_arkode_aa.ctrl | tee -a frs_strain_arkode_aa.log  )
   ( $exe -d $dat frs_strain_arkode_nk.ctrl | tee -a frs_strain_arkode_nk.log  )
   ( $exe -d $dat frs_strain_kinsol_aa.ctrl | tee -a frs_strain_kinsol_aa.log  )
   ( $exe -d $dat frs_strain_kinsol_nk.ctrl | tee -a frs_strain_kinsol_nk.log  )
fi

if [ $mode = "parallel" ]; then
   echo "executing parallel tests..."
   ( mpirun -n 8 $exe -d $dat frs_strain_trap.ctrl      | tee -a frs_strain_trap.log       )
   ( mpirun -n 8 $exe -d $dat frs_strain_arkode_aa.ctrl | tee -a frs_strain_arkode_aa.log  )
   ( mpirun -n 8 $exe -d $dat frs_strain_arkode_nk.ctrl | tee -a frs_strain_arkode_nk.log  )
   ( mpirun -n 8 $exe -d $dat frs_strain_kinsol_aa.ctrl | tee -a frs_strain_kinsol_aa.log  )
   ( mpirun -n 8 $exe -d $dat frs_strain_kinsol_nk.ctrl | tee -a frs_strain_kinsol_nk.log  )
fi

