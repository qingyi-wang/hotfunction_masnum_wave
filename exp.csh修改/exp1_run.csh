#!/bin/tcsh -f
#-------------------------------------------------------------------------------
# ************ Part I. Set parameters for the MASNUM_WAVE model. ***************
#-------------------------------------------------------------------------------

# --- If need to make the executable file, set needmake as "YES", else "NO".
set needmake = "YES"

# --- Set number of processors for MPI run.
set nproc = 24

# --- Path for wave model running.
set masnum_home = $HOME/online1/masnum_wave

# --- File of depth for wave model.
set depfile = $masnum_home/inputdata/topo/topo_pacific_p15.nc

# --- File of ice mask for wave model.
#set icefile = $masnum_home/inputdata/topo/ice_clim_mask.nc

# --- Path for code, w/subpath: scripts, wave_cor & pre_time
set src_path = $masnum_home/source

# --- Path for wind & model setting.
#     NOTE: end with slash; keep it agree with windtype.
set wind_path = ./forcing/

#-------------------------------------------------------------------------------

set title      = "pac_ncep"    # --- Symbal for model output.
set istime     = 20090101      # --- Integral start time
set ietime     = 20090121      # --- Integral end time
set cools_days = 0             # --- The time (days) for cool start.
set delttm     = 7.5           # --- Length of integral time step, in minutes.
                               # --- Maximum value is 7.5 in this example

set wndfreq = 6   # --- The frequence of wind data (hours).
set wndtype = 3   # --- The wind type:
                  #  0  for wind in the same grid with model, files by monthly.
                  #  1  for GFS wind (0.5 * 0.5), no interp.
                  #  2  for QuikSCAT BLN wind (0.5 * 0.5), interp.
                  #  3  for NCEP re-anal wind, with interp.
set outflag = 3   #  output wave variables into file multi-records,
                  #  1    : one file every year,
                  #  2    : one file every month,
                  #  3    : one file every day,
                  #  else : one file every run.
set wiofreq = 24  # --- The output frequence for wave results (hour).
set ciofreq = 24  # --- The output frequence for current coef.s (hour).
set rstfreq = 24  # --- The output frequence for model restart (hour).

#-------------------------------------------------------------------------------
# ************ PART II. Prepare work directory, executable files, **************
# @@@  NOTE: The following part is not necessary to change, just keep them.  @@@
#-------------------------------------------------------------------------------

set BIN = $masnum_home/source/bin
set EXP = `pwd`
if ($needmake == "YES")then
  cd $BIN
  make -f makefile clean
  make -f makefile
endif
cd $EXP

#-------------------------------------------------------------------------------

rm -rf masnum.wam.mpi wamyyz.nc ice_clim_mask.nc wave_rest.nc
cp $BIN/masnum.wam.mpi masnum.wam.mpi
cp -rf $masnum_home/inputdata/forcing ./
cp $depfile            wamyyz.nc
#ln -s $depfile            wamyyz.nc
#ln -s $depfile            ice_clim_mask.nc

#-------------------------------------------------------------------------------

cat > ctlparams << EOF
&CTLPARAMS
data_path   = ""          ,
wind_path   = "$wind_path",
TITLE       = "$title"    ,
CISTIME     = "$istime"    ,
CIETIME     = "$ietime"    ,
COOLS_DAYS  = $cools_days ,
DELTTM      = $delttm     ,
WNDFREQ     = $wndfreq    ,
WNDTYPE     = $wndtype    ,
OUTFLAG     = $outflag    ,
WIOFREQ     = $wiofreq    ,
CIOFREQ     = $ciofreq    ,
RSTFREQ     = $rstfreq
/
EOF
#CISTIME     = $istime     ,

#-------------------------------------------------------------------------------

#mpirun -np $nproc ./masnum.wam.mpi > out.qrunout

#-------------------------------------------------------------------------------

exit

#-------------------------------------------------------------------------------
# ***************************** THE END ****************************************
IETIME
