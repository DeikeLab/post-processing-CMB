#!/bin/bash
#SBATCH --job-name=stats          # create a short name for your job
#SBATCH --nodes=1                   # node count
#SBATCH --ntasks-per-node=112       # number of tasks per node
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G            # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00             # total run time limit (HH:MM:SS)
##SBATCH --time=00:30:00             # total run time limit (HH:MM:SS)
#SBATCH --account=deike             # total run time limit (HH:MM:SS)
####SBATCH --mail-type=begin        # send email when job begins
####SBATCH --mail-type=end          # send email when job ends
####SBATCH --mail-user=ns8802@princeton.edu
#
module purge
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load intel-tbb/2021.13 intel-rt/2024.2
module load anaconda3/2024.10

ulimit -s unlimited
#
# Input parameters
#
Re_ast=720.0;
BO=200.0;
Re_wave=2.5597571943e+04;
UstarRATIO=0.5;
ak=0.3;
r_L0lam=4.0;
rho_r=0.001225;           
mu_r=2.2471881948940954e-06; 
MAXLEVEL=10;
MINLEVEL=6;
RELEASETIME=1629.28;
uemaxRATIO=0.3;
end_sim=100000000000000.0;
do_eta_loc=1;
do_profile=1;
do_fields=0;
do_slices=0;
do_tagging=0;
prt_res=9; 
from_pr=0;
st_wave=0;
rand_num=2;
N_mod=64;
do_eta_stats=1;
eta_stats_stride=1;
eta_stats_window_Tp=2.0;
eta_pdf_bins=161;
eta_stats_keep_last_windows=50;
tout_eta_plot_my=8;   # partial checkpoint for live plots: every T0_/8 ~ 56 min wall time

# Live plot refresh cadence (seconds) while the solver is running
plot_refresh_sec=300;
#
srun windwave_turb -m 23:59:00 $Re_ast $BO $Re_wave $UstarRATIO $ak $r_L0lam \
	                       $rho_r $mu_r $MAXLEVEL $MINLEVEL $RELEASETIME \
	                       $uemaxRATIO $end_sim \
			       $do_eta_loc $do_profile $do_fields $do_slices $do_tagging \
	                       $prt_res $from_pr \
	                       $st_wave $rand_num $N_mod \
	                       $do_eta_stats $eta_stats_stride $eta_stats_window_Tp \
	                       $eta_pdf_bins $eta_stats_keep_last_windows \
	                       $tout_eta_plot_my > out.log 2>&1 &
solver_pid=$!

# Keep refreshing plots from accumulated window files while solver is alive.
(
	while kill -0 "$solver_pid" 2>/dev/null; do
		if ! python3 generate_eta_plots.py statistics >> out.log 2>&1; then
			echo "[plot] WARN: live plotting failed; will retry." | tee -a out.log >&2
		fi
		sleep "$plot_refresh_sec"
	done
) &
plotter_pid=$!

wait "$solver_pid"
solver_rc=$?

# Stop live plot loop after solver exit.
kill "$plotter_pid" 2>/dev/null || true
wait "$plotter_pid" 2>/dev/null || true

# Final pass to ensure newest windows are plotted.
if ! python3 generate_eta_plots.py statistics >> out.log 2>&1; then
	echo "[plot] ERROR: final generate_eta_plots.py failed. Check Python dependencies (e.g., pyparsing)." | tee -a out.log >&2
fi

exit "$solver_rc"
#
