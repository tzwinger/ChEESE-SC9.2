#!/bin/bash -f
#SBATCH --ntasks-per-node=128
#### nodes are set PER RESOLUTION by submit_sweep.sh (weak scaling); this is only
#### the fallback when the file is sbatch'ed directly (superfine == 16 nodes).
#SBATCH --nodes=16
#SBATCH --exclusive
#SBATCH --mem=0
#### for how long (hh:mm:ss) ##########
#SBATCH -t 10:30:00
##### name of job, output file, error file ##########
#SBATCH -J bowdoin_sweep
#SBATCH -o Logs/%x_%A_%a.out
#SBATCH -e Logs/%x_%A_%a.err
##### project to bill from##########
#SBATCH --account=project_462000932
##### choose your queue #####
#SBATCH --partition=standard
##SBATCH --partition=debug
################## OpenMP Stuff ##########
#SBATCH --cpus-per-task=1

# ===========================================================================
# Parameter study: footprint resolution x iceberg attachment bandwidth
#                  ( x cloud-interpolation parameters, off by default ).
#
# TWO STAGES (set via STAGE env):
#   STAGE=init : run Bowdoin_footprint_init_incmp_lumi.sif ONCE for RES. The init
#                is independent of area/interp, so it is a single job per mesh;
#                it drops footprint_bc_<RES>/.init_<restart>.done and is skipped
#                if that marker already exists.
#   STAGE=run  : (default) one Slurm array index == one (area,idwpow,nn) combo
#                of RES; refuses to run unless the init marker is present.
#
# WEAK SCALING: each resolution runs on its own node count, so there is one job
# (init) + one array (run) PER RESOLUTION. Submit everything with:
#
#       make                 # build the shared objects ONCE, on a login node
#       bash submit_sweep.sh # init job (if needed) + dependent transient array, per res
#
# Drive one resolution by hand:
#       RES=fine bash runelmer_multiple.sh --list     # index -> combo map
#       RES=fine STAGE=init sbatch --nodes=4 --export=ALL,RES=fine,STAGE=init runelmer_multiple.sh
#       RES=fine sbatch --nodes=4 --array=0-2 --export=ALL,RES=fine runelmer_multiple.sh
# ===========================================================================

set -u

# ---- parameter grid -------------------------------------------------------
RES_ALL=(coarse fine superfine superfine_vert ultrafine)

# attachment-band widths [m] per resolution. Keep the smallest value >= ~1x
# the vertical element height at the front, else the band collapses onto a
# single node layer (mesh-limited). coarse ~30 m/elem, fine|superfine ~15 m.
# superfine_vert = the superfine footprint mesh with a FINER vertical extrusion
# (runresolution_superfine_vert.sif: neLo=21 neUp=7), A=30 only. It needs its own
# mesh dir footprint_bc_superfine_vert -- a copy of footprint_bc_superfine (the
# horizontal footprint is identical; the extrusion is applied at runtime):
#     cp -r footprint_bc_superfine footprint_bc_superfine_vert
declare -A AREA_LIST=(
  [coarse]="30 40 50"
  [fine]="10 20 30"
  [superfine]="10 20 30"
  [superfine_vert]="30"
  [ultrafine]="30"
)
# ultrafine = a finer HORIZONTAL footprint mesh (footprint_bc_ultrafine, its own 4096-part
# mesh / 32 nodes), with superfine_vert's vertical extrusion (runresolution_ultrafine.sif:
# neLo=21 neUp=7), A=30 -- for horizontal mesh convergence.

# cloud interpolation -- single value each == NOT swept (recommended default).
# Set e.g. IDW_LIST="1.0 2.0 3.0" / NN_LIST="4 6 12" to turn it into a dimension.
IDW_LIST="2.0"
NN_LIST="6"

# stock SIFs (unchanged templates; generated copies get a _<tag>/_<res> suffix)
RUN_SIF=Bowdoin_footprint_run_incmp_lumi.sif
INIT_SIF=Bowdoin_footprint_init_incmp_lumi.sif
RESTART=Bowdoin_v1_init_incmp          # == $restartfile in the run SIF (no .result)

marker_of () { echo "footprint_bc_${1}/.init_${RESTART}.done"; }

# ---- completion detection: a combo is "done" if a prior-run marker exists, OR a
#      complete output .dat (reaching FINAL_TS) is found. Lets submit_sweep.sh launch
#      ONLY the missing combos. Searches the live results dir and any archive roots.
FINAL_TS=${FINAL_TS:-470}                    # last timestep of a COMPLETE run
RESULT_ROOTS=${RESULT_ROOTS:-". New_VTUs"}   # dirs searched for existing outputs

combo_done () {   # $1=res $2=area $3=idwpow $4=nn  -> 0 (true) if already complete
  local res=$1 area=$2 tag="${1}_A${2}_p${3}_k${4}" root pref dat last
  [[ -e ".sweep_locks/run_${tag}.done" ]] && return 0
  for root in $RESULT_ROOTS; do
    for pref in VTU_footprint Output_footprint; do
      dat=$(ls "$root/${pref}_${res}_A${area}"/*_"${res}_A${area}"_*_front.dat 2>/dev/null | head -1)
      [[ -n "$dat" ]] || continue
      last=$(tail -n1 "$dat" | awk '{print $1}')
      [[ "${last:-0}" =~ ^[0-9]+$ && "$last" -ge "$FINAL_TS" ]] && return 0
    done
  done
  return 1
}

# ---- flatten the grid (only RES if set, else the whole grid for --list) ---
build_combos () {
  COMBOS=()
  local reslist
  if [[ -n "${RES:-}" ]]; then reslist=("$RES"); else reslist=("${RES_ALL[@]}"); fi
  local res area p k
  for res in "${reslist[@]}"; do
    for area in ${AREA_LIST[$res]}; do
      for p in $IDW_LIST; do
        for k in $NN_LIST; do
          COMBOS+=("$res $area $p $k")
        done
      done
    done
  done
}
build_combos

case "${1:-}" in
  --count)  echo "${#COMBOS[@]}"; exit 0 ;;
  --list)   for i in "${!COMBOS[@]}"; do printf '%4d  %s\n' "$i" "${COMBOS[$i]}"; done; exit 0 ;;
  --marker) : "${RES:?--marker needs RES}"; marker_of "$RES"; exit 0 ;;
  --missing) : "${RES:?--missing needs RES}"
             m=(); for i in "${!COMBOS[@]}"; do
               read -r _r _a _p _k <<<"${COMBOS[$i]}"
               combo_done "$_r" "$_a" "$_p" "$_k" || m+=("$i")
             done
             (IFS=,; echo "${m[*]}"); exit 0 ;;
esac

: "${RES:?set RES=<coarse|fine|superfine> (submit via submit_sweep.sh)}"
STAGE=${STAGE:-run}
marker=$(marker_of "$RES")
mkdir -p Logs .sweep_locks

# ---- environment (common to both stages) ----------------------------------
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
echo "running OpenMP on $OMP_NUM_THREADS"
export FI_CXI_RX_MATCH_MODE=hybrid
export MPICH_MAX_THREAD_SAFETY=multiple
module load elmer/latest

# rebuild SOs if needed (serialized so concurrent array tasks don't race object files)
( flock 9; make; ) 9>.sweep_locks/make.lock
echo "checking modules and executable:"; module list; which ElmerSolver_mpi

# ---- run_stage: a stage succeeds iff rank 0 printed "ALL DONE" -------------
run_stage () {
  local sif=$1
  local ntasks=${SLURM_NTASKS:-$(( SLURM_JOB_NUM_NODES * SLURM_NTASKS_PER_NODE ))}
  local log=Logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}_${sif%.sif}_np${ntasks}.log
  echo "$(date '+%F %T') starting ${sif}  (log: ${log})"
  srun ElmerSolver_mpi "${sif}" > "${log}" 2>&1
  local rc=$?
  if grep -q "ALL DONE" "${log}"; then
    (( rc != 0 )) && echo "WARNING: ${sif} reached ALL DONE but srun rc=${rc} (crash in finalize?) -> continuing"
    echo "$(date '+%F %T') finished ${sif}"; return 0
  fi
  echo "ERROR: ${sif} failed (srun rc=${rc}, no ALL DONE) -> aborting"; tail -n 30 "${log}"; return 1
}

# ===========================================================================
# STAGE: init -- once per mesh, independent of area/interp
# ===========================================================================
if [[ "$STAGE" == "init" ]]; then
  # The init SIF includes "params.sif" (like the run SIF). Generate a FLAT params file
  # for it -- inline runresolution_<RES>.sif ($res/$neLo/$neUp), never a nested include --
  # and repoint. The init only needs res/neLo/neUp; area/idwpow/nn are added so the params
  # set is complete/self-consistent (unused by the init SIF).
  [[ -f "runresolution_${RES}.sif" ]] || {
    echo "ERROR: runresolution_${RES}.sif not found in $(pwd) -> sync it to the run dir"; exit 1; }
  { cat "runresolution_${RES}.sif"
    echo "\$area=${AREA_LIST[$RES]%% *}"
    echo "\$idwpow=${IDW_LIST%% *}"
    echo "\$nn=${NN_LIST%% *}"
  } > "params_init_${RES}.sif"
  sed -E 's#include[[:space:]]+"?params\.sif"?#include "params_init_'"${RES}"'.sif"#' "$INIT_SIF" > "init_${RES}.sif"
  grep -q "params_init_${RES}.sif" "init_${RES}.sif" || {
    echo "ERROR: failed to repoint the params include in $INIT_SIF (check its 'include' line) -> aborting"; exit 1; }
  (
    flock 9
    if [[ -e "$marker" ]]; then
      echo "init for res=$RES already present ($marker) -> skipping"
    else
      echo "starting initialization run for res=$RES on ${SLURM_JOB_NUM_NODES} nodes"
      run_stage "init_${RES}.sif" && touch "$marker"
    fi
  ) 9>".sweep_locks/init_${RES}.lock"
  [[ -e "$marker" ]] || { echo "ERROR: init for res=$RES did not complete"; exit 1; }
  echo "init stage done for res=$RES"
  exit 0
fi

# ===========================================================================
# STAGE: run -- one array index == one transient combo (needs init marker)
# ===========================================================================
idx=${SLURM_ARRAY_TASK_ID:?submit with --array (see: RES=$RES bash $0 --list)}
read -r res area idwpow nn <<<"${COMBOS[$idx]}"
tag="${res}_A${area}_p${idwpow}_k${nn}"
echo "combo $idx -> res=$res area=$area idwpow=$idwpow nn=$nn  (tag=$tag)  nodes=${SLURM_JOB_NUM_NODES}"

if combo_done "$res" "$area" "$idwpow" "$nn"; then
  echo "combo $tag already complete -> skipping"; exit 0
fi

[[ -e "$marker" ]] || { echo "ERROR: init marker $marker missing -> run 'STAGE=init' for $res first"; exit 1; }

# Inline runresolution_<res>.sif (do NOT `include` it): Elmer's ReadAndTrim reads the
# first line of an included file without re-checking it for `include`, so a nested include
# (run -> params_<tag> -> runresolution) fails with "Unknown input field in header section".
# Keeping params_<tag>.sif flat (res/neLo/neUp/area/idwpow/nn) leaves a single include level.
[[ -f "runresolution_${res}.sif" ]] || {
  echo "ERROR: runresolution_${res}.sif not found in $(pwd) -> sync it to the run dir"; exit 1; }
{
  cat "runresolution_${res}.sif"
  echo "\$area=${area}"
  echo "\$idwpow=${idwpow}"
  echo "\$nn=${nn}"
} > "params_${tag}.sif"
sed -E 's#include[[:space:]]+"?params\.sif"?#include "params_'"${tag}"'.sif"#' "$RUN_SIF" > "run_${tag}.sif"
grep -q "params_${tag}.sif" "run_${tag}.sif" || {
  echo "ERROR: failed to repoint the params include in $RUN_SIF (check its 'include' line) -> aborting"; exit 1; }

echo "starting transient run for $tag"
run_stage "run_${tag}.sif" || exit 1
touch ".sweep_locks/run_${tag}.done"
echo "combo $tag completed"
