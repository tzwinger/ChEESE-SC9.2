#!/bin/bash
# Submit the parameter study, weak-scaled per footprint resolution:
#   coarse -> 1 node (128 ranks), fine -> 4 (512), superfine|superfine_vert -> 16 (2048)
#
# RESUMABLE: only combos whose results are missing are launched. The worker's
# --missing reports them (a combo counts as done if a .sweep_locks/run_<tag>.done
# marker exists OR a complete output .dat reaching FINAL_TS is found under
# RESULT_ROOTS, default ". New_VTUs"). Override: RESULT_ROOTS="." FINAL_TS=470 ...
# NOTE completeness != correctness: an old run that reached FINAL_TS with garbage
# would be treated as done, so clean stale outputs (clean_sweep.sh) first if needed.
#
# superfine_vert = superfine footprint with finer vertical extrusion (21/7), A=30;
# provide its mesh once:  cp -r footprint_bc_superfine footprint_bc_superfine_vert
#
# For each resolution:
#   * if the init marker is missing, submit a single init job (STAGE=init) and
#     make the transient array depend on it (afterok) -- the init runs exactly
#     ONCE per mesh, independent of area/interp;
#   * otherwise submit the transient array straight away (init already done).
#
# Run from the (login-node) working directory after `make`.
# The per-res mesh footprint_bc_<res> must be pre-partitioned to nodes*128 ranks.
# Optional: CAP=<n> caps simultaneously-running array tasks.
set -eu

declare -A NODES=( [coarse]=1 [fine]=4 [superfine]=16 [superfine_vert]=16 )
SCRIPT=runelmer_multiple.sh
ACCOUNT=${ACCOUNT:-project_462000932}   # billed project (overridable via env)
CAP=${CAP:-}
MANIFEST=${MANIFEST:-sweep_jobs.txt}    # job-id log for post-run sacct (sweep_stats.sh)
# partition / walltime still come from the #SBATCH header of $SCRIPT

[[ -e "$MANIFEST" ]] || echo "# submit_time res nodes combos init_jobid array_jobid" > "$MANIFEST"

for res in coarse fine superfine superfine_vert; do
  # resumable: ask the worker which combo indices still need running (missing output
  # or no completion marker). Override the search with RESULT_ROOTS=... / FINAL_TS=...
  miss=$(RES=$res bash "$SCRIPT" --missing)
  [[ -n "$miss" ]] || { echo "$res: all combos already present -> skip"; continue; }
  n=$(awk -F, '{print NF}' <<<"$miss")
  marker=$(RES=$res bash "$SCRIPT" --marker)
  arr="$miss"; [[ -n "$CAP" ]] && arr="${arr}%${CAP}"

  dep=""; jid=none
  if [[ -e "$marker" ]]; then
    echo "$res: init present ($marker) -> transient array only"
  else
    jid=$(sbatch --parsable \
                 --account=${ACCOUNT} \
                 --nodes=${NODES[$res]} --ntasks-per-node=128 \
                 --job-name=bowdoin_init_${res} \
                 --export=ALL,RES=${res},STAGE=init \
                 "$SCRIPT")
    echo "$res: init job ${jid} (${NODES[$res]} nodes)"
    dep="--dependency=afterok:${jid}"
  fi

  aid=$(sbatch --parsable $dep \
               --account=${ACCOUNT} \
               --nodes=${NODES[$res]} --ntasks-per-node=128 \
               --job-name=bowdoin_sweep_${res} \
               --array="${arr}" \
               --export=ALL,RES=${res} \
               "$SCRIPT")
  echo "$res: transient array ${aid} (${NODES[$res]} nodes, ${n} missing combo(s), array ${arr})${dep:+ depends on init}"
  printf '%(%F)T_%(%T)T %s %s %s %s %s\n' -1 -1 "$res" "${NODES[$res]}" "$n" "$jid" "$aid" >> "$MANIFEST"
done
