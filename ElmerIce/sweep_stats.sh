#!/bin/bash
# Runtime + core-hour statistics for the sweep jobs.
#
#   bash sweep_stats.sh                 # read job ids from $MANIFEST (sweep_jobs.txt)
#   bash sweep_stats.sh 123456 123460   # or query given (array-parent / init) job ids
#
# Per task: wall time, allocated cores, allocated core-hours (= AllocCPUS*Elapsed,
# the LUMI billing proxy), resolved back to its (res,area,idwpow,nn) combo.
# sacct reads the Slurm accounting DB, so it works after logout and for finished
# *and* still-running jobs (running tasks show time accrued so far).
set -u

MANIFEST=${MANIFEST:-sweep_jobs.txt}
WORKER=${WORKER:-runelmer_multiple.sh}

# ---- gather job ids (+ res label per id) ----------------------------------
declare -A RESOF
ids=()
if (($#)); then
  ids=("$@")
else
  [[ -f "$MANIFEST" ]] || { echo "no manifest '$MANIFEST' and no ids given"; exit 1; }
  # columns: submit_time res nodes combos init_jobid array_jobid
  while read -r _t res _nodes _nc initjid arrayjid _rest; do
    [[ "${res:-}" == \#* || -z "${res:-}" ]] && continue
    [[ "$arrayjid" =~ ^[0-9]+$ ]] && { ids+=("$arrayjid"); RESOF[$arrayjid]="$res"; }
    [[ "$initjid"  =~ ^[0-9]+$ ]] && { ids+=("$initjid");  RESOF[$initjid]="init:$res"; }
  done < "$MANIFEST"
fi
((${#ids[@]})) || { echo "no job ids found"; exit 1; }
idcsv=$(IFS=,; echo "${ids[*]}")

# ---- combo lookup (res + array index -> "area idwpow nn"), cached ----------
declare -A LIST
combo_of () {  # $1=res  $2=array-index
  [[ -n "${LIST[$1]:-}" ]] || LIST[$1]=$(RES="$1" bash "$WORKER" --list 2>/dev/null)
  awk -v i="$2" '$1==i{$1="";sub(/^[ \t]+/,"");print;exit}' <<<"${LIST[$1]}"
}

# ---- query + report -------------------------------------------------------
fmt=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,NNodes,ExitCode
printf '%-14s %-11s %9s %6s %10s  %s\n' JobID State wall[s] cores core-h combo
tot_ch=0; tot_task=0
declare -A CH_BY_RES

while IFS='|' read -r jid state elraw alloc cpuraw nnodes ecode; do
  [[ "$jid" == *.* ]] && continue                     # skip steps (shouldn't occur with -X)
  [[ "$cpuraw" =~ ^[0-9]+$ ]] || cpuraw=0
  ch=$(awk -v c="$cpuraw" 'BEGIN{printf "%.2f", c/3600}')

  label="-"; reskey="-"
  if [[ "$jid" == *_* ]]; then                         # array task  base_idx
    base=${jid%_*}; idx=${jid##*_}
    reskey=${RESOF[$base]:-}
    [[ "$idx" =~ ^[0-9]+$ && -n "$reskey" ]] && label="$(combo_of "$reskey" "$idx")"  # "res area idwpow nn"
  else                                                 # single (init) job
    reskey=${RESOF[$jid]:-}; label="$reskey"
  fi

  printf '%-14s %-11s %9s %6s %10s  %s\n' "$jid" "$state" "$elraw" "$alloc" "$ch" "${label:--}"
  tot_ch=$(awk -v a="$tot_ch" -v b="$ch" 'BEGIN{printf "%.2f", a+b}')
  CH_BY_RES[$reskey]=$(awk -v a="${CH_BY_RES[$reskey]:-0}" -v b="$ch" 'BEGIN{printf "%.2f", a+b}')
  ((tot_task++))
done < <(sacct -X -n -P -j "$idcsv" --format="$fmt")

echo "--------------------------------------------------------------------"
for k in "${!CH_BY_RES[@]}"; do printf '  %-14s %10s core-h\n' "$k" "${CH_BY_RES[$k]}"; done
printf 'tasks: %d    TOTAL allocated core-hours: %s\n' "$tot_task" "$tot_ch"
echo "(core-h = AllocCPUS*Elapsed; running tasks counted as accrued so far)"
