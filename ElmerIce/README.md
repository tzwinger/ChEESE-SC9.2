# Instructions how to run the case

## Directories

**Geometry** contains the mesh-files (in Gmsh geo and msh format). To obtain the correct mesh for Elmer, one has to (here for the coarsest footprint mesh, which runs on a normal size workstation):

    ElmerGrid 14 2 footprint_bc_coarse.msh -autoclean

and then copy/link the new created folder `footprint_bc_coarse` to to the directory above (e.g. `ln -s footprint_bc_coarse ../footprint_bc_coarse`)

The partitioning of the mesh (make one partition for each MPI-rank) is done by (again with the coarsest mesh example, here for 12 ranks):

    ElmerGrid 2 2 footprint_bc_coarse -partdual -metiskway 12

To do all of this **systematically for every resolution** in one go — the Gmsh step (only where just a `.geo` exists, i.e. coarse), the ElmerGrid conversion, and the partitioning to each resolution's `nodes×128` rank count, plus the `superfine_vert` copy — run the helper script from this directory:

    bash make_meshes.sh                 # all resolutions
    bash make_meshes.sh fine superfine  # or only the named ones

It writes the `footprint_bc_<res>/` directories straight into this run directory (no separate copy/link step needed) and partitions to: coarse 128, fine 512, superfine / superfine_vert 2048, ultrafine 4096. It needs `gmsh` on `PATH` for the coarse `.geo` (the finer meshes already ship a `.msh`); if your ElmerGrid build uses `-metis N` rather than `-metiskway N`, adjust that one line in the script.

**SRC** contains the source code of the additional solvers and userfunctions. Compilation of these (the shared objects will thereafter reside in **BIN**) is done with the provided Makefile. To makes sure that a compatible set of shared objects is created, do a 

    make clean && make all
    
## Running jobs and arrays on LUMI-C

The parameter sweep (footprint resolution × iceberg attachment-band width) is driven by three
scripts:

| script | role |
|---|---|
| `submit_sweep.sh` | orchestrator — submits init + dependent transient array per resolution |
| `runelmer_multiple.sh` | the Slurm worker (two stages: `STAGE=init`, then `STAGE=run`) |
| `clean_sweep.sh` | removes only the files the scripts generate |

It is **weak-scaled**: each resolution runs on its own node count.

| resolution | nodes | ranks (`nodes×128`) |
|---|---|---|
| coarse | 1 | 128 |
| fine | 4 | 512 |
| superfine / superfine_vert | 16 | 2048 |
| ultrafine | 32 | 4096 |

### Prerequisites

```bash
git pull                     # sync the SIFs/scripts
make                         # build BIN/*.so on a login node  (make clean && make after a .F90 change)
```

- each mesh `footprint_bc_<res>` must be pre-partitioned to `nodes×128` ranks;
- the billing project goes in `--account` (default `project_462000932`, override with `ACCOUNT=...`);
- walltime / partition come from the `#SBATCH` header of `runelmer_multiple.sh`.

### Submit the whole sweep

```bash
bash submit_sweep.sh
```

For each resolution this submits one **init** job (if its marker is missing) and a **transient
array** that depends on it (`afterok`), launching only the combos that are still missing. Job ids
are logged to `sweep_jobs.txt`. It is **resumable** — rerun it and it skips completed combos.

```bash
ACCOUNT=project_XXXXXXXXX bash submit_sweep.sh   # bill a different project
CAP=4 bash submit_sweep.sh                        # cap simultaneously-running array tasks
```

### Inspect the grid / what is still missing

```bash
RES=fine bash runelmer_multiple.sh --list        # index -> (area,idwpow,nn) combo map
RES=fine bash runelmer_multiple.sh --missing     # comma-list of indices still to run
```
### Test / drive a single job by hand

Bypass the orchestrator and `sbatch` the worker directly, setting the node count for that
resolution:

```bash
# init once per mesh (ocean pressure; writes the restart the transient needs)
RES=fine STAGE=init sbatch --nodes=4 \
         --export=ALL,RES=fine,STAGE=init runelmer_multiple.sh

# one transient combo -- pick the array index from --list (here A30)
idx=$(RES=fine bash runelmer_multiple.sh --list | awk '$3=="30"{print $1; exit}')
RES=fine sbatch --nodes=4 --array=$idx \
         --export=ALL,RES=fine runelmer_multiple.sh
```

Force an already-complete combo to re-run (e.g. after a solver change) — remove its marker and
appended `.dat`, and point the completeness check at the live dir only:

```bash
rm -f .sweep_locks/run_fine_A30_*.done VTU_footprint_fine_A30/*_front.dat*
RES=fine sbatch --nodes=4 --array=$idx \
         --export=ALL,RES=fine,RESULT_ROOTS=. runelmer_multiple.sh
```
### Monitor, cancel, statistics

```bash
squeue --me                                  # queue state
scancel <jobid>                              # one job or array

# cancel only this sweep's jobs (init + transient), leaving anything else alone:
squeue --me -h -o "%i %j" | awk '$2 ~ /^bowdoin_/{print $1}' | xargs -r scancel

bash sweep_stats.sh                          # per-task runtime + allocated core-hours (from sweep_jobs.txt)
bash sweep_stats.sh <jobid> <jobid> ...      # or query given job ids
```

### Clean up

`clean_sweep.sh` removes **only** files the scripts generate — never sources, meshes, loads or
`.so`. Always dry-run with `-n` first.

```bash
bash clean_sweep.sh -n            # dry-run: list what WOULD be removed
bash clean_sweep.sh               # generated SIFs + transient .dat(.names) + .sweep_locks
                                  #   (keeps .bp / VTU_footprint_* dirs and the init restarts)
bash clean_sweep.sh --outputs     # + VTU_footprint_* result dirs and .bp files
bash clean_sweep.sh --init        # + init restarts + markers  (forces a re-init next submit)
bash clean_sweep.sh --logs        # + Logs/
bash clean_sweep.sh --manifest    # + sweep_jobs.txt
bash clean_sweep.sh --all         # everything the scripts produced
RES=fine bash clean_sweep.sh ...  # restrict any of the above to one resolution
```
## Geometry data
The file `bowdoin_v1_geometry_filled.nc` is derived from the DEM used in (van Dongen et al., 2019, doi:10.1017/jog.2019.89). Any further usage of the data should mention this publication. 


