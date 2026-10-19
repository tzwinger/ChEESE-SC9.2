#!/bin/bash
# =============================================================================
# make_meshes.sh -- build + partition the footprint meshes for the sweep.
#
# Run from the ElmerIce/ run directory (where the SIFs live). It reads the Gmsh
# sources in Geometry/ and writes the Elmer mesh directories footprint_bc_<res>/
# (serial mesh + partitioning.<ranks>/) that the SIFs' `Mesh DB "." "footprint_bc_<res>"`
# expects.
#
# Per resolution it:
#   1. runs Gmsh to make the 2D .msh from the .geo, ONLY if the .msh is missing
#      (coarse ships only a .geo; fine/superfine/ultrafine already have a .msh);
#   2. converts the .msh (Gmsh fmt 14) to an Elmer mesh dir (fmt 2) with ElmerGrid;
#   3. partitions it to  nodes*128  ranks (the weak-scaling layout) with METIS.
# superfine_vert reuses the superfine FOOTPRINT (the finer vertical extrusion is
# applied at runtime), so it is just a copy of footprint_bc_superfine.
#
# Requires ElmerGrid (module load elmer/latest) and, for any missing .msh, gmsh.
# The 2D footprint is what gets partitioned; the runtime extrusion is per-partition.
#
#   bash make_meshes.sh                 # all resolutions
#   bash make_meshes.sh fine superfine  # only the named ones
# =============================================================================
set -eu

GEO=Geometry
declare -A RANKS=( [coarse]=128 [fine]=512 [superfine]=2048 [ultrafine]=4096 )
ALL=(coarse fine superfine ultrafine)       # superfine_vert handled separately (copy)

RES_LIST=("$@"); [[ ${#RES_LIST[@]} -gt 0 ]] || RES_LIST=("${ALL[@]}" superfine_vert)

build_one () {
  local res=$1 np=${RANKS[$1]:-}
  local msh="$GEO/footprint_bc_${res}.msh" geo="$GEO/footprint_bc_${res}.geo"
  local dir="footprint_bc_${res}"
  [[ -n "$np" ]] || { echo "ERROR: no rank count defined for '$res'"; exit 1; }

  # 1. Gmsh: .geo -> .msh only if the .msh is absent
  if [[ ! -f "$msh" ]]; then
    [[ -f "$geo" ]] || { echo "ERROR: neither $msh nor $geo exists for $res"; exit 1; }
    command -v gmsh >/dev/null || { echo "ERROR: $msh missing and gmsh not in PATH"; exit 1; }
    echo "== [$res] gmsh: $geo -> $msh"
    gmsh -2 "$geo" -o "$msh"
  fi

  # 2. ElmerGrid: Gmsh mesh (14) -> Elmer mesh directory (2)
  echo "== [$res] ElmerGrid convert: $msh -> $dir/"
  rm -rf "$dir"
  ElmerGrid 14 2 "$msh" -out "$dir"
  # some ElmerGrid builds write next to the input; move it into the run dir if so
  [[ -d "$dir" ]] || { [[ -d "$GEO/$dir" ]] && mv "$GEO/$dir" "$dir"; }
  [[ -f "$dir/mesh.header" ]] || { echo "ERROR: $dir/mesh.header not produced -- check ElmerGrid output path"; exit 1; }

  # 3. ElmerGrid: partition the Elmer mesh to np ranks (METIS k-way on the dual graph)
  echo "== [$res] ElmerGrid partition: $dir -> partitioning.$np"
  ElmerGrid 2 2 "$dir" -partdual -metiskway "$np"
  [[ -d "$dir/partitioning.$np" ]] || { echo "ERROR: $dir/partitioning.$np not produced"; exit 1; }
}

for res in "${RES_LIST[@]}"; do
  if [[ "$res" == "superfine_vert" ]]; then
    # reuses the superfine footprint (serial mesh + its 2048-partitioning)
    [[ -d footprint_bc_superfine ]] || { echo "ERROR: build 'superfine' before 'superfine_vert'"; exit 1; }
    echo "== [superfine_vert] copy of footprint_bc_superfine"
    rm -rf footprint_bc_superfine_vert
    cp -r footprint_bc_superfine footprint_bc_superfine_vert
  else
    build_one "$res"
  fi
done

echo
echo "done. meshes + partitionings in $(pwd):"
for res in coarse fine superfine superfine_vert ultrafine; do
  [[ -d "footprint_bc_${res}" ]] || continue
  printf "  %-16s %s\n" "$res" "$(ls -d footprint_bc_${res}/partitioning.* 2>/dev/null | tr '\n' ' ')"
done
