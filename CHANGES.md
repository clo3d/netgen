# Changes from Netgen 6.2.2604

Applied in this order, one commit each (see `git log v6.2.2604..main`). `git diff -w v6.2.2604 main`
shows the changes without the line-ending conversion in item 0.

0. **Line endings (no content change):** `libsrc/occ/occgeom.cpp`, `libsrc/occ/occgeom.hpp` and
   `libsrc/occ/occgenmesh.cpp` converted from CRLF to LF.

## Meshing robustness on degenerate geometry
1. `libsrc/meshing/meshing2.cpp`: the advancing front gives up on a face only after 2000 consecutive
   trials without creating a surface element. Fixes an infinite loop on degenerate sliver faces.
2. `libsrc/occ/occgeom.cpp`, `libsrc/occ/occgeom.hpp`:
   - A face that fails to mesh is retried: parameter space, then plane space, then both again at
     `maxh * 0.3`, then meshing from the face boundary. Keeping a partial mesh as the last resort is
     now behind `OCCParameters::keep_partial_mesh` (default off).
   - `BuildFMap`: an edge that does not expose two vertices is skipped instead of read out of bounds.
3. `libsrc/occ/occgenmesh.cpp`:
   - New `OCCMeshFaceFromBoundary`: ear-clips a face's own boundary segments for otherwise
     unmeshable single-loop faces. `OCCMeshFace` takes a `maxh_scale` retry parameter.
   - `OCCSetLocalMeshSize`: skips edges with no 3D curve, skips a face that cannot be triangulated
     (instead of throwing away the whole mesh), and skips faces with zero or negative area (which
     made the mesher run forever).
4. `libsrc/occ/occ_edge.cpp`: edges without a 3D curve (seam/pole edges) no longer crash the
   `OCCEdge` constructor or `GetPoint`; `CalcStep` is implemented (it threw "not implemented").

## Correctness
5. `libsrc/meshing/improve2.cpp`: `CombineImprove` allocates `candidate_edges` for both
   orientations of every edge (`2 * edges.Size()`), fixing a buffer overflow.

## Build (Windows, MSVC Debug; Release output unchanged)
6. `CMakeLists.txt`: `CMAKE_MSVC_RUNTIME_LIBRARY` uses the debug runtime in Debug builds.
7. `libsrc/core/bitarray.hpp`: `NGCORE_API` removed from the inline template method `TBitArray::Or`.
8. `cmake/SuperBuild.cmake`: the static zlib file name follows the build configuration with
   `BUILD_ZLIB=ON`.

## Performance
9. Parallel per-face surface meshing: `libsrc/meshing/basegeom.cpp`, `libsrc/meshing/basegeom.hpp`,
   `libsrc/meshing/meshing2.cpp`, `libsrc/occ/occgeom.cpp`, `libsrc/occ/occgeom.hpp`,
   `libsrc/gprim/adtree.cpp`, `libsrc/gprim/geomtest3d.cpp`, `libsrc/meshing/smoothing2.cpp`.
   - Surface meshing honors `MeshingParameters::parallel_meshing` / `nthreads`. Faces are meshed on
     worker threads into per-thread mesh copies and merged in face order, so the result is the same
     for any thread count of 2 or more.
   - `NG_PARALLEL_MESHING=0` restores the previous sequential behavior; `NG_NTHREADS=N` overrides
     the thread count.
   - Thread-safety fixes in `adtree` and `geomtest3d` (thread-local scratch) and in OCC property
     caches, and a guard around worker threads on Windows.
   - `smoothing2.cpp`: no longer reads the empty `lochs` array when `uselocalh` is off.
