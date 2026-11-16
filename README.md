# netgen (CLO modified)

CLO Virtual Fashion's modified Netgen 6.2.2604, as built into CLO products. Published to make the
modified source available under the GNU Lesser General Public License v2.1 (`LICENSE`, unchanged
from the original).

Netgen's own README is `README.netgen.md` (renamed from `README.md` so this file is the repository's
front page). All original copyright notices and license headers are kept.

## Specs
- Original: Netgen 6.2.2604, https://github.com/NGSolve/netgen/tree/v6.2.2604
  (tag `v6.2.2604`, commit `3ee489c7d58fdbc2a6708cca3cbaefaae506dc17`)
- License: LGPL 2.1 (`LICENSE`)
- This repository holds the full upstream history up to that tag, with CLO's changes as separate
  commits on top. To see every change: `git diff v6.2.2604 main`.
- `CHANGES.md` lists the modified files and what each change does.

## Build notes
- Built as shared libraries (`nglib`, `ngcore`) and loaded dynamically.
- Windows: built with `-DBUILD_ZLIB=ON`, which links zlib 1.2.11 (zlib License) statically into
  `nglib`. macOS: zlib is the system library.
