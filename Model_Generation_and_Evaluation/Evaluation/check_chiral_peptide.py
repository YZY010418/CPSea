#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import shutil
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Optional, Tuple


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vcross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vdot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def ca_chirality(N, CA, C, CB) -> str:
    normal = vcross(vsub(N, CA), vsub(C, CA))
    chi = vdot(normal, vsub(CB, CA))
    return "D" if chi < 0 else "L"


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def parse_atoms(path: Path) -> Dict[str, Dict[int, Dict[str, Optional[Tuple[float, float, float]]]]]:
    want_chains = {"L", "R"}
    data: Dict[str, Dict[int, Dict[str, Optional[Tuple[float, float, float]]]]] = {}
    with open_text(path) as f:
        for line in f:
            if not line.startswith("ATOM") or len(line) < 54:
                continue

            chain = (line[21] or "").strip()
            if chain not in want_chains:
                continue

            atom = line[12:16].strip()
            if atom not in ("N", "CA", "C", "CB"):
                continue

            try:
                res_id = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            chain_map = data.setdefault(chain, {})
            res_map = chain_map.setdefault(res_id, {"N": None, "CA": None, "C": None, "CB": None})
            res_map[atom] = (x, y, z)

    return data


def has_d_residue(path: Path) -> bool:
    try:
        data = parse_atoms(path)
    except Exception:
        return True

    for residues in data.values():
        if not residues:
            continue
        ids = sorted(residues.keys())
        for rid in ids:  # default: check all residues
            atoms = residues.get(rid, {})
            if None in (atoms.get("N"), atoms.get("CA"), atoms.get("C"), atoms.get("CB")):
                continue
            if ca_chirality(atoms["N"], atoms["CA"], atoms["C"], atoms["CB"]) == "D":
                return True
    return False


def safe_move(src: Path, dst_dir: Path) -> Path:
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    if not dst.exists():
        shutil.move(str(src), str(dst))
        return dst

    name = src.name
    i = 1
    while True:
        cand = dst_dir / f"{name}.{i}"
        if not cand.exists():
            shutil.move(str(src), str(cand))
            return cand
        i += 1


def worker(path_str: str) -> Tuple[str, bool]:
    p = Path(path_str)
    return path_str, has_d_residue(p)


def main():
    ap = argparse.ArgumentParser(description="Compute fully-correct ratio and move problematic PDB files.")
    ap.add_argument("-i", required=True, help="Input folder containing .pdb / .pdb.gz")
    ap.add_argument(
        "-d",
        default=None,
        help='Destination folder for problematic files (default: "<input>/problematic_pdb_files")',
    )
    ap.add_argument("--no_move", action="store_true", help="Only compute ratios; do not move files")
    args = ap.parse_args()

    folder = Path(args.i)
    if not folder.is_dir():
        print(f"ERROR: not a folder: {folder}", file=sys.stderr)
        sys.exit(1)

    files = sorted(
        p for p in folder.iterdir()
        if p.is_file() and (p.suffix == ".pdb" or p.name.endswith(".pdb.gz"))
    )
    if not files:
        print("ERROR: no .pdb or .pdb.gz files found.", file=sys.stderr)
        sys.exit(1)

    dest = Path(args.d) if args.d else (folder / "problematic_pdb_files")

    t0 = time.time()
    total = len(files)
    problematic = 0
    moved_paths = []

    try:
        from tqdm import tqdm  # type: ignore
        it = tqdm(total=total, unit="file")
        use_tqdm = True
    except Exception:
        it = None
        use_tqdm = False

    with Pool() as pool:
        for path_str, is_bad in pool.imap_unordered(worker, (str(p) for p in files), chunksize=200):
            if is_bad:
                problematic += 1
                if not args.no_move:
                    moved_paths.append(str(safe_move(Path(path_str), dest)))
            if use_tqdm:
                it.update(1)  # type: ignore

    if use_tqdm:
        it.close()  # type: ignore

    correct = total - problematic
    print(f"total_files={total}", file=sys.stderr)
    print(f"correct_files={correct}", file=sys.stderr)
    print(f"correct_ratio={correct/total*100:.2f}%", file=sys.stderr)
    print(f"problematic_files={problematic}", file=sys.stderr)
    print(f"problematic_ratio={problematic/total*100:.2f}%", file=sys.stderr)
    if not args.no_move:
        print(f"moved_to={dest}", file=sys.stderr)

    if moved_paths:
        for p in moved_paths:
            print(p)

    print(f"elapsed_sec={time.time()-t0:.1f}", file=sys.stderr)


if __name__ == "__main__":
    main()
