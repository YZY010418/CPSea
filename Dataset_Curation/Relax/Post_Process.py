import os, gzip, shutil, csv, time, argparse, gc, sys
from multiprocessing import Pool, Manager
from functools import partial
from threading import Thread

def load_energy_data(tsv_path):
    """Load filename->energy map from a TSV without header (filename in col2, energy in col6)."""
    energy_map, t0 = {}, time.time()
    with open(tsv_path, 'r', newline='') as f:
        for i, row in enumerate(csv.reader(f, delimiter='\t'), 1):
            if len(row) < 6: 
                continue
            fn = row[1].strip()
            try:
                energy_map[fn] = float(row[5].strip())
            except ValueError:
                continue
    print(f"Loaded {len(energy_map)} energy records in {time.time()-t0:.2f}s")
    return energy_map

def create_output_directories(base_dir):
    d = {'energy_bad': os.path.join(base_dir, 'energy_bad'),
         'chiral_bad': os.path.join(base_dir, 'chiral_bad')}
    for p in d.values(): os.makedirs(p, exist_ok=True)
    return d

def filter_by_energy(pdb_files_with_dir, energy_map, output_dirs, thr):
    """Move high-energy files to energy_bad; return passed list [(filename, dir)]."""
    t0 = time.time()
    moved, passed = 0, []
    n = len(pdb_files_with_dir)
    print(f"Found {n} PDB files")
    for k, (fn, d) in enumerate(pdb_files_with_dir, 1):
        if fn not in energy_map: 
            continue
        if energy_map[fn] >= thr:
            src = os.path.join(d, fn)
            dst = os.path.join(output_dirs['energy_bad'], fn)
            if os.path.exists(dst):
                name, ext = os.path.splitext(fn)
                dst = os.path.join(output_dirs['energy_bad'], f"{name}_dup{ext}")
            shutil.move(src, dst); moved += 1
        else:
            passed.append((fn, d))
        if k % 10000 == 0:
            print(f"Energy filter: {k}/{n} ({k/n*100:.2f}%), moved {moved}")
    print(f"Energy filter done in {time.time()-t0:.2f}s | moved: {moved}, kept: {len(passed)}")
    return passed

def subtract(v1, v2): return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])
def cross(v1, v2): return (v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])
def dot(v1, v2): return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def ca_chirality(N, CA, C, CB):
    """Return 'D' if negative triple product else 'L'."""
    n = cross(subtract(N, CA), subtract(C, CA))
    return 'D' if dot(n, subtract(CB, CA)) < 0 else 'L'

def parse_pdb_gz(path, include_r=True):
    """Return {res_id: {'N':(x,y,z),'CA':...,'C':...,'CB':...}} for chain L and (optionally) R."""
    data = {}
    with gzip.open(path, 'rt', errors='ignore') as f:
        for line in f:
            if not line.startswith('ATOM'): 
                continue
            chain = line[21] if len(line) > 21 else ''
            if include_r:
                if chain not in ('L','R'): 
                    continue
            else:
                if chain != 'L': 
                    continue
            atom = line[12:16].strip() if len(line) > 16 else ''
            if atom not in ('N','CA','C','CB'): 
                continue
            try:
                rid = int(line[22:26].strip())
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            if rid not in data: data[rid] = {'N':None,'CA':None,'C':None,'CB':None}
            data[rid][atom] = (x,y,z)
    return data

def process_single(file_info, q, check_all=True, include_r=True):
    fn, d = file_info
    has_problem = False
    data = parse_pdb_gz(os.path.join(d, fn), include_r)
    if data:
        ids = sorted(data.keys())
        if check_all:
            for rid in ids:
                a = data[rid]
                if None not in (a['N'], a['CA'], a['C'], a['CB']) and ca_chirality(a['N'], a['CA'], a['C'], a['CB']) == 'D':
                    has_problem = True; break
        else:
            a = data[ids[0]]
            if None not in (a['N'], a['CA'], a['C'], a['CB']) and ca_chirality(a['N'], a['CA'], a['C'], a['CB']) == 'D':
                has_problem = True
            if not has_problem and ids[0] != ids[-1]:
                b = data[ids[-1]]
                if None not in (b['N'], b['CA'], b['C'], b['CB']) and ca_chirality(b['N'], b['CA'], b['C'], b['CB']) == 'D':
                    has_problem = True
    q.put(1); gc.collect()
    return (fn, d, has_problem)

def progress_worker(q, total, every=1000):
    done = 0
    while done < total:
        done += q.get()
        if done % every == 0 or done == total:
            print(f"Chirality check: {done}/{total} ({done/total*100:.2f}%)", file=sys.stderr)

def filter_by_chirality(passed_energy, output_dirs, nproc):
    t0 = time.time()
    total = len(passed_energy)
    print(f"Chirality check on {total} files using {nproc or 'all'} processes (chains L+R, all residues)")
    m = Manager(); q = m.Queue()
    t = Thread(target=progress_worker, args=(q, total), daemon=True); t.start()
    func = partial(process_single, q=q, check_all=True, include_r=True)
    with Pool(processes=nproc) as pool:
        results = pool.map(func, passed_energy)
    t.join()
    bad, ok = [], []
    for fn, d, badflag in results:
        if badflag: bad.append((fn, d))
        else: ok.append(fn)
    for fn, d in bad:
        src = os.path.join(d, fn)
        dst = os.path.join(output_dirs['chiral_bad'], fn)
        if os.path.exists(dst):
            name, ext = os.path.splitext(fn)
            dst = os.path.join(output_dirs['chiral_bad'], f"{name}_dup{ext}")
        shutil.move(src, dst)
    print(f"Chirality check done in {time.time()-t0:.2f}s | moved: {len(bad)}, passed: {len(ok)}")
    return ok

def update_tsv_file(src_tsv, dst_tsv, passed_files):
    t0 = time.time()
    keep = set(passed_files)
    total = kept = 0
    with open(src_tsv, 'r', newline='') as inf, open(dst_tsv, 'w', newline='') as outf:
        r = csv.reader(inf, delimiter='\t'); w = csv.writer(outf, delimiter='\t')
        for row in r:
            total += 1
            if len(row) >= 2 and row[1].strip() in keep:
                w.writerow(row); kept += 1
            if total % 100000 == 0:
                print(f"TSV update: {total} rows, kept {kept}")
    print(f"TSV updated in {time.time()-t0:.2f}s | total: {total}, kept: {kept}")

def collect_pdb_files_from_dirs(dirs):
    out = []
    for d in dirs:
        if not os.path.isdir(d): 
            continue
        for f in os.listdir(d):
            if f.endswith('.pdb.gz'):
                out.append((f, d))
    return out

def main():
    p = argparse.ArgumentParser(description='PDB energy & chirality filter (multi-dir)')
    p.add_argument('pdb_dir1'); p.add_argument('pdb_dir2'); p.add_argument('pdb_dir3')
    p.add_argument('tsv_file')
    p.add_argument('-p','--processes', type=int)
    p.add_argument('-t','--threshold', type=float, default=0.0)
    args = p.parse_args()

    pdb_dirs = [args.pdb_dir1, args.pdb_dir2, args.pdb_dir3]
    for d in pdb_dirs:
        if not os.path.isdir(d):
            print(f"Error: {d} is not a directory"); return
    if not os.path.isfile(args.tsv_file):
        print(f"Error: {args.tsv_file} is not a file"); return

    T0 = time.time()
    out_dirs = create_output_directories(os.path.dirname(args.tsv_file))
    files = collect_pdb_files_from_dirs(pdb_dirs)
    if not files:
        print("Error: no .pdb.gz files found"); return

    energy_map = load_energy_data(args.tsv_file)
    kept_energy = filter_by_energy(files, energy_map, out_dirs, args.threshold)
    if not kept_energy:
        print("No files passed energy filter"); return

    kept_chiral = filter_by_chirality(kept_energy, out_dirs, args.processes)

    new_tsv = os.path.join(os.path.dirname(args.tsv_file), 'energy_data_new.tsv')
    update_tsv_file(args.tsv_file, new_tsv, kept_chiral)

    total_time = int(time.time() - T0)
    h, rem = divmod(total_time, 3600)
    m, s = divmod(rem, 60)
    total_original = len(files)
    total_energy_bad = len(os.listdir(out_dirs['energy_bad']))
    total_chiral_bad = len(os.listdir(out_dirs['chiral_bad']))

    print("\n" + "="*60)
    print("Done")
    print(f"Elapsed: {h}h {m}m {s}s")
    print(f"Original files: {total_original}")
    print(f"High-energy moved: {total_energy_bad}")
    print(f"Chirality moved: {total_chiral_bad}")
    print(f"Passed: {len(kept_chiral)}")
    print(f"New TSV: {new_tsv}")
    print("="*60)

if __name__ == "__main__":
    main()
