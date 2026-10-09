#!/usr/bin/env python3
"""Regenerate preparations and run the four-condition factorial for the VS benchmark.

Conditions (all with the substructure filter OFF, so all four act on the full library):
  A  naive preparation (OpenBabel gen3d)      x  naive parameters (20 A box, exhaustiveness 8)
  C  skill preparation (RDKit ETKDGv3+Meeko)  x  skill parameters (25 A box, exhaustiveness 32)
  D  skill preparation                        x  naive parameters
  E  naive preparation                        x  skill parameters

The archived deposit has no PDBQT files, so ligands are regenerated from the deposited SMILES
and receptors from the PDB entries used in the paper. Every capability that cannot be
reproduced exactly is reported in capabilities.json rather than silently substituted.

Subcommands
  caps                     report tool versions and available preparation entry points
  prep-receptor TARGET     fetch the PDB entry and prepare it by both routes
  prep-ligands  TARGET     regenerate ligand PDBQT by both routes from the deposited SMILES
  dock TARGET COND         dock one condition
  all TARGET               prep-receptor, prep-ligands, then dock A, D, C, E
"""
import argparse, concurrent.futures as cf, csv, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path

PDB = {"fpr2":("7T6S","R"),"egfr":("1M17","A"),"cdk2":("1H1Q","A"),"cox2":("3LN1","A"),
       "esr1":("1SJ0","A"),"dpp4":("2RGU","A"),"ache":("4EY7","A"),"bace1":("4IVS","A"),
       "hsp90":("2WI7","A"),"p38":("1KV2","A"),"thrombin":("2ZDM","A")}
PARAMS = {"naive":{"box":20.0,"exh":8},"skill":{"box":25.0,"exh":32}}
COND = {"A":("naive","naive"),"C":("skill","skill"),"D":("skill","naive"),"E":("naive","skill")}
NUM_MODES = 9
ROOT = Path(os.environ.get("VSBENCH_ROOT", Path.cwd()))
DATA = ROOT / "data"           # unpacked deposit: data/targets/<t>/...
WORK = ROOT / "work"           # regenerated inputs
OUT  = ROOT / "out"            # score tables and manifests
for d in (WORK, OUT): d.mkdir(parents=True, exist_ok=True)


def run(cmd, timeout=None, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True,
                          timeout=timeout, **kw)


def caps():
    c = {"python": sys.version.split()[0], "executable": sys.executable,
         "host": os.uname().nodename, "time": time.ctime()}
    for tool in ("obabel", "unidock", "mk_prepare_ligand.py", "mk_prepare_receptor.py",
                 "prepare_receptor4.py", "reduce", "pdb2pqr30"):
        c[tool] = shutil.which(tool)
    for mod in ("rdkit", "meeko"):
        try:
            m = __import__(mod); c[mod] = getattr(m, "__version__", "unknown")
        except Exception as e:
            c[mod] = f"ERROR {e}"
    r = run("obabel -V"); c["obabel_version"] = (r.stdout or r.stderr).strip()[:80]
    r = run("unidock --version"); c["unidock_version"] = ((r.stdout or "") + (r.stderr or "")).strip()[:200]
    r = run("nvidia-smi --query-gpu=index,name,memory.used --format=csv,noheader")
    c["gpus"] = (r.stdout or "").strip().splitlines()
    try:
        from meeko import MoleculePreparation  # noqa
        c["meeko_api"] = "MoleculePreparation importable"
        try:
            from meeko import PDBQTWriterLegacy  # noqa
            c["meeko_writer"] = "PDBQTWriterLegacy"
        except Exception:
            c["meeko_writer"] = "legacy write_pdbqt_string"
    except Exception as e:
        c["meeko_api"] = f"ERROR {e}"
    (OUT / "capabilities.json").write_text(json.dumps(c, indent=1))
    print(json.dumps(c, indent=1))
    return c


def fetch_pdb(code, dest):
    if dest.exists(): return dest
    for url in (f"https://files.rcsb.org/download/{code}.pdb",
                f"https://files.rcsb.org/download/{code.lower()}.pdb"):
        r = run(f"curl -sL -m 120 -o {dest} '{url}'")
        if dest.exists() and dest.stat().st_size > 10000: return dest
    raise SystemExit(f"could not fetch {code}")


def prep_receptor(target):
    code, chain = PDB[target]
    wd = WORK / target; wd.mkdir(parents=True, exist_ok=True)
    raw = fetch_pdb(code, wd / f"{code}.pdb")
    # keep one protein chain, drop waters and heteroatoms
    keep = []
    for line in raw.read_text(errors="replace").splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            keep.append(line)
        elif line.startswith("TER") and len(line) > 21 and line[21] == chain:
            keep.append(line)
    clean = wd / "receptor_clean.pdb"
    clean.write_text("\n".join(keep) + "\nEND\n")
    rep = {"target": target, "pdb": code, "chain": chain, "n_atom_records": len(keep)}
    # naive route: OpenBabel, rigid receptor
    nai = wd / "receptor_naive.pdbqt"
    r = run(f"obabel {clean} -O {nai} -xr -p 7.4", timeout=900)
    rep["naive_ok"] = nai.exists() and nai.stat().st_size > 1000
    rep["naive_log"] = ((r.stdout or "") + (r.stderr or ""))[-300:]
    # skill route: Meeko receptor prep if available, else declare the substitution
    skl = wd / "receptor_skill.pdbqt"
    mk = shutil.which("mk_prepare_receptor.py")
    if mk:
        # installed Meeko CLI: --read_pdb, -o basename, -p writes <basename>_rigid.pdbqt,
        # -x deletes residues that fail templating, charges computed with gasteiger
        base = wd / "receptor_skill_mk"
        cmd = (f"{mk} --read_pdb {clean} -o {base} -p -x "
               f"--default_altloc A --compute_charges --charge_model gasteiger")
        r = run(cmd, timeout=1800)
        rep["skill_entrypoint"] = "mk_prepare_receptor.py"
        rep["skill_cmd"] = cmd
        rep["skill_log"] = ((r.stdout or "") + (r.stderr or ""))[-600:]
        cand = sorted(wd.glob("receptor_skill_mk*rigid*.pdbqt")) or sorted(wd.glob("receptor_skill_mk*.pdbqt"))
        if cand: shutil.copy(cand[0], skl)
    if not (skl.exists() and skl.stat().st_size > 1000):
        pr4 = shutil.which("prepare_receptor4.py")
        if pr4:
            r = run(f"{pr4} -r {clean} -o {skl} -A hydrogens -U nphs_lps", timeout=1800)
            rep["skill_entrypoint"] = "prepare_receptor4.py"
            rep["skill_log"] = ((r.stdout or "") + (r.stderr or ""))[-500:]
    if not (skl.exists() and skl.stat().st_size > 1000):
        shutil.copy(nai, skl)
        rep["skill_entrypoint"] = "FALLBACK OpenBabel (Meeko/ADT receptor prep unavailable)"
        rep["receptor_routes_identical"] = True
    rep["skill_ok"] = skl.exists() and skl.stat().st_size > 1000
    for k, p in (("naive_sha256", nai), ("skill_sha256", skl)):
        if p.exists(): rep[k] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    (OUT / f"receptor_{target}.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
    return rep


def library(target):
    p = DATA / "targets" / target / "03_library_preparation" / "library_labels.csv"
    rows = []
    with open(p, newline="") as f:
        for row in csv.DictReader(f):
            rows.append((row["compound_id"], row["smiles"], 1 if row["label"].strip().lower().startswith("act") else 0))
    return rows


def _lig_naive(args):
    cid, smi, outdir, timeout = args
    p = Path(outdir) / f"{cid}.pdbqt"
    if p.exists() and p.stat().st_size > 50: return (cid, "ok", "")
    try:
        r = subprocess.run(["obabel", f"-:{smi}", "-opdbqt", "-O", str(p), "--gen3d", "-p", "7.4"],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return (cid, "timeout", "gen3d timeout")
    if p.exists() and p.stat().st_size > 50: return (cid, "ok", "")
    return (cid, "fail", ((r.stdout or "") + (r.stderr or ""))[-120:])


def _lig_skill(args):
    cid, smi, outdir, timeout = args
    p = Path(outdir) / f"{cid}.pdbqt"
    if p.exists() and p.stat().st_size > 50: return (cid, "ok", "")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation
        m = Chem.MolFromSmiles(smi)
        if m is None: return (cid, "fail", "rdkit parse")
        m = Chem.AddHs(m)
        ps = AllChem.ETKDGv3(); ps.randomSeed = 0xf00d
        if AllChem.EmbedMolecule(m, ps) != 0:
            ps.useRandomCoords = True
            if AllChem.EmbedMolecule(m, ps) != 0: return (cid, "fail", "embed")
        try: AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception: pass
        prep = MoleculePreparation()
        setups = prep(m)
        setup = setups[0] if isinstance(setups, (list, tuple)) else setups
        try:
            from meeko import PDBQTWriterLegacy
            s = PDBQTWriterLegacy.write_string(setup)
            s = s[0] if isinstance(s, tuple) else s
        except Exception:
            s = prep.write_pdbqt_string()
        if not s or len(s) < 50: return (cid, "fail", "empty pdbqt")
        p.write_text(s)
        return (cid, "ok", "")
    except Exception as e:
        return (cid, "fail", f"{type(e).__name__}: {e}"[:120])


def prep_ligands(target, limit=None, nproc=None, timeout=60):
    rows = library(target)
    if limit: rows = rows[:limit]
    rep = {"target": target, "n_requested": len(rows), "timeout_s": timeout, "routes": {}}
    for route, fn in (("naive", _lig_naive), ("skill", _lig_skill)):
        od = WORK / target / f"ligands_{route}"; od.mkdir(parents=True, exist_ok=True)
        args = [(cid, smi, str(od), timeout) for cid, smi, _ in rows]
        t0 = time.time(); res = []
        with cf.ProcessPoolExecutor(max_workers=nproc or max(1, os.cpu_count() - 4)) as ex:
            for i, out in enumerate(ex.map(fn, args, chunksize=8), 1):
                res.append(out)
                if i % 500 == 0: print(f"  {target} {route}: {i}/{len(args)}", flush=True)
        st = {}
        for _, s, _ in res: st[s] = st.get(s, 0) + 1
        rep["routes"][route] = {"status_counts": st, "wall_s": round(time.time() - t0, 1),
                                "n_files": len(list(od.glob('*.pdbqt')))}
        fails = [(c, m) for c, s, m in res if s != "ok"]
        (OUT / f"ligprep_fail_{target}_{route}.json").write_text(json.dumps(fails[:500], indent=1))
    (OUT / f"ligprep_{target}.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
    return rep


def parse_scores(outdir):
    best = {}
    for p in Path(outdir).glob("*.pdbqt"):
        v = None
        for line in p.read_text(errors="replace").splitlines():
            if line.startswith("REMARK VINA RESULT"):
                try: x = float(line.split()[3])
                except Exception: continue
                v = x if v is None else min(v, x)
        if v is not None: best[p.stem.replace("_out", "")] = v
    return best


def dock(target, cond, limit=None, gpu=None, batch=500, retry=20, dry=False, mode="gpu_batch"):
    prep_route, par_name = COND[cond]; par = PARAMS[par_name]
    wd = WORK / target
    rec = wd / f"receptor_{prep_route}.pdbqt"
    ligdir = wd / f"ligands_{prep_route}"
    bc = json.loads((DATA / "targets" / target / "box_center.json").read_text())
    ctr = bc.get("center") or [bc["center_x"], bc["center_y"], bc["center_z"]]
    ligs = sorted(p for p in ligdir.glob("*.pdbqt") if p.stat().st_size > 50)
    if limit: ligs = ligs[:limit]
    od = wd / f"dock_{cond}"; od.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    if gpu is not None: env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    man = {"target": target, "condition": cond, "preparation_route": prep_route,
           "parameters": par_name, "box": par["box"], "exhaustiveness": par["exh"],
           "num_modes": NUM_MODES, "center": ctr, "receptor": str(rec),
           "receptor_sha256": hashlib.sha256(rec.read_bytes()).hexdigest()[:16],
           "n_ligands": len(ligs), "gpu": gpu, "invocation": mode, "batch": batch,
           "commands": [], "started": time.ctime()}
    def one(chunk, tag):
        if mode == "gpu_batch":
            head = ["unidock", "--receptor", str(rec), "--gpu_batch"] + [str(p) for p in chunk]
        else:
            idx = od / f"idx_{tag}.txt"; idx.write_text("\n".join(str(p) for p in chunk) + "\n")
            head = ["unidock", "--receptor", str(rec), "--ligand_index", str(idx)]
        cmd = head + [
               "--center_x", f"{ctr[0]:.3f}", "--center_y", f"{ctr[1]:.3f}", "--center_z", f"{ctr[2]:.3f}",
               "--size_x", f"{par['box']:.1f}", "--size_y", f"{par['box']:.1f}", "--size_z", f"{par['box']:.1f}",
               "--exhaustiveness", str(par["exh"]), "--num_modes", str(NUM_MODES), "--dir", str(od)]
        man["commands"].append(" ".join(cmd[:8]) + f" ... ({len(chunk)} ligands)")
        if dry: print("  DRY", " ".join(cmd)); return 0
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        (od / f"log_{tag}.txt").write_text((r.stdout or "") + "\n--stderr--\n" + (r.stderr or ""))
        return r.returncode
    t0 = time.time()
    for i in range(0, len(ligs), batch):
        ch = ligs[i:i + batch]
        if one(ch, f"b{i//batch:04d}") != 0:
            for j in range(0, len(ch), retry): one(ch[j:j + retry], f"b{i//batch:04d}r{j//retry}")
        if (i // batch) % 5 == 0: print(f"  {target} {cond}: {min(i+batch,len(ligs))}/{len(ligs)} {time.time()-t0:.0f}s", flush=True)
    man["wall_s"] = round(time.time() - t0, 1)
    sc = {} if dry else parse_scores(od)
    lab = {c: l for c, _, l in library(target)}
    rows, miss = [], 0
    for p in ligs:
        cid = p.stem; s = sc.get(cid)
        if s is None: miss += 1
        rows.append({"mol_id": cid, "label": lab.get(cid, ""), "score": "" if s is None else f"{s:.4f}",
                     "status": "no_score" if s is None else "ok"})
    man["n_scored"] = len(ligs) - miss; man["n_no_score"] = miss; man["finished"] = time.ctime()
    f = OUT / f"scores_{target}_{cond}.csv"
    with open(f, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["mol_id", "label", "score", "status"]); w.writeheader(); w.writerows(rows)
    (OUT / f"manifest_{target}_{cond}.json").write_text(json.dumps(man, indent=1))
    print(f"  {target} {cond}: {man['n_scored']}/{len(ligs)} scored in {man['wall_s']}s")
    return man


def bench(target, n=200, gpu=None):
    out = {}
    for mode in ("ligand_index", "gpu_batch"):
        for b in (100, 500):
            if mode == "ligand_index" and b == 500: continue
            t0 = time.time()
            m = dock(target, "A", limit=n, gpu=gpu, batch=b, mode=mode)
            dt = time.time() - t0
            out[f"{mode}_b{b}"] = {"wall_s": round(dt, 1), "per_ligand_s": round(dt / max(1, n), 3),
                                   "scored": m["n_scored"]}
            print(f"  {mode} batch={b}: {dt:.0f}s for {n} ligands "
                  f"({dt/max(1,n):.2f} s/ligand), scored {m['n_scored']}", flush=True)
    (OUT / f"bench_{target}.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("target", nargs="?"); ap.add_argument("cond", nargs="?")
    ap.add_argument("--limit", type=int); ap.add_argument("--gpu", type=int)
    ap.add_argument("--nproc", type=int); ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--mode", default="gpu_batch", choices=["gpu_batch", "ligand_index"])
    ap.add_argument("--batch", type=int, default=500)
    ap.add_argument("-n", type=int, default=200)
    a = ap.parse_args()
    if a.cmd == "caps": caps()
    elif a.cmd == "prep-receptor": prep_receptor(a.target)
    elif a.cmd == "prep-ligands": prep_ligands(a.target, a.limit, a.nproc, a.timeout)
    elif a.cmd == "dock": dock(a.target, a.cond.upper(), a.limit, a.gpu, batch=a.batch, dry=a.dry, mode=a.mode)
    elif a.cmd == "bench": bench(a.target, a.n, a.gpu)
    elif a.cmd == "all":
        prep_receptor(a.target); prep_ligands(a.target, a.limit, a.nproc, a.timeout)
        for cnd in ("A", "D", "C", "E"): dock(a.target, cnd, a.limit, a.gpu, batch=a.batch, mode=a.mode)
    else: ap.error("unknown command")
