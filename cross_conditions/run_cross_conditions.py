#!/usr/bin/env python3
"""Cross-condition docking runner for the skill-file virtual-screening benchmark.

Conditions (filter off in both):
  D = skill preparation (ETKDGv3 + Meeko PDBQT)  x  naive parameters (20 A box, exhaustiveness 8)
  E = naive preparation (OpenBabel gen3d)        x  skill parameters (25 A box, exhaustiveness 32)

Only the docking step is new; ligand and receptor PDBQT files are reused from the archive.

Usage
  python run_cross_conditions.py --check
  python run_cross_conditions.py --smoke egfr
  python run_cross_conditions.py --positive-control egfr
  python run_cross_conditions.py --run --targets all --conditions D,E
  python run_cross_conditions.py --run --targets egfr,cdk2 --conditions D --gpu 3

Design decisions and the flags that must be verified against the installed Uni-Dock build
are documented in PROTOCOL.md. The script prints every command it issues.
"""
import argparse, csv, hashlib, json, os, shlex, shutil, subprocess, sys, time
from pathlib import Path

TARGETS = ["egfr","cdk2","cox2","esr1","dpp4","ache","bace1","hsp90","p38","thrombin"]

# Published protocol parameters. Do not edit without updating PROTOCOL.md and the manuscript.
PARAMS = {
    "naive": {"box": 20.0, "exhaustiveness": 8},
    "skill": {"box": 25.0, "exhaustiveness": 32},
}
CONDITIONS = {
    "A": {"prep": "naive", "params": "naive"},   # published
    "C": {"prep": "skill", "params": "skill"},   # published, filter ablated
    "D": {"prep": "skill", "params": "naive"},   # new
    "E": {"prep": "naive", "params": "skill"},   # new
}
BATCH, RETRY_BATCH = 100, 20


def sha256(p, limit=None):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Layout:
    """Resolves archive paths. Override with --root and, if the archive differs, edit here only."""

    def __init__(self, root):
        self.root = Path(root)

    def target_dir(self, t):
        for cand in (self.root / "targets" / t, self.root / t):
            if cand.is_dir():
                return cand
        raise SystemExit(f"target directory not found for {t} under {self.root}")

    def receptor(self, t, prep):
        d = self.target_dir(t)
        pats = [f"receptor_{prep}.pdbqt", f"{prep}/receptor.pdbqt", f"{prep}_receptor.pdbqt",
                f"receptor/{prep}.pdbqt"]
        for p in pats:
            if (d / p).is_file():
                return d / p
        hits = sorted(d.rglob("*.pdbqt"))
        hits = [h for h in hits if "receptor" in h.name.lower() and prep in str(h).lower()]
        if len(hits) == 1:
            return hits[0]
        raise SystemExit(f"receptor pdbqt for {t}/{prep} not resolved; candidates: {hits[:5]}")

    def ligand_dir(self, t, prep):
        d = self.target_dir(t)
        for p in (f"ligands_{prep}", f"{prep}/ligands", f"pdbqt_{prep}", f"{prep}_pdbqt"):
            if (d / p).is_dir():
                return d / p
        raise SystemExit(f"ligand directory for {t}/{prep} not resolved under {d}")

    def box_center(self, t):
        d = self.target_dir(t)
        for p in ("box.json", "docking_box.json", "box_center.json"):
            if (d / p).is_file():
                j = json.loads((d / p).read_text())
                for k in ("center", "box_center"):
                    if k in j:
                        return [float(x) for x in j[k]]
                if all(k in j for k in ("center_x", "center_y", "center_z")):
                    return [float(j["center_x"]), float(j["center_y"]), float(j["center_z"])]
        xtal = None
        for p in ("crystal_ligand.pdb", "ref_ligand.pdb", "cocrystal.pdb"):
            if (d / p).is_file():
                xtal = d / p
                break
        if xtal is None:
            raise SystemExit(f"no box file and no co-crystal ligand for {t}; cannot define the box")
        xs, ys, zs = [], [], []
        for line in xtal.read_text().splitlines():
            if line.startswith("HETATM"):
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
        if not xs:
            raise SystemExit(f"no HETATM records in {xtal}")
        return [sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)]

    def labels(self, t):
        """Returns {mol_id: (chembl_id, label)} with label 1 for actives, 0 for decoys."""
        d = self.target_dir(t)
        out = {}
        for fn, lab in (("actives.csv", 1), ("decoys.csv", 0)):
            for cand in (d / fn, d / "library" / fn):
                if cand.is_file():
                    with open(cand, newline="") as f:
                        for row in csv.DictReader(f):
                            mid = row.get("mol_id") or row.get("id") or row.get("name") or row.get("chembl_id")
                            out[str(mid)] = (row.get("chembl_id", ""), lab)
                    break
        if not out:
            raise SystemExit(f"active/decoy identifier files not found for {t} under {d}")
        return out


def unidock_cmd(binary, receptor, ligand_index, center, box, exhaustiveness, outdir, extra):
    cmd = [binary,
           "--receptor", str(receptor),
           "--ligand_index", str(ligand_index),
           "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
           "--size_x", f"{box:.1f}", "--size_y", f"{box:.1f}", "--size_z", f"{box:.1f}",
           "--exhaustiveness", str(exhaustiveness),
           "--dir", str(outdir)]
    if extra:
        cmd += shlex.split(extra)
    return cmd


def parse_scores(outdir):
    """Best (most negative) Vina score per output pdbqt."""
    scores = {}
    for p in Path(outdir).glob("*.pdbqt"):
        best = None
        for line in p.read_text(errors="replace").splitlines():
            if line.startswith("REMARK VINA RESULT"):
                try:
                    v = float(line.split()[3])
                except (IndexError, ValueError):
                    continue
                best = v if best is None else min(best, v)
            elif line.startswith("REMARK minimizedAffinity"):
                try:
                    v = float(line.split()[-1])
                except ValueError:
                    continue
                best = v if best is None else min(best, v)
        name = p.stem.replace("_out", "")
        if best is not None:
            scores[name] = best
    return scores


def run_condition(lay, binary, target, cond, workdir, gpu=None, limit=None, dry=False, extra=""):
    spec = CONDITIONS[cond]
    prep, par = spec["prep"], PARAMS[spec["params"]]
    receptor = lay.receptor(target, prep)
    ligdir = lay.ligand_dir(target, prep)
    center = lay.box_center(target)
    labels = lay.labels(target)

    ligands = sorted(ligdir.glob("*.pdbqt"))
    if limit:
        ligands = ligands[:limit]
    if not ligands:
        raise SystemExit(f"no ligand pdbqt files in {ligdir}")

    wd = Path(workdir) / f"{target}_{cond}"
    (wd / "out").mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)

    manifest = {
        "target": target, "condition": cond, "preparation_route": prep,
        "box_size": par["box"], "exhaustiveness": par["exhaustiveness"],
        "box_center": center, "receptor": str(receptor),
        "receptor_sha256": sha256(receptor), "ligand_dir": str(ligdir),
        "n_ligands": len(ligands), "batch": BATCH, "retry_batch": RETRY_BATCH,
        "engine_binary": binary, "gpu": gpu, "cuda_visible_devices": env.get("CUDA_VISIBLE_DEVICES"),
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"), "commands": [], "failures": {},
        "ligand_sha256": {},
    }
    for lp in ligands:
        manifest["ligand_sha256"][lp.stem] = sha256(lp)

    def dock(chunk, tag):
        idx = wd / f"index_{tag}.txt"
        idx.write_text("\n".join(str(p) for p in chunk) + "\n")
        cmd = unidock_cmd(binary, receptor, idx, center, par["box"], par["exhaustiveness"], wd / "out", extra)
        manifest["commands"].append(" ".join(cmd))
        print("  $", " ".join(cmd), flush=True)
        if dry:
            return 0
        t0 = time.time()
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        (wd / f"log_{tag}.txt").write_text((r.stdout or "") + "\n--- stderr ---\n" + (r.stderr or ""))
        print(f"    exit {r.returncode} in {time.time()-t0:.0f}s", flush=True)
        return r.returncode

    t_start = time.time()
    for i in range(0, len(ligands), BATCH):
        chunk = ligands[i:i + BATCH]
        if dock(chunk, f"b{i//BATCH:04d}") != 0:
            for j in range(0, len(chunk), RETRY_BATCH):
                dock(chunk[j:j + RETRY_BATCH], f"b{i//BATCH:04d}r{j//RETRY_BATCH}")
    manifest["wall_seconds"] = round(time.time() - t_start, 1)

    scores = {} if dry else parse_scores(wd / "out")
    rows, missing = [], 0
    for lp in ligands:
        mid = lp.stem
        chembl, lab = labels.get(mid, ("", ""))
        if lab == "":
            # identifier mismatch between filenames and label tables is a hard error, not a silent drop
            manifest["failures"].setdefault("unlabelled", []).append(mid)
        s = scores.get(mid)
        if s is None:
            missing += 1
            rows.append({"mol_id": mid, "chembl_id": chembl, "label": lab, "score": "", "status": "no_score"})
        else:
            rows.append({"mol_id": mid, "chembl_id": chembl, "label": lab, "score": f"{s:.4f}", "status": "ok"})
    manifest["n_scored"] = len(ligands) - missing
    manifest["n_no_score"] = missing
    manifest["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    out_csv = Path(workdir) / f"scores_{target}_{cond}.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mol_id", "chembl_id", "label", "score", "status"])
        w.writeheader(); w.writerows(rows)
    (Path(workdir) / f"manifest_{target}_{cond}.json").write_text(json.dumps(manifest, indent=1))

    print(f"  {target} {cond}: {manifest['n_scored']}/{len(ligands)} scored, "
          f"{missing} without a score, {manifest.get('wall_seconds','-')}s")
    if manifest["failures"].get("unlabelled"):
        print(f"  WARNING {len(manifest['failures']['unlabelled'])} ligands have no label table entry")
    return out_csv


def cmd_check(lay, binary):
    print(f"archive root: {lay.root}")
    ok = True
    print(f"{'target':10s} {'naive lig':>10s} {'skill lig':>10s} {'shared':>8s} {'box centre':>34s}")
    for t in TARGETS:
        try:
            n = {p.stem for p in lay.ligand_dir(t, "naive").glob("*.pdbqt")}
            s = {p.stem for p in lay.ligand_dir(t, "skill").glob("*.pdbqt")}
            c = lay.box_center(t)
            lay.receptor(t, "naive"); lay.receptor(t, "skill")
            lay.labels(t)
            print(f"{t:10s} {len(n):10d} {len(s):10d} {len(n & s):8d} "
                  f"[{c[0]:9.3f},{c[1]:9.3f},{c[2]:9.3f}]")
            if len(n & s) < 0.5 * max(len(n), len(s)):
                print(f"  NOTE {t}: the two preparation routes share few ligand files; "
                      f"condition E will be limited to the naive-attempted set (PROTOCOL.md section 2)")
        except SystemExit as e:
            ok = False
            print(f"{t:10s} FAILED: {e}")
    print("\nengine:", shutil.which(binary) or "NOT FOUND")
    if shutil.which(binary):
        r = subprocess.run([binary, "--help"], capture_output=True, text=True)
        txt = (r.stdout or "") + (r.stderr or "")
        for flag in ("--ligand_index", "--gpu_batch", "--exhaustiveness", "--size_x", "--dir", "--search_mode"):
            print(f"  {flag:16s} {'present' if flag in txt else 'ABSENT - adjust unidock_cmd()'}")
    r = subprocess.run(["nvidia-smi", "--query-gpu=index,name,memory.used,memory.total",
                        "--format=csv,noheader"], capture_output=True, text=True)
    print("\nGPUs:\n" + (r.stdout or r.stderr))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("VSBENCH_ROOT", "."),
                    help="benchmark archive root containing targets/")
    ap.add_argument("--workdir", default="cross_conditions_out")
    ap.add_argument("--binary", default=os.environ.get("UNIDOCK_BIN", "unidock"))
    ap.add_argument("--targets", default="all")
    ap.add_argument("--conditions", default="D,E")
    ap.add_argument("--gpu", type=int, default=None)
    ap.add_argument("--extra", default="", help="extra flags passed verbatim to the engine")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--smoke", metavar="TARGET")
    ap.add_argument("--positive-control", metavar="TARGET")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    lay = Layout(a.root)
    Path(a.workdir).mkdir(parents=True, exist_ok=True)

    if a.check:
        sys.exit(cmd_check(lay, a.binary))
    if a.smoke:
        for c in ("D", "E"):
            run_condition(lay, a.binary, a.smoke, c, a.workdir, gpu=a.gpu, limit=300,
                          dry=a.dry_run, extra=a.extra)
        return
    if a.positive_control:
        run_condition(lay, a.binary, a.positive_control, "C", a.workdir, gpu=a.gpu, limit=500,
                      dry=a.dry_run, extra=a.extra)
        print("compare scores_*_C.csv against the archived condition C table for this target")
        return
    if a.run:
        ts = TARGETS if a.targets == "all" else [x.strip() for x in a.targets.split(",")]
        cs = [x.strip().upper() for x in a.conditions.split(",")]
        for t in ts:
            for c in cs:
                print(f"== {t} condition {c}", flush=True)
                run_condition(lay, a.binary, t, c, a.workdir, gpu=a.gpu, dry=a.dry_run, extra=a.extra)
        return
    ap.print_help()


if __name__ == "__main__":
    main()
