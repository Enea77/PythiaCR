#!/usr/bin/env python3
"""Production summary (for slides): per sample events, jets > 60 GeV, V0 counts, CPU time and memory per job, disk.

  python3 macros/production_summary.py   -> prints a table, writes compare/production_summary.txt

Events / jets / V0s: compare/v2/<sample>.root (genCompare.C v2): hPthatRaw (generated events), hJetEta (jets > 60 GeV,
|eta_WTA| < 2), hV0PtAll_* (all V0s with |eta| < 2, any pT), hInCone_* (V0s with dr < 0.4 of a jet > 60 GeV,
1 < pT < 10 GeV, counted per jet-V0 pair). All counts are generated (unweighted) entries.
CPU / memory: condor job logs (condor/logs/<cluster>.log), jobs that terminated normally.
"""
import os, re, subprocess
import numpy as np
import uproot

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(D, "condor", "logs")
EOS = "/eos/user/e/eprifti/PythiaCR"
# cluster -> (process, [tunes in queue order], jobs per tune, events per job)
CLUSTERS = {
    9724871: ("dijet", ["cuetp8m1", "crblc_mode0", "crblc_mode2", "crblc_mode3"], 100, 100000),
    9726939: ("dijet", ["crblc_mode0"], 2, 100000),
    9726940: ("dijet", ["crblc_mode3"], 1, 100000),
    9726942: ("dijet", ["cp5", "cp5cr1", "cp5cr2"], 100, 100000),
    9747148: ("dijet", ["cuetp8m1_py8212"], 100, 100000),
    9747155: ("dijetppb", ["cuetp8m1", "crblc_mode2", "crblc_mode3"], 200, 50000),
    9747318: ("dijetppbcentral", ["cuetp8m1", "crblc_mode2"], 200, 50000),
    9756651: ("minbiasppb", ["cuetp8m1", "crblc_mode2", "crblc_mode3"], 100, 100000),
    9756652: ("minbiasppbcentral", ["cuetp8m1", "crblc_mode2"], 100, 100000),
}


def parse_log(path):
    """{proc: (cpu seconds, memory MB)} for normally terminated jobs."""
    out, cur = {}, None
    txt = open(path, errors="ignore").read()
    for block in txt.split("...\n"):
        m = re.search(r"^005 \(\d+\.(\d+)\.\d+\)", block, re.M)
        if not m or "Normal termination (return value 0)" not in block:
            continue
        cpu = re.search(r"Usr (\d+) (\d+):(\d+):(\d+), Sys (\d+) (\d+):(\d+):(\d+)\s+-\s+Total Remote Usage", block)
        mem = re.search(r"Memory \(MB\)\s+:\s+(\d+)", block)
        if cpu:
            g = list(map(int, cpu.groups()))
            sec = g[0] * 86400 + g[1] * 3600 + g[2] * 60 + g[3] + g[4] * 86400 + g[5] * 3600 + g[6] * 60 + g[7]
            out[int(m.group(1))] = (sec, int(mem.group(1)) if mem else np.nan)
    return out


usage = {}   # sample -> list of (cpu s, mem MB, events)
for cl, (proc, tunes, nj, nev) in CLUSTERS.items():
    p = os.path.join(LOGS, f"{cl}.log")
    if not os.path.exists(p):
        continue
    for pid, (cpu, mem) in parse_log(p).items():
        t = tunes[pid // nj] if pid // nj < len(tunes) else None
        if t is None:
            continue
        s = f"{proc}_cuetp8m1_py8212" if t == "cuetp8m1_py8212" else f"{proc}_{t}"
        usage.setdefault(s, []).append((cpu, mem, nev))


def disk(sample):
    try:
        r = subprocess.run(["du", "-sb", os.path.join(EOS, sample)], capture_output=True, text=True, timeout=120)
        return int(r.stdout.split()[0]) / 1e9
    except Exception:
        return np.nan


def entries(f, name):
    return f[name].member("fEntries") if name in f else np.nan


rows = []
samples = sorted(set(usage) | {os.path.basename(p)[:-5] for p in
                               __import__("glob").glob(os.path.join(D, "compare", "v2", "dijet*.root"))})
for s in samples:
    p = os.path.join(D, "compare", "v2", s + ".root")
    f = uproot.open(p) if os.path.exists(p) else None
    ev = f["hPthatRaw"].values().sum() if f else np.nan
    if f is None and s in usage:
        ev = sum(u[2] for u in usage[s])
    u = np.array(usage.get(s, []), dtype=float)
    cpu_ev = (u[:, 0].sum() / u[:, 2].sum() * 1e3) if len(u) else np.nan
    rows.append(dict(sample=s, events=ev,
                     jets=entries(f, "hJetEta") if f else np.nan,
                     k0s=entries(f, "hV0PtAll_K0s") if f else np.nan,
                     lam=entries(f, "hV0PtAll_Lam") if f else np.nan,
                     k0s_c=entries(f, "hInCone_K0s") if f else np.nan,
                     lam_c=entries(f, "hInCone_Lam") if f else np.nan,
                     jobs=len(u), cpu_h=u[:, 0].sum() / 3600 if len(u) else np.nan, ms_ev=cpu_ev,
                     mem=np.nanmax(u[:, 1]) if len(u) else np.nan, disk=disk(s)))


def fmt(x, kind="M"):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    return f"{x / 1e6:.1f}M" if kind == "M" else f"{x:.0f}"


hdr = ("| sample | events | jets > 60 GeV | K0S (all) | Λ (all) | K0S in cone | Λ in cone | jobs | CPU (h) | "
       "CPU/event (ms) | max memory (MB) | disk (GB) |")
lines = [hdr, "|" + "---|" * 12]
tot = {k: 0. for k in ("events", "jets", "k0s", "lam", "k0s_c", "lam_c", "cpu_h", "disk")}
for r in rows:
    lines.append(f"| {r['sample']} | {fmt(r['events'])} | {fmt(r['jets'])} | {fmt(r['k0s'])} | {fmt(r['lam'])} | "
                 f"{fmt(r['k0s_c'])} | {fmt(r['lam_c'])} | {r['jobs']} | {fmt(r['cpu_h'], 'n')} | {fmt(r['ms_ev'], 'n')} | "
                 f"{fmt(r['mem'], 'n')} | {fmt(r['disk'], 'n')} |")
    for k in tot:
        if not np.isnan(r[k]):
            tot[k] += r[k]
lines.append(f"| **total** | {fmt(tot['events'])} | {fmt(tot['jets'])} | {fmt(tot['k0s'])} | {fmt(tot['lam'])} | "
             f"{fmt(tot['k0s_c'])} | {fmt(tot['lam_c'])} | | {fmt(tot['cpu_h'], 'n')} | | | {fmt(tot['disk'], 'n')} |")
txt = "\n".join(lines)
print(txt)
open(os.path.join(D, "compare", "production_summary.txt"), "w").write(txt + "\n")
