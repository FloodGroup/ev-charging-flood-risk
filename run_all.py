#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run the analysis pipeline, module by module.

Usage:
  python run_all.py standard          recompute every statistic from the stored
                                      intermediate tables and run all
                                      verification scripts (about 20 minutes)
  python run_all.py full              also run the steps that read the raw
                                      inputs. Checkpointed steps resume where
                                      they stopped, so with all checkpoints in
                                      place this takes hours, not days
  python run_all.py full --fresh      first move every stored checkpoint aside
                                      (into data/_previous_<time>), so that
                                      everything is recomputed from the raw
                                      inputs (the complete rebuild)
  python run_all.py figures           redraw composites from their figure data
  python run_all.py standard fig3     only the named module(s): stage0, fig1,
                                      fig2, fig3, fig4, fig5, si, reanalysis

A step marked "loop" is called again until its completion marker appears in
the output. The analysis stops at the first failing step. Every step writes
its full output to a log file; a summary table is printed at the end.
"""
import glob
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "shared"))
import paths as P  # noqa: E402


LOGS = os.path.join(HERE, "run_logs")
FIGS = P.FIGURES_ROOT
RPS = ["RP010", "RP020", "RP050", "RP075", "RP100", "RP200", "RP500"]


def S(script, mode="standard", variants=((),), loop=None, max_loops=1, full_extra=()):
    return dict(script=script, mode=mode, variants=[list(v) for v in variants],
                loop=loop, max_loops=max_loops, full_extra=list(full_extra))


# Checkpoints of the heavy steps, relative to <module>/data. With --fresh they
# are moved aside so that the steps start again from the raw inputs.
CHECKPOINTS = {
    "fig2_flood_exposure": ["sample_index.npz", "depth_RP*.npy", "metro_p*.npy"],
    "fig3_opening_year_cohorts": ["checkpoints"],
    "fig4_grid_cascade": ["checkpoints"],
    "fig5_inequality_protection": ["checkpoints", "fig05_station_covariates.parquet"],
}


PIPELINE = [
    ("stage0_station_database", [
        S("01_convert_amap_coordinates.py", "full"),
        S("02_merge_charging_networks.py", "full"),
        S("03_assign_countries.py", "full"),
        S("04_export_geopackage.py", "full"),
        S("05_verify_database.py"),
    ]),
    ("fig2_flood_exposure", [
        S("00_build_sampling_index.py", "full"),
        S("01_sample_flood_depths.py", "full", loop="ALL RPS DONE", max_loops=200),
        S("02_assemble_station_depths.py", "full"),
        S("03_assign_urban_centres.py", "full", loop="ALL METRO PARTS DONE", max_loops=8),
        S("04_compute_core_statistics.py"),
        S("05_compute_city_kitagawa.py"),
        S("06_compute_stock_tables.py"),
        S("08_verify_results.py"),
    ]),
    ("fig3_opening_year_cohorts", [
        S("00_ingest_ocm_archive.py", "full"),
        S("01_sample_flood_depths.py", "full", variants=[[rp] for rp in RPS] + [["assemble"]],
          loop="DONE", max_loops=60),
        S("02_assign_urban_centres.py", "full"),
        S("03_compute_cohort_statistics.py"),
        S("04_extract_city_maps.py", "full", variants=[["2046", "28"], ["1677", "24"]]),
        S("07_leave_one_operator_out.py"),
        S("06_verify_results.py"),
    ]),
    ("fig4_grid_cascade", [
        S("00_extract_power_nodes.py", "full",
          variants=[["sub_pt"], ["sub_poly"], ["plant"], ["wri"], ["merge"]]),
        S("01_sample_substation_depths.py", "full", variants=[[rp] for rp in RPS] + [["assemble"]],
          loop="DONE", max_loops=120),
        S("02_assign_substation_countries.py", "full", loop="DONE", max_loops=60),
        S("03_assign_chargers_to_substations.py", "full", loop="DONE", max_loops=60),
        S("04_build_grid_graph.py", "full"),
        S("05_attach_assets.py", "full"),
        S("06_network_feeders.py", "full"),
        S("07_probabilistic_cascade.py"),
        S("08_islanding_check.py"),
        S("09_extract_city_vignettes.py", full_extra=["full"],
          variants=[["2288", "12"], ["12703", "15.5"], ["10912", "18"], ["1048", "13"]]),
        S("11_export_grid_model.py", "full"),
        S("12_verify_results.py"),
        S("13_test_site_failure_topology.py"),
    ]),
    ("fig5_inequality_protection", [
        S("01_build_station_covariates.py", "full", loop="DONE", max_loops=8),
        S("02_compute_statistics.py"),
        S("04_verify_results.py"),
    ]),
    ("fig1_global_infrastructure", [
        S("01_compute_statistics.py"),
        S("03_verify_results.py"),
    ]),
    ("si_tables", [
        S("01_build_si_tables.py"),
        S("02_verify_si_tables.py"),
    ]),
    ("reanalysis_products", [
        S("01_build_station_density_grid.py"),
    ]),
]


def run(cmd, cwd, log, loop=None, max_loops=1):
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", MPLBACKEND="Agg")
    out = ""
    for k in range(max(1, max_loops)):
        p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = p.stdout.decode("utf-8", "replace")
        with open(log, "a", encoding="utf-8") as f:
            f.write(f"\n===== {' '.join(os.path.basename(c) if os.sep in c else c for c in cmd[1:])}"
                    f" (call {k + 1}) =====\n{out}")
        if p.returncode != 0:
            return False, out
        if not loop or loop in out:
            return True, out
    return False, out + f"\n[marker '{loop}' not reached after {max_loops} calls]"


def last_line(out):
    lines = [ln for ln in out.strip().splitlines() if ln.strip()]
    return lines[-1][:70] if lines else ""


def missing_map_inputs(script):
    """List local-only layers required by this composite."""
    folder = os.path.dirname(script)
    figure = os.path.basename(folder)
    if figure == "Figure3":
        city_ids = (2046, 1677)
        layer = "citymap"
    elif figure == "Figure4":
        city_ids = (2288, 12703, 10912, 1048)
        layer = "feeder_map"
    else:
        return []
    private = os.path.join(folder, "data", "private")
    names = [f"{layer}_{city}.npz" for city in city_ids]
    names += [f"{kind}_{city}.npz" for city in city_ids for kind in ("relief", "basemap")]
    return [item for item in names if not os.path.isfile(os.path.join(private, item))]


def selected(module, only):
    return not only or module in only or module.split("_")[0] in only


def set_checkpoints_aside(stamp, only):
    for module, patterns in CHECKPOINTS.items():
        if not selected(module, only):
            continue
        data = P.get(module.split("_", 1)[0].upper() + "_DATA")
        hits = sorted(p for pat in patterns for p in glob.glob(os.path.join(data, pat)))
        if not hits:
            continue
        dest = os.path.join(data, f"_previous_{stamp}")
        os.makedirs(dest, exist_ok=True)
        for h in hits:
            shutil.move(h, os.path.join(dest, os.path.basename(h)))
        print(f"[fresh] {module}: {len(hits)} checkpoint item(s) moved to {dest}", flush=True)


def main():
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except Exception:
            pass
    args = [a for a in sys.argv[1:] if a != "--fresh"]
    fresh = "--fresh" in sys.argv[1:]
    mode = args[0] if args else "standard"
    only = set(args[1:])
    if mode not in ("standard", "full", "figures") or (fresh and mode != "full"):
        raise SystemExit(__doc__)
    os.makedirs(LOGS, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    rows, ok_all = [], True
    print(f"mode        : {mode}")
    print(f"package     : {HERE}")
    print(f"input data  : {P.SOURCE_ROOT}")
    print(f"figures     : {FIGS}")
    print(f"logs        : {LOGS}\n", flush=True)

    if fresh:
        set_checkpoints_aside(stamp, only)

    if mode in ("standard", "full"):
        for module, steps in PIPELINE:
            if not selected(module, only):
                continue
            code = os.path.join(HERE, "analysis", module)
            log = os.path.join(LOGS, f"{stamp}_{module}.log")
            for st in steps:
                if st["mode"] == "full" and mode != "full":
                    continue
                for v in st["variants"]:
                    v = v + (st["full_extra"] if mode == "full" else [])
                    label = f"{module}/{st['script']} {' '.join(v)}".strip()
                    print(f"[run ] {label}", flush=True)
                    t0 = time.time()
                    ok, out = run([sys.executable, st["script"]] + v, code, log, st["loop"], st["max_loops"])
                    rows.append((label, "ok" if ok else "FAILED", time.time() - t0))
                    print(f"[{'ok  ' if ok else 'FAIL'}] {label}  ({time.time() - t0:.0f} s)  {last_line(out)}", flush=True)
                    if not ok:
                        ok_all = False
                        break
                if not ok_all:
                    break
            if not ok_all:
                break

    draw = mode == "figures" or (mode == "full" and not only)
    if ok_all and draw:
        log = os.path.join(LOGS, f"{stamp}_figures.log")
        if not FIGS or not os.path.isdir(FIGS):
            rows.append(("locate figures directory", "FAILED", 0.0))
            print("Figures directory missing. Check EVFR_FIGURES_ROOT.")
            ok_all = False
        else:
            prep = []
            if mode == "full":
                prep.append(("map backgrounds", [sys.executable, os.path.join("analysis", "map_backgrounds", "01_build_map_backgrounds.py"), FIGS]))
                prep.append(("export figure source data", [sys.executable, "export_figure_source_data.py", FIGS, "--local-maps"]))
            for label, cmd in prep:
                print(f"[run ] {label}", flush=True)
                t0 = time.time()
                ok, out = run(cmd, HERE, log)
                rows.append((label, "ok" if ok else "FAILED", time.time() - t0))
                print(f"[{'ok  ' if ok else 'FAIL'}] {label}  ({time.time() - t0:.0f} s)  {last_line(out)}", flush=True)
                ok_all = ok_all and ok
            scripts = sorted(glob.glob(os.path.join(FIGS, "*", "make_*.py")))
            for script in scripts if ok_all else []:
                folder = os.path.dirname(script)
                label = "draw " + os.path.relpath(script, FIGS)
                missing = missing_map_inputs(script)
                if missing:
                    status = "FAILED" if mode == "full" else "skipped"
                    rows.append((label, status, 0.0))
                    if mode == "full":
                        ok_all = False
                    print(f"[{status}] {label}  (local map inputs: {', '.join(missing)})", flush=True)
                    continue
                print(f"[run ] {label}", flush=True)
                t0 = time.time()
                ok, out = run([sys.executable, os.path.basename(script)], folder, log)
                rows.append((label, "ok" if ok else "FAILED", time.time() - t0))
                print(f"[{'ok  ' if ok else 'FAIL'}] {label}  ({time.time() - t0:.0f} s)", flush=True)
                ok_all = ok_all and ok

    print("\n" + "=" * 78)
    for label, status, sec in rows:
        print(f"{status:7s} {sec:7.0f} s  {label}")
    print("=" * 78)
    print("RESULT: " + ("ALL STEPS COMPLETED" if ok_all else
                        f"AT LEAST ONE STEP FAILED -- see the log files in {LOGS}"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
