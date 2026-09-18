# C2K-GRO — Coupled Cotton2K–GroIMP Cotton Growth Model

Daily bidirectional coupling of the process-based crop model **Cotton2K** and the
functional–structural plant modelling platform **GroIMP**, developed to simulate the
growth of cotton (cultivar **Tahe-2 / 塔河2号**) under film-mulched drip irrigation in
the Alar reclamation area, southern Xinjiang, China.

This repository accompanies the manuscript *"Cotton Growth Simulation in the Alar
Reclamation Area Based on Dynamic Bidirectional Process–Structure Coupling"*
(submitted to *Frontiers in Plant Science*; citation details will be added once the
paper is published — see [Citation](#citation) below).

## What C2K-GRO does

Standalone Cotton2K estimates canopy light interception with the Beer–Lambert law,
which does not capture the self-shading created by dense, wide–narrow-row canopies
and, in this system, underestimates lint yield by 44–50%. C2K-GRO replaces that
estimate with **FluxLM**, GroIMP's Monte-Carlo 3-D ray tracer, and exchanges state
with Cotton2K once per simulated day over a local TCP socket:

- **GroIMP → Cotton2K**: 3-D canopy light interception, structural LAI, boll count
- **Cotton2K → GroIMP**: phenology, carbon allocation, leaf weight (drives geometry update)

Coupling runs in three stages: Cotton2K alone from sowing to emergence (Stage A),
full bidirectional exchange from emergence to boll-opening (Stage B, ~142 days), and
Cotton2K alone again from boll-opening to harvest (Stage C). A carbon "Bridge" damps
the cold-start of the 3-D canopy during early seedling growth, and an LAI-consistency
"Nudge" softly reconciles Cotton2K's carbon-derived LAI with GroIMP's structural LAI
whenever they diverge by more than 15%.

The model was evaluated by **bidirectional cross-validation** between two contrasting
growing seasons — 2022 (low radiation) and 2025 (high radiation) — each serving in
turn as the calibration year for the other's validation.

## Repository structure

```
2022_calibration/   2022 data, model calibrated on 2022
2022_validation/    2022 data, model calibrated on 2025 (independent, out-of-sample)
2025_calibration/   2025 data, model calibrated on 2025
2025_validation/    2025 data, model calibrated on 2022 (independent, out-of-sample)
```

**Important — same file names hold different content across folders.** In
particular, `INERTIA_COEF` inside `cotton2k_daily.py` is `0.001843` in the two 2022
folders and `0.00062` in the two 2025 folders (the two independently calibrated
values discussed in the paper, §3.4 and Supplementary Table S2); the `_validation`
folders intentionally *freeze* the other year's calibrated value to test cross-year
transfer. Do not assume a file is interchangeable between folders just because it
has the same name.

### Files inside each folder

| File | Contents |
|---|---|
| `weather20{22,25}.csv` | Daily observed weather forcing for that season |
| `cotton2k_daily.py` | Cotton2K process-model core (phenology, carbon allocation, yield formation) |
| `cotton2k_standalone*.py` | Standalone Cotton2K driver (Beer–Lambert only, no GroIMP) — the baseline used for comparison in the paper |
| `daily_coupling_server2_20{22,25}.py` | Coupling server: runs Cotton2K day-by-day and exchanges state with GroIMP over a TCP socket |
| `CottonModel_20{22,25}_Coupled.rgg` | GroIMP/XL script defining the 3-D cotton canopy and FluxLM light calculation; connects to the Python server above |
| `output_standalone.csv` | Daily output of the standalone Cotton2K run |
| `output_nudge_ON.csv` / `output_nudge_OFF.csv` | Daily output of the coupled model with/without the LAI-Nudge mechanism (ablation pair) |
| `coupling_log_ON.csv` / `coupling_log_OFF.csv` | Daily exchanged variables between Cotton2K and GroIMP |
| `groimp_feedback.csv` | Raw daily feedback sent from GroIMP to Cotton2K |
| `nudge_detail_log_ON.csv` | Per-trigger detail whenever the Nudge mechanism fired |
| `nudge_compare_summary_20{22,25}.csv` | Season-end accuracy summary (NSE, RMSE, nRMSE, MBE, yield error, harvest efficiency) for NudgeON/OFF |

`2025_calibration/` additionally contains two files not present in the other three
folders:

| File | Contents |
|---|---|
| `calibrate_inertia_300mm.py` | Standalone recalibration utility: reads a coupled run's lint yield, compares it to the measured target, and (if outside ±2%) proposes/writes a new `INERTIA_COEF` into `cotton2k_daily.py`, logging every check to `calibration_history.csv` |
| `calibration_history.csv` | Log produced by the script above. Currently holds the single 2025 convergence check cited in Supplementary Table S2 (2026-09-17T10:03:16, error +1.06%, converged, no change needed). The corresponding 2022 check (cited in the same table, 2026-09-17T10:03:38) is not yet included here — see the note in the paper/response letter if that record has not been located. |

## Requirements

- Python 3.9+ (standard library only for the coupling server; `pandas`/`numpy` for
  post-processing the CSV outputs)
- [GroIMP](https://long.grogra.de/) (third-party, GPL-licensed software; not
  redistributed here) with the FluxLM ray-tracing module, to execute the `.rgg`
  scripts

## Running a simulation

The Python side runs as the coupling server; GroIMP connects to it as a client
while executing the matching `.rgg` script. From inside one of the four folders:

```bash
python daily_coupling_server2_20{22,25}.py \
    --weather weather20{22,25}.csv \
    --port 9527 \
    --output runs/my_run
```

Then open the corresponding `CottonModel_20{22,25}_Coupled.rgg` in GroIMP and run
it — it connects to `localhost:9527` and steps forward one simulated day per
exchange. Pass `--nudge-off` to the Python server to reproduce the NudgeOFF
ablation (GroIMP receives a `nudge_enabled=false` flag each day and the two models
evolve without the soft LAI constraint).

Field weather and management inputs used to force the model are documented in
Supplementary Table S3 of the paper; note that the simplified irrigation forcing
built into the model differs somewhat from the raw field irrigation record (both
are reported transparently in the Supplementary Material).

## Reproducing the paper's numbers

`nudge_compare_summary_*.csv` in each folder contains the season-end NSE/RMSE/yield
figures cited in the paper's Tables 2 and Figures 2–5. `coupling_log_ON.csv` and
`output_nudge_ON.csv` contain the daily boll count, boll weight and lint-yield
trajectories used to derive the harvest-efficiency figures in §3.4.

## License

- **Code** (`*.py`, `*.rgg`) is released under the [MIT License](LICENSE).
- **Data** (`*.csv` weather, management and output files) is released under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — please cite the paper
  above when reusing.

## Citation

If you use this code or data, please cite the accompanying paper (details to be
added on publication) and/or this repository via its Zenodo DOI (see the badge on
this page, or `CITATION.cff`, once archived).

## Contact

Corresponding author: Zhenqi Fan (tarimfanzq@163.com), College of Information
Engineering / Key Laboratory of Oasis Agriculture (Ministry of Education), Tarim
University, Alar, China.
