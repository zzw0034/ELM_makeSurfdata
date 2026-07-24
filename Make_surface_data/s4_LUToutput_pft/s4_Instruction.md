# Instruction for s4_LUToutput_pft

```bash
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft
NLCD2ELM and 
LUToutput
```

The `NLCD2ELM` and `LUToutput` folders contain the workflows I used to generate `4km surfdata.nc` and `landuse.timeseries.nc`.
For details, refer to

```bash
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/NLCD2ELM/README_NLCD2ELM_workflow.md
```

I keep these two folders as references and use them as the basis for preparing the `s4_LUToutput_pft` procedure.

## 1. Convert NLCD to ELM PFT

**Notes**: Dan generated AI landuse from NLCD in 
>/projects/hpcl-cli185/proj-shared/zdr/hires_data/AI_nlcd/nlcd_frac_pred_1850-2023_1_24deg.nc

Run **s4_1_convert_nlcd_to_PFT.py**.

```bash
/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/s4_1_convert_nlcd_to_PFT.py
```

> **Key output:** `elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc`

> output_file = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"

## 2. Generate LUToutput

**mksurfdata** requires input from Land Use Translator (LUT)，details in 
```text
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/doc/ED-Creating land surface datasets for E3SM - Internal-151025-234757.pdf

/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/AIprompt_s4_2_donwscale_LUH2harvest.md
```
**In /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft**
**run** s4_2_donwscale_LUH2harvest.py
> sbatch run_s4_2_downscale_LUH2harvest.slurm

**Generated output:** 
>/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest

LUT output example:
>/projects/hpcl-cli185/world-shared/e3sm/inputdata/lnd/clm2/rawdata/LUT_LUH2_SSP2_RCP45_LUH1f_07222020. 


- I don't use LUT, instead I make python script generate the file with LUT format from NLCD-ELM-PFT data.

- **1/24° 全球或北美** 的 LUT results 给 mksurfdata，用 **Step 1 的输出** `elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc` 作为 LUT 生成脚本（`s4_2_donwscale_LUH2harvest.py`）的输入，再经 mksurfdata 生成目标网格的 surfdata 与 landuse.timeseries。

```bash
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/AIprompt_s4_2_donwscale_LUH2harvest.md
```

## 3. Generate LUToutput single text file

All input files (which are the LUT annual outputs) for a given time series need to be listed in a single text file with one filename per line and the year of the file starting at exactly character 197 (the code expects the first 195 characters to be the filename with full path, then a space, then the year)

**run**  
> /s4_LUToutput_pft/make_fdynuse_list.py

## 4. Recent period 2015–2023: continue harvest with LUH2 v2f SSP2-4.5

LUH2 v2h historical harvest (`transitions.nc`) ends at calendar year **2014**.
For the recent period we no longer freeze harvest at 2014; instead we continue
the series with the LUH2 **v2f SSP2-4.5** scenario:

> `/projects/hpcl-cli185/proj-shared/zw5/luh/ssp2rcp45_transitions.nc`
> (0.25°, 720×1440, same grid / lat-order / harvest variables as v2h; calendar 2015–2099)

**Year mapping** — s4_2 uses the E3SM LUT convention that the file labelled year
`Y` carries the harvest rate of calendar year `Y−1`:

| output file year | harvest calendar year | source |
| --- | --- | --- |
| 1850 … 2015 | 1849 … 2014 | LUH2 v2h `transitions.nc` |
| 2016 … 2023 | 2015 … 2022 | LUH2 v2f SSP2-4.5 `ssp2rcp45_transitions.nc` |

**Rationale**: SSP scenarios diverge appreciably only after ~2030, so over
2015–2022 wood harvest is nearly scenario-independent — which scenario fills these
years barely matters. SSP2-4.5 is chosen as the "near-reality" middle path. This
replaces the previous behaviour of freezing harvest at v2h 2014 for 2015–2023.

> **Implementation status**: decided, not yet coded. As of this note
> `s4_2_donwscale_LUH2harvest.py` still reads only v2h and freezes at 2014 for the
> recent years; the two-file splice above is pending.

### Harvest downscaling corrections (2026-07-23) — the year mapping depends on these

`s4_2_donwscale_LUH2harvest.py` was corrected in three ways:

- **Area conservation**: harvest is distributed conserving harvested *area* within
  each 0.25° coarse cell (`frac_i = F · w_i · ΣA / Σ(w·A)`). The previous version
  conserved the sum of fractions, underestimating harvest by the coarse/fine
  cell-count ratio (~36×). Validated: `LUH2_allocatable / ours` went 35.99 → 1.00.
- **Annual vegetated fraction**: the divisor/weight `PCT_NATVEG` is read per year
  (was frozen at 1850). `LANDMASK` stays static at 1850 on purpose.
- **Time label (k=−1)**: file labelled `Y` carries LUH2 calendar `Y−1`, per the
  E3SM doc convention — this is what makes the year-mapping table above well-defined.

Old (buggy) outputs are kept under `output_downscaled_luh2_harvest/prev_buggy_c0711/`.
