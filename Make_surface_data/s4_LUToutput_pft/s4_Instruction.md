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

## 5. 1 km SEUS PFT 图（s4_3，2026-07-29）

**run** `s4_3_nlcd_to_1km_SEUS.py`

```bash
sbatch --export=NONE run_s4_3_nlcd_1km_SEUS.slurm            # 默认 2023
sbatch --export=NONE,YEAR=2015 run_s4_3_nlcd_1km_SEUS.slurm
```

Step 1 的产物是 1/24°（~4 km），**画不出 1 km 的图**。1 km 产品在这里直接从
Dan 重投影后的 Annual NLCD 生成：

> 输入：`/projects/hpcl-cli185/proj-shared/zdr/hires_data/landcover/NLCD/nlcd_<year>_CONUS1.1800deg.nc`
> （变量 `nlcd`，45000 × 108000，dy ≈ 0.00068°、dx ≈ 0.00061°，即 ~70 m）

> 输出：`scr_out/nlcd1km_SEUS/nlcd_elmpft_<year>_SEUS_0.01deg.nc` + 四张 PNG

流程：把 ~70 m 的 NLCD 类别按格点聚合到 **0.01°（~1 km）**，得到每格各 NLCD
类的百分比；再用与 `s4_1_convert_nlcd_to_PFT.py` **完全相同**的映射
（Dan 的 `nlcd_to_elmpft.py`，`fc4 = 0.2`）转成 ELM PFT。

### 与 Dan 的 `regrid_nlcd.py` 的区别

Dan 的脚本按整数像元数分箱，再用最近邻 `interp` 落到目标网格（`_corr` 文件），
因为 `target_res / pixel_size` 不是整数。这里改成**把每个源像元直接归到包含它
中心的目标格子**，一步落在精确的 0.01° 网格上，没有最近邻重采样，也没有被截断
的边缘格子。

### 约定

- **“1 km” = 0.01°**，与项目已有的 `0.01x0.01` 网格命名一致
  （见 `map_0.01x0.01_*_to_SEUS_1_24deg_*`）。SEUS 域内 0.01° 南北 ~1.11 km、
  东西 ~0.9 km。
- **域**取 SEUS_1_24deg 的外边界（`s1_Gridfiles/scr_out/SCRIPgrid_SEUS_1_24deg.nc`，
  504 × 324 @ 1/24°）：lon [-95, -74]，lat [24, 37.5] → 2100 × 1350 @ 0.01°。
  注意 step 1 的 1/24° 产物只到 25°N，而这里的 NLCD 源覆盖到 22.3°N，所以
  24–25°N 那一条在 1 km 产品里是有数据的。

### 交叉校验（2023，lon [-85,-84] × lat [35,36]）

| | 1 km（本脚本） | 1/24°（s4_1 产物） |
|---|---|---|
| PCT_NATVEG 均值 | 90.65 | 91.05 |
| PFT 1 / 7 / 13 / 14 | 15.2 / 52.2 / 25.0 / 4.8 | 15.7 / 50.6 / 25.7 / 5.2 |

差异来自源数据不同（这里是 NLCD 观测直接分箱，s4_1 用的是 Dan 的 AI 预测
NLCD 分数），量级合理。

### 图

1. `fig1_dominant_landcover_*` — 聚合 NLCD 优势类型（Dan 在 `animate_landuse.py`
   里的 8 类分组）。
2. `fig2_dominant_pft_*` — 优势 ELM 自然 PFT，每个 PFT 单独一个图例项。

   **并列（tie）的真实构成**（SEUS 2023，47006 格 = 植被格点的 3.5%）：

   | 并列的 PFT 对 | 占并列 | 来源 |
   |---|---|---|
   | **1 ↔ 7** | **69.9%** | NLCD 43 混交林 + 90 木本湿地，各 50/50 拆给 PFT 1 和 7 |
   | 7 ↔ 13 | 21.4% | NLCD 21 城市开阔地，50/50 拆给 PFT 7 和 13 |
   | 9 ↔ 10 | 7.4% | NLCD 51/52 灌木，50/50 拆 |
   | 其余 | 1.3% | |

   主导来源是**混交林/木本湿地**，不是灌木（早前把整个 tie 类标成 "shrub 9/10"
   是错的）。`argmax` 一律判给低索引，所以 PFT 1 的占比含 2.5% 的 1/7 并列、
   PFT 7 的占比含 0.8% 的 7/13 并列——这两个数当作上界看。

   **PFT 9 与 PFT 10 按构造恒等**，所以 10 在图上永远赢不了；它仍被强制留在图例里
   （`ALWAYS_IN_LEGEND`），让这个假象看得见而不是被吞进 9。`14: c4_grass`
   同理永远赢不了 PFT 13，但它不在 `PFT_MAP_COLORS` 里，所以不出现。

   配色取自 fig1 的 `AGG_COLORS`（去掉 water、shrub，再补回 water 的蓝给 PFT 10），
   颜色与 fig1 的地表类型含义**不对应**。这组色过不了 `--pairs all` 校验：
   `#1b5e20` ↔ `#c62828` 在 protan 下 ΔE **2.0**、正常视觉最差 ΔE **13.8**
   （硬门槛 ≥15），且 `#8d6e63` 过不了 chroma floor。深绿/深红那对撞这六个色
   怎么分配都躲不掉，已放在针叶林（32%）与 PFT 9（0.3%）之间。
   通过校验的一版配色留在 commit `74ef6f6`。
3. `fig3_natveg_urban_*` — PCT_NATVEG、PCT_URBAN（22/23/24 及合计）、有效像元占比。
4. `fig4_pct_nat_pft_panels_*` — 17 个 PFT 分面板。

`make_surfdata_pf` 环境**没有 cartopy**，所以图用纯 matplotlib 的经纬度坐标轴
（`set_aspect(1/cos φ)` 校正长宽比），没有州界。
