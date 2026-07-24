# 我的 ELM make surface data workflow 指南

## Top-level layout

- `Make_surface_data/`: my_workflow_Make_surface_data.md, organized by processing step s1-s8.
- `Docs/`: reference PDFs, Word documents
- `pyproject.toml`, `uv.lock`, `.venv/`: local Python environment in Make_surface_data.

---
**THIS whole workflow on Pathfinder.**

**/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data**

---

## `Make_surface_data/` layout

- `s1_Gridfiles/`: build and inspect grid / SCRIP inputs.
- `s2_Mappingfiles/`: mapping-generation scripts and related notes.
- `s3_Domainfiles/`: domain-file generation and resolution checks.
- `s4_LUToutput_pft/`: NLCD-to-PFT, LUT generation, LUH2 harvest downscaling,
and related references.
- `s5_mksurfdata_map/`: notes and assets for `mksurfdata_map`.
- `s6_updateSoildata/`: update soil order, soc, soil texture.
- `s7_updataPdata/`: update soil P.
- `s8_bilinearLargeSquare_hdm/`: update hdm population density.

---

## 1. High-level workflow (from LAND summary)

The ELM surface dataset pipeline has these stages in order:

| Step                     | Description                                                                                                 |
| ------------------------ | ----------------------------------------------------------------------------------------------------------- |
| **1. Grid files**        | SCRIP-format grid defining the target simulation grid (e.g. `small_SEUS_1_24deg_grids.nc`)                  |
| **2. Mapping files**     | Maps from source data grids (e.g. 0.5°, 3×3min, …) to your target grid; used by mksurfdata to regrid inputs |
| **3. Domain files**      | Land domain for the target grid (e.g. `domain.lnd.small_SEUS_1_24deg.nc`) with landfrac/mask                |
| **4. PFT/CFT data**      | Plant functional type (and optionally crop) fractions. Either standard LUH2→LUT or custom (e.g. NLCD-based) |
| **5. ELM surface files** | Final `surfdata_*.nc` produced by **mksurfdata** (optionally with topounits)                                |

---

### working directory:

```bash
Mac:
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data

Pathfinder:
/gpfs/wolf2/cades/cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data
```

### enviroment

uv on Mac
conda on Pathfinder

#### conda environments:

conda activate make_surfdata_pf

### Step 1 — Grid files

- **脚本**：`s1_Gridfiles/make_scrip_1_24_deg.py`
- **Pathfinder**: sbatch run_make_scrip_1_24_deg.slurm
- **输出目录**：`s1_Gridfiles/scr_out`
- **产物**：SCRIP 格式网格，例如 `SCRIPgrid_SEUS_1_24deg.nc`

---

### Step 2 — Mapping files（Make 17 mapping files + 2 extra mapping files）

- **内容**：
Make 17 mapping files + 2 extra mapping files
- **用法**：
  > s2_Mappingfiles/s2_Generate_Mappingfiles.md

---

### Step 3 — Domain files (需要先完成S4)

从 PFT/植被 NetCDF 与 SCRIP 目标网格生成 ELM 陆面 domain (landfrac + mask).
- **脚本目录**：`/ELM_makeSurfdata/Make_surface_data/s3_Domainfiles/`
- **主要脚本**：`make_land_domain_from_pft.py`（读 SCRIP + PFT/surfdata，写出 domain）
- **产出**：例如 `domain.lnd.SEUS_1_24deg.nc`
- **note**: use vegetation pft data to quantify landfrac 0/1

---

### Step 4 — PFT only（No crop），harvest downscaling from LUH2

```bash
s4_Instruction: 

/s4_LUToutput_pft/

s4_Instruction.md
LUT_mapping_files_s2s3s4.md
```

- **4a. NLCD → ELM PFT**：
draft coding in  
`/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/NLCD2ELM` 
`/Users/zw5/ORNL_workplace/ELM_makeSurfdata/NLCD2ELM`
完成（Dan 的 NLCD 转 ELM PFT）。

Standard code: 
>/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/s4_1_convert_nlcd_to_PFT.py


- **4b. Downscale harvest data from LUH2**
**做成 LUT output 格式**：  
  - `/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/s4_2_donwscale_LUH2harvest.py`  
  - **脚本**：`s4_2_donwscale_LUH2harvest.py` 生成 **LUT 格式的年度文件**。
  - 这些 LUT output 才是 mksurfdata 直接读的输入（单年文件给 `mksrf_fvegtyp` / `-dynpft`，多年则列成列表给 `mksrf_fdynuse`）。
- **给 mksurfdata 的输入**：即 LUToutput 目录产出的 LUT 格式文件（单年或年度列表）。
- **download LUH2 data**：
- /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/LUH2
- /projects/hpcl-cli185/proj-shared/zw5

- **Recent period 2015–2023（harvest）**：LUH2 v2h `transitions.nc` 的采伐止于日历年 2014。
  近期不再冻结在 2014，而是用 LUH2 **v2f SSP2-4.5** 续接：
  `/projects/hpcl-cli185/proj-shared/zw5/luh/ssp2rcp45_transitions.nc`
  （0.25°，2015–2099，网格 / lat 顺序 / 采伐变量与 v2h 一致）。
  按 LUT 约定（标签年 Y 装日历年 Y−1 的采伐）：**输出 2016–2023 取 SSP2-4.5 日历 2015–2022，其余仍取 v2h**。
  **原因**：SSP 情景 2030 前几乎不分岔，2015–2022 各情景采伐近乎重合，拿哪个都差不多；SSP2-4.5 为"近现实"中间路径。

- **s4_2 采伐降尺度已修正（2026-07-23）**：① 面积守恒（`frac_i = F·w_i·ΣA/Σ(w·A)`，旧版守恒的是分数之和，低估约 36×，已验证 `LUH2可分配/ours` 35.99→1.00）；② `PCT_NATVEG` 除数/权重改逐年（LANDMASK 仍用 1850 静态）；③ 时间标签 k=−1（标签年 Y 装日历年 Y−1）。`PCT_PFT` 不受影响（174 年逐位不变）。旧文件存于 `output_downscaled_luh2_harvest/prev_buggy_c0711/`。
  > **状态**：上面 SSP2-4.5 续接为已决定、**待实现**；当前 s4_2 仍冻结在 v2h 2014。
---

## Step 5 — 上手跑 mksurfdata

mksurfdata 位于：`E3SM/components/elm/tools/mksurfdata_map`。  

- 本地：`/Users/zw5/ORNL_workplace/models/E3SM/components/elm/tools/mksurfdata_map`  
- Pathfinder：`/projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map`

### 5.0 Set up the input files（mksrf_fdynuse 列表）— 已完成

在跑 mksurfdata 前，需要把**年度 LUT 文件**列在一个文本文件里，namelist 的 `mksrf_fdynuse` 指向该文件。格式：每行一个完整路径，路径占前 195 字符（不足用空格补齐），**第 197 列起为 4 位年份**。

> /s4_LUToutput_pft/make_fdynuse_list.py

- **已生成列表文件**：
  - 本地 `Make_surface_data/s4_LUToutput_pft/LUT_nlcd2elm_luh2_historical_list.txt`；
  - Pathfinder 上为 `/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest/LUT_nlcd2elm_luh2_historical_list.txt`。  
  - 年份范围：1850–2023（174 行）  
  - 每行格式：`<完整路径>/LUT_nlcd2elm_historical_YYYY_03032026.nc`（补足 195 字符）+ 空格 + 年份  

之后在 5.6 里把 namelist 的 **mksrf_fdynuse** 设为该 `.txt` 的**完整路径**即可。

### 5.1 确认手头已有

- **Step 5.0**：mksrf_fdynuse 列表文件已生成（`s4_LUToutput_pft/LUT_nlcd2elm_luh2_historical_list.txt`，1850–2023）；脚本 `make_fdynuse_list.py` 可复用于换路径/年份
- **Step 5.2、5.3**：5.2 加载环境（含 USER_FFLAGS 与 E3SM 两处源码修改）、5.3 编译 mksurfdata_map 已完成（详见《Load the relevant modules and configure the environment》）

```bash
'/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s5_mksurfdata_map/Load the relevant modules and configure the environment.md'
```

- Step 1：`scr_out/<你的 SCRIP>.nc`（如 `small_SEUS_1_24deg_grids.nc`）
- Step 2：17 个 mapping 文件在同一目录（MAP_DIR）；**若用自有 NLCD LUT**，MAP_DIR 内还须有 LUT 源网格→small_SEUS 的 mapping（`map_nlcd_NA_1_24deg_to_small_SEUS_1_24deg_*.nc`），见 Step 2 补充与 `s4_LUToutput_pft/README_LUT_mapping_usage.md`。
- Step 3：`domain.lnd.<网格名>.nc`
- Step 4：LUToutput / s4_LUToutput_pft 产出的 LUT 格式文件（单年 + 年度列表）

### 5.2 加载环境 — 已完成

在 Pathfinder 上按以下步骤操作（详见 **s5_mksurfdata_map/Load_modules_configure_environment.md**）：

>/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s5_mksurfdata_map/Load_modules_configure_environment.md


### 5.3 编译 mksurfdata_map — 编译已通过

在完成 5.2 后，在 E3SM 的 `mksurfdata_map/src` 下执行：

```bash
cd /projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map/src
gmake clean
gmake -j 8
```

```text
之前baseline的proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map
I cp it to proj-shared/zw5/software 作为backup。
里面有之前生成surfdata
可以删除。
```

**当前状态**：编译已通过，可执行文件在 `mksurfdata_map/mksurfdata_map`（即 Pathfinder：`/projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map/mksurfdata_map`）。

### 5.4 生成 namelist（先 debug 看会生成什么）— 已跑通。

- **mksurfdata.pl** 根据参数写出 `**namelist`**（默认文件名 `namelist`），供 `**mksurfdata_map`** 读取。
- 常见情况是：**无论是否 `-d`，都会写 namelist 并调用 `mksurfdata_map`**；`-d` 多用来**多打印/多检查**，便于发现 mapping 缺失等。**若你已在 `-d` 下看到 “Successfully created fsurdat files”，surfdata 通常已经生成，不必仅为「正式」再去掉 `-d` 重跑一遍**（除非你想换输出名或确认无 `-d` 时行为一致）。
- 改 LUT、改列表等：可**编辑 `namelist`** 后执行 `./mksurfdata_map < namelist`，不必重复改 `mksurfdata.pl` 的长命令。

**用户自定义网格（usrspec）**。下面为当前使用的值（pathfinder），可直接复制运行：

```bash
cd /projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map

./mksurfdata.pl -res usrspec \
  -usr_gname SEUS_1_24deg \
  -usr_gdate 260710 \
  -y 1850 \
  -dinlc /projects/hpcl-cli185/world-shared/e3sm/inputdata \
  -usr_mapdir /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/SEUS_1_24deg \
  -d

===========================================
Successfully created fsurdat files
```

**当前状态**：已用上述命令在 Pathfinder 上跑通：namelist 已生成，`mksurfdata_map < namelist` 已成功执行并生成 surfdata 文件（如 `surfdata_small_SEUS_1_24deg_simyr1850_c260304.nc`、`.log`）。若要用**自己的 NLCD LUT** 做 PFT/transient，需按 5.5、5.6 编辑 namelist 中的 `mksrf_fvegtyp`、`mksrf_fdynuse` 等，再重新执行 mksurfdata.pl（可加 `-dynpft` 或直接改 namelist 后 `./mksurfdata_map < namelist`）。


### 5.5 用你自己的 PFT（NLCD→LUT）替代默认

若要用 s4_LUToutput_pft / LUToutput 产出的 LUT 作为 PFT 来源：**直接在 5.6 编辑 namelist** 填写 `mksrf_fvegtyp`、`mksrf_fdynuse` 等 5 项即可。

### 5.6 编辑 namelist 里 5 + 2 个关键量

按你的 workflow（自定义 LUT、small_SEUS 1/24°、1850–2023）逐项说明，并在最后给出可直接参考的 namelist 片段。

**1) mksrf_fvegtyp**

- **含义**：生成单一年份 surfdata 时使用的 PFT 源文件。
- **在本流程中**：即某一年的 LUT 输出，通常选 1850。
- **示例**：`mksrf_fvegtyp = '/path/LUT_SE/LUT_SE_1850.nc'`
- **要点**：该文件必须与 fsurdat 的年份一致（例如 simyr1850）；只用于构建单年 surfdata。

**2) mksrf_fdynuse**

- **含义**：列出所有年份 LUT 文件的文本列表，用于生成动态 landuse 时间序列。
- **在本流程中**：即 5.0 里生成的列表文件（如 `LUT_nlcd2elm_historical_list.txt`，1850–2023）。
- **示例**：`mksrf_fdynuse = '/path/fdynuse_1850_2023.txt'`
- **要点**：程序按列表顺序读取 LUT 文件，生成 time dimension。

**3) fdyndat**

- **含义**：输出的 landuse 时间序列文件名（多年的动态 land cover）。
- **推荐格式**：`landuse.timeseries_<res>_<scenario>_simyr<yyyy>-<yyyy>_c<yymmdd>.nc`
  - `<res>`：你的网格名，如 `small_SEUS_1_24deg`
  - `<scenario>`：自定义标签，如 `nlcd2elm` 或 `custom`
  - `<yyyy-yyyy>`：时间范围
- **示例**：`fdyndat = 'landuse.timeseries_small_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260305.nc'`

**4) fsurdat**

- **含义**：输出的单年 surfdata 文件名，包含 soil、slope、glacier、wetlands、landunit 以及来自 mksrf_fvegtyp 的 PFT。
- **推荐格式**：`surfdata_<res>_<scenario>_simyr<yyyy>_c<yymmdd>.nc`
- **示例**：`fsurdat = 'surfdata_small_SEUS_1_24deg_nlcd2elm_simyr1850_c260305.nc'`
- **要点**：simyr1850 必须与 mksrf_fvegtyp 对应的年份一致。

**5) fsurlog**

- **含义**：surfdata 生成过程的 log 文件。
- **示例**：`fsurlog = 'surfdata_small_SEUS_1_24deg_nlcd2elm_simyr1850_c260305.log'`

**完整 namelist 片段示例**（在 `namelist` 中找到 `&elmexp` 或 `&clmexp`，修改或确认以下 5 项）：

```text
&elmexp
 mksrf_fvegtyp  = '/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest/LUT_nlcd2elm_luh2_historical_1850_07112026.nc'
 mksrf_fdynuse  = '/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/output_downscaled_luh2_harvest/LUT_nlcd2elm_luh2_historical_list.txt'
 fdyndat        = 'landuse.timeseries_small_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260712.nc'
 fsurdat        = 'surfdata_SEUS_1_24deg_simyr1850_c260712.nc'
 fsurlog        = 'surfdata_SEUS_1_24deg_simyr1850_c260712.log'
 ...

  map_fpft          = '/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/SEUS_1_24deg/map_nlcd_NA_1_24deg_to_SEUS_1_24deg_nomask_aave_da_c260710.nc'
  map_fharvest      = '/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/SEUS_1_24deg/map_nlcd_NA_1_24deg_to_SEUS_1_24deg_nomask_aave_da_c260710.nc'
 
/
```

以上为 Pathfinder 上当前使用的路径；日期后缀 `c260409` 可按需要改为当前日期（yymmdd）。

### 5.6A 运行 mksurfdata_map 时如何更改原始数据（以 soil texture 为例）

**NOTE: 目前我才用s6-s8，用python script post-process surfdata。**
```test
用## Step 6 update soil data
生成surfdata之后，update surfdata里的soil data

```
可以。`mksurfdata_map` 支持替换原始输入数据。  
对 soil texture 来说，最关键的是 namelist 里的两项：

- `mksrf_fsoitex`：soil texture 原始数据文件路径（你要替换成自己的）
- `map_fsoitex`：该原始数据源网格 -> 目标网格的 mapping 文件

**推荐做法（最稳妥）**：

1. 先用 `mksurfdata.pl -d` 生成 `namelist`（先检查配置、文件名与路径）。
2. 在 `namelist` 中手动修改：
  - `mksrf_fsoitex = '/path/to/your/new_soil_texture.nc'`
  - `map_fsoitex   = '/path/to/map_newsoil_to_small_SEUS_1_24deg_...nc'`
3. 运行：
  ```bash
   ./mksurfdata_map < namelist
  ```
4. 在输出日志中确认程序确实读了新文件（会打印 `soil texture from` 和 `mapping for soil texture`）。

**注意**：

- 若新 soil texture 与原始数据使用**同一源网格**，通常只改 `mksrf_fsoitex` 即可。
- 若新 soil texture 的**源网格改变**，必须重做并替换 `map_fsoitex`，否则可能报错或结果不正确。
- `mksurfdata.pl` 还有常数覆盖选项（不是换数据集）：`-soil_cly XX -soil_snd YY`，用于全域统一 clay/sand 百分比。

### 5.7 运行 mksurfdata_map

```bash
在 login 节点跑会被 Killed
所以在 Pathfinder 上启动 Interactive Batch Job

interactive -N 1 -n 8 -c 1 --mem=128gb --time=1-00:00:00


module purge
module load gcc/12.4.0
module load openmpi/5.0.5
module load hdf5/1.14.5-mpi
module load netcdf-c/4.9.2-mpi-h5f
module load netcdf-fortran/4.6.1-mpi-h5f

cd /projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map
./mksurfdata_map < namelist

结束时直接输入：
exit
```

### Output

```text
 Surface data output file = surfdata_SEUS_1_24deg_simyr1850_c260712.nc
    This file contains the land model surface data
 Diagnostic log file      = surfdata_SEUS_1_24deg_simyr1850_c260712.log
    See this file for a summary of the dataset

 Successfully created surface dataset
```
### copy sufdata.nc to Make_surface_data

1. > cd /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data
   > mkdir surfdata_results
   > mkdir surfdata_results
2. > cp /projects/hpcl-cli185/proj-shared/zw5/E3SM/components/elm/tools/mksurfdata_map/{surfdata_SEUS_1_24deg_simyr1850_c260712.nc,landuse.timeseries_SEUS_1_24deg_nlcd2elm_simyr1850-2023_c260712.nc} /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/surfdata_results/


## Step 6 update soil data

**s6 this step to update the soil data in surfdata.**
See

```bash
s6_updateSoildata/

Explain_s6_updateSoildata_short.md
Explain_s6_updateSoildata.md
```
> /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata
> Explain_s6_updateSoildata_short.md

> sbatch run_update_soil.slurm


## Step 7 update soil P data

**Update soil P data.**

1. generate P data with Xiaojuan method.
2. updata surfdata.nc

```bash
s7_updataPdata
Explain_s7_updataPdata.md
Explain_s7_updataPdata-PF.md
```

## Step 8 - s8_update_hdm

**Update hdm population density (0.5) data to avoid large square.**

原因：粗分辨率人口密度（HDM）在0.5°网格边界上的真实台阶，经 Li 火灾模型 + spin-up 死亡率放大，触发火-植被双稳态，造成 GPP 锐利分界。

> /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s8_update_hdm

> sbatch run_update_hdm.slurm
---
**!!!生成最终surfdata.nc!!!**
---
