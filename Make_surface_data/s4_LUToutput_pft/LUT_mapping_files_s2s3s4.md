# LUT（NA 1/24°）→ domain（small_SEUS）mapping：哪里用、怎么用、放哪、mksurfdata 怎么找

## 1. 这个 mapping 是干什么的

- **LUT 数据**（如 `LUT_nlcd2elm_historical_1850_*.nc`）在 **北美 1/24°** 网格上（与 NLCD/elmpft 文件同一套 lat/lon）。
- **目标域**是 **SEUS_1_24deg**（东南部 US 子集，也是 1/24° 但范围更小、格点不同）。
- 需要一张 **权重文件（mapping file）**：把 LUT 源网格上每个格点的 PFT/harvest 按权重插值/聚合到目标域格点上，即 **LUT 源网格（NA 1/24°）→ small_SEUS_1_24deg** 的 mapping。

---

## 5. 怎么生成这个 mapping 文件（本目录下两步完成）

### 步骤 1：生成 LUT 源网格 SCRIP（s4_LUToutput_pft/make_scrip_grid_from_nlcd.py）

在本目录执行：

```bash
cd /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft
conda activate  make_surfdata_pf
python make_scrip_grid_from_nlcd.py
```

- 输入：NLCD/elmpft 文件
```bash
（默认 `/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc`）

```
- 输出：
> /s4_LUToutput_pft/scr_out/
> SCRPgrid_NLCD_PFTs_NA_1_24deg.nc
### 步骤 2：在 Pathfinder 上用 ESMF 生成权重并写入 MAP_DIR

map_domain_19.run need this in s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/
```bash
已经整合到
ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/map_domain_19.run

refer to
/s2_Mappingfiles/s2_Generate_Mappingfiles.md
```

---

## 2. 在哪里被用到

**在 Step 5 运行 `mksurfdata_map` 时，由程序内部使用。**

流程是：

1. 你在 namelist 里设置 `mksrf_fvegtyp`（单年 LUT）、`mksrf_fdynuse`（年度 LUT 列表）。
2. 运行 `./mksurfdata_map < namelist` 时，程序会：
   - 读入你的 LUT 文件（PFT、harvest 等），这些数据在 **NA 1/24°** 网格上；
   - 需要把这类数据 **regrid** 到 **目标网格**（small_SEUS_1_24deg）；
   - 此时会在 mapping 目录里查找 **从 LUT 源网格到目标网格的 mapping（权重）文件**；
   - 用该权重把 LUT 从 NA 1/24° 插值/聚合到 small_SEUS 格点，再写入 `surfdata_*.nc` 和 `landuse.timeseries_*.nc`。

也就是说：**mapping 不是在命令行里显式指定的，而是 mksurfdata_map 读 LUT 并做 regrid 时，按“源网格名 → 目标网格名”的约定在 mapping 目录里查找并使用的。**

---

## 3. 文件放在哪里

- 与 Step 2 那 **17 个** mapping 文件放在 **同一目录**，即 **MAP_DIR**。
- MAP_DIR 是你运行 `mksurfdata.pl` 时用 **`-usr_mapdir`** 指定的目录；mksurfdata.pl 会把它写进 namelist 里每条 `map_*` 的路径前缀，所以 mksurfdata_map 读 namelist 时得到的已经是“每个 mapping 文件的完整路径”，这些路径都指向该目录下的文件。

**当前流程中的 MAP_DIR（Pathfinder）：**

```text
/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/SEUS_1_24deg
```

把 **LUT 源网格 → small_SEUS_1_24deg** 的 mapping 文件生成后，放到这个目录即可。

---

## 4. mksurfdata 怎么找到这个路径

- **mapping 目录（MAP_DIR）** 不是由 mksurfdata_map 单独读一个“MAP_DIR 变量”得到的，而是：你运行 **mksurfdata.pl** 时传入了 **`-usr_mapdir <路径>`**，mksurfdata.pl 在生成 namelist 时，把所有 `map_*` 变量都写成「该路径 + 文件名」，例如  
  `map_fpft = '/gpfs/.../small_SEUS_1_24deg/map_0.5x0.5_MODIS_to_small_SEUS_1_24deg_...'`  
  因此 mksurfdata_map 读 namelist 时，拿到的就是每个 mapping 的**完整路径**。
- **LUT 用的 mapping**：当 mksurfdata_map 读入 LUT 并要做 regrid 时，会根据 LUT 的源网格名（或 E3SM 内部约定）在**同一目录**（即上述 map_* 所在目录）里查找形如  
  `map_<LUT源网格名>_to_small_SEUS_1_24deg_*.nc`  
  的文件。所以只要把「LUT 源网格 → small_SEUS_1_24deg」的 mapping 文件放在 MAP_DIR 下，且文件名符合 E3SM/ELM 的命名约定，程序就能找到并用来 regrid LUT。

---



若 mksurfdata 查找 LUT 源网格时使用的名字与 `nlcd_NA_1_24deg` 不同，可能需要在 MAP_DIR 内做符号链接或重命名以匹配 E3SM 约定。
