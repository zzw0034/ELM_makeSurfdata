## ELM Surface Dataset Notes
### Step 1: Make ESMF

Build a true MPI-enabled ESMF on Baseline:

> ESMF_APP:
> ESMF_APP=/projects/hpcl-cli185/proj-shared/zw5/software/esmf-8.8.1-openmpi-gcc12/bin/binO/Linux.gfortran.64.openmpi.default/ESMF_RegridWeightGen

references: esmf_pathfinder_setup.md 
> /Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/doc
### Step 2: Make 17 mapping files + 2 extra mapping files

Work directory:

```bash
Mac:
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan

Pathfinder:
/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles
```

#### I prepare template:

>/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/

Reference: Bisht, Gautam (`https://github.com/bishtgautam/elm-surface-dataset`).

**steps:**
**1. Create s2_Mappingfiles folder  on Pathfinder:**
(对应/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/)
  ```bash
    mkdir /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles
  ```
**2. Copy my template to Pathfinder:**
```bash
rsync -av \
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/doc \
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF \
pathfinder:/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/
```

**3. Run create_mapping _scripts.sh script:**
  ```bash
   bash create_SEUS_mapping_scripts.sh
  ```
  > This will create a folder SEUS_1_24deg under 
  /s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/

  >  copy and edited map_*.run (1-19)

  **Note: 要手动检查修改 SEUS_1_24deg.map_domain_19.run 里的**
 ```bash
-s /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/nlcd_NA_1_24deg_grids.nc
```
> 应该是
/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/SCRPgrid_NLCD_PFTs_NA_1_24deg.nc

**3.1 生成 LUT 源网格 SCRIP（s4_LUToutput_pft/make_scrip_grid_from_nlcd.py** 
（跳到s4，完成s4，然后回来）
>/s4_LUToutput_pft/s4_Instruction.md
> /s4_LUToutput_pft/LUT_mapping_files_s2s3s4.md

map_domain_19.run need 
- Grids file that defines the grid (same grid as your LUT output)

```bash
  /s4_LUToutput_pft/scr_out/
SCRPgrid_NLCD_PFTs_NA_1_24deg.nc
  ```

   **Jump to**

  ```bash
/s4_LUToutput_pft 
 > s4_Instruction.md
 > s4_LUToutput_pft/LUT_mapping_files_s2s3s4.md
  ```

   run：

  ```bash
    make_scrip_grid_from_nlcd.py
  ```

   get SCRIP：
> /s4_LUToutput_pft/scr_out/
> SCRPgrid_NLCD_PFTs_NA_1_24deg.nc

**4. Change to the  SEUS directory:**

  ```bash
   cd SEUS_1_24deg/
   mkdir done
  ```

**5. Submit all jobs:**

  ```bash
cd /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/elm-surface-dataset-zhuonan-PF/SEUS_1_24deg
bash ../submitall.sh
  ```

**After all the job done will generate all mapping files. about 10-20min**

### List of 19 mapping grids

- **Map_01**: `SCRIPgrid_0.5x0.5_AVHRR_c110228.nc`
- **Map_02**: `SCRIPgrid_0.5x0.5_MODIS_c110228.nc`
- **Map_03**: `SCRIPgrid_3minx3min_LandScan2004_c120517.nc`
- **Map_04**: `SCRIPgrid_3minx3min_MODIS_c110915.nc`
- **Map_05**: `SCRIPgrid_3x3_USGS_c120912.nc`
- **Map_06**: `SCRIPgrid_5x5min_nomask_c110530.nc`
- **Map_07**: `SCRIPgrid_5x5min_IGBP-GSDP_c110228.nc`
- **Map_08**: `SCRIPgrid_5x5min_ISRIC-WISE_c111114.nc`
- **Map_09**: `SCRIPgrid_10x10min_nomask_c110228.nc`
- **Map_10**: `SCRIPgrid_10x10min_IGBPmergeICESatGIS_c110818.nc`
- **Map_11**: `SCRIPgrid_3minx3min_GLOBE-Gardner_c120922.nc`
- **Map_12**: `SCRIPgrid_3minx3min_GLOBE-Gardner-mergeGIS_c120922.nc`
- **Map_13**: `SCRIPgrid_0.9x1.25_GRDC_c130307.nc`
- **Map_14**: `SCRIPgrid_360x720_cruncep_c120830.nc`
- **Map_15**: `UGRID_1km-merge-10min_HYDRO1K-merge-nomask_c130402.nc`
- **Map_16**: `SCRIPgrid_0.5x0.5_GSDTG2000_c240125.nc`
- **Map_17**: `SCRIPgrid_0.1x0.1_nomask_c110712.nc`
- **Map_18**: `SCRIPgrid_0.01x0.01_nomask_c240501.nc`
- **Map_19**: `SCRPgrid_NLCD_PFTs_NA_1_24deg.nc`

---

### Additional 0.01° mapping required for `mksurfdata.pl`

`mksurfdata.pl` 需要 0.01°×0.01° 的权重文件，例如：

`map_0.01x0.01_nomask_to_small_SEUS_1_24deg_nomask_aave_da_c251030.nc`

这是用来把 **0.01° global nomask 源网格上的原始数据**
（例如：

  `/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/lnd/clm2/rawdata/mksrf_toprad_0.01x0.01.c240422.nc`

  ）插值到目标网格 `small_SEUS_1_24deg` 的 mapping。  

如果目录中只有 0.1°×0.1°（如 `map_0.1x0.1_to_..._c251030.nc`），会在 `mksurfdata.pl debug` 时报错缺少 0.01° mapping。  
Bisht 提供的 17 个 run 里没有生成 0.01x0.01，所以需要额外一步：

```text
  1. — 用 `map_0.01x0.01_18.run` 在 Baseline 上生成 0.01° mapping：  
  2. 作业完成后，在同一目录下应出现：  

     `map_0.01x0.01_nomask_to_small_SEUS_1_24deg_nomask_aave_da_c251030.nc`  

     （日期后缀与脚本中 `CDATE` 一致，需与 `mksurfdata.pl` 的 `-usr_gdate` 一致）。
```

### **Baseline 上遇到的问题与日志记录**：

#### -1 `map_0.01x0.01_nomask.run` 在 Baseline 上运行失败，`ESMF_RegridWeightGen` 报告源网格中存在退化单元（degenerate element），典型日志片段为：

```text
20260316 124702.903 ERROR            PET0 ~~~~~~~~~~~~~~~~~~~~ Degenerate Element Detected ~~~~~~~~~~~~~~~~~~~~
20260316 124702.908 ERROR            PET0   degenerate elem. id=268435457
20260316 124702.911 ERROR            PET0
20260316 124702.911 ERROR            PET0   degenerate elem. coords (lon [-180 to 180], lat [-90 to 90]) (x,y,z)
20260316 124702.911 ERROR            PET0   -----------------------------------------------------------------
20260316 124702.911 ERROR            PET0     0  (15.440000,  0.000000)  (0.963910, 0.266229, 0.000000)
20260316 124702.911 ERROR            PET0     1  (15.440000,  0.000000)  (0.963910, 0.266229, 0.000000)
20260316 124702.911 ERROR            PET0     2  (15.450000,  0.000000)  (0.963863, 0.266397, 0.000000)
20260316 124702.911 ERROR            PET0     3  (15.450000,  0.000000)  (0.963863, 0.266397, 0.000000)
20260316 124702.911 ERROR            PET0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
20260316 124702.911 ERROR            PET0 ESMCI_Mesh_Regrid_Glue.C:270 ESMCI_regrid_create() Invalid argument  - - Src contains a cell that has corners close enough that the cell collapses to a line or point
```

  日志副本保存于：  
  `PET0.RegridWeightGen.Log`

- 一个简单的权宜解决方案是在 `map_0.01x0.01_nomask.run` 中对 `ESMF_RegridWeightGen` 增加 `--ignore_degenerate` 选项，只改这一行，其余保持不动，例如：
  ```bash
  # Same number of MPI ranks as nodes (1 rank per node for good memory distribution).
  srun -n 20 ESMF_RegridWeightGen \
    --ignore_unmapped \
    --ignore_degenerate \
    -s "$SRC_SCRIP" \
    -d "$TARGET_SCRIP_PATH" \
    -m conserve \
    -w "$OUT_NAME" \
    --src_type SCRIP \
    --dst_type SCRIP \
    --64bit_offset
  ```

```bash
The 0.01° × 0.01° global SCRIP to small_SEUS_1_24deg mapping failed 
under --64bit_offset 
because the weight file exceeded NetCDF format constraints. 
After switching the output format 
to --netcdf4
the mapping completed successfully.

This has been changed in template map_0.01x0.01_18.run
```

### **Important: `*_regional` flags in `ESMF_RegridWeightGen`**

- ESMF 默认把 source 和 destination 都按 **global grid** 处理。
- 因此：**哪一边不是全球网格，就给哪一边加 `*_regional`**。
- `global -> global`：不加。
- `global -> regional`：加 `--dst_regional`。
- `regional -> global`：加 `--src_regional`。
- `regional -> regional`：同时加 `--src_regional --dst_regional`。
- 本流程里的 `small_SEUS_1_24deg` 是 **regional target**，所以所有 `... -> small_SEUS_1_24deg` 的 mapping 至少都应加 `--dst_regional`。
- 对自有 LUT 的 `nlcd_NA_1_24deg -> small_SEUS_1_24deg`，source 是北美区域、destination 也是区域，因此应加 `--src_regional --dst_regional`。
- 若遗漏这些选项，常见后果是 mapping 文件里的 `frac_a/frac_src > 1`，随后在 `mksurfdata_map` 中报 `(gridmap_map_read) ERROR: frac_src out of bounds`。

### Notes on Map 15 performance

see 

```bash
/Users/zw5/ORNL_workplace/ELM_makeSurfdata/Make_surface_data/s2_Mappingfiles/doc/map_15_issue.md
```

### Notes on sbatch 
On Baseline, this mapping case appears to require at least 8 nodes for stable execution. At 8 nodes, configurations with 88 and 96 tasks have succeeded, while 120 tasks failed, indicating that the job is sensitive not only to node count but also to the specific task layout once the minimum stable node count is reached.

---

## Mapping dataset categories (中文说明)

### 一、Land cover / PFT 相关（0.5° 或 3min）

- **Map_01 — SCRIPgrid_0.5x0.5_AVHRR**
  - 含义：AVHRR 土地覆盖产品的 0.5° 网格  
  - 来源：AVHRR（Advanced Very High Resolution Radiometer）  
  - 用途：早期 PFT/land cover 参考  
  - 分辨率：0.5°
- **Map_02 — SCRIPgrid_0.5x0.5_MODIS**
  - 含义：MODIS 土地覆盖产品的 0.5° 网格  
  - 来源：MODIS（MCD12Q1 等）  
  - 用途：Year-2000 PFT reference 常用  
  - 分辨率：0.5°
- **Map_04 — SCRIPgrid_3minx3min_MODIS**
  - 含义：更高分辨率的 MODIS 网格  
  - 用途：某些精细 PFT 处理  
  - 分辨率：3 arc-min（约 5km）
- **Map_16 — SCRIPgrid_0.5x0.5_GSDTG2000**
  - 含义：GSDTG 2000 土地覆盖网格  
  - 来源：全球 2000 年土地覆盖产品  
  - 用途：某些版本 LUT 使用的 reference  
  - 分辨率：0.5°

### 二、Urban / Population

- **Map_03 — SCRIPgrid_3minx3min_LandScan2004**
  - 含义：LandScan 2004 全球人口网格  
  - 来源：ORNL LandScan  
  - 用途：构建 urban landunit fraction  
  - 分辨率：3 arc-min

### 三、Soil 数据

- **Map_07 — SCRIPgrid_5x5min_IGBP-GSDP**
  - 含义：IGBP 全球土壤数据  
  - 用途：土壤属性（texture 等）  
  - 分辨率：5 arc-min
- **Map_08 — SCRIPgrid_5x5min_ISRIC-WISE**
  - 含义：ISRIC-WISE 全球土壤数据库  
  - 来源：ISRIC  
  - 用途：soil texture / SOC  
  - 分辨率：5 arc-min

### 四、Topography / Elevation

- **Map_11 — SCRIPgrid_3minx3min_GLOBE-Gardner**
  - 含义：GLOBE 高程数据  
  - 来源：GLOBE DEM（Gardner 处理版）  
  - 用途：地形坡度、地表高度  
  - 分辨率：3 arc-min
- **Map_12 — SCRIPgrid_3minx3min_GLOBE-Gardner-mergeGIS**
  - 含义：合并 GIS 的 GLOBE 版本  
  - 用途：改进地形数据

### 五、Hydrology / River / Drainage

- **Map_13 — SCRIPgrid_0.9x1.25_GRDC**
  - 含义：GRDC 河流流域网格  
  - 来源：Global Runoff Data Centre  
  - 用途：river routing / runoff
- **Map_15 — UGRID_1km-merge-10min_HYDRO1K**
  - 含义：HYDRO1K 全球水文数据  
  - 来源：USGS HYDRO1K  
  - 用途：流域、排水网络

### 六、Nomask 通用几何网格

- **Map_06 — 5x5min_nomask**  
- **Map_09 — 10x10min_nomask**  
- **Map_17 — 0.1x0.1_nomask**  
  - 含义：没有陆海 mask 的纯几何网格  
  - 用途：某些变量 remap 时不使用掩膜  
  - 分辨率：分别 5min / 10min / 0.1°

### 七、USGS 地形

- **Map_05 — SCRIPgrid_3x3_USGS**
  - 含义：USGS 地形网格  
  - 用途：地形 / 土壤支持数据

### 八、Climate Forcing Grid

- **Map_14 — SCRIPgrid_360x720_cruncep**
  - 含义：CRUNCEP forcing 网格  
  - 用途：气候强迫相关数据对齐  
  - 分辨率：0.5°（360x720）

