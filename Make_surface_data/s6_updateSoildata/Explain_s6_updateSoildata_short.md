# s6：用 1 km 土壤/地形数据更新 ELM surfdata（Pathfinder 简要说明）

> 本文针对 **Pathfinder（项目 hpcl-cli185）** 上的运行。完整技术细节见同目录 `Explain_s6_updateSoildata.md`。
> 主脚本：`update_Soil_surfdat.py`；HPC 包装：`update_Soil_surfdat_hpc.py`；作业脚本：`run_update_soil.slurm`。

---

## 1. 做什么

`s5_mksurfdata_map` 生成的 ELM `surfdata` 里，部分土壤/地形字段仍是约 **0.5°** 粗分辨率（图上呈大方块）。本步骤用约 **1 km** 源数据，把这些字段重采样到目标 surfdata 网格上，写出新文件，**默认不改动原文件**。

| 项目 | 内容 |
|------|------|
| 输入 | 全球 1 km NetCDF（`global_cf_float`）；北美 1 km 土壤目 GeoTIFF；一份 `surfdata_*.nc` |
| 输出 | `surfdata_Updatedsoil_<原文件名>.nc`（或 `--output` 指定路径） |
| 更新变量 | `TOPO`, `SLOPE`, `ORGANIC`, `PCT_CLAY`, `PCT_SAND`, `SOIL_ORDER`, `STD_ELEV`；可选平滑 `SOIL_COLOR` |
| 保证 | `PFTDATA_MASK == 1` 的陆地格点都有有效值；海洋格点保留原值 |

---

## 2. Pathfinder 上的路径

所有路径均在 `/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/` 下（下文记作 `$MSD`）：

| 内容 | 路径 |
|------|------|
| 工作目录（本步骤） | `$MSD/s6_updateSoildata/` |
| 1 km 全球土壤/地形 NetCDF | `$MSD/Soil_Properties/global_cf_float/`（32 GB，14 个数据文件） |
| 1 km 北美土壤目 GeoTIFF | `$MSD/Soil_Properties/Soil_Order_NA_1km.tif` |
| 输入 surfdata（s5 输出） | `$MSD/surfdata_results/surfdata_SEUS_1_24deg_simyr1850_c260712.nc` |
| 输出 surfdata | `$MSD/s6_updateSoildata/Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc` |
| 诊断图（`--plot`） | `$MSD/s6_updateSoildata/intermediate/merged_surfdat/check_*.png` |
| conda 环境 | `/projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf` |

- `Soil_Properties/` 数据复制自 Yaoping 的共享目录 `/projects/hpcl-cli185/proj-shared/ywo/Soil_Properties/data/`（见其中 `readme.txt`）。**注意：本地副本没有 `data/` 子目录**，文件直接位于 `Soil_Properties/` 下。
- `update_Soil_surfdat_hpc.py` 已把上述 `--global-cf` / `--soil-order-tif` 设为默认值，未显式传参时自动使用。

## 3. 怎么运行（slurm）

```bash
cd /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s6_updateSoildata
sbatch run_update_soil.slurm      # serial 分区，1 节点 2 核，32 GB 内存，2 小时上限
squeue -u $USER                   # 查看状态
# 日志：updateSoilSurfdat-<JOBID>.out / .err（写在提交目录）
```

作业脚本做的事：

1. `source` Pathfinder 系统 miniforge（`/software/baseline/nsp/miniforge3/24.11.3-0`），激活 `make_surfdata_pf` 环境（含 xarray / rioxarray / rasterio / scipy / netCDF4 / matplotlib）。
2. 调 `update_Soil_surfdat_hpc.py`，传入输入/输出路径及 `--smooth-soil-color --soil-color-window 9 --plot`（当前配置**开启** SOIL_COLOR 平滑，窗口 9×9；如需关闭，改回 `--no-smooth-soil-color`，或调窗口大小改 `--soil-color-window N`，N 为奇数）。

也可在计算节点上交互式直接运行（同样先激活 conda 环境）：

```bash
python update_Soil_surfdat_hpc.py \
  --surfdata ../surfdata_results/surfdata_SEUS_1_24deg_simyr1850_c260712.nc \
  --output   Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc \
  --plot
```

---

## 4. 处理流程

1. **复制**输入 surfdata → 输出文件，再在副本上修改。
2. 读 `LATIXY` / `LONGXY`，用地表掩膜 `PFTDATA_MASK` 区分陆地/海洋。
3. 按变量类型处理（只读目标区域附近的 1 km 切片，避免整幅全球数据进内存）：
   - **连续量**（高程、坡度、有机质、黏/砂粒）：双线性重采样
   - **土壤目 `SOIL_ORDER`**：众数（多数投票）
   - **`STD_ELEV`**：用 1 km 高程与子格点标准差，按总方差定律合成
   - **`SOIL_COLOR`**：对已有字段做 N×N 众数滤波（当前 slurm 配置开启，窗口 9×9；出图 `check_SOIL_COLOR_smoothing.png`）
4. 陆地缺测用 EDT 最近邻填补；只写回陆地格点。
5. 可选 `--plot` 出对比图。

---

## 5. 原理要点

| 类型 | 方法 | 原因 |
|------|------|------|
| 连续场 | GDAL 双线性 | 平滑插值到目标格点 |
| 分类场（土壤目） | mode | 类别不能取平均；水体值 0 作 nodata，不参与投票 |
| 子格点高程标准差 | \(\sqrt{\overline{\sigma_{1km}^2} + \mathrm{Var}(z_{1km})}\) | 同时保留 1 km 内变异与格点间变异，不能直接平均 `STDEV_ELEV` |
| 土壤颜色 | 众数滤波 | 消除粗分辨率上采样后的块状伪影（不另读 1 km 源） |
| 陆地填补 | EDT 距离变换 | 保证陆地无空洞，且比点云 KDTree 更快 |

内存上：先对 lat/lon 做整数切片再运算，工作量由**目标陆地区域 + 小幅 padding** 决定，与全球 1 km 全图无关（SEUS 1/24° 网格为 324×504，32 GB 内存足够）。

---

## 6. 结果

- **文件**：新的 `Updatedsoil_surfdata_*.nc`（写在 `s6_updateSoildata/` 下），上述土壤/地形变量已换成 1 km 信息；原 surfdata 不变。
- **质量**：陆地格点无 NaN；海洋保持原值；图上不再是 0.5° 大方块。
- **诊断图**（`--plot`）：`intermediate/merged_surfdat/check_<变量>.png`（左 1 km 源 / 右更新后）。

更细的 CLI 与算法说明见英文长版 `Explain_s6_updateSoildata.md`。
