# s7_updataPdata — 更新 surfdata 中的磷 (P) 数据

> **快速上手：** 下次做 s7，只看下面「快速上手」这一节就够了。
> 详细原理、单位推导、脚本逐行说明、历史决策都放在后面的
> 「背景与原理（详细）」部分。

---

## 快速上手（下次直接看这里）

### 为什么要做 s7（一句话）

ELM 原始 surfdata 里的**磷 (P) 数据来自粗分辨率输入（1° / 0.5°）**，
在北佛罗里达、Florida Panhandle、Chesapeake Bay 等区域 P 值被设成了
**0**，导致模拟结果出现方块状的空间伪影。s7 用 Dr. Xiaojuan Yang 的
原始岩石 P 数据 + Yang et al. (2013) 方法，重新生成高质量的 P 场，
写回 s6 输出的 surfdata，消除这些伪影。

s7 紧接在 `s6_updateSoildata`（更新土壤序 SOIL_ORDER）之后。

### 整个流程分两步

| 步骤 | 脚本 | slurm | 作用 |
|---|---|---|---|
| **01** | `01_Generate_P_data.py` | `run_generate_P.slurm` | 由岩石 P + SOIL_ORDER 生成 TPS 和 5 个 Hedley P 分量（`outputs/P_forms_*.nc`） |
| **02** | `02_Update_Surfdata_Pdata_hpc.py` | `run_update_P.slurm` | 把其中 4 个矿质 P 分量写回 surfdata，得到最终 `surfdata_UpdatedsoilP_*.nc` |

**必须先跑 01，再跑 02**（02 读的是 01 生成的 `outputs/P_forms_SEUS_1_24deg.nc`）。

### 需要什么输入

| 输入文件 | 用途 | 被谁读 |
|---|---|---|
| `/projects/hpcl-cli185/proj-shared/xyk/P_rock_V3/P_rock_den.nc` | 母质岩石 P 面密度（kg P/m², 全球 0.5°） | 01 |
| `../s6_updateSoildata/Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc` | s6 输出的 surfdata，提供 `SOIL_ORDER` + `LATIXY`/`LONGXY`（SEUS 1/24° 网格）；也是 02 要写回的底稿 | 01 和 02 |
| `outputs/P_forms_SEUS_1_24deg.nc` | 01 的产物，02 的输入 | 02 |

> 如果换了新的 s6 输出文件（日期戳变化），需要同步修改两个脚本里的
> 文件名：`01_Generate_P_data.py` 里的 `surfdata_path`，以及
> `02_Update_Surfdata_Pdata_hpc.py` 里的 `SURF_IN`。

### 生成什么输出

**01 生成：**
- `outputs/TPS_SEUS_1_24deg.nc` — 顶部 50 cm 土壤总磷 TPS（g P m⁻²）
- `outputs/P_forms_SEUS_1_24deg.nc` — TPS + 5 个 Hedley 分量（**02 的输入**）
- `intermediate/` — `ppdi_eps_stage_*.nc`、若干 GeoTIFF、3 张 sanity PNG（供 QGIS / 目视检查）

**02 生成：**
- `surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c<今天>.nc` — **最终 surfdata**（这是给 ELM 用的）
- `intermediate/Pform_*_{old,new}.tif`（8 个 GeoTIFF）+ `check_surfdata_P_update.png`

被改写的 4 个变量：`LABILE_P`、`SECONDARY_P`、`APATITE_P`、`OCCLUDED_P`
（单位 g P m⁻²，只改 `PFTDATA_MASK==1` 的陆地格点；其余变量逐比特保留）。

### 怎么在 Pathfinder 上运行

脚本和 slurm 都在
`/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s7_updataPdata/`，
`WORKDIR` 已指向该目录，直接 `sbatch` 即可：

```bash
cd /projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s7_updataPdata

# 第一步：生成 P 场
sbatch run_generate_P.slurm
# 等它跑完，确认 outputs/P_forms_SEUS_1_24deg.nc 已生成

# 第二步：写回 surfdata
sbatch run_update_P.slurm
```

两个 slurm 都：
- 复用 s6 建好的 conda 环境
  `/projects/hpcl-cli185/proj-shared/zw5/conda_envs/make_surfdata_pf`
  （已含 xarray / scipy / netcdf4 / matplotlib / rasterio）；
- `export MPLBACKEND=Agg`（计算节点无显示，强制 Agg 后端让诊断 PNG 正常写出）；
- 单节点串行作业（`-N 1 -n 1 -c 2`，`--mem=32gb`，30 分钟墙钟上限）。

日志：`generatePdata-<JOBID>.out/.err` 和 `updatePsurfdat-<JOBID>.out/.err`。
跑完检查 `intermediate/` 下的 sanity PNG 确认无异常，再把最终
`surfdata_UpdatedsoilP_*.nc` 交给下游 ELM。

---

# 背景与原理（详细）

## 为什么创建这个文件夹

跑完 ELM 模拟后，结果出现**异常空间格局**——大块方形伪影，位于：

- 北佛罗里达
- Florida Panhandle
- Chesapeake Bay 区域

### 根本原因

排查后确认来源：**surfdata 里的磷 (P) 数据来自粗分辨率输入**
（1° 或 0.5°）。在受影响区域 P 值被设为**零**，进而在模拟输出中产生
上述伪影。

参考诊断：
- `checkSurfdata/` — surfdata QC 检查
- `/Users/zw5/ORNL_workplace/ELM_output_read/outputs/SE_hires_20TRC_pathfinder` — 显示伪影的模拟输出

本步骤紧接 `s6_updateSoildata`，用更高质量的源数据填补 P 数据缺口。

---

## 与 Dr. Xiaojuan Yang 的讨论

P 数据问题与 **Dr. Xiaojuan Yang** 讨论过，她：

1. **建议用她自己创建的原始源数据来更新 P**。
2. **方法已发表在同行评审论文**：
   > Yang, X. et al. (2013). *The distribution of soil phosphorus for global biogeochemical modeling.* Biogeosciences, 10, 2525–2537.
   > [https://bg.copernicus.org/articles/10/2525/2013/bg-10-2525-2013.pdf](https://bg.copernicus.org/articles/10/2525/2013/bg-10-2525-2013.pdf#page=4.93)
   > *（论文给出了构建 P 数据的分步方法。）*
3. **原始数据位置**：
   ```
   /projects/hpcl-cli185/proj-shared/xyk/P_rock_V3/
   ```

参考链接：
- https://bg.copernicus.org/articles/10/2525/2013/
- https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-global-phosphorus-dist-map-1223-1#toc-product-summary

---

## `P_rock_V3` 里有什么

> **一句话总结：** `P_rock_V3` 产出一张**全球 0.5° 母质 P 浓度图**
> (`P_rock.bin` / `rockP.nc`)，以及配套的 0.5 m 母质层 P 面密度
> (`P_rock_den.*`)，对应论文 **Sect. 3.1 / Fig. 1 Step 1**。

### 最终输出文件（论文 Step 1）

| 文件 | 说明 | 单位 / 备注 |
|---|---|---|
| `P_rock.bin` | 母质/基岩 P 浓度 (C_P) | ppm (= mg/kg) |
| `rockP.nc` | `P_rock.bin` 的 NetCDF 版 | 0.5°×0.5°, 720×360, 变量 `rockP` (mg/kg) |
| `P_rock_den.bin` | 0.5 m 母质层的 P 储量 | kg P/m²（**不是** 50 cm 土壤 TPS） |
| `P_rock_den.nc` | `P_rock_den.bin` 的 NetCDF 版 | 同上 |
| `P_rock.eps` | `P_rock.bin` 的全球图 | 仅图 |
| `P_rock_den.eps` | `P_rock_den.bin` 的全球图 | 仅图 |

**重要澄清：**
- 这**不是**顶部 50 cm 的土壤总磷 (TPS)。
- 这**不是**各 P 分量（labile / occluded / organic 等）的拆分。
- 它是 0.5° 全球分辨率的**母质（基岩来源）P 浓度**——推导生态系统可用 P 的第一步。

### 两个 NetCDF 输出的单位

**`rockP.nc` — 变量 `rockP`**
文件里标注正确：`rockP:units = "mg/kg"`（即质量 ppm）。这是 $C_P$，**母质岩石 P 浓度**。
文件内经验范围：~264 – 2288 mg kg⁻¹（均值 ≈ 646），与 `P_rock.f90` 的岩性分类表一致。

**`P_rock_den.nc` — 变量 `rockP_den`**
⚠️ **文件里的元数据是错的。** NetCDF 头写着 `rockP_den:units = "ppm"`，
但那是 `rockP_den.ncl`（第 34 行）从 `rockP_ncout.ncl` 克隆时忘了改单位串的复制粘贴 bug。

**真实单位是 kg P m⁻²**（0.5 m 母质层每平方米的 P 储量）。可直接从 `P_rock.f90` 读出：

```fortran
where(p_rock > 0 .and.bd_rock >0.)
     P_rock_den = p_rock*1.0e-6*bd_rock*1.0e3*0.5
elsewhere
     p_rock_den = -9999.
endwhere
```

**逐项量纲检查：**
- $p\_rock$ [mg P / kg rock] $\times 1.0\text{e}{-6}$ → kg P / kg rock
- $bd\_rock$ [g cm$^{-3}$] $\times 1.0\text{e}3$ → kg m$^{-3}$（岩石体密度）
- $0.5$ → m（母质层深度 0.5 m）
- **乘积 = kg P m⁻²**

文件内经验范围：0.229 – 2.974（均值 ≈ 0.78）。作为 kg P m⁻² 物理上合理，
作为 ppm 则不合理——进一步证明 `units = "ppm"` 是过时标签。

**换算到 Yang et al. (2013) Eq. 2 用的 g P m⁻²：**
$$
\text{g P m}^{-2} = \text{rockP\_den (kg P m}^{-2}) \times 1000
$$

---

## Fig. 1 — 生成土壤 P 图的三步（Yang et al., 2013）

完整土壤 P 数据集分三步构建（引自论文 Sect. 3.3）：

**Step 1 — 母质 P 浓度图**
把全球**母质（岩性）图**与**岩石 P 浓度数据库**结合，生成母质 P 浓度 C_P 的网格图。
→ 这就是 `P_rock_V3` 产出的 `P_rock.bin` / `rockP.nc`（0.5°，单位 ppm = mg/kg）。

**Step 2 — 顶部 50 cm 土壤总磷 (TPS)**
把土壤序图叠加到母质 P 图上，应用 **Eq. 2**：

$$
\mathrm{TPS} = 0.01\,D\;\frac{\rho_{P}\,C_{P}\,\left(1-\mathrm{PPDI}\right)}{\varepsilon + 1} \tag{2}
$$

即 `TPS = 0.01 · D · ρ_P · C_P · (1 − PPDI) / (ε + 1)`

| 符号 | 含义 | 单位 |
|---|---|---|
| TPS | 顶部 50 cm 土壤总 P | g P m⁻² |
| D | 土壤深度 | 50 cm |
| ρp | 母质体密度 | g cm⁻³ |
| CP | 母质 P 浓度 | ppm |
| PPDI | 磷成土耗损指数（风化过程累积 P 损失） | — |
| ε | 顶部 50 cm 体积应变 | — |

PPDI 按风化阶段（Table 2）：
- 轻度风化 (Slight)：~12%
- 中度风化 (Intermediate)：~54%
- 强风化 (Strong)：Spodosols 65%、Ultisols 70%、Oxisols 90%

**Step 3 — 各 P 分量图**
用 Hedley 分级比例（来自 178 条文献测量、按土壤序分组）把 TPS 拆成不同 P 池：
- Labile P（resin-P + NaHCO₃-P）
- Moderately labile / Secondary P（NaOH-P）
- Occluded / mineral P（HCl-P）+ Apatite P
- Organic P

> **对本工作的结论：** `P_rock_V3` 只覆盖 **Step 1**。Step 2、3 由本文件夹的
> `01_Generate_P_data.py` 完成，得到 ELM 需要的 TPS 和 P 分量。

---

## 决策日志（历史）

**Step A — 关于 1 km 岩性数据（2026-05-19 解决）**
询问 Xiaojuan 能否找到原始 **Dürr et al. (2005)** 岩性矢量数据及其派生的 1 km 栅格。
**→ Xiaojuan 确认找不到。** 因此放弃 1 km 岩性方案，改用 `P_rock_V3` 的 0.5° `rockP.nc` 作为 C_P 源。

**Step B — 用 `P_rock_V3` (0.5°) 作为 C_P 输入（2026-05-28）**
- 源：`/projects/hpcl-cli185/proj-shared/xyk/P_rock_V3/P_rock_den.nc`（C_P，kg P/m²，0.5°）
- 与原论文使用的 C_P 一致；如需 g P m⁻²：`× 1000`。

**Step C — 用 Eq. 2 计算 TPS**
- `rockP_den_resampled` 已经编码了 `0.01 · D · ρ_P · C_P`（g P m⁻²），见 `01_Generate_P_data.py` 的 TPS 计算段。
- 用 SEUS 1/24° 的 `SOIL_ORDER` 推 PPDI（Table 2）和应变 ε（按风化类别）。

**Step D — 推导各 P 分量（Step 3）**
按土壤序对 TPS 应用 Hedley 分级比例，得到 labile / secondary / apatite / occluded / organic 各分量图。

关于 **Histosol**（沿海有分布）：Yang 2013 Table S2 不含 Histosol / Andisol / Gelisol。
Histosol 的 P 分量从 Yang & Post (2011, Biogeosciences 8, 2907–2916) Fig. 3 按柱高读取：

| P fraction | 占总土壤 P (%) |
|---|---|
| Labile Pi | 10.1% |
| Secondary Pi | 5.8% |
| Apatite P | 8.3% |
| Occluded P | 40.8% |
| Organic P | 35.0% |

（总和 ≈ 100%。）ELM soilorder 2:Gelisol、1:Andisol 不在 Table S2，且不在 SE US，故不考虑。

**Step E — 回写 surfdata**
把 1/24° 的 P 分量图写回 ELM surfdata，替换掉粗分辨率/为零的 P 值。→ 由 `02_Update_Surfdata_Pdata_hpc.py` 完成。

---

## `01_Generate_P_data.py` 详解

端到端流程：把 0.5° 母质 P 面密度 (`P_rock_den.nc`) 和 1/24° SEUS `SOIL_ORDER` 栅格，
按 Yang et al. (2013) Eq. 2 转成逐像元的 TPS 和 5 个 Hedley P 池图。

### 输入
- `/projects/hpcl-cli185/proj-shared/xyk/P_rock_V3/P_rock_den.nc` — 全球 0.5° 母质 P 面密度（kg P/m²）
- `../s6_updateSoildata/Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc` — s6 输出，提供 1/24° `SOIL_ORDER` 栅格及 `LATIXY`/`LONGXY`

> 这两条路径在脚本里是**硬编码的绝对/相对路径**；换输入时改脚本顶部对应行。

### 参考字典（"物理"查找表）
- `SOIL_ORDER_NAMES` / `SOIL_ORDER_CODES` — 整数 ↔ USDA 土壤序名（代码 0–15）
- `SOIL_ORDER_NONLAND_CODES = {0, 13, 14, 15}` — Water / shiftingsand / rockland / iceglacier（排除出陆地掩膜）
- `SOIL_ORDERS_BY_STAGE` — 9 个风化阶段土壤序分成 Slight / Intermediate / Strong（Yang 2013 Sect. 3.1）
- `PPDI_PCT` — Yang 2013 Table 2；嵌套格式 `{stage: {order: %}}`
- `STRAIN_EPSILON` — Yang 2013 Table 3；同嵌套格式（Slight 按母质类型键控，其余按土壤序）
- `DEFAULT_SLIGHT_PARENT_MATERIAL = "Non_carbonate"` — SEUS 全域假设（Option 3；后续可用 GLiM 岩性细化）
- `HEDLEY_P_FRACTIONS_PCT` — Yang 2013 Table S2（9 序）+ Yang & Post 2011 Fig. 3（Histosols）；每值是 5 分量字典

### Step 1 — 构建 TPS 图（Eq. 2）
1. **读 `rockP_den`**，kg P/m² → g P/m²（换算因子已吸收 Eq. 2 的 `0.01 · D · ρ_P · C_P`）。
2. **重采样到 1/24° SEUS 网格**：用 `xarray.DataArray.interp` 以 2-D `LATIXY`/`LONGXY` 作为目标（逐点插值）；线性插值 + `fill_value=None`（在海岸格点外推）。
3. **海岸最近邻填补**：用 `cKDTree` 对「陆地 + 有效 rockP」格点的 `(lat, lon)` 建树，填补 `SOIL_ORDER` 为陆地但 `rockP_den` 为 NaN 的格点（Florida Keys、Outer Banks、Mississippi Delta 等）。
4. **Sanity 图**：3 联图 `(全球 rockP | 重采样 rockP | SOIL_ORDER)` 写入 `intermediate/`。
5. **逐像元 PPDI / ε / stage 图**：3 张 length-16 numpy 查找表 + 一次 fancy-indexing。非阶段土壤序默认 PPDI = 0、ε = 0，Eq. 2 退化为 `TPS = rockP_den`。
6. **计算 TPS** `tps = rockP_den · (1 − PPDI/100) / (ε + 1)`；断言 `(ε + 1) > 0`；非陆地格点 NaN 自然传播。
7. **保存**：`outputs/TPS_SEUS_1_24deg.nc` + `intermediate/TPS_*.tif` + 2 联 sanity PNG。

### Step 2 — 把 TPS 拆成 5 个 Hedley P 池
1. **建 `(16, 5)` 查找表** `code_to_fracs_pct`，按 `SOIL_ORDER` 代码索引；无 Hedley 数据的代码保持全零行。
2. **向量化乘法** `p_form_arr = tps[..., None] · fracs / 100` → 形状 `(324, 504, 5)`。
3. **校验** `sum(5 池) == TPS`。
4. **每池一个 xarray DataArray**，带 `long_name` / `units` / `formula` 属性。
5. **保存**：`outputs/P_forms_SEUS_1_24deg.nc`（TPS + 5 池）+ 5 个 GeoTIFF + 2×3 sanity PNG。

### 目录约定
```
s7_updataPdata/
├── outputs/                                  ← 下游消费的最终产品
│   ├── TPS_SEUS_1_24deg.nc
│   └── P_forms_SEUS_1_24deg.nc               ← 02 的输入
└── intermediate/                             ← 临时 / QGIS / sanity
    ├── ppdi_eps_stage_SEUS_1_24deg.nc
    ├── *_SEUS_1_24deg.tif                    (8 个栅格)
    └── check_*.png                           (3 张 sanity 图)
```

---

## `02_Update_Surfdata_Pdata_hpc.py` 详解（Step E）

把 `outputs/P_forms_SEUS_1_24deg.nc` 里的 4 个矿质 P 池写回 s6 输出的 surfdata。
脚本用 `THIS_DIR = Path(__file__).resolve().parent` 锚定自身目录，所以在
`s7_updataPdata/` 下无论从哪里 `sbatch` 启动都能解析到正确路径。

### 输入
- `SURF_IN` — `../s6_updateSoildata/Updatedsoil_surfdata_SEUS_1_24deg_simyr1850_c260712.nc`（s6 输出）
- `P_SRC` — `outputs/P_forms_SEUS_1_24deg.nc`（01 产物；含 `LABILE_P` / `SECONDARY_P` / `APATITE_P` / `OCCLUDED_P`，g P m⁻²，SEUS 1/24°）

### 输出
- `SURF_OUT` — `surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c<today>.nc`，日期戳用今天，避免覆盖旧结果。

### 被改写的变量
| 变量 | 含义 | 单位 |
|---|---|---|
| `LABILE_P` | Labile inorganic phosphorus | g P m⁻² |
| `SECONDARY_P` | Secondary mineral phosphorus | g P m⁻² |
| `APATITE_P` | Apatite phosphorus | g P m⁻² |
| `OCCLUDED_P` | Occluded phosphorus | g P m⁻² |

文件里其它一切（维度、属性、`_FillValue`、压缩）逐比特保留：脚本先
`shutil.copyfile(SURF_IN → SURF_OUT)`，再用 `netCDF4` 以 `r+` 打开副本，只覆盖这 4 个池。

### 算法
1. **读 `SURF_IN`**：加载 `LATIXY` / `LONGXY` / `PFTDATA_MASK`（==1 ⇒ 陆地）及 4 个 P 池作为"原始"参照。
2. **读 `P_SRC`** 并验证与目标在同一 `(lsmlat, lsmlon)` 网格；若 P-forms 带 `LATIXY`/`LONGXY` 则逐元素比对。
3. **最近邻填补**：对 `PFTDATA_MASK==1` 但 P 池为 NaN 的格点，用 `cKDTree` 从有数据的陆地格点借值（每池独立）；非陆地格点最终置 0，下游不会见到 NaN。
4. **`shutil.copyfile` 后 `netCDF4.Dataset(r+)`**：对每池 `out[pft_mask] = filled[v][pft_mask]`（只更新陆地格点），并写 `var.history` 及文件级 `history`。
5. **Sanity 统计**：每池打印陆地格点上的 `orig mean / new mean / orig max / new max`。
6. **GeoTIFF 导出**：每池 `*_old.tif` / `*_new.tif` 到 `intermediate/`（EPSG:4326, 北向上, `nodata=-9999`）。
7. **3 联 sanity PNG**：`(original | new | new − original)` 到 `intermediate/check_surfdata_P_update.png`。

### landuse.timeseries 传播（默认关闭）
脚本含第二块 `UPDATE_LANDUSE_TIMESERIES`，可把同样 4 池写进
`../landuse.timeseries_..._c*.nc`。**默认 `= False` 不执行**：Pathfinder 上
标准 CLM/ELM 输入里 landuse.timeseries 的 P 池全为零，例如
- `/projects/hpcl-cli185/proj-shared/xyk/inputdata/lnd/clm2/surfdata_map/landuse.timeseries_360x720cru_hist_simyr1850-2015_c180220.nc`
- `/projects/hpcl-cli185/proj-shared/zdr/trendy_2024/landuse/landuse.timeseries_360x720cru_hist_TRENDY_simyr1700-2024_c240826.nc`

即 ELM 只从 surfdata 读矿质 P 池。代码保留完好，将来若假设改变，只需把 flag 翻成 `True`。

### 输出目录
```
s7_updataPdata/
├── surfdata_UpdatedsoilP_SEUS_1_24deg_simyr1850_c<today>.nc   ← 最终 surfdata
└── intermediate/
    ├── Pform_{LABILE,SECONDARY,APATITE,OCCLUDED}_P_{old,new}.tif   ← 8 个 GeoTIFF
    └── check_surfdata_P_update.png                                 ← sanity 图
```

---

*文档结构最后整理：2026-07-12（改为「快速上手在前、原理在后」，路径统一到 `Make_surface_data/s7_updataPdata`）*
*方法与决策日志原始日期：2026-05-19 / 05-28 / 06-01 / 06-05*
