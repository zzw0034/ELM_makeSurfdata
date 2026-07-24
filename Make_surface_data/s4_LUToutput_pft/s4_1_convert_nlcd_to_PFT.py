# convert North America 1850-2023 NLCD land cover percentage to ELM PFT land cover percentage
# Input: convert nlcd_frac_pred_1850-2023_1_24deg.nc
# Output: elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc

# note: carafully deal with NaN values in the input data. replace NaN with 0.0 before processing.
# PCT_NATVEG: sum of all vegetation classes change by year 1850 to year 2023

# %%
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
import os

# %%
# read in data
ncfile = Dataset(
    "/projects/hpcl-cli185/proj-shared/zdr/hires_data/AI_nlcd/nlcd_frac_pred_1850-2023_1_24deg.nc",
    "r",
)
# copy from Pathfinder Dan folder to local folder
# Pathfinder:
# /projects/hpcl-cli185/proj-shared/zdr/hires_data/AI_nlcd/nlcd_frac_pred_1850-2023_1_24deg.nc
# %%
# create a new file as ouput result with the same dimensions and ELM PFT variables
output_file = "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/scr_out/elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc"

# Ensure output directory exists
output_dir = os.path.dirname(output_file)
if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    print(f"Created output directory: {output_dir}")

# Remove file if it exists to avoid permission errors
if os.path.exists(output_file):
    try:
        os.remove(output_file)
        print(f"Removed existing file: {output_file}")
    except PermissionError:
        print(f"Error: Cannot remove existing file {output_file}.")
        print(
            "The file may be open in another program (e.g., another Python session, xarray, or a viewer)."
        )
        print("Please close the file and try again.")
        raise
    except Exception as e:
        print(f"Warning: Error removing existing file: {e}")
        raise

new_ncfile = Dataset(
    output_file,
    "w",
    format="NETCDF4",
)

# Get dimensions from input file (more robust than hardcoding)
ntime = len(ncfile.dimensions["time"])
nlat = len(ncfile.dimensions["lat"])
nlon = len(ncfile.dimensions["lon"])

time_dim = new_ncfile.createDimension("time", ntime)
lat_dim = new_ncfile.createDimension("lat", nlat)
lon_dim = new_ncfile.createDimension("lon", nlon)
pft_dim = new_ncfile.createDimension("natpft", 17)  # Land cover indices
urb_dim = new_ncfile.createDimension("numurbl", 3)


# Helper function to create variable with optional _FillValue
def create_var_with_fill(new_file, var_name, source_var, dims):
    """Create variable, copying _FillValue if it exists in source."""
    fill_val = getattr(source_var, "_FillValue", None)
    if fill_val is not None:
        var = new_file.createVariable(
            var_name, source_var.dtype, dims, fill_value=fill_val
        )
    else:
        var = new_file.createVariable(var_name, source_var.dtype, dims)
    return var


# set the actual values of the dimensions from the original ncfile (nlcd_frac_pred_1850-2023_1_24deg.nc) to
# the new ncfile (elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc)
# create lat variable
lat_var = create_var_with_fill(new_ncfile, "lat", ncfile.variables["lat"], ("lat",))
# assign the values of the lat variable
lat_var[:] = ncfile.variables["lat"][:]
# create lon variable
lon_var = create_var_with_fill(new_ncfile, "lon", ncfile.variables["lon"], ("lon",))
# assign the values of the lon variable
lon_var[:] = ncfile.variables["lon"][:]
# create time variable
time_var = create_var_with_fill(new_ncfile, "time", ncfile.variables["time"], ("time",))
# assign the values of the time variable
time_var[:] = ncfile.variables["time"][:]

# copy the attributes of the original variables (skip _FillValue which must be set at creation)
for var_name in ["lat", "lon", "time"]:
    for attr in ncfile.variables[var_name].ncattrs():
        if attr != "_FillValue":  # Skip _FillValue - must be set during createVariable
            setattr(
                new_ncfile.variables[var_name],
                attr,
                getattr(ncfile.variables[var_name], attr),
            )

# Output file:
# create output variables (ELM PFT variables) in the output ncfile (elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc)
# PCT_NAT_PFT, PCT_NATVEG, PCT_URBAN
pct_nat_pft_var = new_ncfile.createVariable(
    "PCT_NAT_PFT", np.float64, ("time", "natpft", "lat", "lon")
)
pct_natveg_var = new_ncfile.createVariable(
    "PCT_NATVEG", np.float64, ("time", "lat", "lon")
)

pct_urban_var = new_ncfile.createVariable(
    "PCT_URBAN", np.float64, ("time", "numurbl", "lat", "lon")
)
# %%
# translate nlcd land cover percentage to ELM PFT land cover percentage
# this part is from Dan's code (nlcd_to_elmpft.py)
fc4 = 0.2  # Do this later
veg_indices = [
    # 11,  # 11, is water
    # 12,  # 12, is perennial ice/snow
    21,
    31,
    41,
    42,
    43,
    51,
    52,
    71,
    72,
    73,
    74,
    81,
    82,  # cultivated crops
    90,
    95,
]  # 11, is water
veg_index_names = [f"lc_{idx}" for idx in veg_indices]
veg_vars = [v for v in veg_index_names if v in ncfile.variables]
# veg_vars = ['lc_21', 'lc_31', 'lc_41', 'lc_42', 'lc_43', 'lc_52', 'lc_71', 'lc_81', 'lc_82', 'lc_90', 'lc_95']


def _lc(code: int, t: int) -> np.ndarray:
    """Return lc_<code>[t,:,:] if present; otherwise zeros (time,lat,lon) compatible.
    NaN values are replaced with 0.0."""
    name = f"lc_{code}"
    if name in ncfile.variables:
        data = ncfile.variables[name][t, :, :]
        # Replace NaN with 0.0
        data = np.where(np.isnan(data), 0.0, data)
        return data
    return np.zeros((nlat, nlon), dtype=np.float64)


# %%
for t in range(ntime):
    pft = np.zeros((17, nlat, nlon), dtype=np.float64)
    urb = np.zeros((3, nlat, nlon), dtype=np.float64)

    pft[0, :, :] = _lc(31, t)  # 0, is bare ground

    pft[1, :, :] = (
        _lc(42, t) + (_lc(43, t) + _lc(90, t)) * 0.5
    )  # 1 Temperature evergreen
    pft[7, :, :] = (
        _lc(41, t) + (_lc(43, t) + _lc(90, t)) * 0.5
    )  # 7 Temperature deciduous
    pft[9, :, :] = (_lc(52, t) + _lc(51, t)) * 0.5
    pft[10, :, :] = (_lc(52, t) + _lc(51, t)) * 0.5
    pft[13, :, :] = (1.0 - fc4) * (_lc(71, t) + _lc(72, t) + _lc(81, t)) + _lc(
        95, t
    )  # C3 Grass
    pft[14, :, :] = (fc4) * (_lc(71, t) + _lc(72, t) + _lc(81, t))  # C4 Grass
    pft[15, :, :] = _lc(82, t)  # Cultivated crop
    # Assume urban open space is grass/tree mix
    pft[7, :, :] = pft[7, :, :] + _lc(21, t) * 0.5  # C21 Open Space
    pft[13, :, :] = pft[13, :, :] + _lc(21, t) * 0.5

    # Normalize to 100%
    sums = np.sum(pft, axis=0)  # (lat, lon)
    mask = sums > 1e-6

    # Normalize where sum is valid
    # pft shape: (17, lat, lon), mask shape: (lat, lon)
    # Use broadcasting: pft[:, mask] selects all PFTs at masked grid points
    pft[:, mask] = pft[:, mask] / sums[mask] * 100.0

    # Where the sum is <= threshold: set all PFTs to 0
    # Note: Setting PFT[0]=100% here may cause total > 100% if other types (urban, lake) exist
    # Consider removing the next line if you have other land cover types
    pft[:, ~mask] = 0.0
    pft[0, ~mask] = (
        100.0  # Safe as long as PCT_NAT_PFT is only used where PCT_NATVEG > 0  # (i.e., gated by PCT_NATVEG; ignored when PCT_NATVEG == 0).
    )

    # Urban classes
    urb[0, :, :] = _lc(22, 0)
    urb[1, :, :] = _lc(23, 0)
    urb[2, :, :] = _lc(24, 0)
    # ELM Urban donot change with time。
    # Separate urban data are used by mksurfdata to reduce the vegetated land unit to make room for the urban land unit.
    # LUH urban is not used by mksurfdata.
    # So, we keep the urban data unchanged with time.
    # lc(22, 0)  0 means the urban data in year 1850, lc(22, t) means the urban data in year t.

    # Handle NaN in urban data: replace NaN with 0
    urb = np.where(np.isnan(urb), 0.0, urb)

    # Calculate PCT_NATVEG: sum of all vegetation classes
    s = np.zeros((nlat, nlon), dtype=np.float64)
    for v in veg_vars:
        veg_data = ncfile.variables[v][t, :, :]
        # Replace NaN with 0 before adding
        veg_data = np.where(np.isnan(veg_data), 0.0, veg_data)
        s += veg_data
    pct_natveg_var[t, :, :] = s
    # PCT_NATVEG: sum of all vegetation classes change by year 1850 to year 2023
    # Write all variables for this time step
    new_ncfile.variables["PCT_NAT_PFT"][t, :, :, :] = pft
    new_ncfile.variables["PCT_NATVEG"][t, :, :] = pct_natveg_var[t, :, :]
    new_ncfile.variables["PCT_URBAN"][t, :, :, :] = urb

"""
#%% test------------------------------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.pcolormesh(pct_natveg_var[0, :, :], cmap="viridis", shading="auto")
plt.colorbar(label="PCT_NATVEG (%)")
plt.title("PCT_NATVEG at t=0")
plt.xlabel("Longitude Index")
plt.ylabel("Latitude Index")
plt.tight_layout()
plt.show()

#%% test PCT_NATVEG------------------------------------------------------------------------------
# %% test------------------------------------------------------------------------------
plt.figure(figsize=(8, 6))
# PCT_NAT_PFT has dims (time, natpft, lat, lon)
plt.pcolormesh(pct_nat_pft_var[0, 0, :, :], cmap="viridis", shading="auto")
plt.colorbar(label="PCT_NAT_PFT (%)")
plt.title("PCT_NAT_PFT at t=0, natpft=0")
plt.xlabel("Longitude Index")
plt.ylabel("Latitude Index")
plt.tight_layout()
plt.show()
# %% test end------------------------------------------------------------------------------
import xarray as xr
ds = xr.open_dataset(output_file)
ds["PCT_NAT_PFT"].isel(time=0).plot(x="lon", y="lat", col="natpft", col_wrap=3)
# %% test end------------------------------------------------------------------------------
"""
# %%
# Check for NaN values in the output variables
print("\n" + "=" * 60)
print("Checking for NaN values in output variables")
print("=" * 60)

# Check PCT_NAT_PFT
pct_nat_pft_data = new_ncfile.variables["PCT_NAT_PFT"][:]
nan_count_pct_nat_pft = np.isnan(pct_nat_pft_data).sum()
total_elements_pct_nat_pft = pct_nat_pft_data.size
pct_nan_pct_nat_pft = (
    (nan_count_pct_nat_pft / total_elements_pct_nat_pft) * 100
    if total_elements_pct_nat_pft > 0
    else 0
)

print(f"\nPCT_NAT_PFT:")
print(f"  Total elements: {total_elements_pct_nat_pft}")
print(f"  NaN count: {nan_count_pct_nat_pft}")
print(f"  NaN percentage: {pct_nan_pct_nat_pft:.4f}%")
if nan_count_pct_nat_pft > 0:
    # Find where NaNs are located
    nan_indices = np.where(np.isnan(pct_nat_pft_data))
    print(f"  NaN locations: {len(nan_indices[0])} positions")
    if len(nan_indices[0]) > 0:
        print(f"    First few NaN positions (time, natpft, lat, lon):")
        for i in range(min(5, len(nan_indices[0]))):
            print(
                f"      ({nan_indices[0][i]}, {nan_indices[1][i]}, {nan_indices[2][i]}, {nan_indices[3][i]})"
            )
else:
    print("  ✓ No NaN values found")

# Check PCT_NATVEG
pct_natveg_data = new_ncfile.variables["PCT_NATVEG"][:]
nan_count_pct_natveg = np.isnan(pct_natveg_data).sum()
total_elements_pct_natveg = pct_natveg_data.size
pct_nan_pct_natveg = (
    (nan_count_pct_natveg / total_elements_pct_natveg) * 100
    if total_elements_pct_natveg > 0
    else 0
)

print(f"\nPCT_NATVEG:")
print(f"  Total elements: {total_elements_pct_natveg}")
print(f"  NaN count: {nan_count_pct_natveg}")
print(f"  NaN percentage: {pct_nan_pct_natveg:.4f}%")
if nan_count_pct_natveg > 0:
    nan_indices = np.where(np.isnan(pct_natveg_data))
    print(f"  NaN locations: {len(nan_indices[0])} positions")
    if len(nan_indices[0]) > 0:
        print(f"    First few NaN positions (time, lat, lon):")
        for i in range(min(5, len(nan_indices[0]))):
            print(
                f"      ({nan_indices[0][i]}, {nan_indices[1][i]}, {nan_indices[2][i]})"
            )
else:
    print("  ✓ No NaN values found")

# Check PCT_URBAN
pct_urban_data = new_ncfile.variables["PCT_URBAN"][:]
nan_count_pct_urban = np.isnan(pct_urban_data).sum()
total_elements_pct_urban = pct_urban_data.size
pct_nan_pct_urban = (
    (nan_count_pct_urban / total_elements_pct_urban) * 100
    if total_elements_pct_urban > 0
    else 0
)

print(f"\nPCT_URBAN:")
print(f"  Total elements: {total_elements_pct_urban}")
print(f"  NaN count: {nan_count_pct_urban}")
print(f"  NaN percentage: {pct_nan_pct_urban:.4f}%")
if nan_count_pct_urban > 0:
    nan_indices = np.where(np.isnan(pct_urban_data))
    print(f"  NaN locations: {len(nan_indices[0])} positions")
    if len(nan_indices[0]) > 0:
        print(f"    First few NaN positions (time, numurbl, lat, lon):")
        for i in range(min(5, len(nan_indices[0]))):
            print(
                f"      ({nan_indices[0][i]}, {nan_indices[1][i]}, {nan_indices[2][i]}, {nan_indices[3][i]})"
            )
else:
    print("  ✓ No NaN values found")

print("\n" + "=" * 60)
if (
    nan_count_pct_nat_pft == 0
    and nan_count_pct_natveg == 0
    and nan_count_pct_urban == 0
):
    print("✓ All variables are free of NaN values")
else:
    print("⚠ Some variables contain NaN values - please review")
print("=" * 60)

ncfile.close()
new_ncfile.close()

print("Done!")
print("Output file: ")
print("elmpft_from_nlcd_frac_pred_1850-2023_1_24deg.nc ")
print(
    "/projects/hpcl-cli185/proj-shared/zw5/ELM_makeSurfdata/Make_surface_data/s4_LUToutput_pft/"
)

# %%
