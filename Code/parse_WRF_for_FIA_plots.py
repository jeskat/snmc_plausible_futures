import netCDF4
import pandas as pd
import numpy as np
import glob
from pathlib import Path
import concurrent.futures
import os
import seaborn as sns

os.chdir('../')

# Path and pattern to find .nc files.
# If needed, change and run for different GCMs or projections
FILE_PATTERN = "/pscratch/sd/j/jkatz/WRF/wrf-hist/*.nc"

## Load FIA plot locations
plots_raw = pd.read_csv('Outputs/plots_coordinates_for_prism.csv')
plots_df = plots_raw[['pltID', 'LAT', 'LON']].copy()
plots_df.rename(columns = {
    'LAT':'lat',
    'LON':'lon'
}, inplace = True)
plots_df.set_index('pltID', inplace = True)

## Convert to dictionary
TARGET_POINTS = plots_df.to_dict('index')
TARGET_POINTS

# List of variable names to extract from the .nc files
VARIABLE_NAMES = ['TBOT', 'PRECTmms']

# Define variables and their processing method ('mean' or 'sum')
VARIABLE_PROCESSING = {
    'TBOT': 'mean',        # This will be averaged
    'PRECTmms': 'sum' # This will be summed
}

# Define timestep in seconds (for 'sum' variables). This is pretty good approximation.
TIMESTEP_SECONDS = 3600

# Names of coordinate and time variables in .nc files
LAT_VAR_NAME = 'LATIXY'
LON_VAR_NAME = 'LONGXY'
TIME_VAR_NAME = 'time'

def find_nearest_grid_index(nc_file, target_lat, target_lon):
    """
    Finds the (y, x) indices of the grid cell closest to the target lat/lon.
    
    This function reads the 2D LONGXY and LATIXY variables and uses
    numpy to find the index of the minimum squared distance.
    """
    print(f"    Finding nearest grid cell for (lat={target_lat}, lon={target_lon})...")
    with netCDF4.Dataset(nc_file, 'r') as nc:
        # Load the 2D coordinate arrays
        lats_2d = nc.variables[LAT_VAR_NAME][:]
        lons_2d = nc.variables[LON_VAR_NAME][:]

        # Calculate squared Euclidean distance in lat/lon space
        # This is a good approximation for finding the nearest neighbor
        dist_sq = (lats_2d - target_lat)**2 + (lons_2d - target_lon)**2

        # Find the 1D index of the minimum distance
        min_index_1d = np.argmin(dist_sq)

        # Convert the 1D index to 2D (y, x) indices
        # np.unravel_index gives (y_index, x_index)
        y_idx, x_idx = np.unravel_index(min_index_1d, lats_2d.shape)

        found_lat = lats_2d[y_idx, x_idx]
        found_lon = lons_2d[y_idx, x_idx]
        print(f"    Target: (lat={target_lat}, lon={target_lon})")
        print(f"    Found:  (lat={found_lat:.4f}, lon={found_lon:.4f}) at index [y={y_idx}, x={x_idx}]")

        return y_idx, x_idx

def process_single_file(nc_file_path, unique_cell_list, variable_processing, time_var_name):
    """
    WORKER FUNCTION:
    Opens one .nc file, loops through all unique grid cells,
    and extracts/processes the data for each.
    
    Returns a tuple: (month_year_str, data_for_all_cells)
    
    'data_for_all_cells' is a dict like:
    {(y,x): {'T2_avg': 10, 'PRECIP_total': 3}, (y,x): ...}
    """
    
    # Get the PID for logging
    pid = os.getpid()
    print(f"  (PID {pid}) Processing file: {Path(nc_file_path).name}")

    try:
        with netCDF4.Dataset(nc_file_path, 'r') as nc:
            
            # --- 1. Get Month-Year Identifier ---
            time_var = nc.variables[time_var_name]
            first_time_val = time_var[0]
            time_units = time_var.units
            calendar = getattr(time_var, 'calendar', 'standard')
            dt_obj = netCDF4.num2date(first_time_val, 
                                      units=time_units, 
                                      calendar=calendar)
            month_year_str = dt_obj.strftime('%Y-%m')

            # This dict will hold results for all cells *for this month*
            # Key: (y, x) tuple
            # Value: {'T2_avg': ..., 'PRECIP_total_mm': ...}
            data_for_all_cells = {}

            # --- 2. Loop through all cells we need from this file ---
            for y_idx, x_idx in unique_cell_list:
                
                cell_key = (y_idx, x_idx)
                cell_data = {} # Results for this single cell
                
                for var_name, method in variable_processing.items():
                    if var_name not in nc.variables:
                        # Log warning once per file
                        print(f"  (PID {pid}) Warning: Var '{var_name}' not in {Path(nc_file_path).name}")
                        continue
                        
                    # Extract the time series for this (y, x) point
                    var_data_series = nc.variables[var_name][:, y_idx, x_idx]
                    
                    if method == 'mean':
                        result = np.mean(var_data_series)
                        col_name = f"{var_name}_avg"
                        
                    elif method == 'sum':
                        mm_per_timestep = var_data_series * TIMESTEP_SECONDS
                        result = np.sum(mm_per_timestep)
                        col_name = f"{var_name}_total_mm"
                    
                    cell_data[col_name] = result
                
                # Add this cell's data to the main dictionary
                data_for_all_cells[cell_key] = cell_data

            return (month_year_str, data_for_all_cells)

    except Exception as e:
        print(f"  **ERROR** (PID {pid}) Failed to process file {nc_file_path}: {e}")
        return None # Return None on failure

def main():
    """
    Main processing workflow, parallelized by file.
    """
    
    # --- Find all files (Same as before) ---
    file_list = sorted(glob.glob(FILE_PATTERN, recursive=True))
    if not file_list:
        print(f"Error: No files found matching pattern: {FILE_PATTERN}")
        return
    print(f"Found {len(file_list)} files to process.")

    
    # --- PHASE 1: Map points to grid cells  ---
    print("\n--- Phase 1: Mapping points to grid cells ---")
    
    grid_cell_to_points_map = {} # (y, x) -> [point_name_list]
    
    first_file = file_list[0]
    for point_name, coords in TARGET_POINTS.items():
        print(f"  Mapping '{point_name}'...")
        try:
            y_idx, x_idx = find_nearest_grid_index(first_file, coords['lat'], coords['lon'])
            grid_cell_key = (y_idx, x_idx)
            if grid_cell_key not in grid_cell_to_points_map:
                grid_cell_to_points_map[grid_cell_key] = []
            grid_cell_to_points_map[grid_cell_key].append(point_name)
        except Exception as e:
            print(f"  Error finding grid index for {point_name}: {e}")

    unique_cells_list = list(grid_cell_to_points_map.keys())
    
    print("\n...Point mapping complete.")
    print(f"Total points to process: {len(TARGET_POINTS)}")
    print(f"Total unique grid cells to extract: {len(unique_cells_list)}")


    # --- PHASE 2: Parallel data extraction (by file) ---
    print("\n--- Phase 2: Submitting all file jobs in parallel ---")
    
    # This is our intermediate "collector"
    # It will store a list of results for each cell
    # { (y,x): [ {'month_year':..., 'T2_avg':...}, ... ], ... }
    temp_results = {cell_key: [] for cell_key in unique_cells_list}

    with concurrent.futures.ProcessPoolExecutor() as executor:
        
        # Submit one job for each file
        futures = []
        for nc_file in file_list:
            future = executor.submit(process_single_file, 
                                     nc_file, 
                                     unique_cells_list, 
                                     VARIABLE_PROCESSING, 
                                     TIME_VAR_NAME)
            futures.append(future)

        print(f"  ...All {len(futures)} file-jobs submitted. Waiting for results...")
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(futures):
            result = future.result() # (month_year_str, data_for_all_cells)
            
            if result is None:
                continue # Skip failed jobs

            month_year, data_for_all_cells = result
            
            # --- This is the re-assembly step ---
            # Distribute the results to the correct cell "bucket"
            for cell_key, data in data_for_all_cells.items():
                data['month_year'] = month_year # Add the month
                temp_results[cell_key].append(data)
                
    print("\n...All file-jobs complete. Re-assembly finished.")


    # --- \Finalize DataFrames ---
    # Now, convert the lists of dicts into final, sorted DataFrames
    
    # Build it *after* the parallel step, not during.
    grid_cell_data_cache = {} 
    
    for cell_key, data_list in temp_results.items():
        if not data_list:
            print(f"Warning: No data found for cell {cell_key}")
            continue
            
        df = pd.DataFrame(data_list)
        df = df.set_index('month_year').sort_index()
        grid_cell_data_cache[cell_key] = df
        

    # --- PHASE 3: Assign cached data  ---
    print("\n--- Phase 3: Assigning data to final output ---")
    final_results = {}
    for grid_cell_key, point_names_list in grid_cell_to_points_map.items():
        if grid_cell_key in grid_cell_data_cache:
            cached_df = grid_cell_data_cache[grid_cell_key]
            for point_name in point_names_list:
                final_results[point_name] = cached_df
        else:
            print(f"Warning: No data in cache for {grid_cell_key} (points: {point_names_list})")


    # --- PHASE 4: Save all results  ---
    print("\n--- Workflow Complete ---")
    
    if not final_results:
        print("No results to save.")
        return

    output_dir = "Outputs/Climate_data_extraction"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nSaving results to '{output_dir}/' directory...")

    # --- 6a. Save the Grid Cell-to-Point Map ---
    map_file_path = os.path.join(output_dir, "grid_cell_to_point_map.csv")
    print(f"  Saving grid map to: {map_file_path}")
    map_records = []
    for (y, x), point_names_list in grid_cell_to_points_map.items():
        map_records.append({
            'cell_y_index': y,
            'cell_x_index': x,
            'point_names': ', '.join(point_names_list) 
        })
    map_df = pd.DataFrame(map_records)
    map_df.to_csv(map_file_path, index=False)

    # --- 6b. Save the Data as a SINGLE combined CSV file ---
    print("  Combining all point data into a single DataFrame...")
    all_data_list = []
    for point_name, df in final_results.items():
        df_copy = df.copy() 
        df_copy['point_name'] = point_name
        all_data_list.append(df_copy)

    master_df = pd.concat(all_data_list)
    master_df.reset_index(inplace=True) 
    cols = list(master_df.columns)
    cols = ['point_name', 'month_year'] + [c for c in cols if c not in ['point_name', 'month_year']]
    master_df = master_df[cols]
    master_csv_path = os.path.join(output_dir, "all_points_monthly_data.csv")
    print(f"  Saving combined data to: {master_csv_path}")
    master_df.to_csv(master_csv_path, index=False) 

    print("... All data saved successfully.")


# --- Don't forget this! ---
if __name__ == "__main__":
    main()