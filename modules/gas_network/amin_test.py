#%%
import geopandas as gpd
import matplotlib.pyplot as plt
import json
import pandas as pd
#%%

#%%

shape = gpd.read_file("workflow/resources/user/shapes_e_highways.geojson")
# %%
pipes = gpd.read_file (
    r"C:\Users\TAHAVORM\Downloads\GitRepos\ec_modules\modules\gas_network\workflow\resources\automatic\gas_grid\gas_grid\data\IGGIELGN_PipeSegments.geojson"
)
#%%

#%%
pipes["max_cap_M_m3_per_d"] = pipes.param.apply(json.loads).apply(pd.Series)[["max_cap_M_m3_per_d"]]
#
# %%
gas_grid_eu = gpd.sjoin(pipes, shape, how="inner", predicate="intersects")

#%%
fig, axs = plt.subplots(1,1)
shape.plot(ax=axs,facecolor="white",edgecolor= "#c2c2c0")
gas_grid_eu.plot(ax=axs,column='max_cap_M_m3_per_d', cmap='viridis')

# %%
# Perform a spatial join to filter the gas grid within EU borders
gas_grid_eu = gpd.sjoin(pipes, shape, how="inner", predicate="intersects")

# Drop unnecessary columns from the join, if any
gas_grid_eu = gas_grid_eu[pipes.columns]


# %%
