import matplotlib.pyplot as plt
import geopandas
import contextily as cx
from shapely.geometry import Point, Polygon
from matplotlib.patches import Circle
from matplotlib_scalebar.scalebar import ScaleBar


file = 'Coastal_Confined_Gravel_Aquifer_System.shp'
map = geopandas.read_file(file)
map = map.to_crs(epsg=4326)
ax = map.boundary.plot(figsize=(5, 3), alpha=0.5, edgecolor='blue', facecolor='lightblue', label="Confined Aquifer System")
file2 = 'location-and-extent-of-nzs-aquifers-2015.shp'
map2 = geopandas.read_file(file2)
map2 = map2.to_crs(epsg=4326)
map2.boundary.plot(alpha=1.0, edgecolor='black', ax=ax, label="Central Plains")
file3 = 'Circles.shp'
# map3 = geopandas.read_file(file3)
# map3 = map3.to_crs(epsg=4326)
# map3.boundary.plot(alpha=0.5, edgecolor='red', facecolor="lightsalmon", ax=ax)
west_melton = [-43.52217, 172.36920]
chch = [-43.53024, 172.63346]
m361159 = [-43.56130, 172.69571]
ax.set_xlim(172.2, 173.2)
ax.set_ylim(-43.75, -43.25)
cx.add_basemap(ax, crs=4326, zoom=10)
s1 = ax.scatter(chch[1], chch[0], s=100, marker="*", c="yellow", edgecolors="black", label="Christchurch", zorder=10)
s2 = ax.scatter(west_melton[1], west_melton[0], s=100, marker="*", c="lightgreen", edgecolors="black", label="West Melton", zorder=10)
s3 = ax.scatter(m361159[1], m361159[0], s=50, marker="o", c="red", edgecolors="black", label="Well M36/1159", zorder=10)
ax.legend()
plt.show()
pass



