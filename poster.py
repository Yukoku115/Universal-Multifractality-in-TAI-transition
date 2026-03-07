
import osmnx as ox
import matplotlib.pyplot as plt
import pandas as pd

# 1. Define the specific location
# You can use "Manchester, United Kingdom" for just the city center, 
# or "Greater Manchester, United Kingdom" for the wider sprawl.
place_name = "Greater Manchester, United Kingdom"

# 2. Download the street network from OpenStreetMap
# We use network_type='drive' to get the main physical roads.
# 'simplify=True' is crucial: it removes the curved points along a road 
# and ONLY keeps the actual junctions/intersections!
print(f"Downloading street network for {place_name}. This might take a minute...")
G = ox.graph_from_place(place_name, network_type='drive', simplify=True)

# 3. Convert the graph into a GeoDataFrame to easily access the nodes
nodes, edges = ox.graph_to_gdfs(G)

# 4. Isolate the coordinate data (The SIPP)
# We drop all the street data and just keep the X (Longitude) and Y (Latitude) of the nodes
intersections = nodes[['x', 'y']]

# 5. Save the data to a CSV for your physics math later!
filename = "manchester_sipp.csv"
intersections.to_csv(filename, index=False)
print(f"Success! Extracted {len(intersections)} intersections and saved to {filename}.")

# 6. Visualize the data (This will be Figure 1 on your poster!)
plt.figure(figsize=(10, 10), facecolor='white')
# We plot them as tiny black dots (s=0.1) to see the density clouds
plt.scatter(intersections['x'], intersections['y'], s=0.1, c='black', alpha=0.5)

plt.title("Street Intersection Point Pattern (SIPP) - Greater Manchester", fontsize=16)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
# This ensures the map doesn't look stretched out
plt.axis('equal') 

# Remove the box border for a cleaner, "physics" look
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.show()