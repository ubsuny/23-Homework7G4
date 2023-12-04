import numpy as np

# Define the edge length of the tetrahedron based on expected ion separation
edge_length = 0.3  # This is an arbitrary starting point, in nanometers

# Define the vertices of a regular tetrahedron
# These are normalized coordinates that need to be scaled by the edge length
vertices = np.array([[1, 1, 1],
                     [-1, -1, 1],
                     [-1, 1, -1],
                     [1, -1, -1]])

# Scale the vertices to the edge length
tetrahedron_positions = vertices * (edge_length / np.sqrt(3))

# Assign Na+ and Cl- to the vertices alternately
ideal_geometry = []
for i, pos in enumerate(tetrahedron_positions):
    ion = 'Na' if i % 2 == 0 else 'Cl'
    ideal_geometry.append((ion, pos))

ideal_geometry = np.array(ideal_geometry)
