import matplotlib.pyplot as plt

# List of pairs of coordinates, where each pair is a line segment
# Each element is ((x1, y1), (x2, y2))
line_segments = [
    ((80.000000, 41.000000), (80.000000, 25.000000)),
    ((80.000000, 25.000000), (72.000000, 25.000000)),
]

# Create a plot
plt.figure(figsize=(8, 6))

# Plot each line segment
for segment in line_segments:
    start_point, end_point = segment
    x_values = [start_point[0], end_point[0]]
    y_values = [start_point[1], end_point[1]]

    plt.plot(x_values, y_values, marker='o', linestyle='-', markersize=5)  # 'o' for circle markers

# Add titles and labels
plt.title('Line Segments')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.grid(True)  # Optionally add a grid

# Show the plot
plt.show()

