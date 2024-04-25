import matplotlib.pyplot as plt

points = [
    (418.920013, 200.149994), (176.919998, 328.149994),
    (582.419983, 321.149994), (155.419998, 150.649994),
]


plt.figure(figsize=(8, 6))

points_num = len(points)
for idx in range(points_num):
    point1 = points[idx]
    if idx == points_num - 1:
        point2 = points[0]
    else:
        point2 = points[idx + 1]
    x_values = [point1[0], point2[0]]
    y_values = [point1[1], point2[1]]
    plt.plot(x_values, y_values, marker='o', color="red", linestyle='--', markersize=5)  # 'o' for circle markers

# Add titles and labels
plt.title('Line Segments')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.grid(True)  # Optionally add a grid

# Show the plot
plt.show()

