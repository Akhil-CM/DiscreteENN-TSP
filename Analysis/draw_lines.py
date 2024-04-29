import matplotlib.pyplot as plt

def plotPoints(points, color):
    points_num = len(points)
    for idx in range(points_num - 1):
        point1 = points[idx]
        if idx == points_num - 1:
            point2 = points[0]
        else:
            point2 = points[idx + 1]
        if not point1 or not point2:
            continue
        x_values = [point1[0], point2[0]]
        y_values = [point1[1], point2[1]]
        plt.plot(x_values, y_values, marker='o', color=color, linestyle='--', markersize=5)  # 'o' for circle markers

plt.figure(figsize=(8, 6))

points = [
    (850.000000, 700.000000), (900.000000, 600.000000),
    (800.000000, 600.000000), (850.000000, 520.000000),
    (750.000000, 490.000000)
]
plotPoints(points, "red")

points = [
    (850.000000, 700.000000), (800.000000, 600.000000),
    (900.000000, 600.000000), (850.000000, 520.000000),
    (750.000000, 490.000000)
]
plotPoints(points, "blue")

# Add titles and labels
plt.title('Line Segments')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.grid(True)  # Optionally add a grid

# Show the plot
plt.show()

