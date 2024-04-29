import matplotlib.pyplot as plt

points = [
    (254.919998, 271.649994), (418.920013, 200.149994),
    (582.419983, 321.149994), (155.419998, 150.649994),
]
points_ref2 = [
    (418.920013, 200.149994), (176.919998, 328.149994),
    (582.419983, 321.149994), (155.419998, 150.649994),
]
points = [
    (418.920013, 200.149994), (176.919998, 328.149994),
    (582.419983, 321.149994), (155.419998, 150.649994),
]


plt.figure(figsize=(8, 6))

point1, point, point2, point3 = points

x_values = [point1[0], point[0]]
y_values = [point1[1], point[1]]
plt.plot(x_values, y_values, marker='o', color="red", linestyle='--', markersize=5)  # 'o' for circle markers

x_values = [point[0], point2[0]]
y_values = [point[1], point2[1]]
plt.plot(x_values, y_values, marker='o', color="red", linestyle='-', markersize=5)  # 'o' for circle markers

x_values = [point2[0], point3[0]]
y_values = [point2[1], point3[1]]
plt.plot(x_values, y_values, marker='o', color="violet", linestyle='-', markersize=5)  # 'o' for circle markers

plt.plot(point[0], point[1], marker='x', color="green", markersize=10)  # 'o' for circle markers

x_values = [point1[0], point2[0]]
y_values = [point1[1], point2[1]]
plt.plot(x_values, y_values, marker='o', color="blue", linestyle='--', markersize=5)  # 'o' for circle markers

# Add titles and labels
plt.title('Line Segments')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.grid(True)  # Optionally add a grid

# Show the plot
plt.show()

