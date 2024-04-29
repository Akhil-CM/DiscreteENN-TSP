import math
import matplotlib.pyplot as plt

def getArea(coords):
    a, b, c = coords
    return abs(a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))/2.0

def isInside1(coords):
    p, a, b, c = coords
    tolerance = 1e-8
    area_total = getArea([a, b, c])

    area1 = getArea([p, a, b])
    area2 = getArea([p, b, c])
    area3 = getArea([p, c, a])

    if abs(area1) < tolerance or abs(area2) < tolerance or abs(area3) < tolerance:
        print(f"area1 : {area1}, area2 : {area2}, area3 : {area3}")
        return True

    area_diff = abs(area_total - (area1 + area2 + area3))
    print(f"area_diff : {area_diff}")
    return area_diff < tolerance

def isCollinearAndBetween(coords):
    c, a, b = coords
    return min(a[0], b[0]) <= c[0] and c[0] <= max(a[0], b[0]) and\
            min(a[1], b[1]) <= c[1] and c[1] <= max(a[1], b[1])

def isInside2(coords):
    p, a, b, c = coords
    tolerance = 1e-8
    area_total = getArea([a, b, c])

    area1 = getArea([p, a, b])
    area2 = getArea([p, b, c])
    area3 = getArea([p, c, a])

    if abs(area1) < tolerance:
        print(f"area1 : {area1}")
        return isCollinearAndBetween([p, a, b])
    if abs(area2) < tolerance:
        print(f"area2 : {area2}")
        return isCollinearAndBetween([p, b, c])
    if abs(area3) < tolerance:
        print(f"area3 : {area3}")
        return isCollinearAndBetween([p, c, a])

    area_diff = abs(area_total - (area1 + area2 + area3))
    print(f"area_diff : {area_diff}")
    return area_diff < tolerance


def plotTriangle(coords):
    separate_point = coords[0]
    triangle_points = coords[1:]

    plt.plot(separate_point[0], separate_point[1], 'bo')  # 'bo' for blue dot

    triangle_points.append(triangle_points[0])
    xs, ys = zip(*triangle_points)  # This creates two lists: one with all x coordinates, and one with all y coordinates

    plt.plot(xs, ys, 'r-')  # 'r-' for red lines

    for point in triangle_points[:-1]:  # Exclude the last point since it's a repeat of the first
        plt.plot([separate_point[0], point[0]], [separate_point[1], point[1]], 'g--')  # 'g--' for green dotted lines

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Triangle with Dotted Lines to Separate Point')
    plt.grid(True)
    plt.axis('equal')  # Ensures equal scaling for both axes

    plt.show()

if __name__ == "__main__":
    coords = [[945, 685], [700, 580], [770, 610], [700, 500]]
    plotTriangle(coords)
    is_inside = isInside1(coords)
    print(f"coords : {coords}\nisInside : {is_inside}")
    is_inside = isInside2(coords)
    print(f"coords : {coords}\nisInside : {is_inside}")
    coords = [[710, 600], [700, 580], [770, 610], [700, 500]]
    plotTriangle(coords)
    is_inside = isInside1(coords)
    print(f"coords : {coords}\nisInside : {is_inside}")
    is_inside = isInside2(coords)
    print(f"coords : {coords}\nisInside : {is_inside}")
