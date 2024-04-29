import math

def getDistance(a, b):
    x_sqr = (a[0] - b[0]) * (a[0] - b[0])
    y_sqr = (a[1] - b[1]) * (a[1] - b[1])
    return math.sqrt(x_sqr + y_sqr)

def getCost(c, a, b):
    return getDistance(c, a) + getDistance(c, b) - getDistance(a, b)

if __name__ == "__main__":
    point = (800.000000, 600.000000)
    pointA = (900.000000, 600.000000)
    pointB = (850.000000, 520.000000)
    print(f"cost for ({pointA}, {point}, {pointB}):\n{getCost(point, pointA, pointB)}")
