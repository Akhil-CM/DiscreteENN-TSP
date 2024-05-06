import matplotlib.pyplot as plt
import numpy as np

def parse_coordinates(filename):
    """
    Parse a TSP file to extract node coordinates.

    :param filename: Path to the TSP file.
    :return: List of tuples, each containing coordinates (x, y).
    """
    coordinates = []
    start_parsing = False
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "NODE_COORD_SECTION":
                start_parsing = True
                continue
            if line == "EOF":
                break
            if start_parsing:
                parts = line.split()
                # The coordinates are expected to be in parts[1] (x) and parts[2] (y)
                x, y = float(parts[1]), float(parts[2])
                coordinates.append((x, y))
    return coordinates

def plot_coordinates(coordinates):
    """
    Plot the coordinates on a scatter plot with red points and a white background.

    :param coordinates: List of tuples, each containing coordinates (x, y).
    """
    # Unpack the list of tuples into x and y coordinates
    x, y = zip(*coordinates)
    max_x = max(x, default=0)
    max_y = max(y, default=0)
    max_coord = max(max_x, max_y)
    min_x = min(x, default=0)
    min_y = min(y, default=0)
    min_coord = min(min_x, min_y)
    x, y = np.array(x), max_y - np.array(y)

    plt.figure(figsize=(10, 8))  # Set the figure size
    plt.scatter(x, y, color='red')  # Plot with red points
    plt.title('TSP Node Locations')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.grid(True)  # Enable grid for better visualization
    plt.gca().set_facecolor('white')  # Set background to white
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)
    plt.show()
    print(f"min and max are : ({min_coord}, {max_coord})")

# Example usage
filename = './Data/ALL_tsp/berlin52.tsp'  # Update this with the path to your TSP file
coordinates = parse_coordinates(filename)
plot_coordinates(coordinates)
