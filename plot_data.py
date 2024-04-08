import numpy as np
import os
from os import path as osPath
from matplotlib import pyplot as plt

DATA_DIR = os.getcwd()
DATA_FILENAME = "./enn_tsp_data.csv"

DATA_FILE = osPath.join(DATA_DIR, DATA_FILENAME)

def plotData(data, load_data: bool, sequence: bool):
    data = np.loadtxt(data, delimiter=',') if load_data else data
    # data = np.genfromtxt(data, delimiter=',') if load_data else data

    if data.shape[1] < 2:
        print(f"[Error] (plotData): Data without two columns cannot be used for plotting x and y.")
        exit(1)
    x = data[:, 0]  # First column
    y = data[:, 1]  # Second column

    print(f"{x}\t{y}")
    plt.figure(figsize=(10, 6))  # Optional: Adjusts the figure size
    plt.plot(x, y, marker='o', linestyle='-')  # Adjust markers and lines as needed
    plt.title('Number of Cites vs Time (seconds)')
    plt.xlabel('Cities')  # Adjust with actual column name or description
    plt.ylabel('Time (s)')  # Adjust with actual column name or description
    if sequence:
        x_labels = data[:, 0]  # First column
        x = np.arange(1, len(x_labels) + 1)
        axis = plt.gca()
        plt.xticks(x, labels = x_labels)
    plt.grid(True)  # Optional: Adds a grid
    plt.show()

if __name__ == "__main__":
    plotData(DATA_FILE, load_data = True, sequence = False)
