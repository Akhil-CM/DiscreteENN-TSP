import numpy as np
import math
import os
from os import path as osPath
from matplotlib import pyplot as plt
import pandas as pd

DATA_DIR = os.getcwd()
DATA_FILENAME1 = "./enn_tsp_data.csv"
DATA_FILENAME2 = "./DiscreteENN_TSP_table.csv"

DATA_FILE1 = osPath.join(DATA_DIR, DATA_FILENAME1)
DATA_FILE2 = osPath.join(DATA_DIR, DATA_FILENAME2)

def plotData(data_file, use_pandas: bool, sequence: bool):
    time_str = "$\mu$s"
    show_grid = False

    if use_pandas:
        color = "red"
        data = pd.read_csv(data_file)
        # print(f"Pandas dataframe is:\n{data}")
        # exit(1)
        x = data.iloc[:, 1].to_numpy()
        y = data.iloc[:, 3].to_numpy()
        x_labels = (data.iloc[:, 1].to_numpy()).astype(int)  # First column
        # xy_pd = pd.DataFrame({"x" : x, "y": y})
        # print(f"Pandas dataframe is:\n{data}\nx and y are:\n{xy_pd}")
        # exit(1)
    else:
        color = "blue"
        data = np.loadtxt(data_file, delimiter=',')
        # data = np.genfromtxt(data_file, delimiter=',')
        if data.shape[1] < 2:
            print(f"[Error] (plotData): Data without two columns cannot be used for plotting x and y.")
            exit(1)
        x = data[:, 0]  # First column
        y = data[:, 1]  * 1e6 # Second column
        x_labels = data[:, 0].astype(int)  # First column

    y1 = y/x
    y2 = y/(x*x)
    y3 = y/(x*np.log(x))
    if sequence:
        x = np.arange(1, len(x_labels) + 1)

    plt.figure(figsize=(10, 10))  # Optional: Adjusts the figure size
    plt.subplot(131)
    plt.plot(x, y1, color=color, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.title(f"Time ({time_str}) / N")
    plt.xlabel("Cities")  # Adjust with actual column name or description
    # plt.ylabel(f"Time ({time_str})")  # Adjust with actual column name or description
    plt.grid(show_grid)  # Optional: Adds a grid
    if sequence:
        axis = plt.gca()
        plt.xticks(x, labels = x_labels, rotation='vertical')
    plt.subplot(132)
    plt.plot(x, y2, color=color, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.title(f"Time ({time_str}) / $N^2$")
    plt.xlabel('Cities')  # Adjust with actual column name or description
    # plt.ylabel(f"Time ({time_str})")  # Adjust with actual column name or description
    plt.grid(show_grid)  # Optional: Adds a grid
    if sequence:
        axis = plt.gca()
        plt.xticks(x, labels = x_labels, rotation='vertical')
    plt.subplot(133)
    plt.plot(x, y3, color=color, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.title(f"Time ({time_str}) / $N \cdot log(N)$")
    plt.xlabel('Cities')  # Adjust with actual column name or description
    # plt.ylabel(f"Time ({time_str})")  # Adjust with actual column name or description
    plt.grid(show_grid)  # Optional: Adds a grid
    if sequence:
        axis = plt.gca()
        plt.xticks(x, labels = x_labels, rotation='vertical')
    if use_pandas:
        plt.suptitle("Data from current implementation")
    else:
        plt.suptitle("Data from paper")
    plt.show()

def plotDataBoth(data_file1, data_file2, sequence: bool):
    time_str = "$\mu$s"
    show_grid = False

    color1 = "blue"
    data1 = np.loadtxt(data_file1, delimiter=',')
    # data1 = np.genfromtxt(data_file1, delimiter=',')
    if data1.shape[1] < 2:
        print(f"[Error] (plotData): Data without two columns cannot be used for plotting x and y.")
        exit(1)
    x1 = data1[:, 0]  # First column
    y1 = data1[:, 1]  * 1e6 # Second column
    x1_labels = data1[:, 0].astype(int)  # First column

    color2 = "red"
    data2 = pd.read_csv(data_file2)
    x2 = data2.iloc[:, 1].to_numpy()
    y2 = data2.iloc[:, 3].to_numpy()
    x2_labels = (data2.iloc[:, 1].to_numpy()).astype(int)  # First column

    if sequence:
        x_labels = x1_labels
        x = np.arange(1, len(x_labels) + 1)
    else:
        x = x1

    y1_1 = y1/x1
    y1_2 = y1/(x1*x1)
    y1_3 = y1/(x1*np.log(x1))

    y2_1 = y2/x2
    y2_2 = y2/(x2*x2)
    y2_3 = y2/(x2*np.log(x2))

    plt.figure(figsize=(10, 10))  # Optional: Adjusts the figure size
    print(f"{x}\n{y1_1}")
    print(f"{x}\n{y2_1}")
    plt.subplot(131)
    plt.plot(x, y1_1, color=color1, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.plot(x, y2_1, color=color2, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.title(f"Time ({time_str}) / N")
    plt.xlabel("Cities")  # Adjust with actual column name or description
    # plt.ylabel(f"Time ({time_str})")  # Adjust with actual column name or description
    plt.grid(show_grid)  # Optional: Adds a grid
    if sequence:
        axis = plt.gca()
        plt.xticks(x, labels = x_labels, rotation='vertical')
    print(f"{x}\n{y1_2}")
    print(f"{x}\n{y2_2}")
    plt.subplot(132)
    plt.plot(x, y1_2, color=color1, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.plot(x, y2_2, color=color2, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.title(f"Time ({time_str}) / $N^2$")
    plt.xlabel('Cities')  # Adjust with actual column name or description
    # plt.ylabel(f"Time ({time_str})")  # Adjust with actual column name or description
    plt.grid(show_grid)  # Optional: Adds a grid
    if sequence:
        axis = plt.gca()
        plt.xticks(x, labels = x_labels, rotation='vertical')
    print(f"{x}\n{y1_3}")
    print(f"{x}\n{y2_3}")
    plt.subplot(133)
    plt.plot(x, y1_3, color=color1, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.plot(x, y2_3, color=color2, marker='o', linestyle='-.')  # Adjust markers and lines as needed
    plt.title(f"Time ({time_str}) / $N \cdot log(N)$")
    plt.xlabel('Cities')  # Adjust with actual column name or description
    # plt.ylabel(f"Time ({time_str})")  # Adjust with actual column name or description
    plt.grid(show_grid)  # Optional: Adds a grid
    if sequence:
        axis = plt.gca()
        plt.xticks(x, labels = x_labels, rotation='vertical')
    plt.show()

if __name__ == "__main__":
    plotData(DATA_FILE1, False, sequence = True)
    plotData(DATA_FILE2, True, sequence = True)
    plotDataBoth(DATA_FILE1, DATA_FILE2, sequence = True)
