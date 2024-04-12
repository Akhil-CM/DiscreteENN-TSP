

Data_File = "output.log"

if __name__ == "__main__":
    distance_min = 10e8
    distance_max = -1.0
    distance_avg = 0.0
    time_min = -1.0
    time_max = -1.0
    time_avg = 0.0
    time_curr = -1.0
    line_min = 0
    line_max = 0
    line_curr = 0
    iter_min = 0
    iter_max = 0
    iter_curr = 0
    with open(Data_File, "r") as input_file:
        for line in input_file:
            line = line.rstrip()
            if not line: continue
            line_curr += 1
            if "Algorithm finished in " in line:
                # time_curr = line[line.find(" in ") + 3:]
                time_curr = float(line[line.find(" in ") + 3:line.find("ms")].strip())
            if "Total distance" not in line: continue
            iter_curr += 1
            distance = line.split(":")[-1]
            distance = float(distance)
            distance_avg += distance
            time_avg += time_curr
            if distance < distance_min:
                distance_min = distance
                time_min = time_curr
                line_min = line_curr
                iter_min = iter_curr
            if distance > distance_max:
                distance_max = distance
                time_max = time_curr
                line_max = line_curr
                iter_max = iter_curr

    distance_avg /= iter_curr
    time_avg /= iter_curr
    print(f"Out of {iter_curr} iterations")
    print(f"Minimum distance {distance_min} found in time {time_min} ms from the {iter_min}th iteration [line number: {line_min}]")
    print(f"Maximum distance {distance_max} found in time {time_max} ms from the {iter_max}th iteration [line number: {line_max}]")
    print(f"Average distance {distance_avg} and average time {time_avg} ms")

