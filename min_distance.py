

Data_File = "output.log"

if __name__ == "__main__":
    distance_min = 10e8
    time_final = -1.0
    time_curr = -1.0
    line_num = 0
    iter_num = 0
    line_curr = 0
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
            if distance < distance_min:
                distance_min = distance
                time_final = time_curr
                line_num = line_curr
                iter_num = iter_curr

    print(f"Minimum distance {distance_min} found in time {time_final} ms at line number {line_num} from the {iter_num} iteration")

