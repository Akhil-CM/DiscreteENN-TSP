

Data_File = "output.log"

if __name__ == "__main__":
    distance_min = 10e8
    line_num = 0
    iter_num = 0
    line_curr = 0
    iter_curr = 0
    with open(Data_File, "r") as input_file:
        for line in input_file:
            line = line.rstrip()
            if not line: continue
            line_curr += 1
            if "Total distance" not in line: continue
            iter_curr += 1
            distance = line.split(":")[-1]
            distance = float(distance)
            if distance < distance_min:
                distance_min = distance
                line_num = line_curr
                iter_num = iter_curr

    print(f"Minimum distance {distance_min} at line number {line_num} from the {iter_num} iteration")

