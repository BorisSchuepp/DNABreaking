

if __name__ == '__main__':
    fig_num = "57"
    file = open(f"C:/Users/Boris/Desktop/Data/Supplement/NEW/FigS{fig_num}/FigS{fig_num}.csv", "r")
    out_file = open(f"C:/Users/Boris/Desktop/Data/Supplement/NEW/FigS{fig_num}/FigS{fig_num}R.csv", "w")
    start = True
    good_indices = []
    names = []
    for line in file:
        line = line.strip("\n")
        line = line.split(",")
        if start:
            for i in range(0, len(line)):
                if line[i].find("BOT") != -1:
                    good_indices.append(i)
                    names.append(line[i].split("#")[0].split("_")[0] + " Sample " + line[i].split("#")[0].split("_")[1])
            start = False
            out_file.write(line[0] + "," + ",".join(names) + "\n")
        else:
            out_file.write(line[0] + "," + ",".join([line[j] for j in good_indices]) + "\n")