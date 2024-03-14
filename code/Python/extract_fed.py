import sys

def extract(input_file):
    with open(input_file) as f:
        content = f.read()
    lines1 = [line for line in content.split('\n') if "AUC:  " in line]
    lines2 = [line for line in content.split('\n') if "generating auc" in line]
    res = "AUC, CI, Site, Model"
    for i in range(len(lines1)):
        line = lines1[i]
        splitted = line.split(" ")
        auc = splitted[2]
        ci = splitted[-2]

        line = lines2[i]
        splitted = line.split(" ")
        site = splitted[4]
        model = splitted[8][:-1]
        res += "\n{}, {}, {}, {}".format(auc, ci, site, model)
    return res

if __name__ == "__main__":
    input_file, output_file = sys.argv[1], sys.argv[2]
    with open(output_file, "w") as out:
        res = extract(input_file)
        out.write(res)
