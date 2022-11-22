import sys

def extract(input_file):
    f = open(input_file)
    content = f.read()
    f.close()

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
        # print(splitted)
        sgh = splitted[4]
        model = splitted[8][:-1]
        res += "\n{}, {}, {}, {}".format(auc, ci, sgh, model)
        
    return res

if __name__ == "__main__":
    input_file, output_file = sys.argv[1], sys.argv[2]
    res = extract(input_file)
    print(res)
    out = open(output_file, "w")
    out.write(res)
    out.close()