f_path = "7l6v-A-B.pdb"

with open(f_path, "r") as pf:
    lines = pf.readlines()
with open(f_path[:-4] + "-unbound.pdb", "w") as pf:
    for line in lines:
        if line.startswith("ATOM  "):
            if line[21] == "A":
                new_x = str(round(float(line[31:38]) + 122, 3))
                if len(new_x) == 5:
                    new_x += "0"
                elif len(new_x) == 4:
                    new_x += "00"
                elif len(new_x) == 3:
                    new_x += "000"
                elif len(new_x) == 2:
                    new_x += ".000"
                pf.write(line[:31] + " " + new_x + line[38:])
            else:
                pf.write(line)
