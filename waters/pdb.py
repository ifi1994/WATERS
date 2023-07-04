from Bio.PDB import PDBParser, Structure

file1 = open( "/Users/pantelispanka/Downloads/conf-protein.pdb", 'r')
Lines = file1.readlines()

count = 0
# Strips the newline character

# model = None
# new = True
# i = 0
# for line in Lines:
#     if new:
#         file = open("./" + str(i) + "_a.pdb", "w")
#     file.write(line)
#     new = False
#     if "ENDMDL" in line:
#         model = None
#         i += 1
#         new = True
    # print("Line{}: {}".format(count, line.strip()))




sloppyparser = PDBParser()
structure = sloppyparser.get_structure("0", "/Users/pantelispanka/Jaqpot/WATERS/waters/1_a.pdb")



# structure = sloppyparser.get_structure("0", "/Users/pantelispanka/Downloads/conf-protein.pdb")
#
#
for model in structure.get_models():
    print(model)
