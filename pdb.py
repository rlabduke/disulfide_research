import os, sys

pdbsetoffiles = open("top2018_pdbs_with_disulfides.txt", "r"); 
def checkdiffchains(disulfide): 
    if disulfide[15] != disulfide[29]: 
        diffchainsarray.append("This SSBOND has disulfides on different chains: " + str(line.strip()) + " " + str(disulfide[1:6].strip()) + " " + str(disulfide[16].strip()) + " " + str(disulfide[18:21].strip()) + " " + str(disulfide[30].strip()) + " " + str(disulfide[32:35].strip()))
        # above is file name, chain name, residues name

def checkdistances(disulfide): 
    #print ("|" + disulfide[17:20] + "|")
    #print ("|" + disulfide[31:34] + "|")
    distance = abs(int(disulfide[17:21].strip())-int(disulfide[31:35].strip()))
    if distance == 1 and disulfide[15] == disulfide[29]: 
        distancesarray1.append("From this file: " + str(line.strip()) + ", This SSBOND " + str(disulfide[11:36].strip()) + " has a  sequence distance of 1")
    if distance == 2 and disulfide[15] == disulfide[29]: 
        distancesarray2.append("From this file: " + str(line.strip()) + ", This SSBOND " + str(disulfide[11:36].strip()) + " has a  sequence distance of 2")
    if distance == 3 and disulfide[15] == disulfide[29]: 
        distancesarray3.append("From this file: " + str(line.strip()) + ", This SSBOND " + str(disulfide[11:36].strip()) + " has a  sequence distance of 3")
    if distance == 4 and disulfide[15] == disulfide[29]: 
        distancesarray4.append("From this file: " + str(line.strip()) + ", This SSBOND " + str(disulfide[11:36].strip()) + " has a  sequence distance of 4")
    if distance == 5 and disulfide[15] == disulfide[29]:
        distancesarray5.append("From this file: " + str(line.strip()) + ", This SSBOND " + str(disulfide[11:36].strip()) + " has a  sequence distance of 5")
    if distance == 6 and disulfide[15] == disulfide[29]: 
        distancesarray6.append("From this file: " + str(line.strip()) + ", This SSBOND " + str(disulfide[11:36].strip()) + " has a  sequence distance of 6")

#for distance array put pdb, chain, both res numbers, insertion code
#distances only applies for things with same chain

    # col 1-6 is Record Name, Field: SSBOND
    #col 8-10 is serial Number
    #col 12-14 is residue name
    #col 16 is chain identifier
    # col 18-21 is residue sequence number
    # col 22 is insertion code 
    # col 26-28 is residue name 
    #col 30 is chain identifier
    # col 32-35 is residue sequence number 
    

def filesetter(pdbfile): 
  pdb = pdbfile[0:4]
  short = pdb[0:2]
  top_dir = "top2018_disulfide_pdbs"
  filepathlist = [top_dir,short,pdb,pdb+"_FH.pdb"]
  filepath = os.path.join(*filepathlist)
  #os.path.join is a safe way to join file paths that will work on any system.
  currentfile = open(filepath)
  return currentfile

diffchainsarray = [] 
distancesarray1 = [] 
distancesarray2 = []
distancesarray3 = [] 
distancesarray4 = []
distancesarray5 = []
distancesarray6 = []

for line in pdbsetoffiles: 
   print(line.strip(),file=sys.stderr)
   currentfile = filesetter(line)
   for fileline in currentfile: 
       if fileline.startswith("TITLE     "):
           print("The name of the structure that the following information will be under is this: " + fileline.strip())
       if str(fileline[0:6]) == "SSBOND":
          #  print(fileline)
            checkdiffchains(fileline)
            checkdistances(fileline)
    

for element in diffchainsarray: 
    print(element)

for element in distancesarray1: 
    print (element)

for element in distancesarray2: 
    print (element)   

for element in distancesarray3: 
    print (element)

for element in distancesarray4: 
    print (element)

for element in distancesarray5: 
    print (element)

for element in distancesarray6: 
    print (element)





