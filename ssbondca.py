import os, sys

pdbsetoffiles = open("top2018_pdbs_with_disulfides.txt", "r"); 
ssbondres1list = []
ssbondres2list = []
c_alphas = {}
sscount = 0
c_alpharesidlist = []
matchcount = -1 
def filesetter(pdbfile): 
  pdb = pdbfile[0:4]
  short = pdb[0:2]
  top_dir = "top2018_disulfide_pdbs"
  filepathlist = [top_dir,short,pdb,pdb+"_FH.pdb"]
  filepath = os.path.join(*filepathlist)
  #os.path.join is a safe way to join file paths that will work on any system.
  currentfile = open(filepath)
  return currentfile

def parseatom(fileline): 
    res_id = fileline[22].strip() + ':' + fileline[23:26] + ':' + fileline[27].strip() 
    c_alpharesidlist.append(res_id)
    c_alphas[res_id] = {"X": float(fileline[31:38].strip()), "Y": float(fileline[39:46].strip()), "Z": float(fileline[47:54].strip())}
    return res_id

def filesetter(pdbfile): 
  pdb = pdbfile[0:4]
  short = pdb[0:2]
  top_dir = "top2018_disulfide_pdbs"
  filepathlist = [top_dir,short,pdb,pdb+"_FH.pdb"]
  filepath = os.path.join(*filepathlist)
  #os.path.join is a safe way to join file paths that will work on any system.
  currentfile = open(filepath)
  return currentfile

def ssbondreslists(fileline): 
    ssbondres1list.append(fileline[16].strip() + ':'+ fileline[18:21] + ':' + fileline[22].strip())
    ssbondres2list.append(fileline[16].strip() + ':'+ fileline[32:35] + ':' + fileline[22].strip())


for line in pdbsetoffiles: 
   print(line.strip(),file=sys.stderr)
   currentfile = filesetter(line)
   for fileline in currentfile:
       if fileline.startswith("TITLE     "): 
            print("The name of the structure that the following information will be under is this: " + fileline.strip())
       if str(fileline[0:6]) == "SSBOND": #fileline.startswith("SSBOND")
       #12345678901234567890
       #SSBOND *** CYS A   23    CYS A   88                          1555   1555  2.04
          
          ssbondreslists(fileline)
          print(str(fileline[7:10].strip()) + "," + str(fileline[11:14].strip()) + "," + str(fileline[15].strip()) + "," + str(fileline[17:21].strip()) + str(fileline[22].strip()) + "," + str(fileline[25:28].strip()) + "," + str(fileline[29].strip()) + "," + str(fileline[31:35].strip()) + "," + str(fileline[36].strip()))
       
 #for fileline in currentfile:

        #this is where I am calculating the distances     
  #      if fileline.startswith("ATOM"): 
   #         res_id = parseatom(fileline)
    #        if fileline[13:16] == " CA ": 
     #                 for x in ssbondres1list:
      #                    for y in c_alpharesidlist: 
       #                       if x == y: 
        #                           for z in ssbondres2list: 
         #                              if y == z: 
          #                                  distance = abs(((c_alphas.y("X")) - (c_alphas.z("X"))))
                         #d = ((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5

