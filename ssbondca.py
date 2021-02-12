import os, sys

pdbsetoffiles = open("top2018_pdbs_with_disulfides.txt", "r")
ssbondres1list = []
ssbondres2list = []

def filesetter(pdbfile): 
  pdb = pdbfile[0:4]
  short = pdb[0:2]
  top_dir = "top2018_disulfide_pdbs"
  filepathlist = [top_dir,short,pdb,pdb+"_FH.pdb"]
  filepath = os.path.join(*filepathlist)
  #os.path.join is a safe way to join file paths that will work on any system.
  currentfile = open(filepath)
  return currentfile

def filesetter(pdbfile): 
  pdb = pdbfile[0:4]
  short = pdb[0:2]
  top_dir = "top2018_disulfide_pdbs"
  filepathlist = [top_dir,short,pdb,pdb+"_FH.pdb"]
  filepath = os.path.join(*filepathlist)
  currentfile = open(filepath)
  return currentfile



def gothroughCArecords(currentfile): 
  for fileline in currentfile: 
    if fileline.startswith("ATOM") and fileline[12:16] == " CA ": 
      #explicit parsing of the line helps with readability and debugging
      chainid = fileline[22]
      resseq = fileline[23:26]
      inscode = fileline[27]
      ### Debug parsing ###
      #print("|"+chainid+"|"+resseq+"|"+inscode+"|")
      #sys.exit()
      ### End debug ###
      cares_id = chainid + resseq + inscode
      #cares_id = fileline[22].strip() + fileline[23:26].strip() + fileline[27].strip()
      for element in ssbondres1list: 
        if cares_id == element: 
          for fileline1 in currentfile: 
            if fileline1.startswith("ATOM") and fileline1[12:16] == " CA ":
              Xcares_id = fileline1[22].strip() + fileline1[23:26].strip() + fileline1[27].strip()
              if Xcares_id == ssbondres2list.index(element): 
                distance = ((fileline1[31:38]-fileline[31:38])**2+(fileline1[31:38]-fileline[31:38])**2 + (fileline1[31:38]-fileline[31:38])**2)**0.5
                print(distance)
        elif cares_id in ssbondres2list: 
          for fileline2 in currentfile: 
            if fileline2.startswith("ATOM") and fileline2[12:16] == " CA ":
              Ycares_id = fileline1[22].strip() + fileline1[23:26].strip() + fileline1[27].strip()
              if Ycares_id == ssbondres2list(ssbondres2list.index(cares_id)): 
                distance = ((fileline1[31:38]-fileline[31:38])**2+(fileline1[31:38]-fileline[31:38])**2 + (fileline1[31:38]-fileline[31:38])**2)**0.5

for line in pdbsetoffiles: 
  print(line.strip(),file=sys.stderr)
  currentfile = filesetter(line)
  for fileline in currentfile: 
    if str(fileline[0:6]) == "SSBOND": #fileline.startswith("SSBOND")
       #12345678901234567890
       #SSBOND *** CYS A   23    CYS A   88                          1555   1555  2.04
          
          #res_id is chain, num, and icode
      res_id1 =  fileline[16].strip() + "," + fileline[17:21].strip() + "," + fileline[22].strip()
      res_id2 = fileline[30].strip() + "," + fileline[31:35].strip() + "," + fileline[36].strip()
      ssbondres1list.append(res_id1)
      ssbondres2list.append(res_id2)
         # print(str(fileline[7:10].strip()) + "," + str(fileline[11:14].strip()) + "," + str(fileline[15].strip()) + "," + str(fileline[17:21].strip()) + str(fileline[22].strip()) + "," + str(fileline[25:28].strip()) + "," + str(fileline[29].strip()) + "," + str(fileline[31:35].strip()) + "," + str(fileline[36].strip()))

  currentfile = filesetter(line)
  gothroughCArecords(currentfile)
