import os, sys

#below are the arrays for the distance intervals
firstinterval = []
secondinterval = []
thirdinterval = []
fourthinterval = []
fifthinterval = []
sixthinterval = []
seventhinterval = []
eighthinterval = []
everythingelse = []




pdbsetoffiles = open("top2018_pdbs_with_disulfides.txt", "r")
#ssbondres1list = []
#ssbondres2list = []

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

def findCAcoords(pdbfile, ssbondres1list, ssbondres2list):
  ca_coords = {}
  currentfile= filesetter(pdbfile)
  for line in currentfile:
    if line.startswith("ATOM") or line.startswith("HETATM"):
      resname = line[17:20]
      if resname != "CYS":
        continue
      atomname = line[12:16]
      if atomname != " CA ":
        continue
      chainid = line[21]
      resseq = line[22:26]
      inscode = line[26]
      resid = chainid + resseq + inscode
      #print("****",resid,"****")
      if resid in ssbondres1list or resid in ssbondres2list:
        #print(resid)
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        altloc = line[16]
        if altloc == " ":
          ca_coords[resid] = [x,y,z, altloc]
        #ca_coords[resid] = {"x":x,"y":y,"z":z}
        else: 
          ca_coords[resid]= None
  return ca_coords
      
def printCAdistance(pdbid, ssbondres1list, ssbondres2list, ca_coords):
  i = 0
  for res1 in ssbondres1list:
    res2 = ssbondres2list[i]
    try:
      ca1 = ca_coords[res1]
    except KeyError:
      print(pdbid.strip(), ca_coords)
      continue
    try:
      ca2 = ca_coords[res2]
    except KeyError:
      print(pdbid.strip(), res2, "\n", ca_coords)
      continue 
    if ca1 is None or ca2 is None: 
      continue   
    distance = (((ca1[0]-ca2[0])**2) + ((ca1[1]-ca2[1])**2) + ((ca1[2]-ca2[2])**2))**0.5
    #print(pdbid.strip(),res1,res2,"%.3f" % distance)
    if distance > 0.001 and distance <= 0.99: 
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      firstinterval.append(outlist)
    elif distance > 0.99 and distance <= 1.99:
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      secondinterval.append(outlist)
    elif distance > 1.99 and distance <= 2.99: 
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      thirdinterval.append(outlist)
    elif distance >2.99 and distance <= 3.99:
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      fourthinterval.append(outlist)
    elif distance > 3.99 and distance <=4.99: 
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      fifthinterval.append(outlist)
    elif distance > 4.99 and distance <= 5.99: 
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      sixthinterval.append(outlist)
    elif distance > 5.99 and distance <= 6.99: 
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      seventhinterval.append(outlist)
    elif distance > 6.99 and distance <=7.99: 
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      eighthinterval.append(outlist)
    #Add another bin or two here,
    #trying to find a bin that is almost unpopulated
    #trying to find the cutoff of real possible lengths
    else:
      outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
      everythingelse.append(outlist)
    i += 1
"""
    outlist = [pdbid.strip(),res1, ca1[3], res2, ca2[3],"%.3f" % distance]
    print(",".join(outlist))
    """
    
"""
def gothroughCArecords(pdbfile):
  currentfile = filesetter(pdbfile)
  for fileline in currentfile: 
    if fileline.startswith("ATOM") and fileline[12:16] == " CA ": 
      #explicit parsing of the line helps with readability and debugging
      chainid = fileline[21]
      resseq = fileline[22:26]
      inscode = fileline[26]
      ### Debug parsing ###
      #print("|"+chainid+"|"+resseq+"|"+inscode+"|")
      #sys.exit()
      ### End debug ###
      cares_id = chainid + resseq + inscode
      #cares_id = fileline[22].strip() + fileline[23:26].strip() + fileline[27].strip()
      for element in ssbondres1list: 
        if cares_id == element:
          loop2file = filesetter(pdbfile)
          for fileline1 in loop2file: 
            if fileline1.startswith("ATOM") and fileline1[12:16] == " CA ":
              chainid = fileline[21]
              resseq = fileline[22:26]
              inscode = fileline[26]
              #Xcares_id = fileline1[21].strip() + fileline1[22:26].strip() + fileline1[26].strip()
              Xcares_id = chainid + resseq + inscode
             # print(ssbondres2list)
              if Xcares_id == ssbondres2list[ssbondres1list.index(element)]:
                x1 = float(fileline1[30:38].strip())
                y1 = float(fileline1[38:46].strip())
                z1 = float(fileline1[46:54].strip())
                x2 = float(fileline[30:38].strip())
                y2 = float(fileline[38:46].strip())
                z2 = float(fileline[46:54].strip())
               
                distance = (((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))**0.5
                print(distance, pdbfile)
                sys.exit()
                #add else statement
      #  elif cares_id in ssbondres2list: 
       #   for fileline2 in currentfile: 
        #    if fileline2.startswith("ATOM") and fileline2[12:16] == " CA ":
         #     Ycares_id = fileline1[21].strip() + fileline1[22:26].strip() + fileline1[26].strip()
          #    if Ycares_id == ssbondres2list(ssbondres2list.index(cares_id)): 
           #     distance = ((fileline1[31:38]-fileline[31:38])**2+(fileline1[31:38]-fileline[31:38])**2 + (fileline1[31:38]-fileline[31:38])**2)**0.5
"""
for pdbfile in pdbsetoffiles: 
  #print(line.strip(),file=sys.stderr)
  currentfile = filesetter(pdbfile)
  ssbondres1list = []
  ssbondres2list = []
  for fileline in currentfile: 
    if str(fileline[0:6]) == "SSBOND": #fileline.startswith("SSBOND")
       #12345678901234567890
       #SSBOND *** CYS A   23    CYS A   88                          1555   1555  2.04
      #TODO: add a check for non-CYS residues
      #Allow only CYS
      #skip others and do not add those bonds to the list
      #One time, make a list of all the different residue types involved in SSBOND records
      # give pdbid of anything not CYS
      chainid1 = fileline[15]
      resseq1 = fileline[17:21]
      inscode1 = fileline[21]    
          #res_id is chain, num, and icode
      chainid2 = fileline[29]
      resseq2 = fileline[31:35]
      inscode2 = fileline[35]

      res_id1 =  chainid1 + resseq1 + inscode1
      res_id2 = chainid2 + resseq2 + inscode2 
      ssbondres1list.append(res_id1)

      ssbondres2list.append(res_id2)
         # print(str(fileline[7:10].strip()) + "," + str(fileline[11:14].strip()) + "," + str(fileline[15].strip()) + "," + str(fileline[17:21].strip()) + str(fileline[22].strip()) + "," + str(fileline[25:28].strip()) + "," + str(fileline[29].strip()) + "," + str(fileline[31:35].strip()) + "," + str(fileline[36].strip()))
    
  #currentfile = filesetter(pdbfile)
  #currentfile.seek(0)
  #gothroughCArecords(pdbfile)
  ca_coords = findCAcoords(pdbfile, ssbondres1list, ssbondres2list)
  printCAdistance(pdbfile, ssbondres1list, ssbondres2list, ca_coords)

sortedfirstinterval = sorted(firstinterval)
sortedsecondinterval = sorted(secondinterval)
sortedthirdinterval = sorted(thirdinterval)
sortedfourthinterval = sorted(fourthinterval)
sortedfifthinterval = sorted(fifthinterval)
sortedsixthinterval = sorted(sixthinterval)
sortedseventhinterval = sorted(seventhinterval)
sortedeighthinterval = sorted(eighthinterval)
sortedeverythingelse = sorted(everythingelse)
print("CA-CA dist 1")
for i in sortedfirstinterval:
  print(i)
  #print(sortedfirstinterval[i])
print("CA-CA dist 2")
for i in sortedsecondinterval:
  print(i)
  #print(sortedsecondinterval[i])
print("CA-CA dist 3")
for i in sortedthirdinterval:
  #print(sortedthirdinterval[i])
  print(i)
print("CA-CA dist 4")
for i in sortedfourthinterval:
  #print(sortedfourthinterval[i])
  print(i)
print("CA-CA dist 5")
for i in sortedfifthinterval:
  #print(sortedfifthinterval[i])
  print(i)
print("CA-CA dist 6")
for i in sortedsixthinterval:
  #print(sortedsixthinterval[i])
  print(i)
print("CA-CA dist 7")
for i in sortedseventhinterval:
  #print(sortedseventhinterval[i])
  print(i)
print("CA-CA dist 8")
for i in sortedeighthinterval:
  #print(sortedeighthinterval[i])
  print(i)
print("CA-CA dist other")
for i in sortedeverythingelse:
  #print(sortedeverythingelse[i])
  print(i)

