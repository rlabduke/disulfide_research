#This is pseudocode and may not exactly work as written

def parseatom(line):
  #parse the line and return a dictionary of the contents
  #"chain"
  #"resseq"
  #"icode"

#disulfide bonds could be stored as a list of paired res_ids
#can look up stored CAs in their dictionary using these as keys
ssbonds = [[res1a, res2a], [res1b, res2b], [res1c, res2c]]

#store CA coordinates (or whatever) as a dictionary, keyed by res_id (chain+num+icode)
c_alphas = {}

for line in pdbfile:
  if line.startswith("ATOM"):
    atomline = parseatom(line)
    res_id = atomline["chain"]+":"+atomline["resseq"]+":"+atomline["icode"]
    #construct a unique residue identifier, this one will look something like:
    #" A:  14: "
    if res_id in [list,of,ssbond,ids]:
      #keep the line and do more with it
      #check if atomname is " CA "
      if atomname == " CA ":
        c_alphas[res_id] = atom_xyz
      #  or " CB ", later
