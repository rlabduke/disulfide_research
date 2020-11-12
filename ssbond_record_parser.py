import sys

#Official formatting from https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
#         1         2          3        4         5         6         7         8
#12345678901234567890123456789012345678901234567890123456789012345678901234567890
#SSBOND   1 CYS A    6    CYS A  127                          1555   1555  2.03 
#SSBOND   2 CYS A   30    CYS A  115                          1555   1555  2.07 
#SSBOND   3 CYS A   64    CYS A   80                          1555   1555  2.06 
#SSBOND   4 CYS A   76    CYS A   94                          1555   1555  2.04 

#COLUMNS        DATA  TYPE     FIELD            DEFINITION
#--------------------------------------------------------------------------------
# 1 -  6        Record name    "SSBOND"
# 8 - 10        Integer        serNum           Serial number.
#12 - 14        LString(3)     "CYS"            Residue name.
#16             Character      chainID1         Chain identifier.
#18 - 21        Integer        seqNum1          Residue sequence number.
#22             AChar          icode1           Insertion code.
#26 - 28        LString(3)     "CYS"            Residue name.
#30             Character      chainID2         Chain identifier.
#32 - 35        Integer        seqNum2          Residue sequence number.
#36             AChar          icode2           Insertion code.
#60 - 65        SymOP          sym1             Symmetry operator for residue 1.
#67 - 72        SymOP          sym2             Symmetry operator for residue 2.
#74 - 78        Real(5.2)      Length           Disulfide bond distance

class ssbond_record():
  #This class accepts an SSBOND record line from a pdb and stores the important information
  def __init__(self, line):
    #first, do a safety check that the input is an SSBOND record
    if not line.startswith("SSBOND"):
      sys.stderr.write("This line was incorrectly passed as an SSBOND record:\n")
      sys.stderr.write(line)
      return None
    #Parse the SSBOND record and store the contents
    #self.record = line[0:6]
    self.serNum = line[7:10]
    self.chainid1 = line[14:16].strip()
      #the official format has 1 char for chain ids, but some large files use 2 chars
    self.resseq1 = line[17:21]
      #I've kept residue number as 4 white-space padded characters
      #This format is very useful for printing and easier to convert out of than back into
    self.icode1 = line[21:22]

    self.chainid2 = line[29:31].strip()
    self.resseq2 = line[31:35]
    self.icode2 = line[35:36]

    self.ss_dist = float(line[73:78].strip())

if __name__ == "__main__":
  #This file is meant to be a library, imported for use by other scripts but not run on its own
  #However, it can be useful to run it from commandline for testing purposes
  #'if __name__ == "__main__":' checks to see if this script is being run on commandline
  pdbfile = open(sys.argv[1])
  for line in pdbfile:
    if not line.startswith("SSBOND"):
      continue
    ss = ssbond_record(line)
    sys.stdout.write(line)
    sys.stdout.write("serNum |%s|, chain1 |%s|, num1 |%s|, icode1 |%s|, chain2 |%s|, num2 |%s|, icode2 |%s|, ss_dist |%.2f|\n\n"
      % (ss.serNum, ss.chainid1, ss.resseq1, ss.icode1, ss.chainid2, ss.resseq2, ss.icode2, ss.ss_dist))
      # "%s" and "%.2f" are what's called 'string formatting', a way to pass an ordered list of variables into a string

    
