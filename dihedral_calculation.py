#get ssbond proxies, as shown in ssbond_proxies_from_cctbx.py
#fetch the SG atom
from __future__ import absolute_import, division, print_function
import sys
import os 

pdbfiles = open("top2018_pdbs_with_disulfides.txt", "r")
outputlist = []
for line in pdbfiles: 
    pdbfile = line[0:4]
    outputlist.append(pdbfile)

def create_model_mmtbx(path_to_pdb_file):
  import mmtbx.model
  import iotbx.pdb
  from libtbx.utils import null_out
  pdb_inp = iotbx.pdb.input(file_name = path_to_pdb_file)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    pdb_interpretation_params = params,
    build_grm   = True,
    stop_for_unknowns = False,
    log         = null_out())
  return model

def create_hierarchy_from_model(model):
  return model.get_hierarchy()

def get_ssbond_proxies(grm):
  from cctbx.geometry_restraints.linking_class import linking_class
  origin_ids = linking_class()
  specific_origin_id = origin_ids.get_origin_id('SS BOND')
  pair_proxies = grm.pair_proxies()
  return(pair_proxies.bond_proxies.simple.proxy_select(origin_id=specific_origin_id))
  #origin id selects for a particular kind of bond in this case ssbond
  #see modules/cctbx_project/cctbx/geometry_restraints/auto_linking_types.py for supported interaction types
  #see modules/cctbx_project/cctbx/geometry_restraints/manager.py for additional sample code

def find_cys_atoms(sg_atom):
  atom_group = sg_atom.parent()
  ca = None
  cb = None
  n  = None
  #already have sg
  for atom in atom_group.atoms():
    if atom.name == " CA ":
      ca = atom
    elif atom.name == " CB ":
      cb = atom 
    elif atom.name == " N  ":
      n = atom
    #elif ()find other atoms
  return ca,cb,n

def calculate_dihedral(atom1xyz, atom2xyz, atom3xyz, atom4xyz):
  #this is expecting atom.xyx as arguments
  from scitbx.matrix import dihedral_angle
  d = dihedral_angle(sites=[atom1xyz, atom2xyz, atom3xyz, atom4xyz], deg=True)
  return d

def get_residue_identity_from_atom(atom):
  ag = atom.parent()
  rg = ag.parent()
  chain = rg.parent()
  resname = ag.resname
  resseq = rg.resseq
  icode = rg.icode
  altloc = ag.altloc
  chainid = chain.id
  return ":".join([chainid,resname,altloc,resseq,icode])

#calculate dihedrals for:
#N1-CA1-CB1-SG1    (X1)
#CA1-CB1-SG1--SG2  (X2)
#CB1-SG1--SG2-CB2  (X3)
#SG1--SG2-CB2-CA2  (X2')
#SG2-CB2-CA2-N2    (X1')
#these are the 5 dihedrals for a disulfide

#print each of these out with sufficient other information (PDBid, residue identifiers) that we can look it up

#Probably format as one line per disulfide, listing each dihedral as its own column

#at some point, we will add the CA-CA, CB-CB, and SG--SG distances

try:
  from scitbx.matrix import dihedral_angle
except ImportError:
  sys.stderr.write("CCTBX environment not sourced.\n")
  sys.stderr.write("source ~/MolProbity/build/setpaths.sh\n")
  sys.exit()
pdbdir = '/Users/sushritp1/Documents/GitHub/disulfide_research/top2018_disulfide_pdbs/'
for element in outputlist:
  path_to_pdb_file = pdbdir + element[0:2] + "/" + element + "/" + element + "_FH.pdb"
  if not os.path.isfile(path_to_pdb_file):
    print("Could not find file" + path_to_pdb_file,file=sys.stderr)
    continue
  try: 
   
    model = create_model_mmtbx(path_to_pdb_file)
  
    hierarchy = create_hierarchy_from_model(model)
    all_atoms=hierarchy.atoms()
    grm = model.get_restraints_manager().geometry
   
    test = get_ssbond_proxies(grm)
   
    for proxy in test:
      
      sg1 = all_atoms[proxy.i_seqs[0]]
      sg2 = all_atoms[proxy.i_seqs[1]]
      ca1,cb1,n1=find_cys_atoms(sg1)
      ca2,cb2,n2=find_cys_atoms(sg2)
      if None in [n1,ca1,cb1,sg1]:
        x1 = 9999
      else:
        x1=calculate_dihedral(n1.xyz, ca1.xyz, cb1.xyz, sg1.xyz)
      
      if None in [ca1,cb1,sg1,sg2]:
        x2 = 9999
      else:
        x2=calculate_dihedral(ca1.xyz,cb1.xyz , sg1.xyz, sg2.xyz)
      
      if None in [cb1,sg1,sg2,cb2]:
        x3 = 9999
      else:
        x3=calculate_dihedral(cb1.xyz, sg1.xyz, sg2.xyz, cb2.xyz)
  
      if None in [sg1,sg2,cb2,ca2]:
        x2prime = 9999
      else:
        x2prime=calculate_dihedral(sg1.xyz, sg2.xyz, cb2.xyz, ca2.xyz)
      
      if None in [sg2,cb2,ca2,n2]:
        x1prime = 9999
      else:
        x1prime=calculate_dihedral(sg2.xyz, cb2.xyz, ca2.xyz, n2.xyz)
      
      print("%s,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f" % (element,get_residue_identity_from_atom(sg1),get_residue_identity_from_atom(sg2),x1,x2,x3,x2prime,x1prime))
      #N1-CA1-CB1-SG1    (X1)
#CA1-CB1-SG1--SG2  (X2)
#CB1-SG1--SG2-CB2  (X3)
#SG1--SG2-CB2-CA2  (X2')
#SG2-CB2-CA2-N2    (X1')
  except ImportError:
    print("Source the phenix environment")
    print("source ~/MolProbity/build/setpaths.sh")
    #sys.exit(0)
   #print("Model creation failed: " + path_to_pdb_file,file=sys.stderr)
    continue



"""def get_ssbond_proxies_asu(grm):
  #*a*s*ymmetric *u*nit bond proxies. Should find bonds between symmetry mates, etc.
  #I don't have this quite working yet
  from cctbx.geometry_restraints.linking_class import linking_class
  origin_ids = linking_class()
  specific_origin_id = origin_ids.get_origin_id('SS BOND')

  pair_proxies = grm.pair_proxies()
  return(pair_proxies.bond_proxies.asu.proxy_select(origin_id=specific_origin_id))"""




"""
def calculate_dihedral(atom1xyz, atom2xyz, atom3xyz, atom4xyz):
  #this is expecting atom.xyx as arguments
  from scitbx.matrix import dihedral_angle
  d = dihedral_angle(sites=[atom1xyz, atom2xyz, atom3xyz, atom4xyz], deg=True)
  return d
"""
"""
"""
"""
x1 = calculate_dihedral(1,2,3,4)
x2 = calculate_dihedral(1,2,3,4)
x3 = calculate_dihedral(1,2,3,4)
x2prime = calculate_dihedral(1,2,3,4)
x1prime = calculate_dihedral(1,2,3,4)

print(element + " " + x1)
print(element + " " + x2)
print(element + " " + x3)
print(element + " " + x4)
print(element + " " + x5)
"""
