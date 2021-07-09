from __future__ import absolute_import, division, print_function
import sys, os

def create_model_from_file(path_to_pdb_file):
  from iotbx.data_manager import DataManager    # Load in the DataManager
  dm = DataManager()             # Initialize the DataManager and call it dm
  model = dm.get_model(path_to_pdb_file)
  return model

def create_hierarchy_from_model(model):
  return model.get_hierarchy()

def calculate_dihedral(atom1xyz, atom2xyz, atom3xyz, atom4xyz):
  from scitbx.matrix import dihedral_angle
  d = dihedral_angle(sites=[atom1xyz, atom2xyz, atom3xyz, atom4xyz], deg=True)
  return d


pdbfiles = open("top2018_pdbs_with_disulfides.txt", "r")
print(pdbfiles)
outputlist = []
for line in pdbfiles: 
    pdbfile = line[0:4]
    outputlist.append(pdbfile)



pdbdir = '/Users/sushritp1/Documents/GitHub/disulfide_research/top2018_disulfide_pdbs/'
for element in outputlist:
  
  path_to_pdb_file = pdbdir + element[0:2] + "/" + element + "/" + element + "_FH.pdb"
  if not os.path.isfile(path_to_pdb_file):
    print("Could not find file" + path_to_pdb_file,file=sys.stderr)
    continue
  try: 

    model = create_model_from_file(path_to_pdb_file)
    hierarchy = create_hierarchy_from_model(model)
  except:
    print("Model creation failed: " + path_to_pdb_file,file=sys.stderr)
    continue


  for chain in hierarchy.chains(): 
    for rg in chain.residue_groups():
        for ag in rg.atom_groups():
            if ag.resname == "CYS":
                for atom in ag.atoms():
                    if atom.name == " SG ":
                       # print(os.path.basename(path_to_pdb_file), chain.id, rg.resseq+rg.icode,ag.resname, atom.name, atom.b, atom.occ)
                        print(",".join([os.path.basename(path_to_pdb_file),   atom.name, ag.altloc, ag.resname, chain.id, rg.resseq+rg.icode, str(atom.b), str(atom.occ)]))