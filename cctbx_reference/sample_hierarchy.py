from __future__ import absolute_import, division, print_function
import sys

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

path_to_pdb_file = sys.argv[1]

model = create_model_from_file(path_to_pdb_file)
hierarchy = create_hierarchy_from_model(model)

#very simple walk over the hierarchy and its contents
for chain in hierarchy.chains():
  for rg in chain.residue_groups():
    for ag in rg.atom_groups(): #alternate conformations are split into separate atom groups
      for atom in ag.atoms():
        if atom.name == " CA ":
          print(ag.resname, rg.resseq+rg.icode, atom.name, atom.xyz)
