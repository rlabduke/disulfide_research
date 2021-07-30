from __future__ import absolute_import, division, print_function
import sys

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
  #origin id selects for a particular kind of bond in this case ssbond
  #see modules/cctbx_project/cctbx/geometry_restraints/auto_linking_types.py for supported interaction types
  #see modules/cctbx_project/cctbx/geometry_restraints/manager.py for additional sample code

  pair_proxies = grm.pair_proxies()
  return(pair_proxies.bond_proxies.simple.proxy_select(origin_id=specific_origin_id))

def get_ssbond_proxies_asu(grm):
  #*a*s*ymmetric *u*nit bond proxies. Should find bonds between symmetry mates, etc.
  #I don't have this quite working yet
  from cctbx.geometry_restraints.linking_class import linking_class
  origin_ids = linking_class()
  specific_origin_id = origin_ids.get_origin_id('SS BOND')

  pair_proxies = grm.pair_proxies()
  return(pair_proxies.bond_proxies.asu.proxy_select(origin_id=specific_origin_id))

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

def get_bfactor_from_atom(atom):
  return atom.b

def get_occupancy_from_atom(atom):
  return atom.occ

def get_atom_coordinates(atom):
  return {"x":atom.xyz[0],
          "y":atom.xyz[1],
          "z":atom.xyz[2]}

#accept a pdb file as a commandline argument
path_to_pdb_file = sys.argv[1]
#get a high level model object from which other objects are derived
#get a hierachy object for most of our needs
#get a geometry restraints manager (grm for short)
#  The grm is aware of what atoms are bonded to each other
model = create_model_mmtbx(path_to_pdb_file)
hierarchy = create_hierarchy_from_model(model)
grm = model.get_restraints_manager().geometry


ssbond_proxies_simple = get_ssbond_proxies(grm)
###printing for proxy contents
#for simple_proxy in ssbond_proxies_simple:
#  for thing in dir(simple_proxy):
#    print(thing)
#  sys.exit()
ssbond_proxies_asu = get_ssbond_proxies_asu(grm)
###printing for proxy contents
#for asu_proxy in ssbond_proxies_asu:
#  for thing in dir(asu_proxy):
#    print(thing)
#  sys.exit()

#bond proxies are objects that represent covalent bonds (or other) in a structure
#For us, the useful part is bond_proxy.i_seqs, which will let us look up the atoms in each bond
#  proxy.i_seqs = (89, 2027), for example
#bond_proxy.i_seqs is a tuple containing two integer indices.  The indices map to a list of atoms in the hierarchy.atoms()


all_atoms = hierarchy.atoms() #a list of all atoms in the hierarchy, its indices are the proxy i_seqs

print("SIMPLE")
for proxy in ssbond_proxies_simple:
  atom1 = all_atoms[proxy.i_seqs[0]]
  atom2 = all_atoms[proxy.i_seqs[1]]
  print(get_residue_identity_from_atom(atom1), atom1.name, atom1.xyz)
  print(get_residue_identity_from_atom(atom2), atom2.name, atom2.xyz)
  print("--------------------------------------------------------------------------------")

### asu proxies (that's *a*s*ymmetric *u*nit) find ssbonds across crystal contacts or symmetry mates
###this is where the 0-length and 10+length ssbonds are proberly recorded
###however, looking them up by i_seqs doesn't work
###We'll revisit this once I learn a bit more
#print("\nASU")
#for proxy in ssbond_proxies_asu:
#  atom1 = all_atoms[proxy.i_seqs[0]]
#  atom2 = all_atoms[proxy.i_seqs[1]]
#  print(get_residue_identity_from_atom(atom1), atom1.name, atom1.xyz)
#  print(get_residue_identity_from_atom(atom2), atom2.name, atom2.xyz)
#  print("--------------------------------------------------------------------------------")
