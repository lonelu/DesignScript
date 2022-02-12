import prody as pr, numpy as np
import freesasa

def calc_sasa(pdb_path, probe_radius, residues_list):
    ''' @probe_radius determines probe size used to calcualte sasa.
    probe_radius should be around 2 or 3 for design. however, for more realistic
    sasa (non-design purposes), select probe_radius of 1.4, which is approximately 
    the radius of a water molecule. 
        @residues_list should be a list of lists storing chains and resnums 
    of residues you want to probe, in this format: 
    [ ['A', 5,], ['A', 16'], ['A', 25] ]'''

    sasa_dict = {}

    # load PDB as prody object
    prody_obj = pr.parsePDB(pdb_path)

    # should remove solvent, hetatms, etc. when calculating sasa
    dry = prody_obj.select('protein and not element H D')

    # define atoms to calculate SASA over. for design, make sure you convert 
    # the residue to alanine 
    atom_sele = 'name N CA C O CB'
    
    # set up freesasa
    fs_results =  get_fs_result(dry, probe_radius)

    # iterate through residues to calculate SASA for
    for residue in residues_list: 

        # define residue chain, number, and atoms using prody syntax
        res_sele = 'chain {} and resnum {} and {}'.format(
            residue[0], residue[1], atom_sele)

        # store sasa in dict 
        sasa_dict['{}{}'.format(residue[0], residue[1])] = get_sasa(
            prody_obj, dry, res_sele, fs_results)

    return sasa_dict

def get_fs_result(dry, probe_radius):
    coords = list(x for y in dry.getCoords() for x in y)
    radii = list(freesasa.Classifier().radius(x, y) for x, y in zip(
        dry.getResnames(), dry.getNames()))
    return freesasa.calcCoord(coords, radii, freesasa.Parameters(
        {'probe-radius': probe_radius}))

def get_sasa(par_olig, dry, sel, fs_results):
    sel_inds = par_olig.select(sel).getIndices()
    
    # Map index of atom in selection to its position in @dry indices.
    dryindices = [x.getIndex() for x in dry]
    mapped_ind = np.where(np.in1d(dryindices, sel_inds))[0]
    
    sasa = sum(fs_results.atomArea(i) for i in mapped_ind)
    return sasa
