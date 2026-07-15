import biotite
import biotite.structure as struc
from StrucTools import *

def extract_binder_contig(contigs_str: str):
    """ Extract Binder contig from contigs string to know where the motifs are for given design

        Parameters
        ----------
        contigs_str: str
            Contigs string of the overall design (includes target and binder)

        Returns
        -------
        contigs_binder: str
            Contigs string of the binder
    """
    contigs = contigs_str.split(';')
    # Store in a list initially if having multiple targets and binders
    contigs_target = []
    contigs_binder = []
    for contig in contigs:
        if "/" in contig:
            contigs_binder.append(contig)
        else:
            contigs_target.append(contig)
    # Check first if there are multiple binders
    if len(contigs_binder) > 1:
        contigs_binder = ';'.join(contigs_binder)
    else:
        contigs_binder = contigs_binder[0]
    return contigs_binder

def extract_input_motif_coords(path_input_pdb: str, contig_binder: str, desired_motif_chain_id: str = ""):
    """ Extract coordinates of the motif atoms from the input PDB providing the motif contig

        Parameters
        ----------
        path_input_pdb: str
            Path to the input PDB file
        motif_contig: str
            Contig of the motif in the input PDB file. Example: B1-21/20-20/C1-21/20-20/D1-21
        motif_chain_id: str
            Chain ID of the specific motif desired for alignment and rmsd calculation. By default, since "" will extract all motif atoms
        Returns
            all_coords_motif_input: np.array
            Coordinates of the motif atoms
    """
    # 1. Load in the pdb to get atom_array
    atom_array_input = extract_atom_array(path_input_pdb)
    all_coords_motif_input = []
    # 2. Iterate through all the contigs within the binder, may specify motif coordinates or regions to be designed
    for motif_contig in contig_binder.split('/'):
        
        # 2.1 Handle case where desired_motif_chain_id is not provided so goal is to extract all coordinates of the motif atoms from the design PDB
        if desired_motif_chain_id == "":
            motif_extraction_condition = (motif_contig[0].isalpha())
        else:
            motif_extraction_condition = (motif_contig[0] == desired_motif_chain_id)
        
        # 2.2 Now iterate through with the motif_extraction_condition
        if motif_extraction_condition:
            motif_chain_start, motif_chain_end = motif_contig.split('-')
            motif_chain_id = motif_chain_start[0]
            motif_res_start, motif_res_end = int(motif_chain_start[1:]), int(motif_chain_end)
    
            # Create input atom mask of motif atoms in the input
            atom_array_input_motif = atom_array_input[((atom_array_input.chain_id == motif_chain_id) & (atom_array_input.res_id >= motif_res_start) & (atom_array_input.res_id <= motif_res_end))]
            coords_motif_input = atom_array_input_motif.coord
            all_coords_motif_input.append(coords_motif_input)
        else:
            continue
    all_coords_motif_input = np.concatenate(all_coords_motif_input, axis = 0)
    return all_coords_motif_input

def extract_design_motif_coords(path_pdb_design: str, contig_binder: str, desired_motif_chain_id: str = "", binder_chain_id: str = "B", return_coords_mask: bool = False):
    """ Extract coordinates of the motif atoms from the design PDB providing the motif contig
        By default: Extracts all coordinates of the motif atoms from the design PDB unless a specific motif_chain_id is provided

        Parameters
        ----------
        contig_binder: str
            Contig of the motif in the input PDB file. Example: B1-21/20-20/C1-21/20-20/D1-21
        path_pdb_design: str
            Path to the design PDB file
        desired_motif_chain_id: str
            Chain ID of the specific motif desired for alignment and rmsd calculation. By default, since "" will extract all motif atoms
        binder_chain_id: str
            Chain ID of the binder in the design PDB. By default, "B" since typically "A" corresponds to target
    
    """
    # 1. Create residue-specific mask of the motif residues in the design binder sequence
    res_motif_mask = []
    for region in contig_binder.split('/'):
        
        # 2. Handle case where motif_chain_id is not provided so goal is to extract all coordinates of the motif atoms from the design PDB
        if desired_motif_chain_id == "":
            motif_extraction_condition = (region[0].isalpha())
        else:
            motif_extraction_condition = (region[0] == desired_motif_chain_id)
        
        # 3. Handle case where motif_chain_id is provided so goal is to extract coordinates of the motif atoms from the design PDB
        if motif_extraction_condition:
            #print("Region: ", region)
            preserved_region = region[1:]
            start, end = preserved_region.split('-')
            motif_len = int(end) - int(start) + 1
            motif_mask = [True] * motif_len
        else:
            design_len1, design_len2 = region.split('-')
            # Handle case of chains corresponding to those preserved from initial input
            if design_len1 != design_len2:
                design_strand_size = int(design_len2) - int(design_len1[1:]) + 1
            # Handle case of chains corresponding to those designed via RFDiffusion
            elif design_len1 == design_len2:
                design_strand_size = design_len2
            motif_mask = [False] * int(design_strand_size)
        
        res_motif_mask += motif_mask
    #print(len(res_motif_mask))
    
    atom_array_design_complex = extract_atom_array(path_pdb_design)
    atom_array_design_binder = atom_array_design_complex[atom_array_design_complex.chain_id == binder_chain_id]
    mask_atom_motif_design = struc.spread_residue_wise(atom_array_design_binder, res_motif_mask)
    atom_array_design_motif = atom_array_design_binder[mask_atom_motif_design]
    coords_motif_design = atom_array_design_motif.coord
    # Condition to return mask if desired for alignment calculations in other settings
    if return_coords_mask:
        return coords_motif_design, res_motif_mask # Since alignment occurs on alpha-carbons, residue mask = alpha-carbon mask
    else:
        return coords_motif_design