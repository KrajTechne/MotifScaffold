import sys
sys.path.append("/Workspace/Users/karthik.raj@bio-techne.com/ScoringPhysics")
from StrucTools import *
import biotite
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import random

class MotifSampler:
    def __init__(self, path_input_structure: str, target_chain_id: str, binder_chain_dict: dict, mergable_tolerance: int = 4):
        self.path_input_structure = path_input_structure
        self.path_clean_structure = self.clean_and_save_structure()
        self.target_chain_id = target_chain_id
        self.binder_chain_dict = binder_chain_dict
        self.mergable_tolerance = mergable_tolerance
    
    def clean_and_save_structure(self):
        """ Clean the input structure and save it to a file with _clean appended to the filename """
        
        atom_array_clean = clean_structure(self.path_input_structure)
        
        # Save the cleaned structure to same location as original, but with _clean appended to the filename
        file = pdb.PDBFile()
        file.set_structure(atom_array_clean)

        # Create the output file path
        filename, filename_extension = self.path_input_structure.split(".")
        output_file_path = filename + "_clean." + filename_extension
        print(f"Saving cleaned structure to {output_file_path}")

        # Write the cleaned structure to the output_file_path location
        file.write(output_file_path)
        return output_file_path
    
    def motif_indices_to_rfdcontig(self, paratope_indices: list[int], binder_chain_id: str, mergable_tolerance: int = 4) -> str:
        """ For a respective chain, create a RFDiffusion contig string that specifies preserved residues and linker lengths between motifs. """
        
        # 0. Handle edge cases (empty lists)
        if not paratope_indices:
            return ""

        # Ensure indices are sorted to prevent negative differences/logic errors
        sorted_indices = sorted(paratope_indices)
    
        # 1. Group indices based on tolerance
        groups = []
        current_group = [sorted_indices[0]]

        for pos in sorted_indices[1:]:
            # If the gap between the current residue and the last residue in our group is within tolerance
            if (pos - current_group[-1]) <= mergable_tolerance:
                current_group.append(pos)
            else:
                # Tolerance broken: save the current group and start a new one
                groups.append(current_group)
                current_group = [pos]
            
        # Append the final group
        groups.append(current_group)

        # 2. Build the contig string
        contig_parts = []
        for group in groups:
            min_res = group[0]
            max_res = group[-1] 
        
            linker_length = random.randint(mergable_tolerance, 30)
        
            if len(group) > 1:
                contig_parts.append(f"{binder_chain_id}{min_res}-{max_res}/{mergable_tolerance}-{linker_length}")
            else:
                contig_parts.append(f"{binder_chain_id}{min_res}/{mergable_tolerance}-{linker_length}")

        return "/".join(contig_parts)
    
    def sample_motifs_from_binder_chain(self, binder_chain_id: str, approach_motif_scaffold: str = "rfdiffusion", mergable_tolerance: int = 4, cutoff: int = 6):
        """ Sample motif indices from a binder chain and save as a list of 1-indexed paratope indices """
        # Initialize default variables for contact check 
        desired_epitope_residues = []
        dict_contact_binder_target = determine_binding_interface(pdb_file_path= self.path_clean_structure, desired_epitope_residues= desired_epitope_residues, 
                            binder_chain_id= binder_chain_id, target_chain_id= self.target_chain_id, cutoff= cutoff)
        
        # Extract paratope indices
        paratope_indices = dict_contact_binder_target['paratope_indices'].split(',')
        paratope_indices = [int(x) for x in paratope_indices]
        
        # Filter paratope indices for those are truly desirable as contacts
        final_paratope_indices = list(set(paratope_indices) - set(self.binder_chain_dict[binder_chain_id]['undesirable_paratope_indices']))
        sampled_paratope_indices = sorted(random.sample(final_paratope_indices, k = min(len(final_paratope_indices), self.binder_chain_dict[binder_chain_id]['num_motifs'])))
        
        if approach_motif_scaffold == "rfdiffusion":
            contig_string = self.motif_indices_to_rfdcontig(sampled_paratope_indices, binder_chain_id, mergable_tolerance= mergable_tolerance)
            return contig_string
        else:
            return sampled_paratope_indices
    
    def sample_motifs_from_all_binder_chains(self, approach_motif_scaffold: str = "rfdiffusion", mergable_tolerance: int = 4, cutoff: int = 6):
        """ Sample motifs from all binder chains and return a whole contig string if using rfdiffusion or dictionary of {chain_id : paratope_indices} """
        chain_contigs = {}
        for binder_chain_id in self.binder_chain_dict.keys():
            if approach_motif_scaffold == "rfdiffusion":
                contig_string = self.sample_motifs_from_binder_chain(binder_chain_id, approach_motif_scaffold= approach_motif_scaffold, mergable_tolerance= mergable_tolerance,
                                                                    cutoff= cutoff)
                chain_contigs[binder_chain_id] = contig_string
            else:
                paratope_indices = self.sample_motifs_from_binder_chain(binder_chain_id, approach_motif_scaffold= approach_motif_scaffold, mergable_tolerance= mergable_tolerance,
                                                                        cutoff= cutoff)
                chain_contigs[binder_chain_id] = paratope_indices
        
        # Return the full joined/merged contig string or dictionary of chain IDs to paratope indices
        if approach_motif_scaffold == "rfdiffusion":
            chain_contigs = "/".join(chain_contigs.values())
            target_chain_contigs = f"{self.target_chain_id}:{chain_contigs}"
            return target_chain_contigs
        else:
            return chain_contigs


        