import os
import yaml

class CreateRFDYaml:
    def __init__(self, path_yaml_template: str = ""):
        
        self.dict_template = self.create_template_dict(path_yaml_template = path_yaml_template)

    def create_rfd_dict(self):
        """ Create a dictionary corresponding to RFDiffusion Hyperparameters for final YAML file"""
        dict_rfd = {'rfdiffusion' : {'contigs' : "A:40-50/B5-13/20-30/C5-13/3-8",
                                     'pdb' : '/Volumes/sandbox/karthik/motif_scaffolding/alpha2_beta1_double/input/collagen_triple_helix_alpha2_model_0.pdb',
                                     'iterations' : 50,
                                     'hotspot' : 'A12,A13,A14,A15,A74,A75,A76,A77,A78,A79,A80,A115,A116,A117',
                                     'num_designs' : 12,
                                     'visual' : 'image',
                                     'symmetry': 'none',
                                     'order' : 1,
                                     'chains' : '',
                                     'add_potential' : True,
                                     'partial_T' : 'auto',
                                     'use_beta_model' : False,
                                     }
                    }
        return dict_rfd
    
    def create_mpnn_dict(self):
        """ Create a dictionary corresponding to MPNN Hyperparameters for final YAML file"""
        dict_mpnn = {'mpnn' :{'num_seqs' : 8,
                              'mpnn_sampling_temp' : 0.1,
                              'rm_aa' : 'C',
                              'use_solubleMPNN' : True,
                              }
                     }
        return dict_mpnn
    
    def create_af2_dict(self):
        """ Create a dictionary corresponding to AlphaFold2 Hyperparameters for final YAML file"""
        dict_af2 = {'alphafold2' : {'initial_guess' : True, 
                                    'num_recycles' : 3,
                                    'use_multimer' : False
                                    }
                    }
        return dict_af2
    
    def create_template_dict(self, path_yaml_template: str = ""):
        """ Create overarching dictionary for running Sergey's RFDiffusion -> MPNN -> AF2 pipeline
            Args:
                path_yaml_template (str): Path to YAML template file. Allows for modification of pre-made yaml files if provided as input or just create a new one if empty
            Returns:
                dict_template (dict): Dictionary containing all tunable hyperparameters for running Sergey's RFDiffusion -> MPNN -> AF2 pipeline

        """
        if path_yaml_template != "":
            with open(path_yaml_template, 'r') as file:
                dict_template = yaml.safe_load(file)
        else:
            dict_rfd = self.create_rfd_dict()
            dict_mpnn = self.create_mpnn_dict()
            dict_af2 = self.create_af2_dict()
            dict_template = dict_rfd | dict_mpnn | dict_af2
        
        return dict_template
    
    def update_template_dict(self, pipeline_step: str,  hyperparameter: str, value):
        """ Update template dictionary for running Sergey's RFDiffusion -> MPNN -> AF2 pipeline
            Args:
                pipeline_step (str): Name of pipeline step to update
                hyperparameter (str): Name of hyperparameter to update
                value: Value to update hyperparameter with (Value could be int, float, str or bool) varies based on pipeline step and hyperparameter chosen to be modified
        """
        # Check if pipeline step chosen to be updated is a valid option
        if pipeline_step not in self.dict_template.keys():
            raise ValueError(f"Pipeline step: {pipeline_step} not found in template dictionary")
            print(f"Modifiable pipeline steps are: {self.dict_template.keys()}")
        # If prior check passes, verify if hyperparameter chosen to be updated is a valid option
        if hyperparameter not in self.dict_template[pipeline_step].keys():
            raise ValueError(f"Hyperparameter {hyperparameter} not found in pipeline step {pipeline_step} tunable options")
            print(f"Modifiable hyperparameters for pipeline step {pipeline_step} are: {self.dict_template[pipeline_step].keys()}")
        # If prior check passes, check if value is the same data type as the current value so that the update is valid
        if type(value) != type(self.dict_template[pipeline_step][hyperparameter]):
            raise ValueError(f"Value {value} is not the same data type as the current value {self.dict_template[pipeline_step][hyperparameter]}")
            print(f"Value {value} should be of type {type(self.dict_template[pipeline_step][hyperparameter])}")
        # Update template dictionary
        self.dict_template[pipeline_step][hyperparameter] = value

    def write_yaml(self, path_save_folder: str, yaml_name: str):
        """ Create YAML file for running Sergey's RFDiffusion -> MPNN -> AF2 pipeline & save to specified save_path
        Args:
            path_save_folder (str): Path to folder where YAML file should be saved
            yaml_name (str): Name of YAML file

        """
        # Check if save folder exists, if not create it
        if not os.path.exists(path_save_folder):
            print(f"Folder: {path_save_folder} to save YAML file does not exist, so creating folder instead")
            os.makedirs(path_save_folder, exist_ok=True)
        
        # Create YAML file save-path
        path_save_yaml = os.path.join(path_save_folder, yaml_name)
        print(f"Saving YAML file to: {path_save_yaml}")
        
        # Write YAML file
        with open(path_save_yaml, 'w') as file:
            documents = yaml.dump(self.dict_template, file)