<h1 align="center">
  TAMie force field
</h1>
<p align="center">This repository enables users to perform molecular dynamics simulations utilizing LAMMPS with the Transferable Anisotropic Mie (TAMie) force field. The process begins with a SMILES and graph representation of each component, where PLAYMOL constructs initial systems at specific densities. The moleculegraph and pyLAMMPs software are then used to generate a LAMMPS data and input file via jinja2 templates. </p>


## 🚀 Getting Started

Get started by running the following command to install:

1. moleculegraph
```
git clone https://github.com/maxfleck/moleculegraph
cd moleculegraph
pip install -I .
```
2. pyLAMMPS
```
git clone https://github.com/samirdarouich/pyLAMMPS.git
cd pyLAMMPS
pip install -I .
```
3. PLAYMOL
```
git clone https://github.com/atoms-ufrj/playmol
cd playmol
make
```

## 🐍 Example program

## Simulation setup

pyLAMMPS enanbles the structured and FAIR setup of simulation systems. Together with moleculegraph and PLAYMOL it enables the building of an initial configuration and to write the necessary LAMMPS data and input files

1) Read in the YAML files to define the system and simulation/sampling settings. Also provide the simulated ensembles and further settings

```python
# Read in YAML files
lammps_setup = LAMMPS_setup( system_setup = "input/setup.yaml", 
                             simulation_default = "input/defaults.yaml",
                             simulation_ensemble = "input/ensemble.yaml",
                             simulation_sampling = "input/sampling.yaml",
                             submission_command = "qsub"
                            )

# Define simulation folder
sim_folder = f'{lammps_setup.system_setup["folder"]}/{lammps_setup.system_setup["name"]}'

# Define the ensembles that should be simulated (definition what each ensemble means is provided in yaml file)
ensembles = [ "em", "npt" ] 

# Define the simulation time per ensemble in nano seconds (for em the number of iterations is provided in the ensemble yaml)
simulation_times = [ 0, 60.0 ]

# Define initial systems. This can be .data or .restart files. For each temperature & pressure state define one system
initial_systems = [ "/home/st/st_st/st_ac137577/workspace/software/TAMie-force-field/example/1.2-ethanediol/temp_298_pres_1/build/system.data" ]

# Define number of copies
copies = 2

# Define if the initial system (if not provided) should be build on a cluster or local machine
on_cluster = True

# Define any further input kwargs for the input template
input_kwargs = {}

# Define the starting number for the first ensemble ( 0{off_set}_ensemble )
off_set = 0
```

2) Utilize moleculegraph and the LAMMPS_molecules class provided by pyLAMMPS to map the force field parameters with the molecules. <br>
   Apply the charge group approach and write LAMMPS force field input


```python
# Call the LAMMPS molecule class
lammps_molecules = LAMMPS_molecules( mol_str = [ mol["graph"] for mol in lammps_setup.system_setup["molecules"] ],
                                     force_field_path = lammps_setup.system_setup["paths"]["force_field_path"] 
                                    ) 

# Use the charge group approach (and overwrite the bonded information of the moleculegraph objects)
lammps_molecules.bonds = apply_charge_group_approach( mol_list = lammps_molecules.mol_list, 
                                                      force_field = lammps_molecules.ff, 
                                                      table_path = f"{sim_folder}/bonded_interactions.table" 
                                                    )

# Prepare the LAMMPS force field (this gathers all the necessary input from the force field toml file)
lammps_molecules.prepare_lammps_force_field()

# Get shake dictionary
shake_dict = lammps_molecules.get_shake_indices( lammps_setup.simulation_default["shake_dict"] )

# Write LAMMPS force field file
# Define the number of spline interpolations LAMMPS stores of the tabled potential (https://docs.lammps.org/bond_table.html)
lammps_ff_file = write_lammps_ff( ff_template = lammps_setup.system_setup["paths"]["template"]["lammps_ff_file"], 
                                lammps_ff_path = f"{sim_folder}/force_field.params", 
                                potential_kwargs = { **lammps_setup.simulation_default["non_bonded"]["vdw_style"], 
                                                     **lammps_setup.simulation_default["non_bonded"]["coulomb_style"] },
                                atom_numbers_ges = lammps_molecules.atom_numbers_ges, 
                                nonbonded = lammps_molecules.nonbonded, 
                                bond_numbers_ges = lammps_molecules.bond_numbers_ges, 
                                bonds = lammps_molecules.bonds,
                                angle_numbers_ges = lammps_molecules.angle_numbers_ges, 
                                angles = lammps_molecules.angles,
                                torsion_numbers_ges = lammps_molecules.torsion_numbers_ges, 
                                torsions = lammps_molecules.torsions,
                                only_self_interactions = lammps_setup.simulation_default["non_bonded"]["lammps_mixing"], 
                                mixing_rule = lammps_setup.simulation_default["non_bonded"]["mixing"],
                                ff_kwargs = lammps_setup.simulation_default["non_bonded"],
                                n_eval = 1000 
                                )

# Write tabled bond interactions (at the same destination as the LAMMPS ff file)
lammps_table_file = write_tabled_bond( mol_list = lammps_molecules.mol_list, 
                                       force_field = lammps_molecules.ff, 
                                       table_path = f"{sim_folder}/bonded_interactions.table", 
                                       table_template = lammps_setup.system_setup["paths"]["template"]["lammps_table_file"],
                                       evaluation_steps = 1000
                                      )
```

3) Write input (and data) file for every temperature & pressure state

```python
lammps_setup.job_files = []

for i, (temperature, pressure, density) in enumerate( zip( lammps_setup.system_setup["temperature"], 
                                                           lammps_setup.system_setup["pressure"], 
                                                           lammps_setup.system_setup["density"] ) ):
    
    job_files = []
    
    # Define folder for specific temp and pressure state
    state_folder = f"{sim_folder}/temp_{temperature:.0f}_pres_{pressure:.0f}"

    # Build system with PLAYMOL and write LAMMPS data if no initial system is provided
    if not initial_systems:
        
        lammps_data_file = generate_initial_configuration( lammps_molecules = lammps_molecules,
                                                            destination_folder = state_folder,
                                                            molecules_dict_list = lammps_setup.system_setup["molecules"],
                                                            density = density,
                                                            template_xyz = lammps_setup.system_setup["paths"]["template"]["xyz_file"],
                                                            playmol_ff_template = lammps_setup.system_setup["paths"]["template"]["playmol_ff_file"],
                                                            playmol_input_template = lammps_setup.system_setup["paths"]["template"]["playmol_input_file"],
                                                            playmol_bash_file = lammps_setup.system_setup["paths"]["template"]["playmol_bash_file"],
                                                            lammps_data_template = lammps_setup.system_setup["paths"]["template"]["lammps_data_file"],
                                                            submission_command = lammps_setup.submission_command, 
                                                            on_cluster = on_cluster
                                                        )
    
        flag_restart = False
    else:
        lammps_data_file = initial_systems[i]
        print(f"\nIntial system provided for at: {lammps_data_file}\n")
        flag_restart = ".restart" in lammps_data_file
        if flag_restart: 
            print("Restart file is provided. Continue simulation from there!\n")

    # Define folder for each copy
    for copy in range( copies + 1 ):
        copy_folder = f"{state_folder}/copy_{copy}"

        # Produce input files (for each ensemble an own folder 0x_ensemble)
        input_files = generate_input_files( destination_folder = copy_folder, 
                                            input_template = lammps_setup.system_setup["paths"]["template"]["lammps_input_file"],
                                            ensembles = ensembles, 
                                            temperature = temperature, 
                                            pressure = pressure,
                                            data_file = lammps_data_file, 
                                            ff_file = lammps_ff_file,
                                            simulation_times = simulation_times,
                                            dt = lammps_setup.simulation_default["system"]["dt"], 
                                            kwargs = { **lammps_setup.simulation_default,
                                                        **lammps_setup.simulation_sampling, 
                                                        **input_kwargs,
                                                        "style": style_dict,
                                                        "shake_dict": shake_dict, 
                                                        "restart_flag": flag_restart }, 
                                            ensemble_definition = lammps_setup.simulation_ensemble,
                                            off_set = off_set
                                            )
        
        # Create job file
        job_files.append( generate_job_file( destination_folder = copy_folder, 
                                             job_template = lammps_setup.system_setup["paths"]["template"]["job_file"], 
                                             input_files = input_files, 
                                             ensembles = ensembles,
                                             job_name = f'{lammps_setup.system_setup["name"]}_{temperature:.0f}_{pressure:.0f}',
                                             job_out = f"job_{temperature:.0f}_{pressure:.0f}.sh", 
                                             off_set = off_set 
                                            ) 
                        )
        
    lammps_setup.job_files.append( job_files )

```
4) Submit the simulations
```
lammps_setup.submit_simulation()
```

## 🚑 Help

Help will arrive soon ...

## 👫 Authors

Maximilian Fleck - University of Stuttgart, Samir Darouich - University of Stuttgart

## 📄 License

This project is licensed under the MIT License - see the LICENSE.md file for details
