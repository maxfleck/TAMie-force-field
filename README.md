<h1 align="center">
  TAMie force field
</h1>
<p align="center">This repository enables users to perform molecular dynamics simulations utilizing LAMMPS with the Transferable Anisotropic Mie (TAMie) force field. The process begins with a SMILES and graph representation of each component, where PLAYMOL constructs initial systems at specific densities. The moleculegraph software and supplementary Python code are then used to generate a LAMMPS data and input file via jinja2 templates. </p>


## 🚀 Getting Started

Get started by running the following command 
```
git clone https://github.com/maxfleck/TAMie-force-field.git
```

## 🐍 Example program

The following example will demonstrate how to setup MD simulations (this is also demonstrated in the example.ipynb). Essentially, the workflow can be summarized as follows:
1. Get intial molecule coordinates (e.g.: using SMILES and PubChem) and provide a graph representation.
2. Use PLAYMOL to construct a system of any mixture.
3. Use the moleculegraph software and pyLAMMPS tools to generate LAMMPS data and input files via jinja2 templates.

## 1. Define general settings ##

1. Define path to force field toml (in a moleculegraph understandable format),
2. Define name, SMILES, and graph strings of the molecules under investigation.
3. If the charge group approach is utilized, define per component the two atoms in the torsion that was reparametrized (e.g.: 1,2-ethanediol: ["OH_alcohol","OH_alcohol"] )
4. Define the relative path from the final simulation folder to the tabled bond interactions 
5. Call the python LAMMPS input class

```python
# 1: Path tho force field toml
force_field_path     = "force-fields/forcefield_lammps.toml"

# 2: Names, SMILES, and graphs of molecules (further examples on constructing molecule graphs available at https://github.com/maxfleck/moleculegraph)
molecule_name1       = "ethanediol"
molecule_graph1      = "[cH_alcohol][OH_alcohol][CH2_alcohol][CH2_alcohol][OH_alcohol][cH_alcohol]"
molecule_smiles1     = "OCCO"

molecule_name_list   = [ molecule_name1 ]
molecule_graph_list  = [ molecule_graph1 ]
molecule_smiles_list = [ molecule_smiles1 ]

# Path to xyz template
template_xyz         = "tools/templates/template_write_xyz.xyz"

# Define output path for final xyz files
xyz_destinations     = [ "example/%s/%s.xyz"%(name,name) for name in molecule_name_list ]

# Get the single molecule coordinates for each component
get_molecule_coordinates( molecule_name_list = molecule_name_list, molecule_graph_list = molecule_graph_list, molecule_smiles_list = molecule_smiles_list,
                          xyz_destinations = xyz_destinations, template_xyz = template_xyz, verbose = False )


# 3: Define the two atoms in each reparametrized torsion per component
torsion_pairs       = [ [ "OH_alcohol", "OH_alcohol" ] ]

# 4: Relative path to the tabled bond interactions
table_path          = "../bonded_interactions.table"

# 5: Call the LAMMPS input class
LAMMPS_class        = LAMMPS_input( mol_str = molecule_graph_list, ff_path = force_field_path, charge_group_approach = True, table_path = table_path )
```

## 2. Write system size independent files ##

1. The PLAYMOL force field file, which is used to build the initial configurations
2. The tabled bond interactions for the system

```python
# 1: PLAYMOL force field file
playmol_force_field_template    = "tools/templates/template_playmol_forcefield.playmol"
playmol_force_field_destination = "example/ethandiol/playmol_ff.playmol"

LAMMPS_class.prepare_playmol_input( playmol_template = playmol_force_field_template, playmol_ff_path = playmol_force_field_destination )

# 2: Tabled bond interactions
LAMMPS_table_template          = "tools/templates/template_lammps_bond_table.table"
LAMMPS_table_destination       = "example/ethandiol/bonded_interactions.table"

LAMMPS_class.write_tabled_bond( torsion_pairs = torsion_pairs, table_path = LAMMPS_table_destination, table_template = LAMMPS_table_template )
```

## 3. Write system size dependent files ##

1. Define general (thermodynamic) settings for each system.
2. Write PLAYMOL input file and execute it to build the system (if wanted).
3. Write LAMMPS data file using the from PLAYMOL generated xyz file
4. Write LAMMPS input file for each system

```python
# 1: Define general settings for each system
# Temperatures [K], pressures [bar] (if wanted, otherwise use 0.0) and initial denisties [kg/m^3] for each system. Also define the number of molecules per component.

temperatures     = [ 273.15, 298.15, 323.15, 348.15, 373.15 ]
pressures        = [ 1.01325, 1.01325, 1.01325, 1.01325, 1.01325 ]
densities        = [ 1092.677, 1077.871, 1063.222, 1048.52, 1033.567 ]
molecule_numbers = [ 500 ]

# Simulation path
simulation_path  = "example/ethandiol/sim_%d"

# Define if PLAYMOL should be executed and further settings
build_playmol             = True
molecule_xyz_files        = [ "../1.2-Ethanediol_playmol.xyz" ]
relative_playmol_ff_path  = "../playmol_ff.playmol"
playmol_input_template    = "tools/templates/template_playmol_input.mol"

# Define LAMMPS template paths
LAMMPS_data_template      = "tools/templates/template_lammps_data.data"
LAMMPS_input_template     = "tools/templates/template_lammps_input.in"

for i, (temp, press, dens) in enumerate( zip( temperatures, pressures, densities ) ):

    # Create the simulation folder (if not already done)
    os.makedirs( simulation_path%i, exist_ok = True )

    # Prepare LAMMPS with molecules numbers and density of the system
    LAMMPS_class.prepare_lammps_data (nmol_list = molecule_numbers, densitiy = dens )

    # 2: Write PLAYMOL input and execute (if wanted.)
    if build_playmol:
        playmol_input_destination = simulation_path%i + f"/ethandiol_{i}.mol"
        
        LAMMPS_class.write_playmol_input( playmol_template = playmol_input_template, playmol_path = playmol_input_destination, 
                                          playmol_ff_path = relative_playmol_ff_path, xyz_paths = molecule_xyz_files )

    # 3: Write LAMMPS data file from generated xyz file
    system_xyz              = simulation_path%i + f"/ethandiol_{i}.xyz"
    LAMMPS_data_destination = simulation_path%i + "/lammps.data"

    LAMMPS_class.write_lammps_data( xyz_path = system_xyz, data_template = LAMMPS_data_template, data_path = LAMMPS_data_destination )

    # 4: Write LAMMPS input file
    LAMMPS_input_destination  = simulation_path%i + "/lammps.input"
    relative_LAMMPS_data_path = LAMMPS_data_destination[ LAMMPS_data_destination.rfind(simulation_path%i) + len(simulation_path%i) + 1 : ]
   
    LAMMPS_class.prepare_lammps_input( )
    LAMMPS_class.write_lammps_input( input_path = LAMMPS_input_destination, template_path = LAMMPS_input_template, data_file = relative_LAMMPS_data_path,
                                     temperature = temp, pressure = press, equilibration_time = 2e6, production_time = 4e6 )
```

## 🚑 Help

Help will arrive soon ...

## 👫 Authors

Maximilian Fleck - University of Stuttgart, Samir Darouich - University of Stuttgart

## 📄 License

This project is licensed under the MIT License - see the LICENSE.md file for details