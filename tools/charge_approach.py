import os
import re
import numpy as np
import moleculegraph
from jinja2 import Template
from typing import List, Dict
from moleculegraph import molecule
from .numeric_utils import calc_mie, calc_coulomb, calc_bond

def get_local_dipol( molecule: molecule, force_field: Dict[str,Dict]) -> List[List]:
    """
    Function to get the local dipols of a molecule.

    Args:
        molecule (molecule): Moleculegraph object of the molecule under investigation.
        force_field (dict[dict]): Dictionary contains the force field types of the molecule under investigation.

    Returns:
        dipol_list (List[List]): List with sublists for each local dipol with the corresponding atom types in the dipol.
    """

    flag       = 0
    dipol_list = []
    dipol      = []

    # Search through all bonds and check if an atom is charged. If thats the case, start a local dipol list.
    # Check the neighboring atom (in the bond), and if its also charged, append it to the list and add up the local dipol charge. (if one atom is already in the local dipol list, skip the charge evaluation)
    # Once the charge is zero, a local dipol is identified and appended to the overall dipol list

    for bond_types in molecule.bond_list:

        for atom_type in bond_types:
            if atom_type in dipol: 
                continue

            elif force_field["atoms"][molecule.atom_names[atom_type]]["charge"] != 0 and flag == 0:
                dipol = [ atom_type ]
                chrg  = force_field["atoms"][molecule.atom_names[atom_type]]["charge"]
                flag  = 1
            
            elif force_field["atoms"][molecule.atom_names[atom_type]]["charge"] != 0 and flag == 1 :
                chrg += force_field["atoms"][molecule.atom_names[atom_type]]["charge"]
                dipol.append(atom_type)

                if np.round(chrg,3) == 0:
                    dipol_list.append(dipol)
                    flag = 0
            
    return dipol_list


def apply_charge_group_approach( mol_list: List[molecule], force_field: Dict[str,Dict], table_path: str):
        """
        Function that applies the charge group approach to provided moleculegraph components. This function alters the bond entries in the moleculegraph objects, 
        as well as the force field passed.

        Args:
            mol_list (List[molecule]): Moleculegraph objects of the molecules under investigation.
            force_field (dict[dict]): Dictionary contains the force field types of the molecule under investigation.
            table_path (str): Relative path from input file to the bond table
        """

        dipol_lists = [ get_local_dipol(mol,force_field) for mol in mol_list ]

        for j,(dipol_list, mol) in enumerate(zip( dipol_lists, mol_list )):
            
            # If dipol list is empty or only contains one local dipol, skip !
            if len(dipol_list) < 2: continue

            print("\nCharge group approach is applied to molecule: %s"%mol_list[j].molstring)

            dipol_names = [ "Local dipol n°%d: "%i + " ".join(mol.atom_names[np.array(dipol)]) for i,dipol in enumerate(dipol_list)] 
            
            print("\nThese are the local dipols identified:\n%s\n"%("\n".join(dipol_names)))
            
            # Use the dipol list, to create a special bond list, as well as the corresponding moleculegraph representation of it.
            # This means matching the correct force field types to each of the atoms in the special bonds list.
            special_bond_indexes = np.array( np.meshgrid( dipol_list[0], dipol_list[1] ) ).T.reshape(-1, 2)
            special_bond_names   = mol.atom_names[special_bond_indexes]
            special_bond_keys    = [ "[special]"+ moleculegraph.make_graph( moleculegraph.sort_force_fields(x) ) for x in special_bond_names ]

            # Delete unnecessary standard bonds. Checks which standard bonds also exist in the special bonds and therefore remove them
            idx             = np.array( [ i for i,entry in enumerate(mol.bond_list) if tuple(entry) not in set(map(tuple,special_bond_indexes)) ] )

            # Overwrite the bonding information of the moleculegraph object with the new bonds including the special bonds!
            mol.bond_list   = np.concatenate( [mol.bond_list[idx], special_bond_indexes], axis=0 )
            mol.bond_names  = np.concatenate( [mol.bond_names[idx], special_bond_names], axis=0 )
            mol.bond_keys   = np.concatenate( [mol.bond_keys[idx], special_bond_keys], axis=0 )

            mol.unique_bond_keys, mol.unique_bond_indexes, mol.unique_bond_inverse = moleculegraph.molecule_utils.unique_sort( mol.bond_keys, return_inverse=True )
            mol.unique_bond_names   = mol.bond_names[ mol.unique_bond_indexes ]
            mol.unique_bond_numbers = mol.bond_list[ mol.unique_bond_indexes ]    

            # Add the force field information for special bonds
            special_bond_keys_unique, uidx = np.unique(special_bond_keys, return_index=True)
            special_bond_names_unique      = special_bond_names[uidx]

            print("\nThese are the unique special bonds that are added for intramolecular interaction:\n%s\n"%("\n".join( " ".join(sb) for sb in special_bond_names_unique )))
            
            # Jinja2 template uses p[1] as first entry and p[0] as second --> thats why at first the key then the table name.
            for ( bkey, bname ) in zip(special_bond_keys_unique, special_bond_names_unique):
                force_field["bonds"][bkey] = { "list": list(bname), "p":  [ bkey, table_path ], "style": "table", "type": -1}

            # Write the new bonded force field and return it
            bonded_force_field = [j for sub in [molecule.map_molecule( molecule.unique_bond_keys, force_field["bonds"] ) for molecule in mol_list] for j in sub]

        return bonded_force_field



def write_tabled_bond( mol_list: List[molecule], force_field: Dict[str,Dict], 
                       table_path: str, table_template: str, evaluation_steps: int = 1000 ):
    """
    This function writes an input table for the LAMMPS bond style "table". The charge group approach is used, to include the 1.4-Mie and Coulombic interactions between local dipols.
    In the case two atoms of different local dipols are bonded, also evaluate the bonded potential.

    Args:
        mol_list (List[molecule]): Moleculegraph objects of the molecules under investigation.
        force_field (dict[dict]): Dictionary contains the force field types of the molecule under investigation.
        table_path (str): Path where the table should be written to.
        table_template (str): Path to the Jinja2 template to write the table.
        steps_nonbonded (int, optional): Number of evaluation points for all special bonds. Defaults to 1000.
    """

    tabled_dict_list = []

    # As the gradient of the Coulomb and Mie potential is quiete high at low distances, use more 50% of the intermediates in the first 30% of the cut off range. 
    # And less intermediates for the distances higher than 30% of the cut off.
    n1              = int( 0.5 * evaluation_steps )
    n2              = int( evaluation_steps - n1 )

    # Get overall cutoff used in the simulation
    cut_off         = max( p["cut"] for _,p in force_field["atoms"].items() )

    for mol in mol_list:
        # Loop through every molecule and check if there are special bonds
        for sbi,sbk in zip( mol.bond_list, mol.bond_keys ):

            # Check if it is a special bond and if the special bond is already evaluted. In that case, skip the reevaluation.
            if "special" in sbk and not sbk in [ td["list"] for td in tabled_dict_list ]:

                # Extract the force field keys of the atoms in this bond.
                ff_keys = re.findall(r'\[(.*?)\]', sbk)[1:]

                # Evalute the bonded interaction for every bonded special bond and the Coulomb interaction
                if mol.get_distance(*sbi) == 1:

                    # K is given in Kcal/mol/Angstrom^2, r0 in Angstrom
                    r0, K    = force_field["bonds"][moleculegraph.make_graph( ff_keys )]["p"]

                    # Get Coulomb parameter (charges and cut off)
                    charges  = [ force_field["atoms"][atom_key]["charge"] for atom_key in ff_keys ]
                    
                    # Distances evaluated ±30% around r0.
                    r_eval   = np.linspace( 0.7*r0, 1.3*r0, evaluation_steps )
                    
                    # Compute bonded interaction 
                    u_bond, f_bond  = calc_bond( r_eval, r0, K, "energy"), calc_bond( r_eval, r0, K, "force")

                    # Compute Coulomb interaction
                    u_coulomb, f_coulomb = calc_coulomb( r_eval, *charges, "energy"), calc_coulomb( r_eval, *charges, "force")

                    # Get total energy and force
                    u_total, f_total     = u_bond + u_coulomb, f_bond + f_coulomb


                # Evaluate the 1-4 Mie and Coulomb interaction for every special bond molecule.
                # Furthermore, include all intramolecular interactions with a higher bond distance than 3. (Since LAMMPS tabled bonds alter the
                # understanding of the molecule. Thus, if a tabled bond of a 1-4 pair is introduced, the 1-5 interaction do not longer have a bond distance of 4,
                # rather it now only has a distance of 2. Hence, LAMMPS will not compute non bonded interactions, and therefore they needed to be presented as tabled
                # bonds as well)
                elif mol.get_distance(*sbi) >= 3 :

                    # Get Coulomb and vdW parameter (charges, epsilon, sigma and repulsive exponent)
                    charges  = [ force_field["atoms"][atom_key]["charge"] for atom_key in ff_keys ]
                    epsilons = [ force_field["atoms"][atom_key]["epsilon"] for atom_key in ff_keys ]
                    sigmas   = [ force_field["atoms"][atom_key]["sigma"] for atom_key in ff_keys ]
                    ns       = [ force_field["atoms"][atom_key]["m"] for atom_key in ff_keys ]

                    # Use mixing Lorentz-Berthelot mixing rule for sigma and epsilon, as well as an arithmetic mean for the repulsive exponent
                    n        = np.mean( ns )
                    sigma    = np.mean( sigmas )
                    epsilon  = np.sqrt( np.prod( epsilons ) )

                    # Evaluated distances 
                    r_eval   = np.concatenate( [np.linspace( 0.01*cut_off, 0.3*cut_off, n1), np.linspace( 0.301*cut_off, cut_off, n2)] )
                    
                    # Compute Mie interaction up to the cut off radius
                    u_mie, f_mie         = calc_mie( r_eval, sigma, epsilon, n, 6, "energy"), calc_mie( r_eval, sigma, epsilon, n, 6, "force")

                    # Compute Coulomb interaction up to the cut off radius. 
                    u_coulomb, f_coulomb = calc_coulomb( r_eval, *charges, "energy"), calc_coulomb( r_eval, *charges, "force")

                    # Get total energy and force
                    u_total, f_total     = u_mie + u_coulomb, f_mie + f_coulomb
                

                # Evaluate only the Coulomb interaction (in the case for all 1-3 interactions)
                else:
                    
                    # Get Coulomb parameter (charges and cut off)
                    charges  = [ force_field["atoms"][atom_key]["charge"] for atom_key in ff_keys ]

                    # Evaluated distances 
                    r_eval   = np.concatenate( [np.linspace( 0.01*cut_off, 0.3*cut_off, n1), np.linspace( 0.301*cut_off, cut_off, n2)] )

                    # Compute Coulomb interaction up to the cut off radius. 
                    u_coulomb, f_coulomb = calc_coulomb( r_eval, *charges, "energy"), calc_coulomb( r_eval, *charges, "force")

                    # Get total energy and force
                    u_total, f_total     = u_coulomb, f_coulomb

                # Save the evaluated distance, energy and force into the table dict
                tabled_dict_list.append( { "list": sbk, 
                                            "N": len(r_eval), 
                                            "p": list( zip( range(1, evaluation_steps + 1), r_eval, u_total, f_total ) )
                                            } )

    # Write the table using Jinja2 template
    with open(table_template) as file_: 
        template = Template(file_.read())

    rendered = template.render( rd = tabled_dict_list )

    # Create folder (if necessary) where table should be writen to
    os.makedirs( os.path.dirname(table_path), exist_ok=True )

    with open(table_path, "w") as fh:
        fh.write(rendered) 
    
    return
