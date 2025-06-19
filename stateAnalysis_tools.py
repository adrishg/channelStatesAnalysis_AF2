import math
import re
import os # Added for os.path.exists check

# Function to extract B-factors (which are pLDDT scores in AlphaFold models) from a PDB file
# This function already correctly uses the chainID parameter.
def extractBFactor(pdbFile, chainID):
    patternATOM = re.compile(
        r'^ATOM\s+'
        r'\d+\s+'
        r'CA\s+'
        r'([A-Z]{3})\s+'
        + re.escape(chainID) +  # Chain ID is explicitly used here
        r'\s?(\d+)\s+'
        r'[-+]?\d*\.\d+\s+'
        r'[-+]?\d*\.\d+\s+'
        r'[-+]?\d*\.\d+\s+'
        r'\d*\.\d+\s+'
        r'(\d+\.\d+)'
    )
    bFactors = {}

    with open(pdbFile, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and " CA " in line:
                atomInfo = re.search(patternATOM, line)
                if atomInfo:
                    resName = atomInfo.group(1)
                    resNumber = int(atomInfo.group(2))
                    bFactor = float(atomInfo.group(3))
                    bFactors[resName + ' ' + str(resNumber)] = bFactor

    return bFactors

# Function to calculate the average pLDDT from the extracted B-factors
def calculate_average_plddt(bFactors):
    if not bFactors:
        return None
    return sum(bFactors.values()) / len(bFactors)


# Define the side chain atoms for each amino acid, including CA
sidechain_atoms = {
    'ALA': ['CA', 'CB'], 'ARG': ['CA', 'NH2', 'NH1', 'CZ', 'NE', 'CD', 'CG', 'CB'],
    'ASN': ['CA', 'ND2', 'OD1', 'CG', 'CB'], 'ASP': ['CA', 'OD2', 'OD1', 'CG', 'CB'],
    'CYS': ['CA', 'SG', 'CB'], 'GLN': ['CA', 'NE2', 'OE1', 'CD', 'CG', 'CB'],
    'GLU': ['CA', 'OE2', 'OE1', 'CD', 'CG', 'CB'], 'GLY': ['CA'],
    'HIS': ['CA', 'NE2', 'CE1', 'CD2', 'ND1', 'CG', 'CB'], 'ILE': ['CA', 'CD1', 'CG2', 'CG1', 'CB'],
    'LEU': ['CA', 'CD2', 'CD1', 'CG', 'CB'], 'LYS': ['CA', 'NZ', 'CE', 'CD', 'CG', 'CB'],
    'MET': ['CA', 'CE', 'SD', 'CG', 'CB'], 'PHE': ['CA', 'CZ', 'CE2', 'CE1', 'CD2', 'CD1', 'CG', 'CB'],
    'PRO': ['CA', 'CD', 'CG', 'CB'], 'SER': ['CA', 'OG', 'CB'],
    'THR': ['CA', 'OG1', 'CG2', 'CB'], 'TRP': ['CA', 'CH2', 'CZ3', 'CZ2', 'CE3', 'CE2', 'NE1', 'CD2', 'CD1', 'CG', 'CB'],
    'TYR': ['CA', 'OH', 'CZ', 'CE2', 'CE1', 'CD2', 'CD1', 'CG', 'CB'], 'VAL': ['CA', 'CG2', 'CG1', 'CB']
}

# Define the extreme side chain atoms for each amino acid
extreme_sidechain_atoms = {
    'ALA': ['CB'], 'ARG': ['NH2', 'NH1'], 'ASN': ['ND2', 'OD1'], 'ASP': ['OD2', 'OD1'],
    'CYS': ['SG'], 'GLN': ['NE2', 'OE1'], 'GLU': ['OE2', 'OE1'], 'GLY': ['CA'],
    'HIS': ['NE2', 'CE1'], 'ILE': ['CD1', 'CG2'], 'LEU': ['CD2', 'CD1'], 'LYS': ['NZ'],
    'MET': ['CE'], 'PHE': ['CZ'], 'PRO': ['CD'], 'SER': ['OG'], 'THR': ['OG1'],
    'TRP': ['CH2'], 'TYR': ['OH'], 'VAL': ['CG2', 'CG1']
}

def parse_pdb(file_path):
    """
    Parses a PDB file and returns its lines.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def get_atom_coordinates(lines, residue_type, residue_number, atom_name, chain=None):
    """
    Extracts the coordinates of a specific atom from a PDB file's lines,
    optionally filtering by chain ID.
    This function already correctly handles the 'chain' parameter.
    """
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_residue_type = line[17:20].strip()
            atom_residue_number = int(line[22:26].strip())
            atom_chain = line[21].strip()  # Chain identifier
            atom_name_field = line[12:16].strip()

            # Ensure we're only working with the specified chain, if provided
            if chain and atom_chain != chain:
                continue
            
            # Check if the residue type, number, and atom name match
            if (atom_residue_type == residue_type and atom_residue_number == residue_number 
                and atom_name_field == atom_name):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    return (x, y, z)
                except ValueError:
                    print(f"Error parsing coordinates for {residue_type}{residue_number} {atom_name}")
                    return None
    print(f"Could not find {atom_name} atom for {residue_type} {residue_number}" + (f" in chain {chain}" if chain else ""))
    return None


def calculate_distance(coord1, coord2):
    """
    Calculates the Euclidean distance between two 3D coordinates.
    """
    return math.sqrt(sum([(a - b) ** 2 for a, b in zip(coord1, coord2)]))

def measure_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number, chain=None):
    """
    Measures all atom-to-atom distances between two specified residues,
    optionally within a given chain.
    """
    lines = parse_pdb(pdb_file)

    distances = []

    for atom1 in sidechain_atoms[residue1_type]:
        # Pass the chain argument to get_atom_coordinates
        coord1 = get_atom_coordinates(lines, residue1_type, residue1_number, atom1, chain)
        if coord1:
            for atom2 in sidechain_atoms[residue2_type]:
                # Pass the chain argument to get_atom_coordinates
                coord2 = get_atom_coordinates(lines, residue2_type, residue2_number, atom2, chain)
                if coord2:
                    distance = calculate_distance(coord1, coord2)
                    distances.append((distance, residue1_type, residue1_number, atom1, residue2_type, residue2_number, atom2))
                    print(f"Distance between {residue1_type}{residue1_number} ({atom1})" + # Corrected variable name here: residue1_number
                          (f" in chain {chain}" if chain else "") +
                          f" and {residue2_type}{residue2_number} ({atom2})" +
                          (f" in chain {chain}" if chain else "") +
                          f": {distance:.2f} Å")

    return distances

def measure_ca_distances(pdb_file, residue_pairs, chain=None):
    """
    Measures the distances between the CA atoms of the provided residue pairs, optionally by chain.
    This function already correctly handles the 'chain' parameter.
    
    Parameters:
    pdb_file (str): The path to the PDB file.
    residue_pairs (list of tuples): A list of tuples where each tuple contains two residue types and numbers.
                                    Example: [('ARG', 10, 'GLU', 20), ('ARG', 126, 'PHE', 59)]
    chain (str): The chain identifier (optional). If not provided, it will assume all residues are in the same chain.
    
    Returns:
    dict: A dictionary of residue pairs and their CA distances in the format "RESNUM_CA".
    """

    lines = parse_pdb(pdb_file)
    ca_coords = []

    # Dictionary to store residue pairs and their distances
    residue_pairs_distances = {}

    # Iterate over the residue pairs and calculate distances
    for residue1_type, residue1_number, residue2_type, residue2_number in residue_pairs:
        # Pass the chain argument to get_atom_coordinates
        coord1 = get_atom_coordinates(lines, residue1_type, residue1_number, 'CA', chain)
        coord2 = get_atom_coordinates(lines, residue2_type, residue2_number, 'CA', chain)
        
        if coord1 and coord2:
            dist = calculate_distance(coord1, coord2)
            pair_key = f"{residue1_type}{residue1_number}_CA-{residue2_type}{residue2_number}_CA"
            residue_pairs_distances[pair_key] = dist
            print(f"{pair_key}: Distance = {dist:.2f} Å")
        else:
            print(f"Could not find CA atom for one or both residues in pair: {residue1_type}{residue1_number}, {residue2_type}{residue2_number}" +
                  (f" in chain {chain}" if chain else ""))

    return residue_pairs_distances


# Very similar to measure_(all)distances but in theory should be faster because only iterate on the tips of aas
def measure_extreme_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number, chain=None):
    """
    Measures distances between extreme side chain atoms of two residues,
    optionally within a given chain.
    """
    lines = parse_pdb(pdb_file)

    distances = []

    for atom1 in extreme_sidechain_atoms[residue1_type]:
        # Pass the chain argument to get_atom_coordinates
        coord1 = get_atom_coordinates(lines, residue1_type, residue1_number, atom1, chain)
        if coord1:
            for atom2 in extreme_sidechain_atoms[residue2_type]:
                # Pass the chain argument to get_atom_coordinates
                coord2 = get_atom_coordinates(lines, residue2_type, residue2_number, atom2, chain)
                if coord2:
                    distance = calculate_distance(coord1, coord2)
                    distances.append((distance, residue1_type, residue1_number, atom1, residue2_type, residue2_number, atom2))

    return distances

def measure_shortest_distance(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number, chain=None):
    """
    Measures the shortest distance between extreme atoms of two residues and returns the formatted string
    with the residue and atom information, optionally within a given chain.
    
    Parameters:
    pdb_file (str): The path to the PDB file.
    residue1_type (str): The residue type of the first residue (e.g., 'GLU').
    residue1_number (int): The residue number of the first residue.
    residue2_type (str): The residue type of the second residue (e.g., 'ARG').
    residue2_number (int): The residue number of the second residue.
    chain (str): The chain identifier (optional). If not provided, it will consider all chains.
    
    Returns:
    str: A formatted string with the shortest distance in the format "RESNUM_ATOM".
    """
    
    # Pass the chain argument to measure_extreme_distances
    distances = measure_extreme_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number, chain)
    if distances:
        # Find the shortest distance
        shortest_distance = min(distances, key=lambda x: x[0])
        distance, res1_type, res1_num, atom1, res2_type, res2_num, atom2 = shortest_distance
        
        # Format the output to include both residue and atom info, e.g., GLU20_OE1 or ARG126_NH1
        res1_atom = f"{res1_type}{res1_num}_{atom1}"
        res2_atom = f"{res2_type}{res2_num}_{atom2}"
        
        # Round the distance to 1 decimal point
        distance = round(distance, 1)
        
        # Print and return the formatted result
        print(f"\nShortest distance between {res1_atom} and {res2_atom}: {distance:.1f} Å" +
              (f" in chain {chain}" if chain else ""))
        
        return f"{res1_atom}-{res2_atom}: {distance:.1f} Å"
    else:
        print("Could not find specified atoms in the PDB file." + (f" for chain {chain}" if chain else ""))
        return None  # Return None if no distance is found

def calculate_area_from_residues_with_chain(pdb_file, residues, chain=None):
    """
    Calculates the area formed by connecting the C-alpha atoms of the provided residues, optionally by chain.
    This function already correctly handles the 'chain' parameter.
    
    Parameters:
    pdb_file (str): The path to the PDB file.
    residues (list of tuples): A list of tuples where each tuple contains the residue type and residue number.
                                Example: [('ARG', 10), ('GLU', 20), ('SER', 30)]
    chain (str): The chain identifier (optional). If not provided, it will assume all residues are in the same chain.
    
    Returns:
    tuple: A dictionary of residue pairs and distances, and the calculated area of the polygon.
    """

    lines = parse_pdb(pdb_file)
    ca_coords = []

    # Extract the CA coordinates for each residue
    for residue_type, residue_number in residues:
        # Pass the chain argument to get_atom_coordinates
        coord = get_atom_coordinates(lines, residue_type, residue_number, 'CA', chain)
        if coord:
            ca_coords.append((residue_type, residue_number, coord))
        else:
            print(f"Could not find CA atom for {residue_type}{residue_number} in chain {chain if chain else 'any chain'}")
            return None

    if len(ca_coords) < 3:
        print("At least 3 residues are required to calculate an area.")
        return None

    # Dictionary to store residue pairs and their distances
    residue_pairs_distances = {}

    # print the residue pairs and their coordinates #DEBUGGG
    print("\nResidue pairs and coordinates:")
    for i in range(len(ca_coords)):
        res1 = ca_coords[i]
        res2 = ca_coords[(i + 1) % len(ca_coords)]  # This will loop back to the first residue to close the polygon
        dist = calculate_distance(res1[2], res2[2])

        # Create pair key in the format resThreeLetterCodeResNum-resThreeLetterCodeResNum
        pair_key = f"{res1[0]}{res1[1]}-{res2[0]}{res2[1]}"
        residue_pairs_distances[pair_key] = dist

        print(f"{res1[0]}{res1[1]} (CA: {res1[2]}) to {res2[0]}{res2[1]} (CA: {res2[2]}): Distance = {dist:.2f} Å")

    # Project the 3D coordinates onto a 2D plane (ignoring z-coordinate)
    projected_coords = [(coord[0], coord[1]) for _, _, coord in ca_coords]

    # Calculate the area using the Shoelace's formula
    n = len(projected_coords)
    area = 0

    for i in range(n):
        x1, y1 = projected_coords[i]
        x2, y2 = projected_coords[(i + 1) % n]  # Connects the last vertex back to the first
        area += x1 * y2 - x2 * y1

    area = abs(area) / 2.0

    print(f"\nArea formed by the residues: {area:.2f} square Å" + (f" in chain {chain}" if chain else ""))
    
    return residue_pairs_distances, area

# CORRECTION: previous version of area vas considering tips
def calculate_area_from_ca_distances(pdb_file, residues, chain=None):
    """
    Calculates the area formed by connecting the C-alpha atoms of the provided residues, optionally by chain.
    This function already correctly handles the 'chain' parameter.
    
    Parameters:
    pdb_file (str): The path to the PDB file.
    residues (list of tuples): A list of tuples where each tuple contains the residue type and residue number.
                                Example: [('ARG', 10), ('GLU', 20), ('SER', 30)]
    chain (str): The chain identifier (optional). If not provided, it will assume all residues are in the same chain.
    
    Returns:
    tuple: A dictionary of residue pairs and distances, and the calculated area of the polygon in the format RESNUM_CA.
    """

    lines = parse_pdb(pdb_file)
    ca_coords = []

    # Extract the CA coordinates for each residue
    for residue_type, residue_number in residues:
        # Pass the chain argument to get_atom_coordinates
        coord = get_atom_coordinates(lines, residue_type, residue_number, 'CA', chain)
        if coord:
            ca_coords.append((residue_type, residue_number, coord))
        else:
            print(f"Could not find CA atom for {residue_type}{residue_number} in chain {chain if chain else 'any chain'}")
            return None

    if len(ca_coords) < 3:
        print("At least 3 residues are required to calculate an area.")
        return None

    # Dictionary to store residue pairs and their distances
    residue_pairs_distances = {}

    # Print the residue pairs and their CA distances
    print("\nResidue pairs and CA distances:")
    for i in range(len(ca_coords)):
        res1 = ca_coords[i]
        res2 = ca_coords[(i + 1) % len(ca_coords)]  # This will loop back to the first residue to close the polygon
        dist = calculate_distance(res1[2], res2[2])

        # Create pair key in the format RESNUM_CA (e.g., GLU20_CA)
        pair_key = f"{res1[0]}{res1[1]}_CA-{res2[0]}{res2[1]}_CA"
        residue_pairs_distances[pair_key] = dist

        print(f"{pair_key}: Distance = {dist:.2f} Å")

    # Project the 3D coordinates onto a 2D plane (ignoring z-coordinate)
    projected_coords = [(coord[0], coord[1]) for _, _, coord in ca_coords]

    # Calculate the area using Shoelace's formula
    n = len(projected_coords)
    area = 0

    for i in range(n):
        x1, y1 = projected_coords[i]
        x2, y2 = projected_coords[(i + 1) % n]  # Connects the last vertex back to the first
        area += x1 * y2 - x2 * y1

    area = abs(area) / 2.0

    print(f"\nArea formed by the residues: {area:.2f} square Å" + (f" in chain {chain}" if chain else ""))
    
    return residue_pairs_distances, area


def check_residue_atom_existence(pdb_file, distance_pairs, chain=None):
    """
    Checks if all residues and their specified atoms (CA for CA distances, extreme side chain
    atoms for shortest distances) exist in the PDB file for the given chain.

    Parameters:
    pdb_file (str): The path to the PDB file.
    distance_pairs (list of tuples): A list of tuples, where each tuple represents a pair
                                      for distance calculation.
                                      Format for CA distance: ('RES1_TYPE', RES1_NUM, 'RES2_TYPE', RES2_NUM)
                                      Format for shortest distance: ('RES1_TYPE', RES1_NUM, 'RES2_TYPE', RES2_NUM)
    chain (str): The chain identifier (optional).

    Returns:
    bool: True if all residues and atoms exist, False otherwise.
    """
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found.")
        return False

    lines = parse_pdb(pdb_file)
    all_ok = True

    print(f"\n--- Running PDB existence check for {pdb_file} (Chain: {chain if chain else 'All'}) ---")

    # Keep track of checked residues/atoms to avoid redundant messages
    checked_items = set()

    for pair in distance_pairs:
        # Assuming the input format for distance_pairs is consistent with
        # measure_ca_distances and measure_shortest_distance, which take
        # (res1_type, res1_num, res2_type, res2_num).
        # We need to distinguish if it's a CA distance check or a shortest distance check.
        # For simplicity, we'll check both CA and extreme atoms if the function is generic.
        # A more robust solution might require passing a 'type' for each pair.
        # For now, let's check for CA and extreme atoms for each residue.

        res1_type, res1_num, res2_type, res2_num = pair

        # Check for first residue's CA atom
        item_key_1_ca = (res1_type, res1_num, 'CA', chain)
        if item_key_1_ca not in checked_items:
            coord1_ca = get_atom_coordinates(lines, res1_type, res1_num, 'CA', chain)
            if coord1_ca is None:
                print(f"Error: Residue {res1_type}{res1_num} CA atom not found in PDB file '{pdb_file}'" +
                      (f" for chain '{chain}'" if chain else ""))
                all_ok = False
            checked_items.add(item_key_1_ca)

        # Check for second residue's CA atom
        item_key_2_ca = (res2_type, res2_num, 'CA', chain)
        if item_key_2_ca not in checked_items:
            coord2_ca = get_atom_coordinates(lines, res2_type, res2_num, 'CA', chain)
            if coord2_ca is None:
                print(f"Error: Residue {res2_type}{res2_num} CA atom not found in PDB file '{pdb_file}'" +
                      (f" for chain '{chain}'" if chain else ""))
                all_ok = False
            checked_items.add(item_key_2_ca)

        # Check for first residue's extreme side chain atoms
        if res1_type in extreme_sidechain_atoms:
            for atom1 in extreme_sidechain_atoms[res1_type]:
                item_key_1_extreme = (res1_type, res1_num, atom1, chain)
                if item_key_1_extreme not in checked_items:
                    coord1_extreme = get_atom_coordinates(lines, res1_type, res1_num, atom1, chain)
                    if coord1_extreme is None:
                        print(f"Warning: Residue {res1_type}{res1_num} extreme atom '{atom1}' not found in PDB file '{pdb_file}'" +
                              (f" for chain '{chain}'" if chain else ""))
                        # Note: This is a warning, as shortest distance might still work with other extreme atoms
                    checked_items.add(item_key_1_extreme)
        else:
            print(f"Warning: Extreme side chain atoms not defined for residue type {res1_type}.")


        # Check for second residue's extreme side chain atoms
        if res2_type in extreme_sidechain_atoms:
            for atom2 in extreme_sidechain_atoms[res2_type]:
                item_key_2_extreme = (res2_type, res2_num, atom2, chain)
                if item_key_2_extreme not in checked_items:
                    coord2_extreme = get_atom_coordinates(lines, res2_type, res2_num, atom2, chain)
                    if coord2_extreme is None:
                        print(f"Warning: Residue {res2_type}{res2_num} extreme atom '{atom2}' not found in PDB file '{pdb_file}'" +
                              (f" for chain '{chain}'" if chain else ""))
                        # Note: This is a warning
                    checked_items.add(item_key_2_extreme)
        else:
            print(f"Warning: Extreme side chain atoms not defined for residue type {res2_type}.")


    if all_ok:
        print(f"All specified residues and CA atoms found in '{pdb_file}'. Proceeding with calculations.")
    else:
        print(f"Errors detected in '{pdb_file}'. Please revise your input or PDB file.")

    return all_ok