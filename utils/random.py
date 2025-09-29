'''Random functions that assist with using MD to study polymer-enzyme systems'''

__author__ = 'Joe Laforet Jr.'
__email__ = 'jola3134@colorado.edu'

########################################################################################################################################################

# Takes in the smiles string of a molecule and converts it into
# a volume estimate, gives estimate in units of Nm^3
# Also prints the MW and the volume

def smiles_to_vol(smiles):
    """Calculate the volume of a molecule given its SMILES using the VdW approximation.
    Parameters:
    ----------------------------------------
    smiles: (str) SMILES from online resource
    
    Returns: (float) Volume of the molecule in cubic Nm
    """
    
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    formula = CalcMolFormula(mol)

    atomic_count = defaultdict(lambda : 0)
    
    for atom in mol.GetAtoms():
        atomic_count[atom.GetSymbol()] += 1
            
    total_vol = 0
    
    #Dictionary contains Atomic Symbol and VdW radius retrieved from Wikipedia
    atoms_to_parse = {'C':1.7 , 'H':1.2, 'N':1.55, 'O':1.52, 'P':1.8, 'Cl':1.75, 'S':1.8, 'F':1.47, 'I':1.98}
    
    for x in atoms_to_parse.keys():
        if(x not in atoms_to_parse.keys()):
            raise Exception("Error, Atom not found in VdW dictionary. Find VdW for {} and add it to dictionary.".format(x))
        
        total_vol += atomic_count[x] * (4*np.pi* (atoms_to_parse[x]**3))/3
    
    total_vol/=1000
    print("Volume of molecule is: {:.4f} Cubic Nm".format(total_vol))
    print("MW of molecule is: {:.4f} g/mol".format(Chem.Descriptors.ExactMolWt(mol)))
    
    return total_vol
