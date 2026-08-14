"""
Chemistry Utilities Module.
Contains helper functions leveraging the RDKit library to compute
molecular properties, functional groups, and formulas from SMILES strings.
"""

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Disable unnecessary RDKit warnings (e.g., unusual valences) to keep the console clean for the user.
RDLogger.DisableLog('rdApp.*')

# OPTIMIZATION: Pre-compiling SMARTS patterns at the module level.
# This prevents RDKit from recompiling the regex-like patterns thousands of times during the database loop.
SMARTS_G1 = [Chem.MolFromSmarts(s) for s in ['[CX3]=[OX1]', '[P]=[OX1]', '[S]=[OX1]', 'N=C=O']]
SMARTS_G2 = [Chem.MolFromSmarts(s) for s in ['[NX3;!$(NC=O)]', '[NX2]=C', '[CX2]#N']]
SMARTS_G3 = [Chem.MolFromSmarts(s) for s in ['[OX2H]', '[OX2]([#6])[#6]', '[Si]-[OX2]']]

def get_mol_props(smiles):
    """
    Analyzes a SMILES string and extracts molecular properties.
    
    Args:
        smiles (str): The Canonical SMILES representation of the molecule.
        
    Returns:
        tuple: Contains the RDKit Mol object, specific atom counts (C, Cl, Br, Si, S, N),
               Group 1 presence flag, DART chemical family classification, and the chemical formula.
    """
    if pd.isna(smiles) or str(smiles).strip() == "":
        return None, 0, 0, 0, 0, 0, 0, False, "Alkane/Aromatic/Other", ""
    
    try:
        mol = Chem.MolFromSmiles(str(smiles).strip())
        if not mol:
            return None, 0, 0, 0, 0, 0, 0, False, "Alkane/Aromatic/Other", ""
        
        # We manually count atoms critical for isotopic pattern verification (M+1, M+2)
        nc = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
        ncl = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'Cl')
        nbr = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'Br')
        nsi = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'Si')
        ns = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'S')
        nn = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'N')
        
        # Check for specific functional groups using the pre-compiled SMARTS
        has_g1 = any(mol.HasSubstructMatch(pat) for pat in SMARTS_G1)
        has_g2 = any(mol.HasSubstructMatch(pat) for pat in SMARTS_G2)
        has_g3 = any(mol.HasSubstructMatch(pat) for pat in SMARTS_G3)
                  
        # Assign broad chemical families used for DART-MS ionization predictions
        fams = []
        if has_g1: fams.append("Carbonyl/P=O/S=O")
        if has_g2: fams.append("Amine/Nitrile")
        if has_g3: fams.append("Alcohol/Ether/Si-O")
        if not fams: fams.append("Alkane/Aromatic/Other")
        
        form = CalcMolFormula(mol)
        return mol, nc, ncl, nbr, nsi, ns, nn, has_g1, " + ".join(fams), form
    except:
        # Fallback to prevent the entire processing batch from crashing on a single bad SMILES
        return None, 0, 0, 0, 0, 0, 0, False, "Alkane/Aromatic/Other", ""