from dataclasses import dataclass
from typing import Optional
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover, InputFormat
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass
class Structure:
    S0: str
    mol: Optional[Chem.Mol]
    S1: str = None  # canonicalised
    S2: str = None  # desalted
    S3: str = None  # fix charges
    S4: str = None  # tautomerized
    # S5: str = None  # removed isotopes
    # S6: str = None  # removed chirality
    # S7: str = None  # removed bond information
    # S8: str = None  # removed atom type information


import os

SALT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conservativeSaltLibrary.sdf')
SALT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'aggressiveSaltLibrary.sdf')


class StructurePipeline:
    def __init__(self) -> None:
        # salts = [mol for mol in Chem.SDMolSupplier(str(SALT_FILE))]

        self.remover = SaltRemover(defnFilename=str(SALT_FILE), defnFormat=InputFormat.MOL)  # defnData="[Cl,Br]")
        self.tautomerizer = rdMolStandardize.TautomerEnumerator()
        self.tautomerizer.SetMaxTautomers(256)
        self.tautomerizer.SetRemoveSp3Stereo(False)  # Keep stereo information of keto/enol tautomerization
        self.largestFragment = rdMolStandardize.LargestFragmentChooser()

    def process(self, smiles: str) -> Structure:
        structure = Structure(smiles, mol=Chem.MolFromSmiles(smiles))
        if structure.mol is None:
            raise Exception(f'SMILES {smiles} cannot be processed')
        structure = self.transformS0S1(structure)
        structure = self.transformS1S2(structure)
        structure = self.transformS2S3(structure)
        structure = self.transformS3S4(structure)
        # structure = self.transformS4S5(structure)
        # structure = self.transformS5S6(structure)
        # structure = self.transformS6S7(structure)
        # structure = self.transformS7S8(structure)
        return structure

    # Return canonical SMILES and mols from canonical SMILES
    def transformS0S1(self, structure: Structure) -> Structure:
        structure.S1 = Chem.MolToSmiles(structure.mol)
        structure.mol = Chem.MolFromSmiles(structure.S1)
        return structure

    # Remove salt
    def transformS1S2(self, structure: Structure) -> Structure:
        mol = self.remover.StripMol(structure.mol)
        s2 = Chem.MolToSmiles(mol)
        if len(s2) == 0 or '.' in s2:
            mol = self.largestFragment.choose(structure.mol)
            s2 = Chem.MolToSmiles(mol)
        structure.mol = mol
        structure.S2 = s2
        return structure

    # Neutralize atoms (ionization state)
    def transformS2S3(self, structure: Structure) -> Structure:
        structure.mol = neutralize_atoms(structure.mol)
        structure.S3 = Chem.MolToSmiles(structure.mol)
        structure.mol = Chem.MolFromSmiles(structure.S3)
        structure.S3 = Chem.MolToSmiles(structure.mol)
        return structure

    # tautomer standardizing and canonicalizing (Returns the canonical tautomer for a molecule)
    def transformS3S4(self, structure: Structure) -> Structure:
        mol = self.tautomerizer.Canonicalize(structure.mol)
        structure.S4 = Chem.MolToSmiles(mol)
        structure.mol = Chem.MolFromSmiles(structure.S4)
        structure.S4 = Chem.MolToSmiles(structure.mol)
        return structure


# ensures that molecules are properly neutralized, adjusting charges and hydrogen counts while preserving molecular integrity
def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol
