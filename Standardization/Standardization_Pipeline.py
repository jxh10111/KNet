from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from pathlib import Path
from AutoMID_classifier_standarization_pipeline import StructurePipeline

def standardization_pipeline():
    # Initialize the pipeline
    pipeline = StructurePipeline()

    # SMILES strings
    smiles_list = [
        "CC(=O)O.Cl",  # Acetic acid with a chloride salt
        "C[N+](C)(C)CC(=O)[O-]",  # Betaine zwitterion
        "c1ccccc1O",  # Phenol
        "C1=CC=CN=C1",  # Pyridine
        "[Na+].[Cl-]"  # Sodium chloride
    ]

    print("Testing StructurePipeline...\n")

    for smiles in smiles_list:
        print(f"Processing SMILES: {smiles}")
        try:
            structure = pipeline.process(smiles)
            print(f"  Original: {structure.S0}")
            print(f"  Canonicalized (S1): {structure.S1}")
            print(f"  Desalted (S2): {structure.S2}")
            print(f"  Neutralized (S3): {structure.S3}")
            print(f"  Tautomer Standardized (S4): {structure.S4}")
            print("-" * 50)
        except Exception as e:
            print(f"Error processing {smiles}: {e}")
            print("-" * 50)

if __name__ == "__main__":
    standardization_pipeline()
