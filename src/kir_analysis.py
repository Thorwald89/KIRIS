"""
Modulo per l'analisi dei ligandi KIR basato sulle sequenze amminoacidiche HLA.
"""

class KIRLigandProfiler:
    def __init__(self):
        # Motivi amminoacidici di riferimento alla posizione 80 (HLA-C)
        self.c1_residue = 'N'  # Asparagina -> C1
        self.c2_residue = 'K'  # Lisina -> C2

    def determine_hla_c_group(self, residue_80: str) -> str:
        """
        Determina il gruppo di ligando KIR (C1 o C2) in base al residuo in posizione 80.
        """
        res = residue_80.upper().strip()
        if res == self.c1_residue:
            return "C1"
        elif res == self.c2_residue:
            return "C2"
        return "Unknown"

    def evaluate_missing_self(self, kir_genotype: list, hla_ligands: list) -> dict:
        """
        Valuta il potenziale 'Missing Self' (GvL / Reattività NK).
        """
        status = {}
        
        # Interazione KIR2DL1 -> C2
        if "KIR2DL1" in kir_genotype:
            status["KIR2DL1_C2"] = "Inhibited" if "C2" in hla_ligands else "Missing Self (GvL)"
            
        # Interazione KIR2DL2/3 -> C1
        if "KIR2DL2" in kir_genotype or "KIR2DL3" in kir_genotype:
            status["KIR2DL2/3_C1"] = "Inhibited" if "C1" in hla_ligands else "Missing Self (GvL)"
            
        return status
