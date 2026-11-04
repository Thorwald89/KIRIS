"""
KIRIS - Deterministic Clinical Rules Engine for KIR-HLA Alloreactivity
Module: kir_analysis.py
Compliance: ISO 13485 (Deterministic calculations, no stochastic/probabilistic models)
"""

import math
import re
import pandas as pd

# Standard ISO 13485 e Pattern IMGT/HLA
HLA_REGEX = re.compile(r'^[A-C]\*[0-9]{2,3}:[0-9]{2,3}$')

C1_ALLELES_PREFIXES = {'01', '03', '07', '08', '12', '14', '16'}
C2_ALLELES_PREFIXES = {'02', '04', '05', '06', '09', '10', '15', '17', '18'}
A3_A11_ALLELES = {'03', '11'}

# Affinità Rheostat Model per HLA-Bw4 (Classificazione dettagliata dei residui in posizione 80)
Bw4_I80_PREFIXES = {'05', '17', '27', '37', '47', '49', '51', '52', '53', '57', '58', '59', '63', '77', '23', '24', '25', '32', '09'}
Bw4_T80_PREFIXES = {'38', '44'}

CEN_GENES = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR2DS3"]
TEL_GENES = ["KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DL5B", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"]
FRAMEWORK_GENES = {"KIR3DL3", "KIR2DL4", "KIR3DL2", "KIR3DP1", "KIR2DP1"}

CEN_B_MARKERS = {"KIR2DL2", "KIR2DS2", "KIR2DS3", "KIR2DL5B"}
TEL_B_MARKERS = {"KIR2DS1", "KIR3DS1", "KIR2DL5A", "KIR2DS5"}

def validate_and_extract_alleles(input_string, locus_name):
    """
    Validazione e normalizzazione deterministica degli alleli HLA secondo standard IMGT/HLA (ISO 13485).
    Supporta formati standard (es. A*03:01) e normalizza formati frequenti (es. 03:01, A03:01).
    """
    if not input_string or not str(input_string).strip():
        raise ValueError(f"ERRORE DI VALIDAZIONE ISO 13485 ({locus_name}): Campo vuoto.")

    locus_char = locus_name.strip()[-1].upper()
    raw_alleles = [x.strip().upper() for x in str(input_string).split(',') if x.strip()]
    validated_alleles = []

    for allele in raw_alleles:
        clean = allele.replace(" ", "")
        
        # Auto-normalizzazione se manca il prefisso del locus (es: '03:01' -> 'A*03:01')
        if not clean.startswith(('A', 'B', 'C')):
            clean = f"{locus_char}*{clean}"
        # Se manca l'asterisco (es: 'A03:01' -> 'A*03:01')
        elif len(clean) > 1 and clean[0] in ('A', 'B', 'C') and clean[1] != '*':
            clean = f"{clean[0]}*{clean[1:]}"

        if not HLA_REGEX.match(clean):
            raise ValueError(
                f"ERRORE DI VALIDAZIONE ISO 13485 ({locus_name}): "
                f"L'allele '{allele}' (normalizzato in '{clean}') non rispetta lo standard IMGT/HLA (es: {locus_char}*03:01)."
            )
        if not clean.startswith(locus_char):
            raise ValueError(
                f"ERRORE DI COERENZA GENICA ({locus_name}): "
                f"L'allele '{clean}' non appartiene al Locus {locus_char}."
            )
        validated_alleles.append(clean)
        
    if not validated_alleles:
        raise ValueError(f"ERRORE DI VALIDAZIONE ISO 13485 ({locus_name}): Nessun allele valido inserito.")
        
    return validated_alleles

def check_kir_biological_coherence(kir_set):
    """
    Controllo di coerenza biologica dei geni KIR selezionati.
    Verifica la presenza dei geni framework e la presenza degli ancoraggi per alleli dell'aplotipo B.
    """
    warnings = []
    mancanti_fw = FRAMEWORK_GENES.difference(kir_set)
    if mancanti_fw:
        warnings.append(f"⚠️ ANOMALIA GENE FRAMEWORK: Deselezionati geni costitutivi: {', '.join(sorted(mancanti_fw))}.")

    attivatori_b_specifici = {"KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR3DS1"}
    if attivatori_b_specifici.intersection(kir_set) and not ({"KIR2DL2", "KIR2DL5A", "KIR2DL5B"} & set(kir_set)):
        warnings.append("⚠️ INCOERENZA DI APLOTIPO: Geni attivatori B senza geni di ancoraggio dell'aplotipo B.")
    return warnings

def extract_allele_number(allele):
    """Estrae il gruppo sierologico/due campi numerici da una stringa allelica."""
    try:
        return allele.split('*')[1].split(':')[0].zfill(2)
    except:
        return None

def map_hla_to_ligands(hla_c, hla_b, hla_a):
    """
    Mappatura deterministica dagli alleli HLA ai rispettivi ligandi ed epitopi KIR.
    Classifica C1/C2 (incluso C*16:01/04 come C1 e C*16:02 come C2),
    Bw4 con residuo in posizione 80 (I80, T80, L80) e Bw6,
    nonché gli alleli HLA-A con epitopi Bw4 o A3/A11.
    """
    ligs = set()
    
    # Aggiornamento immunogenetico C*16: C*16:01 e C*16:04 sono C1, C*16:02 è C2
    C2_SPECIFIC = {'C*16:02', 'C*12:04'}
    C1_SPECIFIC = {'C*16:01', 'C*16:04', 'C*12:02', 'C*12:03'}

    for a in hla_c:
        if any(exc in a for exc in C2_SPECIFIC):
            ligs.add('C2')
        elif any(exc in a for exc in C1_SPECIFIC):
            ligs.add('C1')
        else:
            p = extract_allele_number(a)
            if p in C1_ALLELES_PREFIXES:
                ligs.add('C1')
            elif p in C2_ALLELES_PREFIXES:
                ligs.add('C2')

    Bw4_EXCEPTIONS = {'B*08:02', 'B*08:03', 'B*15:13', 'B*15:16', 'B*15:17', 'B*35:87'}
    Bw6_EXCEPTIONS = {'B*27:08', 'B*44:06'}
    B13_ANERGIC_EXCEPTIONS = {'B*13:01', 'B*13:02'}

    for a in hla_b:
        if any(exc in a for exc in B13_ANERGIC_EXCEPTIONS): 
            ligs.add('Bw4-L80')  # Eccezione Leucina-80 (Anergica per KIR3DL1)
        elif any(exc in a for exc in Bw4_EXCEPTIONS): 
            ligs.add('Bw4-I80')  # Isoleucina-80 (Alta affinità)
        elif any(exc in a for exc in Bw6_EXCEPTIONS): 
            ligs.add('Bw6')
        else:
            p = extract_allele_number(a)
            if p in Bw4_I80_PREFIXES:
                ligs.add('Bw4-I80')
            elif p in Bw4_T80_PREFIXES:
                ligs.add('Bw4-T80')  # Treonina-80 (Media affinità)
            else:
                ligs.add('Bw6')

    BW4_A_ALLELES = {'A*23:01', 'A*24:02', 'A*24:03', 'A*25:01', 'A*32:01'}
    for a in hla_a:
        p = extract_allele_number(a)
        if p in A3_A11_ALLELES:
            ligs.add('A3/A11')
        if any(b_allele in a for b_allele in BW4_A_ALLELES):
            ligs.add('Bw4-I80')

    return ligs

def calculate_b_content_logic(kir_set):
    """
    Calcolo dell'assetto aplotipico KIR e del B-Content Score (Cooley Model).
    Valuta indipendentemente i marker Centromerici (CEN) e Telomerici (TEL).
    """
    has_cen_b = any(g in kir_set for g in CEN_B_MARKERS)
    has_tel_b = any(g in kir_set for g in TEL_B_MARKERS)

    if not has_cen_b:
        cen_h, cen_pts = "A/A", 0
    elif "KIR2DL3" in kir_set:
        cen_h, cen_pts = "A/B", 1
    else:
        cen_h, cen_pts = "B/B", 2

    if not has_tel_b:
        tel_h, tel_pts = "A/A", 0
    elif "KIR3DL1" in kir_set:
        tel_h, tel_pts = "A/B", 1
    else:
        tel_h, tel_pts = "B/B", 2

    b_score = cen_pts + tel_pts
    return ("A/A" if b_score == 0 else "B/X"), b_score, cen_h, tel_h

def assess_donor_education(don_kir, don_ligs, don_b_alleles, don_c_alleles, kir3dl1_null_variant=False):
    """
    Calcolo dell'educazione funzionale NK (Licensing) basato sul Rheostat Model.
    Valuta con precisione:
    - KIR3DL1/Bw4 (I80 alta affinità, T80 media affinità, L80 anergica, variante Null *004)
    - KIR2DL1/C2 (alta affinità)
    - KIR2DL2/C1 (media affinità per C1, cross-reattività debole per C2)
    - KIR2DL3/C1 (bassa affinità)
    - KIR2DS1 (Attivatorio):
        * C1/C1: Educato (Forte - Hyperexcitable)
        * C1/C2: Educato (Intermedio - Responsivo)
        * C2/C2: Non Educato (Anergico / Tollerizzato)
    - KIR3DL2/A3-A11
    """
    res = []
    has_b13_anergic = any(x in don_b_alleles for x in ['B*13:01', 'B*13:02'])

    # 1. KIR3DL1 / Bw4 System
    if "KIR3DL1" in don_kir:
        if kir3dl1_null_variant:
            edu = "Non Espresso (ER-trapped / Null - *004)"
            lig_detail = "Bw4 (Null Variant)"
        elif "Bw4-I80" in don_ligs:
            edu = "Educato (Forte - High-Affinity Isoleucina-80)"
            lig_detail = "Bw4 (Isoleucina-80 / I80)"
        elif "Bw4-T80" in don_ligs:
            edu = "Educato (Moderato - Interm-Affinity Treonina-80)"
            lig_detail = "Bw4 (Treonina-80 / T80)"
        elif has_b13_anergic or "Bw4-L80" in don_ligs:
            edu = "Non Educato (Anergico - B*13 Leucina-80)"
            lig_detail = "Bw4 (Leucina-80 / L80 - B*13)"
        else:
            edu = "Non Educato (Anergico)"
            lig_detail = "Assente / Bw6"
        res.append({'KIR': 'KIR3DL1', 'Lig': lig_detail, 'Status': 'Pos (*004 Null)' if kir3dl1_null_variant else 'Pos', 'Edu': edu})
    else:
        res.append({'KIR': 'KIR3DL1', 'Lig': 'Bw4', 'Status': 'Neg', 'Edu': 'N/A'})

    # 2. KIR2DL1 / C2 System
    if "KIR2DL1" in don_kir:
        edu = "Educato (Forte - High-Affinity C2)" if "C2" in don_ligs else "Non Educato (Anergico)"
        res.append({'KIR': 'KIR2DL1', 'Lig': 'C2', 'Status': 'Pos', 'Edu': edu})
    else:
        res.append({'KIR': 'KIR2DL1', 'Lig': 'C2', 'Status': 'Neg', 'Edu': 'N/A'})

    # 3. KIR2DL2 / C1 System
    if "KIR2DL2" in don_kir:
        if "C1" in don_ligs:
            edu = "Educato (Moderato - C1)"
        elif "C2" in don_ligs:
            edu = "Educato (Debole - Cross-reactivity C2)"
        else:
            edu = "Non Educato (Anergico)"
        res.append({'KIR': 'KIR2DL2', 'Lig': 'C1', 'Status': 'Pos', 'Edu': edu})
    else:
        res.append({'KIR': 'KIR2DL2', 'Lig': 'C1', 'Status': 'Neg', 'Edu': 'N/A'})

    # 4. KIR2DL3 / C1 System
    if "KIR2DL3" in don_kir:
        edu = "Educato (Debole - Low-Affinity C1)" if "C1" in don_ligs else "Non Educato (Anergico)"
        res.append({'KIR': 'KIR2DL3', 'Lig': 'C1', 'Status': 'Pos', 'Edu': edu})
    else:
        res.append({'KIR': 'KIR2DL3', 'Lig': 'C1', 'Status': 'Neg', 'Edu': 'N/A'})

    # 5. Bw6 (Neutro)
    has_bw6 = "Bw6" in don_ligs
    res.append({'KIR': '— (No KIR)', 'Lig': 'Bw6', 'Status': "Presente" if has_bw6 else "Assente", 'Edu': "N/A (Neutro)"})

    # 6. KIR2DS1 (Attivatorio - Regole di Educazione basate su dosaggio C1/C2)
    k2ds1_pos = "KIR2DS1" in don_kir
    if not k2ds1_pos:
        res.append({'KIR': 'KIR2DS1', 'Lig': 'C1/C2', 'Status': "Neg", 'Edu': "N/A"})
    else:
        don_c_ligs = set()
        for a in don_c_alleles:
            if any(exc in a for exc in ['C*16:02', 'C*12:04']):
                don_c_ligs.add('C2')
            elif any(exc in a for exc in ['C*16:01', 'C*16:04', 'C*12:02', 'C*12:03']):
                don_c_ligs.add('C1')
            else:
                p = extract_allele_number(a)
                if p in C1_ALLELES_PREFIXES:
                    don_c_ligs.add('C1')
                elif p in C2_ALLELES_PREFIXES:
                    don_c_ligs.add('C2')

        has_c1 = 'C1' in don_c_ligs
        has_c2 = 'C2' in don_c_ligs

        # Regola Immunogenetica: C1/C1 -> Hyperexcitable; C1/C2 -> Responsivo; C2/C2 -> Anergico/Tollerizzato
        if has_c1 and not has_c2:
            k2ds1_label = "C1/C1"
            k2ds1_edu = "Educato (Forte - Hyperexcitable C1/C1)"
        elif has_c1 and has_c2:
            k2ds1_label = "C1/C2"
            k2ds1_edu = "Educato (Intermedio - Responsivo C1/C2)"
        elif has_c2 and not has_c1:
            k2ds1_label = "C2/C2"
            k2ds1_edu = "Non Educato (Anergico / Tollerizzato C2/C2)"
        else:
            k2ds1_label = "Indeterminato"
            k2ds1_edu = "Non Educato (Anergico)"
            
        res.append({'KIR': 'KIR2DS1', 'Lig': f"Locus C ({k2ds1_label})", 'Status': "Pos", 'Edu': k2ds1_edu})

    # 7. KIR3DL2 (Presenza di ligando)
    k3dl2_pos = "KIR3DL2" in don_kir
    has_a3_a11 = "A3/A11" in don_ligs
    if not k3dl2_pos:
        k3dl2_edu = "N/A"
    elif has_a3_a11:
        k3dl2_edu = "Presenza di ligando (A*03/A*11)"
    else:
        k3dl2_edu = "Assenza di ligando specifico"
        
    res.append({'KIR': 'KIR3DL2', 'Lig': 'A3/A11', 'Status': "Pos" if k3dl2_pos else "Neg", 'Edu': k3dl2_edu})

    return res

def calculate_vectors_quantitative(edu_results, don_ligs, rec_ligs, haplo, b_score, k2ds1_edu):
    """
    Calcolo quantitativo dei vettori clinici (Score 1-5):
    - GvH (Alloreattività Donatore -> Ricevente basata su Missing-Self inibitorio ponderato + KIR2DS1 attivatorio)
    - HvG (Rischio Rigetto Graft da parte delle cellule NK dell'ospite)
    - GvL (Effetto Anti-Leucemico guidato da missing ligand, KIR2DS1 attivo su ricevente C2+, e B-score)
    """
    rec_has_c1 = 'C1' in rec_ligs
    rec_has_c2 = 'C2' in rec_ligs
    rec_has_bw4 = ('Bw4-I80' in rec_ligs) or ('Bw4-T80' in rec_ligs)
    rec_has_a3_a11 = 'A3/A11' in rec_ligs

    don_has_c1 = 'C1' in don_ligs
    don_has_c2 = 'C2' in don_ligs
    don_has_bw4 = ('Bw4-I80' in don_ligs) or ('Bw4-T80' in don_ligs)

    edu_map = {r['KIR']: r for r in edu_results}

    # ----------------------------------------------------
    # 1. VETTORE GvH (Alloreattività NK Donatore-vs-Ospite)
    # ----------------------------------------------------
    gvh_tags = []
    gvh_weight = 0.0

    r_3dl1 = edu_map.get('KIR3DL1')
    if r_3dl1 and "Educato" in r_3dl1['Edu']:
        if not rec_has_bw4:
            if "Isoleucina-80" in r_3dl1['Edu']:
                gvh_weight += 2.0
                gvh_tags.append("KIR3DL1(Bw4-I80) ML")
            else:
                gvh_weight += 1.2
                gvh_tags.append("KIR3DL1(Bw4-T80) ML")

    r_2dl1 = edu_map.get('KIR2DL1')
    if r_2dl1 and "Educato" in r_2dl1['Edu']:
        if not rec_has_c2:
            gvh_weight += 2.0
            gvh_tags.append("KIR2DL1(C2) ML")

    r_2dl2 = edu_map.get('KIR2DL2')
    r_2dl3 = edu_map.get('KIR2DL3')
    c1_missing_in_rec = not rec_has_c1

    if c1_missing_in_rec:
        if r_2dl2 and "Educato (Moderato" in r_2dl2['Edu']:
            gvh_weight += 1.2
            gvh_tags.append("KIR2DL2(C1) ML")
        elif r_2dl3 and "Educato" in r_2dl3['Edu']:
            gvh_weight += 0.8
            gvh_tags.append("KIR2DL3(C1) ML")

    r_3dl2 = edu_map.get('KIR3DL2')
    if r_3dl2 and "Presenza di ligando" in r_3dl2['Edu'] and not rec_has_a3_a11:
        gvh_weight += 0.4
        gvh_tags.append("KIR3DL2(A3/A11) ML")

    if rec_has_c2:
        if "Hyperexcitable" in k2ds1_edu:
            gvh_weight += 1.0
            gvh_tags.append("KIR2DS1(Iperecc. vs C2)")
        elif "Responsivo" in k2ds1_edu:
            gvh_weight += 0.5
            gvh_tags.append("KIR2DS1(Resp. vs C2)")

    if gvh_weight == 0.0:
        gvh_score = 1
        gvh_status = "Nullo / Tollerante (Nessun Mismatch)"
        gvh_color = "#27ae60"
    elif gvh_weight < 1.0:
        gvh_score = 2
        gvh_status = f"Basso ({', '.join(gvh_tags)})"
        gvh_color = "#f1c40f"
    elif gvh_weight <= 2.2:
        gvh_score = 3
        gvh_status = f"Moderato ({', '.join(gvh_tags)})"
        gvh_color = "#e67e22"
    elif gvh_weight <= 3.6:
        gvh_score = 4
        gvh_status = f"Alto Potenziale ({', '.join(gvh_tags)})"
        gvh_color = "#e74c3c"
    else:
        gvh_score = 5
        gvh_status = f"Massima Alloreattività NK ({', '.join(gvh_tags)})"
        gvh_color = "#c0392b"

    # ----------------------------------------------------
    # 2. VETTORE HvG (Host-versus-Graft / Rischio Rigetto)
    # ----------------------------------------------------
    hvg_tags = []
    hvg_weight = 0.0

    if rec_has_c2 and not don_has_c2:
        hvg_weight += 2.0
        hvg_tags.append("Manca C2")
    if rec_has_bw4 and not don_has_bw4:
        hvg_weight += 2.0
        hvg_tags.append("Manca Bw4")
    if rec_has_c1 and not don_has_c1:
        hvg_weight += 1.2
        hvg_tags.append("Manca C1")

    if hvg_weight == 0.0:
        hvg_score = 1
        hvg_status = "Nullo (Graft protetto)"
        hvg_color = "#27ae60"
    elif hvg_weight <= 1.2:
        hvg_score = 2
        hvg_status = f"Basso ({', '.join(hvg_tags)})"
        hvg_color = "#f1c40f"
    elif hvg_weight <= 2.0:
        hvg_score = 3
        hvg_status = f"Moderato ({', '.join(hvg_tags)})"
        hvg_color = "#e67e22"
    elif hvg_weight <= 3.5:
        hvg_score = 4
        hvg_status = f"Alto Rischio ({', '.join(hvg_tags)})"
        hvg_color = "#e74c3c"
    else:
        hvg_score = 5
        hvg_status = f"Massimo Rischio Rigetto ({', '.join(hvg_tags)})"
        hvg_color = "#c0392b"

    # ----------------------------------------------------
    # 3. VETTORE GvL (Graft-versus-Leukemia / Effetto Anti-Leucemico)
    # ----------------------------------------------------
    gvl_points = 1.0
    gvl_reasons = []

    inhibitory_mismatches = [t for t in gvh_tags if "ML" in t]
    has_high_aff_ml = any(("C2" in t or "Bw4-I80" in t) for t in inhibitory_mismatches)

    if has_high_aff_ml:
        gvl_points += 1.5
        gvl_reasons.append("Missing-Self Alta Affinità (+1.5)")
        if len(inhibitory_mismatches) > 1:
            gvl_points += 0.5
            gvl_reasons.append("Multiplo Missing-Self (+0.5)")
    elif inhibitory_mismatches:
        gvl_points += 1.0
        gvl_reasons.append("Missing-Self Media/Bassa Affinità (+1.0)")

    if rec_has_c2:
        if "Hyperexcitable" in k2ds1_edu:
            gvl_points += 1.5
            gvl_reasons.append("KIR2DS1 Ipereccitabile vs C2 (+1.5)")
        elif "Responsivo" in k2ds1_edu:
            gvl_points += 1.0
            gvl_reasons.append("KIR2DS1 Responsivo vs C2 (+1.0)")
    elif "Hyperexcitable" in k2ds1_edu or "Responsivo" in k2ds1_edu:
        gvl_reasons.append("KIR2DS1 non stimolato (Ricevente C2-negativo)")

    if b_score >= 3:
        gvl_points += 1.5
        gvl_reasons.append("Aplotipo B Arricchito (B-Score >= 3) (+1.5)")
    elif b_score == 2:
        gvl_points += 1.0
        gvl_reasons.append("Aplotipo B Favorevole (B-Score 2) (+1.0)")
    elif b_score == 1:
        gvl_points += 0.5
        gvl_reasons.append("Aplotipo B (B-Score 1) (+0.5)")

    gvl_score = min(5, max(1, round(gvl_points)))
    gvl_labels = {
        1: "Minimo",
        2: "Basso",
        3: "Moderato",
        4: "Alto Potenziale",
        5: "Massimo Effetto GvL"
    }
    reasons_str = ", ".join(gvl_reasons) if gvl_reasons else "Effetto NK basale"
    gvl_status = f"{gvl_labels[gvl_score]} [{reasons_str}]"
    gvl_color = "#27ae60" if gvl_score >= 4 else "#2980b9" if gvl_score == 3 else "#7f8c8d"

    return gvh_score, gvh_status, gvh_color, hvg_score, hvg_status, hvg_color, gvl_score, gvl_status, gvl_color

def calculate_donor_recommendation(gvh_score, hvg_score, gvl_score, haplo, b_score, is_dark=True, **kwargs):
    """
    Sintetizzatore decisionale clinico per selezione donatore (Criteri EBMT / GITMO).
    Rilascia una classificazione a 4 livelli con raccomandazione esplicita e palette per dark/light mode.
    """
    if hvg_score >= 4:
        level = "SUBOTTIMALE"
        title = "⚠️ Donatore Subottimale (Rischio Rigetto HvG Elevato)"
        color = "#ef4444"
        badge_bg = "rgba(239, 68, 68, 0.15)" if is_dark else "#fee2e2"
        badge_txt = "#fca5a5" if is_dark else "#991b1b"
        recom = (
            "Attenzione: elevata reattività delle cellule NK residue dell'ospite contro le cellule staminali "
            "del donatore (mismatch inibitorio inverso significativo). Rischio aumentato di mancato attecchimento "
            "(graft failure). Se selezionato, valutare intensificazione del condizionamento o immuno-profilassi rinforzata."
        )
    elif gvl_score >= 4 and hvg_score <= 2:
        level = "PRIORITÀ ELEVATA"
        title = "⭐ Donatore Ottimale (Massimo Beneficio GvL & Innesto Protetto)"
        color = "#10b981"
        badge_bg = "rgba(16, 185, 129, 0.15)" if is_dark else "#d1fae5"
        badge_txt = "#6ee7b7" if is_dark else "#065f46"
        recom = (
            "Profilo immunogenetico ideale: forte potenziale citotossico alloreattivo contro la leucemia "
            "(elevato Missing-Self e/o stimolo attivatorio KIR2DS1 su target C2+) associato a un innesto "
            "altamente protetto da rigetto mediato da NK. Candidato di prima scelta."
        )
    elif gvl_score >= 3 and hvg_score <= 3:
        level = "FAVOREVOLE"
        title = "🟢 Donatore Favorevole (Buon Potenziale GvL)"
        color = "#3b82f6"
        badge_bg = "rgba(59, 130, 246, 0.15)" if is_dark else "#dbeafe"
        badge_txt = "#93c5fd" if is_dark else "#1e40af"
        recom = (
            "Profilo clinico vantaggioso: presenza di fattori alloreattivi e/o aplotipo B favorevole, con "
            "rischio di reiezione contenuto. Donatore raccomandato per la procedura trapiantologica."
        )
    else:
        level = "STANDARD"
        title = "🟡 Donatore Idoneo Standard (Effetto NK Neutro)"
        color = "#f59e0b"
        badge_bg = "rgba(245, 158, 11, 0.15)" if is_dark else "#fef3c7"
        badge_txt = "#fcd34d" if is_dark else "#92400e"
        recom = (
            "Profilo convenzionale neutro: assenza di specifici mismatch alloreattivi inibitori e innesto "
            "in equilibrio immunologico basale. Il trapianto seguirà il decorso standard dei protocolli clinici."
        )

    return {
        "level": level,
        "title": title,
        "color": color,
        "badge_bg": badge_bg,
        "badge_txt": badge_txt,
        "recommendation": recom
    }

def generate_radar_chart_svg(gvh_score, hvg_score, gvl_score, is_dark=True, width=380, height=310, **kwargs):
    """
    Genera un Radar Chart SVG professionale a 3 assi (GvH, HvG, GvL) su scala 1-5.
    Costruito come stringa compatta priva di interruzioni Markdown per garantire rendering perfetto.
    Supporta Dark Mode ad alto contrasto e Light Mode.
    """
    cx = width / 2.0
    cy = 142.0
    max_r = 95.0

    # 3 Assi a 120 gradi: Top (-90°), Bottom-Right (30°), Bottom-Left (150°)
    angles_deg = [-90.0, 30.0, 150.0]
    angles_rad = [math.radians(a) for a in angles_deg]

    def p2c(angle_rad, radius):
        return cx + radius * math.cos(angle_rad), cy + radius * math.sin(angle_rad)

    # Palette tema
    bg_color = "#0f172a" if is_dark else "#f8fafc"
    card_border = "#334155" if is_dark else "#e2e8f0"
    grid_stroke = "#334155" if is_dark else "#cbd5e1"
    axis_stroke = "#475569" if is_dark else "#94a3b8"
    grid_txt_color = "#94a3b8" if is_dark else "#64748b"
    alt_fill = "rgba(30, 41, 59, 0.4)" if is_dark else "rgba(241, 245, 249, 0.8)"

    # Poligoni concentrici di sfondo da 5 a 1
    grid_elements = []
    for step in range(5, 0, -1):
        r = (step / 5.0) * max_r
        pts = [p2c(a, r) for a in angles_rad]
        pts_str = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        fill = alt_fill if step % 2 == 1 else "none"
        grid_elements.append(
            f'<polygon points="{pts_str}" fill="{fill}" stroke="{grid_stroke}" stroke-width="1.2" stroke-dasharray="3,3"/>'
        )
        lx, ly = p2c(angles_rad[0], r)
        grid_elements.append(
            f'<text x="{lx + 6:.1f}" y="{ly + 4:.1f}" font-size="9" font-weight="600" fill="{grid_txt_color}" font-family="sans-serif">{step}</text>'
        )

    # Assi radiali con terminali a pallino
    axis_elements = []
    for a in angles_rad:
        x, y = p2c(a, max_r)
        axis_elements.append(
            f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{x:.1f}" y2="{y:.1f}" stroke="{axis_stroke}" stroke-width="1.5"/>'
        )
        axis_elements.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2" fill="{axis_stroke}"/>'
        )

    # Poligono dati
    scores = [float(gvh_score), float(gvl_score), float(hvg_score)]
    data_pts = []
    for sc, a in zip(scores, angles_rad):
        r = (max(1.0, min(5.0, sc)) / 5.0) * max_r
        data_pts.append(p2c(a, r))

    data_pts_str = " ".join(f"{x:.1f},{y:.1f}" for x, y in data_pts)

    if gvl_score >= 4 and hvg_score <= 2:
        poly_fill = "rgba(16, 185, 129, 0.45)"
        poly_stroke = "#10b981"
        glow_color = "rgba(16, 185, 129, 0.3)"
    elif hvg_score >= 4:
        poly_fill = "rgba(239, 68, 68, 0.45)"
        poly_stroke = "#ef4444"
        glow_color = "rgba(239, 68, 68, 0.3)"
    else:
        poly_fill = "rgba(59, 130, 246, 0.45)"
        poly_stroke = "#3b82f6"
        glow_color = "rgba(59, 130, 246, 0.3)"

    data_polygon = (
        f'<polygon points="{data_pts_str}" fill="{poly_fill}" stroke="{poly_stroke}" '
        f'stroke-width="3" stroke-linejoin="round" style="filter: drop-shadow(0 0 6px {glow_color});"/>'
    )

    data_dots = []
    dot_inner = "#0f172a" if is_dark else "#ffffff"
    for x, y in data_pts:
        data_dots.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="{poly_stroke}" stroke="{dot_inner}" stroke-width="2"/>'
        )

    # Badge ed etichette esterne per i 3 vettori
    # Top: GvH
    gvh_x, gvh_y = p2c(angles_rad[0], max_r + 26)
    # Bottom-Right: GvL
    gvl_x, gvl_y = p2c(angles_rad[1], max_r + 28)
    # Bottom-Left: HvG
    hvg_x, hvg_y = p2c(angles_rad[2], max_r + 28)

    badge_bg = "#1e293b" if is_dark else "#ffffff"
    badge_border = "#475569" if is_dark else "#cbd5e1"

    txt_gvh = "#f8fafc" if is_dark else "#1e293b"
    txt_gvl = "#34d399" if is_dark else "#059669"
    txt_hvg = "#f87171" if is_dark else "#dc2626"

    labels = [
        f'<g transform="translate({gvh_x:.1f},{gvh_y:.1f})">'
        f'<rect x="-65" y="-14" width="130" height="24" rx="12" fill="{badge_bg}" stroke="{badge_border}" stroke-width="1"/>'
        f'<text x="0" y="2" text-anchor="middle" dominant-baseline="middle" font-size="11" font-weight="700" fill="{txt_gvh}" font-family="system-ui, sans-serif">⚡ GvH: {gvh_score}/5</text>'
        f'</g>',
        f'<g transform="translate({gvl_x:.1f},{gvl_y:.1f})">'
        f'<rect x="-65" y="-14" width="130" height="24" rx="12" fill="{badge_bg}" stroke="{badge_border}" stroke-width="1"/>'
        f'<text x="0" y="2" text-anchor="middle" dominant-baseline="middle" font-size="11" font-weight="700" fill="{txt_gvl}" font-family="system-ui, sans-serif">🎯 GvL: {gvl_score}/5</text>'
        f'</g>',
        f'<g transform="translate({hvg_x:.1f},{hvg_y:.1f})">'
        f'<rect x="-65" y="-14" width="130" height="24" rx="12" fill="{badge_bg}" stroke="{badge_border}" stroke-width="1"/>'
        f'<text x="0" y="2" text-anchor="middle" dominant-baseline="middle" font-size="11" font-weight="700" fill="{txt_hvg}" font-family="system-ui, sans-serif">🛡️ HvG: {hvg_score}/5</text>'
        f'</g>'
    ]

    # Assemblaggio compatto senza commenti HTML e senza newline multiple
    svg_body = "".join(grid_elements) + "".join(axis_elements) + data_polygon + "".join(data_dots) + "".join(labels)
    svg_clean = (
        f'<svg width="100%" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'style="background:{bg_color}; border-radius:12px; border:1px solid {card_border}; display:block; margin:auto;">'
        f'{svg_body}</svg>'
    )
    return svg_clean

def format_excel_row(edu_results, haplo, b_score, gvh_sc, hvg_sc, gvl_sc):
    """
    Formatta la riga standard per esportazione verso fogli di calcolo Excel.
    """
    def get_excel_label(kir_name):
        match_kir = next((x for x in edu_results if x['KIR'] == kir_name), None)
        if not match_kir or match_kir['Status'] == 'Neg':
            return "NEG"
        return "POS, E" if ('Educato' in match_kir['Edu'] or match_kir['Edu'] == 'Educato/Responsive') else "POS, NE"

    return (
        f"{haplo}\t{b_score}\t{get_excel_label('KIR3DL1')}\t{get_excel_label('KIR2DL1')}\t"
        f"{get_excel_label('KIR2DL2')}\t{get_excel_label('KIR2DL3')}\t{get_excel_label('KIR2DS1')}\t"
        f"{get_excel_label('KIR3DL2')}\t{gvh_sc}\t{hvg_sc}\t{gvl_sc}\tRUO (Non-IVD)"
    )

def analyze_single_pair(rec_a, rec_b, rec_c, don_a, don_b, don_c, don_kir_set, kir3dl1_null_variant=False, is_dark=True, **kwargs):
    """
    Funzione ad alta affidabilità per analizzare una coppia Ricevente-Donatore.
    Utilizzata sia dalla GUI interattiva che dal motore batch di coorte.
    """
    r_a = validate_and_extract_alleles(rec_a, "Rec HLA-A")
    r_b = validate_and_extract_alleles(rec_b, "Rec HLA-B")
    r_c = validate_and_extract_alleles(rec_c, "Rec HLA-C")
    d_a = validate_and_extract_alleles(don_a, "Don HLA-A")
    d_b = validate_and_extract_alleles(don_b, "Don HLA-B")
    d_c = validate_and_extract_alleles(don_c, "Don HLA-C")

    kir_warnings = check_kir_biological_coherence(don_kir_set)
    r_ligs = map_hla_to_ligands(r_c, r_b, r_a)
    d_ligs = map_hla_to_ligands(d_c, d_b, d_a)

    haplo, b_score, cen_h, tel_h = calculate_b_content_logic(don_kir_set)
    edu_res = assess_donor_education(don_kir_set, d_ligs, d_b, d_c, kir3dl1_null_variant=kir3dl1_null_variant)

    k2ds1_edu_status = next((x['Edu'] for x in edu_res if x['KIR'] == 'KIR2DS1'), "N/A")

    gvh_sc, gvh_s, gvh_c, hvg_sc, hvg_s, hvg_c, gvl_sc, gvl_s, gvl_c = calculate_vectors_quantitative(
        edu_res, d_ligs, r_ligs, haplo, b_score, k2ds1_edu_status
    )

    recom = calculate_donor_recommendation(gvh_sc, hvg_sc, gvl_sc, haplo, b_score, is_dark=is_dark)
    radar_svg = generate_radar_chart_svg(gvh_sc, hvg_sc, gvl_sc, is_dark=is_dark)
    excel_str = format_excel_row(edu_res, haplo, b_score, gvh_sc, hvg_sc, gvl_sc)

    return {
        "r_a": r_a, "r_b": r_b, "r_c": r_c,
        "d_a": d_a, "d_b": d_b, "d_c": d_c,
        "kir_warnings": kir_warnings,
        "r_ligs": r_ligs, "d_ligs": d_ligs,
        "haplo": haplo, "b_score": b_score,
        "cen_h": cen_h, "tel_h": tel_h,
        "edu_res": edu_res,
        "k2ds1_edu_status": k2ds1_edu_status,
        "gvh": {"score": gvh_sc, "status": gvh_s, "color": gvh_c},
        "hvg": {"score": hvg_sc, "status": hvg_s, "color": hvg_c},
        "gvl": {"score": gvl_sc, "status": gvl_s, "color": gvl_c},
        "recommendation": recom,
        "radar_svg": radar_svg,
        "excel_str": excel_str
    }
