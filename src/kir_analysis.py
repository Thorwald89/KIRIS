import re

# Standard ISO e Pattern
HLA_REGEX = re.compile(r'^[A-C]\*[0-9]{2,3}:[0-9]{2,3}$')

C1_ALLELES_PREFIXES = {'01', '03', '07', '08', '12', '14', '16'}
C2_ALLELES_PREFIXES = {'02', '04', '05', '06', '09', '10', '15', '17', '18'}
A3_A11_ALLELES = {'03', '11'}

Bw4_I80_PREFIXES = {'05', '17', '27', '37', '47', '49', '51', '52', '53', '57', '58', '59', '63', '77', '23', '24', '25', '32', '09'}
Bw4_T80_PREFIXES = {'38', '44'}

CEN_GENES = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR2DS3"]
TEL_GENES = ["KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DL5B", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"]
FRAMEWORK_GENES = {"KIR3DL3", "KIR2DL4", "KIR3DL2", "KIR3DP1", "KIR2DP1"}

CEN_B_MARKERS = {"KIR2DL2", "KIR2DS2", "KIR2DS3", "KIR2DL5B"}
TEL_B_MARKERS = {"KIR2DS1", "KIR3DS1", "KIR2DL5A", "KIR2DS5"}

def validate_and_extract_alleles(input_string, locus_name):
    raw_alleles = [x.strip().upper() for x in input_string.split(',') if x.strip()]
    validated_alleles = []

    for allele in raw_alleles:
        allele_clean = allele.replace(" ", "")
        if not HLA_REGEX.match(allele_clean):
            raise ValueError(
                f"ERRORE DI VALIDAZIONE ISO 13485 ({locus_name}): "
                f"L'allele '{allele}' non rispetta lo standard IMGT/HLA (es: A*03:01)."
            )
        if not allele_clean.startswith(locus_name[-1]):
            raise ValueError(
                f"ERRORE DI COERENZA GENICA ({locus_name}): "
                f"L'allele '{allele_clean}' non appartiene al Locus indicato."
            )
        validated_alleles.append(allele_clean)
    return validated_alleles

def check_kir_biological_coherence(kir_set):
    warnings = []
    mancanti_fw = FRAMEWORK_GENES.difference(kir_set)
    if mancanti_fw:
        warnings.append(f"⚠️ ANOMALIA GENE FRAMEWORK: Deselezionati geni costitutivi: {', '.join(mancanti_fw)}.")

    attivatori_b_specifici = {"KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR3DS1"}
    if attivatori_b_specifici.intersection(kir_set) and not ({"KIR2DL2", "KIR2DL5A", "KIR2DL5B"} & set(kir_set)):
        warnings.append("⚠️ INCOERENZA DI APLOTIPO: Geni attivatori B senza geni di ancoraggio dell'aplotipo B.")
    return warnings

def extract_allele_number(allele):
    try: return allele.split('*')[1].split(':')[0].zfill(2)
    except: return None

def map_hla_to_ligands(hla_c, hla_b, hla_a):
    ligs = set()
    C2_SPECIFIC = {'C*16:02', 'C*16:04', 'C*12:04'}
    C1_SPECIFIC = {'C*16:01', 'C*12:02', 'C*12:03'}

    for a in hla_c:
        if any(exc in a for exc in C2_SPECIFIC): ligs.add('C2')
        elif any(exc in a for exc in C1_SPECIFIC): ligs.add('C1')
        else:
            p = extract_allele_number(a)
            if p in C1_ALLELES_PREFIXES: ligs.add('C1')
            elif p in C2_ALLELES_PREFIXES: ligs.add('C2')

    Bw4_EXCEPTIONS = {'B*08:02', 'B*08:03', 'B*15:13', 'B*15:16', 'B*15:17', 'B*35:87'}
    Bw6_EXCEPTIONS = {'B*27:08', 'B*44:06'}
    B13_ANERGIC_EXCEPTIONS = {'B*13:01', 'B*13:02'}

    for a in hla_b:
        if any(exc in a for exc in B13_ANERGIC_EXCEPTIONS): ligs.add('Bw6')
        elif any(exc in a for exc in Bw4_EXCEPTIONS): ligs.add('Bw4-I80')
        elif any(exc in a for exc in Bw6_EXCEPTIONS): ligs.add('Bw6')
        else:
            p = extract_allele_number(a)
            if p in Bw4_I80_PREFIXES: ligs.add('Bw4-I80')
            elif p in Bw4_T80_PREFIXES: ligs.add('Bw4-T80')
            else: ligs.add('Bw6')

    BW4_A_ALLELES = {'A*23:01', 'A*24:02', 'A*24:03', 'A*25:01', 'A*32:01'}
    for a in hla_a:
        p = extract_allele_number(a)
        if p in A3_A11_ALLELES: ligs.add('A3/A11')
        if any(b_allele in a for b_allele in BW4_A_ALLELES): ligs.add('Bw4-I80')

    return ligs

def calculate_b_content_logic(kir_set):
    has_cen_b = any(g in kir_set for g in CEN_B_MARKERS)
    has_tel_b = any(g in kir_set for g in TEL_B_MARKERS)

    if not has_cen_b: cen_h, cen_pts = "A/A", 0
    elif "KIR2DL3" in kir_set: cen_h, cen_pts = "A/B", 1
    else: cen_h, cen_pts = "B/B", 2

    if not has_tel_b: tel_h, tel_pts = "A/A", 0
    elif "KIR3DL1" in kir_set: tel_h, tel_pts = "A/B", 1
    else: tel_h, tel_pts = "B/B", 2

    b_score = cen_pts + tel_pts
    return ("A/A" if b_score == 0 else "B/X"), b_score, cen_h, tel_h

def assess_donor_education(don_kir, don_ligs, don_b_alleles, don_c_alleles, codon_86_present):
    KIR_MAP = {"KIR2DL1":"C2", "KIR2DL2":"C1", "KIR2DL3":"C1", "KIR3DL1":"Bw4"}
    res = []
    has_b13_anergic = any(x in don_b_alleles for x in ['B*13:01', 'B*13:02'])

    for k, l in KIR_MAP.items():
        is_pos = k in don_kir
        has_l = (l=="Bw4" and any('Bw4' in x for x in don_ligs)) or (l in don_ligs)
        if k == "KIR3DL1" and is_pos and has_b13_anergic and not any('Bw4' in x for x in don_ligs):
            edu = "Non Educato (Anergico / L80)"
        else:
            edu = "Educato" if (is_pos and has_l) else "Non Educato" if is_pos else "N/A"
        res.append({'KIR':k, 'Lig': l, 'Status': "Pos" if is_pos else "Neg", 'Edu': edu})

    has_bw6 = "Bw6" in don_ligs
    res.append({'KIR': '— (No KIR)', 'Lig': 'Bw6', 'Status': "Presente" if has_bw6 else "Assente", 'Edu': "N/A (Neutro)"})

    k2ds1_pos = "KIR2DS1" in don_kir
    if not k2ds1_pos:
        res.append({'KIR': 'KIR2DS1', 'Lig': 'C1/C2', 'Status': "Neg", 'Edu': "N/A"})
    else:
        don_c_ligs = set()
        for a in don_c_alleles:
            if any(exc in a for exc in ['C*16:02', 'C*16:04', 'C*12:04']): don_c_ligs.add('C2')
            elif any(exc in a for exc in ['C*16:01', 'C*12:02', 'C*12:03']): don_c_ligs.add('C1')
            else:
                p = extract_allele_number(a)
                if p in C1_ALLELES_PREFIXES: don_c_ligs.add('C1')
                elif p in C2_ALLELES_PREFIXES: don_c_ligs.add('C2')

        k2ds1_edu = "Iporesponsive" if 'C2' in don_c_ligs else "Educato"
        k2ds1_label = ("C1/C2" if 'C1' in don_c_ligs else "C2/C2") if 'C2' in don_c_ligs else "C1/C1"
        res.append({'KIR': 'KIR2DS1', 'Lig': f"Locus C ({k2ds1_label})", 'Status': "Pos", 'Edu': k2ds1_edu})

    k3dl2_pos = "KIR3DL2" in don_kir
    has_a3_a11 = "A3/A11" in don_ligs
    k3dl2_edu = "N/A" if not k3dl2_pos else "Non Educato" if not has_a3_a11 else "Non Educato (Intrappolato 86+)" if codon_86_present else "Educato"
    res.append({'KIR': 'KIR3DL2', 'Lig': 'A3/A11', 'Status': "Pos" if k3dl2_pos else "Neg", 'Edu': k3dl2_edu})

    return res

def calculate_vectors_quantitative(edu_results, don_ligs, rec_ligs, haplo, b_score, k2ds1_edu):
    d_clean = {x.split('-')[0] for x in don_ligs}
    r_clean = {x.split('-')[0] for x in rec_ligs}
    gvh_mismatches = [r['KIR'] for r in edu_results if r['Edu'] in ["Educato", "Educato/Responsive"] and (r['Lig'].split('-')[0] if "Locus C" not in r['Lig'] else "C2") in d_clean and (r['Lig'].split('-')[0] if "Locus C" not in r['Lig'] else "C2") not in r_clean]

    if not gvh_mismatches: gvh_score, gvh_status, gvh_color = 1, "Nullo / Trascurabile", "#27ae60"
    elif len(gvh_mismatches) == 1: gvh_score, gvh_status, gvh_color = 2, f"Basso ({', '.join(gvh_mismatches)} ML)", "#f1c40f"
    elif len(gvh_mismatches) == 2 and b_score < 2: gvh_score, gvh_status, gvh_color = 3, f"Moderato ({', '.join(gvh_mismatches)} ML)", "#e67e22"
    else: gvh_score, gvh_status, gvh_color = 4, "Alto Rischio GVHD", "#e74c3c"

    hvg_mismatches = [lig for lig in r_clean if lig not in d_clean and lig != 'Bw6']
    if not hvg_mismatches: hvg_score, hvg_status, hvg_color = 1, "Nullo (Graft protetto)", "#27ae60"
    elif len(hvg_mismatches) == 1: hvg_score, hvg_status, hvg_color = 3, f"Moderato (Manca {', '.join(hvg_mismatches)})", "#e67e22"
    else: hvg_score, hvg_status, hvg_color = 5, "Massimo Rischio Rigetto Graft", "#e74c3c"

    gvl_points = 1.0
    reasons = []
    if gvh_mismatches: gvl_points += 1.5; reasons.append("Missing Ligand (+1.5)")
    if haplo == "B/X": gvl_points += 1.0; reasons.append("Aplotipo B (+1.0)")
    if "Educato" in k2ds1_edu: gvl_points += 1.5; reasons.append("KIR2DS1 Attivo (+1.5)")
    
    gvl_score = min(5, max(1, round(gvl_points)))
    gvl_labels = {1: "Minimo", 2: "Basso", 3: "Moderato", 4: "Alto Potenziale", 5: "Massimo Effetto GvL"}
    gvl_status = f"{gvl_labels[gvl_score]} [{', '.join(reasons) if reasons else 'Effetto NK basale'}]"
    gvl_color = "#27ae60" if gvl_score >= 4 else "#2980b9" if gvl_score == 3 else "#7f8c8d"

    return gvh_score, gvh_status, gvh_color, hvg_score, hvg_status, hvg_color, gvl_score, gvl_status, gvl_color
