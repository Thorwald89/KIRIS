import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import streamlit as st
import pandas as pd
from src.kir_analysis import (
    validate_and_extract_alleles, check_kir_biological_coherence,
    map_hla_to_ligands, calculate_b_content_logic, assess_donor_education,
    calculate_vectors_quantitative, CEN_GENES, TEL_GENES
)

st.set_page_config(page_title="KIRIS PRO - Alloreattività NK & Immunogenetica", page_icon="🧬", layout="wide")

st.title("🧬 KIRis PRO: Alloreactivity & Immunogenetics Toolkit")
st.markdown("**Versione Validata ISO 13485** con eccezione $B*13$ anergica (Leucina-80), B-Content e vettori di alloreattività.")

tab1, tab2 = st.tabs([
    "📊 Analisi Alloreattività & Vettori Clinici", 
    "📁 Analisi Coorte (CSV/Excel)"
])

# --- TAB 1: ANALISI COMPLETA KIRIS PRO ---
with tab1:
    col_rec, col_don = st.columns(2)
    with col_rec:
        st.subheader("1. Ricevente (HLA Data)")
        rec_a = st.text_input("Rec HLA-A", "A*03:01,A*11:01")
        rec_b = st.text_input("Rec HLA-B", "B*08:01,B*44:02")
        rec_c = st.text_input("Rec HLA-C", "C*07:01,C*05:01")

    with col_don:
        st.subheader("2. Donatore (HLA Data)")
        don_a = st.text_input("Don HLA-A", "A*02:01,A*03:01")
        don_b = st.text_input("Don HLA-B", "B*44:02,B*07:02")
        don_c = st.text_input("Don HLA-C", "C*03:01,C*08:01")

    st.markdown("---")
    st.subheader("3. Assetto Genico KIR Donatore (16 Loci)")
    
    col_cen, col_tel = st.columns(2)
    default_selected = {'KIR3DL3', 'KIR2DP1', 'KIR2DL1', 'KIR2DL3', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR2DS4', 'KIR3DL2'}
    
    with col_cen:
        st.markdown("**Regione Centromerica (CEN):**")
        selected_cen = [g for g in CEN_GENES if st.checkbox(g, value=(g in default_selected), key=f"cb_{g}")]
    with col_tel:
        st.markdown("**Regione Telomerica (TEL):**")
        selected_tel = [g for g in TEL_GENES if st.checkbox(g, value=(g in default_selected), key=f"cb_{g}")]

    selected_kir = set(selected_cen + selected_tel)
    codon_86 = st.checkbox("KIR3DL2: Presenza Codone 86 (SSP)", value=False)

    if st.button("🚀 Esegui Analisi Completa KIRis", type="primary"):
        try:
            r_a = validate_and_extract_alleles(rec_a, "Rec HLA-A")
            r_b = validate_and_extract_alleles(rec_b, "Rec HLA-B")
            r_c = validate_and_extract_alleles(rec_c, "Rec HLA-C")
            d_a = validate_and_extract_alleles(don_a, "Don HLA-A")
            d_b = validate_and_extract_alleles(don_b, "Don HLA-B")
            d_c = validate_and_extract_alleles(don_c, "Don HLA-C")
        except ValueError as err:
            st.error(f"⚠️ {str(err)}")
            st.stop()

        warnings = check_kir_biological_coherence(selected_kir)
        for w in warnings:
            st.warning(w)

        r_ligs = map_hla_to_ligands(r_c, r_b, r_a)
        d_ligs = map_hla_to_ligands(d_c, d_b, d_a)
        haplo, b_score, cen_h, tel_h = calculate_b_content_logic(selected_kir)
        edu_res = assess_donor_education(selected_kir, d_ligs, d_b, d_c, codon_86)

        k2ds1_edu_status = next((x['Edu'] for x in edu_res if x['KIR'] == 'KIR2DS1'), "N/A")
        gvh_sc, gvh_s, gvh_c, hvg_sc, hvg_s, hvg_c, gvl_sc, gvl_s, gvl_c = calculate_vectors_quantitative(
            edu_res, d_ligs, r_ligs, haplo, b_score, k2ds1_edu_status
        )

        st.markdown("---")
        
        # Mappa Cromosomica HTML
        all_ordered_genes = [
            ("3DL3", "KIR3DL3", "FW"), ("2DS2", "KIR2DS2", ""), ("2DL2", "KIR2DL2", ""),
            ("2DL3", "KIR2DL3", ""), ("2DS3", "KIR2DS3", ""), ("2DP1", "KIR2DP1", "Ψ"),
            ("2DL1", "KIR2DL1", ""), ("3DP1", "KIR3DP1", "Ψ"), ("2DL4", "KIR2DL4", "FW"),
            ("3DL1", "KIR3DL1", ""), ("3DS1", "KIR3DS1", ""), ("2DL5A", "KIR2DL5A", ""),
            ("2DL5B", "KIR2DL5B", ""), ("2DS5", "KIR2DS5", ""), ("2DS1", "KIR2DS1", ""),
            ("2DS4", "KIR2DS4", ""), ("3DL2", "KIR3DL2", "FW")
        ]

        gene_boxes_list = []
        for short_name, full_name, badge in all_ordered_genes:
            is_present = full_name in selected_kir
            bg_color = "#1e3a8a" if is_present else "#e2e8f0"
            text_color = "#ffffff" if is_present else "#94a3b8"
            border_style = "2px solid #2563eb" if is_present else "1px dashed #cbd5e1"
            badge_html = f"<span style='font-size:0.65rem; background:#3b82f6; color:white; padding:1px 4px; border-radius:3px; margin-left:3px;'>{badge}</span>" if badge else ""
            
            box = f'<div style="flex:1; min-width:52px; background:{bg_color}; color:{text_color}; border:{border_style}; border-radius:6px; padding:8px 4px; text-align:center; font-weight:bold; font-size:0.8rem; box-shadow:0 1px 3px rgba(0,0,0,0.1);">{short_name}{badge_html}</div>'
            gene_boxes_list.append(box)

        boxes_str = "".join(gene_boxes_list)

        st.markdown(
            f'<div style="background-color:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:15px; margin-bottom:25px;">'
            f'<div style="font-weight:bold; color:#1e293b; margin-bottom:10px; display:flex; justify-content:space-between;">'
            f'<span>🧬 MAPPA INTERATTIVA DEL LOCUS CROMOSOMICO KIR</span><span style="font-size:0.85rem; color:#64748b;">TELOMERO ➡️</span>'
            f'</div>'
            f'<div style="display:flex; gap:4px; overflow-x:auto; padding-bottom:5px;">{boxes_str}</div>'
            f'</div>',
            unsafe_allow_html=True
        )

        # Cards Report HTML
        d_ligs_str = ", ".join(d_ligs)
        r_ligs_str = ", ".join(r_ligs)
        
        st.markdown(
            f'<div style="background:linear-gradient(135deg, #1e293b 0%, #0f172a 100%); color:white; border-radius:12px; padding:20px; margin-bottom:25px;">'
            f'<h3 style="margin-top:0; color:#38bdf8;">📊 KIRis: Report Immunogenetico Alloreattività</h3>'
            f'<div style="display:grid; grid-template-columns:1fr 1fr; gap:20px; margin-top:15px;">'
            f'<div style="background:rgba(255,255,255,0.05); border:1px solid rgba(255,255,255,0.1); border-radius:8px; padding:15px;">'
            f'<h4 style="margin-top:0; color:#f1f5f9;">1. Profilo KIR & Struttura Aplotipica (Donatore)</h4>'
            f'<p style="margin:5px 0;"><b>Aplotipo Matematico:</b> <span style="color:#38bdf8;">{haplo}</span> | <b>B-Score:</b> {b_score}</p>'
            f'<p style="margin:5px 0;"><b>Regione Centromerica (CEN):</b> {cen_h} | <b>Regione Telomerica (TEL):</b> {tel_h}</p>'
            f'</div>'
            f'<div style="background:rgba(255,255,255,0.05); border:1px solid rgba(255,255,255,0.1); border-radius:8px; padding:15px;">'
            f'<h4 style="margin-top:0; color:#f1f5f9;">2. Mappatura Reale Epitopi ed Espressione</h4>'
            f'<p style="margin:5px 0;"><b>Epitopi Donatore:</b> {d_ligs_str}</p>'
            f'<p style="margin:5px 0;"><b>Epitopi Ricevente:</b> {r_ligs_str}</p>'
            f'</div>'
            f'</div>'
            f'</div>',
            unsafe_allow_html=True
        )

        # Tabella Educazione NK
        st.markdown("#### 🎓 Modello di Educazione Funzionale NK")
        st.dataframe(pd.DataFrame(edu_res), use_container_width=True)

        # Scoring Vettori Clinici
        st.markdown("#### 🎯 Scored Vettori Clinici Quantitativi (1-5)")
        st.write(f"⚡ **Vettore GvH:** Score {gvh_sc}/5 — {gvh_s}")
        st.progress(gvh_sc * 20)
        
        st.write(f"🛡️ **Vettore HvG:** Score {hvg_sc}/5 — {hvg_s}")
        st.progress(hvg_sc * 20)

        st.write(f"🎯 **Vettore GvL:** Score {gvl_sc}/5 — {gvl_s}")
        st.progress(gvl_sc * 20)

        # --- LEGENDA LOGICA VETTORI CLINICI ---
        with st.expander("📖 Legenda e Logica di Calcolo dei Vettori Clinici"):
            st.markdown("""
            * **⚡ Vettore GvH (Graft-vs-Host):**
              * **1/5 (Nullo):** Il ricevente possiede tutti i ligandi inibitori prescritti dal donatore (assenza di mismatch inibitorio).
              * **2-5 (Incrementale):** Presenza di *Missing-Self* sul ricevente coordinata da KIR inibitori del donatore **educati**, oltre alla potenziale attivazione diretta mediata da aplotipi B/X o KIR2DS1.

            * **🛡️ Vettore HvG (Host-vs-Graft / Rischio Rigetto):**
              * **1/5 (Nullo):** Il donatore possiede tutti i ligandi espressi dal ricevente.
              * **2-5 (Incrementale):** Il donatore sprovvisto di un ligando presente nel ricevente (es. *Manca C2* o *Manca Bw4*) espone l'innesto all'attacco da parte di cellule NK residue del ricevente.

            * **🎯 Vettore GvL (Graft-vs-Leukemia / Effetto Anti-Leucemico):**
              * **Base (Missing-Self):** Assenza nel ricevente di ligandi inibitori riconosciuti da KIR educati nel donatore (+1.5 - +2.0).
              * **Boost Aplotipo B (+1.0):** Presenza di genotipo B/X nel donatore (elevato contenuto di geni attivatori).
              * **Boost KIR2DS1 (+1.5):** Donatore con KIR2DS1 **educato** (presenza C1/C1 o C1/C2 nel donatore) e ricevente **C2+** (massima risposta alloreattiva tumorale).
            """)

        def get_excel_label(kir_name):
            match_kir = next((x for x in edu_res if x['KIR'] == kir_name), None)
            if not match_kir or match_kir['Status'] == 'Neg': return "NEG"
            if kir_name == 'KIR3DL2': return "POS +L" if match_kir['Edu'] == 'Educato' else "POS"
            return "POS, E" if ('Educato' in match_kir['Edu'] or match_kir['Edu'] == 'Educato/Responsive') else "POS, NE"

        excel_row_string = f"{haplo}\t{b_score}\t{get_excel_label('KIR3DL1')}\t{get_excel_label('KIR2DL1')}\t{get_excel_label('KIR2DL2')}\t{get_excel_label('KIR2DL3')}\t{get_excel_label('KIR2DS1')}\t{get_excel_label('KIR3DL2')}\t{gvh_sc}\t{hvg_sc}\t{gvl_sc}"
        
        st.markdown("#### 📋 Stringa per Excel")
        st.code(excel_row_string, language="text")

# --- TAB 2: UPLOAD DATASET ---
with tab2:
    st.subheader("Analisi Coorte Batch")
    f = st.file_uploader("Carica file dataset (.csv / .xlsx)", type=["csv", "xlsx"])
    if f is not None:
        df = pd.read_csv(f) if f.name.endswith(".csv") else pd.read_excel(f)
        st.dataframe(df.head())
