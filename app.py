import streamlit as st
import pandas as pd
from src.kir_analysis import (
    validate_and_extract_alleles, check_kir_biological_coherence,
    map_hla_to_ligands, calculate_b_content_logic, assess_donor_education,
    calculate_vectors_quantitative, CEN_GENES, TEL_GENES
)
from src.hed_engine import calculate_hed

st.set_page_config(page_title="KIRIS PRO - Alloreattività NK & Immunogenetica", page_icon="🧬", layout="wide")

st.title("🧬 KIRis PRO: Alloreactivity & Immunogenetics Toolkit")
st.markdown("**Versione Validata ISO 13485** con eccezione $B*13$ anergica (Leucina-80), B-Content e Divergenza HED.")

tab1, tab2, tab3 = st.tabs([
    "📊 Analisi Alloreattività & Vettori Clinici", 
    "🧬 Divergenza HED (Locus A/B/C)", 
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

        st.markdown("### 📊 Report Immunogenetico")
        m1, m2 = st.columns(2)
        with m1:
            st.info(f"**Aplotipo:** {haplo} | **B-Score:** {b_score}\n\n**CEN:** {cen_h} | **TEL:** {tel_h}")
            st.write(f"**Epitopi Donatore:** {', '.join(d_ligs)}")
        with m2:
            st.write(f"**Epitopi Ricevente:** {', '.join(r_ligs)}")

        st.markdown("#### 🎓 Modello di Educazione Funzionale NK")
        st.dataframe(pd.DataFrame(edu_res), use_container_width=True)

        st.markdown("#### 🎯 Scored Vettori Clinici Quantitativi (1-5)")
        st.write(f"⚡ **Vettore GvH:** Score {gvh_sc}/5 — {gvh_s}")
        st.progress(gvh_sc * 20)
        
        st.write(f"🛡️ **Vettore HvG:** Score {hvg_sc}/5 — {hvg_s}")
        st.progress(hvg_sc * 20)

        st.write(f"🎯 **Vettore GvL:** Score {gvl_sc}/5 — {gvl_s}")
        st.progress(gvl_sc * 20)

        # Matrice Excel
        def get_excel_label(kir_name):
            match_kir = next((x for x in edu_res if x['KIR'] == kir_name), None)
            if not match_kir or match_kir['Status'] == 'Neg': return "NEG"
            if kir_name == 'KIR3DL2': return "POS +L" if match_kir['Edu'] == 'Educato' else "POS"
            return "POS, E" if ('Educato' in match_kir['Edu'] or match_kir['Edu'] == 'Educato/Responsive') else "POS, NE"

        excel_row_string = f"{haplo}\t{b_score}\t{get_excel_label('KIR3DL1')}\t{get_excel_label('KIR2DL1')}\t{get_excel_label('KIR2DL2')}\t{get_excel_label('KIR2DL3')}\t{get_excel_label('KIR2DS1')}\t{get_excel_label('KIR3DL2')}\t{gvh_sc}\t{hvg_sc}\t{gvl_sc}"
        
        st.markdown("#### 📋 Stringa per Excel")
        st.code(excel_row_string, language="text")

# --- TAB 2: HED ENGINE ---
with tab2:
    st.subheader("Calcolo HLA Evolutionary Distance (HED)")
    st.markdown("Analisi della divergenza amminoacidica nel **Locus A** per valutare l'effetto sulla presentazione peptidica e sul rischio di anti-viral failure.")
    c1, c2 = st.columns(2)
    with c1:
        seq_a = st.text_area("Sequenza Allele 1 (Locus A)", "SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDRNTRNVKAQSQTDRVDLGTLRGYYNQSEAGS")
    with c2:
        seq_b = st.text_area("Sequenza Allele 2 (Locus A)", "SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWEEETRNVKAQSQTDRVDLGTLRGYYNQSEAGS")

    if st.button("Calcola HED"):
        try:
            val = calculate_hed(seq_a, seq_b)
            st.metric("Divergenza HED (Distanza Grantham)", f"{val}")
            if val < 5.0:
                st.warning("⚠️ **Bassa Divergenza (LHH/LLL):** Rischio compromissione presentazione peptidica virale.")
            else:
                st.success("✅ **Alta Divergenza:** Ampio repertorio peptidico.")
        except Exception as e:
            st.error(f"Errore: {e}")

# --- TAB 3: UPLOAD DATASET ---
with tab3:
    st.subheader("Analisi Coorte Batch")
    f = st.file_uploader("Carica file dataset (.csv / .xlsx)", type=["csv", "xlsx"])
    if f is not None:
        df = pd.read_csv(f) if f.name.endswith(".csv") else pd.read_excel(f)
        st.dataframe(df.head())
