import sys
import os
import io

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import streamlit as st
import pandas as pd
import altair as alt
from src.kir_analysis import (
    validate_and_extract_alleles, check_kir_biological_coherence,
    map_hla_to_ligands, calculate_b_content_logic, assess_donor_education,
    calculate_vectors_quantitative, calculate_donor_recommendation,
    generate_radar_chart_svg, analyze_single_pair, format_excel_row,
    CEN_GENES, TEL_GENES, FRAMEWORK_GENES
)

st.set_page_config(page_title="KIRIS PRO - Alloreattività NK & Immunogenetica", page_icon="🧬", layout="wide")

# Stile CSS globale personalizzato
st.markdown("""
<style>
    .metric-card {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 14px 18px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        margin-bottom: 12px;
    }
    .badge-pill {
        display: inline-block;
        font-size: 0.75rem;
        font-weight: 700;
        padding: 2px 8px;
        border-radius: 12px;
        margin-left: 6px;
    }
</style>
""", unsafe_allow_html=True)

st.title("🧬 KIRis PRO: Alloreactivity & Immunogenetics Toolkit")
st.markdown("**Versione Conforme a Validazione Deterministica ISO 13485**, Release 0.3.")

# ==============================================================================
# PRESET CLINICI RAPIDI
# ==============================================================================
PRESET_DICT = {
    "Personalizzato (Inserimento manuale)": None,
    "⭐ Candidato Ottimale: Donatore C1/C1 KIR2DS1+ vs Ricevente C2+ (Missing C2 + Attivazione KIR2DS1)": {
        "rec_a": "A*02:01,A*03:01", "rec_b": "B*07:02,B*08:01", "rec_c": "C*02:02,C*04:01",
        "don_a": "A*03:01,A*11:01", "don_b": "B*44:02,B*57:01", "don_c": "C*07:01,C*08:01",
        "kir": {'KIR3DL3', 'KIR2DS2', 'KIR2DL2', 'KIR2DL3', 'KIR2DP1', 'KIR2DL1', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR2DS1', 'KIR2DS4', 'KIR3DL2'},
        "kir_null": False
    },
    "🛡️ Trapianto Standard: Aplotipo A/A (Missing-Self Assente, Graft Protetto)": {
        "rec_a": "A*01:01,A*02:01", "rec_b": "B*07:02,B*44:02", "rec_c": "C*07:01,C*05:01",
        "don_a": "A*01:01,A*03:01", "don_b": "B*07:02,B*44:02", "don_c": "C*07:01,C*05:01",
        "kir": {'KIR3DL3', 'KIR2DL3', 'KIR2DP1', 'KIR2DL1', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR2DS4', 'KIR3DL2'},
        "kir_null": False
    },
    "🔥 High-Affinity Bw4-I80 Missing-Self (Ricevente Bw6/Bw6)": {
        "rec_a": "A*02:01,A*03:01", "rec_b": "B*07:02,B*08:01", "rec_c": "C*07:01,C*03:01",
        "don_a": "A*02:01,A*24:02", "don_b": "B*51:01,B*08:01", "don_c": "C*07:01,C*03:01",
        "kir": {'KIR3DL3', 'KIR2DL3', 'KIR2DP1', 'KIR2DL1', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR2DS4', 'KIR3DL2'},
        "kir_null": False
    },
    "⚠️ Rischio Rigetto Innesto (HvG Alto): Donatore C1/C1 vs Ricevente C2/C2 Bw4": {
        "rec_a": "A*24:02,A*02:01", "rec_b": "B*57:01,B*44:02", "rec_c": "C*02:02,C*04:01",
        "don_a": "A*03:01,A*11:01", "don_b": "B*07:02,B*08:01", "don_c": "C*07:01,C*08:01",
        "kir": {'KIR3DL3', 'KIR2DL3', 'KIR2DP1', 'KIR2DL1', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR2DS4', 'KIR3DL2'},
        "kir_null": False
    }
}

APLOTIPO_A_GENES = {'KIR3DL3', 'KIR2DL3', 'KIR2DP1', 'KIR2DL1', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR2DS4', 'KIR3DL2'}
APLOTIPO_B_FULL = {'KIR3DL3', 'KIR2DS2', 'KIR2DL2', 'KIR2DL3', 'KIR2DP1', 'KIR2DL1', 'KIR3DP1', 'KIR2DL4', 'KIR3DL1', 'KIR3DS1', 'KIR2DL5A', 'KIR2DL5B', 'KIR2DS5', 'KIR2DS1', 'KIR3DL2'}

tab1, tab2 = st.tabs([
    "📊 Analisi Singola Coppia & Vettori Clinici", 
    "📁 Analisi Coorte Batch (CSV/Excel)"
])

# ==============================================================================
# --- TAB 1: ANALISI COMPLETA SINGOLA COPPIA RICEVENTE / DONATORE ---
# ==============================================================================
with tab1:
    # Selettore Preset Rapido
    selected_preset_name = st.selectbox(
        "⚡ Caricamento Rapido Preset Clinico (Quick Load):",
        list(PRESET_DICT.keys()),
        index=0,
        help="Seleziona un caso clinico tipo per precompilare istantaneamente tutti i dati."
    )

    preset_data = PRESET_DICT[selected_preset_name]

    col_rec, col_don = st.columns(2)
    with col_rec:
        st.subheader("1. Ricevente (HLA Data)")
        val_ra = preset_data["rec_a"] if preset_data else "A*03:01,A*11:01"
        val_rb = preset_data["rec_b"] if preset_data else "B*08:01,B*44:02"
        val_rc = preset_data["rec_c"] if preset_data else "C*07:01,C*05:01"

        rec_a = st.text_input("Rec HLA-A", val_ra, key="input_rec_a", help="Separa gli alleli con virgola (es. A*03:01, A*11:01)")
        rec_b = st.text_input("Rec HLA-B", val_rb, key="input_rec_b")
        rec_c = st.text_input("Rec HLA-C", val_rc, key="input_rec_c")

    with col_don:
        st.subheader("2. Donatore (HLA Data)")
        val_da = preset_data["don_a"] if preset_data else "A*02:01,A*03:01"
        val_db = preset_data["don_b"] if preset_data else "B*44:02,B*07:02"
        val_dc = preset_data["don_c"] if preset_data else "C*03:01,C*08:01"

        don_a = st.text_input("Don HLA-A", val_da, key="input_don_a")
        don_b = st.text_input("Don HLA-B", val_db, key="input_don_b")
        don_c = st.text_input("Don HLA-C", val_dc, key="input_don_c")

    st.markdown("---")
    st.subheader("3. Assetto Genico KIR Donatore (16 Loci) & Controllo Qualità")

    # Pulsanti di selezione rapida aplotipo KIR
    st.markdown("**Scorciatoie Assetto KIR Donatore:**")
    btn_c1, btn_c2, btn_c3, btn_c4 = st.columns(4)

    if "custom_kir_selection" not in st.session_state:
        st.session_state.custom_kir_selection = set(preset_data["kir"]) if preset_data else set(APLOTIPO_A_GENES)

    if preset_data and st.session_state.get("last_preset") != selected_preset_name:
        st.session_state.custom_kir_selection = set(preset_data["kir"])
        st.session_state.last_preset = selected_preset_name

    with btn_c1:
        if st.button("🔵 Aplotipo A Standard (9 Geni)", use_container_width=True):
            st.session_state.custom_kir_selection = set(APLOTIPO_A_GENES)
            st.rerun()
    with btn_c2:
        if st.button("🟢 Aplotipo B Arricchito", use_container_width=True):
            st.session_state.custom_kir_selection = set(APLOTIPO_B_FULL)
            st.rerun()
    with btn_c3:
        if st.button("🧬 Seleziona Tutti i Loci", use_container_width=True):
            st.session_state.custom_kir_selection = set(CEN_GENES + TEL_GENES)
            st.rerun()
    with btn_c4:
        if st.button("⚪ Solo Geni Costitutivi (Framework)", use_container_width=True):
            st.session_state.custom_kir_selection = set(FRAMEWORK_GENES)
            st.rerun()

    current_selected = st.session_state.custom_kir_selection

    col_cen, col_tel = st.columns(2)
    with col_cen:
        st.markdown("**Regione Centromerica (CEN):**")
        selected_cen = []
        for g in CEN_GENES:
            is_checked = st.checkbox(g, value=(g in current_selected), key=f"cb_cen_{g}")
            if is_checked: selected_cen.append(g)

    with col_tel:
        st.markdown("**Regione Telomerica (TEL):**")
        selected_tel = []
        for g in TEL_GENES:
            is_checked = st.checkbox(g, value=(g in current_selected), key=f"cb_tel_{g}")
            if is_checked: selected_tel.append(g)

    selected_kir = set(selected_cen + selected_tel)
    
    # Variante Null KIR3DL1
    val_null = preset_data["kir_null"] if preset_data else False
    kir3dl1_null = st.checkbox(
        "KIR3DL1: Variante Null / ER-trapped (*004)", 
        value=val_null, 
        help="Se attivo, esclude KIR3DL1 dal licensing di superficie (segnale Null)."
    )

    if st.button("🚀 Esegui Analisi Completa KIRis", type="primary", use_container_width=True):
        try:
            results = analyze_single_pair(
                rec_a, rec_b, rec_c,
                don_a, don_b, don_c,
                selected_kir,
                kir3dl1_null_variant=kir3dl1_null
            )
        except ValueError as err:
            st.error(f"⚠️ {str(err)}")
            st.stop()

        for w in results["kir_warnings"]:
            st.warning(w)

        d_ligs = results["d_ligs"]
        r_ligs = results["r_ligs"]
        haplo = results["haplo"]
        b_score = results["b_score"]
        cen_h = results["cen_h"]
        tel_h = results["tel_h"]
        edu_res = results["edu_res"]
        gvh_sc = results["gvh"]["score"]
        gvh_s = results["gvh"]["status"]
        gvh_c = results["gvh"]["color"]
        hvg_sc = results["hvg"]["score"]
        hvg_s = results["hvg"]["status"]
        hvg_c = results["hvg"]["color"]
        gvl_sc = results["gvl"]["score"]
        gvl_s = results["gvl"]["status"]
        gvl_c = results["gvl"]["color"]
        recom = results["recommendation"]
        radar_svg = results["radar_svg"]
        excel_row_string = results["excel_str"]

        st.markdown("---")

        # ======================================================================
        # 1. RACCOMANDAZIONE CLINICA SINTETICA DEL DONATORE (DONOR RECOMMENDATION)
        # ======================================================================
        st.markdown(
            f'<div style="background-color:{recom["badge_bg"]}; border:2px solid {recom["color"]}; '
            f'border-radius:12px; padding:18px; margin-bottom:25px;">'
            f'<div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">'
            f'<h3 style="margin:0; color:{recom["color"]}; font-size:1.25rem;">{recom["title"]}</h3>'
            f'<span style="background:{recom["color"]}; color:white; font-weight:bold; font-size:0.8rem; padding:4px 10px; border-radius:20px;">'
            f'{recom["level"]}</span>'
            f'</div>'
            f'<p style="margin:0; color:#1e293b; font-size:0.95rem; line-height:1.5;">{recom["recommendation"]}</p>'
            f'</div>',
            unsafe_allow_html=True
        )

        # ======================================================================
        # 2. MAPPA CROMOSOMICA INTERATTIVA A COLORI FUNZIONALI
        # ======================================================================
        # Classificazione funzionale: INIBITORI, ATTIVATORI, FRAMEWORK, PSEUDOGENI
        gene_annotations = {
            "KIR3DL3": ("3DL3", "FW", "#7c3aed", "Framework Costitutivo"),
            "KIR2DS2": ("2DS2", "Act", "#059669", "Attivatorio"),
            "KIR2DL2": ("2DL2", "Inh", "#2563eb", "Inibitorio"),
            "KIR2DL3": ("2DL3", "Inh", "#2563eb", "Inibitorio"),
            "KIR2DS3": ("2DS3", "Act", "#059669", "Attivatorio"),
            "KIR2DP1": ("2DP1", "Ψ", "#64748b", "Pseudogene"),
            "KIR2DL1": ("2DL1", "Inh", "#2563eb", "Inibitorio"),
            "KIR3DP1": ("3DP1", "Ψ", "#64748b", "Pseudogene"),
            "KIR2DL4": ("2DL4", "FW", "#7c3aed", "Framework Costitutivo"),
            "KIR3DL1": ("3DL1", "Null" if kir3dl1_null else "Inh", "#ef4444" if kir3dl1_null else "#2563eb", "Inibitorio"),
            "KIR3DS1": ("3DS1", "Act", "#059669", "Attivatorio"),
            "KIR2DL5A": ("2DL5A", "Inh", "#2563eb", "Inibitorio"),
            "KIR2DL5B": ("2DL5B", "Inh", "#2563eb", "Inibitorio"),
            "KIR2DS5": ("2DS5", "Act", "#059669", "Attivatorio"),
            "KIR2DS1": ("2DS1", "Act", "#059669", "Attivatorio"),
            "KIR2DS4": ("2DS4", "Act", "#059669", "Attivatorio"),
            "KIR3DL2": ("3DL2", "FW", "#7c3aed", "Framework Costitutivo"),
        }

        cen_order = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DS3", "KIR2DP1", "KIR2DL1", "KIR3DP1"]
        tel_order = ["KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DL5B", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"]

        def render_boxes(gene_keys):
            boxes = []
            for gk in gene_keys:
                short, badge, color_func, _ = gene_annotations[gk]
                is_present = gk in selected_kir
                bg_color = color_func if is_present else "#f1f5f9"
                txt_color = "#ffffff" if is_present else "#94a3b8"
                border_st = f"2px solid {color_func}" if is_present else "1px dashed #cbd5e1"
                badge_html = f"<span style='font-size:0.6rem; background:rgba(0,0,0,0.25); color:white; padding:1px 3px; border-radius:3px; margin-left:2px;'>{badge}</span>" if badge else ""
                boxes.append(
                    f'<div style="flex:1; min-width:48px; background:{bg_color}; color:{txt_color}; '
                    f'border:{border_st}; border-radius:6px; padding:7px 2px; text-align:center; '
                    f'font-weight:bold; font-size:0.75rem; box-shadow:0 1px 2px rgba(0,0,0,0.06);">'
                    f'{short}{badge_html}</div>'
                )
            return "".join(boxes)

        cen_boxes_str = render_boxes(cen_order)
        tel_boxes_str = render_boxes(tel_order)

        st.markdown(
            f'<div style="background-color:#ffffff; border:1px solid #e2e8f0; border-radius:12px; padding:18px; margin-bottom:25px; box-shadow:0 2px 4px rgba(0,0,0,0.02);">'
            f'<div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">'
            f'<span style="font-weight:bold; font-size:0.95rem; color:#1e293b;">🧬 MAPPA CROMOSOMICA COMPARATIVA DEL LOCUS KIR</span>'
            f'<span style="font-size:0.8rem; color:#64748b;">5\' CENTROMERO ➡️ TELOMERO 3\'</span>'
            f'</div>'
            f'<div style="display:flex; gap:6px; align-items:center; overflow-x:auto; padding-bottom:6px;">'
            f'{cen_boxes_str}'
            f'<div style="background:#f8fafc; border:2px solid #cbd5e1; border-radius:4px; padding:6px 4px; writing-mode:vertical-rl; text-orientation:mixed; font-size:0.6rem; font-weight:bold; color:#64748b; text-align:center;">HOTSPOT</div>'
            f'{tel_boxes_str}'
            f'</div>'
            f'<div style="display:flex; gap:16px; margin-top:10px; font-size:0.75rem; color:#475569; justify-content:center;">'
            f'<span><span style="display:inline-block; width:10px; height:10px; background:#2563eb; border-radius:2px; margin-right:4px;"></span>Inibitori</span>'
            f'<span><span style="display:inline-block; width:10px; height:10px; background:#059669; border-radius:2px; margin-right:4px;"></span>Attivatori</span>'
            f'<span><span style="display:inline-block; width:10px; height:10px; background:#7c3aed; border-radius:2px; margin-right:4px;"></span>Framework</span>'
            f'<span><span style="display:inline-block; width:10px; height:10px; background:#64748b; border-radius:2px; margin-right:4px;"></span>Pseudogeni (Ψ)</span>'
            f'</div>'
            f'</div>',
            unsafe_allow_html=True
        )

        # ======================================================================
        # 3. CARDS REPORT PROFILO & EPITOPI
        # ======================================================================
        d_ligs_str = ", ".join(sorted(d_ligs))
        r_ligs_str = ", ".join(sorted(r_ligs))
        
        st.markdown(
            f'<div style="background:linear-gradient(135deg, #1e293b 0%, #0f172a 100%); color:white; border-radius:12px; padding:20px; margin-bottom:25px;">'
            f'<h3 style="margin-top:0; color:#38bdf8; font-size:1.15rem;">📊 KIRis: Profilo Immunogenetico della Coppia</h3>'
            f'<div style="display:grid; grid-template-columns:1fr 1fr; gap:20px; margin-top:15px;">'
            f'<div style="background:rgba(255,255,255,0.05); border:1px solid rgba(255,255,255,0.1); border-radius:8px; padding:15px;">'
            f'<h4 style="margin-top:0; color:#f1f5f9; font-size:0.95rem;">1. Profilo KIR & Struttura Aplotipica (Donatore)</h4>'
            f'<p style="margin:5px 0;"><b>Aplotipo:</b> <span style="color:#38bdf8;">{haplo}</span> | <b>Cooley B-Score:</b> {b_score}/4</p>'
            f'<p style="margin:5px 0;"><b>Regione Centromerica (CEN):</b> {cen_h} | <b>Regione Telomerica (TEL):</b> {tel_h}</p>'
            f'</div>'
            f'<div style="background:rgba(255,255,255,0.05); border:1px solid rgba(255,255,255,0.1); border-radius:8px; padding:15px;">'
            f'<h4 style="margin-top:0; color:#f1f5f9; font-size:0.95rem;">2. Mappatura Ligandi & Bersagli Espressi</h4>'
            f'<p style="margin:5px 0;"><b>Epitopi Donatore:</b> {d_ligs_str}</p>'
            f'<p style="margin:5px 0;"><b>Epitopi Ricevente:</b> {r_ligs_str}</p>'
            f'</div>'
            f'</div>'
            f'</div>',
            unsafe_allow_html=True
        )

        # ======================================================================
        # 4. TABELLA EDUCAZIONE NK (RHEOSTAT MODEL)
        # ======================================================================
        st.markdown("#### 🎓 Modello di Educazione Funzionale NK (Stratificazione Graduata)")
        st.dataframe(pd.DataFrame(edu_res), use_container_width=True)

        # ======================================================================
        # 5. RADAR CHART SVG & CARD VETTORI CLINICI
        # ======================================================================
        st.markdown("#### 🎯 Valutazione Geometrica dei Vettori Clinici (1-5)")
        col_radar, col_cards = st.columns([5, 7])

        with col_radar:
            st.markdown(radar_svg, unsafe_allow_html=True)

        with col_cards:
            # Card GvH
            st.markdown(
                f'<div class="metric-card" style="border-left:5px solid {gvh_c};">'
                f'<div style="display:flex; justify-content:space-between;">'
                f'<b>⚡ Vettore GvH (Alloreattività Donatore ➡️ Ospite)</b>'
                f'<span style="font-weight:bold; color:{gvh_c};">Score {gvh_sc}/5</span>'
                f'</div>'
                f'<div style="font-size:0.85rem; color:#475569; margin-top:4px;">{gvh_s}</div>'
                f'</div>',
                unsafe_allow_html=True
            )
            st.progress(gvh_sc * 20)

            # Card HvG
            st.markdown(
                f'<div class="metric-card" style="border-left:5px solid {hvg_c}; margin-top:12px;">'
                f'<div style="display:flex; justify-content:space-between;">'
                f'<b>🛡️ Vettore HvG (Rischio Rigetto Graft da NK Ospite)</b>'
                f'<span style="font-weight:bold; color:{hvg_c};">Score {hvg_sc}/5</span>'
                f'</div>'
                f'<div style="font-size:0.85rem; color:#475569; margin-top:4px;">{hvg_s}</div>'
                f'</div>',
                unsafe_allow_html=True
            )
            st.progress(hvg_sc * 20)

            # Card GvL
            st.markdown(
                f'<div class="metric-card" style="border-left:5px solid {gvl_c}; margin-top:12px;">'
                f'<div style="display:flex; justify-content:space-between;">'
                f'<b>🎯 Vettore GvL (Efficacia Anti-Leucemica / Anti-Recidiva)</b>'
                f'<span style="font-weight:bold; color:{gvl_c};">Score {gvl_sc}/5</span>'
                f'</div>'
                f'<div style="font-size:0.85rem; color:#475569; margin-top:4px;">{gvl_s}</div>'
                f'</div>',
                unsafe_allow_html=True
            )
            st.progress(gvl_sc * 20)

        # ======================================================================
        # 6. LEGENDA SCIENTIFICA
        # ======================================================================
        with st.expander("📖 Legenda, Stratificazione Funzionale e Logica dei Vettori"):
            st.markdown(r"""
            ### 1. Stratificazione Graduata dell'Educazione NK (*Rheostat Model*)
            L'educazione NK segue il **Rheostat Model** (Joncker et al.), in cui l'affinità chimico-fisica del legame KIR-HLA determina il potenziale citotossico dell'NK in caso di *missing-self*.

            * **🟢 Licensing Forte (High-Affinity):**
              * **Bw4-Ile80 + KIR3DL1:** Presenza di Isoleucina al residuo 80 dell'HLA-B/A. Assicura massima costante di affinità ed elevata degranulazione.
              * **C2 + KIR2DL1:** Interazione ad alta affinità (Lisina-80) che genera la risposta inibitoria/licensing più solida nel sistema HLA-C.
            * **🟡 Licensing Moderato / Intermedio (Intermediate-Affinity):**
              * **Bw4-Thr80 + KIR3DL1:** Presenza di Treonina al residuo 80; legame a stabilità intermedia, educazione parziale.
              * **C1 + KIR2DL2:** Interazione a media affinità per il gruppo C1.
            * **🟠 Licensing Debole (Low-Affinity):**
              * **C1 + KIR2DL3:** Interazione molecolare a più bassa affinità per HLA-C1 rispetto a KIR2DL2.
            * **🔴 Anergico / No Licensing (Null / ER-trapped):**
              * **KIR3DL1\*004 (Null):** Difetto di ripiegamento che intrappola la proteina nel reticolo endoplasmatico; recettore non espresso in superficie.
              * **Bw4-Leu80 ($B*13$):** La sostituzione con Leucina al residuo 80 impedisce il legame inibitorio con KIR3DL1.
              * **Bw6 (No KIR):** L'epitopo Bw6 non possiede alcun recettore KIR inibitorio o attivatorio corrispondente.

            ---

            ### 2. Recettore Attivatorio KIR2DS1 (Modello Molecolare di Educazione & Bersaglio)
            * **Background Donatore $C1/C1$:** Le cellule NK $KIR2DS1^+$ non subiscono desensibilizzazione da parte del ligando cognito $C2$. Sono **non tollerizzate / ipereccitabili (*hyperexcitable*)**.
            * **Background Donatore $C1/C2$:** In presenza di una sola dose genica di $C2$, mantengono una responsività funzionale significativa: **educate/responsive**.
            * **Background Donatore $C2/C2$:** L'esposizione al doppio dosaggio di $C2$ induce un tuning down cronico: le cellule $KIR2DS1^+$ sono **anergiche / tollerizzate**.
            * **Condizione Bersaglio Ricevente per Effetto GvL:** Affinché $KIR2DS1$ (se ipereccitabile o responsivo) esplichi l'effetto anti-leucemico $GvL$, il **Ricevente DEVE possedere il ligando $C2$**. Se il ricevente è $C1/C1$, il target tumorale manca del ligando attivatore e non vi è induzione di lisi mediata da $KIR2DS1$.

            ---

            ### 3. Logica di Ponderazione dei Vettori Clinici (Score 1-5)
            * **⚡ Vettore GvH (Graft-vs-Host Alloreactivity):**
              * Valuta il mismatch inibitorio (*missing-self*) ponderato con le affinità del Rheostat Model ($Bw4\text{-}I80$ e $C2$ peso alto, $Bw4\text{-}T80$ e $C1$ con 2DL2 peso medio, $C1$ con 2DL3 peso debole) unito all'attivazione additiva da $KIR2DS1$ su bersaglio $C2^+$.
            * **🛡️ Vettore HvG (Host-vs-Graft / Rischio Rigetto):**
              * Valuta la mancata espressione nel donatore dei ligandi $C2$, $Bw4$ o $C1$ presenti nel ricevente.
            * **🎯 Vettore GvL (Graft-versus-Leukemia / Effetto Anti-Leucemico):**
              * Combinazione cooperativa tra Missing-Self inibitorio ad alta/media affinità, stimolazione attivatoria di $KIR2DS1$ (se ricevente è $C2^+$) e $B$-Score di Cooley ($B\text{-}Score \ge 2$ e arricchimento centromerico/telomerico).
            """)

        st.markdown("#### 📋 Stringa per Fogli Elettronici (Excel/Sheets)")
        st.code(excel_row_string, language="text")

# ==============================================================================
# --- TAB 2: ANALISI COORTE BATCH (CSV / EXCEL) ---
# ==============================================================================
with tab2:
    st.subheader("📁 Analisi Coorte Batch ad Alto Rendimento")
    st.markdown(
        "Carica un file **CSV** o **Excel (.xlsx)** contenente una coorte di coppie Donatore-Ricevente. "
        "Il motore eseguirà in modalità deterministica l'analisi immunogenetica per ciascun caso secondo i criteri ISO 13485."
    )

    # Template CSV scaricabile
    sample_template = (
        "ID,Rec_HLA_A,Rec_HLA_B,Rec_HLA_C,Don_HLA_A,Don_HLA_B,Don_HLA_C,Don_KIR,KIR3DL1_Null\n"
        "PAZ-001,A*03:01;A*11:01,B*08:01;B*44:02,C*07:01;C*05:01,A*02:01;A*03:01,B*44:02;B*07:02,C*03:01;C*08:01,\"KIR3DL3,KIR2DS2,KIR2DL2,KIR2DL3,KIR2DP1,KIR2DL1,KIR3DP1,KIR2DL4,KIR3DL1,KIR2DS1,KIR2DS4,KIR3DL2\",False\n"
        "PAZ-002,A*02:01;A*24:02,B*51:01;B*07:02,C*02:02;C*15:02,A*01:01;A*03:01,B*08:01;B*18:01,C*07:01;C*07:02,\"KIR3DL3,KIR2DL3,KIR2DP1,KIR2DL1,KIR3DP1,KIR2DL4,KIR3DL1,KIR2DS4,KIR3DL2\",False\n"
        "PAZ-003,A*03:01;A*03:01,B*07:02;B*08:01,C*07:01;C*07:01,A*24:02;A*02:01,B*57:01;B*44:02,C*06:02;C*07:01,\"KIR3DL3,KIR2DS2,KIR2DL2,KIR2DP1,KIR2DL1,KIR3DP1,KIR2DL4,KIR3DS1,KIR2DS1,KIR3DL2\",False\n"
        "PAZ-004,A*24:02;A*02:01,B*57:01;B*44:02,C*02:02;C*04:01,A*03:01;A*11:01,B*07:02;B*08:01,C*07:01;C*08:01,\"KIR3DL3,KIR2DL3,KIR2DP1,KIR2DL1,KIR3DP1,KIR2DL4,KIR3DL1,KIR2DS4,KIR3DL2\",False\n"
    )

    col_btn_tpl, col_info_tpl = st.columns([1, 3])
    with col_btn_tpl:
        st.download_button(
            label="📥 Scarica Template CSV di Esempio",
            data=sample_template,
            file_name="template_coorte_kiris.csv",
            mime="text/csv"
        )
    with col_info_tpl:
        st.caption("Il template include 4 casi studio clinici pronti: aplotipi A/A e B/X, varianti C1/C2, Bw4 e status KIR2DS1.")

    f = st.file_uploader("Carica file dataset (.csv / .xlsx)", type=["csv", "xlsx"])

    if f is not None:
        try:
            if f.name.endswith(".csv"):
                df_raw = pd.read_csv(f)
            else:
                try:
                    df_raw = pd.read_excel(f)
                except Exception as ex:
                    st.error(f"Errore nella lettura del file Excel: {str(ex)}. Installa openpyxl oppure carica il file in formato .csv.")
                    st.stop()
        except Exception as e:
            st.error(f"Impossibile leggere il file: {str(e)}")
            st.stop()

        st.markdown(f"**Dataset caricato:** `{f.name}` ({len(df_raw)} righe)")
        with st.expander("👁️ Anteprima dati grezzi caricati (prime 5 righe)"):
            st.dataframe(df_raw.head(), use_container_width=True)

        cols_lower = {str(c).strip().lower(): c for c in df_raw.columns}

        def find_col(candidates):
            for cand in candidates:
                if cand in cols_lower:
                    return cols_lower[cand]
            return None

        col_id = find_col(['id', 'id_paziente', 'patient_id', 'paziente', 'codice', 'case_id'])
        col_rec_a = find_col(['rec_hla_a', 'rec_a', 'hla_a_rec', 'rec_a_alleles', 'paz_a', 'recipient_a'])
        col_rec_b = find_col(['rec_hla_b', 'rec_b', 'hla_b_rec', 'rec_b_alleles', 'paz_b', 'recipient_b'])
        col_rec_c = find_col(['rec_hla_c', 'rec_c', 'hla_c_rec', 'rec_c_alleles', 'paz_c', 'recipient_c'])
        col_don_a = find_col(['don_hla_a', 'don_a', 'hla_a_don', 'don_a_alleles', 'donor_a'])
        col_don_b = find_col(['don_hla_b', 'don_b', 'hla_b_don', 'don_b_alleles', 'donor_b'])
        col_don_c = find_col(['don_hla_c', 'don_c', 'hla_c_don', 'don_c_alleles', 'donor_c'])
        col_don_kir = find_col(['don_kir', 'kir_donatore', 'kir_genes', 'kir_don', 'donor_kir', 'kir'])
        col_kir_null = find_col(['kir3dl1_null', 'null_3dl1', 'kir3dl1_004', 'is_3dl1_null'])

        with st.expander("⚙️ Verifica e Associazione Colonne Dataset", expanded=not all([col_rec_a, col_rec_b, col_rec_c, col_don_a, col_don_b, col_don_c])):
            all_cols = ["-- Nessuna --"] + list(df_raw.columns)
            c1, c2, c3 = st.columns(3)
            with c1:
                col_id_sel = st.selectbox("Colonna ID Paziente", all_cols, index=all_cols.index(col_id) if col_id in all_cols else 0)
                col_rec_a_sel = st.selectbox("Colonna Ricevente HLA-A", all_cols, index=all_cols.index(col_rec_a) if col_rec_a in all_cols else 0)
                col_rec_b_sel = st.selectbox("Colonna Ricevente HLA-B", all_cols, index=all_cols.index(col_rec_b) if col_rec_b in all_cols else 0)
                col_rec_c_sel = st.selectbox("Colonna Ricevente HLA-C", all_cols, index=all_cols.index(col_rec_c) if col_rec_c in all_cols else 0)
            with c2:
                col_don_a_sel = st.selectbox("Colonna Donatore HLA-A", all_cols, index=all_cols.index(col_don_a) if col_don_a in all_cols else 0)
                col_don_b_sel = st.selectbox("Colonna Donatore HLA-B", all_cols, index=all_cols.index(col_don_b) if col_don_b in all_cols else 0)
                col_don_c_sel = st.selectbox("Colonna Donatore HLA-C", all_cols, index=all_cols.index(col_don_c) if col_don_c in all_cols else 0)
            with c3:
                col_don_kir_sel = st.selectbox("Colonna Geni KIR Donatore", all_cols, index=all_cols.index(col_don_kir) if col_don_kir in all_cols else 0)
                col_kir_null_sel = st.selectbox("Colonna KIR3DL1 Null (*004)", all_cols, index=all_cols.index(col_kir_null) if col_kir_null in all_cols else 0)

            col_id = None if col_id_sel == "-- Nessuna --" else col_id_sel
            col_rec_a = None if col_rec_a_sel == "-- Nessuna --" else col_rec_a_sel
            col_rec_b = None if col_rec_b_sel == "-- Nessuna --" else col_rec_b_sel
            col_rec_c = None if col_rec_c_sel == "-- Nessuna --" else col_rec_c_sel
            col_don_a = None if col_don_a_sel == "-- Nessuna --" else col_don_a_sel
            col_don_b = None if col_don_b_sel == "-- Nessuna --" else col_don_b_sel
            col_don_c = None if col_don_c_sel == "-- Nessuna --" else col_don_c_sel
            col_don_kir = None if col_don_kir_sel == "-- Nessuna --" else col_don_kir_sel
            col_kir_null = None if col_kir_null_sel == "-- Nessuna --" else col_kir_null_sel

        known_kir_genes = CEN_GENES + TEL_GENES
        binary_kir_cols = [c for c in df_raw.columns if c.strip().upper() in [g.upper() for g in known_kir_genes]]

        missing_required = []
        if not col_rec_a: missing_required.append("Ricevente HLA-A")
        if not col_rec_b: missing_required.append("Ricevente HLA-B")
        if not col_rec_c: missing_required.append("Ricevente HLA-C")
        if not col_don_a: missing_required.append("Donatore HLA-A")
        if not col_don_b: missing_required.append("Donatore HLA-B")
        if not col_don_c: missing_required.append("Donatore HLA-C")
        if not col_don_kir and len(binary_kir_cols) < 5:
            missing_required.append("Geni KIR Donatore")

        if missing_required:
            st.warning(f"⚠️ Mappatura incompleta. Seleziona le colonne mancanti: {', '.join(missing_required)}.")
        else:
            if st.button("🚀 Elabora Coorte Completa", type="primary"):
                results_batch = []
                progress_bar = st.progress(0.0)
                status_text = st.empty()
                tot_rows = len(df_raw)

                for idx, row in df_raw.iterrows():
                    current_id = str(row[col_id]) if col_id else f"Caso_{idx+1}"
                    status_text.text(f"Elaborazione riga {idx+1}/{tot_rows}: {current_id}...")

                    def clean_hla_cell(val):
                        if pd.isna(val): return ""
                        s = str(val).strip()
                        return s.replace(";", ",").replace("/", ",").replace(" ", "")

                    r_a_raw = clean_hla_cell(row[col_rec_a])
                    r_b_raw = clean_hla_cell(row[col_rec_b])
                    r_c_raw = clean_hla_cell(row[col_rec_c])
                    d_a_raw = clean_hla_cell(row[col_don_a])
                    d_b_raw = clean_hla_cell(row[col_don_b])
                    d_c_raw = clean_hla_cell(row[col_don_c])

                    don_kir_set = set()
                    if col_don_kir and pd.notna(row[col_don_kir]):
                        kir_str = str(row[col_don_kir]).replace(";", ",").replace("/", ",")
                        items = [k.strip().upper() for k in kir_str.split(",") if k.strip()]
                        for it in items:
                            for canonical in known_kir_genes:
                                if it == canonical.upper():
                                    don_kir_set.add(canonical)
                    elif binary_kir_cols:
                        for bcol in binary_kir_cols:
                            bval = str(row[bcol]).strip().upper()
                            if bval in ['1', 'TRUE', 'POS', 'POSITIVE', '+', 'SI', 'YES']:
                                for canonical in known_kir_genes:
                                    if bcol.strip().upper() == canonical.upper():
                                        don_kir_set.add(canonical)

                    kir_null_val = False
                    if col_kir_null and pd.notna(row[col_kir_null]):
                        n_str = str(row[col_kir_null]).strip().upper()
                        kir_null_val = n_str in ['1', 'TRUE', 'SI', 'YES', 'NULL', 'POS']

                    try:
                        res = analyze_single_pair(
                            r_a_raw, r_b_raw, r_c_raw,
                            d_a_raw, d_b_raw, d_c_raw,
                            don_kir_set,
                            kir3dl1_null_variant=kir_null_val
                        )

                        results_batch.append({
                            "ID": current_id,
                            "Stato_Analisi": "SUCCESS",
                            "Note_Errore": "",
                            "Idoneità_Donatore": res["recommendation"]["level"],
                            "Aplotipo_KIR": res["haplo"],
                            "B_Score": res["b_score"],
                            "Cen_Haplo": res["cen_h"],
                            "Tel_Haplo": res["tel_h"],
                            "Epitopi_Donatore": ", ".join(sorted(res["d_ligs"])),
                            "Epitopi_Ricevente": ", ".join(sorted(res["r_ligs"])),
                            "KIR2DS1_Status": res["k2ds1_edu_status"],
                            "GvH_Score": res["gvh"]["score"],
                            "GvH_Dettaglio": res["gvh"]["status"],
                            "HvG_Score": res["hvg"]["score"],
                            "HvG_Dettaglio": res["hvg"]["status"],
                            "GvL_Score": res["gvl"]["score"],
                            "GvL_Dettaglio": res["gvl"]["status"],
                            "Stringa_Excel": res["excel_str"]
                        })
                    except Exception as err_row:
                        results_batch.append({
                            "ID": current_id,
                            "Stato_Analisi": "ERROR",
                            "Note_Errore": str(err_row),
                            "Idoneità_Donatore": "ERRORE",
                            "Aplotipo_KIR": "N/D",
                            "B_Score": -1,
                            "Cen_Haplo": "N/D",
                            "Tel_Haplo": "N/D",
                            "Epitopi_Donatore": "",
                            "Epitopi_Ricevente": "",
                            "KIR2DS1_Status": "N/D",
                            "GvH_Score": -1,
                            "GvH_Dettaglio": "",
                            "HvG_Score": -1,
                            "HvG_Dettaglio": "",
                            "GvL_Score": -1,
                            "GvL_Dettaglio": "",
                            "Stringa_Excel": ""
                        })

                    progress_bar.progress((idx + 1) / tot_rows)

                status_text.success(f"Analisi batch completata su {tot_rows} casi!")
                df_results = pd.DataFrame(results_batch)

                # --- METRICHE DI SINTESI DELLA COORTE ---
                st.markdown("### 📈 Indicatori e Metriche di Sintesi della Coorte")
                df_valid = df_results[df_results["Stato_Analisi"] == "SUCCESS"]
                n_valid = len(df_valid)
                n_err = len(df_results) - n_valid

                kpi1, kpi2, kpi3, kpi4 = st.columns(4)
                with kpi1:
                    st.metric("Casi Elaborati con Successo", f"{n_valid} / {tot_rows}")
                with kpi2:
                    pct_high_gvl = round((df_valid['GvL_Score'] >= 4).sum() / max(1, n_valid) * 100, 1)
                    st.metric("Alto Potenziale GvL (Score ≥ 4)", f"{pct_high_gvl}%")
                with kpi3:
                    pct_bx = round((df_valid['Aplotipo_KIR'] == 'B/X').sum() / max(1, n_valid) * 100, 1)
                    st.metric("Donatori Aplotipo B/X", f"{pct_bx}%")
                with kpi4:
                    k2ds1_active_count = df_valid['GvL_Dettaglio'].str.contains("KIR2DS1", na=False).sum()
                    pct_k2ds1 = round(k2ds1_active_count / max(1, n_valid) * 100, 1)
                    st.metric("KIR2DS1 Attivo vs Target C2", f"{pct_k2ds1}%")

                if n_err > 0:
                    st.warning(f"⚠️ {n_err} casi hanno riscontrato errori di sintassi allelica. Consulta la colonna 'Note_Errore'.")

                # ==============================================================
                # GRAFICI INTERATTIVI DI COORTE ALTAIR
                # ==============================================================
                if n_valid > 0:
                    st.markdown("### 📊 Visualizzazioni Interattive di Coorte")
                    chart_col1, chart_col2 = st.columns(2)

                    with chart_col1:
                        st.markdown("**Distribuzione del Potenziale Anti-Leucemico (Score GvL)**")
                        hist_data = df_valid['GvL_Score'].value_counts().reset_index()
                        hist_data.columns = ['Score_GvL', 'Numero_Casi']
                        hist_data['Score_GvL'] = hist_data['Score_GvL'].astype(str)

                        bar_chart = alt.Chart(hist_data).mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4).encode(
                            x=alt.X('Score_GvL:O', title='Score GvL (1-5)'),
                            y=alt.Y('Numero_Casi:Q', title='Numero di Pazienti'),
                            color=alt.Color('Score_GvL:N', scale=alt.Scale(
                                domain=['1', '2', '3', '4', '5'],
                                range=['#94a3b8', '#60a5fa', '#3b82f6', '#10b981', '#059669']
                            ), legend=None),
                            tooltip=['Score_GvL', 'Numero_Casi']
                        ).properties(height=260)
                        st.altair_chart(bar_chart, use_container_width=True)

                    with chart_col2:
                        st.markdown("**Matrice Decisionale: Potenziale GvL vs Rischio Rigetto HvG**")
                        scatter_chart = alt.Chart(df_valid).mark_circle(size=120, opacity=0.8).encode(
                            x=alt.X('HvG_Score:Q', title='Rischio Rigetto HvG (1 = Minimo, 5 = Massimo)', scale=alt.Scale(domain=[0.5, 5.5])),
                            y=alt.Y('GvL_Score:Q', title='Potenziale Anti-Leucemico GvL (1-5)', scale=alt.Scale(domain=[0.5, 5.5])),
                            color=alt.Color('Idoneità_Donatore:N', scale=alt.Scale(
                                domain=['PRIORITÀ ELEVATA', 'FAVOREVOLE', 'STANDARD', 'SUBOTTIMALE'],
                                range=['#10b981', '#3b82f6', '#f59e0b', '#ef4444']
                            )),
                            tooltip=['ID', 'Idoneità_Donatore', 'GvL_Score', 'HvG_Score', 'Aplotipo_KIR']
                        ).properties(height=260)
                        st.altair_chart(scatter_chart, use_container_width=True)

                # --- TABELLA RISULTATI COMPLETI ---
                st.markdown("### 📋 Tabella Risultati Dettagliati")
                st.dataframe(df_results, use_container_width=True)

                # --- PULSANTI DI DOWNLOAD ---
                st.markdown("### 💾 Esportazione Risultati")
                col_csv, col_xlsx = st.columns(2)
                
                with col_csv:
                    csv_data = df_results.to_csv(index=False, sep=';', encoding='utf-8-sig')
                    st.download_button(
                        label="📥 Scarica Risultati Coorte (CSV per Excel)",
                        data=csv_data,
                        file_name=f"kiris_coorte_risultati_{f.name.split('.')[0]}.csv",
                        mime="text/csv"
                    )

                with col_xlsx:
                    try:
                        buffer = io.BytesIO()
                        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                            df_results.to_excel(writer, index=False, sheet_name="Analisi_Coorte_KIRis")
                        st.download_button(
                            label="📥 Scarica Risultati Coorte (.xlsx Excel)",
                            data=buffer.getvalue(),
                            file_name=f"kiris_coorte_risultati_{f.name.split('.')[0]}.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )
                    except Exception:
                        st.caption("Per l'export diretto in .xlsx è consigliato installare openpyxl.")
