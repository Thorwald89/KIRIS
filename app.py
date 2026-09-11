import streamlit as st
import pandas as pd
from src.kir_analysis import KIRLigandProfiler

# Configurazione della pagina Streamlit
st.set_page_config(
    page_title="KIRIS - KIR & HLA Immunogenetics Platform",
    page_icon="🧬",
    layout="wide"
)

# Inizializzazione del profiler
profiler = KIRLigandProfiler()

# Titolo e intestazione
st.title("🧬 KIRIS: KIR-HLA Immunogenetics Analysis Toolkit")
st.markdown("Piattaforma per l'analisi dei ligandi HLA, profiling KIR e simulazione dell'effetto Missing Self / GvL.")

# Barra laterale (Sidebar) per la navigazione
st.sidebar.header("Pannello di Controllo")
mode = st.sidebar.radio("Seleziona Modalità", ["Simulatore Singolo Paziente", "Analisi Dataset (CSV/Excel)"])

if mode == "Simulatore Singolo Paziente":
    st.subheader("Simulazione Interattiva Coppia KIR - Ligando HLA")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 1. Alleli HLA-C e Ligandi")
        res_allele1 = st.selectbox("Residuo Posizione 80 (Allele 1)", ["N (Asparagina - C1)", "K (Lisina - C2)", "Altro"])
        res_allele2 = st.selectbox("Residuo Posizione 80 (Allele 2)", ["N (Asparagina - C1)", "K (Lisina - C2)", "Altro"])
        
        # Estrazione della lettera singola
        code1 = res_allele1[0]
        code2 = res_allele2[0]
        
        group1 = profiler.determine_hla_c_group(code1)
        group2 = profiler.determine_hla_c_group(code2)
        
        st.info(f"**Profilo Ligandi HLA-C identificato:** {group1} / {group2}")
        
    with col2:
        st.markdown("### 2. Genotipo Recettori KIR")
        selected_kirs = st.multiselect(
            "Seleziona i recettori KIR espressi:",
            ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR3DL1", "KIR2DS1"],
            default=["KIR2DL1", "KIR2DL2"]
        )
        
    st.markdown("---")
    st.subheader("Risultato Valutazione Missing Self (Potenziale GvL)")
    
    if st.button("Esegui Analisi Interazione"):
        hla_ligands = [group1, group2]
        results = profiler.evaluate_missing_self(selected_kirs, hla_ligands)
        
        for interaction, status in results.items():
            if "Missing Self" in status:
                st.warning(f"**{interaction}:** {status} ⚡ (Elevata reattività NK / Effetto GvL)")
            else:
                st.success(f"**{interaction}:** {status} ✅ (Inibizione presente)")

elif mode == "Analisi Dataset (CSV/Excel)":
    st.subheader("Caricamento Dataset per Analisi di Coorte")
    uploaded_file = st.file_uploader("Carica un file di dati (.csv o .xlsx)", type=["csv", "xlsx"])
    
    if uploaded_file is not None:
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        else:
            df = pd.read_excel(uploaded_file)
            
        st.write("### Anteprima Dati Caricati")
        st.dataframe(df.head())
    else:
        st.info("Carica un file contenente i dati dei pazienti per avviare l'analisi collettiva.")
