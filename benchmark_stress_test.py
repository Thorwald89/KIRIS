"""
KIRIS Benchmark & Stress Test Suite
Misura la latenza, il throughput, il consumo di memoria e la scalabilità del motore deterministico KIRIS.
"""

import time
import random
import statistics
import tracemalloc
import sys
import os

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from src.kir_analysis import (
    validate_and_extract_alleles,
    map_hla_to_ligands,
    calculate_b_content_logic,
    assess_donor_education,
    calculate_vectors_quantitative,
    calculate_donor_recommendation,
    generate_radar_chart_svg,
    format_excel_row,
    analyze_single_pair,
    CEN_GENES,
    TEL_GENES,
    FRAMEWORK_GENES
)

# Pool di alleli reali per generazione dati sintetici ad alta variabilità
HLA_A_POOL = [
    "A*01:01", "A*02:01", "A*03:01", "A*11:01", "A*24:02", 
    "A*23:01", "A*25:01", "A*29:02", "A*30:01", "A*32:01", "A*68:01"
]
HLA_B_POOL = [
    "B*07:02", "B*08:01", "B*13:01", "B*13:02", "B*15:01", "B*15:13",
    "B*18:01", "B*27:05", "B*35:01", "B*38:01", "B*40:01", "B*44:02", 
    "B*44:03", "B*51:01", "B*57:01", "B*58:01"
]
HLA_C_POOL = [
    "C*01:02", "C*02:02", "C*03:01", "C*04:01", "C*05:01", "C*06:02",
    "C*07:01", "C*07:02", "C*08:01", "C*12:02", "C*12:04", "C*15:02",
    "C*16:01", "C*16:02", "C*16:04", "C*17:01"
]

ALL_KIR = CEN_GENES + TEL_GENES

def generate_random_patient_pair():
    r_a = f"{random.choice(HLA_A_POOL)}, {random.choice(HLA_A_POOL)}"
    r_b = f"{random.choice(HLA_B_POOL)}, {random.choice(HLA_B_POOL)}"
    r_c = f"{random.choice(HLA_C_POOL)}, {random.choice(HLA_C_POOL)}"

    d_a = f"{random.choice(HLA_A_POOL)}, {random.choice(HLA_A_POOL)}"
    d_b = f"{random.choice(HLA_B_POOL)}, {random.choice(HLA_B_POOL)}"
    d_c = f"{random.choice(HLA_C_POOL)}, {random.choice(HLA_C_POOL)}"

    # Assetto KIR: framework sempre presenti + geni variabili randomici
    variable_kir = [g for g in ALL_KIR if g not in FRAMEWORK_GENES]
    selected_variable = [g for g in variable_kir if random.random() > 0.4]
    don_kir = set(list(FRAMEWORK_GENES) + selected_variable)
    kir_null = random.random() < 0.15

    return r_a, r_b, r_c, d_a, d_b, d_c, don_kir, kir_null

def benchmark_single_execution_micro(iterations=1000):
    print("=" * 70)
    print("1. BENCHMARK MICRO-LATENZA (Analisi Singola Coppia)")
    print("=" * 70)

    # Genera 100 coppie pre-istanziate per escludere il tempo del random generator
    pairs = [generate_random_patient_pair() for _ in range(100)]

    latencies_ms = []
    # Warm-up (100 iterazioni)
    for i in range(100):
        p = pairs[i % 100]
        analyze_single_pair(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], is_dark=True)

    # Benchmark effettivo
    for i in range(iterations):
        p = pairs[i % 100]
        t0 = time.perf_counter()
        analyze_single_pair(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], is_dark=True)
        t1 = time.perf_counter()
        latencies_ms.append((t1 - t0) * 1000.0)

    avg_lat = statistics.mean(latencies_ms)
    med_lat = statistics.median(latencies_ms)
    std_lat = statistics.stdev(latencies_ms)
    min_lat = min(latencies_ms)
    max_lat = max(latencies_ms)
    p95_lat = statistics.quantiles(latencies_ms, n=100)[94]
    p99_lat = statistics.quantiles(latencies_ms, n=100)[98]
    throughput = 1000.0 / avg_lat

    print(f"Iterazioni misurate:   {iterations}")
    print(f"Tempo Medio:           {avg_lat:.4f} ms  ({avg_lat * 1000:.1f} microsecondi)")
    print(f"Tempo Mediano:         {med_lat:.4f} ms")
    print(f"Min / Max:             {min_lat:.4f} ms / {max_lat:.4f} ms")
    print(f"Percentile 95 (P95):   {p95_lat:.4f} ms")
    print(f"Percentile 99 (P99):   {p99_lat:.4f} ms")
    print(f"Deviazione Standard:   {std_lat:.4f} ms")
    print(f"Throughput Stimato:    {throughput:,.0f} coppie analizzate / secondo")
    print()

def benchmark_component_breakdown(iterations=1000):
    print("=" * 70)
    print("2. PROFILING DETTAGLIATO DEI SOTTOSISTEMI COMPUTAZIONALI")
    print("=" * 70)

    p = generate_random_patient_pair()
    rec_a, rec_b, rec_c, don_a, don_b, don_c, don_kir, kir_null = p

    # Profiling individuale delle fasi interne
    times = {
        "1. Validazione & Normalizzazione HLA (6 Loci)": [],
        "2. Mapping Epitopi & Ligandi (Donatore + Ricevente)": [],
        "3. Calcolo Aplotipo & Cooley B-Score": [],
        "4. Educazione NK Rheostat (Licensing + 2DS1)": [],
        "5. Calcolo Vettori Quantitativi (GvH, HvG, GvL)": [],
        "6. Sintetizzatore Decisionale EBMT/GITMO": [],
        "7. Generazione Radar Geopolare SVG ad Alta Risoluzione": [],
        "8. Formattazione Stringa Tabellare Excel": []
    }

    for _ in range(iterations):
        # 1. Validazione
        t0 = time.perf_counter()
        r_a = validate_and_extract_alleles(rec_a, "Rec HLA-A")
        r_b = validate_and_extract_alleles(rec_b, "Rec HLA-B")
        r_c = validate_and_extract_alleles(rec_c, "Rec HLA-C")
        d_a = validate_and_extract_alleles(don_a, "Don HLA-A")
        d_b = validate_and_extract_alleles(don_b, "Don HLA-B")
        d_c = validate_and_extract_alleles(don_c, "Don HLA-C")
        times["1. Validazione & Normalizzazione HLA (6 Loci)"].append((time.perf_counter() - t0) * 1000)

        # 2. Mapping
        t0 = time.perf_counter()
        r_ligs = map_hla_to_ligands(r_c, r_b, r_a)
        d_ligs = map_hla_to_ligands(d_c, d_b, d_a)
        times["2. Mapping Epitopi & Ligandi (Donatore + Ricevente)"].append((time.perf_counter() - t0) * 1000)

        # 3. Aplotipo
        t0 = time.perf_counter()
        haplo, b_score, cen_h, tel_h = calculate_b_content_logic(don_kir)
        times["3. Calcolo Aplotipo & Cooley B-Score"].append((time.perf_counter() - t0) * 1000)

        # 4. Educazione
        t0 = time.perf_counter()
        edu_res = assess_donor_education(don_kir, d_ligs, d_b, d_c, kir3dl1_null_variant=kir_null)
        times["4. Educazione NK Rheostat (Licensing + 2DS1)"].append((time.perf_counter() - t0) * 1000)

        k2ds1_edu_status = next((x['Edu'] for x in edu_res if x['KIR'] == 'KIR2DS1'), "N/A")

        # 5. Vettori
        t0 = time.perf_counter()
        gvh_sc, gvh_s, gvh_c, hvg_sc, hvg_s, hvg_c, gvl_sc, gvl_s, gvl_c = calculate_vectors_quantitative(
            edu_res, d_ligs, r_ligs, haplo, b_score, k2ds1_edu_status
        )
        times["5. Calcolo Vettori Quantitativi (GvH, HvG, GvL)"].append((time.perf_counter() - t0) * 1000)

        # 6. Recommendation
        t0 = time.perf_counter()
        recom = calculate_donor_recommendation(gvh_sc, hvg_sc, gvl_sc, haplo, b_score, is_dark=True)
        times["6. Sintetizzatore Decisionale EBMT/GITMO"].append((time.perf_counter() - t0) * 1000)

        # 7. Radar SVG
        t0 = time.perf_counter()
        radar_svg = generate_radar_chart_svg(gvh_sc, hvg_sc, gvl_sc, is_dark=True)
        times["7. Generazione Radar Geopolare SVG ad Alta Risoluzione"].append((time.perf_counter() - t0) * 1000)

        # 8. Excel string
        t0 = time.perf_counter()
        excel_str = format_excel_row(edu_res, haplo, b_score, gvh_sc, hvg_sc, gvl_sc)
        times["8. Formattazione Stringa Tabellare Excel"].append((time.perf_counter() - t0) * 1000)

    total_sum = sum(statistics.mean(v) for v in times.values())

    print(f"{'Componente Computazionale':<55} | {'Media (ms)':<10} | {'Quota %':<8}")
    print("-" * 78)
    for comp, vals in times.items():
        m = statistics.mean(vals)
        pct = (m / total_sum) * 100
        print(f"{comp:<55} | {m:8.4f} ms | {pct:6.1f} %")
    print("-" * 78)
    print(f"{'TEMPO TOTALE PURO SOMMATO':<55} | {total_sum:8.4f} ms | 100.0 %")
    print()

def benchmark_stress_scalability():
    print("=" * 70)
    print("3. STRESS-TEST DI SCALABILITÀ COORTE (Batch Processing)")
    print("=" * 70)

    cohort_sizes = [100, 500, 1000, 5000, 10000, 25000]

    print(f"{'Dimensione Coorte':<20} | {'Tempo Totale':<14} | {'Throughput':<18} | {'Memoria Picco':<14}")
    print("-" * 72)

    for n in cohort_sizes:
        # Pre-generazione coppie
        dataset = [generate_random_patient_pair() for _ in range(n)]

        tracemalloc.start()
        t0 = time.perf_counter()

        # Elaborazione batch sequenziale simulata
        for row in dataset:
            analyze_single_pair(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], is_dark=True)

        t1 = time.perf_counter()
        current_mem, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        elapsed = t1 - t0
        tp = n / elapsed
        peak_mb = peak_mem / (1024 * 1024)

        if elapsed < 1.0:
            elapsed_str = f"{elapsed * 1000:.1f} ms"
        else:
            elapsed_str = f"{elapsed:.2f} s"

        print(f"{n:<20,d} | {elapsed_str:<14} | {tp:14,.0f} casi/s | {peak_mb:10.2f} MB")
    print("-" * 72)
    print()

def benchmark_edge_cases_and_robustness():
    print("=" * 70)
    print("4. STRESS-TEST DI ROBUSTEZZA & CASI LIMITE")
    print("=" * 70)

    edge_cases = [
        ("Omozigosi Totale Ricevente & Donatore (A*02:01, B*08:01, C*07:01)",
         ("A*02:01,A*02:01", "B*08:01,B*08:01", "C*07:01,C*07:01", "A*02:01,A*02:01", "B*08:01,B*08:01", "C*07:01,C*07:01", set(FRAMEWORK_GENES), False)),
        
        ("Tutti i 17 Loci KIR Presenti (Donatore B/B Massimo)",
         ("A*03:01,A*11:01", "B*51:01,B*57:01", "C*02:02,C*04:01", "A*01:01,A*02:01", "B*44:02,B*13:01", "C*16:01,C*07:01", set(ALL_KIR), False)),
        
        ("Solo i 5 Geni Framework (Assetto Minimo Costitutivo)",
         ("A*03:01,A*11:01", "B*08:01,B*44:02", "C*07:01,C*05:01", "A*02:01,A*03:01", "B*44:02,B*07:02", "C*03:01,C*08:01", set(FRAMEWORK_GENES), False)),
        
        ("Input 'Sporco' da Normalizzare (Spazi multipli, prefissi mancanti, separatori misti)",
         ("03:01, 11:01", "B08:01 , 44:02", "C07:01, 05:01", "A02:01 , A*03:01", "B*44:02, 07:02", "03:01 , C*08:01", set(FRAMEWORK_GENES), True)),
        
        ("Caso HLA B*13 Anergico Leucina-80 + KIR3DL1*004 Null",
         ("A*02:01,A*03:01", "B*13:01,B*13:02", "C*07:01,C*08:01", "A*24:02,A*02:01", "B*13:01,B*07:02", "C*03:01,C*07:01", set(ALL_KIR), True))
    ]

    for name, args in edge_cases:
        t0 = time.perf_counter()
        try:
            res = analyze_single_pair(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], is_dark=True)
            t_exec = (time.perf_counter() - t0) * 1000
            print(f"[OK] PASSED: {name}")
            print(f"  Tempo: {t_exec:.4f} ms | GvH: {res['gvh']['score']}/5 | HvG: {res['hvg']['score']}/5 | GvL: {res['gvl']['score']}/5 | Idoneità: {res['recommendation']['level']}")
        except Exception as e:
            print(f"[FAIL] FAILED: {name} -> {str(e)}")
    print()

if __name__ == "__main__":
    print("\nAVVIO BENCHMARK UFFICIALE & STRESS-TEST KIRIS")
    print(f"Interprete Python: {sys.version.split()[0]} | Piattaforma: {sys.platform}\n")
    
    benchmark_single_execution_micro(iterations=2000)
    benchmark_component_breakdown(iterations=1000)
    benchmark_stress_scalability()
    benchmark_edge_cases_and_robustness()
    print("=" * 70)
    print("BENCHMARK COMPLETATO CON SUCCESSO.")
    print("=" * 70)
