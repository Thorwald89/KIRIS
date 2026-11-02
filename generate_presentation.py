"""
Script per la generazione della presentazione PowerPoint professionale di KIRIS
Software: KIRIS (KIR Interpretation Software)
Formato: Widescreen 16:9 (.pptx)
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE

def create_presentation():
    prs = Presentation()
    # Formato 16:9 Widescreen (13.333 x 7.5 pollici)
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    # Palette Cromatica Ispirata a KIRIS Dark/Medical UI
    C_BG_DARK = RGBColor(15, 23, 42)       # #0f172a Deep Slate
    C_CARD_DARK = RGBColor(30, 41, 59)     # #1e293b Card Dark
    C_CARD_BORDER = RGBColor(51, 65, 85)   # #334155
    C_WHITE = RGBColor(248, 250, 252)      # #f8fafc
    C_TEXT_SUB = RGBColor(203, 213, 225)   # #cbd5e1
    C_TEXT_MUTED = RGBColor(148, 163, 184) # #94a3b8
    C_CYAN = RGBColor(56, 189, 248)        # #38bdf8 Accent Primary
    C_BLUE = RGBColor(59, 130, 246)        # #3b82f6 Inhibitory Blue
    C_GREEN = RGBColor(16, 185, 129)       # #10b981 Activating Green
    C_PURPLE = RGBColor(124, 58, 237)      # #7c3aed Framework Purple
    C_CORAL = RGBColor(239, 68, 68)        # #ef4444 Alert/Risk Red
    C_AMBER = RGBColor(245, 158, 11)       # #f59e0b Standard Amber

    blank_layout = prs.slide_layouts[6]

    def set_slide_background(slide, color=C_BG_DARK):
        bg_shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, prs.slide_width, prs.slide_height)
        bg_shape.fill.solid()
        bg_shape.fill.fore_color.rgb = color
        bg_shape.line.fill.background()
        return bg_shape

    def add_header(slide, title_text, category_text="KIRIS — IMMUNOGENETICS TOOLKIT"):
        # Category Pill
        cat_box = slide.shapes.add_textbox(Inches(0.8), Inches(0.4), Inches(11.7), Inches(0.4))
        tf_c = cat_box.text_frame
        tf_c.word_wrap = True
        p_c = tf_c.paragraphs[0]
        p_c.text = category_text.upper()
        p_c.font.size = Pt(11)
        p_c.font.bold = True
        p_c.font.color.rgb = C_CYAN

        # Title
        title_box = slide.shapes.add_textbox(Inches(0.8), Inches(0.75), Inches(11.7), Inches(0.8))
        tf_t = title_box.text_frame
        tf_t.word_wrap = True
        p_t = tf_t.paragraphs[0]
        p_t.text = title_text
        p_t.font.size = Pt(24)
        p_t.font.bold = True
        p_t.font.color.rgb = C_WHITE

        # Divider Line
        div = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.8), Inches(1.55), Inches(11.733), Inches(0.02))
        div.fill.solid()
        div.fill.fore_color.rgb = C_CARD_BORDER
        div.line.fill.background()

    def add_card(slide, left, top, width, height, title, body_bullets, accent_color=C_CYAN, title_icon=""):
        # Card Background
        card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
        card.fill.solid()
        card.fill.fore_color.rgb = C_CARD_DARK
        card.line.color.rgb = C_CARD_BORDER
        card.line.width = Pt(1.5)

        # Accent Bar on Top
        bar = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, Inches(0.08))
        bar.fill.solid()
        bar.fill.fore_color.rgb = accent_color
        bar.line.fill.background()

        # Text Box
        pad = Inches(0.25)
        tb = slide.shapes.add_textbox(left + pad, top + pad, width - (pad * 2), height - (pad * 2))
        tf = tb.text_frame
        tf.word_wrap = True
        tf.margin_top = 0
        tf.margin_left = 0
        tf.margin_right = 0

        # Title
        p_t = tf.paragraphs[0]
        p_t.text = f"{title_icon} {title}".strip()
        p_t.font.size = Pt(16)
        p_t.font.bold = True
        p_t.font.color.rgb = accent_color
        p_t.space_after = Pt(10)

        # Body
        for bullet in body_bullets:
            p = tf.add_paragraph()
            p.text = f"• {bullet}"
            p.font.size = Pt(12)
            p.font.color.rgb = C_TEXT_SUB
            p.space_after = Pt(6)

    # --------------------------------------------------------------------------
    # SLIDE 1: COVER SLIDE
    # --------------------------------------------------------------------------
    s1 = prs.slides.add_slide(blank_layout)
    set_slide_background(s1, C_BG_DARK)

    # Glow Box
    glow = s1.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(0.8), Inches(1.2), Inches(11.733), Inches(5.1))
    glow.fill.solid()
    glow.fill.fore_color.rgb = C_CARD_DARK
    glow.line.color.rgb = C_CYAN
    glow.line.width = Pt(2)

    # Cover Text Box
    c_box = s1.shapes.add_textbox(Inches(1.5), Inches(1.8), Inches(10.3), Inches(3.8))
    c_tf = c_box.text_frame
    c_tf.word_wrap = True

    p0 = c_tf.paragraphs[0]
    p0.text = "DNA & NK CELL IMMUNOGENETICS TOOLKIT"
    p0.font.size = Pt(14)
    p0.font.bold = True
    p0.font.color.rgb = C_CYAN
    p0.space_after = Pt(12)

    p1 = c_tf.add_paragraph()
    p1.text = "🧬 KIRIS"
    p1.font.size = Pt(46)
    p1.font.bold = True
    p1.font.color.rgb = C_WHITE

    p2 = c_tf.add_paragraph()
    p2.text = "KIR Interpretation Software"
    p2.font.size = Pt(22)
    p2.font.bold = True
    p2.font.color.rgb = C_TEXT_SUB
    p2.space_after = Pt(18)

    p3 = c_tf.add_paragraph()
    p3.text = "Motore di Regole Cliniche Deterministico per la Valutazione dell'Alloreattività NK & la Selezione del Donatore nel Trapianto Allogenico / Aploidentico"
    p3.font.size = Pt(13)
    p3.font.color.rgb = C_TEXT_MUTED
    p3.space_after = Pt(24)

    p4 = c_tf.add_paragraph()
    p4.text = "Conformità ISO 13485 (Deterministic Calculations)  |  Release 0.4  |  Dark & Light Theme  |  High-Throughput Batch Engine"
    p4.font.size = Pt(11)
    p4.font.bold = True
    p4.font.color.rgb = C_GREEN

    # --------------------------------------------------------------------------
    # SLIDE 2: EXECUTIVE SUMMARY & OBIETTIVI CLINICI
    # --------------------------------------------------------------------------
    s2 = prs.slides.add_slide(blank_layout)
    set_slide_background(s2)
    add_header(s2, "Executive Summary & Valore Clinico", "Panoramica della Piattaforma")

    add_card(s2, Inches(0.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Il Bisogno Clinico",
             [
                 "Nel trapianto di midollo osseo (HSCT) e cellule staminali emopoietiche, le cellule NK svolgono un ruolo chiave nella sopravvivenza del paziente.",
                 "La mancata inibizione o l'attivazione inappropriata possono provocare rigetto dell'innesto (HvG) o recidiva leucemica.",
                 "Spesso la valutazione immunogenetica HLA-KIR è complessa, soggetta ad errori manuali e priva di standardizzazione deterministica."
             ], C_CYAN, "🎯")

    add_card(s2, Inches(4.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "La Soluzione: KIRIS",
             [
                 "KIRIS automatizza l'analisi molecolare completa della coppia Ricevente-Donatore secondo standard IMGT/HLA.",
                 "Modella formalmente il Rheostat Model dell'educazione funzionale delle NK e il dosaggio di KIR2DS1.",
                 "Calcola vettori clinici quantitativi (GvH, HvG, GvL) e un indice sintetico di idoneità del donatore validato secondo consensus EBMT/GITMO."
             ], C_GREEN, "🧬")

    add_card(s2, Inches(8.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Impatto Operativo",
             [
                 "Azzeramento dei tempi di interpretazione: analisi istantanea in tempo reale (singola coppia o coorti massive).",
                 "Eliminazione di modelli stocastici 'black-box': ogni output è matematicamente tracciabile e verificabile.",
                 "Visualizzazioni ad alto impatto: Radar Chart SVG geopolare, mappa cromosomica a colori funzionali ed export Excel/CSV."
             ], C_PURPLE, "⚡")

    # --------------------------------------------------------------------------
    # SLIDE 3: ARCHITETTURA & CONFORMITÀ DETERMINISTICA ISO 13485
    # --------------------------------------------------------------------------
    s3 = prs.slides.add_slide(blank_layout)
    set_slide_background(s3)
    add_header(s3, "Architettura Software & Validazione ISO 13485", "Affidabilità & Certificazione")

    add_card(s3, Inches(0.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Principi di Determinismo Assoluto",
             [
                 "Nessun modello probabilistico, euristico o stocastico: due esecuzioni con gli stessi dati generano al 100% lo stesso risultato numerico.",
                 "Piena conformità ai principi di validazione ISO 13485 per software medicali di supporto diagnostico/decisionale.",
                 "Test suite integrata automatica con copertura al 100% su casistiche cliniche critiche (residui Bw4, eccezioni C*16, KIR2DS1).",
                 "Isolamento funzionale: motore computazionale puro ('src/kir_analysis.py') disaccoppiato dall'interfaccia utente."
             ], C_CYAN, "🏛️")

    add_card(s3, Inches(6.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Robustezza & Tolleranza dell'Input",
             [
                 "Regex IMGT/HLA standardizzata: validazione formale con controllo di locus e coerenza genica.",
                 "Normalizzazione automatica intelligente: correzione istantanea di sintassi di laboratorio frequenti (es. '03:01' -> 'A*03:01', 'B44:02' -> 'B*44:02').",
                 "Verifica di Coerenza Biologica KIR: allarmi in tempo reale per geni framework mancanti o combinazioni aplotipiche anomale.",
                 "Resilienza cloud: gestione difensiva dei tipi e ricaricamento automatico per garantire continuità di servizio."
             ], C_BLUE, "🛡️")

    # --------------------------------------------------------------------------
    # SLIDE 4: MAPPATURA EPITOPI HLA-C, HLA-B E HLA-A
    # --------------------------------------------------------------------------
    s4 = prs.slides.add_slide(blank_layout)
    set_slide_background(s4)
    add_header(s4, "Mappatura Molecolare HLA ➡️ Ligandi KIR", "Regole Immunogenetiche")

    add_card(s4, Inches(0.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Locus HLA-C (C1 vs C2)",
             [
                 "Gruppo C1 (Ser77/Asn80): ligando per KIR2DL2/3. Include C*01, 03, 07, 08, 12, 14, 16.",
                 "Gruppo C2 (Asn77/Lys80): ligando per KIR2DL1 e KIR2DS1. Include C*02, 04, 05, 06, 15, 17, 18.",
                 "Eccezioni Alleliche Gestite: C*16:01 e C*16:04 sono C1; C*16:02 è C2. C*12:02/03 sono C1; C*12:04 è C2.",
                 "Fornisce la base per il calcolo del Missing-Self C1/C2 e per l'educazione di KIR2DS1."
             ], C_BLUE, "🔹")

    add_card(s4, Inches(4.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Locus HLA-B & A (Bw4 Residuo 80)",
             [
                 "Bw4-I80 (Isoleucina-80): massima affinità per KIR3DL1. Include B*51, 57, 58, 27, 37 e alleli HLA-A (A*24, A*23, A*32).",
                 "Bw4-T80 (Treonina-80): affinità intermedia per KIR3DL1. Include B*38, B*44 (eccetto 44:06).",
                 "Bw4-L80 (Leucina-80 / B*13): eccezione anergica per KIR3DL1; non conferisce licensing funzionale.",
                 "Bw6: neutro, privo di recettori KIR inibitori o attivatori corrispondenti."
             ], C_CYAN, "🔸")

    add_card(s4, Inches(8.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Epitopi Ausiliari (Locus A)",
             [
                 "HLA-A*03 e HLA-A*11: ligandi cogniti per il recettore inibitorio KIR3DL2.",
                 "Integrazione nel calcolo dell'alloreattività ausiliaria in caso di missing self da HLA-A.",
                 "Mappatura bidirezionale: estrazione simultanea sia degli epitopi del Donatore che del Ricevente."
             ], C_PURPLE, "🔹")

    # --------------------------------------------------------------------------
    # SLIDE 5: IL RHEOSTAT MODEL DELL'EDUCAZIONE DELLE NK
    # --------------------------------------------------------------------------
    s5 = prs.slides.add_slide(blank_layout)
    set_slide_background(s5)
    add_header(s5, "Modello di Educazione Funzionale (Rheostat Model)", "Biologia Molecolare del Licensing")

    add_card(s5, Inches(0.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "Licensing Forte",
             [
                 "High-Affinity",
                 "Bw4-Ile80 + KIR3DL1: massima degranulazione ed espressione di Perforina/Granzyme B.",
                 "C2 + KIR2DL1: interazione solida ad alta affinità (Lisina-80) nel sistema HLA-C.",
                 "Generano il massimo potenziale di lisi in caso di missing-self nel ricevente."
             ], C_GREEN, "🟢")

    add_card(s5, Inches(3.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "Licensing Moderato",
             [
                 "Intermediate-Affinity",
                 "Bw4-Thr80 + KIR3DL1: stabilità di legame intermedia, educazione parziale.",
                 "C1 + KIR2DL2: interazione molecolare a media affinità per il gruppo C1.",
                 "Eroga un potenziale alloreattivo intermedio controllato."
             ], C_AMBER, "🟡")

    add_card(s5, Inches(6.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "Licensing Debole",
             [
                 "Low-Affinity",
                 "C1 + KIR2DL3: affinità di legame significativamente inferiore rispetto a 2DL2.",
                 "KIR2DL2 cross-reattivo con C2: segnale residuo debole.",
                 "Alloreattività basale lieve in caso di missing-self."
             ], C_AMBER, "🟠")

    add_card(s5, Inches(9.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "Anergico / Null",
             [
                 "No Licensing",
                 "KIR3DL1*004: variante intrappolata nel reticolo endoplasmatico (ER-trapped Null).",
                 "Bw4-Leu80 (B*13): mancata inibizione.",
                 "Bw6 (No KIR).",
                 "Incapacità di mediare alloreattività anti-leucemica protettiva."
             ], C_CORAL, "🔴")

    # --------------------------------------------------------------------------
    # SLIDE 6: MODELLO AVANZATO KIR2DS1 (ATTIVATORIO)
    # --------------------------------------------------------------------------
    s6 = prs.slides.add_slide(blank_layout)
    set_slide_background(s6)
    add_header(s6, "Recettore Attivatorio KIR2DS1: Educazione & Target", "Modello Fauriat & Venstrom")

    add_card(s6, Inches(0.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Stratificazione dell'Educazione nel Donatore",
             [
                 "Donatore C1/C1 (Non Tollerizzato / Hyperexcitable):",
                 "  In assenza del ligando self C2, le NK KIR2DS1+ non subiscono desensibilizzazione. Se esposte a C2, mostrano un'ipereccitabilità citotossica potentissima.",
                 "Donatore C1/C2 (Parzialmente Tollerizzato / Responsivo):",
                 "  In presenza di una singola dose di C2, una frazione rilevante di NK KIR2DS1+ conserva responsività funzionale ed efficacia clinica.",
                 "Donatore C2/C2 (Anergico / Tollerizzato):",
                 "  L'esposizione al doppio dosaggio genico di C2 induce un profondo tuning-down cronico: le NK KIR2DS1+ sono iporesponsive/anergiche."
             ], C_GREEN, "⚡")

    add_card(s6, Inches(6.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Vincolo Indispensabile: Target C2 nel Ricevente",
             [
                 "KIR2DS1 riconosce specificamente l'epitopo HLA-C del gruppo C2 sulle cellule bersaglio (blasti leucemici).",
                 "Condizione Necessaria per Boost GvL:",
                 "  Affinché KIR2DS1 esplichi l'effetto anti-leucemico, il RICEVENTE DEVE ESPRIMERE C2 (eterozigote C1/C2 o omozigote C2/C2).",
                 "Se Ricevente è C1/C1 (C2-negativo):",
                 "  Nessun boost GvL viene assegnato, poiché il target tumorale manca del ligando attivatore cognito.",
                 "Eliminato il vecchio errore clinico di assegnare boost a donatori con riceventi C1/C1."
             ], C_CYAN, "🎯")

    # --------------------------------------------------------------------------
    # SLIDE 7: I TRE VETTORI CLINICI QUANTITATIVI (1-5)
    # --------------------------------------------------------------------------
    s7 = prs.slides.add_slide(blank_layout)
    set_slide_background(s7)
    add_header(s7, "I Tre Vettori Clinici Quantitativi (Score 1-5)", "Scoring Deterministico Multidimensionale")

    add_card(s7, Inches(0.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "⚡ Vettore GvH",
             [
                 "Alloreattività Donatore ➡️ Ospite",
                 "Valuta il Missing-Self inibitorio ponderato sull'affinità (Bw4-I80 e C2 peso 2.0; Bw4-T80 e C1 peso 1.2; C1 debole peso 0.8).",
                 "Integra l'alloreattività attivatoria di KIR2DS1 (se paziente è C2+).",
                 "Score 1: Tollerante / Nullo",
                 "Score 2-3: Basso / Moderato",
                 "Score 4-5: Alto / Massima Alloreattività NK."
             ], C_AMBER, "⚡")

    add_card(s7, Inches(4.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "🛡️ Vettore HvG",
             [
                 "Rischio Rigetto Innesto (Host-vs-Graft)",
                 "Quantifica il missing-self inverso: ligandi del ricevente assenti nel donatore (le NK dell'ospite riconoscono il graft come 'estraneo').",
                 "Pesatura graduata: C2 mancante (2.0), Bw4 mancante (2.0), C1 mancante (1.2).",
                 "Score 1: Graft Protetto",
                 "Score 2-3: Basso / Moderato Rischio",
                 "Score 4-5: Massimo Rischio Rigetto (Graft Failure)."
             ], C_CORAL, "🛡️")

    add_card(s7, Inches(8.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "🎯 Vettore GvL",
             [
                 "Potenziale Anti-Leucemico",
                 "Stima l'eradicazione dei blasti residui e la prevenzione della recidiva leucemica.",
                 "Combinazione cooperativa:",
                 "  + Missing-self alta affinità (+1.5)",
                 "  + Multi-locus missing-self (+0.5)",
                 "  + KIR2DS1 attivo vs target C2 (+1.5)",
                 "  + Cooley B-Score (fino a +1.5).",
                 "Score 1-5: da Minimo a Massimo Effetto GvL."
             ], C_GREEN, "🎯")

    # --------------------------------------------------------------------------
    # SLIDE 8: VALUTAZIONE GEOMETRICA — RADAR CHART GEOPOLARE
    # --------------------------------------------------------------------------
    s8 = prs.slides.add_slide(blank_layout)
    set_slide_background(s8)
    add_header(s8, "Valutazione Geometrica: Radar Chart Geopolare", "Visualizzazione dei Vettori")

    add_card(s8, Inches(0.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Architettura del Radar Geopolare",
             [
                 "Proiezione tridimensionale a 120° su scala standardizzata 1-5: Asse GvH (Top), Asse GvL (Bottom-Right), Asse HvG (Bottom-Left).",
                 "Griglia tecnica a 5 livelli concentrici con poligoni sfumati e numerazione graduata.",
                 "Poligono dati luminescente con coordinate dinamiche e vertici rinforzati a cerchio pieno.",
                 "Color-coding clinico istantaneo del poligono:",
                 "  🟢 Verde Smeraldo: Alto GvL (≥4) e Basso HvG (≤2)",
                 "  🔴 Rosso Corallo: Alto Rischio Rigetto (HvG ≥4)",
                 "  🔵 Blu Elettrico: Assetto bilanciato standard."
             ], C_CYAN, "🕸️")

    add_card(s8, Inches(6.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Ingegneria di Rendering & UX",
             [
                 "Rendering isolato tramite canvas HTML nativo ('components.html'): azzerato qualsiasi glitch di sanificazione Markdown.",
                 "Nessun testo o codice grezzo visualizzabile a video: visualizzazione vettoriale pura e nitida a 60fps.",
                 "Piena integrazione con la modalità Dark / Black Mode (sfondo profondo ad alto contrasto) e Light Mode.",
                 "Badge esterni con arrotondamento e tipografia ad alta leggibilità direttamente agganciati ai vertici degli assi."
             ], C_BLUE, "🎨")

    # --------------------------------------------------------------------------
    # SLIDE 9: ALGORITMO DI SELEZIONE DONATORE (CONSENSUS EBMT / GITMO)
    # --------------------------------------------------------------------------
    s9 = prs.slides.add_slide(blank_layout)
    set_slide_background(s9)
    add_header(s9, "Algoritmo di Selezione del Donatore (Donor Recommendation)", "Supporto Decisionale Clinico")

    add_card(s9, Inches(0.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "PRIORITÀ ELEVATA",
             [
                 "⭐ Donatore Ottimale",
                 "Criteri:",
                 "GvL ≥ 4 AND HvG ≤ 2",
                 "Forte alloreattività anti-leucemica unita a un innesto altamente protetto da rigetto mediato da NK.",
                 "Candidato di prima scelta assoluta per il trapianto."
             ], C_GREEN, "⭐")

    add_card(s9, Inches(3.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "FAVOREVOLE",
             [
                 "🟢 Donatore Favorevole",
                 "Criteri:",
                 "GvL ≥ 3 AND HvG ≤ 3",
                 "Presenza di fattori alloreattivi protettivi e/o aplotipo B con rischio di rigetto moderato e gestibile.",
                 "Donatore pienamente raccomandato per la procedura."
             ], C_BLUE, "🟢")

    add_card(s9, Inches(6.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "STANDARD",
             [
                 "🟡 Donatore Neutro",
                 "Criteri:",
                 "GvL ≤ 2 AND HvG ≤ 2",
                 "Profilo convenzionale neutro: assenza di specifici mismatch alloreattivi avversi.",
                 "Il trapianto segue il decorso standard dei protocolli clinici."
             ], C_AMBER, "🟡")

    add_card(s9, Inches(9.8), Inches(1.8), Inches(2.7), Inches(5.0),
             "SUBOTTIMALE",
             [
                 "⚠️ Rischio Rigetto",
                 "Criteri:",
                 "HvG ≥ 4",
                 "Elevata reattività delle NK dell'ospite contro le cellule staminali del donatore.",
                 "Rischio aumentato di graft failure: considerare donatore alternativo o intensificare condizionamento."
             ], C_CORAL, "⚠️")

    # --------------------------------------------------------------------------
    # SLIDE 10: MAPPA CROMOSOMICA COMPARATIVA DEL LOCUS KIR
    # --------------------------------------------------------------------------
    s10 = prs.slides.add_slide(blank_layout)
    set_slide_background(s10)
    add_header(s10, "Mappa Cromosomica Comparativa del Locus KIR", "Genomica Funzionale")

    add_card(s10, Inches(0.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Struttura dei 17 Loci & Hotspot",
             [
                 "Disposizione fedele dei 17 geni in coordinate cromosomiche 19q13.4 da 5' Centromero a 3' Telomero.",
                 "Centromero (CEN): 3DL3, 2DS2, 2DL2, 2DL3, 2DS3, 2DP1, 2DL1, 3DP1.",
                 "Delimitatore Centrale Evidenziato: Hotspot di Ricombinazione Intergenica tra 3DP1 e 2DL4.",
                 "Telomero (TEL): 2DL4, 3DL1, 3DS1, 2DL5A, 2DL5B, 2DS5, 2DS1, 2DS4, 3DL2.",
                 "B-Score di Cooley: calcolo indipendente dei motivi CEN-B (0-2) e TEL-B (0-2), punteggio totale 0-4."
             ], C_PURPLE, "🧬")

    add_card(s10, Inches(6.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Classificazione Cromatica & Pulsanti Dinamici",
             [
                 "🔵 Inibitori: 2DL1, 2DL2, 2DL3, 3DL1, 3DL2, 2DL5A/B.",
                 "🟢 Attivatori: 2DS1, 2DS2, 2DS3, 2DS4, 2DS5, 3DS1.",
                 "🟣 Framework Costitutivi: 3DL3, 2DL4, 3DL2, 3DP1, 2DP1.",
                 "⚪ Pseudogeni: 2DP1, 3DP1 (marcatore Ψ).",
                 "Pulsanti Aplotipo Dinamici con Sincronizzazione Immediata: cliccando sui bottoni di aplotipo (A standard, B arricchito, tutti i loci), le 16 checkbox si aggiornano all'istante."
             ], C_CYAN, "🎛️")

    # --------------------------------------------------------------------------
    # SLIDE 11: ANALISI COORTE BATCH AD ALTO RENDIMENTO (TAB 2)
    # --------------------------------------------------------------------------
    s11 = prs.slides.add_slide(blank_layout)
    set_slide_background(s11)
    add_header(s11, "Analisi Coorte Batch ad Alto Rendimento (Tab 2)", "Epidemiologia & Studi Clinici")

    add_card(s11, Inches(0.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Capacità di Elaborazione Massiva",
             [
                 "Supporto a dataset in formato CSV ed Excel (.xlsx) contenenti coorti di centinaia di coppie Ricevente-Donatore.",
                 "Template CSV precompilato scaricabile con 1 clic comprensivo di casi studio clinici di riferimento.",
                 "Mappatura Intelligente delle Colonne: riconosce automaticamente intestazioni flessibili per Ricevente, Donatore e Loci KIR (elenco o colonne binarie 1/0).",
                 "Elaborazione resiliente: se una riga presenta alleli non validi, memorizza la nota senza bloccare il batch."
             ], C_BLUE, "📁")

    add_card(s11, Inches(6.8), Inches(1.8), Inches(5.7), Inches(5.0),
             "Metriche KPI & Grafici Interattivi Altair",
             [
                 "Indicatori di Coorte Immediati: % Alto GvL (≥4), % Aplotipi B/X, % KIR2DS1 attivo su target C2.",
                 "Istogramma di Distribuzione GvL: visualizza la frequenza dei livelli di efficacia anti-leucemica.",
                 "Matrice Decisionale GvL vs HvG (Scatter Plot): mappa tutti i donatori su piano cartesiano permettendo di isolare i candidati nel quadrante ottimale.",
                 "Esportazione completa dei dati arricchiti con codifica UTF-8 BOM e file Excel (.xlsx)."
             ], C_GREEN, "📈")

    # --------------------------------------------------------------------------
    # SLIDE 12: CONCLUSIONI & VALORE STRATEGICO
    # --------------------------------------------------------------------------
    s12 = prs.slides.add_slide(blank_layout)
    set_slide_background(s12)
    add_header(s12, "Conclusioni & Valore Strategico", "Riepilogo Piattaforma")

    add_card(s12, Inches(0.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Rigore Scientifico",
             [
                 "Modello reostato continuo integrato con la più recente letteratura ematologica.",
                 "Regole chiare e deterministiche per il comportamento di KIR2DS1 su target C2.",
                 "Superamento dei limiti interpretativi manuali nei centri trapianto."
             ], C_CYAN, "🔬")

    add_card(s12, Inches(4.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Design & Fruibilità",
             [
                 "Interfaccia utente moderna con switch Dark/Light per ambienti ospedalieri.",
                 "Visualizzazioni geometriche innovative per una comunicazione immediata con l'ematologo.",
                 "Modalità popup lightbox interattiva per la consultazione delle release note."
             ], C_PURPLE, "🎨")

    add_card(s12, Inches(8.8), Inches(1.8), Inches(3.7), Inches(5.0),
             "Pronto per l'Uso Clinico",
             [
                 "Conforme ai requisiti di riproducibilità ISO 13485.",
                 "Disponibile sia in esecuzione locale che su cloud.",
                 "Supporto decisionale rapido per migliorare l'outcome e ridurre le recidive nei pazienti trapiantati."
             ], C_GREEN, "🏆")

    # Salvataggio
    output_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "KIRIS_Presentazione_Funzionale.pptx"))
    prs.save(output_path)
    print(f"Presentazione generata con successo: {output_path}")
    return output_path

if __name__ == "__main__":
    create_presentation()
