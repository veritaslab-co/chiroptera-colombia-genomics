# Script 02: Descarga y cruce con BOLD Systems (API v5)
# Artículo 1 - Vacíos genómicos murciélagos colombianos

import requests
import pandas as pd
import os

# --- CONFIGURACIÓN ---
TAB_DIR = os.path.expanduser("~/chiroptera_colombia/results/tables/")
RAW_DIR = os.path.expanduser("~/chiroptera_colombia/data/raw/")
NCBI_LIMPIO = TAB_DIR + "chiroptera_colombia_limpio.csv"
BOLD_RAW = RAW_DIR + "bold_chiroptera_colombia_raw.tsv"

# --- CARGAR DATOS ---
df_ncbi = pd.read_csv(NCBI_LIMPIO)
df_bold = pd.read_csv(BOLD_RAW, sep='\t', low_memory=False)

print(f"Registros NCBI: {len(df_ncbi)}")
print(f"Registros BOLD: {len(df_bold)}")
print(f"\nColumnas BOLD disponibles:")
print(list(df_bold.columns))

# --- EXPLORAR CAMPOS RELEVANTES DE BOLD ---
campos_interes = ['processid', 'species', 'genus', 'family',
                  'country/ocean', 'province/state', 'coord',
                  'marker_code', 'nuc', 'inst']

campos_disponibles = [c for c in campos_interes if c in df_bold.columns]
print(f"\nCampos encontrados: {campos_disponibles}")

df_bold_limpio = df_bold[campos_disponibles].copy()

# Longitud de secuencia
if 'nuc' in df_bold_limpio.columns:
    df_bold_limpio['longitud_bp'] = df_bold_limpio['nuc'].astype(str).apply(
        lambda x: len(x.replace('-','').replace('N','').replace('nan',''))
    )

# --- RESUMEN BOLD ---
print(f"\n{'='*50}")
print("RESUMEN BOLD")
print(f"{'='*50}")

col_especie = 'species' if 'species' in df_bold_limpio.columns else None
col_marcador = 'marker_code' if 'marker_code' in df_bold_limpio.columns else None

if col_especie:
    print(f"Especies únicas:  {df_bold_limpio[col_especie].nunique()}")
    print(f"\nTop 15 especies:")
    print(df_bold_limpio[col_especie].value_counts().head(15).to_string())

if col_marcador:
    print(f"\nMarcadores:")
    print(df_bold_limpio[col_marcador].value_counts().to_string())

# --- CRUCE NCBI vs BOLD ---
if col_especie:
    especies_ncbi = set(df_ncbi['organismo'].str.strip().dropna().unique())
    especies_bold = set(df_bold_limpio[col_especie].str.strip().dropna().unique())

    # Eliminar entradas vacías o ambiguas
    especies_bold = {e for e in especies_bold if len(str(e).split()) >= 2}

    solo_ncbi  = especies_ncbi - especies_bold
    solo_bold  = especies_bold - especies_ncbi
    en_ambas   = especies_ncbi & especies_bold

    print(f"\n{'='*50}")
    print("CRUCE NCBI vs BOLD")
    print(f"{'='*50}")
    print(f"Especies en NCBI:   {len(especies_ncbi)}")
    print(f"Especies en BOLD:   {len(especies_bold)}")
    print(f"En AMBAS bases:     {len(en_ambas)}")
    print(f"Solo en NCBI:       {len(solo_ncbi)}")
    print(f"Solo en BOLD:       {len(solo_bold)}")

    print(f"\n--- SOLO EN NCBI (sin barcoding) ---")
    for e in sorted(solo_ncbi):
        print(f"  {e}")

    print(f"\n--- SOLO EN BOLD (sin secuencia GenBank) ---")
    for e in sorted(solo_bold):
        print(f"  {e}")

    # Guardar tablas
    os.makedirs(TAB_DIR, exist_ok=True)
    pd.DataFrame(sorted(solo_ncbi), columns=['especie']).to_csv(
        TAB_DIR + "especies_solo_ncbi.csv", index=False)
    pd.DataFrame(sorted(solo_bold), columns=['especie']).to_csv(
        TAB_DIR + "especies_solo_bold.csv", index=False)
    pd.DataFrame(sorted(en_ambas), columns=['especie']).to_csv(
        TAB_DIR + "especies_en_ambas.csv", index=False)
    df_bold_limpio.to_csv(TAB_DIR + "bold_chiroptera_colombia_limpio.csv", index=False)
    print(f"\nTablas guardadas en {TAB_DIR}")

# --- TABLA MAESTRA DE COBERTURA ---
todas_ncbi = set(df_ncbi['organismo'].str.strip().dropna())
todas_bold = set(df_bold_limpio['species'].str.strip().dropna())
todas_bold = {e for e in todas_bold if len(str(e).split()) >= 2}

todas_especies = sorted(todas_ncbi | todas_bold)

resumen = []
for sp in todas_especies:
    resumen.append({
        'especie': sp,
        'en_ncbi': sp in todas_ncbi,
        'en_bold': sp in todas_bold,
        'cobertura': 'ambas' if (sp in todas_ncbi and sp in todas_bold)
                     else ('solo_ncbi' if sp in todas_ncbi else 'solo_bold')
    })

df_resumen = pd.DataFrame(resumen)
df_resumen.to_csv(TAB_DIR + "tabla_maestra_cobertura.csv", index=False)
print(f"\nTabla maestra guardada: {len(df_resumen)} especies totales")
print(df_resumen['cobertura'].value_counts().to_string())

print("\n=== BOLD COMPLETO ===")

