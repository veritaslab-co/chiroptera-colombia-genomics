# Script 03: Análisis de vacíos taxonómicos
# Artículo 1 - Vacíos genómicos murciélagos colombianos

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import os

# --- CONFIGURACIÓN ---
INPUT = os.path.expanduser("~/chiroptera_colombia/data/raw/chiroptera_colombia_metadata.csv")
FIG_DIR = os.path.expanduser("~/chiroptera_colombia/results/figures/")
TAB_DIR = os.path.expanduser("~/chiroptera_colombia/results/tables/")

df = pd.read_csv(INPUT)

# --- NORMALIZACIÓN DE MARCADORES ---
marcador_map = {
    'COX1': 'COI', 'cox1': 'COI', 'CO1': 'COI',
    'Cyt-b': 'cytb', 'CYTB': 'cytb', 'CYB': 'cytb',
    'ND1': 'ND1', 'ND2': 'ND2', 'Fgb': 'Fgb'
}
df['marcador_norm'] = df['marcador'].replace(marcador_map).replace('', 'Sin_anotar').fillna('Sin_anotar')

# --- CLASIFICACIÓN DE TIPO DE SECUENCIA ---
def clasificar_tipo(row):
    if len(str(row['accession'])) > 12:
        return 'WGS/Ensamblaje'
    elif row['longitud_bp'] > 10000:
        return 'Mitogenoma'
    elif row['marcador_norm'] in ['COI','cytb','ND1','ND2']:
        return 'Marcador_mitocondrial'
    elif row['marcador_norm'] == 'Fgb':
        return 'Marcador_nuclear'
    else:
        return 'Sin_clasificar'

df['tipo_secuencia'] = df.apply(clasificar_tipo, axis=1)

print("\n--- MARCADORES NORMALIZADOS ---")
print(df['marcador_norm'].value_counts().to_string())

print("\n--- TIPOS DE SECUENCIA ---")
print(df['tipo_secuencia'].value_counts().to_string())

print(f"\n{'='*50}")
print(f"RESUMEN GENERAL")
print(f"{'='*50}")
print(f"Total registros:        {len(df)}")
print(f"Especies únicas:        {df['organismo'].nunique()}")
print(f"Marcadores únicos:      {df['marcador'].nunique()}")
print(f"Registros sin país:     {(df['pais'] == '').sum() + df['pais'].isna().sum()}")
print(f"Registros sin marcador: {(df['marcador'] == '').sum() + df['marcador'].isna().sum()}")

# --- TABLA 1: Especies con más registros ---
print(f"\n{'='*50}")
print("TOP 15 ESPECIES POR NÚMERO DE SECUENCIAS")
print(f"{'='*50}")
top_especies = df['organismo'].value_counts().head(15)
print(top_especies.to_string())

# --- TABLA 2: Distribución por marcador ---
print(f"\n{'='*50}")
print("DISTRIBUCIÓN POR MARCADOR MOLECULAR")
print(f"{'='*50}")
marcadores = df['marcador'].replace('', 'Sin_anotar').fillna('Sin_anotar').value_counts()
print(marcadores.to_string())

# --- TABLA 3: Extraer género ---
df['genero'] = df['organismo'].str.split().str[0]
print(f"\n{'='*50}")
print("TOP 15 GÉNEROS POR NÚMERO DE SECUENCIAS")
print(f"{'='*50}")
top_generos = df['genero'].value_counts().head(15)
print(top_generos.to_string())

# --- INSPECCIÓN DE SIN_CLASIFICAR ---
sin_clasificar = df[df['tipo_secuencia'] == 'Sin_clasificar'][
    ['accession','organismo','longitud_bp','marcador_norm','tipo_secuencia']
]
print("\n--- REGISTROS SIN CLASIFICAR ---")
print(sin_clasificar.to_string())
sin_clasificar.to_csv(TAB_DIR + "sin_clasificar.csv", index=False)

# --- RECLASIFICACIÓN DE SIN_CLASIFICAR ---

# Eliminar contaminante europeo
df = df[df['accession'] != 'OX031246.1'].copy()
print(f"Registro OX031246.1 (Plecotus austriacus - Europa) eliminado.")

# Eliminar contaminante euro asiatico
df = df[df['organismo'] != 'Rhinolophus ferrumequinum'].copy()
print("Rhinolophus ferrumequinum (Paleártico) eliminado.")

# Reclasificar PX, PP, OM, AY sin marcador con longitud típica de cytb
mascara_cytb = (
    (df['tipo_secuencia'] == 'Sin_clasificar') &
    (df['marcador_norm'] == 'Sin_anotar') &
    (df['longitud_bp'] >= 400)  # ← cambiado de 700 a 400
)
df.loc[mascara_cytb, 'marcador_norm'] = 'cytb_probable'
df.loc[mascara_cytb, 'tipo_secuencia'] = 'Marcador_mitocondrial'

print(f"Registros reclasificados como cytb_probable: {mascara_cytb.sum()}")
print(f"\nDataset final: {len(df)} registros")
print(df['tipo_secuencia'].value_counts().to_string())

# Guardar dataset limpio
df.to_csv(TAB_DIR + "chiroptera_colombia_limpio.csv", index=False)
print(f"\nDataset limpio guardado.")

print(df[df['tipo_secuencia'] == 'Sin_clasificar'][
    ['accession','organismo','longitud_bp','marcador_norm']
].to_string())

# --- GUARDAR TABLAS ---
os.makedirs(TAB_DIR, exist_ok=True)
df['organismo'].value_counts().reset_index().rename(
    columns={'index':'especie','organismo':'n_secuencias'}).to_csv(
    TAB_DIR + "especies_por_n_secuencias.csv", index=False)
marcadores.reset_index().rename(
    columns={'index':'marcador','marcador':'n_registros'}).to_csv(
    TAB_DIR + "distribucion_marcadores.csv", index=False)
print(f"\nTablas guardadas en {TAB_DIR}")

# --- FIGURA 1: Top 20 especies ---
fig, ax = plt.subplots(figsize=(10, 8))
top20 = df['organismo'].value_counts().head(20)
sns.barplot(x=top20.values, y=top20.index, palette="viridis", ax=ax)
ax.set_xlabel("Número de secuencias", fontsize=12)
ax.set_ylabel("")
ax.set_title("Top 20 especies de murciélagos colombianos\npor número de secuencias en NCBI", fontsize=13)
plt.tight_layout()
os.makedirs(FIG_DIR, exist_ok=True)
plt.savefig(FIG_DIR + "fig1_top20_especies.png", dpi=300)
plt.close()
print("Figura 1 guardada.")

# --- FIGURA 2: Distribución de marcadores ---
fig, ax = plt.subplots(figsize=(8, 5))
marcadores_plot = df['marcador'].replace('', 'Sin anotar').fillna('Sin anotar').value_counts()
sns.barplot(x=marcadores_plot.values, y=marcadores_plot.index, palette="mako", ax=ax)
ax.set_xlabel("Número de registros", fontsize=12)
ax.set_ylabel("")
ax.set_title("Distribución de marcadores moleculares\nChiroptera Colombia — NCBI", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_DIR + "fig2_marcadores.png", dpi=300)
plt.close()
print("Figura 2 guardada.")

# --- FIGURA 3: Distribución de longitudes ---
fig, ax = plt.subplots(figsize=(8, 5))
df['longitud_bp'].hist(bins=40, color='steelblue', edgecolor='white', ax=ax)
ax.set_xlabel("Longitud (bp)", fontsize=12)
ax.set_ylabel("Número de secuencias", fontsize=12)
ax.set_title("Distribución de longitudes de secuencia\nChiroptera Colombia — NCBI", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_DIR + "fig3_longitudes.png", dpi=300)
plt.close()
print("Figura 3 guardada.")

print(f"\n{'='*50}")
print("ANÁLISIS COMPLETO")
print(f"{'='*50}")