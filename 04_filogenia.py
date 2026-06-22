# Script 04: Pipeline filogenético
# Artículo 1 - Vacíos genómicos murciélagos colombianos

import pandas as pd
import subprocess
import os
from Bio import Entrez, SeqIO
import time

# --- CONFIGURACIÓN ---
Entrez.email = "veritaslab.latam@gmail.com"  # ← tu email real
TAB_DIR = os.path.expanduser("~/chiroptera_colombia/results/tables/")
SEQ_DIR = os.path.expanduser("~/chiroptera_colombia/data/processed/")
RES_DIR = os.path.expanduser("~/chiroptera_colombia/results/")
os.makedirs(SEQ_DIR, exist_ok=True)

# --- CARGAR DATASET ---
df = pd.read_csv(TAB_DIR + "chiroptera_colombia_limpio.csv")

# --- FILTRAR POR MARCADOR Y CALIDAD ---
def filtrar_secuencias(df, marcadores, longitud_min, longitud_max, nombre):
    subset = df[
        (df['marcador_norm'].isin(marcadores)) &
        (df['longitud_bp'] >= longitud_min) &
        (df['longitud_bp'] <= longitud_max)
    ].copy()
    print(f"\n{nombre}: {len(subset)} secuencias, {subset['organismo'].nunique()} especies")
    return subset

# Árbol A — cytb (fragmentos >= 400 bp, <= 1200 bp)
df_cytb = filtrar_secuencias(df,
    marcadores=['cytb', 'cytb_probable'],
    longitud_min=400, longitud_max=1200,
    nombre="Árbol A (cytb)")

# Árbol B — COI (fragmentos >= 400 bp, <= 800 bp)
df_coi = filtrar_secuencias(df,
    marcadores=['COI'],
    longitud_min=400, longitud_max=800,
    nombre="Árbol B (COI)")

# --- DESCARGAR SECUENCIAS DE NCBI ---
def descargar_fasta(ids, output_file, batch_size=50):
    if os.path.exists(output_file):
        print(f"  Ya existe: {output_file} — omitiendo descarga")
        return
    with open(output_file, "w") as out:
        for i in range(0, len(ids), batch_size):
            lote = ids[i:i+batch_size]
            intentos = 0
            while intentos < 5:
                try:
                    handle = Entrez.efetch(db="nucleotide", id=lote,
                                          rettype="fasta", retmode="text")
                    out.write(handle.read())
                    handle.close()
                    print(f"  Descargados {min(i+batch_size, len(ids))}/{len(ids)}")
                    time.sleep(0.5)
                    break
                except Exception as e:
                    intentos += 1
                    print(f"  Reintento {intentos}/5: {e}")
                    time.sleep(3 * intentos)

print("\nDescargando secuencias FASTA...")
descargar_fasta(df_cytb['accession'].tolist(),
               SEQ_DIR + "cytb_raw.fasta")
descargar_fasta(df_coi['accession'].tolist(),
               SEQ_DIR + "coi_raw.fasta")

# --- LIMPIAR HEADERS FASTA ---
def limpiar_headers(input_fasta, output_fasta, df_meta):
    """Reemplaza headers largos de NCBI por: accession_genero_especie"""
    acc_to_org = dict(zip(df_meta['accession'], df_meta['organismo']))
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            acc = rec.id.split('.')[0]
            # buscar accession con versión también
            org = acc_to_org.get(rec.id, acc_to_org.get(acc, "Unknown"))
            org_clean = org.replace(' ', '_').replace('/', '_')
            rec.id = f"{rec.id}_{org_clean}"
            rec.description = ""
            out.write(f">{rec.id}\n{str(rec.seq)}\n")
    print(f"  Headers limpios: {output_fasta}")

print("\nLimpiando headers...")
limpiar_headers(SEQ_DIR + "cytb_raw.fasta",
               SEQ_DIR + "cytb_clean.fasta", df_cytb)
limpiar_headers(SEQ_DIR + "coi_raw.fasta",
               SEQ_DIR + "coi_clean.fasta", df_coi)

# --- ALINEAMIENTO CON MAFFT ---
def correr_mafft(input_fasta, output_fasta):
    print(f"\nAlineando con MAFFT: {os.path.basename(input_fasta)}")
    cmd = f"mafft --auto --thread -1 {input_fasta} > {output_fasta}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        n = sum(1 for line in open(output_fasta) if line.startswith('>'))
        print(f"  Alineamiento completo: {n} secuencias")
    else:
        print(f"  ERROR MAFFT: {result.stderr[:200]}")

correr_mafft(SEQ_DIR + "cytb_clean.fasta",
            SEQ_DIR + "cytb_aligned.fasta")
correr_mafft(SEQ_DIR + "coi_clean.fasta",
            SEQ_DIR + "coi_aligned.fasta")

# --- ÁRBOL CON IQ-TREE ---
def correr_iqtree(input_fasta, prefix):
    print(f"\nConstruyendo árbol IQ-TREE: {os.path.basename(input_fasta)}")
    cmd = (f"iqtree2 -s {input_fasta} "
           f"--prefix {prefix} "
           f"-B 1000 "
           f"-T AUTO "
           f"--redo")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  Árbol generado: {prefix}.treefile")
    else:
        print(f"  ERROR IQ-TREE: {result.stderr[:300]}")

correr_iqtree(SEQ_DIR + "cytb_aligned.fasta",
             RES_DIR + "arbol_cytb")
correr_iqtree(SEQ_DIR + "coi_aligned.fasta",
             RES_DIR + "arbol_coi")

print("\n=== PIPELINE FILOGENÉTICO COMPLETO ===")
print(f"Árboles en: {RES_DIR}")
print("Archivos .treefile listos para visualizar en R con ggtree")
EOF