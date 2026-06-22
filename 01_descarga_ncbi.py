# Script 01: Descarga de secuencias de Chiroptera Colombia desde NCBI
# Artículo 1 - Vacíos genómicos murciélagos colombianos
# Autor: Diego Alfonso Melo Escobar

from Bio import Entrez, SeqIO
import pandas as pd
import time
import os

# --- CONFIGURACIÓN ---
Entrez.email = "veritaslab.latam@gmail.com"  # ← cambia esto
OUTPUT_DIR = os.path.expanduser("~/chiroptera_colombia/data/raw/")
QUERY = "Chiroptera[Organism] AND Colombia[Country]"

# --- PASO 1: Contar registros disponibles ---
def contar_registros(query):
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

# --- PASO 2: Obtener IDs ---
def obtener_ids(query, total):
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=total)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# --- PASO 3: Descargar registros en lotes ---
def descargar_registros(ids, output_file, batch_size=100):
    import time
    with open(output_file, "w") as out:
        for i in range(0, len(ids), batch_size):
            lote = ids[i:i+batch_size]
            intentos = 0
            while intentos < 5:
                try:
                    handle = Entrez.efetch(db="nucleotide", id=lote,
                                           rettype="gb", retmode="text")
                    out.write(handle.read())
                    handle.close()
                    print(f"  Descargados {min(i+batch_size, len(ids))}/{len(ids)}")
                    time.sleep(1)
                    break
                except Exception as e:
                    intentos += 1
                    print(f"  Error en lote {i} — reintento {intentos}/5: {e}")
                    time.sleep(3 * intentos)
            else:
                print(f"  LOTE FALLIDO sin recuperación: registros {i} a {i+batch_size}")

# --- PASO 4: Extraer metadatos a tabla ---
def extraer_metadatos(genbank_file):
    registros = []
    for rec in SeqIO.parse(genbank_file, "genbank"):
        organismo = rec.annotations.get("organism", "")
        longitud = len(rec.seq)
        fecha = rec.annotations.get("date", "")
        pais = ""
        marcador = ""
        for feature in rec.features:
            if feature.type == "source":
                pais = feature.qualifiers.get("country", [""])[0]
            if feature.type == "gene":
                marcador = feature.qualifiers.get("gene", [""])[0]
        registros.append({
            "accession": rec.id,
            "organismo": organismo,
            "longitud_bp": longitud,
            "marcador": marcador,
            "pais": pais,
            "fecha": fecha
        })
    return pd.DataFrame(registros)

# --- EJECUCIÓN PRINCIPAL ---
if __name__ == "__main__":
    print("Contando registros en NCBI...")
    total = contar_registros(QUERY)
    print(f"Total encontrado: {total} registros")

    print("Obteniendo IDs...")
    ids = obtener_ids(QUERY, total)

    print("Descargando registros GenBank...")
    gb_file = OUTPUT_DIR + "chiroptera_colombia_raw.gb"
    descargar_registros(ids, gb_file)

    print("Extrayendo metadatos...")
    df = extraer_metadatos(gb_file)
    df.to_csv(OUTPUT_DIR + "chiroptera_colombia_metadata.csv", index=False)
    print(f"Tabla guardada: {len(df)} registros")
    print(df.head())