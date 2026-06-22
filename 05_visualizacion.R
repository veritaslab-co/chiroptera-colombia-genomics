# Script 05: Visualización de árboles filogenéticos
# Artículo 1 - Vacíos genómicos murciélagos colombianos

library(ggtree)
library(ape)
library(ggplot2)

RES_DIR <- path.expand("~/chiroptera_colombia/results/")

# --- FUNCIÓN DE VISUALIZACIÓN ---
visualizar_arbol <- function(treefile, titulo, output_png) {
  
  tree <- read.tree(treefile)
  
  # Limpiar etiquetas: extraer género y especie del header
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  tree$tip.label <- sub("^[^ ]+ ", "", tree$tip.label)  # quitar accession
  # Acortar a Género especie solamente
  tree$tip.label <- sapply(strsplit(tree$tip.label, " "), function(x) {
    if(length(x) >= 2) paste(x[1], x[2]) else paste(x, collapse=" ")
  })
  
  # Convertir bootstrap a numérico
  tree$node.label <- as.numeric(tree$node.label)
  
  p <- ggtree(tree, layout="rectangular", 
              color="gray30", linewidth=0.4) +
    
    # Etiquetas de tips
    geom_tiplab(size=2.5, fontface="italic", 
                color="black", offset=0.001) +
    
    # Soporte bootstrap — solo nodos con soporte >= 70
    geom_nodelab(aes(label=ifelse(!is.na(as.numeric(label)) & 
                                    as.numeric(label) >= 70,
                                  as.numeric(label), "")),
                 size=2, color="firebrick", hjust=1.2, vjust=-0.3) +
    
    # Escala
    geom_treescale(x=0, y=-1, fontsize=3, linesize=0.5) +
    
    # Título
    ggtitle(titulo) +
    
    theme_tree2() +
    theme(
      plot.title = element_text(size=13, face="bold", hjust=0.5),
      plot.margin = margin(10, 120, 10, 10)
    )
  
  # Guardar
  ggsave(output_png, plot=p, width=12, height=14, dpi=300)
  cat("Guardado:", output_png, "\n")
  
  return(p)
}

# --- ÁRBOL CYTB ---
p_cytb <- visualizar_arbol(
  treefile  = paste0(RES_DIR, "arbol_cytb.treefile"),
  titulo    = "Filogenia de murciélagos colombianos (citocromo b)\nN = 118 secuencias, 40 especies — NCBI 2026",
  output_png = paste0(RES_DIR, "figures/fig4_arbol_cytb.png")
)

# --- ÁRBOL COI ---
p_coi <- visualizar_arbol(
  treefile  = paste0(RES_DIR, "arbol_coi.treefile"),
  titulo    = "Filogenia de murciélagos colombianos (COI)\nN = 110 secuencias, 16 especies — NCBI 2026",
  output_png = paste0(RES_DIR, "figures/fig4_arbol_coi.png")
)

cat("\n=== VISUALIZACIÓN COMPLETA ===\n")
cat("Figuras en:", paste0(RES_DIR, "figures/"), "\n")