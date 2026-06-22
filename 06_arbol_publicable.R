# Script 06: Árboles filogenéticos publicables
# Artículo 1 - Vacíos genómicos murciélagos colombianos

library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)

RES_DIR <- path.expand("~/chiroptera_colombia/results/")
FIG_DIR <- paste0(RES_DIR, "figures/")

# --- PALETA POR FAMILIA ---
colores_familia <- c(
  "Phyllostomidae" = "#2166AC",
  "Vespertilionidae" = "#D6604D",
  "Molossidae" = "#4DAC26",
  "Thyropteridae" = "#8B008B",
  "Desmodontidae" = "#FF8C00",
  "Unknown" = "gray50"
)

# --- TABLA DE FAMILIAS ---
familias <- data.frame(
  genero = c("Artibeus","Dermanura","Carollia","Glossophaga","Micronycteris",
             "Phyllostomus","Phylloderma","Sturnira","Vampyressa","Uroderma",
             "Ametrida","Gardnerycteris","Mesophylla","Platyrrhinus","Trachops",
             "Lophostoma","Chiroderma",
             "Myotis","Eptesicus","Neoeptesicus","Histiotus","Lasiurus",
             "Molossus","Molossops","Cynomops","Eumops",
             "Thyroptera",
             "Desmodus","Vampyrum","Tonatia"),
  familia = c(rep("Phyllostomidae", 17),
              rep("Vespertilionidae", 5),
              rep("Molossidae", 4),
              "Thyropteridae",
              rep("Desmodontidae", 3))
)

# --- FUNCIÓN PRINCIPAL ---
arbol_publicable <- function(treefile, titulo, output_png, ancho=14, alto=16) {
  
  tree <- read.tree(treefile)
  
  # Limpiar labels originales
  labels_limpios <- gsub("_", " ", tree$tip.label)
  labels_limpios <- sub("^[^ ]+ ", "", labels_limpios)
  labels_limpios <- sapply(strsplit(labels_limpios, " "), function(x) {
    if(length(x) >= 2) paste(x[1], x[2]) else paste(x, collapse=" ")
  })
  
  # Contar por especie
  conteo <- table(labels_limpios)
  
  # Podar a una secuencia por especie
  especies_unicas <- unique(labels_limpios)
  tips_mantener <- sapply(especies_unicas, function(sp) which(labels_limpios == sp)[1])
  tips_eliminar <- setdiff(seq_along(tree$tip.label), tips_mantener)
  tree_podado <- if(length(tips_eliminar) > 0) drop.tip(tree, tree$tip.label[tips_eliminar]) else tree
  
  # Labels limpios del árbol podado
  labs <- gsub("_", " ", tree_podado$tip.label)
  labs <- sub("^[^ ]+ ", "", labs)
  labs <- sapply(strsplit(labs, " "), function(x) {
    if(length(x) >= 2) paste(x[1], x[2]) else paste(x, collapse=" ")
  })
  
  # Labels con conteo
  labs_n <- sapply(labs, function(sp) {
    n <- conteo[sp]
    if(!is.na(n) && n > 1) paste0(sp, " (n=", n, ")") else sp
  })
  
  # Asignar familia por género ANTES de cambiar tip.label
  generos <- sapply(strsplit(labs, " "), `[`, 1)
  familias_asignadas <- sapply(generos, function(g) {
    idx <- match(g, familias$genero)
    if(!is.na(idx)) familias$familia[idx] else "Unknown"
  })
  
  # Ahora sí cambiar tip.label
  tree_podado$tip.label <- unname(labs_n)
  
  # Bootstrap
  tree_podado$node.label <- suppressWarnings(as.numeric(tree_podado$node.label))
  
  # Construir ggtree con labels ya limpios
  p <- ggtree(tree_podado, layout="rectangular", color="gray40", linewidth=0.5)
  
  # Extraer coordenadas AHORA — labels ya son los finales
  tip_coords <- p$data[p$data$isTip == TRUE, ]
  
  # Unir familia por label
  familia_df <- data.frame(
    label = unname(labs_n),
    familia = unname(familias_asignadas),
    stringsAsFactors = FALSE
  )
  tip_coords <- merge(tip_coords, familia_df, by="label", all.x=TRUE)
  tip_coords$familia[is.na(tip_coords$familia)] <- "Unknown"
  
  # Verificar
  cat("Familias asignadas:\n")
  print(table(tip_coords$familia))
  
  p <- p +
    geom_text(data=tip_coords,
              aes(x=x + 0.002, y=y, label=label, color=familia),
              size=3, fontface="italic", hjust=0) +
    geom_nodelab(aes(label=ifelse(!is.na(label) & label >= 70, label, "")),
                 size=2.2, color="gray20", hjust=1.3, vjust=-0.3) +
    geom_nodepoint(aes(subset=!is.na(label) & label >= 95),
                   size=1.5, color="black", alpha=0.7) +
    scale_color_manual(values=colores_familia, name="Familia",
                       guide=guide_legend(override.aes=list(label="■", size=5))) +
    geom_treescale(x=0, y=-1, fontsize=3, linesize=0.5) +
    ggtitle(titulo) +
    theme_tree2() +
    theme(
      plot.title = element_text(size=13, face="bold", hjust=0.5),
      legend.position = "bottom",
      legend.text = element_text(size=10),
      legend.title = element_text(size=11, face="bold"),
      plot.margin = margin(10, 220, 10, 10)
    )
  
  ggsave(output_png, plot=p, width=ancho, height=alto, dpi=300)
  cat("Guardado:", output_png, "\n")
  return(p)
}

# --- EJECUTAR ---
p_cytb <- arbol_publicable(
  treefile   = paste0(RES_DIR, "arbol_cytb.treefile"),
  titulo     = "Filogenia de murciélagos colombianos (citocromo b)\nN = 118 secuencias, 40 especies — NCBI 2026",
  output_png = paste0(FIG_DIR, "fig5_arbol_cytb_publicable.png"),
  ancho=16, alto=18
)

p_coi <- arbol_publicable(
  treefile   = paste0(RES_DIR, "arbol_coi.treefile"),
  titulo     = "Filogenia de murciélagos colombianos (COI)\nN = 110 secuencias, 16 especies — NCBI 2026",
  output_png = paste0(FIG_DIR, "fig5_arbol_coi_publicable.png"),
  ancho=16, alto=10
)

# --- ÁRBOLES ENRAIZADOS ---
p_cytb_enraizado <- arbol_publicable(
  treefile   = paste0(RES_DIR, "arbol_cytb_enraizado.treefile"),
  titulo     = "Filogenia de murciélagos colombianos (citocromo b)\nN = 118 secuencias, 40 especies — Outgroup: Pteropus vampyrus",
  output_png = paste0(FIG_DIR, "fig6_arbol_cytb_enraizado.png"),
  ancho=16, alto=18
)

p_coi_enraizado <- arbol_publicable(
  treefile   = paste0(RES_DIR, "arbol_coi_enraizado.treefile"),
  titulo     = "Filogenia de murciélagos colombianos (COI)\nN = 110 secuencias, 16 especies — Outgroup: Pteropus vampyrus",
  output_png = paste0(FIG_DIR, "fig6_arbol_coi_enraizado.png"),
  ancho=16, alto=10
)

cat("\n=== FIGURAS PUBLICABLES COMPLETAS ===\n")
