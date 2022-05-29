#' Calculate phylogenetic diversity for a raster
#'
#' @param pres_bin_stack a raster of presence-absence. It can be an object of class 'raster' or 'SpatRaster'
#' @param tree an object of class 'phylo'
#'
#' @return SpatRaster
#' @export
#'
#' @examples
rast_pd <- function(pres_bin_stack, tree){
  {
    spatial_names <- as.character(names(pres_bin_stack))   # nomes espaciais e nomes filogeneticos
    tree_phy4 <- phylobase::phylo4(tree) # criei essa arvore de classe phylo4 para poder extrair as tipLabels
    labels <- as.character(phylobase::tipLabels(tree_phy4))
    on_tree <- intersect(spatial_names,labels)   # especies do raster que estao na arvore
    subtree <- ape::keep.tip(tree, on_tree) # arvore com o subconjunto de especies
    stack_reord <- pres_bin_stack[[subtree[["tip.label"]]]] # reordenar o stack de especies de acordo com a ordem das especies na arvore
    stack_reord_t <- terra::rast(stack_reord) # lendo como um spatraster para se adequar ao pacote Terra
    species_name <- names(stack_reord)
    subtree <- phylobase::phylo4(subtree)
    branch_length <- as.numeric(phylobase::edgeLength(subtree, 1:phylobase::nTips(subtree)))
  }
  vec_pd <- function(pres_bin_stack, branch_length, species_names){
    pres_bin_stack[is.na(pres_bin_stack)] <- 0 # atribui 0 a tudo que eh NA
    # para poder somar as presencas
    if(sum(pres_bin_stack)== 0) { # retorna NA quando a soma de stack_rast eh 0
      return(c(NA,NA))
    }
    pres_spp <- pres_bin_stack == 1 # seleciona as especies presentes
    n_species <- length(species_name) # riqueza
    blt <- sum(branch_length[pres_spp]) # div filogenetica: soma do branch_length de especies presentes
    blr <- blt/sum(pres_spp) # div filogenetica relativa a riqueza
    return(c(blt = blt, blr = blr)) # retornando os objetos
  }

  rpd <- terra::app(terra::rast(stack_reord), fun = vec_pd, branch_length = branch_length, species_name = species_name)
  return(rpd)
}
