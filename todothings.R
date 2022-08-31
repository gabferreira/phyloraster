# to do
# checar em todas as funçoes o inherit parameters
# usar inherit parameters com a funcao do pacote phylogrid
# na funcao range size colocar "x" no lugar de pres.rast
# o filename pode dar problemas com arquivos temporarios

# if (!is.null(filename)){ # to save the rasters when the path is provide
#   rend <- terra::writeRaster(rend, filename, ...)
# }

# ver funcao  pd ses e rodar pra ver como sao os objetos
# rowrise aleatoriza pres e ausencia das especies em cada sitios
# para raster devo aleatorizar a presenca/ausencia das especies no
# pixel ate o numero de celulas
# r[1] <- sample(r[1]) fazer isso de 1:ncell

### colwise sample para aleatorizar aonde a especie esta
## somo quantos zeros e faco uma probabilidade
# total = (soma(0)+soma(1))
# soma(1)/total ### prob1 prob0
# vamos pegar cada pixel de cada layer separadamente
# ao inves de retornar uma lista retornar um stack (c(rasts))

# hyperlink para funcoes


library(devtools)
