#para crear
# usethis::use_r("name.R")

## Load recount3 R package
library("recount3")

## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()
