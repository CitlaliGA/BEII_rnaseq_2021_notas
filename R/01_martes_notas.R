

## Para poder conectar tu compu con GitHub
usethis::create_github_token() ## Abrirá una página web, escoje un nombre único
## y luego da click en el botón verde al final. Después copia el token
## (son 40 caracteres)
gitcreds::gitcreds_set() ## Ojo, copia el token, no tu password de git!
## Si no, terminaras en la situación descrita en
## https://github.com/r-lib/usethis/issues/1347
## Tambien puedes usar
 usethis::edit_r_environ()
## y luego agregar las siguientes dos líneas (la línea en blanco es importante)
# GITHUB_PAT=TU_CLAVE_DE_40_LETRAS
#

## Configura tu usuario de GitHub
usethis::edit_git_config()
# [user]
#   name = Leonardo Collado Torres
#   email = lcolladotor@gmail.com

## Para inicializar el repositorio de Git
usethis::use_git()

## Para conectar tu repositorio local de Git con los servidores de GitHub
usethis::use_github()


#PUSH!!!!!!!!!!!!!!!!!!!!!!!!
gert::git_push()
#escribimos un nuevo archivo
writeLines("hola", "R/prueba.R")
#Por ejemplo podríamos probar añadir algo nuevo
gert::git_add("R/prueba.R")

#añadimos commit de lo que se hizo
gert::git_commit("prueba gert")

gert::git_log() #nos da info

gert::git_push() #sube tus cambios del repo local
