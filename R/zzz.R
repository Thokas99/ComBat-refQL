.onLoad <- function(...) S7::methods_register()

.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")
  packageStartupMessage(paste0(
    "   ___          ___       __           ____  __\n",
    "  / __|___ _ __| _ ) __ _| |_  ___ _ _/ __ \\/ /\n",
    " | (__/ _ \\ '  \\ _ \\/ _` |  _|/ -_) '_/ /_/ / /__\n",
    "  \\___\\___/_|_|_|___/\\__,_|\\__|\\___|_| \\___\\_\\___/\n\n",
    "                    ComBat-refQL ", version
  ))
}
