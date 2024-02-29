.onAttach <- function(libname, pkgname) {

  version <- "1.0.0"


  loadmsg <- paste0(
    "\nSince the companion paper of the 'metaConvert' ",
    "package is under review",
    ", this version is still considered as a beta.",
    "\n\n -> Email us at cgosling@parisnanterre.fr if you notice any bug. <-")

  packageStartupMessage(loadmsg)

}
