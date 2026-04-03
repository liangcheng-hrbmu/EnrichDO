.onAttach <- function(libname, pkgname) {
  .initial()
}
.initial <- function() {
  pos <- 1
  envir <- as.environment(pos)
  assign(".EnrichDOenv", new.env(), envir = envir)
  .EnrichDOenv <- get(".EnrichDOenv", envir = envir)

  tryCatch(utils::data(list = "doterms", package = "EnrichDO"))
  doterms <- get("doterms")
  .EnrichDOenv$doterms <- doterms
  rm(doterms, envir = .GlobalEnv)

  tryCatch(utils::data(list = "dotermgenes", package = "EnrichDO"))
  dotermgenes <- get("dotermgenes")
  .EnrichDOenv$dotermgenes <- dotermgenes
  rm(dotermgenes, envir = .GlobalEnv)

  .EnrichDOenv$enrich <- NULL
  .EnrichDOenv$doidCount <- NULL
}
