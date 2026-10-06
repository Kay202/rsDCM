.onAttach <- function(libname, pkgname) {
  message <- c("\n Welcome to rsDCM: robust and sparse group dynamic causal modeling.",
               "\n \n Website: https://github.com/Kay202/rsDCM/",
               "\n Bug report: https://github.com/Kay202/rsDCM/issues")
  packageStartupMessage(message)
}
