.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    msg <- paste0(
      "Please cite: Wang Z (2026). ",
      crayon::italic(
        "scPairs: Go Beyond Marker Genes - Discover Synergistic Gene Pairs in scRNA-seq and Spatial Maps."
      ),
      " R package version", crayon::bold(" 0.1.0."),
      " Available at: https://github.com/zhaoqing-wang/scPairs"
    )
  } else {
    msg <- paste0(
      "Please cite: Wang Z (2026). ",
      "scPairs: Go Beyond Marker Genes - Discover Synergistic Gene Pairs in scRNA-seq and Spatial Maps. ",
      "R package version 0.1.0. Available at: https://github.com/zhaoqing-wang/scPairs"
    )
  }

  packageStartupMessage(msg)
}
