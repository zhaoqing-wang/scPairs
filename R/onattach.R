.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion("scPairs")
  if (requireNamespace("crayon", quietly = TRUE)) {
    msg <- paste0(
      "Please cite: Wang Z (2026). ",
      crayon::italic(
        "scPairs: Identifying Synergistic Gene Pairs in single-Cell and Spatial Transcriptomics."
      ),
      " R package version", crayon::bold(paste0(" ", pkg_version, ".")),
      " Available at: https://github.com/zhaoqing-wang/scPairs"
    )
  } else {
    msg <- paste0(
      "Please cite: Wang Z (2026). ",
      "scPairs: Identifying Synergistic Gene Pairs in single-Cell and Spatial Transcriptomics. ",
      "R package version ", pkg_version, ". Available at: https://github.com/zhaoqing-wang/scPairs"
    )
  }

  packageStartupMessage(msg)
}
