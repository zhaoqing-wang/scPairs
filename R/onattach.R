.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion("scPairs")

  banner <- paste(
    c(
      "           ____       _           ",
      "  ___  ___|  _ \\ __ _(_)_ __ ___  ",
      " / __|/ __| |_) / _` | | '__/ __| ",
      " \\__ \\ (__|  __/ (_| | | |  \\__ \\ ",
      " |___/\\___|_|   \\__,_|_|_|  |___/ "
    ),
    collapse = "\n"
  )

  title <- "scPairs: Identifying Synergistic Gene Pairs in single-Cell and Spatial Transcriptomics."
  url <- "https://github.com/zhaoqing-wang/scPairs"
  use_crayon <- requireNamespace("crayon", quietly = TRUE)

  if (use_crayon) {
    citation <- paste0(
      "Please cite: Wang Z (2026). ",
      crayon::italic(title),
      " R package version",
      crayon::bold(paste0(" ", pkg_version, ".")),
      " Available at: ",
      url
    )
    msg <- paste0(crayon::cyan$bold(banner), "\n", citation)
  } else {
    citation <- paste0(
      "Please cite: Wang Z (2026). ",
      title,
      " R package version ",
      pkg_version,
      ". Available at: ",
      url
    )
    msg <- paste0(banner, "\n", citation)
  }

  packageStartupMessage(msg)
}
