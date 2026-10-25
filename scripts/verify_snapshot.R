# Authenticate the compressed raw inputs before they are read.
# Hashes establish byte integrity, not scientific accuracy or origin authenticity.
verify_power_snapshot <- function(root = ".", snapshot_dir = NULL, quiet = FALSE) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("Install digest first: source('install_dependencies.R').")
  manifest <- read.csv(file.path(root, "data", "snapshot_manifest.csv"),
                       stringsAsFactors = FALSE)
  stopifnot(all(c("path", "sha256") %in% names(manifest)))
  if (!is.null(snapshot_dir))
    manifest <- manifest[dirname(manifest$path) == snapshot_dir, , drop = FALSE]
  if (!nrow(manifest)) stop("No snapshot manifest entries match the requested data.")
  paths <- file.path(root, manifest$path)
  if (any(!file.exists(paths))) stop("A file listed in the snapshot manifest is missing.")
  actual <- vapply(paths, function(p)
    digest::digest(file = p, algo = "sha256"), character(1))
  bad <- actual != manifest$sha256
  if (any(bad)) stop("Snapshot SHA-256 mismatch: ", paste(manifest$path[bad], collapse = ", "))
  if (!quiet) message("Snapshot SHA-256 verified: ", nrow(manifest), " files.")
  invisible(nrow(manifest))
}
