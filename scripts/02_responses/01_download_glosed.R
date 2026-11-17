rm(list = ls()); gc()

# download of the GloSED files used by the responses
#
# GloSED v1.0.0, the Global Standardised Soil Eukaryome Dataset (Mikryukov
# et al. 2026), Zenodo record 17827890, DOI 10.5281/zenodo.17827890, open
# access under CC-BY-4.0. Only the three files read by 03_responses.R are
# downloaded (about 230 MB), through the public Zenodo REST API (no token).
# Each file is checked against the md5 checksum published in the record; a
# file already present with the right checksum is not downloaded again.
#
# inputs: none (internet access)
# outputs (paths$glosed/): GloSED__OTU_table.parquet,
#   GloSED__Taxonomy.parquet, GloSED__SampleMetadata.xlsx
# supports: Methods (fungal data), Text S1

library(dplyr)
library(jsonlite)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

# settings ---------------------------------------------------------------------

zenodo_record <- "17827890"
zenodo_api <- "https://zenodo.org/api/records/"
glosed_files <- c("GloSED__OTU_table.parquet", "GloSED__Taxonomy.parquet", "GloSED__SampleMetadata.xlsx")
options(timeout = max(3600, getOption("timeout")))

ensure_dir(paths$glosed)

# record -----------------------------------------------------------------------

rec <- jsonlite::fromJSON(paste0(zenodo_api, zenodo_record), simplifyVector = FALSE)
message(rec$metadata$title, " | version ", rec$metadata$version, " | DOI ", rec$metadata$doi,
        " | license ", rec$metadata$license$id)
files <- dplyr::bind_rows(lapply(rec$files, function(f) {
  tibble::tibble(
    key = f$key,
    size = as.numeric(f$size),
    checksum = f$checksum,
    url = if (is.null(f$links$self)) {
      paste0(zenodo_api, zenodo_record, "/files/", f$key, "/content")
    } else {
      f$links$self
    }
  )
}))
miss <- setdiff(glosed_files, files$key)
if (length(miss) > 0) stop("not in the Zenodo record: ", paste(miss, collapse = ", "))
files <- files %>%
  dplyr::filter(key %in% glosed_files)
if (!all(grepl("^md5:", files$checksum))) stop("the record does not give md5 checksums")
files$md5 <- sub("^md5:", "", files$checksum)

# download and checksum --------------------------------------------------------

for (i in seq_len(nrow(files))) {
  dest <- file.path(paths$glosed, files$key[i])
  if (file.exists(dest) && identical(unname(tools::md5sum(dest)), files$md5[i])) {
    message(files$key[i], ": present, md5 ok")
    next
  }
  message(format(Sys.time(), "%H:%M:%S"), " ", files$key[i], ": downloading ",
          round(files$size[i] / 2^20, 1), " MB")
  tmp <- paste0(dest, ".part")
  utils::download.file(files$url[i], tmp, mode = "wb", quiet = TRUE)
  if (!identical(unname(tools::md5sum(tmp)), files$md5[i])) {
    file.remove(tmp)
    stop(files$key[i], ": md5 of the download differs from the record")
  }
  if (file.exists(dest)) file.remove(dest)
  file.rename(tmp, dest)
  message(files$key[i], ": md5 ok")
}

message("GloSED files: ", paths$glosed)
