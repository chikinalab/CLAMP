# Fetches GMT gene set libraries from Enrichr and saves them as .rds files
# in inst/extdata/ for use in vignettes.
#
# Usage: Rscript inst/scripts/fetch_vignette_gmts.R
library(CLAMP)

out_dir <- file.path("inst", "extdata")
enrichr_url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary"
libs <- c(
    "CellMarker_2024",
    "KEGG_2021_Human",
    "GO_Biological_Process_2025",
    "GTEx_Tissues_V8_2023",
    "MSigDB_Hallmark_2020",
    "Diabetes_Perturbations_GEO_2022"
)

for (name in libs) {
    message("Fetching ", name, "...")
    url <- paste0(enrichr_url, "?mode=text&libraryName=", name)
    gmt <- getGMT(url, name = name)
    out_path <- file.path(out_dir, paste0(name, ".rds"))
    saveRDS(gmt, out_path)
    message("Saved to ", out_path, " (", length(gmt), " gene sets)")
}