# Intro to genomics

## Updating the website

GitHub does not run the R code. Results are computed locally and saved in `_freeze/`.

After editing an `.Rmd`:

1. Render it locally (from this folder, with Quarto on your PATH):
   `quarto render vignettes/DESeq2_workflow.Rmd`
2. Commit the `.Rmd` **and** the updated `_freeze/` folder.
3. Push to `main`. GitHub assembles the pages and publishes them.

Each page's first code chunk should set `Sys.setenv(LANGUAGE = "en")`, otherwise R messages follow the system language (Polish on this Mac when rendered from a terminal).

To add a page, list the new file in the sidebar in `_quarto.yml`.
