# Intro to genomics

## DESeq2

## Seurat

## Updating the website

GitHub does not run the R code. Results are computed locally and saved in `_freeze/`.

After editing a page in `vignettes/` (`.Rmd` or `.qmd`):

1. Render it locally from the repo root, with Quarto on your PATH:
   `quarto render vignettes/<file>`
2. Commit the edited file **and** the updated `_freeze/` folder.
3. Push to `main`. GitHub assembles the pages and publishes them.

New pages in `vignettes/` appear automatically as tiles on the home page. Give each one a `description:` in its header, and add it to the sidebar in `_quarto.yml`.

Set `Sys.setenv(LANGUAGE = "en")` in a page's first chunk, otherwise R messages follow the system language.
