# Intro to genomics

## Updating the website

GitHub does not run the R code. Results are computed locally and saved in `_freeze/`.

After editing an `.Rmd`:

1. Render it locally (from this folder, with Quarto on your PATH):
   `quarto render DESeq2_workflow.Rmd`
2. Commit the `.Rmd` **and** the updated `_freeze/` folder.
3. Push to `main`. GitHub assembles the pages and publishes them.

To add a page, list the new file in the sidebar in `_quarto.yml`.
