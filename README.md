# ufg-transport-docs

Source for the UFG-Transport documentation site, built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/).

**Rendered site: <https://hglee531.github.io/ufg-transport-docs/>**

Code repository: <https://github.com/hglee531/ufg-transport>

## Local preview

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

Then browse http://127.0.0.1:8000.

## Deployment

Every push to `main` triggers `.github/workflows/docs.yml`, which runs `mkdocs gh-deploy --force` and publishes the built site to the `gh-pages` branch. GitHub Pages serves that branch at the URL above.

## Editing

All documentation content lives under `docs_site/`. The navigation is defined in `mkdocs.yml`. The render pipeline uses MathJax for equations, the MkDocs Material `social` plugin for OpenGraph preview cards, and the Material theme's built-in sitemap generator for search-engine indexing.
