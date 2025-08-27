# ARG/ARM Risk Lookup — Pro v3

New:
- Autocomplete select for single-gene lookup (type to filter).
- Continuous gradient badge for `Final_Risk_score` (0=green → 1=red).
- Keeps display order: all *_level columns first, then *_score columns, then extras.
- Risk Index Calculator retained with CSV upload and fuzzy matching.

Deploy:
1. Put `app.py`, `requirements.txt`, `genes_risk.csv` and optional `assets/` in a GitHub repo root.
2. Streamlit Community Cloud → New app → set `app.py`.
3. Deploy and share the URL.