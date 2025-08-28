
import streamlit as st
import pandas as pd
from rapidfuzz import process, fuzz
from pathlib import Path

st.set_page_config(page_title="ARG/ARM Risk Lookup", layout="wide", page_icon="🧬")

# ----------------- Color utilities -----------------
def _interp(c1, c2, t):
    c1 = c1.lstrip("#"); c2 = c2.lstrip("#")
    r1,g1,b1 = int(c1[0:2],16), int(c1[2:4],16), int(c1[4:6],16)
    r2,g2,b2 = int(c2[0:2],16), int(c2[2:4],16), int(c2[4:6],16)
    r = int(r1 + (r2 - r1) * t)
    g = int(g1 + (g2 - g1) * t)
    b = int(b1 + (b2 - b1) * t)
    return f"#{r:02x}{g:02x}{b:02x}"

def _to_percent(x):
    try:
        v = float(x)
    except:
        return None
    if v < 0: v = 0.0
    if v <= 1.0: v = v * 100.0
    elif v > 100.0: v = 100.0
    return v

def risk_color_hex_percent(val):
    p = _to_percent(val)
    if p is None:
        return "#455a64"
    if p <= 50.0:
        t = p / 50.0
        return _interp("#2e7d32", "#ffffff", t)
    else:
        t = (p - 50.0) / 50.0
        return _interp("#ffffff", "#c62828", t)

def risk_badge_html(val):
    p = _to_percent(val)
    if p is None:
        return '<span style="display:inline-block;padding:4px 10px;border-radius:999px;background:#455a64;color:#fff;">N/A</span>'
    color = risk_color_hex_percent(p)
    label = f"{p:.1f}"
    return f'<span style="display:inline-block;padding:4px 10px;border-radius:999px;background:{color};color:#000;font-weight:700;border:1px solid #ddd;">{label}</span>'

# ----------------- Data -----------------
@st.cache_data(show_spinner=False)
def load_data(path: str) -> pd.DataFrame:
    try:
        d = pd.read_csv(path, low_memory=False)
    except Exception:
        try:
            d = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception:
            d = pd.read_csv(path, sep=None, engine="python", low_memory=False)
    if "Genes" not in d.columns:
        raise ValueError("Your file must include a 'Genes' column.")
    d["gene_key"] = d["Genes"].astype(str).str.strip().str.lower()
    return d

df = load_data("genes_risk.csv")

LEVEL_ORDER = ["Clinical_Importance_level","Transmissibility_level","Mobility_level","Pathogenic_level"]
SCORE_ORDER = ["Clinial_Importance_score","Transmissbilitty_score","Mobility_score","Pathogenic_score","Final_Risk_score"]
present_levels = [c for c in LEVEL_ORDER if c in df.columns]
present_scores = [c for c in SCORE_ORDER if c in df.columns]
base_cols = ["Genes"] + present_levels + present_scores
extra_cols = [c for c in df.columns if c not in set(base_cols + ["gene_key"])]
DISPLAY_COLS = base_cols + extra_cols

with st.sidebar:
    st.markdown("# ARG / ARM Risk Portal")
    st.caption("Levels first → Scores. Calculator computes Σ(Abundance × Risk score).")
    st.write("**Columns in display order:**")
    st.code("\\n".join(DISPLAY_COLS))
    st.markdown("---")
    st.markdown("**Color scale (Final_Risk_score):** 0=green → 50=white → 100=red")
    st.markdown("---")
    st.info("**Data note**: Not all genes have levels defined for every criterion. To maximize coverage and enable ranking, **missing values are recorded as 'Not Defined' and assigned a numeric value of 0** for scoring and risk calculations.")

tab1, tab2, tab3, tab4 = st.tabs(["🔎 Lookup", "📥 Bulk Lookup", "🧮 Risk Index Calculator", "ℹ️ About"])

# ----------------- Lookup -----------------
with tab1:
    st.subheader("Single Gene Lookup")
    c1, c2, c3 = st.columns([3,1,1])

    gene_options = df["Genes"].astype(str).tolist()
    sel = c1.selectbox("Autocomplete", options=["— Select a gene —"] + gene_options, index=0, key="sel_single")
    q_free = c1.text_input("Or type a gene", placeholder="e.g., dfra24", key="free_single").strip()
    fuzzy = c2.checkbox("Fuzzy match", value=True, key="fuzzy_single")
    do_search = c3.button("Search", key="search_btn")

    q = None
    if do_search:
        q = sel if sel and sel != "— Select a gene —" else (q_free if q_free else None)
    elif sel and sel != "— Select a gene —":
        q = sel

    if q:
        q_key = q.lower()
        if fuzzy and (not sel or sel == "— Select a gene —"):
            choices = df["gene_key"].tolist()
            top = process.extract(q_key, choices, scorer=fuzz.WRatio, limit=5)
            if not top:
                st.warning("No fuzzy match found.")
            else:
                best_key = top[0][0]
                best = df[df["gene_key"] == best_key].iloc[0]
                st.success(f"Best match: **{best['Genes']}**")
                if "Final_Risk_score" in df.columns:
                    st.markdown(risk_badge_html(best["Final_Risk_score"]), unsafe_allow_html=True)
                st.dataframe(pd.DataFrame([best])[DISPLAY_COLS], use_container_width=True, hide_index=True)
                with st.expander("Similar matches"):
                    sim_rows = []
                    for choice, score, _ in top:
                        r = df[df["gene_key"] == choice].iloc[0]
                        sim_rows.append({"Match": r["Genes"], "Score": score})
                    st.dataframe(pd.DataFrame(sim_rows), use_container_width=True, hide_index=True)
        else:
            hit = df[df["gene_key"] == q_key]
            if hit.empty:
                st.warning("No exact match found.")
            else:
                row = hit.iloc[0]
                if "Final_Risk_score" in hit.columns:
                    st.markdown(risk_badge_html(row["Final_Risk_score"]), unsafe_allow_html=True)
                st.dataframe(hit[DISPLAY_COLS], use_container_width=True, hide_index=True)

# ----------------- Bulk Lookup (with button) -----------------
with tab2:
    st.subheader("Bulk Lookup")
    st.write("Paste a list of gene names (one per line) and press **Search** to run.")
    bulk_text = st.text_area("Genes list", height=160, placeholder="dfra24\\nblaTEM\\nmecA", key="bulk_text")
    c1, c2, c3 = st.columns([1,1,1])
    fuzzy_bulk = c1.checkbox("Fuzzy match", value=True, key="fuzzy_bulk")
    cutoff = c2.slider("Fuzzy cutoff", min_value=50, max_value=95, value=70, step=1, key="cutoff_bulk")
    do_bulk = c3.button("Search", key="bulk_search_btn")

    if bulk_text and do_bulk:
        queries = [x.strip() for x in bulk_text.splitlines() if x.strip()]
        out_rows = []
        choices = df["gene_key"].tolist()
        for q in queries:
            q_key = q.lower()
            if fuzzy_bulk:
                match = process.extractOne(q_key, choices, scorer=fuzz.WRatio, score_cutoff=cutoff)
                if not match:
                    out_rows.append({"Query": q, "Match": "", "Note": f"No fuzzy match ≥{cutoff}"})
                else:
                    best_key = match[0]
                    row = df[df["gene_key"] == best_key].iloc[0].to_dict()
                    row["Query"] = q
                    row["Match"] = row.get("Genes", "")
                    row["Note"] = "Fuzzy"
                    out_rows.append(row)
            else:
                hit = df[df["gene_key"] == q_key]
                if hit.empty:
                    out_rows.append({"Query": q, "Match": "", "Note": "No exact match"})
                else:
                    row = hit.iloc[0].to_dict()
                    row["Query"] = q
                    row["Match"] = row.get("Genes", "")
                    row["Note"] = "Exact"
                    out_rows.append(row)
        bulk_df = pd.DataFrame(out_rows)
        st.dataframe(bulk_df, use_container_width=True, hide_index=True)
        st.download_button("Download results (CSV)", bulk_df.to_csv(index=False).encode("utf-8"), "bulk_lookup.csv", "text/csv", key="bulk_download")

# ----------------- Risk Index Calculator (with button) -----------------
with tab3:
    st.subheader("Risk Index Calculator")
    eq_img = Path("assets") / "risk_index_equation.png"
    if eq_img.exists():
        st.image(str(eq_img), caption="Risk Index formula", use_container_width=False)
    st.latex(r"\text{Risk Index}_{\text{sample}} = \sum_{i=1}^{n} \left(\text{Abundance}_i \times \text{Risk score}_i\right)")
    st.markdown("Paste or upload genes with abundances; press **Compute** to calculate.")

    c1, c2, c3 = st.columns([2,2,1])
    score_cols = [c for c in ["Final_Risk_score","Pathogenic_score","Mobility_score","Clinial_Importance_score","Transmissbilitty_score"] if c in df.columns]
    if not score_cols:
        st.error("No score columns found. Expect one of: Final_Risk_score, Pathogenic_score, Mobility_score, Clinial_Importance_score, Transmissbilitty_score.")
    else:
        default_index = score_cols.index("Final_Risk_score") if "Final_Risk_score" in score_cols else 0
        chosen_score = c1.selectbox("Risk score column to use", options=score_cols, index=default_index, key="score_choice_calc")
        fuzzy_calc = c2.checkbox("Fuzzy match", value=True, key="fuzzy_calc")
        cutoff_calc = c3.slider("Cutoff", 50, 95, 70, 1, key="cutoff_calc")

        st.write("### Input")
        st.code("dfra24, 12.5\nmecA, 0.8")
        calc_text = st.text_area("Genes and Abundances", height=140, placeholder="geneA, 10\ngeneB, 3.5", key="calc_text")

        uploaded = st.file_uploader("...or upload a CSV with columns: Genes, Abundance", type=["csv"], key="calc_upload")
        compute = st.button("Compute", key="compute_btn")
        input_rows = []

        if compute:
            if calc_text:
                for line in calc_text.splitlines():
                    if "," in line:
                        g, a = line.split(",", 1)
                        g = g.strip()
                        try:
                            a = float(a)
                        except:
                            a = None
                        input_rows.append({"Genes": g, "Abundance": a})
            elif uploaded is not None:
                tmp = pd.read_csv(uploaded)
                colmap = {c.lower(): c for c in tmp.columns}
                gcol = colmap.get("genes") or colmap.get("gene")
                acol = colmap.get("abundance")
                if gcol and acol:
                    for _, r in tmp.iterrows():
                        try:
                            aval = float(r[acol])
                        except:
                            aval = None
                        input_rows.append({"Genes": str(r[gcol]), "Abundance": aval})
                else:
                    st.error("Uploaded CSV must contain 'Genes' and 'Abundance' columns.")

            if input_rows:
                inp = pd.DataFrame(input_rows)
                inp["query_key"] = inp["Genes"].str.strip().str.lower()
                if fuzzy_calc:
                    choices = df["gene_key"].tolist()
                    matches = []
                    for q in inp["query_key"]:
                        match = process.extractOne(q, choices, scorer=fuzz.WRatio, score_cutoff=cutoff_calc)
                        matches.append(match[0] if match else None)
                    inp["match_key"] = matches
                    joined = inp.merge(df, left_on="match_key", right_on="gene_key", how="left", suffixes=("_q",""))
                else:
                    joined = inp.merge(df, left_on="query_key", right_on="gene_key", how="left", suffixes=("_q",""))

                # Use chosen score directly; convert to numeric; replace NaN with 0 (for 'Not Defined')
                joined["Risk_Score"] = pd.to_numeric(joined[chosen_score], errors="coerce").fillna(0.0)
                joined["Product"] = joined["Abundance"] * joined["Risk_Score"]
                total = joined["Product"].sum(skipna=True)

                show_cols = ["Genes_q","Genes","Abundance","Risk_Score","Product"] + present_levels + present_scores
                show_cols = [c for c in show_cols if c in joined.columns]
                st.write("### Matches & Contributions")
                st.dataframe(joined[show_cols], use_container_width=True, hide_index=True)

                st.metric("Risk Index (sample)", f"{total:.6g}")
                st.download_button("Download matched table (CSV)", joined.to_csv(index=False).encode("utf-8"), "risk_index_breakdown.csv", "text/csv", key="calc_download")
            else:
                st.info("Enter some data above, then press **Compute**.")

# ----------------- About -----------------
with tab4:
    st.subheader("About")
    st.markdown("""
This portal supports Antimicrobial Resistance (AMR) risk exploration by mapping gene-level attributes
(Clinical Importance, Transmissibility, Mobility, Pathogenic level/scores) to a final risk score, and providing a calculator for sample-level risk indices based on abundance-weighted risk scores.
""")
    st.markdown("For questions or collaboration, visit our research group: **[Food Safety Risk Analysis (FSRA) Lab, UNL](https://fsra.unl.edu/)**.")
