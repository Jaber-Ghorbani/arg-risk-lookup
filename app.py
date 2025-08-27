
import streamlit as st
import pandas as pd
from rapidfuzz import process, fuzz
from pathlib import Path
import colorsys

st.set_page_config(page_title="ARG/ARM Risk Lookup", layout="wide", page_icon="🧬")

def clamp01(x):
    try:
        x = float(x)
    except:
        return None
    if x < 0:
        x = 0.0
    if x > 1:
        x = 1.0
    return x

def risk_color_hex(val):
    v = clamp01(val)
    if v is None:
        return "#455a64"
    hue = (1.0 - v) * 120.0 / 360.0
    l, s = 0.45, 0.7
    r, g, b = colorsys.hls_to_rgb(hue, l, s)
    return "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))

def risk_badge_html(val):
    v = clamp01(val)
    if v is None:
        return '<span style="display:inline-block;padding:4px 10px;border-radius:999px;background:#455a64;color:#fff;">N/A</span>'
    color = risk_color_hex(v)
    return f'<span style="display:inline-block;padding:4px 10px;border-radius:999px;background:{color};color:white;font-weight:600;">{v:.3f}</span>'

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
    st.caption("Levels first, then scores. Use the calculator to estimate sample-level risk indices.")
    st.write("**Columns in display order:**")
    st.code("\n".join(DISPLAY_COLS))

tab1, tab2, tab3, tab4 = st.tabs(["🔎 Lookup", "📥 Bulk Lookup", "🧮 Risk Index Calculator", "ℹ️ About & Media"])

with tab1:
    st.subheader("Single Gene Lookup")
    c1, c2 = st.columns([3,1])

    gene_options = df["Genes"].astype(str).tolist()
    sel = c1.selectbox("Autocomplete", options=["— Select a gene —"] + gene_options, index=0)
    q_free = c1.text_input("Or type a gene", placeholder="e.g., dfra24").strip()
    fuzzy = c2.checkbox("Fuzzy match", value=True)

    q = None
    if sel and sel != "— Select a gene —":
        q = sel
    elif q_free:
        q = q_free

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
                    badge = risk_badge_html(best["Final_Risk_score"])
                    st.markdown(badge, unsafe_allow_html=True)
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
                    badge = risk_badge_html(row["Final_Risk_score"])
                    st.markdown(badge, unsafe_allow_html=True)
                st.dataframe(hit[DISPLAY_COLS], use_container_width=True, hide_index=True)

with tab2:
    st.subheader("Bulk Lookup")
    st.write("Paste a list of gene names (one per line) and download the results.")
    bulk_text = st.text_area("Genes list", height=160, placeholder="dfra24\nblaTEM\nmecA")
    c1, c2 = st.columns([1,1])
    fuzzy_bulk = c1.checkbox("Fuzzy match", value=True)
    cutoff = c2.slider("Fuzzy cutoff", min_value=50, max_value=95, value=70, step=1)
    if bulk_text:
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
        st.download_button("Download results (CSV)", bulk_df.to_csv(index=False).encode("utf-8"), "bulk_lookup.csv", "text/csv")

with tab3:
    st.subheader("Risk Index Calculator")
    st.markdown("""
    **Equation**:
    \n
    \\[ \\text{Risk Index}_{\\text{sample}} = \\sum_{i=1}^{n} (\\text{Abundance}_i \\times \\text{Risk Score}_i) \\]
    """)
    c1, c2, c3 = st.columns([2,2,1])
    score_cols = [c for c in ["Final_Risk_score","Pathogenic_score","Mobility_score","Clinial_Importance_score","Transmissbilitty_score"] if c in df.columns]
    if not score_cols:
        st.error("No score columns found. Expect one of: Final_Risk_score, Pathogenic_score, Mobility_score, Clinial_Importance_score, Transmissbilitty_score.")
    else:
        chosen_score = c1.selectbox("Risk score column to use", options=score_cols, index=0 if "Final_Risk_score" in score_cols else 0)
        fuzzy_calc = c2.checkbox("Fuzzy match", value=True)
        cutoff_calc = c3.slider("Cutoff", 50, 95, 70, 1)

        st.write("### Input")
        st.write("Paste CSV-like lines with `Gene, Abundance`. Example:")
        st.code("dfra24, 12.5\nmecA, 0.8")
        calc_text = st.text_area("Genes and Abundances", height=140, placeholder="geneA, 10\ngeneB, 3.5")

        uploaded = st.file_uploader("...or upload a CSV with columns: Genes, Abundance", type=["csv"])
        input_rows = []

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
                    input_rows.append({"Genes": str(r[gcol]), "Abundance": float(r[acol])})
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

            if chosen_score not in joined.columns:
                st.error(f"Selected score column '{chosen_score}' is not present in your dataset.")
            else:
                joined["Risk_Score"] = pd.to_numeric(joined[chosen_score], errors="coerce")
                joined["Product"] = joined["Abundance"] * joined["Risk_Score"]
                total = joined["Product"].sum(skipna=True)

                show_cols = ["Genes_q","Genes","Abundance","Risk_Score","Product"] + present_levels + present_scores
                show_cols = [c for c in show_cols if c in joined.columns]
                st.write("### Matches & Contributions")
                st.dataframe(joined[show_cols], use_container_width=True, hide_index=True)

                if chosen_score == "Final_Risk_score":
                    st.markdown("**Sample Risk Index (Σ Abundance × Final_Risk_score):**")
                    st.markdown(risk_badge_html(min(max(total, 0.0), 1.0)), unsafe_allow_html=True)
                else:
                    st.metric("Risk Index", f"{total:.6g}")
                st.download_button("Download matched table (CSV)", joined.to_csv(index=False).encode("utf-8"), "risk_index_breakdown.csv", "text/csv")
        else:
            st.info("Enter some data above to compute the Risk Index.")

with tab4:
    st.subheader("About & Media")
    st.markdown("""
    **ARM context:** This portal supports Antimicrobial Resistance (AMR/ARM) risk exploration by mapping gene-level attributes to a final risk score and enabling sample-level risk index computation.
    """)
    st.write("Upload images (e.g., risk assessment diagram, *Staphylococcus*):")
    imgs = st.file_uploader("Upload images (PNG/JPG/SVG)", type=["png","jpg","jpeg","svg"], accept_multiple_files=True)
    if imgs:
        cols = st.columns(3)
        i = 0
        for im in imgs:
            with cols[i % 3]:
                st.image(im, caption=im.name, use_container_width=True)
            i += 1

    asset_dir = Path("assets")
    if asset_dir.exists():
        files = [p for p in asset_dir.iterdir() if p.suffix.lower() in {".png",".jpg",".jpeg",".svg"}]
        if files:
            st.write("Images in repository assets:")
            cols2 = st.columns(3)
            i = 0
            for p in files:
                with cols2[i % 3]:
                    st.image(str(p), caption=p.name, use_container_width=True)
                i += 1
