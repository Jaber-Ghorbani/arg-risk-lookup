import streamlit as st
import pandas as pd
from rapidfuzz import process, fuzz
from io import StringIO

st.set_page_config(page_title="ARG Risk Lookup", layout="centered", page_icon="🔎")
st.title("ARG Risk Lookup")
st.caption("Search any gene to view Clinical Importance, Transmissibility, Mobility, Pathogenicity score, and Final Risk.")

ALIASES = {
    "genes": "Genes",
    "gene": "Genes",
    "clinical importance": "Clinical_Importance",
    "clinial importance": "Clinical_Importance",
    "clinical_importance": "Clinical_Importance",
    "transmissibility": "Transmissibility",
    "transmissbilitty": "Transmissibility",
    "mobility": "Mobility",
    "pathogenic score": "Pathogenicity_Score",
    "pathogenic_score": "Pathogenicity_Score",
    "pathogenicity score": "Pathogenicity_Score",
    "clinical importance value": "Clinical_Importance_Value",
    "clinical_importance_value": "Clinical_Importance_Value",
    "transmissibility score": "Transmissibility_Score",
    "transmissbilitty score": "Transmissibility_Score",
    "transmissibility_score": "Transmissibility_Score",
    "mobility score": "Mobility_Score",
    "mobility_score": "Mobility_Score",
    "pathogenic score value": "Pathogenicity_Score_Value",
    "pathogenicity score value": "Pathogenicity_Score_Value",
    "pathogenicity_score_value": "Pathogenicity_Score_Value",
    "final risk": "Final_Risk",
    "final_risk": "Final_Risk",
}

CANONICAL_COLS = [
    "Genes",
    "Clinical_Importance",
    "Transmissibility",
    "Mobility",
    "Pathogenicity_Score",
    "Clinical_Importance_Value",
    "Transmissibility_Score",
    "Mobility_Score",
    "Pathogenicity_Score_Value",
    "Final_Risk",
]

def sanitize(name: str) -> str:
    n = (name or "").strip().replace("\xa0"," ").replace("\t"," ").replace("  "," ")
    n = n.replace("-", " ").replace("/", " ").replace("\\", " ")
    n = "_".join([p for p in n.split() if p])
    return n

def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    new_cols = []
    seen = {}
    for c in df.columns:
        key = sanitize(c).lower()
        target = ALIASES.get(key, sanitize(c))
        cnt = seen.get(target, 0)
        if cnt:
            target = f"{target}_{cnt+1}"
        seen[target] = cnt + 1
        new_cols.append(target)
    df = df.copy()
    df.columns = new_cols
    if "Pathogenicity_Score_Value" not in df.columns and "Pathogenicity_Score_2" in df.columns:
        df = df.rename(columns={"Pathogenicity_Score_2":"Pathogenicity_Score_Value"})
    return df

@st.cache_data(show_spinner=False)
def load_data(path: str) -> pd.DataFrame:
    try:
        d = pd.read_csv(path, low_memory=False)
    except Exception:
        try:
            d = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception:
            d = pd.read_csv(path, sep=None, engine="python", low_memory=False)
    d = canonicalize_columns(d)
    if "Genes" not in d.columns:
        raise ValueError("Missing 'Genes' column after normalization.")
    d["gene_key"] = d["Genes"].astype(str).str.strip().str.lower()
    return d.fillna("Unknown")

def pick_columns(d: pd.DataFrame) -> list:
    cols = [c for c in CANONICAL_COLS if c in d.columns]
    extras = [c for c in d.columns if c not in set(cols + ["gene_key"])]
    return cols + extras

def color_badge(text: str) -> str:
    t = str(text).strip().lower()
    if t in {"high","very high","critical"}:
        bg = "#c62828"
    elif t in {"medium","moderate"}:
        bg = "#ef6c00"
    elif t in {"low","marginally important","marginal"}:
        bg = "#2e7d32"
    else:
        bg = "#455a64"
    return f'<span style="display:inline-block;padding:2px 8px;border-radius:999px;background:{bg};color:white;font-size:0.85rem;">{text}</span>'

with st.sidebar:
    st.subheader("Dataset")
    st.write("This app uses **genes_risk.csv** at the app root. Replace it to update data.")
    st.markdown("—")
    st.write("Tips: one row per gene; consistent headers. Misspellings auto-corrected.")
    st.markdown("—")
    st.write("Version: 1.0.0")

df = load_data("genes_risk.csv")
cols_for_display = pick_columns(df)

tab1, tab2 = st.tabs(["🔍 Single search", "📥 Bulk search"])

with tab1:
    left, right = st.columns([3,1])
    q = left.text_input("Gene name", placeholder="e.g., dfra24").strip()
    fuzzy = right.checkbox("Fuzzy match", value=True)
    if q:
        q_key = q.lower()
        if not fuzzy:
            hit = df[df["gene_key"] == q_key]
            if hit.empty:
                st.warning("No exact match found.")
            else:
                row = hit.iloc[0]
                st.subheader(row["Genes"])
                if "Final_Risk" in hit.columns:
                    st.markdown(color_badge(str(row["Final_Risk"])), unsafe_allow_html=True)
                st.dataframe(hit[cols_for_display], use_container_width=True, hide_index=True)
        else:
            choices = df["gene_key"].tolist()
            top = process.extract(q_key, choices, scorer=fuzz.WRatio, limit=5)
            if not top:
                st.warning("No fuzzy match found.")
            else:
                best_key = top[0][0]
                best = df[df["gene_key"] == best_key].iloc[0]
                st.subheader(best["Genes"])
                if "Final_Risk" in df.columns:
                    st.markdown(color_badge(str(best["Final_Risk"])), unsafe_allow_html=True)
                st.dataframe(pd.DataFrame([best])[cols_for_display], use_container_width=True, hide_index=True)
                with st.expander("See similar matches"):
                    sim_rows = []
                    for choice, score, _ in top:
                        r = df[df["gene_key"] == choice].iloc[0]
                        sim_rows.append({"Match": r["Genes"], "Score": score})
                    st.dataframe(pd.DataFrame(sim_rows), use_container_width=True, hide_index=True)

with tab2:
    st.write("Paste a list of genes (one per line).")
    bulk_text = st.text_area("Genes list", height=180, placeholder="dfra24\nblaTEM\nmecA")
    fuzzy_bulk = st.checkbox("Fuzzy match for bulk", value=True, key="bulk")
    cutoff = st.slider("Fuzzy score cutoff", min_value=50, max_value=95, value=70, step=1)
    if bulk_text:
        queries = [x.strip() for x in bulk_text.splitlines() if x.strip()]
        out_rows = []
        choices = df["gene_key"].tolist()
        for q in queries:
            q_key = q.lower()
            if not fuzzy_bulk:
                hit = df[df["gene_key"] == q_key]
                if hit.empty:
                    out_rows.append({"Query": q, "Match": "", "Note": "No exact match"})
                else:
                    row = hit.iloc[0].to_dict()
                    row["Query"] = q
                    row["Match"] = row.get("Genes", "")
                    row["Note"] = "Exact"
                    out_rows.append(row)
            else:
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
        bulk_df = pd.DataFrame(out_rows)
        st.subheader("Bulk results")
        st.dataframe(bulk_df, use_container_width=True, hide_index=True)
        st.download_button("Download results (CSV)", bulk_df.to_csv(index=False).encode("utf-8"), "bulk_lookup.csv", "text/csv")