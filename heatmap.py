from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path("/Users/javierdelgadoluis/Desktop/PROYECTO_FINAL_TFG")

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.colors import LinearSegmentedColormap

# Folder where summaries are stored
RES_FOLDER = Path("/Users/javierdelgadoluis/Desktop/PROYECTO_FINAL_TFG/results_prs")

# ============================================
# Localizar automatically the latest summary
# ============================================
def cargar_ultimo_resumen():
    files = sorted(
        Path(RES_FOLDER).glob("resumen_prs_*.csv"),
        key=lambda p: p.stat().st_mtime,
        reverse=True
    )
    if not files:
        raise FileNotFoundError("No se encontró ningún resumen_prs_*.csv en results_prs")
    resumen_path = files[0]
    print(f"Using summary: {resumen_path}")
    df = pd.read_csv(resumen_path)
    return df

df = cargar_ultimo_resumen()

# Acortar names of PGS
df["PGS"] = df["PGS"].astype(str)
df["VCF"] = df["VCF"].astype(str)
df["Status"] = df["Status"].astype(str).str.strip()

# 1) Quitar the parte redundante of annotation if exists (37 o 38)
df["PGS_short"] = (
    df["PGS"]
    .str.replace(r"_GRCh37.*$", "", regex=True)
    .str.replace(r"(_GRCh38).*$", r"\1", regex=True)
)
df["VCF_short"] = df["VCF"].str.replace(r"_anotado.*$", "", regex=True)

print(f"Loaded  {df.shape[0]} rows")
print(f"Samples VCF unique: {df['VCF_short'].nunique()}")
print(f"PGS unique: {df['PGS'].nunique()}")

# ============================================
# Filter rows valid
# ============================================
df = df[df["Status"] == "OK"].copy()
df = df[df["SNPs_matched"] > 0].copy()
df["Percentil"] = pd.to_numeric(df["Percentil"], errors="coerce")

# ============================================
# Diagnosis: Ver which PGS tienes disponibles
# ============================================
print("\n" + "="*60)
print("PGS DISPONIBLES In TUS Data:")
print("="*60)
pgs_unicos = sorted(df["PGS_short"].unique())
for i, pgs in enumerate(pgs_unicos, 1):
    print(f"{i}. {pgs}")
print("="*60)
print("\n⚠️ COPIA And PEGA estos names in ORDEN_PGS_CUSTOM for ordenarlos\n")

# ============================================
# HEATMAP: Create pivot
# ============================================
heatmap_percentil = df.pivot_table(
    index="VCF_short",
    columns="PGS_short",
    values="Percentil",
    aggfunc="first"
)

print("Before Of Sort:")
print(f"Columns: {heatmap_percentil.columns.tolist()}\n")

# ============================================
# ORDEN PERSONALIZADO - MODIFICA This List
# ============================================
ORDEN_PGS_CUSTOM = ['Benafif_S', 'Benafif_S_GRCh38', 
                    'Schumacher_FR', 'Schumacher_FR_GRCh38',
                    'Kim_ES', 'Kim_ES_GRCh38', 
                    'Xin_J', 'Xin_J_GRCh38',
                    'Shi_Z', 'Shi_Z_GRCh38',
                    'Karunamuni_RA', 'Karunamuni_RA_GRCh38',
                    'Sipeky_C', 'Sipeky_C_GRCh38',
                    'Jia_G', 'Jia_G_GRCh38',
                    'Graff_RE', 'Graff_RE_GRCh38', 
                    'Conti_DV_Darst_BF', 'Conti_DV_Darst_BF_GRCh38',
                    'Wang_A', 'Wang_A_GRCh38'
                   ]

# ============================================
# ORDEN POR RIESGO MEDIO (respetando tu lista)
# ============================================
# Media de percentil por PGS (columna)
media_por_pgs = heatmap_percentil.mean(axis=0)  # Serie: PGS_short -> media

# Dividir en PGS que están en tu lista y el resto
presentes = [pgs for pgs in ORDEN_PGS_CUSTOM if pgs in heatmap_percentil.columns]
resto = [pgs for pgs in heatmap_percentil.columns if pgs not in ORDEN_PGS_CUSTOM]

# Dentro de cada grupo, ordenar por percentil medio (de menor a mayor riesgo)
presentes_ordenados = sorted(presentes, key=lambda x: media_por_pgs[x])
resto_ordenados = sorted(resto, key=lambda x: media_por_pgs[x])

orden_final = presentes_ordenados + resto_ordenados

heatmap_percentil = heatmap_percentil[orden_final]

print("After Of Sort (by mean risk):")
print(f"Columns: {heatmap_percentil.columns.tolist()}\n")

# ============================================
# Colormap
# ============================================
cmap_custom = LinearSegmentedColormap.from_list(
    "percentil_riesgo",
    [
        (0/100,   "#d6eaf8"),  # 0   - azul muy pálido (bajo riesgo, discreto)
        (40/100,  "#aed6f1"),  # 40  - azul suave
        (40/100,  "#f9e79f"),  # 40  - amarillo pálido (riesgo moderado)
        (80/100,  "#f39c12"),  # 80  - naranja
        (80/100,  "#e74c3c"),  # 80  - rojo (riesgo elevado)
        (90/100,  "#c0392b"),  # 90  - rojo medio
        (100/100, "#7b241c"),  # 100 - rojo muy oscuro
    ],
    N=256
)

# ============================================
# HEATMAP GLOBAL
# ============================================
plt.figure(figsize=(16, 10))

ax = sns.heatmap(
    heatmap_percentil,
    cmap=cmap_custom,
    vmin=0,
    vmax=100,
    cbar_kws={
        "label": "Percentil de Riesgo",
    },
    linewidths=0.5,
    linecolor="white",
    annot=True,
    fmt=".1f",
    annot_kws={"size": 9, "weight": "bold"},
)

for text in ax.texts:
    val = float(text.get_text())
    text.set_color("white" if val > 50 else "#2c3e50")
    text.set_fontsize(9)
    text.set_fontweight("bold")

colorbar = ax.collections[0].colorbar
colorbar.set_label("Percentil de Riesgo", fontsize=13, weight="bold")
colorbar.set_ticks([20, 60, 85, 95])
colorbar.set_ticklabels([
    "0-40\n(Bajo)",
    "40-80\n(Moderado)",
    "80-90\n(Elevado)",
    "90-100\n(Alto)"
])

plt.title(
    "Riesgo de Cáncer de Próstata por Muestra y PGS\n(Percentiles de Riesgo Poblacional)",
    fontsize=16,
    pad=20,
    weight="bold",
)
plt.xlabel("Polygenic Score (Estudio)", fontsize=12, weight="bold")
plt.ylabel("Muestra VCF", fontsize=12, weight="bold")
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(rotation=0, fontsize=10)
plt.tight_layout()

out_svg = Path(RES_FOLDER) / "heatmap_percentiles_gradiente.svg"
plt.savefig(out_svg, format="svg", bbox_inches="tight")
print(f"\n✓ Heatmap saved: {out_svg}")
plt.show()

# ============================================
# HEATMAP INDIVIDUAL POR MUESTRA (un SVG por VCF)
# ============================================
for sample_id in heatmap_percentil.index:
    data_ind = heatmap_percentil.loc[[sample_id]]  # 1 fila, mantiene 2D

    # Más alto para que haya espacio vertical
    plt.figure(figsize=(16, 4))  # antes 16 x 2.5

    ax = sns.heatmap(
        data_ind,
        cmap=cmap_custom,
        vmin=0,
        vmax=100,
        cbar_kws={
            "label": "Percentil de Riesgo",
        },
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".1f",
        annot_kws={"size": 8, "weight": "bold"},  # un pelín más pequeño
    )

    for text in ax.texts:
        val = float(text.get_text())
        text.set_color("white" if val > 50 else "#2c3e50")
        text.set_fontsize(8)          # bajar un poco fuente
        text.set_fontweight("bold")

    colorbar = ax.collections[0].colorbar
    # Aumentar tamaño del colorbar y espaciado entre etiquetas
    colorbar.ax.tick_params(labelsize=9, pad=6)  # fuente un poco mayor y separada
    colorbar.set_label("Percentil de Riesgo", fontsize=11, weight="bold", labelpad=8)
    colorbar.set_ticks([20, 60, 85, 95])
    colorbar.set_ticklabels([
        "0-40 (Bajo)",
        "40-80 (Moderado)",
        "80-90 (Elevado)",
        "90-100 (Alto)"
    ])


    plt.title(
        f"Riesgo genético de cáncer de próstata\nMuestra: {sample_id}",
        fontsize=14,
        pad=18,
        weight="bold",
    )
    plt.xlabel("Polygenic Score (Estudio)", fontsize=11, weight="bold")
    plt.ylabel("Muestra VCF", fontsize=11, weight="bold")

    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(rotation=0, fontsize=9)

    plt.tight_layout()

    out_ind_svg = Path(RES_FOLDER) / f"heatmap_percentiles_{sample_id}.svg"
    plt.savefig(out_ind_svg, format="svg", bbox_inches="tight")
    print(f"✓ Heatmap individual guardado: {out_ind_svg}")
    plt.close()


# ============================================
# Summary Statistical
# ============================================
print("\n" + "=" * 60)
print("Summary Of RESULTADOS")
print("=" * 60)

def categorizar_percentil(p):
    if pd.isna(p):
        return np.nan
    elif p < 20:
        return "Bajo (0-20)"
    elif p < 40:
        return "Bajo-Medio (20-40)"
    elif p < 60:
        return "Medio (40-60)"
    elif p < 80:
        return "Medio-Alto (60-80)"
    else:
        return "Alto (80-100)"

df["Categoria_Riesgo"] = df["Percentil"].apply(categorizar_percentil)

print("\nPercentil promedio by sample:")
percentil_promedio = df.groupby("VCF_short")["Percentil"].mean().sort_values(ascending=False)
for vcf, perc in percentil_promedio.items():
    print(f" {vcf}: {perc:.2f}")

print("\nZ-Score promedio by sample:")
df["Z_Score"] = pd.to_numeric(df["Z_Score"], errors="coerce")
zscore_promedio = df.groupby("VCF_short")["Z_Score"].mean().sort_values(ascending=False)
for vcf, score in zscore_promedio.items():
    print(f" {vcf}: {score:.3f}")

print("\nDistribution of risk categories by sample:")
categoria_counts = df.groupby(["VCF_short", "Categoria_Riesgo"]).size().unstack(fill_value=0)
print(categoria_counts)

print("\nPercentil promedio by PGS (in orden of the heatmap):")
for pgs in orden_final:
    perc = df[df["PGS_short"] == pgs]["Percentil"].mean()
    if not pd.isna(perc):
        print(f" {pgs}: {perc:.2f}")


print("\n✓ Analysis completed exitosamente!")


# ============================================
# RESUMEN POR INDIVIDUO (para diapositiva)
# ============================================
df_resumen_ind = df.groupby("VCF_short").agg({
    "Percentil": ["mean", "std", "min", "max", lambda x: (x > 80).sum(), lambda x: (x > 90).sum()]
}).round(2)

df_resumen_ind.columns = ["Media", "Desv.Est", "Mín", "Máx", "PGS>80", "PGS>90"]
df_resumen_ind = df_resumen_ind.sort_values("Media", ascending=False)

print("\n" + "="*80)
print("RESUMEN POR INDIVIDUO (ordenado por riesgo medio)")
print("="*80)
print(df_resumen_ind)
print(f"\nGuardando en {RES_FOLDER}/resumen_por_individuo.csv")
df_resumen_ind.to_csv(Path(RES_FOLDER)/"resumen_por_individuo.csv")


# ============================================
# RESUMEN POR ESTUDIO + CORRELACIONES
# ============================================
df_resumen_pgs = df.groupby("PGS_short").agg({
    "Percentil": ["count", "mean", "std", lambda x: (x > 80).sum()]
}).round(2)

df_resumen_pgs.columns = ["N_Ind", "Media", "Desv.Est", "%>80"]
df_resumen_pgs = df_resumen_pgs.sort_values("Media", ascending=False)

print("\n" + "="*80)
print("RESUMEN POR ESTUDIO PGS")
print("="*80)
print(df_resumen_pgs)
print(f"\nGuardando en {RES_FOLDER}/resumen_por_pgs.csv")
df_resumen_pgs.to_csv(Path(RES_FOLDER)/"resumen_por_pgs.csv")

# MATRIZ CORRELACIONES entre PGS
pivot_correl = df.pivot_table(index="VCF_short", columns="PGS_short", values="Percentil")
correl_matrix = pivot_correl.corr().round(2)

print("\n" + "="*80)
print("CORRELACIONES ENTRE PGS (Pearson)")
print("="*80)
print(correl_matrix)


plt.figure(figsize=(12, 10))
sns.heatmap(correl_matrix, annot=True, cmap="RdBu", center=0, vmin=-1, vmax=1,
            cbar_kws={"label": "Correlación Pearson"})
plt.title("Correlación entre PGS\n(Azul=positiva, Rojo=negativa)", fontsize=14, weight="bold")
plt.tight_layout()
plt.savefig(Path(RES_FOLDER)/"correlaciones_pgs.svg", format="svg", bbox_inches="tight")
plt.show()

# ============================================
# DISCREPANCIAS: Casos extremos (mismo ind, PGS muy distintos)
# ============================================
# Encontrar individuos con mayor rango entre PGS
df["rango_percentil"] = df.groupby("VCF_short")["Percentil"].transform(lambda x: x.max() - x.min())
casos_extremos = df.loc[df["rango_percentil"] > 50].groupby("VCF_short").agg({
    "Percentil": ["min", "max"]  # ← QUITADO "PGS_short"
}).round(2)

print("\n" + "="*80)
print("INDIVIDUOS CON MAYOR DISCREPANCIA (>50 puntos entre PGS)")
print("="*80)
print(casos_extremos)

# Para el más extremo, mostrar qué PGS dan bajo/alto
if len(casos_extremos) > 0:
    individuo_extremo = casos_extremos.index[0]
    print(f"\n📊 CASO ESTUDIO: {individuo_extremo}")
    discrepancias = df[(df["VCF_short"] == individuo_extremo)].nlargest(3, "Percentil")[["PGS_short", "Percentil"]]
    discrepancias = pd.concat([discrepancias, 
                              df[(df["VCF_short"] == individuo_extremo)].nsmallest(3, "Percentil")[["PGS_short", "Percentil"]]])
    print(discrepancias.round(1).drop_duplicates())
