"""
================================================================
Day 16 — AMR Co-resistance Network Analysis (REAL DATA)
Author  : Subhadip Jana
Dataset : example_isolates — AMR R package
          2,000 clinical isolates × 40 antibiotics (R/S/I)

What is Co-resistance?
  Two antibiotics are co-resistant when an isolate is resistant
  to BOTH simultaneously. High co-resistance = cross-resistance
  mechanisms (shared genes, efflux pumps, mobile elements).

Network Definition:
  Nodes  : Antibiotics (34 with ≥80 R isolates)
  Edges  : Co-resistance — φ (phi) coefficient between R/S pairs
  Weight : Strength of co-resistance association
  Colour : Antibiotic class

Key Metrics:
  • φ coefficient (Matthews Correlation) — co-resistance strength
  • Odds Ratio — clinical co-resistance risk
  • Network centrality — which ABs are "hubs" of resistance
  • Community detection — resistance clusters (gene families)
  • Clique detection — fully-connected resistance blocks

Methods:
  • Phi coefficient matrix (binary R/S pairwise)
  • Chi-square significance + BH FDR correction
  • NetworkX graph with weighted edges
  • Louvain-style community detection (greedy modularity)
  • Betweenness, degree, eigenvector centrality
================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
import networkx as nx
from scipy.stats import chi2_contingency, fisher_exact
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# SECTION 1: LOAD & PREPARE
# ─────────────────────────────────────────────────────────────

print("🔬 Loading example_isolates dataset...")
df = pd.read_csv("data/isolates.csv")
META = ["date","patient","age","gender","ward","mo"]
ALL_AB = [c for c in df.columns if c not in META]

# Select antibiotics with ≥80 R isolates for meaningful co-resistance
FOCUS_ABS = [ab for ab in ALL_AB if (df[ab]=="R").sum() >= 80]
print(f"✅ {len(df)} isolates | {len(FOCUS_ABS)} antibiotics for network")

# Antibiotic metadata
AB_CLASS = {
    "PEN":"Penicillin","OXA":"Penicillin","FLC":"Penicillin",
    "AMX":"Penicillin","AMC":"Penicillin","AMP":"Penicillin",
    "TZP":"Cephalosporin","CZO":"Cephalosporin","FEP":"Cephalosporin",
    "CXM":"Cephalosporin","FOX":"Cephalosporin","CTX":"Cephalosporin",
    "CAZ":"Cephalosporin","CRO":"Cephalosporin",
    "GEN":"Aminoglycoside","TOB":"Aminoglycoside",
    "AMK":"Aminoglycoside","KAN":"Aminoglycoside",
    "TMP":"Sulfonamide","SXT":"Sulfonamide",
    "CIP":"Fluoroquinolone","MFX":"Fluoroquinolone",
    "VAN":"Glycopeptide","TEC":"Glycopeptide",
    "ERY":"Macrolide","CLI":"Macrolide","AZM":"Macrolide",
    "IPM":"Carbapenem","MEM":"Carbapenem",
    "NIT":"Other","FOS":"Other","LNZ":"Other","TCY":"Tetracycline",
    "TGC":"Tetracycline","DOX":"Tetracycline",
    "MTR":"Other","CHL":"Other","COL":"Polymyxin","MUP":"Other","RIF":"Other",
}

CLASS_COLORS = {
    "Penicillin":"#E74C3C","Cephalosporin":"#E67E22",
    "Aminoglycoside":"#F1C40F","Sulfonamide":"#2ECC71",
    "Fluoroquinolone":"#1ABC9C","Glycopeptide":"#3498DB",
    "Macrolide":"#9B59B6","Carbapenem":"#E91E63",
    "Tetracycline":"#795548","Polymyxin":"#607D8B","Other":"#BDC3C7",
}

# ─────────────────────────────────────────────────────────────
# SECTION 2: COMPUTE PHI COEFFICIENT MATRIX
# ─────────────────────────────────────────────────────────────

print("\n📐 Computing phi (φ) co-resistance matrix...")

def phi_coefficient(x, y):
    """Matthews correlation coefficient for binary co-resistance."""
    valid = x.notna() & y.notna()
    xv = (x[valid] == "R").astype(int)
    yv = (y[valid] == "R").astype(int)
    n  = len(xv)
    if n < 20: return 0, 1.0
    # 2×2 contingency table
    a = ((xv==1) & (yv==1)).sum()   # both R
    b = ((xv==1) & (yv==0)).sum()   # x R, y S
    c = ((xv==0) & (yv==1)).sum()   # x S, y R
    d = ((xv==0) & (yv==0)).sum()   # both S
    denom = np.sqrt((a+b)*(c+d)*(a+c)*(b+d))
    phi   = (a*d - b*c) / denom if denom > 0 else 0
    # Chi-square p-value
    ct = np.array([[a,b],[c,d]])
    if ct.min() >= 5:
        _, p, _, _ = chi2_contingency(ct)
    else:
        _, p = fisher_exact(ct)
    return round(phi, 4), round(p, 6)

n_ab    = len(FOCUS_ABS)
phi_mat = np.zeros((n_ab, n_ab))
p_mat   = np.ones( (n_ab, n_ab))
or_mat  = np.ones( (n_ab, n_ab))

for i, ab1 in enumerate(FOCUS_ABS):
    for j, ab2 in enumerate(FOCUS_ABS):
        if i >= j:
            continue
        phi, p = phi_coefficient(df[ab1], df[ab2])
        phi_mat[i,j] = phi_mat[j,i] = phi
        p_mat[i,j]   = p_mat[j,i]   = p

# BH FDR correction
from itertools import combinations
pairs = list(combinations(range(n_ab), 2))
pvals = [p_mat[i,j] for i,j in pairs]
n_tests = len(pvals)
sorted_idx = np.argsort(pvals)
fdr_thresh = 0.05
qvals = np.ones(n_tests)
for rank, idx in enumerate(sorted_idx):
    qvals[idx] = min(pvals[idx] * n_tests / (rank+1), 1.0)

sig_mask = np.zeros((n_ab, n_ab), dtype=bool)
for k, (i,j) in enumerate(pairs):
    if qvals[k] < fdr_thresh and abs(phi_mat[i,j]) >= 0.2:
        sig_mask[i,j] = sig_mask[j,i] = True

n_sig = sig_mask.sum() // 2
print(f"   Total pairs tested : {n_tests}")
print(f"   Significant (FDR<0.05, |φ|≥0.2): {n_sig}")
print(f"   Positive co-resistance: {((phi_mat > 0) & sig_mask).sum()//2}")
print(f"   Negative (antagonistic): {((phi_mat < 0) & sig_mask).sum()//2}")

# Save matrix
phi_df = pd.DataFrame(phi_mat, index=FOCUS_ABS, columns=FOCUS_ABS)
phi_df.to_csv("outputs/phi_coresistance_matrix.csv")

# ─────────────────────────────────────────────────────────────
# SECTION 3: BUILD NETWORK
# ─────────────────────────────────────────────────────────────

print("\n🕸️  Building co-resistance network...")

G = nx.Graph()

# Add nodes
for ab in FOCUS_ABS:
    n_R = int((df[ab]=="R").sum())
    G.add_node(ab,
               ab_class = AB_CLASS.get(ab,"Other"),
               color    = CLASS_COLORS.get(AB_CLASS.get(ab,"Other"),"#BDC3C7"),
               n_R      = n_R,
               pct_R    = round(n_R/df[ab].notna().sum()*100, 1))

# Add edges (significant co-resistance only)
for i, ab1 in enumerate(FOCUS_ABS):
    for j, ab2 in enumerate(FOCUS_ABS):
        if j <= i: continue
        if sig_mask[i,j]:
            phi = phi_mat[i,j]
            G.add_edge(ab1, ab2,
                       weight   = abs(phi),
                       phi      = phi,
                       edge_type= "positive" if phi > 0 else "negative")

print(f"   Nodes   : {G.number_of_nodes()}")
print(f"   Edges   : {G.number_of_edges()}")
print(f"   Density : {nx.density(G):.4f}")
pos_edges = [(u,v) for u,v,d in G.edges(data=True) if d["phi"]>0]
neg_edges = [(u,v) for u,v,d in G.edges(data=True) if d["phi"]<0]
print(f"   Positive edges (co-resistance): {len(pos_edges)}")
print(f"   Negative edges (antagonistic) : {len(neg_edges)}")

# ─────────────────────────────────────────────────────────────
# SECTION 4: NETWORK METRICS
# ─────────────────────────────────────────────────────────────

print("\n📊 Computing network centrality metrics...")

degree_cent    = nx.degree_centrality(G)
between_cent   = nx.betweenness_centrality(G, weight="weight", normalized=True)
eigenvec_cent  = nx.eigenvector_centrality(G, weight="weight", max_iter=500)
clustering     = nx.clustering(G, weight="weight")

centrality_df = pd.DataFrame({
    "Antibiotic"  : list(G.nodes()),
    "Class"       : [G.nodes[n]["ab_class"] for n in G.nodes()],
    "Degree"      : [G.degree(n) for n in G.nodes()],
    "Degree_cent" : [round(degree_cent[n],4) for n in G.nodes()],
    "Betweenness" : [round(between_cent[n],4) for n in G.nodes()],
    "Eigenvector" : [round(eigenvec_cent[n],4) for n in G.nodes()],
    "Clustering"  : [round(clustering[n],4) for n in G.nodes()],
    "Pct_R"       : [G.nodes[n]["pct_R"] for n in G.nodes()],
}).sort_values("Degree", ascending=False)
centrality_df.to_csv("outputs/network_centrality.csv", index=False)

print("\n   Top 10 hub antibiotics (by degree):")
print(f"   {'AB':5s} {'Class':15s} {'Degree':>7s} {'Betweenness':>12s} {'%R':>7s}")
for _, row in centrality_df.head(10).iterrows():
    print(f"   {row['Antibiotic']:5s} {row['Class']:15s} {int(row['Degree']):>7d} "
          f"{row['Betweenness']:>12.4f} {row['Pct_R']:>6.1f}%")

# ─────────────────────────────────────────────────────────────
# SECTION 5: COMMUNITY DETECTION
# ─────────────────────────────────────────────────────────────

print("\n🔍 Detecting resistance communities...")
communities = nx.community.greedy_modularity_communities(G, weight="weight")
modularity  = nx.community.modularity(G, communities, weight="weight")

community_map = {}
for i, comm in enumerate(communities):
    for node in comm:
        community_map[node] = i

print(f"   Communities detected : {len(communities)}")
print(f"   Modularity           : {modularity:.4f}")
for i, comm in enumerate(communities):
    classes = [G.nodes[n]["ab_class"] for n in comm]
    print(f"   Community {i+1} ({len(comm)} nodes): "
          f"{', '.join(sorted(comm)[:6])}{'...' if len(comm)>6 else ''}")

# ─────────────────────────────────────────────────────────────
# SECTION 6: TOP CO-RESISTANCE PAIRS
# ─────────────────────────────────────────────────────────────

print("\n💊 Top co-resistance pairs...")
edge_data = []
for u, v, d in G.edges(data=True):
    if d["phi"] > 0:
        # Co-resistance rate
        valid = df[u].notna() & df[v].notna()
        both_R = ((df.loc[valid,u]=="R") & (df.loc[valid,v]=="R")).sum()
        n_valid = valid.sum()
        edge_data.append({
            "AB1":u,"AB2":v,"phi":d["phi"],
            "class1":AB_CLASS.get(u,"Other"),
            "class2":AB_CLASS.get(v,"Other"),
            "coresist_pct":round(both_R/n_valid*100,1),
            "n_both_R":int(both_R),
        })
edge_df = pd.DataFrame(edge_data).sort_values("phi", ascending=False)
edge_df.to_csv("outputs/coresistance_pairs.csv", index=False)

print(f"\n   {'AB1':5s} {'AB2':5s} {'φ':>8s} {'Co-R%':>8s} {'Same class':>12s}")
for _, row in edge_df.head(12).iterrows():
    same = "✅" if row["class1"] == row["class2"] else "❌"
    print(f"   {row['AB1']:5s} {row['AB2']:5s} {row['phi']:>8.4f} "
          f"{row['coresist_pct']:>7.1f}% {same:>12s}")

# ─────────────────────────────────────────────────────────────
# SECTION 7: DASHBOARD (9 panels)
# ─────────────────────────────────────────────────────────────

print("\n🎨 Generating dashboard...")

COMM_COLORS = plt.cm.Set2(np.linspace(0, 1, max(len(communities), 4)))

fig = plt.figure(figsize=(26, 20))
fig.suptitle(
    "AMR Co-resistance Network Analysis — REAL CLINICAL DATA\n"
    "Phi (φ) coefficient network | 34 antibiotics | BH-FDR corrected\n"
    "example_isolates (AMR R package) | Block 2 Day 16 — Final Day",
    fontsize=15, fontweight="bold", y=0.99
)

# ── Plot 1: Full co-resistance network ──
ax1 = fig.add_subplot(3, 3, 1)
ax1.set_facecolor("#0D1117")
pos  = nx.spring_layout(G, weight="weight", seed=42, k=2.5)
node_colors = [G.nodes[n]["color"] for n in G.nodes()]
node_sizes  = [200 + G.degree(n)*30 for n in G.nodes()]
# Positive edges (co-resistance)
pos_edge_w = [G[u][v]["phi"]*5 for u,v in pos_edges]
nx.draw_networkx_edges(G, pos, edgelist=pos_edges,
                       width=pos_edge_w, edge_color="#E74C3C",
                       alpha=0.55, ax=ax1)
# Negative edges (antagonistic)
neg_edge_w = [abs(G[u][v]["phi"])*3 for u,v in neg_edges]
nx.draw_networkx_edges(G, pos, edgelist=neg_edges,
                       width=neg_edge_w, edge_color="#3498DB",
                       alpha=0.4, ax=ax1, style="dashed")
nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                       node_size=node_sizes, ax=ax1, alpha=0.92,
                       edgecolors="white", linewidths=0.6)
nx.draw_networkx_labels(G, pos, ax=ax1, font_size=6.5,
                        font_color="white", font_weight="bold")
ax1.set_title("Co-resistance Network\n(red=co-resistant, blue=antagonistic)",
              fontweight="bold", fontsize=10, color="white", pad=5)
patches = [mpatches.Patch(color=c, label=k)
           for k, c in CLASS_COLORS.items()
           if any(AB_CLASS.get(n,"Other")==k for n in G.nodes())]
ax1.legend(handles=patches, fontsize=6, loc="lower left",
           framealpha=0.3, labelcolor="white", facecolor="#1a1a2e")

# ── Plot 2: Community network ──
ax2 = fig.add_subplot(3, 3, 2)
ax2.set_facecolor("#0D1117")
comm_node_colors = [COMM_COLORS[community_map[n]]
                    for n in G.nodes()]
nx.draw_networkx_edges(G, pos, edgelist=pos_edges,
                       width=1.2, edge_color="#888888",
                       alpha=0.3, ax=ax2)
nx.draw_networkx_nodes(G, pos, node_color=comm_node_colors,
                       node_size=node_sizes, ax=ax2, alpha=0.92,
                       edgecolors="white", linewidths=0.6)
nx.draw_networkx_labels(G, pos, ax=ax2, font_size=6.5,
                        font_color="white", font_weight="bold")
ax2.set_title(f"Community Structure\n"
              f"{len(communities)} communities | modularity={modularity:.3f}",
              fontweight="bold", fontsize=10, color="white", pad=5)
comm_patches = [mpatches.Patch(
    color=COMM_COLORS[i],
    label=f"Community {i+1} (n={len(list(communities)[i])})")
    for i in range(len(communities))]
ax2.legend(handles=comm_patches, fontsize=7, loc="lower left",
           framealpha=0.3, labelcolor="white", facecolor="#1a1a2e")

# ── Plot 3: Phi coefficient heatmap (top 18 ABs) ──
ax3 = fig.add_subplot(3, 3, 3)
top18 = centrality_df["Antibiotic"].head(18).tolist()
phi_sub = phi_df.loc[top18, top18]
mask_diag = np.eye(len(top18), dtype=bool)
sns.heatmap(phi_sub, ax=ax3, cmap="RdBu_r", center=0,
            vmin=-1, vmax=1, annot=True, fmt=".2f",
            linewidths=0.3, mask=mask_diag,
            cbar_kws={"label":"φ coefficient","shrink":0.75},
            annot_kws={"size":6})
ax3.tick_params(axis="both", labelsize=7)
ax3.set_title("Phi Co-resistance Heatmap\n(Top 18 antibiotics by degree)",
              fontweight="bold", fontsize=10)

# ── Plot 4: Degree distribution ──
ax4 = fig.add_subplot(3, 3, 4)
degrees = sorted([d for _, d in G.degree()], reverse=True)
degree_counts = pd.Series(degrees).value_counts().sort_index()
ax4.bar(degree_counts.index, degree_counts.values,
        color="#E74C3C", edgecolor="black", linewidth=0.4, alpha=0.85)
ax4.set_xlabel("Degree (number of co-resistant pairs)")
ax4.set_ylabel("Number of Antibiotics")
ax4.set_title(f"Degree Distribution\n(mean={np.mean(degrees):.1f}, "
              f"max={max(degrees)}, density={nx.density(G):.3f})",
              fontweight="bold", fontsize=10)
# Annotate top hubs
for ab, d in sorted(G.degree(), key=lambda x: x[1], reverse=True)[:3]:
    ax4.annotate(ab, (d, degree_counts.get(d,0)),
                 xytext=(5, 5), textcoords="offset points", fontsize=8)

# ── Plot 5: Top hub centrality comparison ──
ax5 = fig.add_subplot(3, 3, 5)
top10 = centrality_df.head(10)
x5    = np.arange(len(top10))
w5    = 0.28
b1 = ax5.bar(x5-w5, top10["Degree_cent"],   w5, label="Degree",
             color="#E74C3C", alpha=0.85, edgecolor="black", linewidth=0.4)
b2 = ax5.bar(x5,    top10["Betweenness"]*3, w5, label="Betweenness×3",
             color="#3498DB", alpha=0.85, edgecolor="black", linewidth=0.4)
b3 = ax5.bar(x5+w5, top10["Eigenvector"],   w5, label="Eigenvector",
             color="#2ECC71", alpha=0.85, edgecolor="black", linewidth=0.4)
ax5.set_xticks(x5)
ax5.set_xticklabels(top10["Antibiotic"], fontsize=8, rotation=30, ha="right")
ax5.set_ylabel("Centrality Score")
ax5.set_title("Hub Antibiotic Centrality\n(Top 10 by degree)",
              fontweight="bold", fontsize=10)
ax5.legend(fontsize=8)

# ── Plot 6: Top co-resistance pairs (φ bar) ──
ax6 = fig.add_subplot(3, 3, 6)
top_pairs = edge_df.head(15)
pair_labels = [f"{r['AB1']}–{r['AB2']}" for _,r in top_pairs.iterrows()]
pair_colors = ["#E74C3C" if r["class1"]==r["class2"]
               else "#F39C12" for _,r in top_pairs.iterrows()]
ax6.barh(range(len(top_pairs)), top_pairs["phi"].values[::-1],
         color=pair_colors[::-1], edgecolor="black", linewidth=0.3, alpha=0.87)
ax6.set_yticks(range(len(top_pairs)))
ax6.set_yticklabels(pair_labels[::-1], fontsize=8)
ax6.set_xlabel("Phi (φ) coefficient")
ax6.set_title("Top 15 Co-resistance Pairs\n(red=same class, orange=cross-class)",
              fontweight="bold", fontsize=10)
ax6.axvline(0.5, color="gray", lw=1, linestyle="--", alpha=0.5)
ax6.legend(handles=[
    mpatches.Patch(color="#E74C3C", label="Same class"),
    mpatches.Patch(color="#F39C12", label="Cross-class"),
], fontsize=8)

# ── Plot 7: Co-resistance % heatmap (selected clinical pairs) ──
ax7 = fig.add_subplot(3, 3, 7)
# Class-level co-resistance summary
classes_present = list(set(AB_CLASS.get(ab,"Other") for ab in FOCUS_ABS
                           if AB_CLASS.get(ab,"Other") != "Other"))[:8]
class_pairs = pd.DataFrame(index=classes_present, columns=classes_present, dtype=float)
for c1 in classes_present:
    for c2 in classes_present:
        abs1 = [ab for ab in FOCUS_ABS if AB_CLASS.get(ab,"Other")==c1]
        abs2 = [ab for ab in FOCUS_ABS if AB_CLASS.get(ab,"Other")==c2]
        phis = []
        for a1 in abs1:
            for a2 in abs2:
                if a1 != a2:
                    i, j = FOCUS_ABS.index(a1), FOCUS_ABS.index(a2)
                    phis.append(phi_mat[i,j])
        class_pairs.loc[c1,c2] = np.mean(phis) if phis else 0

sns.heatmap(class_pairs.astype(float), ax=ax7, cmap="RdBu_r",
            center=0, vmin=-0.5, vmax=0.5, annot=True, fmt=".2f",
            linewidths=0.5, cbar_kws={"label":"Mean φ","shrink":0.8},
            annot_kws={"size":8})
ax7.tick_params(axis="x", rotation=35, labelsize=7)
ax7.tick_params(axis="y", rotation=0,  labelsize=7)
ax7.set_title("Class-level Co-resistance\n(mean φ between antibiotic classes)",
              fontweight="bold", fontsize=10)

# ── Plot 8: Clustering coefficient by class ──
ax8 = fig.add_subplot(3, 3, 8)
clust_by_class = {}
for cls in CLASS_COLORS:
    nodes_cls = [n for n in G.nodes() if G.nodes[n]["ab_class"]==cls]
    if nodes_cls:
        clust_by_class[cls] = np.mean([clustering[n] for n in nodes_cls])
clust_series = pd.Series(clust_by_class).sort_values(ascending=False)
clust_series = clust_series[clust_series > 0]
bar_cols8 = [CLASS_COLORS.get(c,"#BDC3C7") for c in clust_series.index]
bars8 = ax8.bar(range(len(clust_series)), clust_series.values,
                color=bar_cols8, edgecolor="black", linewidth=0.4, alpha=0.87)
for bar, val in zip(bars8, clust_series.values):
    ax8.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.005,
             f"{val:.2f}", ha="center", fontsize=8, fontweight="bold")
ax8.set_xticks(range(len(clust_series)))
ax8.set_xticklabels(clust_series.index, rotation=30, ha="right", fontsize=8)
ax8.set_ylabel("Mean Clustering Coefficient")
ax8.set_title("Network Clustering by Class\n(higher = tighter co-resistance cluster)",
              fontweight="bold", fontsize=10)
ax8.axhline(nx.average_clustering(G, weight="weight"),
            color="black", lw=1.5, linestyle="--",
            label=f"Network avg ({nx.average_clustering(G,weight='weight'):.3f})")
ax8.legend(fontsize=8)

# ── Plot 9: Summary table ──
ax9 = fig.add_subplot(3, 3, 9)
ax9.axis("off")
rows9 = [
    ["Antibiotics (nodes)",   str(G.number_of_nodes())],
    ["Significant edges",     str(G.number_of_edges())],
    ["Network density",       f"{nx.density(G):.4f}"],
    ["Avg clustering coeff",  f"{nx.average_clustering(G, weight='weight'):.4f}"],
    ["Positive co-R edges",   str(len(pos_edges))],
    ["Negative (antagonist)", str(len(neg_edges))],
    ["Communities detected",  str(len(communities))],
    ["Modularity",            f"{modularity:.4f}"],
    ["Top hub AB",            centrality_df.iloc[0]["Antibiotic"]],
    ["Top hub degree",        str(int(centrality_df.iloc[0]["Degree"]))],
    ["Strongest pair φ",      f"{edge_df.iloc[0]['AB1']}–{edge_df.iloc[0]['AB2']} "
                               f"(φ={edge_df.iloc[0]['phi']:.3f})"],
    ["FDR threshold",         "q < 0.05, |φ| ≥ 0.20"],
]
tbl9 = ax9.table(cellText=rows9,
                  colLabels=["Metric","Value"],
                  cellLoc="left", loc="center")
tbl9.auto_set_font_size(False); tbl9.set_fontsize(9); tbl9.scale(1.7, 1.95)
for j in range(2):  tbl9[(0,j)].set_facecolor("#2C3E50")
for j in range(2):  tbl9[(0,j)].set_text_props(color="white", fontweight="bold")
tbl9[(9,0)].set_facecolor("#E74C3C"); tbl9[(9,0)].set_text_props(color="white")
tbl9[(9,1)].set_facecolor("#E74C3C"); tbl9[(9,1)].set_text_props(color="white")
ax9.set_title("Co-resistance Network Summary\n(Block 2 — Final Day)",
              fontweight="bold", fontsize=11, pad=20)

plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig("outputs/coresistance_network_dashboard.png",
            dpi=150, bbox_inches="tight")
plt.close()
print("✅ Dashboard saved → outputs/coresistance_network_dashboard.png")

# ─────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────

print("\n" + "="*62)
print("FINAL SUMMARY — CO-RESISTANCE NETWORK | BLOCK 2 COMPLETE 🎉")
print("="*62)
print(f"\nNetwork topology:")
print(f"  Nodes     : {G.number_of_nodes()} antibiotics")
print(f"  Edges     : {G.number_of_edges()} significant co-resistance links")
print(f"  Density   : {nx.density(G):.4f}")
print(f"  Avg clust : {nx.average_clustering(G,weight='weight'):.4f}")
print(f"  Communities: {len(communities)} (modularity={modularity:.4f})")
print(f"\nTop 5 hub antibiotics:")
for _, row in centrality_df.head(5).iterrows():
    print(f"  {row['Antibiotic']:5s} ({row['Class']:15s}): "
          f"degree={int(row['Degree'])}, betweenness={row['Betweenness']:.4f}")
print(f"\nStrongest co-resistance pairs (φ):")
for _, row in edge_df.head(5).iterrows():
    print(f"  {row['AB1']}–{row['AB2']:5s}: φ={row['phi']:.4f} | "
          f"co-R rate={row['coresist_pct']:.1f}%")
print("\n" + "="*62)
print("✅ Block 2 (Days 09–16) COMPLETE! 🎉")
print("="*62)
