import os
import pubchempy as pcp
from mordred import Calculator
from rdkit import Chem
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from mordred.MoeType import SlogP_VSA
from mordred.HydrogenBond import HBondDonor, HBondAcceptor
from mordred.RingCount import RingCount
from mordred.Weight import Weight
from mordred.SLogP import SLogP
from mordred.TopoPSA import TopoPSA

cs = pcp.get_compounds('Aspirin', 'name')
target = cs[0]
print(f"CID: {target.cid}")

props = pcp.get_properties(['CanonicalSMILES', 'IsomericSMILES'], target.cid, 'cid')

smiles = props[0].get('SMILES')
print(smiles)
similar = pcp.get_compounds(smiles, namespace='smiles', searchtype='similarity', Threshold=90, listkey_count=1000)

important_desc = [
    SlogP_VSA(k=2),  # SlogP_VSA2
    SlogP_VSA(k=8),  # SlogP_VSA8
    HBondDonor(),
    HBondAcceptor(),
    RingCount(),
    Weight(),
    SLogP(),
    TopoPSA()
]


records = list()
calc = Calculator(important_desc, ignore_3D=True)

try:
    for c in similar:
        similar_props = pcp.get_properties(['CanonicalSMILES', 'IsomericSMILES'], c.cid, 'cid')
        similar_smiles = similar_props[0].get('SMILES')

        mol = Chem.MolFromSmiles(similar_smiles)

        if mol is None:
            print(f"[WARN] CID {c.cid}: Invalid SMILES skipped.")
            continue
        
        desc = calc(mol).asdict()
        print(desc)
        desc["CID"] = c.cid
        desc["SMILES"] = similar_smiles
        desc["IUPAC"] = c.iupac_name
        records.append(desc)
        print(f"[INFO] Processed CID {c.cid}")

except KeyError as e:
    print(f"{c.iupac_name} has no SMILES!!")

print(records)
df = pd.DataFrame(records)
print(df.info())
print(df.columns)

print(f"[INFO] Mordred descriptors calculated for {len(df)} compounds")
print(df)

# Descriptor의 중요한 feature들 
desc_cols = ["SlogP_VSA2", "SlogP_VSA8", "nHBDon", "nHBAcc", "nRing", "MW", "SLogP", "TopoPSA(NO)"]

# 원본 유지지
df_clean = df[desc_cols].apply(pd.to_numeric, errors='coerce').dropna()


# Scaled DataFrame
#scaler = StandardScaler()
#scaled = pd.DataFrame(scaler.fit_transform(df_clean), columns=desc_cols)

print(f"[INFO] Descriptor DataFrame shape: {df_clean.shape}")
print(df_clean)

print(df_clean.corr())
print(df_clean.info())


os.makedirs("plots", exist_ok=True)

# 변수 형태에 따라 EDA 다르게
discrete_cols = ["nHBDon", "nHBAcc", "nRing"]
continuous_cols = [c for c in df_clean.columns if c not in discrete_cols]

# ---------------------------
# 연속형 변수: Boxplot + KDE
# ---------------------------
for col in continuous_cols:
    series = df_clean[col].dropna()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sns.boxplot(y=series, ax=axes[0], color='skyblue')
    axes[0].set_title(f"{col} - Boxplot")

    sns.kdeplot(series, fill=True, ax=axes[1])
    axes[1].set_title(f"{col} - KDE Plot (Continuous)")

    plt.tight_layout()
    plt.savefig(f"plots/{col}_box_kde.png", dpi=300)
    plt.close(fig)

# ---------------------------
# 이산형 변수: Barplot
# ---------------------------
for col in discrete_cols:
    series = df_clean[col].dropna()

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.countplot(x=series, hue=series, palette='muted', ax=ax, legend=False)
    ax.set_title(f"{col} - Bar Plot (Discrete)")

    plt.tight_layout()
    plt.savefig(f"plots/{col}_bar.png", dpi=300)
    plt.close(fig)

print("[INFO] All plots saved in the 'plots/' folder successfully.")