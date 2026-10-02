# %% import packages
import os
import pandas as pd

# %% function creation
# input: dataset from clinvar (tsv file)
# output: dataset filtered (csv file)

def explorar_dataset(nombre_archivo):
    ruta_tsv = os.path.join("data", "raw", nombre_archivo)
    nombre_dataset = os.path.splitext(os.path.basename(ruta_tsv))[0]
    ruta_csv = os.path.join("data", "raw", nombre_dataset + ".csv")
    df = pd.read_csv(ruta_tsv, sep="\t", dtype=str, index_col=False, engine="python")
    df.columns = [c.strip().lstrip("#").strip() for c in df.columns]
    df = df.loc[:, [c for c in df.columns if not c.startswith("Unnamed")]]
    df.to_csv(ruta_csv, index=False)
    df = pd.read_csv(ruta_csv, dtype=str, index_col=False)
    valores_unicos = df["Germline classification"].dropna().unique().tolist()
    return nombre_dataset, valores_unicos

def filtrar_por_clasificacion(nombre_archivo, clasificacion):
    ruta_entrada = os.path.join("data", "raw", nombre_archivo)
    df = pd.read_csv(ruta_entrada, sep="\t", dtype=str)
    df_filtrado = df[
        (df["Germline classification"] == clasificacion)
        & (df["Molecular consequence"].fillna("").str.contains("missense variant", case=False))
        & (df["Variant type"].fillna("").str.strip().str.lower() == "single nucleotide variant")
    ]
    os.makedirs(os.path.join("data", "processed"), exist_ok=True)
    nombre_dataset = os.path.splitext(os.path.basename(ruta_entrada))[0]
    ruta_salida = os.path.join("data", "processed", nombre_dataset + ".csv")
    df_filtrado.to_csv(ruta_salida, index=False)
    return df_filtrado

# %% 
nombre_dataset, clasificaciones = explorar_dataset("HMGCR.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("SLCO1B1.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("ABCG2.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("ABCB1.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("CYP2C9.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("CYP3A4.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("UGT1A1.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("LDLR.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("PCSK9.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("SCARB1.tsv")
print(nombre_dataset)
print(clasificaciones)

nombre_dataset, clasificaciones = explorar_dataset("APOE.tsv")
print(nombre_dataset)
print(clasificaciones)

# %%
df_resultado = filtrar_por_clasificacion("HMGCR.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("SLCO1B1.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("ABCG2.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("ABCB1.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("CYP2C9.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("CYP3A4.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("UGT1A1.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("LDLR.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("PCSK9.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("SCARB1.tsv", "Pathogenic")

df_resultado = filtrar_por_clasificacion("APOE.tsv", "Pathogenic")
