"""
This script reads the EECU data from the fixed-width text file, cleans and normalizes it, and then performs analysis 
to understand the distribution of EECU hours across different OpenET models. 
It generates summary statistics and a boxplot to visualize the differences in EECU hours between models. 
The cleaned data and results are saved to CSV files for further use.

EECU Data Link: https://openet-dri.appspot.com/
Authors: Dr. Sayantan Majumdar, Charles Morton (Desert Research Institute)
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


if __name__ == "__main__":
    # Read the fixed-width file, skipping title, blank line, and both header rows
    df = pd.read_fwf(
        "eecu_data/eecu.txt",
        skiprows=[0, 1, 2, 3],
        header=None,
        names=["Date", "Model", "Tool", "Completed", "Time (min)", "EECU (hours)"],
    )


    # Forward-fill Date and Model (blank rows inherit from above)
    df["Date"] = pd.to_datetime(df["Date"].ffill())
    df["Model"] = df["Model"].ffill()

    # Normalize all model name variants to canonical names
    model_map = {
        "openet-eemetric": "eeMETRIC", "EEMETRIC": "eeMETRIC", "eemetric": "eeMETRIC",
        "openet-geesebal": "geeSEBAL", "GEESEBAL": "geeSEBAL", "geesebal": "geeSEBAL",
        "openet-ptjpl": "PT-JPL", "PTJPL": "PT-JPL", "ptjpl": "PT-JPL",
        "openet-ssebop": "SSEBop", "SSEBOP": "SSEBop", "ssebop": "SSEBop",
        "openet-sims": "SIMS", "sims": "SIMS",
        "openet-disalexi": "DisALEXI", "DISALEXI": "DisALEXI", "disalexi": "DisALEXI",
    }
    df["Model"] = df["Model"].replace(model_map)

    valid_openet_models = [
        "eeMETRIC",
        "geeSEBAL",
        "PT-JPL",
        "SSEBop",
        "SIMS",
        "DisALEXI"
    ]

    cutoff_date = pd.to_datetime("2026-12-31")
    cutoff_date_geesebal = pd.to_datetime("2024-03-26") # geesebal has some different runs after this period
    eecu_output_dir = "eecu_output/"

    # Filter to valid models and relevant export tools
    df = df[
        (df.Model.isin(valid_openet_models))
        & (
            ((df.Model == 'geeSEBAL') & (df.Date <= cutoff_date_geesebal))
            | ((df.Model != 'geeSEBAL') & (df.Date <= cutoff_date))
        )
        & (
            (df.Tool == 'scene_image_wrs2_export')
            | ((df.Model == 'DisALEXI') & (df.Tool == 'tair_image_wrs2_export'))
        ) 
    ]

    # Drop rows with NaN EECU hours (pre-Aug 2023 data lacks EECU)
    df = df.dropna(subset=["EECU (hours)"])

    # Aggregate to one row per (Date, Model)
    # Sum EECU hours for DisALEXI (multiple tools), mean for others (duplicate exports)
    disalexi = df[df["Model"] == "DisALEXI"]
    others = df[df["Model"] != "DisALEXI"]

    # Only keep DisALEXI dates where both scene_image and tair_image tools are present
    disalexi_tool_counts = disalexi.groupby("Date")["Tool"].nunique()
    disalexi_both_dates = disalexi_tool_counts[disalexi_tool_counts == 2].index
    disalexi = disalexi[disalexi["Date"].isin(disalexi_both_dates)]

    disalexi_agg = disalexi.groupby(["Date", "Model"], as_index=False).agg(
        {"EECU (hours)": "sum", "Completed": "sum", "Time (min)": "max"}
    )
    others_agg = others.groupby(["Date", "Model"], as_index=False).agg(
        {"EECU (hours)": "mean", "Completed": "sum", "Time (min)": "max"}
    )
    df = pd.concat([others_agg, disalexi_agg], ignore_index=True) \
        .sort_values(by=["Date", "Model"], ascending=[False, True]) \
        .reset_index(drop=True)
    df_plot = df[["Date", "Model", "EECU (hours)"]].copy()
    df_plot = df_plot[df_plot["EECU (hours)"] > 0]
    os.makedirs(eecu_output_dir, exist_ok=True)
    df_plot.to_csv(f"{eecu_output_dir}eecu_data.csv", index=False)
    for model in df_plot["Model"].unique():
        model_data = df_plot[df_plot["Model"] == model]
        model_data.to_csv(f"{eecu_output_dir}eecu_{model.lower()}.csv", index=False)

    # Show start and end dates per model
    date_range = df_plot.groupby("Model")["Date"].agg(["min", "max"])
    date_range.columns = ["Start Date", "End Date"]
    print("\nExport Date Range per Model")
    print("=" * 50)
    print(date_range.to_string())
    print()

    # Compute average EECU hours per model
    avg_eecu = df_plot.groupby("Model")["EECU (hours)"].mean().sort_values(ascending=False)

    print("Average EECU (hours) per OpenET Model")
    print("=" * 45)
    print(avg_eecu.to_string())

    # Descriptive statistics per model
    desc_stats = df_plot.groupby("Model")["EECU (hours)"].describe()
    print("\nDescriptive Statistics per OpenET Model")
    print("=" * 70)
    print(desc_stats.to_string())
    print()

    # Boxplot of EECU hours per model
    model_order = avg_eecu.index.tolist()
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=df_plot, x="Model", y="EECU (hours)", order=model_order, linewidth=1.2, ax=ax)
    ax.set_yscale("log")
    ax.set_xlabel("OpenET Model", fontsize=12)
    ax.set_ylabel("EECU (hours) [log scale]", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"{eecu_output_dir}eecu_boxplot.png", dpi=600)
