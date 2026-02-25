"""
This script reads the EECU data from the fixed-width text file, cleans and normalizes it, and then performs analysis 
to understand the distribution of EECU hours across different OpenET models. 
It generates summary statistics and a boxplot to visualize the differences in EECU hours between models. 
The cleaned data and results are saved to CSV files for further use.

EECU Data Link: https://openet-dri.appspot.com/ (this link may not be active in the future, but the data file is included in the repository for reproducibility)
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
    df_plot = df[["Date", "Model", "EECU (hours)", "Completed"]].copy()
    df_plot = df_plot[(df_plot["EECU (hours)"] > 0) & (df_plot["Completed"] > 0)]
    df_plot["EECU (hours/task)"] = df_plot["EECU (hours)"] / df_plot["Completed"]
    df_plot["EECU (secs/task)"] = df_plot["EECU (hours/task)"] * 3600
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

    # Compute average EECU hours per task per model
    avg_eecu = df_plot.groupby("Model")["EECU (hours/task)"].mean().sort_values(ascending=False)
    avg_eecu_secs = df_plot.groupby("Model")["EECU (secs/task)"].mean().sort_values(ascending=False)

    print("Average EECU (hours/task) per OpenET Model")
    print("=" * 45)
    print(avg_eecu.to_string())
    print("\nAverage EECU (secs/task) per OpenET Model")
    print("=" * 45)
    print(avg_eecu_secs.to_string())
    # plot_col = "EECU (hours/task)" 
    plot_col = "EECU (secs/task)" # uncomment the above and comment this to plot hours instead of seconds

    # Descriptive statistics per model
    desc_stats = df_plot.groupby("Model")[plot_col].describe()
    print("\nDescriptive Statistics per OpenET Model")
    print("=" * 70)
    print(desc_stats.to_string())
    print()

    # Boxplot of EECU seconds per task per model
    model_order = avg_eecu.index.tolist()[::-1]  # Order by average EECU, lowest to highest
    sns.set_theme(style='whitegrid', font_scale=1.2)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=df_plot, x="Model", y=plot_col, 
        order=model_order, linewidth=1.2, ax=ax,
        # color=sns.color_palette("tab10")[0]
    )
    ax.set_yscale("log")    
    ax.set_xlabel("OpenET Model", fontsize=14)
    eecu_label = plot_col.replace("secs/task", "seconds task$^{-1}$")
    ylab = f"{eecu_label}, log scale)"
    ax.set_ylabel(ylab, fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    plt.savefig(f"{eecu_output_dir}eecu_boxplot.png", dpi=600)

    # Barplot of mean EECU hours per task per model with error bars
    summary = df_plot.groupby("Model")[plot_col].agg(["mean", "std"]).loc[model_order]
    # Clip lower error bars to avoid negative values on log scale
    lower_err = summary['std'].clip(upper=summary['mean'] * 0.99)
    upper_err = summary['std']
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(
        summary.index,
        summary['mean'],
        yerr=[lower_err, upper_err],
        capsize=5,
        edgecolor='black',
        linewidth=0.8,
    )
    # Annotate bars with mean values
    for bar, mean_val in zip(bars, summary['mean']):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() / 2,
            f'{int(round(mean_val))}',
            ha='center', va='center', fontweight='bold', color='black', fontsize=14,
        )
    ax.set_yscale("log")
    ax.set_xlabel("OpenET Model", fontsize=14)
    eecu_label = plot_col.replace("secs/task", "seconds task$^{-1}$")
    ylab = f"{eecu_label}, log scale)"
    ax.set_ylabel(ylab, fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{eecu_output_dir}eecu_barplot.png", dpi=600)
