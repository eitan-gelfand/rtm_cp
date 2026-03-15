"""
D-prime Corrected Workflow
Combines correct D-prime calculation from final_Dpirme_CR.ipynb 
with file I/O system from local_dprime_workflow.py

KEY DIFFERENCE: Uses notebook's compute_dprime_cr logic (NO .clip() on Z-scores after 0.99/0.01 correction)
"""

from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


DEFAULT_DATA_FILE = "full_data_cp.csv"


def find_project_root(start: Path | None = None) -> Path:
    """
    Walk upward until a folder containing /data is found.
    This makes the script portable inside the repo.
    """
    start = (start or Path(__file__).resolve()).resolve()
    for candidate in [start.parent, *start.parents]:
        if (candidate / "data").exists():
            return candidate
    # fallback: script folder
    return start.parent


ROOT = find_project_root()
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"
PLOTS_DIR = OUTPUT_DIR / "plots"
OUTPUT_DIR.mkdir(exist_ok=True)
PLOTS_DIR.mkdir(exist_ok=True, parents=True)


def data_path(*parts: str) -> Path:
    return DATA_DIR.joinpath(*parts)


def output_path(*parts: str) -> Path:
    return OUTPUT_DIR.joinpath(*parts)


def plot_path(*parts: str) -> Path:
    return PLOTS_DIR.joinpath(*parts)


def rename_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Renames duplicate columns by appending a suffix."""
    cols = pd.Series(df.columns)
    duplicates = cols[cols.duplicated()].unique()
    if len(duplicates) == 0:
        return df

    new_cols = []
    counts = {}
    for col in df.columns:
        if col not in counts:
            counts[col] = 0
            new_cols.append(col)
        else:
            counts[col] += 1
            new_cols.append(f"{col}_{counts[col]}")
    df = df.copy()
    df.columns = new_cols
    return df


def split_into_blocks(
    df: pd.DataFrame,
    grouping_cols: list[str],
    acc_col: str,
    num_of_blocks: int,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Splits each group in the DataFrame into a specified number of random blocks.
    (From notebook final_Dpirme_CR.ipynb, unchanged)
    """
    print("[INFO] Splitting FA data into random blocks...")
    results = []

    for group_vals, group_df in df.groupby(grouping_cols):
        if not isinstance(group_vals, tuple):
            group_vals = (group_vals,)

        group_size = len(group_df)

        # Shuffle to ensure randomness
        group_df_shuffled = group_df.sample(frac=1, random_state=random_state).reset_index(drop=True)

        # Determine block sizes
        base_block_size = group_size // num_of_blocks
        leftover = group_size % num_of_blocks
        block_sizes = [base_block_size + 1 if i < leftover else base_block_size for i in range(num_of_blocks)]

        start_idx = 0
        for block_num, size in enumerate(block_sizes, 1):
            if size == 0:
                continue  # Skip if no data for this block

            end_idx = start_idx + size
            block_slice = group_df_shuffled.iloc[start_idx:end_idx]

            block_len = len(block_slice)
            block_sum = block_slice[acc_col].sum()

            # Prepare result row with grouping columns
            result_row = {col: group_vals[i] for i, col in enumerate(grouping_cols)}
            result_row['Block_Number'] = block_num
            result_row['ACC_len'] = block_len
            result_row['ACC_sum_FA'] = block_sum

            # Append to results
            results.append(result_row)

            start_idx = end_idx

    results_df = pd.DataFrame(results)
    return results_df


def process_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Processes the Hits data with hardcoded columns.
    (From notebook final_Dpirme_CR.ipynb, adapted for automation)
    """
    print("\n[INFO] Processing Hits...")
    acc_col = "ACC"
    grouping_cols = ["Subject", "ExperimentName", "Regression", "Range", "Age", "Group"]

    # Step 1: Check if 'Regression' column exists
    if 'Regression' not in df.columns:
        raise KeyError("[ERROR] Column 'Regression' not found in the DataFrame for Hits processing.")

    # Step 2: Filter rows where 'Regression' != 'Null'
    df_filtered = df[df['Regression'] != 'Null'].copy()
    print(f"[INFO] After filtering 'Regression' != 'Null', shape={df_filtered.shape}")

    if df_filtered.empty:
        raise ValueError("[ERROR] No data left after filtering 'Regression' != 'Null' for Hits processing.")

    # Step 3: Convert accuracy to binary (1 if 1, else 0)
    df_filtered.loc[:, acc_col] = df_filtered[acc_col].apply(lambda x: 1 if x == 1 else 0)

    # Step 4: Group and aggregate
    grouped = df_filtered.groupby(grouping_cols, dropna=False)[acc_col].agg(['count', 'sum']).reset_index()
    grouped.rename(columns={'count': 'ACC_count', 'sum': 'ACC_sum_Hits'}, inplace=True)

    print(f"[INFO] Hits: Final shape={grouped.shape}")
    print("[INFO] Sample of Hits results:")
    print(grouped.head(10))

    return grouped


def process_fa(df: pd.DataFrame, num_of_blocks: int) -> pd.DataFrame:
    """
    Processes the False Alarms (FA) data with hardcoded columns.
    (From notebook final_Dpirme_CR.ipynb, adapted for automation)
    """
    print("\n[INFO] Processing False Alarms (FA)...")
    acc_col = "ACC"
    grouping_cols = ["Subject", "ExperimentName", "Age", "Group"]

    # Step 1: Check if 'Regression' column exists
    if 'Regression' not in df.columns:
        raise KeyError("[ERROR] Column 'Regression' not found in the DataFrame for FA processing.")

    # Step 2: Filter rows where 'Regression' == 'Null'
    df_filtered = df[df['Regression'] == 'Null'].copy()
    print(f"[INFO] After filtering 'Regression' == 'Null', shape={df_filtered.shape}")

    if df_filtered.empty:
        raise ValueError("[ERROR] No data left after filtering 'Regression' == 'Null' for FA processing.")

    # Step 3: Convert accuracy to binary (0 if 1, else 1)
    df_filtered.loc[:, acc_col] = df_filtered[acc_col].apply(lambda x: 0 if x == 1 else 1)

    # Step 4: Split each group into blocks based on the number_of_blocks
    fa_blocks = split_into_blocks(df_filtered, grouping_cols, acc_col, num_of_blocks, random_state=42)

    print(f"[INFO] FA: Final shape={fa_blocks.shape}")
    print("[INFO] Sample of FA results:")
    print(fa_blocks.head(10))

    return fa_blocks


def merge_dataframes(hits_df: pd.DataFrame, fa_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merges Hits and FA DataFrames on common grouping columns.
    (From notebook final_Dpirme_CR.ipynb, unchanged)
    """
    print("\n[INFO] Merging Hits and FA DataFrames...")
    # Identify aggregation-related columns to exclude from common columns
    aggregation_cols = {'ACC_count', 'ACC_sum_Hits', 'ACC_sum_FA', 'Block_Number', 'ACC_len'}
    common_cols = list((set(hits_df.columns) & set(fa_df.columns)) - aggregation_cols)

    if not common_cols:
        raise ValueError("[ERROR] No common grouping columns found for merging.")

    print(f"[INFO] Common grouping columns for merging: {common_cols}")

    # Merge DataFrames on all common grouping columns
    combined_df = pd.merge(hits_df, fa_df, on=common_cols, how='outer', suffixes=('_Hits','_FA'))
    print(f"[INFO] Merged DataFrame shape: {combined_df.shape}")
    print("[INFO] Sample of Merged DataFrame:")
    print(combined_df.head(10))

    return combined_df


def compute_dprime_cr(
    df: pd.DataFrame,
    hit_sum_col: str = 'ACC_sum_Hits',
    hit_count_col: str = 'ACC_count',
    fa_sum_col: str = 'ACC_sum_FA',
    fa_count_col: str = 'ACC_len'
) -> pd.DataFrame:
    """
    Computes Dprime and CR based on Hits and FA data.
    
    **KEY DIFFERENCE FROM local_dprime_workflow.py:**
    Uses notebook's logic - NO .clip() on Z-scores after 0.99/0.01 replacement.
    This is the mathematically correct approach.
    
    (From notebook final_Dpirme_CR.ipynb, NOT from local_dprime_workflow.py)
    """
    print("\n[INFO] Computing Dprime and CR...")
    Z = norm.ppf

    # Calculate proportions - fix to avoid 1s and 0s
    df['HitProp'] = (
        df[hit_sum_col] / df[hit_count_col].replace(0, np.nan)
    ).replace({1.0: 0.99, 0.0: 0.01})

    df['FaProp'] = (
        df[fa_sum_col] / df[fa_count_col].replace(0, np.nan)
    ).replace({1.0: 0.99, 0.0: 0.01})

    # ⭐ CRITICAL: NO .clip() applied here (unlike local_dprime_workflow.py)
    # Apply Z-score transformation directly to the 0.99/0.01 corrected proportions
    df['ZHitProp'] = Z(df['HitProp'])
    df['ZFaProp'] = Z(df['FaProp'])

    # Compute Dprime and CR
    df['Dprime'] = df['ZHitProp'] - df['ZFaProp']
    df['CR'] = -((df['ZHitProp'] + df['ZFaProp']) / 2)

    print("[INFO] Dprime and CR computed successfully.")
    print("[INFO] Sample of Dprime and CR:")
    print(df[['Dprime', 'CR']].head(10))

    return df


def filter_blocks_by_range(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters the combined DataFrame to include only rows where the 'Block_Number' matches
    the mapped block number for their 'Range'. Removes the 'Block_Number' column.
    (From notebook final_Dpirme_CR.ipynb, unchanged)
    """
    print("\n[INFO] Filtering Blocks by Range...")

    # Verify necessary columns exist
    required_columns = {'Range', 'Block_Number'}
    if not required_columns.issubset(df.columns):
        missing = required_columns - set(df.columns)
        print(f"[ERROR] Missing required columns for filtering: {missing}")
        return None

    # Extract unique non-zero ranges
    unique_ranges = df['Range'].unique()
    unique_ranges = [r for r in unique_ranges if r != 0]
    print(f"[INFO] Unique non-zero Ranges: {unique_ranges}")

    # Number of blocks is the number of unique non-zero ranges
    num_of_blocks = len(unique_ranges)
    print(f"[INFO] Number of blocks set to: {num_of_blocks}")

    # Create a mapping from Range to Block_Number
    sorted_unique_ranges = sorted(unique_ranges)
    range_to_block = {r: i + 1 for i, r in enumerate(sorted_unique_ranges)}
    print(f"[INFO] Mapping 'Range' to 'Block_Number': {range_to_block}")

    # Apply the mapping to create 'Mapped_Block_Number'
    df['Mapped_Block_Number'] = df['Range'].map(range_to_block)

    # Handle unmapped ranges (those that were 0 or any other unexpected values)
    if df['Mapped_Block_Number'].isnull().any():
        unmapped_ranges = df[df['Mapped_Block_Number'].isnull()]['Range'].unique()
        print(f"[WARNING] The following 'Range' values could not be mapped and will be excluded: {unmapped_ranges}")

    # Filter rows where 'Block_Number' matches 'Mapped_Block_Number'
    filtered_df = df[df['Block_Number'] == df['Mapped_Block_Number']].copy()
    print(f"[INFO] Filtered DataFrame shape: {filtered_df.shape}")

    # Drop the 'Block_Number' and 'Mapped_Block_Number' columns
    filtered_df = filtered_df.drop(columns=['Block_Number', 'Mapped_Block_Number'])
    print(f"[INFO] 'Block_Number' column removed from the filtered DataFrame.")

    return filtered_df


def aggregate_over_range(
    df: pd.DataFrame,
    grouping_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregates data over Range dimension.
    (From local_dprime_workflow.py)
    """
    if grouping_cols is None:
        grouping_cols = ["ExperimentName", "Regression", "Subject"]

    required = set(grouping_cols) | {"Range"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns for Range aggregation: {missing}")

    value_cols = [c for c in df.columns if c not in set(grouping_cols) | {"Range"}]
    if not value_cols:
        return df[grouping_cols].drop_duplicates().reset_index(drop=True)

    agg_spec = {
        col: ("mean" if pd.api.types.is_numeric_dtype(df[col]) else "first")
        for col in value_cols
    }

    aggregated = (
        df.groupby(grouping_cols, dropna=False)[value_cols]
        .agg(agg_spec)
        .reset_index()
    )
    return aggregated


def generate_summary_plot(df: pd.DataFrame, metric: str, out_file: Path) -> None:
    """
    Generates summary plots for metrics across Range.
    (Simplified from local_dprime_workflow.py)
    """
    if metric not in df.columns:
        return

    summary = (
        df.groupby(["ExperimentName", "Regression", "Range"], dropna=False)[metric]
        .mean()
        .reset_index()
        .sort_values(["ExperimentName", "Regression", "Range"])
    )

    if summary.empty:
        return

    plt.figure(figsize=(10, 6))
    for (experiment, regression), sub in summary.groupby(["ExperimentName", "Regression"], dropna=False):
        label = f"{experiment} | {regression}"
        plt.plot(sub["Range"], sub[metric], marker="o", label=label)

    plt.xlabel("Range")
    plt.ylabel(metric)
    plt.title(f"{metric} across Range")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


def run_workflow(data_file: str = DEFAULT_DATA_FILE) -> dict[str, pd.DataFrame]:
    """
    Main workflow - combines notebook logic with local I/O.
    """
    input_path = data_path(data_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input data file not found: {input_path}")

    if input_path.suffix.lower() == ".csv":
        df = pd.read_csv(input_path)
    elif input_path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    else:
        raise ValueError(f"Unsupported input type: {input_path.suffix}")

    print(f"[INFO] Root: {ROOT}")
    print(f"[INFO] Input: {input_path}")
    print(f"[INFO] Original DataFrame shape: {df.shape}")

    df = rename_duplicate_columns(df)

    # Step 1: Determine Number of Blocks Based on Unique Non-Zero Ranges
    unique_non_zero_ranges = df['Range'].unique()
    unique_non_zero_ranges = [r for r in unique_non_zero_ranges if r != 0]
    num_of_blocks = len(unique_non_zero_ranges)
    print(f"\n[INFO] Number of unique non-zero Ranges: {num_of_blocks}")

    if num_of_blocks == 0:
        print("[WARNING] No unique non-zero Ranges found. Defaulting to 1 block.")
        num_of_blocks = 1

    # Step 2: Process Hits
    hits_df = process_hits(df)

    # Step 3: Process FA
    fa_df = process_fa(df, num_of_blocks)

    # Step 4: Merge DataFrames
    combined_df = merge_dataframes(hits_df, fa_df)

    # Step 5: Compute Dprime and CR (⭐ USING NOTEBOOK LOGIC)
    combined_df = compute_dprime_cr(
        combined_df,
        hit_sum_col='ACC_sum_Hits',
        hit_count_col='ACC_count',
        fa_sum_col='ACC_sum_FA',
        fa_count_col='ACC_len'
    )

    # Step 6: Filter Blocks by Mapping
    filtered_df = filter_blocks_by_range(combined_df)

    # Step 7: Aggregate Over Range
    grouped_df = aggregate_over_range(filtered_df)

    # Step 8: Export to CSV
    print("\n[INFO] Exporting results to CSV...")
   # hits_df.to_csv(data_path("Hits_Output.csv"), index=False)
    #fa_df.to_csv(data_path("FA_Output.csv"), index=False)
    #combined_df.to_csv(data_path("Combined_Output_WithDprime.csv"), index=False)
    filtered_df.to_csv(data_path("dprime_results_with_range.csv"), index=False)
    grouped_df.to_csv(data_path("dprime_results.csv"), index=False)

    print(f"[INFO] Combined shape: {combined_df.shape}")
    print(f"[INFO] Filtered shape: {filtered_df.shape}")
    print(f"[INFO] Grouped shape: {grouped_df.shape}")

    # Step 9: Generate Summary Plots
    print("\n[INFO] Generating summary plots...")
    generate_summary_plot(filtered_df, "Dprime", plot_path("Dprime_summary.png"))
    generate_summary_plot(filtered_df, "CR", plot_path("CR_summary.png"))

    return {
        "raw": df,
        "hits": hits_df,
        "fa": fa_df,
        "combined": combined_df,
        "filtered": filtered_df,
        "grouped": grouped_df,
    }


if __name__ == "__main__":
    run_workflow()



