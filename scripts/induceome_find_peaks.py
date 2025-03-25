import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# Define function to parse mpileup file
def parse_mpileup(file_path):
    """
    Parse an mpileup file to extract position and coverage.

    Args:
        file_path (str): Path to the mpileup file.

    Returns:
        pd.DataFrame: Dataframe with 'Position' and 'Coverage'.
    """
    data = []
    with open(file_path, "r") as file:
        for line in file:
            cols = line.strip().split("\t")
            if len(cols) >= 4:  # Ensure line has sufficient columns
                position = int(cols[1])
                coverage = int(cols[3])
                data.append((position, coverage))

    print(f"Total positions: {len(data)}")

    return pd.DataFrame(data, columns=["Position", "Coverage"])


def infer_prominence(coverage_data, factor=20):
    """
    Infer the prominence threshold for peak detection.

    Args:
        coverage_data (np.array): Array of coverage values.

    Returns:
        int: Suggested prominence threshold.
    """
    avg_coverage = np.nanmean(coverage_data)
    max_coverage = np.nanmax(coverage_data)
    print(print(f"Average coverage: {avg_coverage}"))
    print(f"Max coverage: {max_coverage}")
    return max_coverage / factor
    return avg_coverage * factor


def find_peak_regions_sliding_window(df, threshold, window_size=500):
    """
    Identifies peak regions using a sliding window average.

    Args:
        df (pd.DataFrame): DataFrame with 'Position' and 'Coverage' columns.
        threshold (float): Minimum average coverage required for a peak.
        window_size (int): Size of the sliding window in base pairs.

    Returns:
        pd.DataFrame: DataFrame with columns ['Start', 'End', 'Average_Coverage'].
    """
    # Ensure the DataFrame is sorted by position
    # df = df.sort_values(by="Position").reset_index(drop=True)

    # List to store detected peaks
    peaks = []

    # Iterate with a sliding window
    for i in range(0, len(df) - window_size, window_size):
        avg_cov = df.loc[i : i + window_size, "Coverage"].mean()
        if avg_cov >= threshold:
            # Add the new peak or extend the previous one
            if not peaks:
                peaks.append(
                    (
                        df.loc[i, "Position"],
                        df.loc[i + window_size, "Position"],
                        avg_cov,
                    )
                )
            elif df.loc[i, "Position"] - peaks[-1][1] <= window_size:
                peaks[-1] = (peaks[-1][0], df.loc[i + window_size, "Position"], avg_cov)
            else:
                peaks.append(
                    (
                        df.loc[i, "Position"],
                        df.loc[i + window_size, "Position"],
                        avg_cov,
                    )
                )

    # Convert list of peaks to DataFrame
    peak_regions_df = pd.DataFrame(peaks, columns=["Start", "End", "Average_Coverage"])

    return peak_regions_df


def smooth_peaks(peaks, coverage, smoothing_range=1000):
    """
    Smooth the detected peak regions by merging nearby peaks.

    Args:
        peaks (pd.DataFrame): DataFrame with 'Start', 'End', 'Average_Coverage'
        coverage (pd.DataFrame): DataFrame with 'Position' and 'Coverage'
        smoothing_range (int): Maximum distance between peaks to merge

    Returns:
        pd.DataFrame: DataFrame with 'Start', 'End', 'Average_Coverage'
    """
    smoothed_peaks = []
    final_peaks = []
    for i, peak in peaks.iterrows():
        if not smoothed_peaks:
            smoothed_peaks.append(peak)
        elif peak["Start"] - smoothed_peaks[-1][1] <= smoothing_range:
            smoothed_peaks[-1] = (
                smoothed_peaks[-1][0],
                peak["End"],
                peak["Average_Coverage"],
            )
        else:
            smoothed_peaks.append(peak)

    # Recalculate average coverage for smoothed peaks
    for i, peak in enumerate(smoothed_peaks):
        start, end = peak[0], peak[1]
        avg_cov = coverage.loc[
            (coverage["Position"] >= start) & (coverage["Position"] <= end), "Coverage"
        ].mean()
        final_peaks.append((start, end, avg_cov))

    return pd.DataFrame(final_peaks, columns=["Start", "End", "Average_Coverage"])


def plot_peaks(df, peaks, threshold, output_image_fp):
    """
    Plot the coverage data with detected peak regions and threshold line.

    Args:
        df (pd.DataFrame): DataFrame with 'Position' and 'Coverage'
        peaks (pd.DataFrame): DataFrame with 'Start', 'End', 'Average_Coverage'
        threshold (float): Minimum coverage threshold for peaks
        output_image_fp (str): Output image file path
    """
    plt.figure(figsize=(12, 6))
    plt.plot(df["Position"], df["Coverage"], label="Coverage")
    plt.axhline(threshold, color="r", linestyle="--", label="Threshold")

    for _, peak in peaks.iterrows():
        plt.axvspan(peak["Start"], peak["End"], color="gray", alpha=0.5)

    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title("Coverage Plot with Peak Regions")
    plt.legend()
    plt.savefig(output_image_fp)
    plt.close()


# Example usage
df = parse_mpileup(snakemake.input.pileup)
if df.empty:
    print("No coverage reported.")
    Path(snakemake.output.peaks_img).touch()
    Path(snakemake.output.peaks_contigs).touch()
    with open(snakemake.output.peaks_csv, "w") as f:
        f.write("Start,End,Average_Coverage\n")
    exit(0)

threshold = infer_prominence(df["Coverage"].values)
print("Threshold: ", threshold)
peaks = find_peak_regions_sliding_window(
    df, threshold=threshold, window_size=snakemake.params.min_width
)
peaks = smooth_peaks(
    peaks, df, snakemake.params.smoothing_factor * snakemake.params.min_width
)
plot_peaks(df, peaks, threshold, snakemake.output.peaks_img)
# plot_peaks(df.head(10000), pd.DataFrame(), threshold, snakemake.output.peaks_img)

with open(snakemake.input.ref) as f:
    ref = ""
    for line in f.readlines():
        if not line.startswith(">"):
            ref += line.strip()

with open(snakemake.output.peaks_contigs, "w") as f:
    for _, peak in peaks.iterrows():
        start = int(peak["Start"])
        end = int(peak["End"])
        f.write(f"> {start} - {end}\n")
        f.write(f"{ref[start:end]}\n")

peaks.to_csv(snakemake.output.peaks_csv, index=False)
