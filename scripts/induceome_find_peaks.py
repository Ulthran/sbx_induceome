import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths


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
    return pd.DataFrame(data, columns=["Position", "Coverage"])


def infer_prominence(coverage_data, factor=5):
    """
    Infer the best prominence for peak detection based on assumptions:
    1. Few peaks of interest.
    2. Peaks rise substantially above the background (factor x baseline).
    3. Peaks have some thickness (measured by width).

    Args:
        coverage_data (np.array): Array of coverage values.
        factor (int): Factor by which a peak must exceed the baseline.

    Returns:
        float: Suggested prominence value.
    """
    # Use the 75th percentile of non-zero values as the baseline
    non_zero_coverage = coverage_data[coverage_data > 0]
    if len(non_zero_coverage) == 0:
        return 1  # Don't error out on empty coverage data

    baseline = np.percentile(non_zero_coverage, 75)  # 75th percentile
    suggested_prominence = factor * baseline
    print(f"Baseline (75th percentile of non-zero values): {baseline}")
    print(f"Suggested prominence (factor {factor}): {suggested_prominence}")

    return suggested_prominence


# Detect peaks
def detect_peaks(coverage_data, prominence, distance, min_width):
    """
    Detect peaks in coverage data.

    Args:
        coverage_data (np.array): Array of coverage values.
        prominence (int): Minimum prominence of peaks.
        distance (int): Minimum distance between peaks.

    Returns:
        list: Indices of detected peaks.
    """
    # Find peaks in the coverage data
    peaks, _ = find_peaks(coverage_data, prominence=prominence, distance=distance)

    # Calculate peak widths at half prominence height
    widths = peak_widths(coverage_data, peaks, rel_height=0.5)[0]

    # Filter peaks based on minimum width
    valid_peaks = [p for p, w in zip(peaks, widths) if w >= min_width]
    print(f"Total peaks detected: {len(peaks)}")
    print(f"Peaks after width filter (≥ {min_width}): {len(valid_peaks)}")

    return valid_peaks


# Main workflow
def extract_peaks(file_path, output_image_fp, distance, min_width):
    # Parse the mpileup file
    df = parse_mpileup(file_path)

    # Plot plain coverage
    plt.figure(figsize=(12, 6))
    plt.plot(df["Position"], df["Coverage"], label="Coverage")
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title("Plain Coverage Plot")
    plt.legend()
    # plt.show()

    # Detect peaks in the coverage data
    cov_values = df["Coverage"].values
    peaks = detect_peaks(
        cov_values,
        prominence=infer_prominence(cov_values),
        distance=distance,
        min_width=min_width,
    )

    # Extract peak regions
    peak_positions = df.iloc[peaks]

    # Plot the results for visualization
    plt.figure(figsize=(10, 6))
    plt.plot(df["Position"], df["Coverage"], label="Coverage")
    plt.scatter(
        peak_positions["Position"],
        peak_positions["Coverage"],
        color="red",
        label="Peaks",
    )
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title("Coverage Peaks")
    plt.legend()
    plt.savefig(output_image_fp, dpi=300, bbox_inches="tight")
    # plt.show()

    return peak_positions


# Example usage
peak_data = extract_peaks(
    snakemake.input.pileup,
    snakemake.output.peaks_img,
    distance=snakemake.params.distance,
    min_width=snakemake.params.min_width,
)
peak_data.to_csv(snakemake.output.peaks_csv, index=False)
