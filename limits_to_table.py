import ROOT
import os

# Define the MX and MY lists
MX = ["1400", "1600", "1800", "2000", "2200", "2600", "3000"]
MY = ["90", "125", "190", "250", "300", "400"]

# Normalization factor
norm_xsec = 5

# Latex-compatible sigma symbols
sigma_symbols = ["-2 &sigma;", "-1 &sigma;", "Expected", "+1 &sigma;", "+2 &sigma;", "Observed"]

# Function to extract data from a .root file
def extract_data(filepath):
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None
    file = ROOT.TFile.Open(filepath)
    tree = file.Get("limit")
    values = [entry.limit for entry in tree]
    file.Close()
    return values

# Function to extract significance
def extract_significance(mx, my):
    filepath = f"SR_run2/MX{mx}_MY{my}-2_area/higgsCombineTest.Significance.mH120.root"
    values = extract_data(filepath)
    return values[0] if values else None

# Function to format a TWiki table row
def format_twiki_row(mx, my, limits, significance, highlight=False):
    signal = f"MX{mx}_MY{my}"
    formatted_limits = [f"{v * norm_xsec:.2f}" for v in limits]
    formatted_significance = f"{significance:.1f}" if significance is not None else "N/A"
    row = f"| {signal} | " + " | ".join(formatted_limits) + f" | {formatted_significance} |"
    if highlight:
        row = f"%RED%{row}%ENDCOLOR%"
    return row

# Find the row with the maximum significance
significance_data = []
for mx in MX:
    for my in MY:
        significance = extract_significance(mx, my)
        if significance is not None:
            significance_data.append((mx, my, significance))

max_significance_entry = max(significance_data, key=lambda x: x[2]) if significance_data else None

# Create the TWiki table
print("| Signal | " + " | ".join(sigma_symbols) + " | Significance |")
for mx in MX:
    for my in MY:
        limits = extract_data(f"SR_run2/MX{mx}_MY{my}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root")
        significance = extract_significance(mx, my)
        highlight = max_significance_entry and (mx, my) == (max_significance_entry[0], max_significance_entry[1])
        if limits:
            print(format_twiki_row(mx, my, limits, significance, highlight))
