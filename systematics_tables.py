import ROOT
import os
import csv

def get_histogram(file, hist_name):
    """Function to retrieve a histogram from the ROOT file."""
    hist = file.Get(hist_name)
    if not hist:
        print(f"Warning: Histogram {hist_name} not found!")
    return hist

def calculate_yield(hist):
    """Function to calculate the integral (yield) of a histogram."""
    return hist.Integral() if hist else 0.0

def format_relative_change(change):
    """Formats the relative change percentage for LaTeX output."""
    if abs(change) < 0.1:
        return "<0.1\%"
    else:
        return f"{change:.2f}\%"

def calculate_vae_sf_variation(csv_file, years, processes, nominal_yields):
    """Calculates the weighted average of the vae_sf systematic variation."""
    total_weighted_up = 0.0
    total_weighted_down = 0.0
    total_nominal_yield = sum(nominal_yields.values())

    if total_nominal_yield == 0:
        return 0.0, 0.0

    with open(csv_file, mode='r') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            csv_year = row[0]
            csv_process = row[1]
            up_variation = float(row[3])
            down_variation = float(row[4])

            # Loop over years and processes, match with the CSV line
            for year in years:
                if csv_year == year and csv_process == processes[year]:
                    nominal_yield = nominal_yields[year]
                    total_weighted_up += nominal_yield * up_variation
                    total_weighted_down += nominal_yield * down_variation

    # Calculate the weighted averages
    average_up = total_weighted_up / total_nominal_yield
    average_down = total_weighted_down / total_nominal_yield

    return average_up * 100, average_down * 100  # Return in percentage

def main(base_folder, region, systematics, csv_file):
    # Define the years we are interested in
    years = ["2016APV", "2016", "2017", "2018"]
    processes = {year: f"MX2200_MY250" for year in years}  # Define process for each year

    # Initialize total yields for nominal and systematic variations
    total_yield_nominal = 0.0
    total_yields_up = {syst: 0.0 for syst in systematics}
    total_yields_down = {syst: 0.0 for syst in systematics}
    nominal_yields_by_year = {}

    # Loop over years and accumulate yields
    for year in years:
        process = f"MX2200_MY250"  # Define process for the current year
        rootfile_path = os.path.join(base_folder, f"templates_{process}_{year}.root")
        
        if not os.path.exists(rootfile_path):
            print(f"File {rootfile_path} does not exist!")
            continue

        rootfile = ROOT.TFile(rootfile_path, "READ")

        # Get the nominal histogram and sum its yield
        hist_nom_name = f"mjj_my_{process}_{region}_nom"
        hist_nom = get_histogram(rootfile, hist_nom_name)
        yield_nominal = calculate_yield(hist_nom)
        total_yield_nominal += yield_nominal
        nominal_yields_by_year[year] = yield_nominal

        # Loop over all systematics and accumulate their yields
        for syst in systematics:
            hist_up_name = f"mjj_my_{process}_{region}_{syst}_up"
            hist_down_name = f"mjj_my_{process}_{region}_{syst}_down"

            hist_up = get_histogram(rootfile, hist_up_name)
            hist_down = get_histogram(rootfile, hist_down_name)

            yield_up = calculate_yield(hist_up)
            yield_down = calculate_yield(hist_down)

            total_yields_up[syst] += yield_up
            total_yields_down[syst] += yield_down

        # Close the ROOT file after processing
        rootfile.Close()

    # LaTeX table header
    print("\\begin{table}[htbp]")
    print("\\centering")
    print("\\begin{tabular}{|c|c|c|}")
    print("\\hline")
    print("Syst. name & Up yield change & Down yield change \\\\")
    print("\\hline")

    # Print out the relative variations for each systematic
    syst_names = {"pdf":"pdf\_1\_signal", "prefire":"CMS\_l1\_ecal\_prefiring", "pileup":"CMS\_pileup\_13TeV", "PS_ISR":"ps\_isr", "PS_FSR":"ps\_fsr", "jes":"CMS\_scale\_j", "jer":"CMS\_res\_j", "pnet":"CMS\_Xbbtag\_PNet","vae_sf":"CMS\_anomaly\_tagging\_signal"}
    for syst in systematics:
        total_yield_up = total_yields_up[syst]
        total_yield_down = total_yields_down[syst]

        if total_yield_nominal != 0:
            rel_change_up = (total_yield_up - total_yield_nominal) / total_yield_nominal * 100
            rel_change_down = (total_yield_down - total_yield_nominal) / total_yield_nominal * 100
        else:
            rel_change_up = rel_change_down = 0.0

        # Format relative changes
        rel_change_up_str = format_relative_change(rel_change_up)
        rel_change_down_str = format_relative_change(rel_change_down)

        # Print LaTeX table row
        syst_name = syst_names[syst]
        print(f"{syst_name} & ${rel_change_up_str}$ & ${rel_change_down_str}$ \\\\")

    # Calculate and print the vae_sf systematic variation
    average_vae_up, average_vae_down = calculate_vae_sf_variation(csv_file, years, processes, nominal_yields_by_year)
    print(f"vae\_sf & ${format_relative_change(average_vae_up)}$ & ${format_relative_change(average_vae_down)}$ \\\\")
    print("\\hline")
    print("\label{tab:syst}")
    # LaTeX table footer
    print("\\end{tabular}")
    print("\\caption{The effect of systematic uncertainties on the yields in SR Pass for $\MX=2200 \GeV$ and $\MY=250 \GeV$ \XHYbbWW signal.}")
    print("\\end{table}")

if __name__ == "__main__":
    # Define the base folder and region
    base_folder = "/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v7/"
    region = "SR_Pass"
    csv_file = "SFs_vae.txt"  # Path to the CSV file containing vae_sf information

    # Define the list of systematics
    systematics = ["pdf", "prefire", "pileup", "PS_ISR", "PS_FSR","jes", "jer", "pnet", "vae_sf"]

    # Call the main function
    main(base_folder, region, systematics, csv_file)
