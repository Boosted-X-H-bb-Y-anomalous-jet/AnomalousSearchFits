import json

def calculate_average_percentage_difference(json_file_path):
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    results = []
    for param in data['params']:
        r_values = param['r']
        min_val, best_val, max_val = r_values[0], r_values[1], r_values[2]
        diff_min = abs(min_val - best_val)
        diff_max = abs(max_val - best_val)
        perc_diff_min = (diff_min / best_val) * 100
        perc_diff_max = (diff_max / best_val) * 100
        avg_percentage_diff = (perc_diff_min + perc_diff_max) / 2
        results.append({
            "name": param["name"],
            "average_percentage_difference": avg_percentage_diff
        })
    
    # Sort the results by average percentage difference in descending order
    results.sort(key=lambda x: x['average_percentage_difference'], reverse=True)
    
    return results

json_file_path = "SR_run2_one_signal/MX2200_MY250-0_area/impacts.json"
averages = calculate_average_percentage_difference(json_file_path)

for result in averages:
    print(f"Uncertainty: {result['name']}, Average Percentage Difference: {result['average_percentage_difference']:.2f}%")
