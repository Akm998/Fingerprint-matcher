import subprocess

def run_bozorth3_with_verbose(probe_file, gallery_file, output_file="comparison_table.txt"):
    """
    Run Bozorth3 algorithm with verbose output to capture comparison table.

    Parameters:
    probe_file (str): Path to the probe .xyt file
    gallery_file (str): Path to the gallery .xyt file
    output_file (str): File to store the comparison table output

    Returns:
    None
    """
    command = ["bozorth3", "-v", probe_file, gallery_file]
    with open(output_file, 'w') as f:
        result = subprocess.run(command, stdout=f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        raise Exception(f"Bozorth3 failed with error: {result.stderr}")

# Example usage
probe_file = "minutiae_file1.xyt"
gallery_file = "minutiae_file2.xyt"
run_bozorth3_with_verbose(probe_file, gallery_file)
print("Comparison table saved to 'comparison_table.txt'")
