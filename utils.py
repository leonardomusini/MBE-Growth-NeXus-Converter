import os
import re
from datetime import datetime, timedelta
from fractions import Fraction
import numpy as np


def file_searcher(folder_path):
    """
    Extracts the sample ID part of filenames for all .wri files in a given folder.

    :param folder_path: Path to the folder containing .wri files.
    :return: List of sample IDs (as integers).
    """
    extracted_ids = []

    if not os.path.isdir(folder_path):
        print("Error: The provided folder path does not exist.")
        return []

    # Iterate over files in the folder
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(".wri"):  
            sample_id = filename[2:-4] 
            extracted_ids.append(sample_id)

    return extracted_ids

# ----------------------------------

def time_calculator(time_str, time_delta_str):

    time_obj = datetime.strptime(time_str, "%H:%M:%S")
    h, m, s = map(int, time_delta_str.split(":"))  
    time_delta = timedelta(hours=h, minutes=m, seconds=s) 

    new_time_obj = time_obj + time_delta
    new_time_str = new_time_obj.strftime("%H:%M:%S")

    return new_time_str

# ----------------------------------

def area_converter(area: str) -> float:
    """
    Converts an area value given as a fraction or integer string into mm².

    Parameters:
        area (str): A string representing the area fraction (e.g., "1/4", "1/2", "1").

    Returns:
        float: The area converted to mm². If area is "piece" return 0.0
    """
    base_area = 2000  # "1" = 2000 mm²

    if area.lower() == "piece":
        return 0.0

    try:
        numeric_value = float(Fraction(area))  
        return numeric_value * base_area
    except ValueError:
        raise ValueError(f"Invalid area format: {area}. Expected fraction or integer as a string.")
    
# ----------------------------------

def time_converter(date_str, time_str):
    """
    Converts date and time strings into an ISO 8601 timestamp (YYYY-MM-DDTHH:MM:SS).
    
    Parameters:
        date_str (str): The date string in MM-DD-YYYY format.
        time_str (str): The time string in HH:MM:SS format.

    Returns:
        str: ISO 8601 formatted timestamp (YYYY-MM-DDTHH:MM:SS).
    """
    try:
        # MM-DD-YYYY
        parsed_date = datetime.strptime(date_str, "%m-%d-%Y")  
        # HH:MM:SS
        parsed_time = datetime.strptime(time_str, "%H:%M:%S").time()

        combined_datetime = datetime.combine(parsed_date.date(), parsed_time)

        return combined_datetime.isoformat()  # Convert to ISO 8601 format

    except ValueError as e:
        raise ValueError(f"Error parsing date or time: {e}")
    
# ----------------------------------

def arsenic_ranges(value):

    # Typo in scientific notation
    value = value.strip().replace(".E", "E")

    # Pattern to detect ranges"
    match = re.match(r"^(\d+(\.\d+)?)-(\d+(\.\d+)?E[-+]?\d+)$", value)
    
    if match:
        first_match = float(match.group(1))  
        second_match = float(match.group(3))  
        
        # Extract the exponent from the second number 
        exponent_match = re.search(r"E([-+]?\d+)", match.group(3))
        if exponent_match:
            exponent = int(exponent_match.group(1))  
            first_match = first_match * (10 ** exponent)  
        
        return (first_match + second_match) / 2 
    
    # If it's just a normal float, parse it directly
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"Invalid number format: {value}")
    
# ----------------------------------

def doping_calculator(dop_element, Si_curr, C_curr, g_rate):
    
    #gr_ratio = 2.8 / g_rate

    if dop_element is None:
        return 0.0
    
    elif dop_element == "C":
        if C_curr is None or np.isnan(C_curr):
            return 0.0
        return 9.8653e8 * np.exp(0.41434 * C_curr) * 2.8 / g_rate
    
    elif dop_element == "Si":
        if Si_curr is None or np.isnan(Si_curr):
            return 0.0
        elif Si_curr < 16:
            return 1.0223e7 * np.exp(1.5319 * Si_curr) * 2.8 / g_rate
        else:
            return 3.0471e12 * np.exp(0.68786 * Si_curr) * 2.8 / g_rate
    
    else:
        raise ValueError(f"Unknown dopant: {dop_element}")



# Extra utils functions

def file_searcher_old(folder_path):
    """
    Extracts the numerical part of filenames for all .wri files in a given folder.

    :param folder_path: Path to the folder containing .wri files.
    :return: List of extracted numbers (as integers).
    """
    extracted_numbers = []

    if not os.path.isdir(folder_path):
        print("Error: The provided folder path does not exist.")
        return []

    # Regex pattern to capture numbers in the filename before .wri
    pattern = re.compile(r"(\d+)", re.IGNORECASE)

    # Iterate over files in the folder
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(".wri"):  
            match = pattern.search(filename)
            if match:
                extracted_numbers.append(int(match.group(1))) 

    return extracted_numbers

   
