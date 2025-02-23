import numpy as np
import matplotlib.pyplot as plt

# Sample data structures for RF components (can be extended or connected to external files)
antennas = [
    # VHF (30 MHz - 300 MHz)
    {"name": "FAKE VHF Yagi Antenna", "frequency_range": (30, 300), "gain": 7, "type": "Yagi-Uda"},
   # {"name": "UWB Deployable - RW", "frequency range": (30,2000), "gain": 15, "type": "Panel Array"},
    
    # UHF (300 MHz - 3 GHz)
    {"name": "FAKE UHF Log-Periodic Dipole Array", "frequency_range": (300, 1000), "gain": 10, "type": "Log-Periodic"},
    {"name": "FAKE UHF Patch Antenna", "frequency_range": (300, 1000), "gain": 8, "type": "Patch"},
    
    
    # L-Band (1 GHz - 2 GHz)
    {"name": "FAKE L-Band Helical Antenna", "frequency_range": (1000, 2000), "gain": 12, "type": "Helical"},
    {"name": "FALCCON - RW", "frequency_range": (1000, 2000), "gain": 21, "type": "Dipole Array"},
    
    # S-Band (2 GHz - 4 GHz)
    {"name": "FAKE S-Band Parabolic Reflector", "frequency_range": (2000, 4000), "gain": 15, "type": "Parabolic"},
    {"name": "AC-2000 - AAC", "frequency_range": (2000, 2300), "gain": 5.2, "type": "Parabolic Dish"},
    
    # C-Band (4 GHz - 8 GHz)
    {"name": "FAKE C-Band Horn Antenna", "frequency_range": (4000, 8000), "gain": 18, "type": "Horn"},
    
    # X-Band (8 GHz - 12 GHz)
    {"name": "FAKE X-Band Parabolic Dish", "frequency_range": (8000, 12000), "gain": 20, "type": "Parabolic Dish"},
    {"name": "X-Band Patch Antenna - Enduro", "frequency_range": (8025, 8400), "gain": 6, "type": "Patch"},
    {"name": "4x4 X-Band Patch Array - Enduro", "frequency_range": (8025, 8400), "gain": 16, "type": "Patch Array"},
    # Ku-Band (12 GHz - 18 GHz)
    {"name": "FAKE Ku-Band Microstrip Antenna", "frequency_range": (12000, 18000), "gain": 23, "type": "Microstrip"},
    
    # K-Band (18 GHz - 27 GHz)
    {"name": "FAKE K-Band Waveguide Antenna", "frequency_range": (18000, 27000), "gain": 25, "type": "Waveguide"},
    
    # Ka-Band (27 GHz - 40 GHz)
    {"name": "FAKE Ka-Band Horn Antenna", "frequency_range": (27000, 40000), "gain": 30, "type": "Horn"},
    
    # Above 40 GHz (EHF - Millimeter Wave Frequencies)
    {"name": "FAKE EHF Lens Antenna", "frequency_range": (40000, 60000), "gain": 35, "type": "Lens"},
    {"name": "FAKE Terahertz Horn Antenna", "frequency_range": (60000, 100000), "gain": 40, "type": "Horn"},
]


mixers = [
    {"name": "FAKE Mixer 1", "input_freq_range": (2000, 4000),  "output_freq": 0.5},
    {"name": "FAKE Mixer 2", "input_freq_range": (4000, 8000), "output_freq": 1.0},
]

downconverters = [
    {"name": "FAKE Downconverter 1", "input_freq_range": (0.5, 2.0), "IF_freq": 0.1},
    {"name": "FAKE Downconverter 2", "input_freq_range": (1.0, 4.0), "IF_freq": 0.2},
]

modems_sdrs = [
    {"name": "FAKE Modem A", "bandwidth": 20, "modulation": "QPSK"},
    {"name": "FAKE SDR X", "bandwidth": 50, "modulation": "QAM"},
    {"name": "RX-2000 - AAC", "bandwidth": 20, "modulation": "GFSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulation": "OOK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulation": "GMSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulation": "2FSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulation": "4FSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulation": "2GFSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulation": "4GFSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 10, "modulation": "OOK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 10, "modulation": "GMSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 10, "modulation": "2FSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 10, "modulation": "4FSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 10, "modulation": "2GFSK"},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 10, "modulation": "4GFSK"},
    {"name": "TRX-U - Clyde", "bandwidth": 100, "modulation": "FSK"},
    {"name": "TRX-U - Clyde", "bandwidth": 100, "modulation": "GFSK"},
    
]

#s band uplink, s band downlink  
ground_segment = [
    {"name": "VIASAT PENDER", "sup_freq": (2025,2110), "uEIRP": 53.2, "sdown_bw": (2200, 2290), "sdown_gt": 17, "xdown_bw":(8025,8400), "xdown_gt": 30, "kadown_bw": None, "kadown_gt": None},
    {"name": "VIASAT GUILDFORD", "sup_freq": (2025,2110), "uEIRP": 53.2, "sdown_bw": (2200, 2290), "sdown_gt": 17, "xdown_bw":(8025,8400), "xdown_gt": 30, "kadown_bw": None, "kadown_gt": None},
    {"name": "VIASAT ALICE", "sup_freq": (2025,2110), "uEIRP": 65.0, "sdown_bw": (2200, 2290), "sdown_gt": 18, "xdown_bw":(8025,8400), "xdown_gt": 32, "kadown_bw":(25500,27000), "kadown_gt":34.5},
    {"name": "VIASAT GHANA", "sup_freq": (2025,2110), "uEIRP": 65.0, "sdown_bw": (2200, 2290), "sdown_gt": 18, "xdown_bw":(8025,8400), "xdown_gt": 32, "kadown_bw":(25500,27000), "kadown_gt":34.5},  
    {"name": "ATLAS PAUMALI", "sup_freq": (2025,2120), "uEIRP":50.0 , "sdown_bw": (2200, 2300), "sdown_gt": 21, "xdown_bw":(7900,8500), "xdown_gt": 31, "kadown_bw": None, "kadown_gt": None},  
]

cables = [
    {"name": "FAKE Cable A", "loss_per_meter": 0.1, "max_length": 10},
    {"name": "FAKE Cable B", "loss_per_meter": 0.05, "max_length": 20},
]



# Orbital parameters and ground station properties
orbital_parameters = {
    "altitude": 2000,  # km
    "frequency": 2050, # MHz
    "inclination": 45,  # degrees
    "data_rate": 100,  # Mbps
    "P_noise_Earth": (1.38*10**-23)*290*3500000000, # =(k)(T)(B)
    "T_Earth": 290 # avg T in K 
    
# P_rx = P_tx + G_tx + G_rx - L_path - L_other 
# SNR = P_rx /  P_noise
# C = B * log2(1 + SNR)
}

# Placeholder function for link budget calculation
def calculate_link_budget(orbital_params, antenna_gain, cable_loss, modulation_type):
    path_loss = 20 * np.log10(orbital_params["altitude"]) + 20 * np.log10(orbital_params["frequency"]) - 147.56
    eirp = antenna_gain - cable_loss   # Simplified EIRP calculation
    received_power = eirp - path_loss  # Simplified received power calculation
    
    return received_power

# Function to select RF components based on orbital parameters
def select_rf_components(orbital_params):
    chosen_components = {}

    # Antenna selection based on frequency
    for antenna in antennas:
        if antenna["frequency_range"][0] <= orbital_params["frequency"] <= antenna["frequency_range"][1]:
            chosen_components["Antenna"] = antenna
            break

    # Mixer selection based on frequency
    for mixer in mixers:
        if mixer["input_freq_range"][0] <= orbital_params["frequency"] <= mixer["input_freq_range"][1]:
            chosen_components["Mixer"] = mixer
            break

    # Downconverter selection
    for downconverter in downconverters:
        if "Mixer" in chosen_components and downconverter["input_freq_range"][0] <= chosen_components["Mixer"]["output_freq"] <= downconverter["input_freq_range"][1]:
            chosen_components["Downconverter"] = downconverter
            break

    # Modem/SDR selection based on data rate
    for modem_sdr in modems_sdrs:
        
        #if modem_sdr["bandwidth"]*np.log2(1+SNR)
        chosen_components["Modem/SDR"] = modem_sdr  # Basic selection logic
        
        break

    # Cable selection (basic choice)
    chosen_components["Cable"] = cables[0]  # Placeholder choice
    return chosen_components

# Generate the configuration and perform link budget analysis
chosen_components = select_rf_components(orbital_parameters)
antenna_gain = chosen_components["Antenna"]["gain"]
cable_loss = chosen_components["Cable"]["loss_per_meter"] * chosen_components["Cable"]["max_length"]

received_power = calculate_link_budget(orbital_parameters, antenna_gain, cable_loss, chosen_components["Modem/SDR"]["modulation"])

print("\nChosen RF System Components:")
for component_type, details in chosen_components.items():
    print(f"{component_type}: {details['name']} - Specs: {details}")

print(f"\nCalculated Received Power: {received_power:.2f} dBm")

# Visualization of radiation patterns
theta = np.linspace(0, 2 * np.pi, 360)
antenna_pattern = antenna_gain * np.cos(theta)**2
ground_station_pattern = (antenna_gain - np.abs(theta - np.pi / 2)) * 1.5  # Simplified pattern

plt.figure(figsize=(10, 5))
plt.polar(theta, antenna_pattern, label="Satellite Antenna Pattern")
plt.polar(theta, ground_station_pattern, label="Ground Station Antenna Pattern")
plt.title("Radiation Pattern Overlap")
plt.legend()
plt.show()
