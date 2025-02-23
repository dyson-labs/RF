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
    {"name": "K-Band 4x4 Patch Array - Enduro", "frequency_range": (17700, 20200), "gain": 16, "type": "Waveguide"},
    
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

modulation_bits_per_symbol = {
    "OOK": 1,
    "BPSK": 1,
    "QPSK": 2,
    "GMSK": 1,  # Gaussian Minimum Shift Keying (similar to BPSK)
    "2FSK": 1,
    "4FSK": 2,
    "2GFSK": 1,
    "4GFSK": 2,
    "8PSK": 3,
    "16QAM": 4,
    "QAM": 4,  # Assuming 16-QAM for this example
    "GFSK": 1,
    "FSK": 1,
}
    

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

# Orbital and satellite parameters
orbital_parameters = {
    "altitude": 200,  # km
    "uplink_frequency": 2025,  # MHz
    "downlink_frequency": 25600,  # MHz
    "T_system": 290,  # System noise temperature in K
}

satellite_parameters = {
    "G_satellite_receiver_dBi": 10,  # Satellite receive antenna gain
    "G_satellite_transmit_dBi": 10,  # Satellite transmit antenna gain
    "P_satellite_transmit_dBm": 30,  # Satellite transmit power in dBm (1 W)
}

# Function to calculate Free Space Path Loss (FSPL)
def calculate_fspl(altitude, frequency):
    """Calculate Free Space Path Loss (FSPL) in dB."""
    # Convert altitude to meters and frequency to Hz
    d = altitude * 1e3  # Altitude in meters
    f = frequency * 1e6  # Frequency in Hz
    c = 3e8  # Speed of light in m/s
    fspl = 20 * np.log10(4 * np.pi * d * f / c)
    return fspl

# Function to calculate link budget
def calculate_link_budget(orbital_params, ground_station, satellite_params, modem, modulation_bits_per_symbol):
    """
    Calculate link budget for both uplink and downlink.
    
    Returns:
        P_r_uplink_dBm (float): Received power at satellite (uplink) in dBm.
        P_r_downlink_dBm (float): Received power at ground station (downlink) in dBm.
        SNR_dB (float): Signal-to-Noise Ratio in dB.
        capacity_shannon_Mbps (float): Shannon capacity in Mbps.
        data_rate_Mbps (float): Achievable data rate based on modulation in Mbps.
    """
    # Extract bandwidth and modulation
    B = modem["bandwidth"]  # Bandwidth in MHz
    B = B * 1000 # Bandwidth in Hz
    modulation = modem["modulation"]
    bits_per_symbol = modulation_bits_per_symbol.get(modulation, 1)
    
    # Uplink FSPL
    uplink_fspl = calculate_fspl(orbital_params["altitude"], orbital_params["uplink_frequency"])
    
    # Uplink calculation (Ground to Satellite)
    EIRP_ground_dBm = ground_station["uEIRP"]  # Ground station EIRP in dBm
    G_satellite_rx_dBi = satellite_parameters["G_satellite_receiver_dBi"]
    L_uplink_dB = 0  # Assume no additional losses
    P_r_uplink_dBm = EIRP_ground_dBm - uplink_fspl + G_satellite_rx_dBi - L_uplink_dB
    
    # Downlink FSPL
    downlink_fspl = calculate_fspl(orbital_params["altitude"], orbital_params["downlink_frequency"])
    
    # Downlink calculation (Satellite to Ground)
    P_satellite_tx_dBm = satellite_parameters["P_satellite_transmit_dBm"]
    G_satellite_tx_dBi = satellite_parameters["G_satellite_transmit_dBi"]
    EIRP_satellite_dBm = P_satellite_tx_dBm + G_satellite_tx_dBi
    G_ground_rx_dBi = ground_station["sdown_gt"]  # Ground station receive antenna gain
    L_downlink_dB = 0  # Assume no additional losses
    P_r_downlink_dBm = EIRP_satellite_dBm - downlink_fspl + G_ground_rx_dBi - L_downlink_dB
    
    # Noise power calculation
    k = 1.38e-23  # Boltzmann's constant in J/K
    T_system = orbital_params["T_system"]  # System noise temperature in Kelvin
    P_noise_W = k * T_system * B  # Noise power in Watts
    P_noise_dBm = 10 * np.log10(P_noise_W) + 30  # Convert to dBm
    
    # SNR calculation in dB
    SNR_dB = P_r_downlink_dBm - P_noise_dBm
    SNR_linear = 10**(SNR_dB / 10)
    
    # Theoretical capacity using Shannon-Hartley theorem
    capacity_shannon_bps = B * np.log2(1 + SNR_linear)
    capacity_shannon_Mbps = capacity_shannon_bps / 1e6  # Convert to Mbps
    
    # Practical achievable data rate based on modulation scheme
    symbol_rate = B  # Assuming Nyquist rate (symbol rate equals bandwidth)
    data_rate_bps = symbol_rate * bits_per_symbol
    data_rate_Mbps = data_rate_bps / 1e6  # Convert to Mbps
    
    # Adjust data rate based on SNR
    max_bits_per_symbol = np.log2(1 + SNR_linear)
    if bits_per_symbol > max_bits_per_symbol:
        print(f"Warning: The modulation scheme '{modulation}' requires higher SNR.")
        data_rate_bps = symbol_rate * max_bits_per_symbol  # Adjust data rate
        data_rate_Mbps = data_rate_bps / 1e6
    
    return P_r_uplink_dBm, P_r_downlink_dBm, SNR_dB, capacity_shannon_Mbps, data_rate_Mbps

# Select ground station
def select_ground_station(orbital_params, ground_segment):
    for station in ground_segment:
        if station["sup_freq"][0] <= orbital_params["uplink_frequency"] <= station["sup_freq"][1]:
            return station
    return None

# Select appropriate modem/SDR
def select_modem(orbital_params, modems_sdrs):
    for modem in modems_sdrs:
        # Check if modem's bandwidth and modulation are acceptable
        return modem  # For simplicity, just return the first modem
    return None

# Run calculations
ground_station = select_ground_station(orbital_parameters, ground_segment)
modem = select_modem(orbital_parameters, modems_sdrs)

if ground_station is None:
    print("No suitable ground station found for the given uplink frequency.")
elif modem is None:
    print("No suitable modem/SDR found for the given parameters.")
else:
    uplink_power_dBm, downlink_power_dBm, SNR_dB, capacity_shannon_Mbps, data_rate_Mbps = calculate_link_budget(
        orbital_parameters, ground_station, satellite_parameters, modem, modulation_bits_per_symbol
    )
    
    print(f"Selected Ground Station: {ground_station['name']}")
    print(f"Selected Modem/SDR: {modem['name']} with Modulation {modem['modulation']}")
    print(f"Uplink Received Power at Satellite: {uplink_power_dBm:.2f} dBm")
    print(f"Downlink Received Power at Ground Station: {downlink_power_dBm:.2f} dBm")
    print(f"SNR: {SNR_dB:.2f} dB")
    print(f"Shannon Capacity: {capacity_shannon_Mbps:.9f} Mbps")
    print(f"Achievable Data Rate based on Modulation: {data_rate_Mbps:.9f} Mbps")
