import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc, erfcinv
from scipy.optimize import bisect

# Antennas
antennas = [
    # VHF (30 MHz - 300 MHz)
    {"name": "FAKE VHF Yagi Antenna", "frequency_range": (30, 300), "gain": 7, "type": "Yagi-Uda"},

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
    "GFSK": 1,
    "FSK": 1,
}

modems_sdrs = [
    {"name": "FAKE Modem A", "bandwidth": 20, "modulations": ["BPSK", "QPSK", "8PSK", "16QAM"]},
    {"name": "FAKE SDR X", "bandwidth": 50, "modulations": ["BPSK", "QPSK", "16QAM"]},
    {"name": "RX-2000 - AAC", "bandwidth": 20, "modulations": ["BPSK", "GFSK"]},
    {"name": "UHF Transceiver II - Enduro", "bandwidth": 3, "modulations": ["OOK", "GMSK", "2FSK", "4FSK", "2GFSK", "4GFSK"]},
    {"name": "TRX-U - Clyde", "bandwidth": 0.048, "modulations": ["FSK", "GFSK"]},
    # (Add more modems if necessary)
]

mixers = [
    {"name": "FAKE Mixer 1", "input_freq_range": (2000, 4000),  "output_freq": 0.5},
    {"name": "FAKE Mixer 2", "input_freq_range": (4000, 8000), "output_freq": 1.0},
]

downconverters = [
    {"name": "FAKE Downconverter 1", "input_freq_range": (0.5, 2.0), "IF_freq": 0.1},
    {"name": "FAKE Downconverter 2", "input_freq_range": (1.0, 4.0), "IF_freq": 0.2},
]

cables = [
    {"name": "FAKE Cable A", "loss_per_meter": 0.1, "max_length": 10},
    {"name": "FAKE Cable B", "loss_per_meter": 0.05, "max_length": 20},
]

# Ground stations
ground_segment = [
    {"name": "VIASAT PENDER", "sup_freq": (2025, 2110), "uEIRP": 53.2, "sdown_bw": (2200, 2290), "sdown_gt": 17, "xdown_bw": (8025, 8400), "xdown_gt": 30, "kadown_bw": None, "kadown_gt": None},
    {"name": "VIASAT GUILDFORD", "sup_freq": (2025, 2110), "uEIRP": 53.2, "sdown_bw": (2200, 2290), "sdown_gt": 17, "xdown_bw": (8025, 8400), "xdown_gt": 30, "kadown_bw": None, "kadown_gt": None},
    {"name": "VIASAT ALICE", "sup_freq": (2025, 2110), "uEIRP": 65.0, "sdown_bw": (2200, 2290), "sdown_gt": 18, "xdown_bw": (8025, 8400), "xdown_gt": 32, "kadown_bw": (25500, 27000), "kadown_gt": 34.5},
    {"name": "VIASAT GHANA", "sup_freq": (2025, 2110), "uEIRP": 65.0, "sdown_bw": (2200, 2290), "sdown_gt": 18, "xdown_bw": (8025, 8400), "xdown_gt": 32, "kadown_bw": (25500, 27000), "kadown_gt": 34.5},
    {"name": "ATLAS PAUMALU", "sup_freq": (2025, 2120), "uEIRP": 50.0, "sdown_bw": (2200, 2300), "sdown_gt": 21, "xdown_bw": (7900, 8500), "xdown_gt": 31, "kadown_bw": None, "kadown_gt": None},
]

# Orbital and satellite parameters
orbital_parameters = {
    "altitude": 2000 ,  # km
    "uplink_frequency": 2025,  # MHz
    "downlink_frequency": 25600,  # MHz
    "T_system": 290,  # System temperature in K
}

satellite_parameters = {
    "G_satellite_receiver_dBi": None,  # To be selected from antennas
    "G_satellite_transmit_dBi": None,  # To be selected from antennas
    "P_satellite_transmit_W": 5,  # Satellite transmit power in Watts
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

# Q-function
def Q(x):
    return 0.5 * erfc(x / np.sqrt(2))

# Function to get required Eb/N0 for a given BER and modulation
def get_required_Eb_N0_dB(modulation_scheme, BER_target):
    if modulation_scheme in ["BPSK", "QPSK", "GMSK", "GFSK", "FSK", "OOK"]:
        # Use BPSK formula
        Eb_N0_linear = (erfcinv(2 * BER_target))**2
        Eb_N0_dB = 10 * np.log10(Eb_N0_linear)
        return Eb_N0_dB
    elif modulation_scheme == "8PSK":
        M = 8
        sin_pi_M = np.sin(np.pi / M)
        def func(Eb_N0_dB):
            Eb_N0_linear = 10 ** (Eb_N0_dB / 10)
            Q_arg = np.sqrt(2 * Eb_N0_linear) * sin_pi_M
            BER = 2 * Q(Q_arg)
            return BER - BER_target
        Eb_N0_dB = bisect(func, 0, 30)
        return Eb_N0_dB
    elif modulation_scheme == "16QAM":
        M = 16
        k = np.log2(M)
        def func(Eb_N0_dB):
            Eb_N0_linear = 10 ** (Eb_N0_dB / 10)
            Q_arg = np.sqrt(3 * Eb_N0_linear * k / (M - 1))
            BER = (4 / k) * (1 - 1 / np.sqrt(M)) * Q(Q_arg)
            return BER - BER_target
        Eb_N0_dB = bisect(func, 0, 30)
        return Eb_N0_dB
    else:
        # For other modulations, return a high required Eb_N0_dB
        return 30

# Function to calculate link budget
def calculate_link_budget(orbital_params, ground_station, ground_station_receive_gain_key, satellite_params, modem, modulation_bits_per_symbol, satellite_tx_antenna, satellite_rx_antenna, required_Eb_N0_dB_dict):
    """
    Calculate link budget for both uplink and downlink.

    Returns a dictionary with all calculated parameters.
    """
    # Extract bandwidth
    B = modem["bandwidth"]  # Bandwidth in MHz
    B = B * 1e6  # Bandwidth in Hz

    # Boltzmann's constant
    k = 1.38e-23  # Boltzmann's constant in J/K
    T_system = orbital_params["T_system"]  # System noise temperature in Kelvin

    # Uplink FSPL
    uplink_fspl = calculate_fspl(orbital_params["altitude"], orbital_params["uplink_frequency"])

    # Uplink calculation (Ground to Satellite)
    EIRP_ground_dBm = ground_station["uEIRP"]  # Ground station EIRP in dBm
    G_satellite_rx_dBi = satellite_rx_antenna["gain"]  # Satellite receive antenna gain
    L_uplink_dB = 0  # Assume no additional losses
    P_r_uplink_dBm = EIRP_ground_dBm - uplink_fspl + G_satellite_rx_dBi - L_uplink_dB

    # Noise power calculation for uplink
    P_noise_uplink_W = k * T_system * B  # Noise power in Watts
    P_noise_uplink_dBm = 10 * np.log10(P_noise_uplink_W) + 30  # Convert to dBm

    # Uplink SNR calculation in dB
    SNR_uplink_dB = P_r_uplink_dBm - P_noise_uplink_dBm
    SNR_uplink_linear = 10**(SNR_uplink_dB / 10)

    # Downlink FSPL
    downlink_fspl = calculate_fspl(orbital_params["altitude"], orbital_params["downlink_frequency"])

    # Downlink calculation (Satellite to Ground)
    P_satellite_tx_W = satellite_params["P_satellite_transmit_W"]
    P_satellite_tx_dBm = 10 * np.log10(P_satellite_tx_W * 1000)  # Convert W to dBm
    G_satellite_tx_dBi = satellite_tx_antenna["gain"]
    EIRP_satellite_dBm = P_satellite_tx_dBm + G_satellite_tx_dBi
    G_ground_rx_dBi = ground_station[ground_station_receive_gain_key]  # Ground station receive antenna gain
    L_downlink_dB = 0  # Assume no additional losses
    P_r_downlink_dBm = EIRP_satellite_dBm - downlink_fspl + G_ground_rx_dBi - L_downlink_dB

    # Noise power calculation for downlink
    P_noise_downlink_W = k * T_system * B  # Noise power in Watts
    P_noise_downlink_dBm = 10 * np.log10(P_noise_downlink_W) + 30  # Convert to dBm

    # Downlink SNR calculation in dB
    SNR_downlink_dB = P_r_downlink_dBm - P_noise_downlink_dBm
    SNR_downlink_linear = 10**(SNR_downlink_dB / 10)

    # Initialize best modulation schemes
    best_uplink_modulation = None
    best_downlink_modulation = None
    max_bits_per_symbol_uplink = 0
    max_bits_per_symbol_downlink = 0

    # Uplink Modulation Selection
    for modulation in modem["modulations"]:
        bits_per_symbol = modulation_bits_per_symbol[modulation]
        Eb_N0_uplink_dB = SNR_uplink_dB - 10 * np.log10(bits_per_symbol)
        required_Eb_N0_dB = required_Eb_N0_dB_dict[modulation]
        if Eb_N0_uplink_dB >= required_Eb_N0_dB and bits_per_symbol > max_bits_per_symbol_uplink:
            max_bits_per_symbol_uplink = bits_per_symbol
            best_uplink_modulation = modulation

    # Downlink Modulation Selection
    for modulation in modem["modulations"]:
        bits_per_symbol = modulation_bits_per_symbol[modulation]
        Eb_N0_downlink_dB = SNR_downlink_dB - 10 * np.log10(bits_per_symbol)
        required_Eb_N0_dB = required_Eb_N0_dB_dict[modulation]
        if Eb_N0_downlink_dB >= required_Eb_N0_dB and bits_per_symbol > max_bits_per_symbol_downlink:
            max_bits_per_symbol_downlink = bits_per_symbol
            best_downlink_modulation = modulation

    # If no modulation scheme meets the BER requirement, return None
    if best_uplink_modulation is None or best_downlink_modulation is None:
        return None

    # Data rate calculations
    symbol_rate = B  # Assuming Nyquist rate (symbol rate equals bandwidth)
    data_rate_uplink_bps = symbol_rate * max_bits_per_symbol_uplink
    data_rate_uplink_Mbps = data_rate_uplink_bps / 1e6  # Convert to Mbps

    data_rate_downlink_bps = symbol_rate * max_bits_per_symbol_downlink
    data_rate_downlink_Mbps = data_rate_downlink_bps / 1e6  # Convert to Mbps

    return {
        "uplink_power_dBm": P_r_uplink_dBm,
        "downlink_power_dBm": P_r_downlink_dBm,
        "SNR_uplink_dB": SNR_uplink_dB,
        "SNR_downlink_dB": SNR_downlink_dB,
        "data_rate_uplink_Mbps": data_rate_uplink_Mbps,
        "data_rate_downlink_Mbps": data_rate_downlink_Mbps,
        "best_uplink_modulation": best_uplink_modulation,
        "best_downlink_modulation": best_downlink_modulation,
        "modem": modem,
        "ground_station": ground_station,
        "ground_station_receive_gain_key": ground_station_receive_gain_key,
        "satellite_tx_antenna": satellite_tx_antenna,
        "satellite_rx_antenna": satellite_rx_antenna,
    }

# Optimization function
def optimize_link_budget(orbital_params, ground_segment, satellite_params, modems_sdrs, modulation_bits_per_symbol, antennas, required_Eb_N0_dB_dict):
    """
    Iterate over all combinations of ground stations, modems, and satellite antennas to find the configuration
    that maximizes the data rate for both uplink and downlink.
    """
    best_configuration = None
    max_total_data_rate = 0

    for ground_station in ground_segment:
        # Check if the ground station supports the uplink frequency
        if not (ground_station["sup_freq"][0] <= orbital_params["uplink_frequency"] <= ground_station["sup_freq"][1]):
            continue

        # Check if the ground station supports the downlink frequency
        downlink_supported = False
        gs_downlink_bands = ["sdown_bw", "xdown_bw", "kadown_bw"]
        gs_downlink_gains = ["sdown_gt", "xdown_gt", "kadown_gt"]
        for band, gain_key in zip(gs_downlink_bands, gs_downlink_gains):
            if ground_station[band] is not None:
                if ground_station[band][0] <= orbital_params["downlink_frequency"] <= ground_station[band][1]:
                    downlink_supported = True
                    ground_station_receive_gain_key = gain_key  # This key will be used to get the gain
                    break
        if not downlink_supported:
            continue

        for modem in modems_sdrs:
            # Check if modem supports any modulation
            if not modem["modulations"]:
                continue

            # Iterate over satellite transmit antennas
            for sat_tx_ant in antennas:
                # Check if satellite transmit antenna supports downlink frequency
                if not (sat_tx_ant["frequency_range"][0] <= orbital_params["downlink_frequency"] <= sat_tx_ant["frequency_range"][1]):
                    continue

                for sat_rx_ant in antennas:
                    # Check if satellite receive antenna supports uplink frequency
                    if not (sat_rx_ant["frequency_range"][0] <= orbital_params["uplink_frequency"] <= sat_rx_ant["frequency_range"][1]):
                        continue

                    # Calculate link budget
                    results = calculate_link_budget(
                        orbital_params,
                        ground_station,
                        ground_station_receive_gain_key,
                        satellite_params,
                        modem,
                        modulation_bits_per_symbol,
                        sat_tx_ant,
                        sat_rx_ant,
                        required_Eb_N0_dB_dict
                    )

                    if results is None:
                        continue

                    # Sum of uplink and downlink data rates
                    total_data_rate = results["data_rate_uplink_Mbps"] + results["data_rate_downlink_Mbps"]

                    if total_data_rate > max_total_data_rate:
                        max_total_data_rate = total_data_rate
                        best_configuration = results

    if best_configuration is not None:
        return best_configuration
    else:
        return None

# Get user input for desired BER
user_BER = float(input("Enter the desired BER (e.g., 1e-6): "))

# Compute required Eb/N0 for each modulation scheme
required_Eb_N0_dB_dict = {}
for modulation in modulation_bits_per_symbol:
    required_Eb_N0_dB = get_required_Eb_N0_dB(modulation, user_BER)
    required_Eb_N0_dB_dict[modulation] = required_Eb_N0_dB

# Run optimization
best_config = optimize_link_budget(orbital_parameters, ground_segment, satellite_parameters, modems_sdrs, modulation_bits_per_symbol, antennas, required_Eb_N0_dB_dict)

if best_config is not None:
    modem = best_config["modem"]
    ground_station = best_config["ground_station"]
    ground_station_receive_gain_key = best_config["ground_station_receive_gain_key"]
    ground_station_receive_gain = ground_station[ground_station_receive_gain_key]
    print(f"Optimal Ground Station: {ground_station['name']}")
    print(f"Optimal Modem/SDR: {modem['name']}")
    print(f"Satellite Transmit Antenna: {best_config['satellite_tx_antenna']['name']} with Gain {best_config['satellite_tx_antenna']['gain']} dBi")
    print(f"Satellite Receive Antenna: {best_config['satellite_rx_antenna']['name']} with Gain {best_config['satellite_rx_antenna']['gain']} dBi")
    print(f"Ground Station Receive Antenna Gain ({ground_station_receive_gain_key}): {ground_station_receive_gain} dBi")
    print("\n=== Uplink (Ground to Satellite) ===")
    print(f"Uplink Received Power at Satellite: {best_config['uplink_power_dBm']:.2f} dBm")
    print(f"Uplink SNR: {best_config['SNR_uplink_dB']:.2f} dB")
    print(f"Best Uplink Modulation: {best_config['best_uplink_modulation']}")
    print(f"Achievable Uplink Data Rate: {best_config['data_rate_uplink_Mbps']:.6f} Mbps")
    print("\n=== Downlink (Satellite to Ground) ===")
    print(f"Downlink Received Power at Ground Station: {best_config['downlink_power_dBm']:.2f} dBm")
    print(f"Downlink SNR: {best_config['SNR_downlink_dB']:.2f} dB")
    print(f"Best Downlink Modulation: {best_config['best_downlink_modulation']}")
    print(f"Achievable Downlink Data Rate: {best_config['data_rate_downlink_Mbps']:.6f} Mbps")
else:
    print("No suitable configuration found.")
