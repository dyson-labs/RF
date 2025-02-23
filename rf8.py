import numpy as np
from scipy.special import erfc, erfcinv
from scipy.optimize import bisect

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
    "GMSK": 1,
    "MSK": 1,
    "2FSK": 1,
    "4FSK": 2,
    "2GFSK": 1,
    "4GFSK": 2,
    "8PSK": 3,
    "16QAM": 4,
    "16APSK": 4,
    "32APSK": 5,
    "256APSK": 8,
    "GFSK": 1,
    "FSK": 1,
}

modems_sdrs = [
    {
        "type": "transceiver",
        "name": "UHF Transceiver II - Enduro",
        "data_rate": (0.1, 19.2),  # kbps
        "modulations": ["OOK", "GMSK", "2FSK", "4FSK", "4GFSK"],
        "transmit_power_W": [1,2],  # Watts
        "receive_sensitivity_dBm": -121,
        "rx_frequency_range": (400, 403),  # MHz
        "tx_frequency_range": (430, 440),  # MHz
    },
    {
        "type": "transceiver",
        "name": "S Band Transceiver - Enduro",
        "data_rate": (0.1, 125),  # kbps
        "modulations": ["FSK", "MSK", "GFSK","GMSK"],
        "transmit_power_W": (0.4,2),  # Watts
        "receive_sensitivity_dBm": -121,
        "rx_frequency_range": (2025, 2110),  # MHz
        "tx_frequency_range": (2200, 2290),  # MHz
    },
    {
        "type": "transceiver",
        "name": "Flexible and Miniaturised Transceiver - GOMSpace",
        "data_rate": (0.1, 38.4),  # kbps
        "modulations": ["GFSK","GMSK"],
        "transmit_power_W": 1,  # Watts
        "receive_sensitivity_dBm": -137,
        "rx_frequency_range": (430, 440),  # MHz
        "tx_frequency_range": (430, 440),  # MHz
    },
    {
        "type": "transceiver",
        "name": "NanoCom AX2150 - GOMSpace",
        "data_rate": (9.6, 96),  # kbps
        "modulations": ["GFSK","GMSK"],
        "transmit_power_W": (0.008,0.5),  # Watts
        "receive_sensitivity_dBm": -113,
        "rx_frequency_range": (2025, 2110),  # MHz
        "tx_frequency_range": (2200,2290),  # MHz
    },
    {
        "type": "transmitter",
        "name": "S Band Transmitter - Enduro",
        "data_rate": 20000,
        "modulations": ["QPSK", "8PSK", "16APSK"],
        "transmit_power_W": (0,2),
        "frequency_range": [(2200, 2290),(2400,2450)],
    },
    {
        "type": "transmitter",
        "name": "X Band Transmitter - Enduro",
        "data_rate": 150000,
        "modulations": ["QPSK", "8PSK", "16APSK", "32APSK"],
        "transmit_power_W": (0,2),
        "frequency_range": (7900, 8400),
    },
    {
        "type": "transmitter",
        "name": "K Band Transmitter - Enduro",
        "data_rate": 1000000,
        "modulations": ["QPSK", "8PSK", "16APSK", "32APSK", "256APSK"],
        "transmit_power_W": (0,2),
        "frequency_range": (25500, 27000),
    },
]

# Ground stations
ground_segment = [
    {"name": "VIASAT PENDER", "sup_freq": (2025, 2110), "uEIRP": 53.2, "sdown_bw": (2200, 2290), "sdown_gt": 17, "xdown_bw": (8025, 8400), "xdown_gt": 30, "kadown_bw": None, "kadown_gt": None},
    {"name": "VIASAT GUILDFORD", "sup_freq": (2025, 2110), "uEIRP": 53.2, "sdown_bw": (2200, 2290), "sdown_gt": 17, "xdown_bw": (8025, 8400), "xdown_gt": 30, "kadown_bw": None, "kadown_gt": None},
    {"name": "VIASAT ALICE", "sup_freq": (2025, 2110), "uEIRP": 65.0, "sdown_bw": (2200, 2290), "sdown_gt": 18, "xdown_bw": (8025, 8400), "xdown_gt": 32, "kadown_bw": (25500, 27000), "kadown_gt": 34.5},
    {"name": "VIASAT GHANA", "sup_freq": (2025, 2110), "uEIRP": 65.0, "sdown_bw": (2200, 2290), "sdown_gt": 18, "xdown_bw": (8025, 8400), "xdown_gt": 32, "kadown_bw": (25500, 27000), "kadown_gt": 34.5},
    {"name": "ATLAS PAUMALU", "sup_freq": (2025, 2120), "uEIRP": 50.0, "sdown_bw": (2200, 2300), "sdown_gt": 21, "xdown_bw": (7900, 8500), "xdown_gt": 31, "kadown_bw": None, "kadown_gt": None},
]

# Function to calculate Free Space Path Loss (FSPL)
def calculate_fspl(altitude, frequency):
    """Calculate Free Space Path Loss (FSPL) in dB."""
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
    if modulation_scheme in ["BPSK", "QPSK", "GMSK", "GFSK", "FSK", "OOK", "MSK"]:
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
    elif modulation_scheme in ["16QAM", "16APSK"]:
        M = 16
        k = np.log2(M)
        def func(Eb_N0_dB):
            Eb_N0_linear = 10 ** (Eb_N0_dB / 10)
            Q_arg = np.sqrt(3 * Eb_N0_linear * k / (M - 1))
            BER = (4 / k) * (1 - 1 / np.sqrt(M)) * Q(Q_arg)
            return BER - BER_target
        Eb_N0_dB = bisect(func, 0, 30)
        return Eb_N0_dB
    elif modulation_scheme == "32APSK":
        # Approximate Eb/N0 value for 32APSK
        return 18  # dB (approximate)
    elif modulation_scheme == "256APSK":
        # Approximate Eb/N0 value for 256APSK
        return 25  # dB (approximate)
    else:
        return 30  # Default high value

# Function to calculate link budget
def calculate_link_budget(orbital_params, ground_station, ground_station_receive_gain_key, satellite_params, device, modulation_bits_per_symbol, satellite_tx_antenna, satellite_rx_antenna, required_Eb_N0_dB_uplink_dict, required_Eb_N0_dB_downlink_dict):
    """
    Calculate link budget for both uplink and downlink.

    Returns a dictionary with all calculated parameters.
    """
    # Extract maximum data rate in bps
    if isinstance(device["data_rate"], (tuple, list)):
        max_data_rate_kbps = device["data_rate"][1]
    else:
        max_data_rate_kbps = device["data_rate"]
    max_data_rate_bps = max_data_rate_kbps * 1e3  # Convert kbps to bps

    # Boltzmann's constant
    k = 1.38e-23  # J/K
    T_system = orbital_params["T_system"]  # K

    # Initialize best modulation schemes
    best_uplink_modulation = None
    best_downlink_modulation = None
    max_bits_per_symbol_uplink = 0
    max_bits_per_symbol_downlink = 0
    best_uplink_data_rate = 0
    best_downlink_data_rate = 0

    # Uplink calculations (Ground to Satellite)
    EIRP_ground_dBm = ground_station["uEIRP"]  # dBm
    G_satellite_rx_dBi = satellite_rx_antenna["gain"]  # dBi
    uplink_fspl = calculate_fspl(orbital_params["altitude"], orbital_params["uplink_frequency"])
    L_uplink_dB = 0  # Losses
    P_r_uplink_dBm = EIRP_ground_dBm - uplink_fspl + G_satellite_rx_dBi - L_uplink_dB

    # Downlink calculations (Satellite to Ground)
    P_satellite_tx_W = satellite_params["P_satellite_transmit_W"]  # W
    P_satellite_tx_dBm = 10 * np.log10(P_satellite_tx_W * 1000)  # dBm
    G_satellite_tx_dBi = satellite_tx_antenna["gain"]  # dBi
    EIRP_satellite_dBm = P_satellite_tx_dBm + G_satellite_tx_dBi  # dBm
    G_ground_rx_dBi = ground_station[ground_station_receive_gain_key]  # dBi
    downlink_fspl = calculate_fspl(orbital_params["altitude"], orbital_params["downlink_frequency"])
    L_downlink_dB = 0  # Losses
    P_r_downlink_dBm = EIRP_satellite_dBm - downlink_fspl + G_ground_rx_dBi - L_downlink_dB

    # Iterate over modulations to find the best one for uplink and downlink
    for modulation in device["modulations"]:
        bits_per_symbol = modulation_bits_per_symbol.get(modulation, 1)

        # Uplink
        symbol_rate_uplink = max_data_rate_bps / bits_per_symbol
        B_uplink = symbol_rate_uplink  # Hz (assuming Nyquist bandwidth)
        P_noise_uplink_W = k * T_system * B_uplink
        P_noise_uplink_dBm = 10 * np.log10(P_noise_uplink_W) + 30  # dBm
        SNR_uplink_dB = P_r_uplink_dBm - P_noise_uplink_dBm
        Eb_N0_uplink_dB = SNR_uplink_dB - 10 * np.log10(bits_per_symbol)
        required_Eb_N0_dB_uplink = required_Eb_N0_dB_uplink_dict.get(modulation, 30)

        if Eb_N0_uplink_dB >= required_Eb_N0_dB_uplink and bits_per_symbol >= max_bits_per_symbol_uplink:
            max_bits_per_symbol_uplink = bits_per_symbol
            best_uplink_modulation = modulation
            best_uplink_data_rate = symbol_rate_uplink * bits_per_symbol  # bps

        # Downlink
        symbol_rate_downlink = max_data_rate_bps / bits_per_symbol
        B_downlink = symbol_rate_downlink  # Hz
        P_noise_downlink_W = k * T_system * B_downlink
        P_noise_downlink_dBm = 10 * np.log10(P_noise_downlink_W) + 30  # dBm
        SNR_downlink_dB = P_r_downlink_dBm - P_noise_downlink_dBm
        Eb_N0_downlink_dB = SNR_downlink_dB - 10 * np.log10(bits_per_symbol)
        required_Eb_N0_dB_downlink = required_Eb_N0_dB_downlink_dict.get(modulation, 30)

        if Eb_N0_downlink_dB >= required_Eb_N0_dB_downlink and bits_per_symbol >= max_bits_per_symbol_downlink:
            max_bits_per_symbol_downlink = bits_per_symbol
            best_downlink_modulation = modulation
            best_downlink_data_rate = symbol_rate_downlink * bits_per_symbol  # bps

    if best_uplink_modulation is None or best_downlink_modulation is None:
        return None  # No suitable modulation found

    # Convert data rates to Mbps
    data_rate_uplink_Mbps = best_uplink_data_rate / 1e6
    data_rate_downlink_Mbps = best_downlink_data_rate / 1e6

    return {
        "uplink_power_dBm": P_r_uplink_dBm,
        "downlink_power_dBm": P_r_downlink_dBm,
        "SNR_uplink_dB": SNR_uplink_dB,
        "SNR_downlink_dB": SNR_downlink_dB,
        "data_rate_uplink_Mbps": data_rate_uplink_Mbps,
        "data_rate_downlink_Mbps": data_rate_downlink_Mbps,
        "best_uplink_modulation": best_uplink_modulation,
        "best_downlink_modulation": best_downlink_modulation,
        "device": device,
        "ground_station": ground_station,
        "ground_station_receive_gain_key": ground_station_receive_gain_key,
        "satellite_tx_antenna": satellite_tx_antenna,
        "satellite_rx_antenna": satellite_rx_antenna,
    }

# Optimization function
def optimize_link_budget(orbital_params, ground_segment, satellite_params, devices, modulation_bits_per_symbol, antennas, required_Eb_N0_dB_uplink_dict, required_Eb_N0_dB_downlink_dict):
    """
    Iterate over all combinations to find the configuration that maximizes the data rate.
    """
    best_configuration = None
    max_total_data_rate = 0

    for ground_station in ground_segment:
        # Ground station frequency checks
        if not (ground_station["sup_freq"][0] <= orbital_params["uplink_frequency"] <= ground_station["sup_freq"][1]):
            continue

        downlink_supported = False
        gs_downlink_bands = ["sdown_bw", "xdown_bw", "kadown_bw"]
        gs_downlink_gains = ["sdown_gt", "xdown_gt", "kadown_gt"]
        for band, gain_key in zip(gs_downlink_bands, gs_downlink_gains):
            if ground_station[band]:
                if ground_station[band][0] <= orbital_params["downlink_frequency"] <= ground_station[band][1]:
                    downlink_supported = True
                    ground_station_receive_gain_key = gain_key
                    break
        if not downlink_supported:
            continue

        for device in devices:
            # Device frequency checks
            tx_freq_range = device.get("tx_frequency_range")
            rx_freq_range = device.get("rx_frequency_range")
            if tx_freq_range and not (tx_freq_range[0] <= orbital_params["uplink_frequency"] <= tx_freq_range[1]):
                continue
            if rx_freq_range and not (rx_freq_range[0] <= orbital_params["downlink_frequency"] <= rx_freq_range[1]):
                continue

            # Modulations check
            if not device["modulations"]:
                continue

            # Iterate over satellite antennas
            for sat_tx_ant in antennas:
                if not (sat_tx_ant["frequency_range"][0] <= orbital_params["downlink_frequency"] <= sat_tx_ant["frequency_range"][1]):
                    continue
                for sat_rx_ant in antennas:
                    if not (sat_rx_ant["frequency_range"][0] <= orbital_params["uplink_frequency"] <= sat_rx_ant["frequency_range"][1]):
                        continue

                    # Calculate link budget
                    results = calculate_link_budget(
                        orbital_params,
                        ground_station,
                        ground_station_receive_gain_key,
                        satellite_params,
                        device,
                        modulation_bits_per_symbol,
                        sat_tx_ant,
                        sat_rx_ant,
                        required_Eb_N0_dB_uplink_dict,
                        required_Eb_N0_dB_downlink_dict
                    )

                    if results is None:
                        continue

                    total_data_rate = results["data_rate_uplink_Mbps"] + results["data_rate_downlink_Mbps"]
                    if total_data_rate > max_total_data_rate:
                        max_total_data_rate = total_data_rate
                        best_configuration = results

    return best_configuration

# Get user input for orbital parameters
print("Enter Orbital Parameters:")
try:
    altitude = float(input("Satellite Altitude above Earth's surface (in km, e.g., 2000): "))
    uplink_frequency = float(input("Uplink Frequency (in MHz, e.g., 2025): "))
    downlink_frequency = float(input("Downlink Frequency (in MHz, e.g., 25600): "))
    T_system = float(input("System Noise Temperature (in Kelvin, e.g., 290): "))
except ValueError:
    print("Invalid input. Please enter numerical values.")
    exit()

orbital_parameters = {
    "altitude": altitude,  # km
    "uplink_frequency": uplink_frequency,  # MHz
    "downlink_frequency": downlink_frequency,  # MHz
    "T_system": T_system,  # K
}

# Satellite parameters
satellite_parameters = {
    "G_satellite_receiver_dBi": None,  # To be selected from antennas
    "G_satellite_transmit_dBi": None,  # To be selected from antennas
    "P_satellite_transmit_W": 5,  # W
}

# Get user input for desired BER
try:
    user_BER_uplink = float(input("\nEnter the desired BER for UPLINK (e.g., 1e-6): "))
    user_BER_downlink = float(input("Enter the desired BER for DOWNLINK (e.g., 1e-6): "))
except ValueError:
    print("Invalid input. Please enter numerical values.")
    exit()

# Compute required Eb/N0
required_Eb_N0_dB_uplink_dict = {}
required_Eb_N0_dB_downlink_dict = {}
for modulation in modulation_bits_per_symbol:
    required_Eb_N0_dB_uplink = get_required_Eb_N0_dB(modulation, user_BER_uplink)
    required_Eb_N0_dB_downlink = get_required_Eb_N0_dB(modulation, user_BER_downlink)
    required_Eb_N0_dB_uplink_dict[modulation] = required_Eb_N0_dB_uplink
    required_Eb_N0_dB_downlink_dict[modulation] = required_Eb_N0_dB_downlink

# Run optimization
best_config = optimize_link_budget(
    orbital_parameters,
    ground_segment,
    satellite_parameters,
    modems_sdrs,
    modulation_bits_per_symbol,
    antennas,
    required_Eb_N0_dB_uplink_dict,
    required_Eb_N0_dB_downlink_dict
)

if best_config:
    device = best_config["device"]
    ground_station = best_config["ground_station"]
    ground_station_receive_gain_key = best_config["ground_station_receive_gain_key"]
    ground_station_receive_gain = ground_station[ground_station_receive_gain_key]
    print(f"\nOptimal Ground Station: {ground_station['name']}")
    print(f"Optimal Device: {device['name']}")
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
