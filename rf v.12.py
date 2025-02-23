import numpy as np
from scipy.special import erfc
from skyfield.api import Topos, EarthSatellite, load, wgs84

antennas = [
    # VHF (30 MHz - 300 MHz)
    {"name": "FAKE VHF Yagi Antenna", "frequency_range": (30, 300), "gain": 7, "type": "Yagi-Uda"},

    # UHF (300 MHz - 3 GHz)
    {"name": "FAKE UHF Log-Periodic Dipole Array", "frequency_range": (300, 1000), "gain": 10, "type": "Log-Periodic"},
    {"name": "FAKE UHF Patch Antenna", "frequency_range": (300, 1000), "gain": 8, "type": "Patch"},

    # L-Band (1 GHz - 2 GHz)
    {"name": "FAKE L-Band Helical Antenna", "frequency_range": (1000, 2000), "gain": 12, "type": "Helical"},
    {"name": "FALCCON - RW", "frequency_range": (1000, 2000), "gain": 21, "type": "Dipole Array"},
    {"name": "L Band Patch Antenna - Printech", "frequency_range": (1563,1587), "gain": 5, "type": "Patch"},

    # S-Band (2 GHz - 4 GHz)
    {"name": "FAKE S-Band Parabolic Reflector", "frequency_range": (2000, 4000), "gain": 15, "type": "Parabolic"},
    {"name": "AC-2000 - AAC", "frequency_range": (2000, 2300), "gain": 5.2, "type": "Parabolic Dish"},
    {"name": "SANT S-band Patch Antenna", "frequency_range": (2200, 2290), "gain": 6, "type": "Parabolic"},

    # C-Band (4 GHz - 8 GHz)
    {"name": "FAKE C-Band Horn Antenna", "frequency_range": (4000, 8000), "gain": 18, "type": "Horn"},

    # X-Band (8 GHz - 12 GHz)
    {"name": "FAKE X-Band Parabolic Dish", "frequency_range": (8000, 12000), "gain": 20, "type": "Parabolic Dish"},
    {"name": "X-Band Patch Antenna - Enduro", "frequency_range": (8025, 8400), "gain": 6, "type": "Patch"},
    {"name": "4x4 X-Band Patch Array - Enduro", "frequency_range": (8025, 8400), "gain": 16, "type": "Patch Array"},
    {"name": "XANT X-Band Patch Antenna - Cubecom", "frequency_range": (8000, 8400), "gain": 8, "type": "Patch"},
    {"name": "XPLANT X-band Payload Antenna", "frequency_range": (8500, 9600), "gain": 8, "type": "Patch Array"},

    # Ku-Band (12 GHz - 18 GHz)
    {"name": "FAKE Ku-Band Microstrip Antenna", "frequency_range": (12000, 18000), "gain": 23, "type": "Microstrip"},

    # K-Band (18 GHz - 27 GHz)
    {"name": "FAKE K-Band Waveguide Antenna", "frequency_range": (18000, 27000), "gain": 25, "type": "Waveguide"},
    {"name": "K-Band 4x4 Patch Array - Enduro", "frequency_range": (17700, 20200), "gain": 16, "type": "Waveguide"},

    # Ka-Band (27 GHz - 40 GHz)
    {"name": "FAKE Ka-Band Horn Antenna", "frequency_range": (27000, 40000), "gain": 30, "type": "Horn"},
    {"name": "PAN-5151-64-KA - ReliaSat", "frequency_range": (27000, 31000), "gain": 20, "type": "Panel Array"},
    {"name": "4x4 X-Band Patch Array - Enduro", "frequency_range": (8025, 8400), "gain": 16, "type": "Patch Array"},
    
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
        "type": "multi-band transceiver",
        "name": "NanoCom Link SX",
        "x_band_tx_data_rate": (1,225000),  # kbps
        "s_band_tx_data_rate": (500, 7500),
        "x_band_tx_modulations": ["QPSK", "8PSK", "16APSK","32APSK"],
        "s_band_tx_modulations": ["BPSK, QPSK"],
        "x_band_transmit_power_dBm": [0.1,33],  # dBm
        "s_band_transmit_power_dBm": [0.1,32], # dBm
        "receive_sensitivity_dBm": -100, # -111 dBm BSPK 1 MBd, -100 dBm QPSK 7 MBd
        "rx_frequency_range": (2025, 2110),  # MHz
        "tx_frequency_range": [(2200,2290),(8000,8400)],  # MHz
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

# Dictionary of TLEs for different orbital lanes
satellite_tles = {
    "LEO": {
        "name": "ISS (ZARYA)",
        "description": "Low Earth Orbit (LEO) satellite - International Space Station.",
        "tle": [
            "1 25544U 98067A   23314.54692130  .00007237  00000-0  13252-3 0  9992",
            "2 25544  51.6425 282.3050 0002927 134.1747  13.9034 15.49925521424794"
        ],
    },
    "MEO": {
        "name": "GPS BIIR-2  (PRN 18)",
        "description": "Medium Earth Orbit (MEO) satellite - Part of the GPS constellation.",
        "tle": [
            "1 24876U 97033A   23314.47420425  .00000025  00000-0  00000-0 0  9997",
            "2 24876  54.8326 305.6921 0152963  58.7624 304.8789  2.00569910172652"
        ],
    },
    "GEO": {
        "name": "GOES-16",
        "description": "Geostationary Orbit (GEO) satellite - Weather monitoring.",
        "tle": [
            "1 41866U 16071A   23314.57030787 -.00000267  00000-0  00000+0 0  9998",
            "2 41866   0.0171 121.1528 0000291 312.5125  47.5398  1.00272067 25134"
        ],
    },
    "Dawn-Dusk Orbit": {
        "name": "Sentinel-2A",
        "description": "Sun-synchronous dawn-dusk orbit satellite for Earth observation.",
        "tle": [
            "1 40697U 15028A   23314.46294037  .00000027  00000-0  23210-4 0  9995",
            "2 40697  98.5672  44.5289 0001275  90.3575 269.7627 14.30883213437250"
        ],
    },
}


# Ground stations with Skyfield Topos and additional properties
ground_segment = [
    {
        "name": "VIASAT PENDER",
        "location": Topos(latitude_degrees=49.1, longitude_degrees=-123.9, elevation_m=30),  # Placeholder coordinates
        "sup_freq": (2025, 2110),
        "uEIRP": 53.2,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 17,
        "xdown_bw": (8025, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "VIASAT GUILDFORD",
        "location": Topos(latitude_degrees=51.2, longitude_degrees=-0.6, elevation_m=70),  # Placeholder coordinates
        "sup_freq": (2025, 2110),
        "uEIRP": 53.2,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 17,
        "xdown_bw": (8025, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "VIASAT ALICE",
        "location": Topos(latitude_degrees=-23.7, longitude_degrees=133.9, elevation_m=600),  # Placeholder coordinates
        "sup_freq": (2025, 2110),
        "uEIRP": 65.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (8025, 8400),
        "xdown_gt": 32,
        "kadown_bw": (25500, 27000),
        "kadown_gt": 34.5,
    },
    {
        "name": "VIASAT GHANA",
        "location": Topos(latitude_degrees=5.6, longitude_degrees=-0.2, elevation_m=50),  # Placeholder coordinates
        "sup_freq": (2025, 2110),
        "uEIRP": 65.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (8025, 8400),
        "xdown_gt": 32,
        "kadown_bw": (25500, 27000),
        "kadown_gt": 34.5,
    },
    {
        "name": "ATLAS PAUMALU",
        "location": Topos(latitude_degrees=21.6, longitude_degrees=-158.0, elevation_m=100),  # Placeholder coordinates
        "sup_freq": (2025, 2120),
        "uEIRP": 50.0,
        "sdown_bw": (2200, 2300),
        "sdown_gt": 21,
        "xdown_bw": (7900, 8500),
        "xdown_gt": 31,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Alaska 1",
        "location": Topos(latitude_degrees=64.2008, longitude_degrees=-149.4937, elevation_m=100),  # Placeholder coordinates for Alaska
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_bw": (2200, 2290),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_bw": (7750, 8400),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Bahrain 1",
        "location": Topos(latitude_degrees=26.0667, longitude_degrees=50.5577, elevation_m=50),  # Placeholder coordinates for Bahrain
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_bw": (2200, 2290),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_bw": (7750, 8400),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Cape Town 1",
        "location": Topos(latitude_degrees=-33.9249, longitude_degrees=18.4241, elevation_m=50),  # Placeholder coordinates for Cape Town
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_bw": (2200, 2290),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_bw": (7750, 8400),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Dubbo 1",
        "location": Topos(latitude_degrees=-32.2569, longitude_degrees=148.6011, elevation_m=50),  # Placeholder coordinates for Dubbo
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_bw": (2200, 2290),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_bw": (7750, 8400),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Hawaii 1",
        "location": Topos(latitude_degrees=19.8968, longitude_degrees=-155.5828, elevation_m=100),  # Placeholder coordinates for Hawaii
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_bw": (2200, 2290),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_bw": (7750, 8400),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Ireland 1",
        "location": Topos(latitude_degrees=53.1424, longitude_degrees=-7.6921, elevation_m=50),  # Placeholder coordinates for Ireland
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Ohio 1",
        "location": Topos(latitude_degrees=40.4173, longitude_degrees=-82.9071, elevation_m=50),  # Placeholder coordinates for Ohio
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Oregon 1",
        "location": Topos(latitude_degrees=43.8041, longitude_degrees=-120.5542, elevation_m=100),  # Placeholder coordinates for Oregon
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Punta Arenas 1",
        "location": Topos(latitude_degrees=-53.1638, longitude_degrees=-70.9171, elevation_m=50),  # Placeholder coordinates for Punta Arenas
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Seoul 1",
        "location": Topos(latitude_degrees=37.5665, longitude_degrees=126.9780, elevation_m=50),  # Placeholder coordinates for Seoul
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Singapore 1",
        "location": Topos(latitude_degrees=1.3521, longitude_degrees=103.8198, elevation_m=50),  # Placeholder coordinates for Singapore
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Stockholm 1",
        "location": Topos(latitude_degrees=59.3293, longitude_degrees=18.0686, elevation_m=50),  # Placeholder coordinates for Stockholm
        "sup_freq": (2025, 2110),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_bw": (2200, 2290),
        "sdown_gt": 18,
        "xdown_bw": (7750, 8400),
        "xdown_gt": 30,
        "kadown_bw": None,
        "kadown_gt": None,
    },
]

# Function to calculate Free Space Path Loss (FSPL)
def calculate_fspl(distance_m, frequency):
    """Calculate Free Space Path Loss (FSPL) in dB."""
    d = distance_m  # Distance in meters
    f = frequency * 1e6  # Frequency in Hz
    c = 3e8  # Speed of light in m/s
    fspl = 20 * np.log10(4 * np.pi * d * f / c)
    return fspl

# Function to compute BER given Eb/N0 and modulation
def compute_BER(Eb_N0_dB, modulation_scheme):
    Eb_N0_linear = 10 ** (Eb_N0_dB / 10)
    if modulation_scheme in ["BPSK", "QPSK", "GMSK", "GFSK", "FSK", "OOK", "MSK"]:
        # BER for BPSK/QPSK
        BER = 0.5 * erfc(np.sqrt(Eb_N0_linear))
    elif modulation_scheme == "8PSK":
        BER = erfc(np.sqrt(1.5 * Eb_N0_linear * np.log2(8) / (8 - 1)))
    elif modulation_scheme in ["16QAM", "16APSK"]:
        BER = 0.75 * erfc(np.sqrt(0.1 * Eb_N0_linear))
    elif modulation_scheme == "32APSK":
        # Approximate BER for 32APSK
        BER = erfc(np.sqrt(0.068 * Eb_N0_linear))
    elif modulation_scheme == "256APSK":
        # Approximate BER for 256APSK
        BER = erfc(np.sqrt(0.0156 * Eb_N0_linear))
    else:
        # Default high BER for unsupported modulations
        BER = 1.0
    return BER

# Function to get possible transmit powers for the satellite transmitter
def get_possible_satellite_transmit_powers(transmit_power_W):
    """
    Given a transmit_power_W field, return a list of possible transmit powers to try for the satellite transmitter.
    """
    if isinstance(transmit_power_W, list):
        # It's a list of values
        possible_powers = transmit_power_W
    elif isinstance(transmit_power_W, tuple):
        # It's a range
        min_power = transmit_power_W[0]
        max_power = transmit_power_W[1]
        num_steps = 10
        if max_power == min_power:
            possible_powers = [min_power]
        else:
            step_size = (max_power - min_power) / (num_steps - 1)
            possible_powers = [min_power + i * step_size for i in range(num_steps)]
    else:
        # It's a single value
        possible_powers = [transmit_power_W]
    return possible_powers

# Function to calculate link budget
def calculate_link_budget(orbital_params, ground_station, satellite_params, uplink_device, satellite_transmitter, modulation_bits_per_symbol, satellite_tx_antenna, satellite_rx_antenna, desired_data_rate_uplink_bps, desired_data_rate_downlink_bps, satellite_transmit_power_W, distance_m):
    """
    Calculate link budget for both uplink and downlink.
    Returns a dictionary with all calculated parameters.
    """
    # Boltzmann's constant
    k = 1.38e-23  # J/K

    # Uplink calculations (Ground to Satellite)
    # Fixed transmit power for ground station
    uplink_transmit_power_W = uplink_device["transmit_power_W"]
    if isinstance(uplink_transmit_power_W, (list, tuple)):
        # Select the maximum power if it's a list or tuple
        uplink_transmit_power_W = max(uplink_transmit_power_W)
    P_ground_tx_dBm = 10 * np.log10(uplink_transmit_power_W * 1000)  # dBm
    G_ground_tx_dBi = ground_station["G_ground_tx_dBi"]  # dBi
    EIRP_ground_dBm = P_ground_tx_dBm + G_ground_tx_dBi  # dBm
    G_satellite_rx_dBi = satellite_rx_antenna["gain"]  # dBi
    uplink_fspl = calculate_fspl(distance_m, orbital_params["uplink_frequency"])
    L_uplink_dB = 0  # Losses
    P_r_uplink_dBm = EIRP_ground_dBm - uplink_fspl + G_satellite_rx_dBi - L_uplink_dB

    # Downlink calculations (Satellite to Ground)
    P_satellite_tx_W = satellite_transmit_power_W  # W
    P_satellite_tx_dBm = 10 * np.log10(P_satellite_tx_W * 1000)  # dBm
    G_satellite_tx_dBi = satellite_tx_antenna["gain"]  # dBi
    EIRP_satellite_dBm = P_satellite_tx_dBm + G_satellite_tx_dBi  # dBm
    G_ground_rx_dBi = ground_station[satellite_params["ground_station_receive_gain_key"]]  # dBi
    downlink_fspl = calculate_fspl(distance_m, orbital_params["downlink_frequency"])
    L_downlink_dB = 0  # Losses
    P_r_downlink_dBm = EIRP_satellite_dBm - downlink_fspl + G_ground_rx_dBi - L_downlink_dB

    # Initialize variables
    selected_uplink_modulation = None
    selected_downlink_modulation = None
    achieved_uplink_BER = None
    achieved_downlink_BER = None

    # Uplink Modulation Selection
    for modulation in uplink_device["modulations"]:
        bits_per_symbol = modulation_bits_per_symbol.get(modulation, 1)
        # Calculate required symbol rate
        required_symbol_rate_uplink = desired_data_rate_uplink_bps / bits_per_symbol
        # Check if the device supports this symbol rate
        max_device_data_rate_bps = uplink_device["data_rate"][1] * 1e3 if isinstance(uplink_device["data_rate"], (tuple, list)) else uplink_device["data_rate"] * 1e3
        max_symbol_rate_device = max_device_data_rate_bps / bits_per_symbol
        if required_symbol_rate_uplink > max_symbol_rate_device:
            continue  # Cannot support desired data rate with this modulation

        # Calculate noise power
        B_uplink = required_symbol_rate_uplink  # Assuming Nyquist bandwidth
        T_system_uplink = satellite_params["T_system_satellite"]  # K
        P_noise_uplink_W = k * T_system_uplink * B_uplink
        P_noise_uplink_dBm = 10 * np.log10(P_noise_uplink_W) + 30  # dBm

        # Calculate SNR
        SNR_uplink_dB = P_r_uplink_dBm - P_noise_uplink_dBm
        Eb_N0_uplink_dB = SNR_uplink_dB - 10 * np.log10(bits_per_symbol)

        # Compute BER
        BER_uplink = compute_BER(Eb_N0_uplink_dB, modulation)

        # Check if BER meets requirement (assuming a threshold, e.g., BER < 1e-5)
        if BER_uplink < 1e-5:
            # Select the modulation
            selected_uplink_modulation = modulation
            achieved_uplink_BER = BER_uplink
            break  # Select the first modulation that meets the data rate and BER requirement

    # Downlink Modulation Selection
    for modulation in satellite_transmitter["modulations"]:
        bits_per_symbol = modulation_bits_per_symbol.get(modulation, 1)
        # Calculate required symbol rate
        required_symbol_rate_downlink = desired_data_rate_downlink_bps / bits_per_symbol
        # Check if the satellite transmitter supports this symbol rate
        max_device_data_rate_bps = satellite_transmitter["data_rate"] * 1e3 if not isinstance(satellite_transmitter["data_rate"], (tuple, list)) else satellite_transmitter["data_rate"][1] * 1e3
        max_symbol_rate_device = max_device_data_rate_bps / bits_per_symbol
        if required_symbol_rate_downlink > max_symbol_rate_device:
            continue  # Cannot support desired data rate with this modulation

        # Calculate noise power
        B_downlink = required_symbol_rate_downlink  # Assuming Nyquist bandwidth
        T_system_downlink = ground_station["T_system_ground"]  # K
        P_noise_downlink_W = k * T_system_downlink * B_downlink
        P_noise_downlink_dBm = 10 * np.log10(P_noise_downlink_W) + 30  # dBm

        # Calculate SNR
        SNR_downlink_dB = P_r_downlink_dBm - P_noise_downlink_dBm
        Eb_N0_downlink_dB = SNR_downlink_dB - 10 * np.log10(bits_per_symbol)

        # Compute BER
        BER_downlink = compute_BER(Eb_N0_downlink_dB, modulation)

        # Check if BER meets requirement (assuming a threshold, e.g., BER < 1e-5)
        if BER_downlink < 1e-5:
            # Select the modulation
            selected_downlink_modulation = modulation
            achieved_downlink_BER = BER_downlink
            break  # Select the first modulation that meets the data rate and BER requirement

    if selected_uplink_modulation is None or selected_downlink_modulation is None:
        return None  # Could not find suitable modulation

    # Convert data rates to Mbps
    data_rate_uplink_Mbps = desired_data_rate_uplink_bps / 1e6
    data_rate_downlink_Mbps = desired_data_rate_downlink_bps / 1e6

    return {
        "uplink_power_dBm": P_r_uplink_dBm,
        "downlink_power_dBm": P_r_downlink_dBm,
        "SNR_uplink_dB": SNR_uplink_dB,
        "SNR_downlink_dB": SNR_downlink_dB,
        "data_rate_uplink_Mbps": data_rate_uplink_Mbps,
        "data_rate_downlink_Mbps": data_rate_downlink_Mbps,
        "selected_uplink_modulation": selected_uplink_modulation,
        "selected_downlink_modulation": selected_downlink_modulation,
        "achieved_uplink_BER": achieved_uplink_BER,
        "achieved_downlink_BER": achieved_downlink_BER,
        "uplink_device": uplink_device,
        "satellite_transmitter": satellite_transmitter,
        "ground_station": ground_station,
        "satellite_tx_antenna": satellite_tx_antenna,
        "satellite_rx_antenna": satellite_rx_antenna,
    }

# Optimization function
def optimize_link_budget(orbital_params, ground_segment, satellite_params, devices, modulation_bits_per_symbol, antennas, desired_data_rate_uplink_bps, desired_data_rate_downlink_bps, ground_station_max_distances):
    """
    Iterate over all combinations to find the configuration that can achieve the desired data rates.
    """
    best_configuration = None

    # Ground station can always receive; no need for a downlink device
    transmit_devices = [device for device in devices if device["type"] in ["transmitter", "transceiver"]]

    for ground_station in ground_segment:
        # Get maximum distance
        max_distance = ground_station_max_distances[ground_station["name"]]

        if max_distance is None:
            continue  # Satellite is not visible from this ground station

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
                    satellite_params["ground_station_receive_gain_key"] = gain_key
                    break
        if not downlink_supported:
            continue

        # Assuming ground station transmit antenna gain is known (added to ground station data)
        if "G_ground_tx_dBi" not in ground_station:
            ground_station["G_ground_tx_dBi"] = 20  # Example value

        for uplink_device in transmit_devices:
            # Uplink device frequency check
            tx_freq_range = uplink_device.get("tx_frequency_range")
            if tx_freq_range and not (tx_freq_range[0] <= orbital_params["uplink_frequency"] <= tx_freq_range[1]):
                continue

            # Modulations check
            if not uplink_device["modulations"]:
                continue

            # Uplink transmit power is fixed (take maximum if multiple)
            if isinstance(uplink_device["transmit_power_W"], (list, tuple)):
                uplink_transmit_power_W = max(uplink_device["transmit_power_W"])
            else:
                uplink_transmit_power_W = uplink_device["transmit_power_W"]

            for satellite_transmitter in devices:
                if satellite_transmitter["type"] != "transmitter":
                    continue
                # Satellite transmitter frequency check
                freq_range = satellite_transmitter.get("frequency_range")
                if freq_range:
                    if isinstance(freq_range[0], tuple):
                        # Multiple frequency ranges
                        if not any(fr[0] <= orbital_params["downlink_frequency"] <= fr[1] for fr in freq_range):
                            continue
                    else:
                        if not (freq_range[0] <= orbital_params["downlink_frequency"] <= freq_range[1]):
                            continue

                # Modulations check
                if not satellite_transmitter["modulations"]:
                    continue

                possible_satellite_transmit_powers = get_possible_satellite_transmit_powers(satellite_transmitter["transmit_power_W"])

                for satellite_transmit_power_W in possible_satellite_transmit_powers:

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
                                satellite_params,
                                uplink_device,
                                satellite_transmitter,
                                modulation_bits_per_symbol,
                                sat_tx_ant,
                                sat_rx_ant,
                                desired_data_rate_uplink_bps,
                                desired_data_rate_downlink_bps,
                                satellite_transmit_power_W,
                                max_distance
                            )

                            if results is not None:
                                # Keep the configuration with the minimum satellite transmit power
                                if (best_configuration is None) or (satellite_transmit_power_W < best_configuration["satellite_transmit_power_W"]):
                                    best_configuration = results
                                    best_configuration["uplink_transmit_power_W"] = uplink_transmit_power_W
                                    best_configuration["satellite_transmit_power_W"] = satellite_transmit_power_W

        if best_configuration:
            break  # Found a suitable configuration with this ground station

    return best_configuration

# Main code
def main():
    # Load timescale and TLE data
    ts = load.timescale()

    # Select a satellite from satellite_tles
    satellite_key = "LEO"  # You can change this to "MEO", "GEO", etc.
    tle_lines = satellite_tles[satellite_key]["tle"]
    satellite = EarthSatellite(tle_lines[0], tle_lines[1], satellite_tles[satellite_key]["name"], ts)

    # Get user input for orbital parameters (excluding altitude)
    print("Enter Orbital Parameters:")
    try:
        uplink_frequency = float(input("Uplink Frequency (in MHz, e.g., 2025): "))
        downlink_frequency = float(input("Downlink Frequency (in MHz, e.g., 25600): "))
        T_system_ground = float(input("Ground Station System Noise Temperature (in Kelvin, e.g., 290): "))
        T_system_satellite = float(input("Satellite System Noise Temperature (in Kelvin, e.g., 290): "))
    except ValueError:
        print("Invalid input. Please enter numerical values.")
        exit()

    orbital_parameters = {
        "uplink_frequency": uplink_frequency,  # MHz
        "downlink_frequency": downlink_frequency,  # MHz
    }

    # Satellite parameters
    satellite_parameters = {
        "T_system_satellite": T_system_satellite,  # K
        "ground_station_receive_gain_key": None,  # To be set during optimization
    }

    # Update ground stations with their system noise temperature
    for gs in ground_segment:
        gs["T_system_ground"] = T_system_ground  # K
        # Assuming ground station transmit antenna gain is known
        if "G_ground_tx_dBi" not in gs:
            gs["G_ground_tx_dBi"] = 20  # Example value

    # Get user input for desired data rates
    try:
        desired_data_rate_uplink_Mbps = float(input("\nEnter the desired UPLINK data rate (in Mbps): "))
        desired_data_rate_downlink_Mbps = float(input("Enter the desired DOWNLINK data rate (in Mbps): "))
    except ValueError:
        print("Invalid input. Please enter numerical values.")
        exit()

    desired_data_rate_uplink_bps = desired_data_rate_uplink_Mbps * 1e6
    desired_data_rate_downlink_bps = desired_data_rate_downlink_Mbps * 1e6

    # Generate times over 24 hours
    time_steps = np.arange(0, 24 * 60, 1)  # every minute over 24 hours
    times = ts.utc(2023, 11, 12, time_steps // 60, time_steps % 60)

    # Compute maximum distances for each ground station
    ground_station_max_distances = {}
    for gs in ground_segment:
        # Get ground station location
        ground_loc = gs["location"]

        # Compute satellite's position as seen from ground station
        difference = satellite - ground_loc

        topocentric = difference.at(times)
        alt, az, distance = topocentric.altaz()

        # alt is the elevation angle in degrees
        # distance is in Astronomical Units

        # Convert distance to meters
        distance_m = distance.m

        # Get times when satellite is above horizon
        visible = alt.degrees > 10  # Use elevation angle > 10 degrees

        # For times when satellite is visible, compute the maximum distance
        if np.any(visible):
            max_distance = np.max(distance_m[visible])
            ground_station_max_distances[gs["name"]] = max_distance
        else:
            # Satellite is never visible from this ground station
            ground_station_max_distances[gs["name"]] = None

    # Run optimization
    best_config = optimize_link_budget(
        orbital_parameters,
        ground_segment,
        satellite_parameters,
        modems_sdrs,
        modulation_bits_per_symbol,
        antennas,
        desired_data_rate_uplink_bps,
        desired_data_rate_downlink_bps,
        ground_station_max_distances
    )

    if best_config:
        uplink_device = best_config["uplink_device"]
        satellite_transmitter = best_config["satellite_transmitter"]
        ground_station = best_config["ground_station"]
        ground_station_receive_gain_key = satellite_parameters["ground_station_receive_gain_key"]
        ground_station_receive_gain = ground_station[ground_station_receive_gain_key]
        print(f"\nOptimal Ground Station: {ground_station['name']}")
        print(f"Optimal Uplink Device: {uplink_device['name']}")
        print(f"Satellite Transmitter: {satellite_transmitter['name']}")
        print(f"Satellite Transmit Antenna: {best_config['satellite_tx_antenna']['name']} with Gain {best_config['satellite_tx_antenna']['gain']} dBi")
        print(f"Satellite Receive Antenna: {best_config['satellite_rx_antenna']['name']} with Gain {best_config['satellite_rx_antenna']['gain']} dBi")
        print(f"Ground Station Receive Antenna Gain ({ground_station_receive_gain_key}): {ground_station_receive_gain} dBi")
        print(f"Ground Station Transmit Antenna Gain: {ground_station['G_ground_tx_dBi']} dBi")
        print("\n=== Uplink (Ground to Satellite) ===")
        print(f"Uplink Transmit Power: {best_config['uplink_transmit_power_W']:.2f} W")
        print(f"Uplink Received Power at Satellite: {best_config['uplink_power_dBm']:.2f} dBm")
        print(f"Uplink SNR: {best_config['SNR_uplink_dB']:.2f} dB")
        print(f"Selected Uplink Modulation: {best_config['selected_uplink_modulation']}")
        print(f"Achieved Uplink Data Rate: {best_config['data_rate_uplink_Mbps']:.6f} Mbps")
        print(f"Achieved Uplink BER: {best_config['achieved_uplink_BER']:.2e}")
        print("\n=== Downlink (Satellite to Ground) ===")
        print(f"Satellite Transmit Power: {best_config['satellite_transmit_power_W']:.2f} W")
        print(f"Downlink Received Power at Ground Station: {best_config['downlink_power_dBm']:.2f} dBm")
        print(f"Downlink SNR: {best_config['SNR_downlink_dB']:.2f} dB")
        print(f"Selected Downlink Modulation: {best_config['selected_downlink_modulation']}")
        print(f"Achieved Downlink Data Rate: {best_config['data_rate_downlink_Mbps']:.6f} Mbps")
        print(f"Achieved Downlink BER: {best_config['achieved_downlink_BER']:.2e}")

        # Now analyze the link budget across 24 hours
        gs = best_config["ground_station"]
        ground_loc = gs["location"]

        # Compute satellite's position as seen from ground station
        difference = satellite - ground_loc

        topocentric = difference.at(times)
        alt, az, distance = topocentric.altaz()

        # Convert distance to meters
        distance_m = distance.m

        # Get times when satellite is above horizon
        visible = alt.degrees > 10  # Use elevation angle > 10 degrees

        link_budget_over_time = []

        for i in range(len(times)):
            if visible[i]:
                distance_m_i = distance_m[i]

                # Compute the link budget using the best configuration
                result = calculate_link_budget(
                    orbital_parameters,
                    gs,
                    satellite_parameters,
                    best_config["uplink_device"],
                    best_config["satellite_transmitter"],
                    modulation_bits_per_symbol,
                    best_config["satellite_tx_antenna"],
                    best_config["satellite_rx_antenna"],
                    desired_data_rate_uplink_bps,
                    desired_data_rate_downlink_bps,
                    best_config["satellite_transmit_power_W"],
                    distance_m_i
                )

                if result:
                    # Store the result along with the time
                    result["time"] = times[i].utc_datetime()
                    link_budget_over_time.append(result)

        # Output the link budget over time
        if link_budget_over_time:
            print("\nLink Budget Analysis Over 24 Hours:")
            for lb in link_budget_over_time:
                time_str = lb["time"].strftime("%Y-%m-%d %H:%M:%S")
                print(f"Time: {time_str}, Uplink SNR: {lb['SNR_uplink_dB']:.2f} dB, Downlink SNR: {lb['SNR_downlink_dB']:.2f} dB")
        else:
            print("\nSatellite is not visible from the ground station during the 24-hour period.")
    else:
        print("No suitable configuration found to meet the desired data rates.")

if __name__ == "__main__":
    main()
