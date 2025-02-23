import numpy as np
from scipy.special import erfc
from skyfield.api import Topos, EarthSatellite, load, wgs84
import datetime

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
        "modulations": ["FSK", "MSK", "GFSK", "GMSK"],
        "transmit_power_W": (0.4,2),  # Watts
        "receive_sensitivity_dBm": -121,
        "rx_frequency_range": (2025, 2110),  # MHz
        "tx_frequency_range": (2200, 2290),  # MHz
    },
    {
        "type": "transceiver",
        "name": "Flexible and Miniaturised Transceiver - GOMSpace",
        "data_rate": (0.1, 38.4),  # kbps
        "modulations": ["GFSK", "GMSK"],
        "transmit_power_W": 1,  # Watts
        "receive_sensitivity_dBm": -137,
        "rx_frequency_range": (430, 440),  # MHz
        "tx_frequency_range": (430, 440),  # MHz
    },
    {
        "type": "transceiver",
        "name": "NanoCom AX2150 - GOMSpace",
        "data_rate": (9.6, 96),  # kbps
        "modulations": ["GFSK", "GMSK"],
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
        "network": "VIASAT",
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
        "network": "VIASAT",
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
        "network": "VIASAT",
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
        "network": "VIASAT",
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
        "network": "ATLAS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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
        "network": "AWS",
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

# Constants
k = 1.38e-23  # Boltzmann's constant (J/K)
c = 3e8       # Speed of light (m/s)

# Function to calculate Free Space Path Loss (FSPL)
def calculate_fspl(distance_m, frequency_hz):
    """Calculate Free Space Path Loss (FSPL) in dB."""
    fspl = 20 * np.log10((4 * np.pi * distance_m * frequency_hz) / c)
    return fspl

# Function to compute BER given Eb/N0 and modulation
def compute_BER(Eb_N0_dB, modulation_scheme):
    Eb_N0_linear = 10 ** (Eb_N0_dB / 10)
    if modulation_scheme in ["BPSK", "OOK", "GMSK", "GFSK", "FSK", "MSK"]:
        BER = 0.5 * erfc(np.sqrt(Eb_N0_linear))
    elif modulation_scheme == "QPSK":
        BER = 0.5 * erfc(np.sqrt(Eb_N0_linear))
    elif modulation_scheme == "8PSK":
        BER = erfc(np.sqrt(1.5 * Eb_N0_linear * np.log2(8) / (8 - 1)))
    elif modulation_scheme in ["16QAM", "16APSK"]:
        BER = (3/8) * erfc(np.sqrt((4/5) * Eb_N0_linear))
    elif modulation_scheme == "32APSK":
        BER = erfc(np.sqrt(0.068 * Eb_N0_linear))
    elif modulation_scheme == "256APSK":
        BER = erfc(np.sqrt(0.0156 * Eb_N0_linear))
    else:
        BER = 1.0  # Unsupported modulation
    return BER

# Function to select the best antenna based on frequency
def select_antenna(frequency_hz):
    """Select the antenna with the highest gain that covers the given frequency."""
    suitable_antennas = [ant for ant in antennas if ant['frequency_range'][0] <= frequency_hz <= ant['frequency_range'][1]]
    if not suitable_antennas:
        return None
    # Select the antenna with the highest gain
    best_antenna = max(suitable_antennas, key=lambda ant: ant['gain'])
    return best_antenna

# Function to calculate link budget
def calculate_link_budget(params):
    """Calculate link budget and check if it meets BER requirements."""
    # Unpack parameters
    distance_m = params['distance_m']
    frequency_hz = params['frequency_hz']
    transmit_power_W = params['transmit_power_W']
    transmit_gain_dBi = params['transmit_gain_dBi']
    receive_gain_dBi = params['receive_gain_dBi']
    data_rate_bps = params['data_rate_bps']
    modulation = params['modulation']
    T_system = params['T_system']
    acceptable_BER = params['acceptable_BER']
    bits_per_symbol = modulation_bits_per_symbol.get(modulation, 1)

    # Calculate FSPL
    fspl_dB = calculate_fspl(distance_m, frequency_hz)

    # Convert transmit power to dBm
    P_tx_dBm = 10 * np.log10(transmit_power_W * 1e3)  # W to dBm

    # Calculate received power
    P_rx_dBm = P_tx_dBm + transmit_gain_dBi + receive_gain_dBi - fspl_dB

    # Calculate noise power
    bandwidth_Hz = data_rate_bps / bits_per_symbol
    P_noise_dBm = 10 * np.log10(k * T_system * bandwidth_Hz) + 30  # in dBm

    # Calculate SNR
    SNR_dB = P_rx_dBm - P_noise_dBm

    # Calculate Eb/N0
    Eb_N0_dB = SNR_dB - 10 * np.log10(bits_per_symbol)

    # Calculate BER
    BER = compute_BER(Eb_N0_dB, modulation)

    # Check if BER meets acceptable threshold
    if BER <= acceptable_BER:
        return {
            "P_rx_dBm": P_rx_dBm,
            "SNR_dB": SNR_dB,
            "Eb_N0_dB": Eb_N0_dB,
            "BER": BER
        }
    else:
        return None

# Function to find the best downlink configuration
def find_best_downlink_configuration(distance_m, gs, devices, constants):
    feasible_configurations = []
    for device in devices:
        if device['type'] in ['transmitter', 'transceiver']:
            device_tx_freq_range = device.get('tx_frequency_range', device.get('frequency_range', None))
            if not device_tx_freq_range:
                continue
            # Handle multiple frequency ranges
            if isinstance(device_tx_freq_range, list):
                freq_ranges = device_tx_freq_range
            else:
                freq_ranges = [device_tx_freq_range]
            for fr in freq_ranges:
                # Convert frequency range to Hz
                fr_hz = (fr[0]*1e6, fr[1]*1e6)
                # Determine which downlink band this frequency falls into
                downlink_band = None
                if gs['sdown_bw'] and gs['sdown_bw'][0]*1e6 <= fr_hz[0] <= gs['sdown_bw'][1]*1e6:
                    downlink_band = 'sdown_bw'
                    ground_receive_gain = gs['sdown_gt']
                elif gs['xdown_bw'] and gs['xdown_bw'][0]*1e6 <= fr_hz[0] <= gs['xdown_bw'][1]*1e6:
                    downlink_band = 'xdown_bw'
                    ground_receive_gain = gs['xdown_gt']
                elif gs['kadown_bw'] and gs['kadown_bw'][0]*1e6 <= fr_hz[0] <= gs['kadown_bw'][1]*1e6:
                    downlink_band = 'kadown_bw'
                    ground_receive_gain = gs['kadown_gt']
                else:
                    continue  # Frequency does not match any downlink band
                # Select the best antenna for the device's tx frequency
                antenna = select_antenna(fr_hz[0])
                if not antenna:
                    continue
                transmit_gain_dBi = antenna['gain']
                # Determine transmit power
                if 'transmit_power_W' in device:
                    if isinstance(device['transmit_power_W'], (list, tuple)):
                        transmit_power_W = max(device['transmit_power_W'])
                    else:
                        transmit_power_W = device['transmit_power_W']
                elif 'x_band_transmit_power_dBm' in device:
                    transmit_power_W = 10 ** (max(device['x_band_transmit_power_dBm']) / 10) / 1e3  # Convert dBm to W
                else:
                    continue
                # Iterate over modulations
                modulations = device.get('modulations', [])
                if not modulations:
                    continue
                for modulation in modulations:
                    # Determine data rate
                    if isinstance(device.get('data_rate', None), tuple):
                        data_rate_bps = device['data_rate'][1] * 1e3  # Convert kbps to bps
                    elif 'x_band_tx_data_rate' in device:
                        data_rate_bps = device['x_band_tx_data_rate'][1] * 1e3
                    elif 'data_rate' in device:
                        data_rate_bps = device['data_rate'] * 1e3
                    else:
                        continue
                    downlink_params = {
                        'distance_m': distance_m,
                        'frequency_hz': fr_hz[0],
                        'transmit_power_W': transmit_power_W,
                        'transmit_gain_dBi': transmit_gain_dBi,
                        'receive_gain_dBi': ground_receive_gain,
                        'data_rate_bps': data_rate_bps,
                        'modulation': modulation,
                        'T_system': constants['T_system_ground'],
                        'acceptable_BER': constants['acceptable_BER_downlink'],
                    }
                    # Compute link budget
                    result = calculate_link_budget(downlink_params)
                    if result:
                        # Link budget acceptable
                        feasible_configurations.append({
                            'device': device,
                            'modulation': modulation,
                            'data_rate_bps': data_rate_bps,
                            'result': result
                        })
    if feasible_configurations:
        # Find the configuration with the highest data rate
        best_config = max(feasible_configurations, key=lambda x: x['data_rate_bps'])
        return best_config
    else:
        return None

# Function to merge overlapping intervals
def merge_intervals(intervals):
    """Merge overlapping intervals."""
    intervals = sorted(intervals, key=lambda x: x[0].utc_datetime())
    merged = []
    for current in intervals:
        if not merged:
            merged.append(current)
        else:
            prev = merged[-1]
            if current[0].utc_datetime() <= prev[1].utc_datetime():
                # Overlapping intervals, merge them
                merged[-1] = (prev[0], max(prev[1], current[1]))
            else:
                merged.append(current)
    return merged

# Main function
def main():
    # Extract unique networks from ground_segment
    networks = sorted(list(set(gs['network'] for gs in ground_segment)))
    
    # User selects ground station network
    print("Available Ground Station Networks:")
    for i, network_name in enumerate(networks, start=1):
        print(f"{i}. {network_name}")
    try:
        network_choice = int(input("Select a Ground Station Network (enter number): "))
        if not 1 <= network_choice <= len(networks):
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter a valid number corresponding to the network.")
        return
    selected_network = networks[network_choice - 1]
    print(f"\nSelected Network: {selected_network}")
    
    # Gather all ground stations in the selected network
    selected_ground_stations = [gs for gs in ground_segment if gs['network'] == selected_network]
    
    # User selects orbit type
    print("\nAvailable Orbit Types:")
    orbit_types = list(satellite_tles.keys())
    for i, orbit_type in enumerate(orbit_types, start=1):
        print(f"{i}. {orbit_type}")
    try:
        orbit_choice = int(input("Select an Orbit Type (enter number): "))
        if not 1 <= orbit_choice <= len(orbit_types):
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter a valid number corresponding to the orbit type.")
        return
    selected_orbit = orbit_types[orbit_choice - 1]
    print(f"\nSelected Orbit Type: {selected_orbit}")
    
    # Load satellite TLE
    tle_data = satellite_tles[selected_orbit]['tle']
    ts = load.timescale()
    satellite = EarthSatellite(tle_data[0], tle_data[1], satellite_tles[selected_orbit]['name'], ts)
    
    # Assume both the ground station and satellite are at 290 K
    T_system_ground = 290  # K
    T_system_satellite = 290  # K
    
    # Acceptable BERs
    acceptable_BER_uplink = 1e-7
    acceptable_BER_downlink = 1e-5
    
    constants = {
        'T_system_ground': T_system_ground,
        'T_system_satellite': T_system_satellite,
        'acceptable_BER_uplink': acceptable_BER_uplink,
        'acceptable_BER_downlink': acceptable_BER_downlink,
    }
    
    # Simulation time over 30 days
    start_time = datetime.datetime.utcnow()
    end_time = start_time + datetime.timedelta(days=30)
    t0 = ts.utc(start_time.year, start_time.month, start_time.day)
    t1 = ts.utc(end_time.year, end_time.month, end_time.day)
    
    # List to store visibility intervals from all ground stations
    all_visibility_intervals = []
    total_data_transferred = 0.0  # in bits
    best_overall_config = None
    max_data_rate = 0.0  # in bps
    
    # Iterate over each ground station in the selected network
    for gs in selected_ground_stations:
        #print(f"\n--- Evaluating for Ground Station: {gs['name']} ---")
        ground_loc = gs['location']
        
        # Find satellite visibility events from the ground station
        t_events, events = satellite.find_events(ground_loc, t0, t1, altitude_degrees=10.0)
        
        # Collect rise and set times
        passes = []
        current_pass = {}
        for ti, event in zip(t_events, events):
            if event == 0:  # Rise
                current_pass = {'start': ti}
            elif event == 1:  # Culmination
                current_pass['culmination': ti]
            elif event == 2:  # Set
                current_pass['end'] = ti
                if 'start' in current_pass:
                    passes.append(current_pass)
                    all_visibility_intervals.append((current_pass['start'], current_pass['end']))
                current_pass = {}
        
        # Process each pass
        for p in passes:
            pass_start = p['start']
            pass_end = p['end']
            pass_duration_s = (pass_end.utc_datetime() - pass_start.utc_datetime()).total_seconds()
            # Get distance at culmination or midpoint
            pass_midpoint = pass_start.utc_datetime() + (pass_end.utc_datetime() - pass_start.utc_datetime()) / 2
            pass_midpoint_ts = ts.utc(pass_midpoint)
            difference = satellite - ground_loc
            topocentric = difference.at(pass_midpoint_ts)
            distance = topocentric.distance()
            distance_m = distance.m
            
            # Find the best downlink configuration
            best_config = find_best_downlink_configuration(distance_m, gs, modems_sdrs, constants)
            if best_config:
                data_transferred = best_config['data_rate_bps'] * pass_duration_s
                total_data_transferred += data_transferred
                # Check if this configuration has the highest data rate
                if best_config['data_rate_bps'] > max_data_rate:
                    max_data_rate = best_config['data_rate_bps']
                    best_overall_config = best_config
    
    # Merge all visibility intervals to compute total visibility time
    merged_intervals = merge_intervals(all_visibility_intervals)
    total_visibility_time_s = sum((interval[1].utc_datetime() - interval[0].utc_datetime()).total_seconds() for interval in merged_intervals)
    
    total_simulation_time_s = (end_time - start_time).total_seconds()
    visibility_percentage = (total_visibility_time_s / total_simulation_time_s) * 100.0
    
    # Convert total data transferred to GB
    total_data_transferred_GB = total_data_transferred / (8 * 1e9)  # Convert bits to GB
    
    # Output results
    print("\n--- Summary ---")
    print(f"Total Data Transferred: {total_data_transferred_GB:.2f} GB over 30 days")
    print(f"Percentage of time satellite was in view: {visibility_percentage:.2f}%")
    if best_overall_config:
        device = best_overall_config['device']
        print("\nBest RF Configuration:")
        print(f"Device Name: {device['name']}")
        print(f"Modulation: {best_overall_config['modulation']}")
        print(f"Data Rate: {best_overall_config['data_rate_bps'] / 1e6:.2f} Mbps")
    else:
        print("No suitable RF configuration found.")

if __name__ == "__main__":
    main()