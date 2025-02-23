import numpy as np
from scipy.special import erfc
from skyfield.api import Topos, EarthSatellite, load
import datetime

antennas = [
    # VHF (30 MHz - 300 MHz)
    {
        "name": "FAKE VHF Yagi Antenna",
        "frequency_range": (30e6, 300e6),
        "gain": 7,
        "type": "Yagi-Uda"
    },

    # UHF (300 MHz - 1 GHz)
    {
        "name": "FAKE UHF Log-Periodic Dipole Array",
        "frequency_range": (300e6, 1e9),
        "gain": 10,
        "type": "Log-Periodic"
    },
    {
        "name": "FAKE UHF Patch Antenna",
        "frequency_range": (300e6, 1e9),
        "gain": 8,
        "type": "Patch"
    },

    # L-Band (1 GHz - 2 GHz)
    {
        "name": "FAKE L-Band Helical Antenna",
        "frequency_range": (1e9, 2e9),
        "gain": 12,
        "type": "Helical"
    },
    {
        "name": "FALCCON - RW",
        "frequency_range": (1e9, 2e9),
        "gain": 21,
        "type": "Dipole Array"
    },
    {
        "name": "L Band Patch Antenna - Printech",
        "frequency_range": (1.563e9, 1.587e9),
        "gain": 5,
        "type": "Patch"
    },

    # S-Band (2 GHz - 4 GHz)
    {
        "name": "AC-2000 - AAC",
        "frequency_range": (2e9, 2.3e9),
        "gain": 5.2,
        "type": "Parabolic Dish",
        "vswr": {
            "frequency": [2e9, 2.3e9],  # Frequencies in Hz
            "vswr_values": [1.5, 1.5]  # VSWR values
            }
    },
    {
        "name": "SANT S-band Patch Antenna - AAC",
        "frequency_range": (2.2e9, 2.29e9),
        "gain": 7,
        "type": "Patch",
        "connector": "SMA_P",
        "s11_dB": -15
    },
    {
        "name": "Quad S Band Antenna - IQ Tech",
        "frequency_range": (1.98e9, 2.5e9),
        "gain": 11,
        "type": "Patch Array",
        "vswr": {
            "frequency": [1.98e9, 2.5e9],  # Frequencies in Hz
            "vswr_values": [1.8, 1.8]  # VSWR values, check data sheet
            }
    },
    {
        "name": "S-Band Antenna Commercial - Enduro",
        "frequency_range": (2.025e9, 2.11e9),
        "gain": 7,
        "type": "Patch array",
        "vswr": {
            "frequency": [2.025e9, 2.11e9],  # Frequencies in Hz
            "vswr_values": [1.8, 1.8]  # VSWR values, check data sheet
            }
    },
    {
        "name": "S-Band Antenna Wideband - Enduro",
        "frequency_range": (2.2e9, 2.29e9),
        "gain": 5,
        "type": "Patch array",
        "vswr": {
            "frequency": [2.025e9, 2.11e9],  # Frequencies in Hz
            "vswr_values": [1.8, 1.8]  # VSWR values, check data sheet
            }
    },
    # C-Band (4 GHz - 8 GHz)
    {
        "name": "FAKE C-Band Horn Antenna",
        "frequency_range": (4e9, 8e9),
        "gain": 18,
        "type": "Horn"
    },

    # X-Band (8 GHz - 12 GHz)
    {
        "name": "X-Band Patch Antenna - Enduro",
        "frequency_range": (8.025e9, 8.4e9),
        "gain": 6,
        "type": "Patch"
    },
    {
        "name": "4x4 X-Band Patch Array - Enduro",
        "frequency_range": (8.025e9, 8.4e9),
        "gain": 16,
        "type": "Patch Array"
    },
    {
         "name": "XANT X-Band Patch Antenna - Cubecom",
         "frequency_range": [8e9, 8.4e9],  # Frequency range in Hz
         "gain": 8,
         "type": "Patch",
         "return_loss": {
             "frequency": [8.0e9, 8.1e9, 8.2e9, 8.3e9, 8.4e9],  # Frequencies in Hz
             "return_loss_dB": [-21, -20, -18, -16, -17]  # Return loss in dB
        }
    },
    {
        "name": "XPLANT X-band Payload Antenna - Cubecom",
        "frequency_range": (8.5e9, 9.6e9),
        "gain": 8,
        "type": "Patch Array",
        "return_loss": {
            "frequency": [8.50e9, 8.75e9, 8.90e9, 9.00e9, 9.25e9, 9.50e9],
            "return_loss_dB": [-10, -20, -15, -25, -30, -10]
        }
    },
    {
        "name": "High Gain X-Band Antenna - Anywaves",
        "frequency_range": (7.9e9, 8.5e9),
        "gain": 15.5,
        "type": "Patch array",
        "return_loss": {
            "frequency": [7.9e9,8.5e9],
            "return_loss_dB": [-10, -10] #made up, no data
        }
    },
    {
    "name": "High-Gain X-Band Patch Array - Printech",
    "frequency_range": [7.5e9, 8.5e9],  # Frequency range in Hz
    "gain": 20.7,  # Gain in dB
    "type": "Patch Array",
    "vswr": {
        "frequency": [7.5e9, 7.6e9, 7.7e9, 7.8e9, 7.9e9, 8.0e9, 8.1e9, 8.2e9, 8.3e9, 8.4e9],  # Frequencies in Hz
        "vswr_values": [1.3, 1.4, 1.2, 1.5, 1.3, 1.6, 1.2, 1.8, 1.4, 2.5]  # VSWR values
        }
    },
    {
        "name": "Lens Horn Antenna - Anteral",
        "frequency_range": (8.2e9, 12.4e9),
        "gain": 30.4,
        "type": "Lens Horn",
        "s11_dB": -18
    },  # not trusted?
    # Ku-Band (12 GHz - 18 GHz)
    {
        "name": "FAKE Ku-Band Microstrip Antenna",
        "frequency_range": (12e9, 18e9),
        "gain": 23,
        "type": "Microstrip"
    },

    # K-Band (18 GHz - 27 GHz)
    {
        "name": "FAKE K-Band Waveguide Antenna",
        "frequency_range": (18e9, 27e9),
        "gain": 25,
        "type": "Waveguide"
    },
    {
        "name": "K-Band 4x4 Patch Array - Enduro",
        "frequency_range": (17.7e9, 20.2e9),
        "gain": 16,
        "type": "Waveguide"
    },

    # Ka-Band (27 GHz - 40 GHz)
    {
        "name": "FAKE Ka-Band Horn Antenna",
        "frequency_range": (27e9, 40e9),
        "gain": 30,
        "type": "Horn"
    },
    {
        "name": "PAN-5151-64-KA - ReliaSat",
        "frequency_range": (27e9, 31e9),
        "gain": 20,
        "type": "Panel Array"
    },
    {
        "name": "4x4 X-Band Patch Array - Enduro",
        "frequency_range": (8.025e9, 8.4e9),
        "gain": 16,
        "type": "Patch Array"
    },

    # Above 40 GHz (EHF - Millimeter Wave Frequencies)
    {
        "name": "FAKE EHF Lens Antenna",
        "frequency_range": (40e9, 60e9),
        "gain": 35,
        "type": "Lens"
    },
    {
        "name": "FAKE Terahertz Horn Antenna",
        "frequency_range": (60e9, 100e9),
        "gain": 40,
        "type": "Horn"
    }
]


def interpolate_return_loss(frequency, antenna_data):
    freq_points = antenna_data["return_loss"]["frequency"]
    loss_values = antenna_data["return_loss"]["values"]
    return np.interp(frequency, freq_points, loss_values)
     
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
        "type": "receiver",
        "name": "RX-2000 S-Band Receiver - AAC",
        "data_rate": (9.6, 153.6),  # kbps
        "modulations": ["FM", "GFSK"],
        "receive_sensitivity_dBm": -110, #-117 @ 9.6kbps
        "rx_frequency_range": (2000e6, 2400e6),  # Hz
        "rx_ant_con": "SMA",
        "interface": "Micro-D",
    },
    {
        "type": "transceiver",
        "name": "UHF Transceiver II - Enduro",
        "data_rate": (0.1, 19.2),  # kbps
        "modulations": ["OOK", "GMSK", "2FSK", "4FSK", "4GFSK"],
        "transmit_power_W": [1, 2],  # Watts
        "receive_sensitivity_dBm": -121,
        "rx_frequency_range": (400e6, 403e6),  # Hz
        "tx_frequency_range": (430e6, 440e6),  # Hz
        "rx_ant_con": "SMA",
        "tx_ant_con": "SMA",
        "interface": ["RS485", "UART", "I2C", "USB-C"]
    },
    {
        "type": "transceiver",
        "name": "S Band Transceiver - Enduro",
        "data_rate": (0.1, 125),  # kbps
        "modulations": ["FSK", "MSK", "GFSK", "GMSK"],
        "transmit_power_W": (0.4, 2),  # Watts
        "receive_sensitivity_dBm": -121,
        "rx_frequency_range": (2025e6, 2110e6),  # Hz
        "tx_frequency_range": (2200e6, 2290e6),  # Hz
        "rx_ant_con": "SMP",
        "tx_ant_con": "SMP",
        "interface": ["RS-485", "RS-485/422", "USB-C"],
    },
    {
        "type": "transceiver",
        "name": "AX100 - GOMSpace",
        "data_rate": (0.1, 38.4),  # kbps
        "modulations": ["GFSK", "GMSK"],
        "transmit_power_W": 1,  # Watts
        "receive_sensitivity_dBm": -137,
        "rx_frequency_range": (430e6, 440e6),  # Hz
        "tx_frequency_range": (430e6, 440e6),  # Hz
        "rx_ant_con": "MCX",
        "interface": "CSP",
    },
    {
        "type": "transceiver",
        "name": "NanoCom AX2150 - GOMSpace",
        "data_rate": (2.4, 90),  # kbps
        "modulations": ["GFSK", "GMSK"],
        "transmit_power_W": (0.008, 0.5),  # Watts
        "receive_sensitivity_dBm": -113,
        "rx_frequency_range": (2025e6, 2110e6),  # Hz
        "tx_frequency_range": (2200e6, 2290e6),  # Hz
        "rx_ant_con": "SMP",
        "interface": "CSP"
    },
    {
        "type": "transceiver",
        "name": "TOTEM SDR - Alen Space",
        "data_rate": (200,56000),  # kbps
        "modulations": ["GFSK", "GMSK"],
        "transmit_power_W": (0.1,3),  # Watts
        "receive_sensitivity_dBm": -89, #AD9364 transceiver data sheet, assumptions
        "rx_frequency_range": (70e6, 60000e6),  # Hz
        "tx_frequency_range": (70e6, 60000e6),  # Hz
        "rx_ant_con": "MMCX",
        "tx_ant_con": "MMCX", 
        "interface": ["UART", "I2C", "JTAG", "ETHERNET", "CAN"],
    },
    {
        "type": "transmitter",
        "name": "S Band Transmitter - Enduro",
        "data_rate": 20000,  # ksps - check data sheet
        "modulations": ["QPSK", "8PSK", "16APSK"],
        "transmit_power_W": (0.5, 2),
        "frequency_range": [(2200e6, 2290e6), (2400e6, 2450e6)],
        "tx_ant_con": "SMA",
        "interface": ["UART", "RS-485", "LVDS"],
    },
    {
        "type": "transmitter",
        "name": "X Band Transmitter - Enduro",
        "data_rate": 150000,  # kbps
        "modulations": ["QPSK", "8PSK", "16APSK", "32APSK"],
        "transmit_power_W": (0.5, 2),
        "frequency_range": [(7900e6, 8400e6)],
        "tx_ant_con": "SMA",
        "interface": ["UART", "RS-485", "LVDS"],
    },
    {
        "type": "transmitter",
        "name": "K Band Transmitter - Enduro",
        "data_rate": 1000000,  # kbps
        "modulations": ["QPSK", "8PSK", "16APSK", "32APSK", "256APSK"],
        "transmit_power_W": (0.5, 2),
        "frequency_range": [(25500e6, 27000e6)],
        "tx_ant_con": "K-connector",
        "interface": ["CAN", "Ethernet", "ESPS-RS-485", "LVDS"],
    },
    {
        "type": "transmitter",
        "name": "XTX X-Band Transmitter - Cubecom",
        "data_rate": (2500,25000),  # kbps - check symbol rate data sheet
        "modulations": ["QPSK", "8PSK", "16APSK"],
        "transmit_power_W": (0, 2),
        "frequency_range": [(8025e6, 8400e6)],
        "tx_ant_con": "SMP",
        "interface": ["CAN", "I2C", "SpaceWire", "LVDS"],
    },
    {
        "type": "transmitter",
        "name": "HDRTX X-Band Gigabit Transmitter - Cubecom",
        "data_rate": (50000,200000),  # kbps - check symbol rate data sheet
        "modulations": ["8PSK", "16APSK", "32APSK"],
        "transmit_power_W": (0, 2),
        "frequency_range": [(8025e6, 8400e6)],
        "tx_ant_con": "SMP",
        "interface": ["CAN", "SpaceWire", "8B10B"],
    },
    {
        "type": "transmitter",
        "name": "TX-2400 S-Band Transmitter - AAC",
        "data_rate": (56,6000),  # kbps 
        "modulations": ["FM", "FSK"],
        "transmit_power_W": (1,10), # Different models
        "frequency_range": [(2000e6, 2400e6)],
        "tx_ant_con": "SMA",
        "interface": "Micro-D",
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

# connectors, p dictates precision.
rf_connectors = {
    "SMA": {
        "frequency_range": (0,18e9), #Hz
        "impedance": 50,
        "gender": {
            "male": {"contact": "pin", "thread_type": "outer"},
            "female": {"contact": "socket", "thread_type": "inner"}
        },
        "power_handling": "0.5 W (average)",
        "compatible_with": ["3.5 mm", "2.92 mm (K-connector)"],
    },
    
    "SMA_P": {
        "frequency_range": (0,26.5e9), #Hz
        "impedance": 50,
        "gender": {
            "male": {"contact": "pin", "thread_type": "outer"},
            "female": {"contact": "socket", "thread_type": "inner"}
        },
        "power_handling": "0.5 W (average)",
        "compatible_with": ["3.5 mm", "2.92 mm (K-connector)"],
    },
    
    "N-Type": {
        "frequency_range": (0,11e9), #Hz
        "impedance": 50,
        "gender": {
            "male": {"contact": "pin", "thread_type": "outer"},
            "female": {"contact": "socket", "thread_type": "inner"}
        },
        "power_handling": "150 W (average)",
        "compatible_with": ["Weatherproof N-Type"],
    },
    
    "N-Type_P": {
        "frequency_range": (0,18e9), #Hz
        "impedance": 50,
        "gender": {
            "male": {"contact": "pin", "thread_type": "outer"},
            "female": {"contact": "socket", "thread_type": "inner"}
        },
        "power_handling": "150 W (average)",
        "compatible_with": ["Weatherproof N-Type"],
    },
    "Micro-D": {
        "frequency_range": (0,3e9), #Hz
        "pins_sockets": {
            "9-pin": {"type": "male/female"},
            "15-pin": {"type": "male/female"},
            "25-pin": {"type": "male/female"}
        },
        "applications": [
            "Aerospace and defense",
            "Satellite communication",
            "High-reliability systems"
        ],
        "mounting": ["Panel mount", "Cable mount"],
        "notes": "Designed for compact, high-reliability connections."
    },
    "2.92 mm (K-Connector)": {
        "frequency_range": "DC to 40 GHz",
        "impedance": 50,
        "gender": {
            "male": {"contact": "pin", "precision": "high"},
            "female": {"contact": "socket", "precision": "high"}
        },
        "applications": [
            "Precision measurements",
            "High-frequency radar",
            "Satellite payload testing"
        ],
        "compatible_with": ["SMA", "3.5 mm"],
        "notes": "Provides excellent performance at high frequencies and is compatible with SMA and 3.5 mm connectors."
    }
}

# Ground stations with Skyfield Topos and additional properties
ground_segment = [
    {
        "name": "VIASAT PENDER",
        "network": "VIASAT",
        "location": Topos(latitude_degrees=49.1, longitude_degrees=-123.9, elevation_m=30),
        "sup_freq": (2025e6, 2110e6),
        "uEIRP": 53.2,  # dBW
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 17,  # dB/K
        "xdown_fr": (8025e6, 8400e6),
        "xdown_gt": 30,  # dB/K
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "VIASAT GUILDFORD",
        "network": "VIASAT",
        "location": Topos(latitude_degrees=51.2, longitude_degrees=-0.6, elevation_m=70),
        "sup_freq": (2025e6, 2110e6),
        "uEIRP": 53.2,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 17,
        "xdown_fr": (8025e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "VIASAT ALICE",
        "network": "VIASAT",
        "location": Topos(latitude_degrees=-23.7, longitude_degrees=133.9, elevation_m=600),
        "sup_freq": (2025e6, 2110e6),
        "uEIRP": 65.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (8025e6, 8400e6),
        "xdown_gt": 32,
        "kadown_fr": (25500e6, 27000e6),
        "kadown_gt": 34.5,
    },
    {
        "name": "VIASAT GHANA",
        "network": "VIASAT",
        "location": Topos(latitude_degrees=5.6, longitude_degrees=-0.2, elevation_m=50),  # Placeholder coordinates
        "sup_freq": (2025e6, 2110e6),
        "uEIRP": 65.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (8025e6, 8400e6),
        "xdown_gt": 32,
        "kadown_fr": (25500e6, 27000e6),
        "kadown_gt": 34.5,
    },
    {
        "name": "ATLAS PAUMALU",
        "network": "ATLAS",
        "location": Topos(latitude_degrees=21.6, longitude_degrees=-158.0, elevation_m=100),  # Placeholder coordinates
        "sup_freq": (2025e6, 2120e6),
        "uEIRP": 50.0,
        "sdown_fr": (2200e6, 2300e6),
        "sdown_gt": 21,
        "xdown_fr": (7900e6, 8500e6),
        "xdown_gt": 31,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Alaska 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=64.2008, longitude_degrees=-149.4937, elevation_m=100),  # Placeholder coordinates for Alaska
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_fr": (2200e6, 2290e6),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_fr": (7750e6, 8400e6),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Bahrain 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=26.0667, longitude_degrees=50.5577, elevation_m=50),  # Placeholder coordinates for Bahrain
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_fr": (2200e6, 2290e6),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_fr": (7750e6, 8400e6),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Cape Town 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=-33.9249, longitude_degrees=18.4241, elevation_m=50),  # Placeholder coordinates for Cape Town
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_fr": (2200e6, 2290e6),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_fr": (7750e6, 8400e6),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Dubbo 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=-32.2569, longitude_degrees=148.6011, elevation_m=50),  # Placeholder coordinates for Dubbo
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_fr": (2200e6, 2290e6),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_fr": (7750e6, 8400e6),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Hawaii 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=19.8968, longitude_degrees=-155.5828, elevation_m=100),  # Placeholder coordinates for Hawaii
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,  # Assumed value based on capabilities
        "sdown_fr": (2200e6, 2290e6),  # S-band downlink
        "sdown_gt": 18,  # Placeholder value
        "xdown_fr": (7750e6, 8400e6),  # X-band downlink
        "xdown_gt": 30,  # Placeholder value
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Ireland 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=53.1424, longitude_degrees=-7.6921, elevation_m=50),  # Placeholder coordinates for Ireland
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Ohio 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=40.4173, longitude_degrees=-82.9071, elevation_m=50),  # Placeholder coordinates for Ohio
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Oregon 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=43.8041, longitude_degrees=-120.5542, elevation_m=100),  # Placeholder coordinates for Oregon
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Punta Arenas 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=-53.1638, longitude_degrees=-70.9171, elevation_m=50),  # Placeholder coordinates for Punta Arenas
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Seoul 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=37.5665, longitude_degrees=126.9780, elevation_m=50),  # Placeholder coordinates for Seoul
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Singapore 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=1.3521, longitude_degrees=103.8198, elevation_m=50),  # Placeholder coordinates for Singapore
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
        "kadown_gt": None,
    },
    {
        "name": "AWS Stockholm 1",
        "network": "AWS",
        "location": Topos(latitude_degrees=59.3293, longitude_degrees=18.0686, elevation_m=50),  # Placeholder coordinates for Stockholm
        "sup_freq": (2025e6, 2110e6),  # S-band uplink
        "uEIRP": 53.0,
        "sdown_fr": (2200e6, 2290e6),
        "sdown_gt": 18,
        "xdown_fr": (7750e6, 8400e6),
        "xdown_gt": 30,
        "kadown_fr": None,
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

    # Calculate Link Margin
    link_margin_dB = SNR_dB - (10 * np.log10(bits_per_symbol) + acceptable_BER)

    # Calculate BER
    BER = compute_BER(Eb_N0_dB, modulation)

    # Calculate Spectral Efficiency
    spectral_efficiency = bits_per_symbol / bandwidth_Hz  # bits/sec/Hz

    # Calculate Energy per Bit
    energy_per_bit = transmit_power_W / data_rate_bps  # Joules/bit

    # Check if BER meets acceptable threshold
    if BER <= acceptable_BER:
        return {
            "P_rx_dBm": P_rx_dBm,
            "SNR_dB": SNR_dB,
            "Eb_N0_dB": Eb_N0_dB,
            "BER": BER,
            "link_margin_dB": link_margin_dB,
            "spectral_efficiency": spectral_efficiency,
            "energy_per_bit": energy_per_bit,
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
            freq_ranges = []
            if isinstance(device_tx_freq_range, list):
                for fr in device_tx_freq_range:
                    if isinstance(fr, tuple):
                        freq_ranges.append(fr)
                    else:
                        freq_ranges.append((fr, fr))
            elif isinstance(device_tx_freq_range, tuple):
                freq_ranges.append(device_tx_freq_range)
            else:
                freq_ranges.append((device_tx_freq_range, device_tx_freq_range))
            for fr in freq_ranges:
                # Determine which downlink band this frequency falls into
                downlink_band = None
                if gs['sdown_fr'] and gs['sdown_fr'][0] <= fr[0] <= gs['sdown_fr'][1]:
                    downlink_band = 'sdown_fr'
                    ground_receive_gain = gs['sdown_gt']
                elif gs['xdown_fr'] and gs['xdown_fr'][0] <= fr[0] <= gs['xdown_fr'][1]:
                    downlink_band = 'xdown_fr'
                    ground_receive_gain = gs['xdown_gt']
                elif gs['kadown_fr'] and gs['kadown_fr'][0] <= fr[0] <= gs['kadown_fr'][1]:
                    downlink_band = 'kadown_fr'
                    ground_receive_gain = gs['kadown_gt']
                else:
                    continue  # Frequency does not match any downlink band
                # Select the best antenna for the device's tx frequency
                tx_antenna = select_antenna(fr[0])
                if not tx_antenna:
                    continue
                transmit_gain_dBi = tx_antenna['gain']
                # Determine transmit power
                if 'transmit_power_W' in device:
                    if isinstance(device['transmit_power_W'], (list, tuple)):
                        transmit_power_W = max(device['transmit_power_W'])
                    else:
                        transmit_power_W = device['transmit_power_W']
                else:
                    continue
                # Iterate over modulations
                modulations = device.get('modulations', [])
                if not modulations:
                    continue
                for modulation in modulations:
                    # Determine data rate
                    if isinstance(device.get('data_rate', None), tuple):
                        data_rate_bps = device['data_rate'][1] * 1e3  # Convert kbps to bps (max data rate)
                    elif 'data_rate' in device:
                        data_rate_bps = device['data_rate'] * 1e3
                    else:
                        continue
                    downlink_params = {
                        'distance_m': distance_m,
                        'frequency_hz': fr[0],
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
                            'result': result,
                            'tx_antenna': tx_antenna,
                            'rx_antenna': None,  # Ground station antenna
                            'distance_m': distance_m,
                            'pass_duration_s': constants.get('pass_duration_s', 0),
                        })
    if feasible_configurations:
        # Find the configuration with the highest data rate
        best_config = max(feasible_configurations, key=lambda x: x['data_rate_bps'])
        return best_config
    else:
        return None

# Function to find the best uplink configuration
def find_best_uplink_configuration(distance_m, gs, devices, constants):
    feasible_configurations = []
    for device in devices:
        if device['type'] in ['receiver', 'transceiver']:
            device_rx_freq_range = device.get('rx_frequency_range', None)
            if not device_rx_freq_range:
                continue
            # Handle multiple frequency ranges
            freq_ranges = []
            if isinstance(device_rx_freq_range, list):
                for fr in device_rx_freq_range:
                    if isinstance(fr, tuple):
                        freq_ranges.append(fr)
                    else:
                        freq_ranges.append((fr, fr))
            elif isinstance(device_rx_freq_range, tuple):
                freq_ranges.append(device_rx_freq_range)
            else:
                freq_ranges.append((device_rx_freq_range, device_rx_freq_range))
            for fr in freq_ranges:
                # Determine which uplink band this frequency falls into
                uplink_band = None
                if gs['sup_freq'] and gs['sup_freq'][0] <= fr[0] <= gs['sup_freq'][1]:
                    uplink_band = 'sup_freq'
                    ground_transmit_EIRP_dBW = gs['uEIRP']
                else:
                    continue  # Frequency does not match any uplink band
                # Select the best antenna for the device's rx frequency
                rx_antenna = select_antenna(fr[0])
                if not rx_antenna:
                    continue
                receive_gain_dBi = rx_antenna['gain']
                # Ground station transmit power and gain
                transmit_gain_dBi = gs['uEIRP'] - 10 * np.log10(1e3)  # Convert dBW to dBi assuming 1 W transmit power
                transmit_power_W = 1  # Assume 1 W for simplicity
                # Iterate over modulations
                modulations = device.get('modulations', [])
                if not modulations:
                    continue
                for modulation in modulations:
                    # Determine data rate
                    if isinstance(device.get('data_rate', None), tuple):
                        data_rate_bps = device['data_rate'][1] * 1e3  # Convert kbps to bps (max data rate)
                    elif 'data_rate' in device:
                        data_rate_bps = device['data_rate'] * 1e3
                    else:
                        continue
                    uplink_params = {
                        'distance_m': distance_m,
                        'frequency_hz': fr[0],
                        'transmit_power_W': transmit_power_W,
                        'transmit_gain_dBi': transmit_gain_dBi,
                        'receive_gain_dBi': receive_gain_dBi,
                        'data_rate_bps': data_rate_bps,
                        'modulation': modulation,
                        'T_system': constants['T_system_satellite'],
                        'acceptable_BER': constants['acceptable_BER_uplink'],
                    }
                    # Compute link budget
                    result = calculate_link_budget(uplink_params)
                    if result:
                        # Link budget acceptable
                        feasible_configurations.append({
                            'device': device,
                            'modulation': modulation,
                            'data_rate_bps': data_rate_bps,
                            'result': result,
                            'rx_antenna': rx_antenna,
                            'tx_antenna': None,  # Ground station antenna
                            'distance_m': distance_m,
                            'pass_duration_s': constants.get('pass_duration_s', 0),
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
    total_data_transferred_downlink = 0.0  # in bits
    total_data_transferred_uplink = 0.0  # in bits
    best_overall_config_downlink = None
    best_overall_config_uplink = None
    max_data_rate_downlink = 0.0  # in bps
    max_data_rate_uplink = 0.0  # in bps
    total_SNR_downlink = []
    total_SNR_uplink = []
    max_link_distance = 0.0  # in meters

    # Iterate over each ground station in the selected network
    for gs in selected_ground_stations:
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
                current_pass['culmination'] = ti
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
            constants['pass_duration_s'] = pass_duration_s

            # Sample distances at multiple points during the pass
            num_samples = 10
            sample_times = [pass_start.utc_datetime() + i * (pass_end.utc_datetime() - pass_start.utc_datetime()) / num_samples for i in range(num_samples + 1)]
            distances_m = []
            for sample_time in sample_times:
                sample_ts = ts.utc(sample_time)
                difference = satellite - ground_loc
                topocentric = difference.at(sample_ts)
                distance = topocentric.distance()
                distances_m.append(distance.m)

            average_distance_m = sum(distances_m) / len(distances_m)
            max_distance_m = max(max_link_distance, max(distances_m))
            max_link_distance = max_distance_m

            # Collect maximum and average distances
            distance_m = average_distance_m

            # Find the best downlink configuration
            best_downlink_config = find_best_downlink_configuration(distance_m, gs, modems_sdrs, constants)
            if best_downlink_config:
                data_transferred = best_downlink_config['data_rate_bps'] * pass_duration_s
                total_data_transferred_downlink += data_transferred
                total_SNR_downlink.append(best_downlink_config['result']['SNR_dB'])
                # Check if this configuration has the highest data rate
                if best_downlink_config['data_rate_bps'] > max_data_rate_downlink:
                    max_data_rate_downlink = best_downlink_config['data_rate_bps']
                    best_overall_config_downlink = best_downlink_config

            # Find the best uplink configuration
            best_uplink_config = find_best_uplink_configuration(distance_m, gs, modems_sdrs, constants)
            if best_uplink_config:
                data_transferred = best_uplink_config['data_rate_bps'] * pass_duration_s
                total_data_transferred_uplink += data_transferred
                total_SNR_uplink.append(best_uplink_config['result']['SNR_dB'])
                # Check if this configuration has the highest data rate
                if best_uplink_config['data_rate_bps'] > max_data_rate_uplink:
                    max_data_rate_uplink = best_uplink_config['data_rate_bps']
                    best_overall_config_uplink = best_uplink_config

    # Merge all visibility intervals to compute total visibility time
    merged_intervals = merge_intervals(all_visibility_intervals)
    total_visibility_time_s = sum((interval[1].utc_datetime() - interval[0].utc_datetime()).total_seconds() for interval in merged_intervals)

    total_simulation_time_s = (end_time - start_time).total_seconds()
    visibility_percentage = (total_visibility_time_s / total_simulation_time_s) * 100.0

    # Convert total data transferred to GB
    total_data_transferred_downlink_GB = total_data_transferred_downlink / (8 * 1e9)  # Convert bits to GB
    total_data_transferred_uplink_GB = total_data_transferred_uplink / (8 * 1e9)  # Convert bits to GB

    # Calculate average SNRs
    average_SNR_downlink = sum(total_SNR_downlink) / len(total_SNR_downlink) if total_SNR_downlink else None
    average_SNR_uplink = sum(total_SNR_uplink) / len(total_SNR_uplink) if total_SNR_uplink else None

    # Calculate average effective data rates
    average_data_rate_downlink = (total_data_transferred_downlink / total_visibility_time_s) if total_visibility_time_s > 0 else 0
    average_data_rate_uplink = (total_data_transferred_uplink / total_visibility_time_s) if total_visibility_time_s > 0 else 0

    # Assume latency is the time it takes for a signal to travel the maximum link distance
    average_latency = max_link_distance / c  # in seconds

    # Output results
    print("\n--- Summary ---")
    print(f"Total Downlink Data Transferred: {total_data_transferred_downlink_GB:.2f} GB over 30 days")
    print(f"Total Uplink Data Transferred: {total_data_transferred_uplink_GB:.2f} GB over 30 days")
    print(f"Percentage of time satellite was in view: {visibility_percentage:.2f}%")
    if average_SNR_downlink is not None:
        print(f"Average Downlink SNR: {average_SNR_downlink:.2f} dB")
    else:
        print("Average Downlink SNR: N/A")
    if average_SNR_uplink is not None:
        print(f"Average Uplink SNR: {average_SNR_uplink:.2f} dB")
    else:
        print("Average Uplink SNR: N/A")
    print(f"Average Downlink Data Rate: {average_data_rate_downlink / 1e6:.2f} Mbps")
    print(f"Average Uplink Data Rate: {average_data_rate_uplink / 1e6:.2f} Mbps")
    print(f"Average Latency: {average_latency * 1e3:.2f} ms")
    print(f"Maximum Link Distance: {max_link_distance / 1e3:.2f} km")

    if best_overall_config_downlink:
        device = best_overall_config_downlink['device']
        tx_antenna = best_overall_config_downlink['tx_antenna']
        result = best_overall_config_downlink['result']
        print("\nBest Downlink RF Configuration:")
        print(f"Device Type: {device['type']}")
        print(f"Device Name: {device['name']}")
        print(f"Transmit Antenna: {tx_antenna['name']} (Gain: {tx_antenna['gain']} dBi)")
        print(f"Modulation: {best_overall_config_downlink['modulation']}")
        print(f"Maximum Achievable Data Rate: {best_overall_config_downlink['data_rate_bps'] / 1e6:.2f} Mbps")
        print(f"SNR: {result['SNR_dB']:.2f} dB")
        print(f"Link Margin: {result['link_margin_dB']:.2f} dB")
        print(f"Energy per Bit: {result['energy_per_bit']:.2e} J/bit")
        print(f"Spectral Efficiency: {result['spectral_efficiency']:.2e} bits/sec/Hz")
    else:
        print("No suitable downlink RF configuration found.")

    if best_overall_config_uplink:
        device = best_overall_config_uplink['device']
        rx_antenna = best_overall_config_uplink['rx_antenna']
        result = best_overall_config_uplink['result']
        print("\nBest Uplink RF Configuration:")
        print(f"Device Type: {device['type']}")
        print(f"Device Name: {device['name']}")
        print(f"Receive Antenna: {rx_antenna['name']} (Gain: {rx_antenna['gain']} dBi)")
        print(f"Modulation: {best_overall_config_uplink['modulation']}")
        print(f"Maximum Achievable Data Rate: {best_overall_config_uplink['data_rate_bps'] / 1e6:.2f} Mbps")
        print(f"SNR: {result['SNR_dB']:.2f} dB")
        print(f"Link Margin: {result['link_margin_dB']:.2f} dB")
        print(f"Energy per Bit: {result['energy_per_bit']:.2e} J/bit")
        print(f"Spectral Efficiency: {result['spectral_efficiency']:.2e} bits/sec/Hz")
    else:
        print("No suitable uplink RF configuration found.")

if __name__ == "__main__":
    main()