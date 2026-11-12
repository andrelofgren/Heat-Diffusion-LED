import sys
import argparse
import json
import numpy as np
import matplotlib as mpl
# Import some additional colormaps
try:
    import cmasher
except ImportError as ie:
    print(f"Package {ie.name} not found, unable to import some colormaps")
try:
    import cmocean
except ImportError as ie:
    print(f"Package {ie.name} not found, unable to import some colormaps")
try:
    import sciviz
except ImportError as ie:
    print(f"Package {ie.name} not found, unable to import some colormaps")


parser = argparse.ArgumentParser(
    prog="set_color_from_cmap",
    description="Generates rgb values from named colormaps"
)

parser.add_argument("-c", "--config", type=str, help="path of json config file")
parser.add_argument("-m", "--cmap", type=str, help="name of colormap")
parser.add_argument("-n", "--nbins", type=int, help="number of rgb samples of the colormap")
parser.add_argument("-l", "--list-cmaps", action="store_true", help="list all available colormaps")

args = parser.parse_args()

if args.list_cmaps:
    print("Available colormaps:")
    for cmap in mpl.colormaps:
        print(cmap)
    sys.exit(0)

if args.config is None:
    print("Missing required argument '--config'")
    sys.exit(-1)

if args.nbins is None:
    print("Missing required argument '--nbins'")
    sys.exit(-1)

if args.cmap is None:
    print("Missing required argument '--cmap'")
    sys.exit(-1)

num_bins = args.nbins
cmap_name = args.cmap
config_path = args.config

try:
    cmap = mpl.colormaps[cmap_name]
except KeyError:
    print(f"Colormap {cmap_name} not found")
    sys.exit(-1)

cmap_rgb_vals = cmap(np.linspace(0, 1, num_bins))[:, :-1]

with open(config_path, "r") as json_file:
    params = json.load(json_file)

params["colormap"]["num_bins"] = num_bins
colormap = {
    "num_bins": num_bins
}
for i, color in enumerate(cmap_rgb_vals):
    color_key = f"color_{i+1}"
    colormap.update({color_key: list(color)})
params.update({"colormap": colormap})

with open(config_path, "w") as json_file:
    json.dump(params, json_file, indent=4)
