import matplotlib

# My colors
nice_colors = ['peru','maroon','darkgreen','indigo', 'plum'] # autumn at IAR

# My plot config
nice_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": 'serif',
    # Use 10pt font in plots, to match 10pt font in the document
    "axes.labelsize": 18,
    "font.size": 20,
    "axes.linewidth": 1,
    "axes.titlesize": 16,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 14,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "xtick.major.size": 5,  # major tick size in points
    "xtick.minor.size": 5,  # minor tick size in points
    "xtick.major.width": 1.4,  # major tick width in points
    "xtick.minor.width": 1.4,  # minor tick width in points
}

matplotlib.rcParams.update(nice_fonts)

#    if ax is not None:
#       ax.tick_params(axis="y", direction="in", length=10, width=1, color="black")
#       ax.tick_params(which='minor', axis="y", direction="in", length=5, width=1, color="black")
#       ax.tick_params(axis="x", direction="in", length=10, width=1, color="black")
#       ax.tick_params(which='minor', axis="x", direction="in", length=5, width=1, color="black")
#       ax.set_axisbelow(False)
