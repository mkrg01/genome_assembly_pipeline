"""Shared typography settings for plots produced by the workflow."""

FONT_FAMILY = "Helvetica"
# Helvetica and Arial are proprietary and are not installed on every Linux
# host. The remaining entries keep plots usable while preserving sans-serif
# typography; installing Helvetica makes it the automatically selected face.
FONT_FALLBACKS = (
    FONT_FAMILY,
    "Arial",
    "Nimbus Sans",
    "Liberation Sans",
    "DejaVu Sans",
)
BASE_FONT_SIZE = 8
SMALL_FONT_SIZE = 6
TITLE_FONT_SIZE = 8


def display_organism_name(organism_name):
    """Convert a filesystem-safe organism label to a figure label."""
    return " ".join(str(organism_name).replace("_", " ").split())


def multiline_organism_name(organism_name):
    """Format an organism label on two lines, breaking after the genus."""
    genus_and_rest = display_organism_name(organism_name).split(maxsplit=1)
    return "\n".join(genus_and_rest)


def _first_available_font(font_manager):
    for family in FONT_FALLBACKS:
        try:
            font_manager.findfont(
                font_manager.FontProperties(family=family),
                fallback_to_default=False,
            )
        except ValueError:
            continue
        return family
    return FONT_FALLBACKS[-1]


def apply_matplotlib_style():
    """Apply the workflow-wide publication typography to Matplotlib."""
    import matplotlib
    from matplotlib import font_manager

    math_font = _first_available_font(font_manager)

    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": list(FONT_FALLBACKS),
            "font.size": BASE_FONT_SIZE,
            "axes.titlesize": TITLE_FONT_SIZE,
            "axes.titleweight": "bold",
            "axes.labelsize": BASE_FONT_SIZE,
            "xtick.labelsize": SMALL_FONT_SIZE,
            "ytick.labelsize": SMALL_FONT_SIZE,
            "legend.fontsize": BASE_FONT_SIZE,
            "figure.titlesize": TITLE_FONT_SIZE,
            "figure.titleweight": "bold",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.grid": False,
            "legend.frameon": False,
            # Match MathText (used for italic scientific names) to body text.
            "mathtext.fontset": "custom",
            "mathtext.rm": math_font,
            "mathtext.it": f"{math_font}:italic",
            "mathtext.bf": f"{math_font}:bold",
            "mathtext.bfit": f"{math_font}:bold:italic",
            "mathtext.cal": math_font,
            "mathtext.tt": math_font,
            "mathtext.sf": math_font,
            # Keep text editable in vector output and embed TrueType fonts.
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
