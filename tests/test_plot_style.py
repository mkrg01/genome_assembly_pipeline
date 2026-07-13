import sys
from pathlib import Path
from types import SimpleNamespace


SCRIPTS_DIR = Path(__file__).parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from plot_style import (  # noqa: E402
    BASE_FONT_SIZE,
    FONT_FAMILY,
    FONT_FALLBACKS,
    SMALL_FONT_SIZE,
    TITLE_FONT_SIZE,
    apply_matplotlib_style,
    display_organism_name,
)


def test_publication_typography_constants():
    assert FONT_FAMILY == "Helvetica"
    assert FONT_FALLBACKS == (
        "Helvetica",
        "Arial",
        "Nimbus Sans",
        "Liberation Sans",
        "DejaVu Sans",
    )
    assert BASE_FONT_SIZE == 8
    assert SMALL_FONT_SIZE == 6
    assert TITLE_FONT_SIZE == 8


def test_display_organism_name_replaces_filename_separators():
    assert display_organism_name("Dioncophyllum_thollonii") == (
        "Dioncophyllum thollonii"
    )


def test_apply_matplotlib_style_sets_sizes_and_editable_vector_text(monkeypatch):
    font_manager = SimpleNamespace(
        FontProperties=lambda family: SimpleNamespace(family=family),
        findfont=lambda properties, fallback_to_default: f"/{properties.family}.ttf",
    )
    matplotlib = SimpleNamespace(rcParams={}, font_manager=font_manager)
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib)

    apply_matplotlib_style()

    assert matplotlib.rcParams["font.sans-serif"][0] == "Helvetica"
    assert matplotlib.rcParams["font.size"] == 8
    assert matplotlib.rcParams["axes.titlesize"] == 8
    assert matplotlib.rcParams["axes.titleweight"] == "bold"
    assert matplotlib.rcParams["xtick.labelsize"] == 6
    assert matplotlib.rcParams["ytick.labelsize"] == 6
    assert matplotlib.rcParams["mathtext.it"] == "Helvetica:italic"
    assert matplotlib.rcParams["pdf.fonttype"] == 42
    assert matplotlib.rcParams["svg.fonttype"] == "none"
