from io import StringIO
import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from atlas7v import Atlas7V
from compare_readin import _fortran_state, _parse_dump, _python_state


FIXTURES = Path(__file__).with_name("fixtures")
MOLECULES = Path(__file__).resolve().parents[2] / "data" / "atlas_files" / "molecules.dat"


def _fixture(name):
    return FIXTURES / name


def _deck(card, rows, ending="BEGIN"):
    return "\n".join(
        [
            "TEFF 5770 GRAVITY 4.44",
            card,
            *(" ".join(str(value) for value in row) for row in rows),
            ending,
        ]
    ) + "\n"


def test_mode20_continues_after_begin_and_end_between_decks():
    first = [
        (1.0, 5000.0, 2.0, 3.0, 4.0, 5.0, 6.0),
        (2.0, 5100.0, 2.1, 3.1, 4.1, 5.1, 6.1),
        (3.0, 5200.0, 2.2, 3.2, 4.2, 5.2, 6.2),
    ]
    second = [(9.0, 6000.0, 9.1, 9.2, 9.3, 9.4, 9.5)]
    stream = StringIO(
        _deck("READ DECK 3", first)
        + "END\n"
        + _deck("READ DECK 1", second)
    )
    atlas = Atlas7V()

    atlas.readin(stream, mode=20)
    np.testing.assert_array_equal(atlas.RHOX[:3], [1.0, 2.0, 3.0])
    assert atlas.NRHOX == 3

    atlas.readin(stream, mode=20)
    assert atlas.NRHOX == 0

    atlas.readin(stream, mode=20)
    assert atlas.NRHOX == 1
    np.testing.assert_array_equal(atlas.RHOX[:1], [9.0])
    np.testing.assert_array_equal(atlas.T[:1], [6000.0])
    np.testing.assert_array_equal(atlas.RHOX[1:], 0.0)
    np.testing.assert_array_equal(atlas.T[1:], 0.0)


def test_read_deck_and_deck6_populate_common_model_arrays():
    rows = [
        (1.0, 5000.0, 2.0, 3.0, 4.0, 5.0, 6.0),
        (2.0, 5100.0, 2.1, 3.1, 4.1, 5.1, 6.1),
    ]
    deck = Atlas7V().readin(StringIO(_deck("READ DECK 2", rows)), mode=20)
    deck6 = Atlas7V().readin(
        StringIO(
            _deck(
                "READ DECK6 2 RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB",
                rows,
            )
        ),
        mode=20,
    )

    for name in ("RHOX", "T", "P", "XNE", "ABROSS", "VTURB"):
        np.testing.assert_array_equal(
            getattr(deck, name)[:2], getattr(deck6, name)[:2]
        )


def test_mode20_consumes_only_one_deck_from_line_iterator():
    lines = iter(
        (
            _deck(
                "READ DECK 1",
                [(1.0, 5000.0, 2.0, 3.0, 4.0, 5.0, 6.0)],
            )
            + "END\n"
            + _deck(
                "READ DECK 1",
                [(9.0, 6000.0, 9.1, 9.2, 9.3, 9.4, 9.5)],
            )
        ).splitlines(keepends=True)
    )

    atlas = Atlas7V()
    atlas.readin(lines, mode=20)
    assert atlas.NRHOX == 1
    atlas.readin(lines, mode=20)
    assert atlas.NRHOX == 0
    atlas.readin(lines, mode=20)
    assert atlas.NRHOX == 1
    np.testing.assert_array_equal(atlas.RHOX[:1], [9.0])


def test_map1_preserves_legacy_quadratic_interpolation():
    values = Atlas7V._map1(
        np.array([0.0, 1.0, 2.0, 3.0]),
        np.array([0.0, 1.0, 4.0, 9.0]),
        np.array([0.5, 1.5, 2.5]),
    )
    np.testing.assert_allclose(values, [0.5, 2.5, 6.25])


def test_scale_model_resamples_and_reconstructs_pressure():
    atlas = Atlas7V()
    atlas.NRHOX = 4
    atlas.RHOX[:4] = [1.0, 2.0, 3.0, 4.0]
    atlas.T[:4] = [5000.0, 5100.0, 5200.0, 5300.0]
    atlas.P[:4] = [1.0, 2.0, 3.0, 4.0]
    atlas.XNE[:4] = 1.0
    atlas.ABROSS[:4] = 0.1
    atlas.PRAD[:4] = [0.0, 0.1, 0.2, 0.3]
    atlas.VTURB[:4] = 1.0
    atlas.BMIN[:4] = 1.0
    atlas.BHYD[:4, :6] = 1.0
    atlas.TEFF = 5000.0
    atlas.GRAV = 1.0e4

    atlas.readin(StringIO("SCALE MODEL 4 -1 .5 5100 4\nBEGIN\n"), mode=2)

    np.testing.assert_allclose(atlas.TAUROS[:4], 10 ** np.array([-1.0, -0.5, 0.0, 0.5]))
    assert atlas.TEFF == 5100.0
    assert atlas.GRAV == 1.0e4
    assert np.all(np.isfinite(atlas.P[:4]))
    np.testing.assert_allclose(
        atlas.RHOX[:4], (atlas.PTOTAL[:4] - atlas.PZERO) / atlas.GRAV
    )


def test_mode_1_calculate_fixture():
    atlas = Atlas7V().readin(_fixture("mode1.atlas"), mode=1)

    assert atlas.TEFF == 5770.0
    np.testing.assert_allclose(atlas.GRAV, 10**4.44)
    assert atlas.NUMITS == 1
    assert atlas.NUMNU == 1
    assert atlas.NRHOX == 3
    np.testing.assert_allclose(atlas.TAUROS[:3], 10 ** np.array([-2.0, -1.5, -1.0]))


def test_mode_2_deck6_fixture():
    atlas = Atlas7V().readin(_fixture("mode2.atlas"), mode=2)

    assert atlas.NRHOX == 4
    np.testing.assert_allclose(atlas.RHOX[:4], [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(atlas.T[:4], [5000.0, 5100.0, 5200.0, 5300.0])
    np.testing.assert_allclose(atlas.PRAD[:4], [5.0, 10.05, 15.2, 20.45])


def test_parse_fortran_reference_scalars(tmp_path):
    dump = tmp_path / "state.txt"
    dump.write_text("S TEFF   5.7700000000000000E+03\nI NRHOX           3\n")

    assert _parse_dump(dump) == {"TEFF": 5770.0, "NRHOX": 3}


def test_mode_20_end():
    state = Atlas7V()
    with _fixture("mode20.atlas").open() as deck:
        state.readin(deck, mode=20)
        assert state.NRHOX == 2
        state.readin(deck, mode=20)
    assert state.NRHOX == 0


@pytest.mark.parametrize(
    ("fixture", "attribute", "expected"),
    [
        ("abundance_table.atlas", "ABUND", 0.922),
        ("frequencies.atlas", "NUMNU", 2),
        ("starting_t_tau.atlas", "NRHOX", 2),
        ("departure_coefficients.atlas", "NLTEON", 1),
    ],
)
def test_read_card_fixtures(fixture, attribute, expected):
    atlas = Atlas7V().readin(_fixture(fixture), mode=2)

    actual = getattr(atlas, attribute)
    if isinstance(actual, np.ndarray):
        actual = actual[0]
    assert actual == pytest.approx(expected)


def test_change_rhox_fixture():
    atlas = Atlas7V().readin(_fixture("change_rhox.atlas"), mode=2)

    assert atlas.NRHOX == 2
    np.testing.assert_allclose(atlas.RHOX[:2], [1.0, 3.0])


def test_scale_model_fixture():
    atlas = Atlas7V().readin(_fixture("scale_model.atlas"), mode=2)

    assert atlas.NRHOX == 4
    assert atlas.TEFF == 5100.0
    assert np.all(np.isfinite(atlas.P[:4]))


@pytest.mark.skipif(
    not (os.environ.get("FC") or shutil.which("gfortran") or shutil.which("ifx") or shutil.which("ifort")),
    reason="Fortran compiler is required for the reference-state comparison",
)
def test_mode_2_matches_fortran_reference_state():
    expected = _fortran_state(_fixture("mode2.atlas"), MOLECULES, os.environ.get("FC"))
    actual = _python_state(_fixture("mode2.atlas"))

    assert actual["NRHOX"] == expected["NRHOX"]
    for name in ("RHOX", "T", "P", "XNE", "ABROSS", "PRAD", "VTURB"):
        np.testing.assert_allclose(
            actual[name][:expected["NRHOX"]],
            expected[name][:expected["NRHOX"]],
            rtol=1e-12,
            atol=1e-12,
        )
