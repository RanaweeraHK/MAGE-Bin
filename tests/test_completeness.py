import numpy as np

from magebin.completeness import (
    find_terminal_repeats,
    isolate_complete_contigs,
    reverse_complement,
)


def test_reverse_complement():
    assert reverse_complement("AACGTN") == "NACGTT"


def test_direct_terminal_repeat_detection():
    repeat = "ACGTACGTACGTACGTACGTA"
    sequence = repeat + "G" * 80 + repeat
    direct, inverted = find_terminal_repeats(sequence, seed_length=21)
    assert direct == len(repeat)
    assert inverted == 0


def test_complete_contigs_become_singletons():
    labels = np.array([0, 0, 1, 1])
    result = isolate_complete_contigs(
        labels, {"locked": np.array([False, True, False, False])}
    )
    assert result[0] != result[1]
    assert result[2] == result[3]

