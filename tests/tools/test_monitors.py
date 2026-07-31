"""Tests for the monitors tools."""

import numpy as np

from aiida_quantumespresso.tools.monitors.utils import is_stuck, parse_accuracy_lines


# ---------------------------------------------------------------------------
# parse_accuracy_lines
# ---------------------------------------------------------------------------


def test_parse_accuracy_lines_basic():
    """Basic multi-line input returns one value per matching line."""
    content = 'Final accuracy: 98.5 %\nTest accuracy: 0.965 units'
    result = parse_accuracy_lines(content)
    np.testing.assert_array_almost_equal(result, [98.5, 0.965])


def test_parse_accuracy_lines_case_insensitive():
    """The word 'accuracy' is matched case-insensitively."""
    content = 'ACCURACY = 1.23 eV'
    result = parse_accuracy_lines(content)
    np.testing.assert_array_almost_equal(result, [1.23])


def test_parse_accuracy_lines_scientific_notation():
    """Scientific-notation values are parsed correctly."""
    content = 'accuracy: 1.5e-3 Ry'
    result = parse_accuracy_lines(content)
    np.testing.assert_array_almost_equal(result, [1.5e-3])


def test_parse_accuracy_lines_no_match_returns_none():
    """Returns None when no line contains the word 'accuracy'."""
    content = 'convergence achieved after 10 iterations'
    assert parse_accuracy_lines(content) is None


def test_parse_accuracy_lines_partial_word_not_matched():
    """Substrings like 'inaccuracy' must NOT be matched (whole-word boundary)."""
    content = 'inaccuracy measure: 0.5 %'
    assert parse_accuracy_lines(content) is None


def test_parse_accuracy_lines_line_without_number_skipped():
    """Lines with 'accuracy' but no parseable number are silently skipped."""
    content = 'accuracy: unknown\naccuracy: 0.1 Ry'
    result = parse_accuracy_lines(content)
    np.testing.assert_array_almost_equal(result, [0.1])


def test_parse_accuracy_lines_thousands_separator():
    """Numbers with comma thousands separators are handled correctly."""
    content = 'accuracy: 1,000 steps'
    result = parse_accuracy_lines(content)
    np.testing.assert_array_almost_equal(result, [1000.0])

    content = 'accuracy: 1000 steps'
    result = parse_accuracy_lines(content)
    np.testing.assert_array_almost_equal(result, [1000.0])


def test_parse_accuracy_lines_empty_string_returns_none():
    """Empty input returns None."""
    assert parse_accuracy_lines('') is None


# ---------------------------------------------------------------------------
# is_stuck
# ---------------------------------------------------------------------------


def test_is_stuck_true_exact_threshold():
    """Returns True when the last `min_repeats` values are identical and consecutive."""
    nums = np.array([1.0, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5])
    assert is_stuck(nums, min_repeats=5) is True


def test_is_stuck_false_fewer_than_threshold():
    """Returns False when fewer than `min_repeats` consecutive identical values exist."""
    nums = np.array([1.0, 0.8, 0.5, 0.5, 0.5])
    assert is_stuck(nums, min_repeats=5) is False


def test_is_stuck_false_not_at_end():
    """Returns False when the plateau is not at the tail of the sequence."""
    nums = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.2])
    assert is_stuck(nums, min_repeats=5) is False


def test_is_stuck_custom_threshold():
    """Custom min_repeats is respected."""
    nums = np.array([1.0, 0.9, 0.9, 0.9])
    assert is_stuck(nums, min_repeats=3) is True
    assert is_stuck(nums, min_repeats=4) is False


def test_is_stuck_single_element():
    """A single-element array with min_repeats=1 is considered stuck."""
    nums = np.array([0.42])
    assert is_stuck(nums, min_repeats=1) is True


def test_is_stuck_all_different():
    """Strictly decreasing sequence is never stuck."""
    nums = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
    assert is_stuck(nums, min_repeats=2) is False


# ---------------------------------------------------------------------------
# monitor (integration-level, using mocks)
# ---------------------------------------------------------------------------


class _FakeNode:
    """Minimal CalcJobNode stand-in for testing."""

    def __init__(self, filename):
        self._filename = filename

    @property
    def attributes(self):
        return {'output_filename': self._filename}


class _FakeTransport:
    """Minimal Transport stand-in that writes a fixed string to a local file."""

    def __init__(self, content):
        self._content = content

    def getfile(self, remote_path, local_path):
        with open(local_path, 'w', encoding='utf-8') as fh:
            fh.write(self._content)


def test_monitor_returns_none_when_not_stuck():
    """monitor() returns None when accuracy is still improving."""
    from aiida_quantumespresso.tools.monitors.accuracy_stuck import monitor

    lines = '\n'.join(f'accuracy: {v:.4f} Ry' for v in [0.5, 0.4, 0.3, 0.2, 0.1])
    node = _FakeNode('output.txt')
    transport = _FakeTransport(lines)
    assert monitor(node, transport) is None


def test_monitor_returns_message_when_stuck():
    """monitor() returns a non-None message when accuracy is stuck."""
    from aiida_quantumespresso.tools.monitors.accuracy_stuck import monitor

    lines = '\n'.join('accuracy: 0.1000 Ry' for _ in range(6))
    node = _FakeNode('output.txt')
    transport = _FakeTransport(lines)
    result = monitor(node, transport)
    assert result is not None
    assert 'stuck' in result.lower()


def test_monitor_returns_none_when_no_output_filename():
    """monitor() returns None when the node has no output_filename attribute."""
    from aiida_quantumespresso.tools.monitors.accuracy_stuck import monitor

    node = _FakeNode(None)
    transport = _FakeTransport('')
    assert monitor(node, transport) is None


def test_monitor_respects_custom_threshold():
    """monitor() honors a custom `min_repeats` argument instead of the hardcoded default."""
    from aiida_quantumespresso.tools.monitors.accuracy_stuck import monitor

    lines = '\n'.join('accuracy: 0.1000 Ry' for _ in range(3))
    node = _FakeNode('output.txt')
    transport = _FakeTransport(lines)

    # Only 3 identical values are present, so the default min_repeats (5) is not met.
    assert monitor(node, transport) is None

    # Lowering min_repeats to 3 should now detect the plateau.
    result = monitor(node, transport, min_repeats=3)
    assert result is not None
    assert '3' in result


def test_monitor_returns_none_when_transport_fails():
    """monitor() returns None gracefully when the transport raises an exception."""
    from aiida_quantumespresso.tools.monitors.accuracy_stuck import monitor

    class _FailingTransport:
        def getfile(self, remote_path, local_path):
            raise OSError('connection lost')

    node = _FakeNode('output.txt')
    assert monitor(node, _FailingTransport()) is None
