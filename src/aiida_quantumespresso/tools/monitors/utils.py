import re
import numpy as np


def parse_accuracy_lines(content):
    """Parse lines containing the word "accuracy" and extract their numerical values.

    Searches for lines containing the word "accuracy" (case-insensitive), then
    extracts the last floating-point number found on each matching line.
    Scientific notation and sign prefixes are supported.

    Args:
        content (str): Multiline string to parse.

    Returns:
        np.ndarray or None: Array of extracted float values, or None if none are found.

    Example:
        >>> parse_accuracy_lines("Final accuracy: 98.5 %\\nTest accuracy: 0.965 units")
        array([98.5 , 0.965])
    """
    word_pat = re.compile(r"\baccuracy\b", flags=re.IGNORECASE)

    # Match a number followed by a non-space token (the unit)
    num_unit_pat = re.compile(
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+(\S+)[\.,;:]*$"
    )

    results = []
    for line in content.splitlines():
        if not word_pat.search(line):
            continue
        m = num_unit_pat.search(line)
        if m:
            num_s = m.group(1)
        else:
            # Fallback: take the last number on the line
            pairs = re.findall(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+(\S+)", line)
            if not pairs:
                continue
            num_s = pairs[-1][0]

        num_s = num_s.replace(",", "")  # allow thousands separators
        try:
            results.append(float(num_s))
        except ValueError:
            continue

    return np.array(results, dtype=float) if results else None


def is_stuck(nums, threshold=5):
    """Determine whether the accuracy metric has plateaued.

    Returns True if the last `threshold` values in `nums` are all equal to the
    final value *and* are consecutive (i.e. the metric has not moved for the
    last `threshold` steps).

    Args:
        nums (array-like): Sequence of numerical values to check.
        threshold (int, optional): Number of consecutive identical values
            required. Defaults to 5.

    Returns:
        bool: True if the metric appears stuck, False otherwise.
    """
    indices = np.where(nums == nums[-1])[0]
    if len(indices) < threshold:
        return False
    return indices[-threshold] == len(nums) - threshold