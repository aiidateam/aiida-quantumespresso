import os
import tempfile
from typing import Optional

from aiida.orm import CalcJobNode
from aiida.transports import Transport

from .utils import parse_accuracy_lines, is_stuck


def _fetch_output_content(node: CalcJobNode, transport: Transport) -> Optional[str]:
    """Fetch the output file of a running calculation and return its text content.

    Returns None if the output filename cannot be determined or the file cannot
    be retrieved.
    """
    outfile = node.attributes.get("output_filename", None)
    if not outfile:
        return None

    fd, tmpname = tempfile.mkstemp()
    os.close(fd)
    content = None
    try:
        transport.getfile(outfile, tmpname)
        with open(tmpname, "r", encoding="utf-8", errors="replace") as fh:
            content = fh.read()
    except Exception:
        pass
    finally:
        try:
            os.remove(tmpname)
        except OSError:
            pass
    return content


def monitor(node: CalcJobNode, transport: Transport) -> Optional[str]:
    """Monitor a running calculation for a stuck accuracy metric.

    Fetches the output file from the remote working directory, extracts all
    numerical values that appear on lines containing the word "accuracy", and
    checks whether the metric has stopped improving.  Returns a human-readable
    message if the job appears stuck, or None to let the calculation continue.
    """
    content = _fetch_output_content(node, transport)
    if content is None:
        return None

    nums = parse_accuracy_lines(content)
    if nums is not None and is_stuck(nums, threshold=5):
        return "Job appears stuck: accuracy has not improved in the last 5 occurrences."

    return None