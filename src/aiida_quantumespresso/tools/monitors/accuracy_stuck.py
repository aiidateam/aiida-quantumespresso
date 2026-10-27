import os
import tempfile
from typing import Optional

from aiida.orm import CalcJobNode
from aiida.transports import Transport

from .utils import parse_accuracy_lines, is_stuck

# Message prefix identifying this monitor as the cause of a `STOPPED_BY_MONITOR` exit code, since that exit code is
# shared by all monitors. Error handlers (e.g. in `PwBaseWorkChain`) match on this prefix, as the full message
# depends on the configurable `min_repeats`.
ACCURACY_STUCK_MESSAGE_PREFIX = 'Job appears stuck: accuracy has not improved'


def _fetch_output_content(node: CalcJobNode, transport: Transport) -> Optional[str]:
    """Fetch the output file of a running calculation and return its text content.

    Returns None if the output filename cannot be determined or the file cannot
    be retrieved.
    """
    outfile = node.attributes.get('output_filename', None)
    if not outfile:
        return None

    fd, tmpname = tempfile.mkstemp()
    os.close(fd)
    content = None
    try:
        transport.getfile(outfile, tmpname)
        with open(tmpname, 'r', encoding='utf-8', errors='replace') as fh:
            content = fh.read()
    except Exception:
        pass
    finally:
        try:
            os.remove(tmpname)
        except OSError:
            pass
    return content


def monitor(node: CalcJobNode, transport: Transport, min_repeats: int = 5) -> Optional[str]:
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
    if nums is not None and is_stuck(nums, min_repeats=min_repeats):
        return f'{ACCURACY_STUCK_MESSAGE_PREFIX} in the last {min_repeats} occurrences.'

    return None
