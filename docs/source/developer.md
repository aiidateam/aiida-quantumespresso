# Developer Guide

## Documentation

We use the [Diátaxis](https://diataxis.fr/) approach for organising the documentation in four sections:

* Tutorials
* How-to's
* Topics (Explanation)
* Reference

All of our documentation should be written in [MyST Markdown](https://myst-parser.readthedocs.io/en/latest).
Below you can find a list of current style guide items:

1. Write **one sentence per line** and otherwise **no manual line wrapping** to make easy to create and review  diffs.
   All standard editors allow for dynamic line wrapping, and the line length is irrelevant for the rendered documentation in, e.g., HTML or PDF format.
1. **File and directory names should be alphanumeric** and all lower-case with underscores as word-separators. Example: `entry_points.rst`
1. **Headers must be set in sentence-case**.
   Example: "Entry points"
1. Separate paragraphs by one empty line, but not more.
1. Use the `-` symbol for itemized lists.

### Notes

#### `run` vs `submit`

Originally discussed in [this issue](https://github.com/aiidateam/aiida-quantumespresso/issues/1127).

In the majority of real use cases, `submit` is the preferred engine command, and the documentation should reflect that.
Hence, we limit using the `run` commando to the quick start documentation, and perhaps showing other use cases such as testing calculation setups or workflows (with caching).
