"""Validation functions for calculation parameters."""

import warnings

from aiida import orm
from aiida.common.warnings import AiidaDeprecationWarning


def validate_parameters(node: orm.Dict, _) -> str | None:
    """Validate the `parameters` input.

    Checks that the `parameters` input adheres to the following conventions:

    1. Namelist keys should be UPPERCASE, and their values must be a dictionary of
       parameters.
    2. Parameter keys should be lowercase.

    Emits deprecation warnings for incorrect casing. In v5.0, unstored nodes will be
    auto-corrected and stored nodes will be rejected.

    :param node: Dict node containing the namelists/parameters.
    :return: Error message if validation fails, None otherwise.
    """
    for namelist_key, namelist_value in node.items():
        if namelist_key != namelist_key.upper():
            warnings.warn(
                f"Namelist '{namelist_key}' should be UPPERCASE. "
                f"In v5.0, this will be auto-corrected to '{namelist_key.upper()}' for unstored nodes "
                'and enforced for stored nodes. Please update your code to use UPPERCASE namelist names.',
                AiidaDeprecationWarning,
            )

        if not isinstance(namelist_value, dict):
            return f"Namelist '{namelist_key}' should contain a dictionary of parameters, got {type(namelist_value).__name__}"

        for param_key, param_value in namelist_value.items():
            if param_key != param_key.lower():
                warnings.warn(
                    f"Parameter '{param_key}' in namelist '{namelist_key}' should be lowercase. "
                    f"In v5.0, this will be auto-corrected to '{param_key.lower()}' for unstored nodes "
                    'and enforced for stored nodes. Please update your code to use lowercase parameter names.',
                    AiidaDeprecationWarning,
                )
