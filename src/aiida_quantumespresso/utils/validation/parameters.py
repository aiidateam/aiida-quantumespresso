"""Validation functions for calculation parameters."""

from aiida import orm
import warnings


def validate_parameters(node: orm.Dict, _) -> str | None:
    """Validate the `parameters` input.

    Makes sure the `parameters` input adheres to the following conventions:

    1. Namelist keys should be UPPERCASE, and their values must be a dictionary of
       parameters.
    2. Parameter keys should be lowercase.

    Auto-corrects unstored nodes and emits warnings, rejects stored nodes with error
    message.

    :param node: Dict node containing the namelists/parameters.
    """
    normalized = {}

    for namelist_key, namelist_value in node.items():
        if namelist_key != namelist_key.upper():
            if node.is_stored:
                return (
                    f"Namelist '{namelist_key}' should be UPPERCASE. "
                    'Since the Dict node is stored, it cannot be auto-corrected. '
                    f"Please use '{namelist_key.upper()}' instead."
                )
            else:
                warnings.warn(
                    f"Namelist '{namelist_key}' should be UPPERCASE. Changed to '{namelist_key.upper()}'.", UserWarning
                )
                namelist_key_normalized = namelist_key.upper()
        else:
            namelist_key_normalized = namelist_key

        if not isinstance(namelist_value, dict):
            return f"Namelist '{namelist_key_normalized}' should contain a dictionary of parameters, got {type(namelist_value).__name__}"

        normalized_params = {}
        for param_key, param_value in namelist_value.items():
            if param_key != param_key.lower():
                if node.is_stored:
                    return (
                        f"Parameter '{param_key}' in namelist '{namelist_key_normalized}' should be lowercase. "
                        'Since the Dict node is stored, it cannot be auto-corrected. '
                        f"Please use '{param_key.lower()}' instead."
                    )
                else:
                    warnings.warn(
                        f"Parameter '{param_key}' in namelist '{namelist_key_normalized}' should be lowercase. "
                        f"Changed to '{namelist_key_normalized}.{param_key.lower()}'.",
                        UserWarning,
                    )
                    param_key_normalized = param_key.lower()
            else:
                param_key_normalized = param_key
            normalized_params[param_key_normalized] = param_value
        normalized[namelist_key_normalized] = normalized_params

    if normalized != node.get_dict():
        node.set_dict(normalized)
