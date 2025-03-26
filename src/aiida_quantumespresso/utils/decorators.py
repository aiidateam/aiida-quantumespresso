# -*- coding: utf-8 -*-
"""Decorators for several purposes."""


def remove_none_overrides(func):
    """Remove namespaces of the returned builder of a `get_builder*` method."""

    def recursively_remove_nones(item):
        """Recursively remove keys with None values from dictionaries."""
        if isinstance(item, dict):
            return {key: recursively_remove_nones(value) for key, value in item.items() if value is not None}
        return item

    def remove_keys_from_builder(builder, keys, path=()):
        """Recursively remove specified keys from the builder based on a path."""
        if not keys:
            return
        current_level = keys.pop(0)
        if hasattr(builder, current_level):
            if keys:
                next_attr = getattr(builder, current_level)
                remove_keys_from_builder(next_attr, keys, path + (current_level,))
            else:
                delattr(builder, current_level)

    def wrapper(*args, **kwargs):
        """Wrap the function."""
        if 'overrides' in kwargs and kwargs['overrides'] is not None:
            original_overrides = kwargs['overrides']

            # Identify paths to keys with None values to be removed
            paths_to_remove = []

            def find_paths(item, path=()):
                """Find the paths to remove."""
                if isinstance(item, dict):
                    for key, value in item.items():
                        if value is None:
                            paths_to_remove.append(path + (key,))
                        else:
                            find_paths(value, path + (key,))

            find_paths(original_overrides)

            # Recursively remove keys with None values from overrides
            cleaned_overrides = recursively_remove_nones(original_overrides)
            kwargs['overrides'] = cleaned_overrides

            # Call the original function to get the builder
            builder = func(*args, **kwargs)

            # Remove specified keys from the builder
            for path in paths_to_remove:
                remove_keys_from_builder(builder, list(path))

            return builder

        return func(*args, **kwargs)

    return wrapper
