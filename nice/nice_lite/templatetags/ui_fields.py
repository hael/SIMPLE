from django import template

register = template.Library()


@register.filter
def bracket_default(value):
    """
    Extract the actual default value from a SIMPLE UI descr_placeholder string.

    `binary`/`multi` keytypes encode their allowed options and default as a
    single string, e.g. "(yes|no){no}" or "(none|flip_auto|...){none}", so the
    raw string is not directly usable as a field value or comparison target.
    Scalar keytypes (int/float/dir/file) store a plain default with no
    braces, e.g. "2.7", and are returned unchanged.
    """
    if not isinstance(value, str):
        return value
    start = value.rfind('{')
    end = value.rfind('}')
    if start != -1 and end != -1 and end > start:
        return value[start + 1:end]
    return value
