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


@register.filter
def strip_trailing_zeros(value):
    """
    Strip trailing zeros (and a trailing decimal point) from a formatted
    numeric string, e.g. "0.9800" -> "0.98", "34.0000" -> "34".

    Intended to be chained after `floatformat`, which always pads to a
    fixed number of decimal places, e.g. `{{ value|floatformat:4|strip_trailing_zeros }}`.
    """
    text = str(value)
    if '.' not in text:
        return text
    return text.rstrip('0').rstrip('.')


@register.filter
def format_float(value, decimals=4):
    """
    Format a numeric value with at most `decimals` decimal places, stripping
    any trailing zeros (and a trailing decimal point), e.g.
    "300.0" -> "300", "0.98000" -> "0.98", "0.12345" -> "0.1235" (4 dp).

    Unlike Django's builtin `floatformat`, this always drops insignificant
    trailing zeros regardless of sign of the `decimals` argument.
    """
    try:
        number = float(value)
    except (TypeError, ValueError):
        return value
    text = f"{number:.{int(decimals)}f}"
    if '.' in text:
        text = text.rstrip('0').rstrip('.')
    return text


@register.filter
def format_int(value):
    """
    Coerce a numeric value to an integer string, e.g. "300.0" -> "300",
    300.0 -> "300". Used for int-keytype fields whose stored value/default
    may be an unformatted float (no decimal point should ever be shown).
    """
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return value
