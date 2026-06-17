import datetime
from django import template

register = template.Library()

@register.filter
def subtract(value, arg):
    try:
        return int(value) - int(arg)
    except (ValueError, TypeError):
        return ''

@register.filter
def print_timestamp(ts_in):
    #//TODO - ensure in local time
    try:
        #assume, that timestamp is given in seconds with decimal point
        ts = float(ts_in)
    except ValueError:
        return None
    return datetime.datetime.fromtimestamp(ts)


@register.filter
def mean_with(value, other):
    """Return the arithmetic mean of two values, or None if either is invalid."""
    try:
        left = float(value)
        right = float(other)
    except (TypeError, ValueError):
        return None
    return 0.5 * (left + right)