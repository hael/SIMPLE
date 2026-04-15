import datetime
from django import template

register = template.Library()

@register.filter
def print_timestamp(ts_in):
    #//TODO - ensure in local time
    try:
        #assume, that timestamp is given in seconds with decimal point
        ts = float(ts_in)
    except ValueError:
        return None
    return datetime.datetime.fromtimestamp(ts)