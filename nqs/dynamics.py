from ._C_nqs.dynamics import *
from . import _core


@_core.deprecated("function has been renamed to `timestepper`")
def create_timestepper(*args, **kwargs):
    return timestepper(*args, **kwargs)
