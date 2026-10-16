"""Selectable analytic adapters; molecular adapters share the same runner."""
from .motion import MotionAdapter
from .ctf import CtfAdapter
from .picking import PickAdapter
from .solve2d import Solve2DAdapter
from .solve3d import Solve3DAdapter

ADAPTERS = {a.name: a for a in (MotionAdapter(), CtfAdapter(), PickAdapter(), Solve2DAdapter(), Solve3DAdapter())}
