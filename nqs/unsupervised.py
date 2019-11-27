from ._C_nqs.unsupervised import *

import itertools


def _Qsr_iter(self, n_iter=None, step_size=1):
    """
    Returns a generator which advances the QSR optimization, yielding
    after every step_size steps up to n_iter.

    Args:
        n_iter (int, optional): The number of steps or None, for no limit.
        step_size (int): The number of steps the simulation is advanced.

    Yields:
        int: The current step.
    """
    self.reset()
    for i in itertools.count(step=step_size):
        if n_iter and i >= n_iter:
            return
        self.advance(step_size)
        yield i


Qsr.iter = _Qsr_iter
