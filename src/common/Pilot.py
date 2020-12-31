from control import series, tf
from src.common.tools import lead_lag, pade_1


class Pilot:
    def __init__(self, gain, delay):
        self.gain = gain
        self.delay = delay

    def delay(self):
        """return acceleration due to gravity density wrt altitude."""
        pilot = series(self.gain, pade_1(self.delay))
        return pilot

    def delay_lead_lag(self, lead, lag):
        """return acceleration due to gravity density wrt altitude."""
        pilot = series(self.gain, series(pade_1(self.delay), lead_lag(lead, lag)))
        return pilot

    def delay_lead_lag_nmd(self, lead, lag, nm_delay):
        """return acceleration due to gravity density wrt altitude."""
        pilot = series(self.gain, series(pade_1(self.delay), series(lead_lag(lead, lag), tf(1, [nm_delay, 1]))))
        return pilot

# Public Methods #######################################################################################################
