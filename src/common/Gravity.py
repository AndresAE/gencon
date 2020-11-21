from src.common import constants


class Gravity:
    def __init__(self, altitude):
        self.altitude = altitude

    def gravity(self):
        """return acceleration due to gravity density wrt altitude."""
        r = constants.radius_earth() + self.altitude / constants.m2ft()
        g = constants.m2ft() * constants.gravitational_constant() * constants.mass_earth() / r ** 2
        return g

# Public Methods #######################################################################################################
