import math
import unittest

import highschool_physics as hp


class TestHighSchoolPhysics(unittest.TestCase):
    def test_kinetic_energy(self):
        self.assertAlmostEqual(hp.kinetic_energy({"m": 2.0, "v": 3.0}), 9.0)

    def test_gravitational_potential_energy_default_g(self):
        # Using default g0=9.81
        self.assertAlmostEqual(hp.gravitational_pe({"m": 2.0, "h": 10.0}), 2.0 * hp.g0 * 10.0)

    def test_gravitational_potential_energy_custom_g(self):
        self.assertAlmostEqual(hp.gravitational_pe({"m": 1.0, "h": 1.0, "g": 10.0}), 10.0)

    def test_final_velocity(self):
        self.assertAlmostEqual(hp.final_velocity({"u": 5.0, "a": 2.0, "t": 3.0}), 11.0)

    def test_displacement_suvat(self):
        # s = u t + 0.5 a t^2 = 5*2 + 0.5*3*4 = 10 + 6 = 16
        self.assertAlmostEqual(hp.displacement_suvat({"u": 5.0, "a": 3.0, "t": 2.0}), 16.0)

    def test_final_velocity_suvat(self):
        # v^2 = u^2 + 2 a s => v = sqrt(4 + 2*3*2) = sqrt(16) = 4
        self.assertAlmostEqual(hp.final_velocity_suvat({"u": 2.0, "a": 3.0, "s": 2.0}), 4.0)
        self.assertTrue(math.isnan(hp.final_velocity_suvat({"u": 0.0, "a": -1.0, "s": 1.0})))

    def test_centripetal(self):
        self.assertAlmostEqual(hp.centripetal_acc({"v": 4.0, "r": 2.0}), 8.0)
        self.assertAlmostEqual(hp.centripetal_force({"m": 1.5, "v": 4.0, "r": 2.0}), 12.0)

    def test_ohms_law_and_power(self):
        self.assertAlmostEqual(hp.ohms_law({"I": 2.0, "R": 5.0}), 10.0)
        self.assertAlmostEqual(hp.electrical_power({"V": 10.0, "I": 2.0}), 20.0)

    def test_pressure_and_density(self):
        self.assertAlmostEqual(hp.pressure({"F": 10.0, "A": 2.0}), 5.0)
        self.assertAlmostEqual(hp.density({"m": 10.0, "V": 2.0}), 5.0)

    def test_hooke_and_periods(self):
        self.assertAlmostEqual(hp.hooke_force({"k": 100.0, "x": 0.1}), 10.0)
        self.assertAlmostEqual(hp.period_mass_spring({"m": 1.0, "k": 4.0}), 2 * math.pi * math.sqrt(0.25))
        self.assertAlmostEqual(hp.period_pendulum({"L": 1.0, "g": 9.81}), 2 * math.pi * math.sqrt(1.0 / 9.81))

    def test_wave_and_photon(self):
        self.assertAlmostEqual(hp.wave_speed({"f": 2.0, "lam": 3.0}), 6.0)
        self.assertAlmostEqual(hp.photon_energy({"f": 1.0e14}), hp.h * 1.0e14)

    def test_coulomb_force(self):
        # F = k q1 q2 / r^2
        self.assertAlmostEqual(hp.coulomb_force({"q1": 1e-6, "q2": 2e-6, "r": 0.1}), hp.k_e * 2e-12 / 0.01)


if __name__ == "__main__":
    unittest.main(verbosity=2)
