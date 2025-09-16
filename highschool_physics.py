import math
import sys

#!/usr/bin/env python3
"""
Simple Physics Calculator (high-school level)

Run this script and pick an equation by number. You'll be prompted for the required
variables and the script will compute the result.
"""


# Physical constants
G = 6.67430e-11          # gravitational constant, m^3 kg^-1 s^-2
g0 = 9.81                # standard gravity, m/s^2
c = 299_792_458          # speed of light m/s
h = 6.62607015e-34       # Planck's constant J*s
k_e = 8.9875517923e9     # Coulomb constant N*m^2/C^2
PI = math.pi

# Helper input functions
def ask_float(prompt):
    """Ask for a float repeatedly until valid input is given or user enters 'q' to quit."""
    while True:
        s = input(prompt).strip()
        if s.lower() in ('q', 'quit', 'exit'):
            print("Exiting.")
            sys.exit(0)
        try:
            return float(s)
        except ValueError:
            print("Please enter a number (or 'q' to quit).")

def show_result(name, value, unit=""):
    """Print a neat formatted result."""
    print(f"\n{name} = {value:.6g} {unit}\n")


# Equation implementations
def kinetic_energy(vars):
    # KE = 1/2 m v^2
    m = vars['m']; v = vars['v']
    return 0.5 * m * v * v

def gravitational_pe(vars):
    # U = m g h
    m = vars['m']; h_ = vars['h']; g = vars.get('g', g0)
    return m * g * h_

def momentum(vars):
    # p = m v
    return vars['m'] * vars['v']

def newton_f(vars):
    # F = m a
    return vars['m'] * vars['a']

def work_force_distance(vars):
    # W = F d cos(theta)
    F = vars['F']; d = vars['d']; theta_deg = vars.get('theta', 0.0)
    theta_rad = math.radians(theta_deg)
    return F * d * math.cos(theta_rad)

def power_from_work(vars):
    # P = W / t
    return vars['W'] / vars['t']

def final_velocity(vars):
    # v = u + a t
    return vars['u'] + vars['a'] * vars['t']

def displacement_suvat(vars):
    # s = u t + 0.5 a t^2
    return vars['u'] * vars['t'] + 0.5 * vars['a'] * vars['t']**2

def final_velocity_suvat(vars):
    # v^2 = u^2 + 2 a s  -> return v (positive root if possible; may return NaN if negative)
    u = vars['u']; a = vars['a']; s = vars['s']
    val = u*u + 2*a*s
    if val < 0:
        return float('nan')
    return math.sqrt(val)

def centripetal_acc(vars):
    # a_c = v^2 / r
    v = vars['v']; r = vars['r']
    return (v*v) / r

def centripetal_force(vars):
    # F = m v^2 / r
    return vars['m'] * vars['v']**2 / vars['r']

def pressure(vars):
    # P = F / A
    return vars['F'] / vars['A']

def density(vars):
    # rho = m / V
    return vars['m'] / vars['V']

def hooke_force(vars):
    # F = k x (magnitude)
    return vars['k'] * vars['x']

def period_mass_spring(vars):
    # T = 2*pi*sqrt(m/k)
    return 2 * PI * math.sqrt(vars['m'] / vars['k'])

def period_pendulum(vars):
    # T = 2*pi*sqrt(L/g)
    return 2 * PI * math.sqrt(vars['L'] / vars.get('g', g0))

def wave_speed(vars):
    # v = f * lambda
    return vars['f'] * vars['lam']

def photon_energy(vars):
    # E = h f
    return h * vars['f']

def ohms_law(vars):
    # V = I R
    return vars['I'] * vars['R']

def electrical_power(vars):
    # P = V I
    return vars['V'] * vars['I']

def coulomb_force(vars):
    # F = k q1 q2 / r^2
    q1 = vars['q1']; q2 = vars['q2']; r = vars['r']
    return k_e * q1 * q2 / (r*r)

# Catalog of equations: id -> {name, vars, func, unit, formula}
EQUATIONS = {
    1:  {"name": "Kinetic energy", "vars": [("m","mass (kg)"), ("v","speed (m/s)")], "func": kinetic_energy, "unit":"J", "formula":"KE = 1/2 m v^2"},
    2:  {"name": "Gravitational potential energy (near Earth)", "vars":[("m","mass (kg)"), ("h","height (m)"), ("g","gravity (m/s^2) [press enter for 9.81]")], "func": gravitational_pe, "unit":"J", "formula":"U = m g h"},
    3:  {"name": "Linear momentum", "vars":[("m","mass (kg)"), ("v","velocity (m/s)")], "func": momentum, "unit":"kg·m/s", "formula":"p = m v"},
    4:  {"name": "Newton's 2nd law (Force)", "vars":[("m","mass (kg)"), ("a","acceleration (m/s^2)")], "func": newton_f, "unit":"N", "formula":"F = m a"},
    5:  {"name": "Work by a constant force", "vars":[("F","force (N)"), ("d","displacement (m)"), ("theta","angle between F and d (deg)")], "func": work_force_distance, "unit":"J", "formula":"W = F d cos(theta)"},
    6:  {"name": "Power from work", "vars":[("W","work (J)"), ("t","time (s)")], "func": power_from_work, "unit":"W", "formula":"P = W / t"},
    7:  {"name": "Final velocity (v = u + a t)", "vars":[("u","initial velocity (m/s)"), ("a","acceleration (m/s^2)"), ("t","time (s)")], "func": final_velocity, "unit":"m/s", "formula":"v = u + a t"},
    8:  {"name": "Displacement (s = u t + 1/2 a t^2)", "vars":[("u","initial velocity (m/s)"), ("a","acceleration (m/s^2)"), ("t","time (s)")], "func": displacement_suvat, "unit":"m", "formula":"s = u t + 1/2 a t^2"},
    9:  {"name": "Final speed from v^2 = u^2 + 2 a s", "vars":[("u","initial speed (m/s)"), ("a","acceleration (m/s^2)"), ("s","displacement (m)")], "func": final_velocity_suvat, "unit":"m/s", "formula":"v^2 = u^2 + 2 a s"},
    10: {"name": "Centripetal acceleration", "vars":[("v","speed (m/s)"), ("r","radius (m)")], "func": centripetal_acc, "unit":"m/s^2", "formula":"a_c = v^2 / r"},
    11: {"name": "Centripetal force", "vars":[("m","mass (kg)"), ("v","speed (m/s)"), ("r","radius (m)")], "func": centripetal_force, "unit":"N", "formula":"F = m v^2 / r"},
    12: {"name": "Pressure", "vars":[("F","force (N)"), ("A","area (m^2)")], "func": pressure, "unit":"Pa", "formula":"P = F / A"},
    13: {"name": "Density", "vars":[("m","mass (kg)"), ("V","volume (m^3)")], "func": density, "unit":"kg/m^3", "formula":"rho = m / V"},
    14: {"name": "Hooke's law (spring force)", "vars":[("k","spring constant (N/m)"), ("x","displacement (m)")], "func": hooke_force, "unit":"N", "formula":"F = k x"},
    15: {"name": "Period of mass-spring (T = 2π√(m/k))", "vars":[("m","mass (kg)"), ("k","spring constant (N/m)")], "func": period_mass_spring, "unit":"s", "formula":"T = 2π√(m/k)"},
    16: {"name": "Period of simple pendulum (small angles)", "vars":[("L","length (m)"), ("g","gravity (m/s^2) [press enter for 9.81]")], "func": period_pendulum, "unit":"s", "formula":"T = 2π√(L/g)"},
    17: {"name": "Wave speed (v = f λ)", "vars":[("f","frequency (Hz)"), ("lam","wavelength (m)")], "func": wave_speed, "unit":"m/s", "formula":"v = f λ"},
    18: {"name": "Photon energy (E = h f)", "vars":[("f","frequency (Hz)")], "func": photon_energy, "unit":"J", "formula":"E = h f"},
    19: {"name": "Ohm's law (V = I R)", "vars":[("I","current (A)"), ("R","resistance (Ω)")], "func": ohms_law, "unit":"V", "formula":"V = I R"},
    20: {"name": "Electrical power (P = V I)", "vars":[("V","voltage (V)"), ("I","current (A)")], "func": electrical_power, "unit":"W", "formula":"P = V I"},
    21: {"name": "Coulomb's law (electrostatic force)", "vars":[("q1","charge 1 (C)"), ("q2","charge 2 (C)"), ("r","distance (m)")], "func": coulomb_force, "unit":"N", "formula":"F = k q1 q2 / r^2"},
}


def prompt_for_vars(var_list):
    """Prompt the user for each variable described in var_list.
    var_list: list of (key, prompt_text) tuples
    Returns dict key -> float
    """
    values = {}
    for key, prompt_text in var_list:
        prompt = f"Enter {prompt_text}: "
        # For defaults like gravity, accept empty to use default in function (we store as not provided)
        s = input(prompt).strip()
        if s == "":
            # leave absent so function can use default if coded
            # but if the function requires the value, try to convert (this will error later)
            continue
        if s.lower() in ('q', 'quit', 'exit'):
            print("Exiting.")
            sys.exit(0)
        try:
            values[key] = float(s)
        except ValueError:
            print("Invalid number entered. Please try again.")
            return prompt_for_vars(var_list)
    return values

def main():
    print("Simple Physics Calculator")
    print("Enter 'q' at any prompt to quit.\n")
    while True:
        print("Equations:")
        for idx in sorted(EQUATIONS.keys()):
            print(f"{idx:2d}. {EQUATIONS[idx]['name']}")
        print(" 0. Quit")
        choice = input("\nPick an equation number and press Enter: ").strip()
        if choice.lower() in ('q', 'quit', 'exit', '0', ''):
            print("Goodbye.")
            break
        try:
            ni = int(choice)
        except ValueError:
            print("Please enter a valid number.")
            continue
        if ni not in EQUATIONS:
            print("Unknown selection. Try again.")
            continue

        item = EQUATIONS[ni]
        print(f"\nSelected: {item['name']}")
        print(f"Formula: {item['formula']}")
        print("Provide the variables (press Enter to use typical default where noted).")
        # prompt variables
        var_list = item['vars']
        vars_provided = prompt_for_vars(var_list)

        # For variables where user left blank but a default exists, do not override:
        # functions expect certain keys; to allow defaults many functions check .get
        try:
            result = item['func'](vars_provided)
        except KeyError as e:
            print(f"Missing required variable: {e}. Please try again.")
            continue
        except ZeroDivisionError:
            print("Math error: division by zero (check inputs).")
            continue
        except Exception as e:
            print(f"Error computing result: {e}")
            continue

        unit = item.get('unit', '')
        name = item['name']
        # special message for NaN or inf
        if isinstance(result, float) and (math.isnan(result) or math.isinf(result)):
            print(f"\n{name} = {result} ({unit}) — result is not a real finite number; check inputs.\n")
        else:
            show_result(name, result, unit)

        # loop to continue or quit
        again = input("Compute another? (Y/n): ").strip().lower()
        if again in ('n', 'no'):
            print("Goodbye.")
            break
        print("\n" + "-"*40 + "\n")

if __name__ == "__main__":
    main()