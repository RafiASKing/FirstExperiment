#!/usr/bin/env python3
"""highschool_physics
----------------------

A small, ready-to-use physics helper for common high-school level
equations. Can run interactively or via a simple CLI for quick
calculations.

Examples
========

- List available equations:
    python highschool_physics.py --list

- Compute directly from the CLI (no prompts):
    python highschool_physics.py --eq 1 --set m=2 --set v=3

- Interactive menu:
    python highschool_physics.py
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import Callable, Dict, List, Mapping, MutableMapping, Sequence, Tuple


# Physical constants
G = 6.67430e-11          # gravitational constant, m^3 kg^-1 s^-2
g0 = 9.81                # standard gravity, m/s^2
c = 299_792_458          # speed of light m/s
h = 6.62607015e-34       # Planck's constant J*s
k_e = 8.9875517923e9     # Coulomb constant N*m^2/C^2
PI = math.pi

Params = Mapping[str, float]


def show_result(name: str, value: float, unit: str = "", precision: int = 6) -> None:
    """Print a neat formatted result with configurable precision."""
    fmt = f"{{value:.{precision}g}}"
    print(f"\n{name} = {fmt.format(value=value)} {unit}\n")


# Equation implementations
def kinetic_energy(params: Params) -> float:
    """Kinetic energy: KE = 1/2 m v^2"""
    m = params["m"]; v = params["v"]
    return 0.5 * m * v * v

def gravitational_pe(params: Params) -> float:
    """Gravitational potential energy near Earth: U = m g h"""
    m = params["m"]; h_ = params["h"]; g = params.get("g", g0)
    return m * g * h_

def momentum(params: Params) -> float:
    """Linear momentum: p = m v"""
    return params["m"] * params["v"]

def newton_f(params: Params) -> float:
    """Newton's 2nd law (Force): F = m a"""
    return params["m"] * params["a"]

def work_force_distance(params: Params) -> float:
    """Work by a constant force: W = F d cos(theta) (theta in degrees)"""
    F = params["F"]; d = params["d"]; theta_deg = params.get("theta", 0.0)
    theta_rad = math.radians(theta_deg)
    return F * d * math.cos(theta_rad)

def power_from_work(params: Params) -> float:
    """Power: P = W / t"""
    return params["W"] / params["t"]

def final_velocity(params: Params) -> float:
    """Final velocity: v = u + a t"""
    return params["u"] + params["a"] * params["t"]

def displacement_suvat(params: Params) -> float:
    """Displacement: s = u t + 1/2 a t^2"""
    return params["u"] * params["t"] + 0.5 * params["a"] * params["t"] ** 2

def final_velocity_suvat(params: Params) -> float:
    """Final speed from v^2 = u^2 + 2 a s (positive root)."""
    u = params["u"]; a = params["a"]; s = params["s"]
    val = u * u + 2 * a * s
    if val < 0:
        return float("nan")
    return math.sqrt(val)

def centripetal_acc(params: Params) -> float:
    """Centripetal acceleration: a_c = v^2 / r"""
    v = params["v"]; r = params["r"]
    return (v * v) / r

def centripetal_force(params: Params) -> float:
    """Centripetal force: F = m v^2 / r"""
    return params["m"] * params["v"] ** 2 / params["r"]

def pressure(params: Params) -> float:
    """Pressure: P = F / A"""
    return params["F"] / params["A"]

def density(params: Params) -> float:
    """Density: rho = m / V"""
    return params["m"] / params["V"]

def hooke_force(params: Params) -> float:
    """Hooke's law (spring force): F = k x"""
    return params["k"] * params["x"]

def period_mass_spring(params: Params) -> float:
    """Period of mass-spring: T = 2π√(m/k)"""
    return 2 * PI * math.sqrt(params["m"] / params["k"])

def period_pendulum(params: Params) -> float:
    """Period of simple pendulum (small angles): T = 2π√(L/g)"""
    return 2 * PI * math.sqrt(params["L"] / params.get("g", g0))

def wave_speed(params: Params) -> float:
    """Wave speed: v = f λ"""
    return params["f"] * params["lam"]

def photon_energy(params: Params) -> float:
    """Photon energy: E = h f"""
    return h * params["f"]

def ohms_law(params: Params) -> float:
    """Ohm's law: V = I R"""
    return params["I"] * params["R"]

def electrical_power(params: Params) -> float:
    """Electrical power: P = V I"""
    return params["V"] * params["I"]

def coulomb_force(params: Params) -> float:
    """Coulomb's law: F = k q1 q2 / r^2"""
    q1 = params["q1"]; q2 = params["q2"]; r = params["r"]
    return k_e * q1 * q2 / (r * r)

@dataclass(frozen=True)
class EquationSpec:
    id: int
    name: str
    variables: List[Tuple[str, str]]
    func: Callable[[Params], float]
    unit: str
    formula: str


# Catalog of equations: id -> EquationSpec
EQUATIONS: Dict[int, EquationSpec] = {
    1: EquationSpec(1, "Kinetic energy", [("m", "mass (kg)"), ("v", "speed (m/s)")], kinetic_energy, "J", "KE = 1/2 m v^2"),
    2: EquationSpec(2, "Gravitational potential energy (near Earth)", [("m", "mass (kg)"), ("h", "height (m)"), ("g", "gravity (m/s^2) [optional; default 9.81]")], gravitational_pe, "J", "U = m g h"),
    3: EquationSpec(3, "Linear momentum", [("m", "mass (kg)"), ("v", "velocity (m/s)")], momentum, "kg·m/s", "p = m v"),
    4: EquationSpec(4, "Newton's 2nd law (Force)", [("m", "mass (kg)"), ("a", "acceleration (m/s^2)")], newton_f, "N", "F = m a"),
    5: EquationSpec(5, "Work by a constant force", [("F", "force (N)"), ("d", "displacement (m)"), ("theta", "angle between F and d (deg) [optional; default 0]")], work_force_distance, "J", "W = F d cos(theta)"),
    6: EquationSpec(6, "Power from work", [("W", "work (J)"), ("t", "time (s)")], power_from_work, "W", "P = W / t"),
    7: EquationSpec(7, "Final velocity (v = u + a t)", [("u", "initial velocity (m/s)"), ("a", "acceleration (m/s^2)"), ("t", "time (s)")], final_velocity, "m/s", "v = u + a t"),
    8: EquationSpec(8, "Displacement (s = u t + 1/2 a t^2)", [("u", "initial velocity (m/s)"), ("a", "acceleration (m/s^2)"), ("t", "time (s)")], displacement_suvat, "m", "s = u t + 1/2 a t^2"),
    9: EquationSpec(9, "Final speed from v^2 = u^2 + 2 a s", [("u", "initial speed (m/s)"), ("a", "acceleration (m/s^2)"), ("s", "displacement (m)")], final_velocity_suvat, "m/s", "v^2 = u^2 + 2 a s"),
    10: EquationSpec(10, "Centripetal acceleration", [("v", "speed (m/s)"), ("r", "radius (m)")], centripetal_acc, "m/s^2", "a_c = v^2 / r"),
    11: EquationSpec(11, "Centripetal force", [("m", "mass (kg)"), ("v", "speed (m/s)"), ("r", "radius (m)")], centripetal_force, "N", "F = m v^2 / r"),
    12: EquationSpec(12, "Pressure", [("F", "force (N)"), ("A", "area (m^2)")], pressure, "Pa", "P = F / A"),
    13: EquationSpec(13, "Density", [("m", "mass (kg)"), ("V", "volume (m^3)")], density, "kg/m^3", "rho = m / V"),
    14: EquationSpec(14, "Hooke's law (spring force)", [("k", "spring constant (N/m)"), ("x", "displacement (m)")], hooke_force, "N", "F = k x"),
    15: EquationSpec(15, "Period of mass-spring (T = 2π√(m/k))", [("m", "mass (kg)"), ("k", "spring constant (N/m)")], period_mass_spring, "s", "T = 2π√(m/k)"),
    16: EquationSpec(16, "Period of simple pendulum (small angles)", [("L", "length (m)"), ("g", "gravity (m/s^2) [optional; default 9.81]")], period_pendulum, "s", "T = 2π√(L/g)"),
    17: EquationSpec(17, "Wave speed (v = f λ)", [("f", "frequency (Hz)"), ("lam", "wavelength (m)")], wave_speed, "m/s", "v = f λ"),
    18: EquationSpec(18, "Photon energy (E = h f)", [("f", "frequency (Hz)")], photon_energy, "J", "E = h f"),
    19: EquationSpec(19, "Ohm's law (V = I R)", [("I", "current (A)"), ("R", "resistance (Ω)")], ohms_law, "V", "V = I R"),
    20: EquationSpec(20, "Electrical power (P = V I)", [("V", "voltage (V)"), ("I", "current (A)")], electrical_power, "W", "P = V I"),
    21: EquationSpec(21, "Coulomb's law (electrostatic force)", [("q1", "charge 1 (C)"), ("q2", "charge 2 (C)"), ("r", "distance (m)")], coulomb_force, "N", "F = k q1 q2 / r^2"),
}

__all__ = [
    "G", "g0", "c", "h", "k_e", "PI",
    "kinetic_energy", "gravitational_pe", "momentum", "newton_f",
    "work_force_distance", "power_from_work", "final_velocity",
    "displacement_suvat", "final_velocity_suvat", "centripetal_acc",
    "centripetal_force", "pressure", "density", "hooke_force",
    "period_mass_spring", "period_pendulum", "wave_speed", "photon_energy",
    "ohms_law", "electrical_power", "coulomb_force", "EQUATIONS",
]


def prompt_for_vars(var_list: Sequence[Tuple[str, str]]) -> Dict[str, float]:
    """Prompt the user for each variable described in var_list.

    Accepts empty string to allow defaults handled inside the equation
    functions (e.g., g = 9.81 if omitted).
    """
    values: Dict[str, float] = {}
    for key, prompt_text in var_list:
        while True:
            s = input(f"Enter {prompt_text}: ").strip()
            if s == "":
                break  # leave absent so the function may use defaults
            if s.lower() in ("q", "quit", "exit"):
                print("Exiting.")
                sys.exit(0)
            try:
                values[key] = float(s)
                break
            except ValueError:
                print("Invalid number entered. Please try again.")
    return values


def list_equations() -> None:
    print("Equations:")
    for idx in sorted(EQUATIONS.keys()):
        item = EQUATIONS[idx]
        print(f"{idx:2d}. {item.name}")


def parse_kv_pairs(pairs: Sequence[str]) -> Dict[str, float]:
    """Parse repeated --set key=value pairs into a {key: float} mapping."""
    result: Dict[str, float] = {}
    for p in pairs:
        if "=" not in p:
            raise ValueError(f"Invalid --set entry '{p}'. Expected key=value.")
        k, v = p.split("=", 1)
        try:
            result[k] = float(v)
        except ValueError:
            raise ValueError(f"Value for '{k}' must be a number; got '{v}'.")
    return result


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="High-school physics calculator")
    parser.add_argument("--list", action="store_true", help="List available equations and exit")
    parser.add_argument("--eq", type=int, help="Equation ID to compute (use --list to see options)")
    parser.add_argument("--set", dest="sets", action="append", default=[], help="Set a variable as key=value (repeatable)")
    parser.add_argument("--precision", type=int, default=6, help="Significant digits in output (default: 6)")
    parser.add_argument("--interactive", action="store_true", help="Force interactive mode (menu and prompts)")
    return parser

def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.list:
        list_equations()
        return 0

    # Non-interactive CLI compute
    if args.eq and not args.interactive:
        if args.eq not in EQUATIONS:
            print(f"Unknown equation id: {args.eq}. Use --list to see options.")
            return 2
        spec = EQUATIONS[args.eq]
        try:
            params = parse_kv_pairs(args.sets)
            result = spec.func(params)
        except KeyError as e:
            print(f"Missing required variable: {e} for '{spec.name}'.")
            return 2
        except ValueError as e:
            print(str(e))
            return 2
        except ZeroDivisionError:
            print("Math error: division by zero (check inputs).")
            return 2
        except Exception as e:
            print(f"Error computing result: {e}")
            return 2

        if isinstance(result, float) and (math.isnan(result) or math.isinf(result)):
            print(f"{spec.name} = {result} {spec.unit} — not a real finite number; check inputs.")
        else:
            show_result(spec.name, result, spec.unit, precision=args.precision)
        return 0

    # Interactive mode
    print("Simple Physics Calculator")
    print("Enter 'q' at any prompt to quit.\n")
    while True:
        list_equations()
        print(" 0. Quit")
        choice = input("\nPick an equation number and press Enter: ").strip()
        if choice.lower() in ("q", "quit", "exit", "0", ""):
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
        print(f"\nSelected: {item.name}")
        print(f"Formula: {item.formula}")
        print("Provide the variables (press Enter to use typical default where noted).")
        vars_provided = prompt_for_vars(item.variables)
        try:
            result = item.func(vars_provided)
        except KeyError as e:
            print(f"Missing required variable: {e}. Please try again.")
            continue
        except ZeroDivisionError:
            print("Math error: division by zero (check inputs).")
            continue
        except Exception as e:
            print(f"Error computing result: {e}")
            continue

        if isinstance(result, float) and (math.isnan(result) or math.isinf(result)):
            print(f"\n{item.name} = {result} ({item.unit}) — result is not a real finite number; check inputs.\n")
        else:
            show_result(item.name, result, item.unit, precision=args.precision)

        again = input("Compute another? (Y/n): ").strip().lower()
        if again in ("n", "no"):
            print("Goodbye.")
            break
        print("\n" + "-" * 40 + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())