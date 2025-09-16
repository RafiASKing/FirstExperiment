# FirstExperiment

Small experiments and utilities.

highschool_physics.py
---------------------

A simple, ready-to-use physics calculator for common high-school equations.
It can run interactively (menu + prompts) or via a command-line interface.

Quick start (Windows PowerShell)
--------------------------------

- List available equations
```
python .\highschool_physics.py --list
```

- Compute directly from CLI (no prompts). Example: Kinetic Energy (id 1) with m=2 kg, v=3 m/s
```
python .\highschool_physics.py --eq 1 --set m=2 --set v=3
```

- Interactive mode (menu and prompts)
```
python .\highschool_physics.py
```

Options
-------

- `--list`: Show equations and exit.
- `--eq <id>`: Select equation by id.
- `--set key=value`: Provide variables for non-interactive mode (repeatable).
- `--precision <n>`: Significant digits in output (default 6).
- `--interactive`: Force interactive prompts.

