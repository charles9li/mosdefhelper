# mosdefhelper

**mosdefhelper** is a Python library that extends the capabilities of the Molecular Simulation Design Framework (MoSDeF).
It provides additional capabilities or building molecular systems, applying force fields, and managing simulation inputs in the MosDeF environment.

## Features

- **mbuildhelper**: Assists in building molecular systems with reusable components and utilities by extending the [mBuild](https://mbuild.mosdef.org/) package.
- **foyerhelper**: Simplifies force field application and management by extending the [Foyer](https://foyer.mosdef.org/) package.
- **mosdefhelper**: Integrates various functionalities for creating and managing simulation systems.

## Directory Overview

- `foyerhelper`: Tools for managing force fields using [Foyer](https://foyer.mosdef.org/).
  - Includes `forcefield.py` for force field-related logic.
  - Customizable force fields located in the `forcefields` directory.
  - Utilities and test cases provided in `utils` and `tests`.

- `mbuildhelper`: Extensions for [mBuild](https://mbuild.mosdef.org/) to facilitate system construction.
  - Includes libraries (`lib`) for reusable components.
  - Utilities and test cases in `utils` and `tests`.

- `mosdefhelper`: High-level tools for system creation and management.
  - Core functionality in `system.py`.
  - Includes test cases for robustness.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/charles9li/mosdefhelper.git
   cd mosdefhelper

2. Add the path of the package to `PYTHONPATH` in your `.bashrc` file:
   ```bash
   PYTHONPATH=${mosdefhelperdirectory}:$PYTHONPATH

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments

This library leverages the power of [mBuild](https://mbuild.mosdef.org/) and [Foyer](https://foyer.mosdef.org/).
