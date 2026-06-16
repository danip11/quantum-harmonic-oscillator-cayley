````markdown
# Quantum Harmonic Oscillator Simulation with the Cayley Method

This repository contains a computational physics project focused on the numerical simulation of the time-dependent Schrödinger equation for one-dimensional quantum systems.

The code is written in Fortran and uses an implicit Cayley scheme together with finite differences to study the time evolution of quantum wave functions.

## Project overview

The main objective of this project is to simulate and analyse the dynamics of confined quantum systems. In particular, the project includes:

- Time evolution of quantum wave functions using the Cayley method.
- Simulation of the quantum harmonic oscillator.
- Evolution of Gaussian wave packets.
- Conservation analysis of the norm and mean energy.
- Calculation of expected values of position and momentum.
- Study of the uncertainty principle.
- Quantum-classical comparison for the harmonic oscillator.
- Simulation of tunnelling in a double-well potential.

## Numerical method

The time evolution is based on the Cayley approximation of the unitary evolution operator. This method provides a stable and norm-conserving scheme for solving the time-dependent Schrödinger equation.

The Hamiltonian is discretised using finite differences, leading to a tridiagonal system that is solved at each time step.

## Repository structure

```text
src/
    gaussiana_cayley.f90
    oscilador_cayley.f90
    autofuncion_caja_cayley.f90

docs/
    informe_oscilador_armonico_cuantico.pdf
````

## Requirements

To compile and run the Fortran codes, a Fortran compiler is required. For example:

```bash
gfortran
```

## Compilation example

```bash
gfortran -O2 -o gaussiana src/gaussiana_cayley.f90
./gaussiana
```

The programs generate output files with numerical data that can be used for plotting the wave function, probability density, norm, energy and expectation values.

## Author

Daniel Pérez Pérez
MSc in Physics, University of Granada
GitHub: https://github.com/danip11

```
```
