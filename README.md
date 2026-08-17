# Resistor Network Analysis

A computational project analyzing **resistive grid networks using linear algebra and circuit theory**. I modeled a 5×5 resistor network as a system of equations and implemented LU factorization to solve for node voltages, currents, power dissipation, and effective resistance.

## How It Works

- Constructed a conductance matrix using Kirchhoff’s Current Law
- Applied fixed-voltage boundary conditions to the resistor network
- Implemented LU factorization with forward and backward substitution
- Calculated node voltages and currents through each resistor
- Analyzed power dissipation and current as network parameters changed
- Used the graph Laplacian pseudoinverse to calculate effective resistance

## Technologies & Concepts

**Linear Algebra • Circuit Analysis • LU Factorization • Graph Laplacians • Numerical Methods • KCL • Effective Resistance**

## Project Report

For a detailed explanation of the methodology, mathematical formulation, and results, see the [Resistor Network Analysis Report](./Resistor_Network_Analysis_Report.pdf).

## What I Learned

This project strengthened my understanding of how **physical systems can be represented and solved computationally using matrix methods**. I gained experience connecting circuit theory with linear algebra, implementing numerical solvers, and analyzing how network properties change with resistance and grid size.

## Future Improvements

- Extend the model to non-uniform and random resistor networks
- Compare the custom LU solver with iterative numerical methods
- Scale the analysis to significantly larger networks
