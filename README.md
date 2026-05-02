NOTE: This is a work in progress—functional, but not yet ready for production use.

`sdc_optimizer` is a preprocessor that takes a description of a superconducting digital circuit (SDC), modifies circuit parameters as directed, and runs a circuit simulator, specifically:

JoSIM: Josephson Junction Superconductive SPICE Circuit Simulator
Copyright (C) 2020 by Johannes Delport ([jdelport@sun.ac.za](mailto:jdelport@sun.ac.za))
v2.6.8.e80ec7c, compiled on Aug 28, 2024 at 12:46:48

The tool then analyzes the simulation results. Depending on the type of analysis, it may modify the circuit and repeat the simulation.

The input file is essentially a SPICE-like circuit description augmented with an application-specific language embedded as pragmas. These pragmas appear as comments to the SPICE interpreter.

Currently, only JoSIM is supported. It has a number of idiosyncrasies (such as capitalizing all characters and offering very limited parameterization). In principle, other SPICE-like simulators (e.g., WRspice) could also be supported.

Supported operations include:

* Plotting complex dependencies—for example, output delay as a function of competing inputs to investigate metastability issues
* Optimizing circuit component values using Bayesian optimization and a simulated annealing subsystem
* Determining and optimizing component value margins, such as finding optimal combinations in a multi-dimensional Schmoo plot
* Defining and verifying correct circuit behavior algorithmically

Many use cases involve thousands to millions of simulations, so analyses must be fully automated and cannot require manual intervention. This requirement is a major source of the framework’s complexity.

Although this work is aimed at SDCs, there is nothing inherently limiting it to that domain, and it could be applied to other types of circuits.

For more information, contact [agn@acm.org](mailto:agn@acm.org).
