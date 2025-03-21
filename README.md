# 2D SPH Shallow Water Equations Simulator

A dam break simulation using Smoothed Particle Hydrodynamics (SPH) for solving the Shallow Water Equations (SWE) built with Rust and the Bevy game engine.

## Overview

This project implements a 2D fluid simulation using the SPH method to model the shallow water equations, which are commonly used to simulate fluid dynamics in situations where the horizontal scale is much larger than the vertical scale (like floods, tsunamis, or dam breaks).

The simulation uses particles to represent the water and calculates interactions between them to simulate fluid behavior, including:
- Water depth calculations
- Adaptive smoothing lengths
- Velocity and position updates based on elevation gradients
- Surface friction

## Features

- Real-time interactive dam break simulation
- Physics-based water flow over arbitrary terrain
- Adjustable simulation parameters
- Visualization tools for water depth and flow

## Requirements

- Rust (2024 Edition)
- Cargo (Rust's package manager)
- Dependencies:
  - Bevy (v0.15.3) - Game engine for rendering and ECS
  - Various Bevy plugins for debugging and asset loading

## Building and Running

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/twod-sph-swe.git
   cd twod-sph-swe
   ```

2. **Build the project:**
   ```bash
   cargo build
   ```

3. **Run the simulation:**
   ```bash
   cargo run
   ```

## Input Files

The simulation requires two image files to be placed in the "assets" directory:

1. **dam.png** - Defines the initial water location (red pixels)
2. **elevation.png** - Defines the terrain elevation (grayscale, where brighter pixels represent higher elevations)

## Utility Tools

The `utils` directory contains Python scripts to help create and manipulate input files:

- **draw_dam.py**: A GUI tool to create and edit dam layouts
  ```bash
  python utils/draw_dam.py input_image.png output_image.png threshold
  ```

- **geotiff_to_png.py**: Converts GeoTIFF elevation data to the required PNG format
  ```bash
  python utils/geotiff_to_png.py scale input.tif output.png
  ```

## How It Works

### SPH Method

Smoothed Particle Hydrodynamics (SPH) is a computational method used to simulate fluid flows. Unlike grid-based methods, SPH is a meshless Lagrangian method where the fluid is represented by a set of particles.

In this implementation:
1. Each particle represents a small volume of water
2. Particles interact with neighboring particles within a smoothing length
3. The smoothing length adapts based on water depth
4. Forces applied include gravity, terrain gradient, and friction

### Shallow Water Equations

The Shallow Water Equations (SWE) are a set of hyperbolic partial differential equations that describe the flow of a fluid under the influence of gravity, assuming that the depth of the fluid is small compared to the wavelengths of disturbances.

Key parameters in the simulation:
- `GRAVITY`: Gravitational constant (9.8 m/s²)
- `CFL`: Courant–Friedrichs–Lewy condition for stability
- `ROUGHNESS_COEFF`: Manning's roughness coefficient

## Simulation Parameters

You can adjust the following parameters in `main.rs`:

- `SCALE`: Controls the resolution of the simulation (default: 30)
- `PARTICLES_PER_PIXEL`: Number of particles per pixel (default: 3)
- `TOTAL_VOLUME`: Total water volume in the simulation
- `DAM_COLOR`: Color used in the dam image to represent water

## Debugging

The simulation includes:
- FPS counter
- World inspector (provided by bevy-inspector-egui)
- Automatic pause when numerical instabilities are detected

## Acknowledgments

- The Bevy open-source community
- The implementation is based on the research paper:
  > Kao, H.M., Chang, T.J. (2012). Numerical modeling of dambreak-induced flood and inundation using smoothed particle hydrodynamics. Journal of Hydrology, 448-449, 232-244.
