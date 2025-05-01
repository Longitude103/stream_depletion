# Stream Depletion Calculator

## Description

The Stream Depletion Calculator is a Rust library designed to estimate the impact of groundwater pumping on streamflow. This library implements various methods for calculating stream depletion, including:

1. Glover Solution for Infinite Aquifers
2. Glover Solution for Alluvial Aquifers
3. Stream Depletion Factor (SDF) Method
4. Unit Response Function (URF) Method

These methods allow users to model the effects of groundwater extraction on nearby streams over time, considering different aquifer characteristics and pumping scenarios.

## Features

- Calculate stream depletion using multiple analytical methods
- Support for both infinite and alluvial aquifer scenarios
- Flexible input of pumping volumes on a monthly basis
- Time-series output of stream depletion estimates
- Conversion utilities for different units and time scales

## Usage

This package is not published yet to crates.io. To use at this point you must clone the repo and then add it to your "cargo.toml" file from your local computer.

Your toml file should look like the following:

```toml
[dependencies]
your-library = { path = "/path/to/your-library" }
```

## Installation

Once the project is added to your new app, you should be able to build and run it with the new library installed:
```bash
cd /home/user/projects/my-app
cargo build
cargo run
```

# Documentation

## Glover Equation with Image Wells for Alluvial Boundaries

### Overview
The Glover equation (Glover and Balmer, 1954) estimates streamflow depletion due to groundwater pumping in alluvial aquifers. For aquifers with finite boundaries, such as impermeable valley walls, the equation is modified using the method of image wells to account for boundary effects. This approach, first introduced by Charles V. Theis in 1941, simulates the hydraulic impact of boundaries by adding contributions from imaginary wells.
Standard Glover Equation
For an infinite aquifer, the streamflow depletion \( Q_s(t) \) is calculated as:

$$
Q_s(t) = Q_w \cdot \text{erfc}\left( \sqrt{\frac{S d^2}{4 T t}} \right)
$$

Where:

\( Q_s(t) \): Streamflow depletion rate (L³/T)
\( Q_w \): Pumping rate (L³/T)
\( S \): Specific yield or storativity (dimensionless)
\( T \): Transmissivity (L²/T)
\( d \): Distance from well to stream (L)
\( t \): Time since pumping began (T)
\( \text{erfc} \): Complementary error function

### Modification for Alluvial Boundaries
In alluvial aquifers bounded by impermeable features (e.g., valley walls), an image well is placed to simulate the boundary. For a single impermeable boundary at distance ( W ) from the well, perpendicular to the stream, the modified equation is:

$$
Q_s(t) = Q_w \cdot \left[ \text{erfc}\left( \sqrt{\frac{S d^2}{4 T t}} \right) + \text{erfc}\left( \sqrt{\frac{S (2W - d)^2}{4 T t}} \right) \right]
$$

Where:

\( W \): Distance from well to the boundary (L)
\( 2W - d \): Distance from image well to the stream (L)

Multiple boundaries require additional image wells, increasing computational complexity.

### Historical Context
Charles V. Theis (1941): First proposed image wells for groundwater depletion problems in The Effect of a Well on the Flow of a Nearby Stream, addressing stream and impermeable boundaries.
Glover and Balmer (1954): Developed the original equation for alluvial aquifers, later adapted with image wells for finite systems.
Jenkins (1968): Introduced the Stream Depletion Factor (SDF) but did not explicitly use image wells, assuming infinite aquifers.

### Applications
Estimating lagged depletions in water rights management (e.g., Colorado Water Court).
Modeling streamflow impacts in bounded alluvial aquifers.
Use in tools like streamDepletr (R package) and AQTESOLV for automated calculations.

### Limitations
Assumes homogeneous, isotropic aquifers and fully penetrating streams.
May underestimate long-term depletion in complex systems; validate with numerical models (e.g., MODFLOW) for critical applications.
Requires accurate aquifer parameters (( T ), ( S ), ( d ), ( W )).

### References
- Theis, C.V. (1941). The Effect of a Well on the Flow of a Nearby Stream. Transactions, American Geophysical Union, 22(3), 734–738.
- Glover, R.E., & Balmer, G.G. (1954). River Depletion Resulting from Pumping a Well near a River. Eos, Transactions American Geophysical Union, 35(3), 468–470.
- Jenkins, C.T. (1968). Techniques for Computing Rate and Volume of Stream Depletion by Wells. Ground Water, 6(2), 37–46.
- USGS Techniques of Water-Resources Investigations

## Jenkins Stream Depletion Factor (SDF) Method

### Overview
The Jenkins Stream Depletion Factor (SDF) method, introduced by Charles T. Jenkins in 1968, simplifies the calculation of streamflow depletion caused by groundwater pumping in alluvial aquifers. It builds on the Glover equation (1954) by introducing the SDF, a time parameter that quantifies the rate of depletion, making it easier to estimate impacts on streams for water management.
Stream Depletion Factor (SDF). The SDF is defined as:

$$ \text{sdf} = \frac{d^2 S}{T} $$

Where:

\( d \): Distance from the well to the stream (L)
\( S \): Specific yield (unconfined) or storativity (confined) (dimensionless)
\( T \): Transmissivity (L²/T)

The SDF (units: time) indicates how quickly pumping affects the stream. A lower SDF implies faster depletion due to proximity, high transmissivity, or low storage.

### Depletion Equation
Using the SDF, the streamflow depletion ( Q_s(t) ) is calculated as:

$$ Q_s(t) = Q_w \cdot \text{erfc}\left( \sqrt{\frac{\text{sdf}}{4 t}} \right) $$

Where:

- \( Q_s(t) \): Streamflow depletion rate (L³/T)
- \( Q_w \): Pumping rate (L³/T)
- \( t \): Time since pumping began (T)
- \( \text{erfc} \): Complementary error function

### Key Assumptions

- Infinite, homogeneous, isotropic aquifer.
- Fully penetrating stream with no streambed resistance.
- Horizontal flow (Dupuit assumptions).
- No external recharge or leakage.

### Applications

Estimating lagged streamflow depletions for water rights management (e.g., Colorado Water Court).
Assessing well pumping impacts in alluvial aquifers.

### Limitations

Assumes infinite aquifer extent; for finite boundaries (e.g., valley walls), image wells are needed (not covered in Jenkins’ original method).
May not account for streambed resistance or complex geology; validate with numerical models (e.g., MODFLOW) for critical cases.
Requires accurate aquifer parameters (( T ), ( S ), ( d )).

### Historical Context

Jenkins (1968): Introduced SDF in Techniques for Computing Rate and Volume of Stream Depletion by Wells.
Builds on Glover and Balmer (1954) and Theis (1941) for stream depletion modeling.
Extended by later works (e.g., streamDepletr) to include image wells for bounded aquifers.

References

- Jenkins, C.T. (1968). Techniques for Computing Rate and Volume of Stream Depletion by Wells. Ground Water, 6(2), 37–46.
- Glover, R.E., & Balmer, G.G. (1954). River Depletion Resulting from Pumping a Well near a River. Eos, Transactions American Geophysical Union, 35(3), 468–470.
- USGS Streamflow Depletion by Wells

## Contributing

You are welcome contributions to stream_depletion Library! Whether you're fixing bugs, adding features, improving documentation, or reporting issues, your help is greatly appreciated. This guide outlines how to contribute to the project.

### Getting Started
1. **Fork the Repository**: Click the "Fork" button on the [repository page](https://github.com/Longitude103/stream_depletion) to create your own copy.
2. **Clone Your Fork**:
   ```bash
   git clone https://github.com/Longitude103/stream_depletion.git
   cd stream_depletion
   ```
3. **Set Up the Environment**:
    - Ensure you have [Rust](https://www.rust-lang.org/tools/install) installed (`rustup` recommended).
    - Run `cargo build` to build the project.
    - Run `cargo test` to verify tests pass.

### How to Contribute
#### Reporting Issues
- Check the [issue tracker](https://github.com/Longitude103/stream_depletion/issues) to avoid duplicates.
- Open a new issue with a clear title, description, and steps to reproduce (if applicable).
- Use provided templates for bug reports or feature requests.

#### Submitting Code
1. **Create a Branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. **Follow Coding Standards**:
    - Adhere to Rust's idiomatic style (use `cargo fmt` for formatting).
    - Write clear, concise code with appropriate documentation (`///` comments for public APIs).
    - Ensure all tests pass (`cargo test`).
    - Use `cargo clippy` to catch common issues.
3. **Write Tests**:
    - Add unit tests in the `tests/` directory or inline with `#[test]`.
    - Ensure new features or bug fixes are covered by tests.
4. **Commit Changes**:
    - Use clear, descriptive commit messages (e.g., `Add feature X to module Y`).
    - Follow the [Conventional Commits](https://www.conventionalcommits.org/) format if applicable.
5. **Push and Create a Pull Request**:
   ```bash
   git push origin feature/your-feature-name
   ```
    - Open a pull request (PR) on the main repository.
    - Reference related issues (e.g., `Fixes #123`).
    - Describe the changes and their purpose in the PR description.
6. **Code Review**:
    - Respond to feedback from maintainers.
    - Make requested changes and push updates to the same branch.

#### Documentation
- Update `README.md` or other documentation for new features or changes.
- Add or improve doc comments for public APIs using `///`.
- Consider contributing to the `examples/` directory to demonstrate usage.

### Community Guidelines
- Be respectful and inclusive, following the [Rust Code of Conduct](https://www.rust-lang.org/policies/code-of-conduct).
- Engage constructively in discussions on issues and PRs.
- Ask questions if you're unsure—new contributors are welcome!

### Development Tools
- **Rust Version**: Use the version specified in `rust-toolchain.toml` or the latest stable release.
- **Formatter**: Run `cargo fmt` before committing.
- **Linter**: Run `cargo clippy --all-targets --all-features` to check for warnings.
- **Tests**: Use `cargo test` to run all tests.
- **Dependency Management**: Update dependencies with `cargo update` cautiously, ensuring compatibility.

### Contact
For questions or guidance, reach out via:
- [GitHub Issues](https://github.com/Longitude103/stream_depletion/issues)

Thank you for contributing to stream_depletion!

# License

stream_depletion is licensed under the MIT License, a permissive open-source license that allows you to use, modify, and distribute the software freely, provided that the license notice is included with all copies or substantial portions of the software.
