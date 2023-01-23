# Geode Finder

High density geode and budding amethyst finder for Minecraft. This tool searches through a minecraft world for areas with large amounts of potential locations for geodes to generate and runs a simulation for budding amethyst generation in each area.

## Installation

### Github Releases
(WIP)

### Cargo
(WIP)

## Usage

0. This guide currently assumes that you have a working Rust tool chain with Cargo and some knowledge of modifying Rust code.

2. Compile the binary with `cargo build --release`. The executable will be located in the `target/releases` directory.

1. Use the `--help` argument to see a list of arguments you can modify such as game version, seed, search radius, geode threshold and budding amethyst threshold. The geode threshold and budding amethyst threshold will need to be set depending on your search radius, but good defaults for 10,000 chunk search radius are a geode threshold of 26-27 and a budding amethyst threshold of 1000.


## Todo

- [x] Add user input and CLI arguments
- [x] Add support for 1.17
- [ ] Publish executables for Linux and Windows
- [x] Create variables for custom geode feature configuration
- [x] Improve search algorithm for geode search to minimize repeated checks
- [ ] Add carpet script for automating world generation

## Credits

- [KaptainWutax](https://github.com/KaptainWutax) for reference implementation of a geode finder
- [MrSpike](https://github.com/MrSpike63) for some help with logic and information about geode generation
