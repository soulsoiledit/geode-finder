# Geode Finder

This program helps find high density areas of geodes and budding amethyst in a given Minecraft world. It runs a search through each chunk within the search radius, and finds areas exceeding a given threshold of geodes within the random tick range. Then, it runs a simulation of geode feature generation, including budding amethysts, to filter out areas failing a budding amethyst threshold. Once complete, the program will return a list of coordinates of the center chunks in each location.

## Installation

### GitHub Releases

1. An executable for Windows is located in GitHub releases.

### Cargo

This guide assumes that you already have Rust and Cargo installed and working.

1. Use `cargo install geode-finder` to compile and install the program.

2. The executable will be located your cargo binary directory.

## Usage

1. To see a list of all available arguments and defaults, run the executable with the `--help` argument. The defaults are helpful to minimize the amount of searching done.

2. Set the variables accordingly and run the program (`./geode-finder`):
- `--seed`: the seed of the world that you want to search in (default: 177013)
- `--search-radius`: the radius of chunks to search through (default: 10000)
- `--geode-threshold`: the minimum number of geodes in random tick range (default: 26)
- `--budding-threshold`: the minimum number of budding amethyst in random tick range (default: 1000)
- `--game-version`: the game version to use (1.18 or 1.17)

3. The program will produce a list of coordinates of the center chunks of each valid location.

4. If you have Carpet mod installed, I've included a helper script to facilitate world pregeneration. Copy the list of valid geode locations to `[worldname]/scripts/shared/geodes.txt` and place `geodegen.sc` inside the `[worldname]/scripts` directory. Load the script and begin the search with `./geodegen`.

5. You can now use the region files from your world in this [Geode AFK Spot Finder](https://russellsprouts.github.io/minecraft-amethyst-tool/) tool to obtain the best locations. 

## Todo

- [x] Add user input and CLI arguments
- [x] Add support for 1.17
- [x] Publish executables for Linux and Windows
- [x] Create variables for custom geode feature configuration
- [x] Improve search algorithm for geode search to minimize repeated checks
- [x] Add carpet script for automating world generation

## Credits

- [KaptainWutax](https://github.com/KaptainWutax) for reference implementation of a geode finder
- [MrSpike](https://github.com/MrSpike63) for some help with logic and information about geode generation
