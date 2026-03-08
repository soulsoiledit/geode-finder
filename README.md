# Geode Finder

This programs helps find high density clusters of geodes and budding amethyst in a Minecraft world without running real world generation. It runs a search through each chunk within the search radius and finds areas exceeding the selected threshold of geodes within the loaded radius. Then, it goes through each area and fully simulates the geode generation to find the areas that also exceed the budding amethyst threshold. The program both prints out a list of coordinates of central chunks of each cluster and saves the results to an output file.

## Installation

### GitHub Releases

1. Executables for Windows and Linux are located in GitHub releases.

### Cargo

1. Install and set up a Rust development environment.

1. Use `cargo install geode-finder` to download, build, and install the program.

1. The executable will be located your cargo binary directory.

### Source

1. Install and set up Rust development environment.

1. Download or clone this repository.

1. Use `cargo build --release` to build the program.

1. The executable will be located in the `./target/release/` directory.

## Usage

1. Open a terminal and navigate to the folder containing the binary. On Windows, you can navigate to the folder with the .exe in File Explorer, right click, and select "Open in Terminal".

1. Run `./geode-finder.exe --help` or `./geode-finder --help` to view the list of options:

    | **Option**                  | **Information**                                            | **Default** |
    | --------------------------- | ---------------------------------------------------------- | ----------- |
    | `-h`, `--help`              | Shows the list of options                                  |
    | `-V`, `--version`           | Prints the version of the program                          |
    | `-m`, `--minecraft-version` | Minecraft version to use (1.17, 1.18, or 1.19)             | 1.19        |
    | `-s`, `--seed`              | Seed of your world                                         | 0           |
    | `-r`, `--search-radius`     | Radius of chunks to search over                            | 1000        |
    | `-g`, `--geode-threshold`   | Minimum number of geodes in a cluster                      | 20          |
    | `-b`, `--budding-threshold` | Minimum number of budding amethyst in a cluster            | 800         |
    | `--loaded-radius`           | Random tickable radius                                     | 6           |
    | `--center-x`                | x coordinate of the center chunk                           | 0           |
    | `--center-z`                | z coordinate of the center chunk                           | 0           |
    | `--output-path`             | Where to save results                                      | output.json |
    | `--estimate`                | Estimate the number of clusters without running the search | false       |

1. Run the program with your selected options.

    Example: `./geode-finder.exe -r 10000 -g 22 -b 850` will search a 10000 chunk radius around (0, 0) to find 13x13 (`--loaded-radius` \* 2 + 1) clusters that have at least 22 geodes at 850 budding amethyst.

1. The program will find all valid geode clusters at the beginning, then check if each cluster will meet the budding amethyst threshold. It will print the clusters out and save them to file selected by `--output-path`.

1. If you have Carpet mod installed, I've included a helper script to facilitate world pregeneration. Copy the printed list of valid geode locations to `[world]/scripts/shared/geodes.txt` and place `geodegen.sc` inside the `[worldname]/scripts` directory. Load the script and begin the search with `./geodegen`. This script currently does not work with the `--output-path` file.

## Credits

- [KaptainWutax](https://github.com/KaptainWutax) for reference implementation of a geode finder
- [MrSpike](https://github.com/MrSpike63) for some help with logic and information about geode generation
