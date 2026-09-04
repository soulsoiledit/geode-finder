{
  description = "A geode finder";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    crane.url = "github:ipetkov/crane";
  };

  outputs =
    {
      self,
      nixpkgs,
      crane,
    }:
    let
      supportedSystems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
      ];
      forAllSystems =
        f: nixpkgs.lib.genAttrs supportedSystems (system: f nixpkgs.legacyPackages.${system});
      mkCrate =
        pkgs:
        let
          craneLib = crane.mkLib pkgs;
          commonArgs = {
            src = craneLib.cleanCargoSource ./.;
            strictDeps = true;
          };
        in
        {
          inherit craneLib;
          crate = craneLib.buildPackage (
            commonArgs
            // {
              cargoArtifacts = craneLib.buildDepsOnly commonArgs;
            }
          );
        };
    in
    {
      checks = forAllSystems (pkgs: {
        inherit (mkCrate pkgs) crate;
      });

      packages = forAllSystems (pkgs: {
        default = (mkCrate pkgs).crate;
      });

      devShells = forAllSystems (pkgs: {
        default = (mkCrate pkgs).craneLib.devShell {
          checks = self.checks.${pkgs.stdenv.hostPlatform.system};
          packages = with pkgs; [
            # please please please benchmark and profile
            hyperfine
            poop

            perf
            cargo-flamegraph
          ];
        };
      });
    };
}
