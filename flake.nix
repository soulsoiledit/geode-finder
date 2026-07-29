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
        f: nixpkgs.lib.genAttrs supportedSystems (system: f nixpkgs.legacyPackages.${system} system);
    in
    {
      packages = forAllSystems (
        pkgs: system:
        let
          craneLib = crane.mkLib pkgs;

          inherit (craneLib.crateNameFromCargoToml { cargoToml = ./Cargo.toml; }) pname version;
          src = craneLib.cleanCargoSource ./.;

          mkRustPackage =
            {
              lib,
              suffix ? "",
              postInstall ? "",
            }:
            let
              commonArgs = {
                inherit pname version src;
                strictDeps = true;
              };
            in
            lib.buildPackage (
              commonArgs
              // {
                inherit postInstall;
                cargoArtifacts = lib.buildDepsOnly commonArgs;
                pname = "${pname}-${suffix}";
              }
            );
        in
        {
          nix = mkRustPackage { lib = craneLib; };

          static = mkRustPackage {
            lib = crane.mkLib pkgs.pkgsStatic;
            suffix = "static";
            postInstall = ''
              mv $out/bin/${pname} $out/bin/${pname}-static
            '';
          };

          windows = mkRustPackage {
            lib = crane.mkLib pkgs.pkgsCross.mingwW64;
            suffix = "windows";
          };

          default = self.packages.${system}.nix;
          all = pkgs.symlinkJoin {
            name = "${pname}-all-${version}";
            paths = with self.packages.${system}; [
              nix
              static
              windows
            ];
          };
        }
      );

      devShells = forAllSystems (
        pkgs: system: {
          default = pkgs.mkShell {
            packages = with pkgs; [
              rustc
              cargo
              rustfmt
              clippy
            ];
          };
        }
      );
    };
}
