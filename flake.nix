{
  description = "C development environment with OpenMP and MPI support";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in
      {
        packages.default = pkgs.stdenv.mkDerivation {
          pname = "heat_parallel";
          version = "1.0";
          src = ./.;

          buildInputs = [ pkgs.gcc pkgs.openmpi ];

          buildPhase = ''
            # mpicc -o heat_parallel heat_parallel.c
            make all
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp heat_parallel $out/bin/
          '';
        };


        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            gcc
            openmpi
            autoconf
            automake
            libtool
            pkg-config
            bear
          ];
          shellHook = ''
            # Use bear to generate compile_commands.json for LSP
            bear -- mpicc -o heat_parallel heat_parallel.c
            echo "Development environment with OpenMP and MPI support is ready."
          '';
        };
      }
    );
}

