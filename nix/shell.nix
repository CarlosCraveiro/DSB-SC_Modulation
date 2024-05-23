{ pkgs ? import <nixpkgs> { } }:
with pkgs;
mkShell {
    buildInputs = [
        octave
    ];
    shellHook = ''
        #...
    '';
}
