{ pkgs ? import <nixpkgs> { } }:
with pkgs;
mkShell {
    buildInputs = [
        ghostscript
        octaveFull 
    ];
    shellHook = ''
        #...
    '';
}
