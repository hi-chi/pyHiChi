Test configuration:
- photons emit from the origin in all directions
- G4_Galactic world
- spherical Au197 target
- spherical G4_Galactic detector behind the target
Photons may spawn other particle types.

Geant4 should be installed.

1. Build pyHiChi in a standart way.
2. Run geant4.bat or geant4.sh if necessary to use Geant4 environment variables.
3. Start Hi-Chi: "python generate_data.py"
   The script generates the "hichi_particles.txt" file
4. Build the test (run CMakeLists.txt in this directory).
5. Start the geant test program. TestGeant supports 2 modes: a visualization mode and a command line mode.
   You can start a single-threaded or a multithreaded launch. To obtain more information use "TestGeant --help".
   Example 1: "TestGeant.exe -r run.mac -p hichi_particles.txt" (cmd mode)
   Example 2: "TestGeant.exe -v init_vis.mac run.mac -p hichi_particles.txt" (visual mode, start vis.mac in Idle after modeling)
