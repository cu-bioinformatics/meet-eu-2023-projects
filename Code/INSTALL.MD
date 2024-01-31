# Install 

## Requirements
- Python (3.10.13)
- Java (min ??, max 20)
- WSL
- Pymol

## 0 - Visualization
- Create a conda environment with Pymol.
- RDKit.

## 1 - POCKET DETECTION: P2RANK
You need:
1. Use gitbash or Install Cygwin (https://www.cygwin.com/)
2. Install Java (jdk20)
3. Download P2Rank (https://github.com/rdk/p2rank)
4. Run make.sh script
   - There's probably a problem with gradlew, which is a Gradle dependency. To solve it, follow these steps:
   A. Check which version corresponds to your Java (see here: https://docs.gradle.org/current/userguide/compatibility.html) and memorize the appropriate one for compilation/testing.
   B. Then go to gradle>wrapper>gradle-wrapper.properties and change the distUrl version number to the one you found (e.g. if you have "distributionUrl=https\://services.gradle.org/distributions/gradle-8.5-bin.zip", change it to "distributionUrl=https\://services.gradle.org/distributions/gradle-8.1-bin.zip").
   C. Then run `./gradlew clean`.
5. Copy the `local_env` variable from `misc` to get the `./prank.sh` command.
6. Example command to see if it works:
   ``bash
   ./prank.sh predict -f ../NSP13/6zsl.cif -o output_here ../NSP13/results/test.ds

# 2 - LIGAND PRESELECTION

- The script is in `ligands>ligands_down...>lig_down_script`.