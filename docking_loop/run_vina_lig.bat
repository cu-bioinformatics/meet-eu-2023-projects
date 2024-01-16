@echo off

for %%a in ("%~dp0\lig\*.pdbqt") do (
"C:\Program Files (x86)\The Scripps Research Institute\Vina\vina.exe" --receptor protein.pdbqt --ligand "%%a" --config config.txt --log "%%~dpna_log.txt"
)

pause
