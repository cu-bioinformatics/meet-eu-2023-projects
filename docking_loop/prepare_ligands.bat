@echo off

for %%a in ("%~dp0\lig\*.pdb") do (
"C:\Program Files (x86)\MGLTools-1.5.6\python.exe" "C:\Program Files (x86)\MGLTools-1.5.6\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py" -l "%%a" -v
)

pause