$env:Path = "C:\ProgramData\miniconda3\envs\qgis;C:\ProgramData\miniconda3\envs\qgis\Scripts;$env:Path"
$env:CONDA_PREFIX = "C:\ProgramData\miniconda3\envs\qgis"
$env:CONDA_DEFAULT_ENV = "qgis"

Write-Host "Using Python: C:\ProgramData\miniconda3\envs\qgis\python.exe"
& "C:\ProgramData\miniconda3\envs\qgis\python.exe" --version
