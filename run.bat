@echo off
SETLOCAL ENABLEDELAYEDEXPANSION

for /f "delims=" %%a in ('where /R "C:\Program Files\R" "Rscript.exe"2^>nul') do set p=%%~fa
if defined p (
echo %p%
) else (
echo File not found
)
SET ROPTS=--no-save --no-environ --no-init-file --no-restore --no-Rconsole
"%p%" runShinyApp.R 2>&1