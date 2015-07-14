@echo off
setlocal enabledelayedexpansion

set ifile=%1
set ifile=!ifile:\=/!

if exist ".\tsplot.sh" (
  rem If run from source directory
  bash -c "./tsplot.sh %ifile%"
) else (
  rem If source directory is contained in the path variable
  bash -c "tsplot.sh %ifile%"
)

