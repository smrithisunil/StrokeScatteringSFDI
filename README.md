## Purpose

This README provides a brief documentation for the code and data provided in this repository that is linked with the following paper submission:
Stroke core revealed by tissue scattering using spatial frequency domain imaging
Smrithi Sunil, Sefik Evren Erdener, Xiaojun Cheng, Sreekanth Kura, Jianbo Tang, John Jiang, Kavon Karrobi, Kivilcim Kiliç, Darren Roblyer, and David A. Boas

## Directory Structure

### …/figuresData
- Contains main script, generateFigures.m, to generate the figures from the manuscript
- Contains all the relevant data needed to generate the figures
- Reproduce figures by running generateFigures.m in MATLAB

### …/expData/analysisCode
- Contains code used to analyze OCT angiogram raw data to obtain OCT attenuation measures (run OCTanalysis.m)
- Contains code used to analyze SFDI scattering maps using a Monte Carlo lookup table (run SFDIanalysis.m)
- Contains code used to analyze SFDI scattering and OCT attenuation overlap (run contourAnalysis.m)

### …/expData/MX
- Data from seven individual mice, one folder for each mouse 
- Data contains OCT attenuation and SFDI scattering maps before stroke and at various time points after stroke (1 hour, 2 hours, 4 hours, 24 hours, and 72 hours)
- Data contains calculated overlap contour maps for each time point between OCT and SFDI
- Data also contains TTC image and corresponding overlap with OCT and SFDI
- All data is processed and registered to the pre-stroke SFDI scattering map
