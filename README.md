# Muon Scattering / Absorption Analysis

This script analyzes the proportion of cosmic muons absorped/scattered from the data generated by Monte Carlo Simulation carried out using Geant4 and Cry.

## Purpose

Determining the proportions of “missing” muons in the MC due to scattering/absorption phenomena caused by the Water Tower

## Methodology

Simulate two MC experiments- one with regular water tower (wtp case) and the second with the water tower being completely made of air (i.e. the control case).

These experiments were conducted using the same random seed, inputs and code files so that the comparison is fair, and we can find the cause of the “missing” muons. I changed the current G4 code so that I can track and create text files of the important variables under consideration:

The data format is as shown below:

!(data)[DATAFORMAT.png]

Screen capture of the text data file output:

!(txt)[pic1.png]

## Information Flow

Raw data from G4 -> process.py —> processed data —> analyze.py —> result!

## Pseudocode

~~~~~~~~~~~~~~~~~~~~~
>Read the two text files.
>Import the data and create relevant data structures using an "indexed muon event” system
>The generated structure would be created using the following filters for each of the data file:
> events that generated mu-
>    out of these events sort events into two categories: detected (4/4) mu- and not detected muons
>    Using this new information the following algorithm is implemented to determine the scattering/absorption effects:
> count number of detected (4/4) mu- events in control case and wtp case
> count number of undetected mu- events in control case and wtp case
> initialize scattered mu counting variable
> initialize absorbed mu counting variable
> Identify detected (4/4) mu- events in control case that are not present in the wtp case and put them into an array
> Count these “missing” events
> Loop over the array:
    > look at location of last position of mu-
    > if last position in tower/tube/water:
            > absorbed mu counting variable += 1
    > else:
            > scattered mu counting variable += 1
> return ( scattered mu counting variable / “missing” events) and ( absorbed mu counting variable / “missing” events)
~~~~~~~~~~~~~~~~~~~~~

Using a similar process the number of muons scattered into the detector due to the presence of the wtp can also be found that. Thus, providing a thorough analysis.

## Results

Muons are scattered 87 percent of the time and absorbed 13 percent of the time.
