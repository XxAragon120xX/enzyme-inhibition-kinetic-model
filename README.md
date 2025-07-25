# Interactive Kinetic Model of Enzyme Inhibition

This project implements an interactive graphical application in Python to simulate the effect of inhibitors on bacterial growth, specifically on *Escherichia coli*, using a competitive inhibition model of the DHFR (dihydrofolate reductase) enzyme.

## Description

The tool allows you to:

- Adjust kinetic and pharmacokinetic parameters such as Ki, ka, ke, and initial inhibitor concentration.
- Visualize in real time:
  - Inhibition curves (μ vs. concentration).
  - Pharmacokinetic profile of the inhibitor.
  - Specific growth rate.
  - Bacterial biomass over time.
  - Comparison between two compounds with different activities.
  - Compare experimental data with model predictions.
- Export simulated data to `.csv` files.

## Requirements

This project was developed in Python 3 and uses the following libraries:

- numpy
- pandas
- matplotlib
- seaborn
- tkinter (included by default in Python installations)
- ttk (part of tkinter)

Install the dependencies with:

```bash
pip install numpy pandas matplotlib seaborn
```

## Usage

To start the application, run:

```bash
python Trabajo2.py
```

A window with the interactive model interface will appear.

## Included Models

- **Competitive inhibition:**  
  μ = μ_max / (1 + [I] / Ki)
- **Pharmacokinetics (elimination):**  
  C(t) = C₀·e^(-(ka+ke)t)
- **Logistic population growth:**  
  N(t) = (K_max · N₀) / [N₀ + (K_max - N₀)·e^(-μ·t)]

## Adjustable Parameters

- μ_max: Maximum growth rate without inhibitor.
- Ki_comp1 / Ki_comp2: Inhibition constants for two compounds.
- ka: Inhibitor absorption constant.
- ke: Inhibitor elimination constant.
- Initial inhibitor concentration.
- Time points for simulation.

## Calculated Results

- Estimated ED50 for both compounds.
- Relative potency between the compounds.
- Plots of concentration, growth rate, and biomass.

## Data Export

The generated kinetic data can be exported as `.csv` files for external analysis.

## Target Audience

This software is intended for:

- Students of biochemistry, pharmacology, and biotechnology.
- Instructors who wish to interactively demonstrate inhibition models.
- Researchers needing to validate or compare experimental data with a simple theoretical model.
