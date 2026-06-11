[![DOI](https://zenodo.org/badge/1242719459.svg)](https://doi.org/10.5281/zenodo.20637315)

# Tropical Dry Forest Tree Growth Forecasting (TDF-Growth)

Welcome to the **TDF-Growth* repository, an integral part of the research at the [Forest Ecophys Lab](https://forestecophys.com/). This repository is the central hub for our team, encompassing our project overview, codebase, field protocols and more...

## Our Project
The goal of our project is to produce seasonal forecast of tree growth for an array of tropical dry forest tress.

## Documentation
- Access detailed documentation on our [GitHub Pages site](https://your-gh-pages-url/).

## Author
- German Vargas G., Oregon State University

## Code Repository Structure
- **Data Processing**: Scripts for cleaning, merging, and managing datasets.
- **Analysis Code**: Scripts for data analysis, statistical modeling, etc.
- **Python Code**: Jupyter workbooks for downloading climate re-analysis data from ERA5.

```
ACG-GROWTH-FORECASTING

├── code
│   └── analysis
│   ├── data-processing
│   ├── functions
│   └── python
│       └── Pynotebooks
│           └── downloadERA5.ipynb
├── ACG-Growth-Forecasting.Rproj
├── requirements.txt
└── README.md
```

## Contributing to This Repository
- Please adhere to these guidelines:
  - Ensure commits have clear and concise messages.
  - Document major changes as annotations in the code.
  - Review and merge changes through pull requests for oversight.

## Setup Instructions

1. **Clone the Repository**: 
   ```bash
   git clone <repository-url>
   cd ACG-Growth-Forecasting
   ```

2. **Install Dependencies**: 
    Recommended quick start (venv + pip):
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade pip
    pip install -r requirements.txt
    ```
3. **Prepare Access to Data**: 
   Make sure you have cloud storage as part of your file system