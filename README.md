[![DOI](https://zenodo.org/badge/761899116.svg)](https://zenodo.org/doi/10.5281/zenodo.11167012)

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

## Meeting Notes and Agendas
- Regular updates to keep all group members informed and engaged with the project's progress and direction.

## Contributing to This Repository
- Please adhere to these guidelines:
  - Ensure commits have clear and concise messages.
  - Document major changes as annotations in the code.
  - Review and merge changes through pull requests for oversight.

## Customize Your Repository
- **Edit This Readme**: Update with information specific to your project.
- **Update Group Member Bios**: Add detailed information about each group member's expertise and role.
- **Organize Your Code**: Use logical structure and clear naming conventions.
- **Document Your Data**: Include a data directory with README files for datasets.
- **Outline Your Methods**: Create a METHODS.md file for methodologies and tools.
- **Set Up Project Management**: Use 'Issues' and 'Projects' for task tracking.
- **Add a License**: Include an appropriate open-source license.
- **Create Contribution Guidelines**: Establish a CONTRIBUTING.md file.
- **Review and Merge Workflow**: Document your process for reviewing and merging changes.
- **Establish Communication Channels**: Set up channels like Slack or Discord for discussions.

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