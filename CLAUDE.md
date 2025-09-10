# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains plugin configurations for CellCraft, a web-based visual programming application for gene regulatory network (GRN) inference. It manages multiple GRN reconstruction tools integrated through a modular plugin system.

## Key Commands

### Plugin Management
- **Generate plugins.csv**: `python generate_plugins_csv.py` - Combines plugin_initialization.csv with metadata.json and dependency files from each plugin folder
- **Test plugins.csv**: `python test_plugins_csv.py` - Validates JSON parsing and data structure compatibility before database initialization
- **Initialize database**: `python init_db.py` - Loads plugins.csv 
data into the database

### Docker Operations
Each plugin has its own Dockerfile for containerized execution:
- Build: `docker build -t cellcraft-[plugin-name] ./[PLUGIN-NAME]/`
- Run: Docker containers use Snakemake workflows with micromamba environments

## Architecture & Structure

### Plugin System
Each plugin folder (e.g., `TENET/`, `FastSCODE/`, `GENIE3/`) contains:
- **metadata.json**: Defines workflow nodes (drawflow) and processing rules
  - `drawflow`: Visual workflow configuration with node connections
  - `rules`: Processing steps with input/output specifications and parameters
- **scripts/**: Python scripts for data processing
  - `tenet_input.py`: Common preprocessing script for h5ad input files
  - `[plugin]_for_cellcraft.py`: Main plugin execution script
- **dependency/**: Requirements files for the plugin environment
- **Dockerfile**: Container definition using micromamba for dependency management
- **Snakefile**: Workflow orchestration (if applicable)

### Data Flow
1. **Input**: Single-cell RNA-seq data in h5ad format
2. **Processing**: Plugin-specific GRN inference algorithms
3. **Output**: Network files (.sif format) and analysis results

### Key Scripts
- **generate_plugins_csv.py**: 
  - Reads plugin_initialization.csv for basic info
  - Loads metadata.json from each plugin folder
  - Applies transformations (FastScode → FastSCODE, plugin → cellcraft-plugin)
  - Generates unified plugins.csv with escaped JSON fields
  
- **test_plugins_csv.py**: 
  - Validates CSV structure and required columns
  - Tests JSON parsing for dependencies, drawflow, and rules
  - Simulates database initialization process
  - Verifies string transformations were applied correctly

## Plugin Development

### Adding a New Plugin
1. Create plugin folder with standard structure
2. Define metadata.json with drawflow and rules
3. Add dependency requirements
4. Create Dockerfile following the template pattern
5. Add entry to plugin_initialization.csv
6. Run `python generate_plugins_csv.py` to update plugins.csv
7. Test with `python test_plugins_csv.py`

### Plugin Metadata Structure
```json
{
  "name": "PluginName",
  "author": "author",
  "description": "Plugin description",
  "drawflow": { /* Visual workflow configuration */ },
  "rules": { /* Processing rules and parameters */ }
}
```

### Common Parameter Types
- `inputFile`: Required input file
- `optionalInputFile`: Optional input file  
- `outputFile`: Generated output file
- `h5adParameter`: Parameter from h5ad file metadata
- `int`, `float`: Numeric parameters with min/max constraints

## Available Plugins

- **TENET**: Transfer Entropy-based Network Reconstruction
- **FastTENET**: Accelerated TENET implementation
- **SCODE/FastSCODE**: ODE-based regulatory network inference
- **SCRIBE**: Causal network inference using restricted directed information
- **LEAP**: Pseudotime-based co-expression network construction
- **GRNBoost2**: Gradient boosting-based scalable GRN inference
- **GENIE3**: Tree-based network inference using Random Forests
- **GRNViz**: Visualization tools for network analysis

## Testing & Validation

Before deploying plugins:
1. Run `python test_plugins_csv.py` to validate data structures
2. Check JSON parsing for all metadata fields
3. Verify string transformations are applied correctly
4. Ensure all required columns are present
5. Test simulated database initialization

## Important Notes

- All plugin paths should use `./cellcraft-plugin/` prefix (not `./plugin/`)
- JSON fields in CSV must be properly escaped for database compatibility
- Each plugin runs in an isolated Docker container with micromamba environment
- Snakemake is used for workflow orchestration within containers