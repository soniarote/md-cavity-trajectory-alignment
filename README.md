# MD Cavity Trajectory Alignment: Dynamic Cavity Visualization

## Overview

This tool performs structural alignment of molecular dynamics trajectories containing cavity representations (CAV residues), enabling visualization and analysis of cavity evolution throughout simulations. It's specifically designed for comparative studies between Wild-Type and Mutant protein systems, providing aligned trajectories ready for visualization in VMD, ChimeraX, or PyMOL.

## Key Features

- **Robust trajectory alignment** using customizable reference atoms (metals, backbone, etc.)
- **CAV residue validation** ensures cavity representations are present
- **Multiple output formats** (XTC, DCD, PDB) for different visualization software
- **Cavity tracking** throughout MD simulation timescale
- **Water-free analysis** optimized for cavity-focused studies

## Scientific Applications

- **Cavity dynamics analysis**: Track how binding pockets change over time
- **Mutation impact studies**: Compare cavity accessibility between WT and mutant
- **Druggability assessment**: Visualize transient binding sites
- **Allosteric analysis**: Monitor cavity communication pathways
- **Enzyme mechanism studies**: Observe active site flexibility

## Requirements

```bash
pip install MDAnalysis numpy
```

## Input Requirements

### Trajectory Files
- **Topology files**: PDB files with protein structure
- **Trajectory files**: XTC, DCD, or TRR format molecular dynamics trajectories
- **CAV residues**: Pre-computed cavity representations as CAV residues

### Cavity Detection Prerequisites
Cavities must be pre-computed and represented as CAV residues using tools such as:
- **CASTp** (Computed Atlas of Surface Topography)
- **fpocket** (Fast protein cavity detection)
- **POVME** (Pocket Volume Measurer)
- **Custom cavity detection workflows**

## Quick Start

### Installation

1. **Clone the repository**:
```bash
git clone https://github.com/your_username/md-cavity-trajectory-alignment.git
cd md-cavity-trajectory-alignment
```

2. **Install dependencies**:
```bash
pip install -r requirements.txt
```

3. **Update file paths** in `md_cavity_alignment.py`:
```python
# Edit the main() function with your file paths
file_paths = {
    'wt_topology': "/path/to/your/WT_topology_with_cavities.pdb",
    'wt_trajectory': "/path/to/your/WT_trajectory.xtc",
    'mutant_topology': "/path/to/your/mutant_topology_with_cavities.pdb",
    'mutant_trajectory': "/path/to/your/mutant_trajectory.xtc"
}
```

4. **Run the analysis**:
```bash
python md_cavity_alignment.py
```

## Usage

### Method 1: Direct Execution (Easiest)

Simply update the file paths in the script and run:
```bash
python md_cavity_alignment.py
```

### Method 2: Import as Module

```python
from md_cavity_alignment import CavityAnalyzer

# Initialize analyzer
analyzer = CavityAnalyzer(
    wt_topology="wt_system_with_cavities.pdb",
    wt_trajectory="wt_trajectory.xtc",
    mutant_topology="mutant_system_with_cavities.pdb", 
    mutant_trajectory="mutant_trajectory.xtc"
)

# Align trajectories (default: copper atoms)
analyzer.align_trajectories()

# Analyze cavity information at specific frames
analyzer.get_cavity_info(frame_idx=0)
analyzer.get_cavity_info(frame_idx=100)

# Save aligned trajectories
analyzer.save_aligned_trajectories(format="xtc")
```

### Advanced Configuration

```python
# Align using different reference atoms
analyzer.align_trajectories(selection="protein and name CA")  # Backbone
analyzer.align_trajectories(selection="resname CU")           # Metal centers
analyzer.align_trajectories(selection="resid 50:150")        # Specific region

# Save in different formats
analyzer.save_aligned_trajectories(format="xtc")  # Compressed (recommended)
analyzer.save_aligned_trajectories(format="dcd")  # CHARMM/NAMD format
analyzer.save_aligned_trajectories(format="pdb")  # Multi-frame PDB
```

## Methodology

### Alignment Strategy
1. **Reference selection**: Choose stable atoms for alignment (metals, backbone, conserved regions)
2. **Structural superposition**: Align all frames to first WT frame using MDAnalysis
3. **Cross-system alignment**: Align mutant trajectory to WT reference frame
4. **Cavity preservation**: CAV residue positions maintained throughout alignment

### CAV Residue System
- **Sphere representation**: Each CAV atom represents a cavity sphere
- **Dynamic positioning**: CAV coordinates update with protein conformational changes
- **Radius information**: Sphere sizes typically 1.4-2.0 Å
- **Volume calculation**: Total cavity volume from sphere ensemble

## Output Files

### Aligned Trajectories
- `wt_aligned.xtc` / `mutant_aligned.xtc`: Aligned trajectory files
- `wt_topology.pdb` / `mutant_topology.pdb`: Reference structure files

### Visualization Ready
All output files are immediately compatible with:
- **VMD**: Load topology + trajectory for cavity animation
- **ChimeraX**: Import aligned structures for comparative analysis  
- **PyMOL**: Visualize cavity evolution over time

## Visualization Workflow

### VMD Visualization
```tcl
# Load WT system
mol new wt_topology.pdb type pdb
mol addfile wt_aligned.xtc type xtc waitfor all

# Load mutant system  
mol new mutant_topology.pdb type pdb
mol addfile mutant_aligned.xtc type xtc waitfor all

# Visualize cavities as spheres
mol selection "resname CAV"
mol representation VDW 1.4 12.0
mol color ColorID 0
mol addrep 0

# Animate to see cavity dynamics
animate forward
```

### ChimeraX Commands
```bash
# Load aligned systems
open wt_topology.pdb
open wt_aligned.xtc coordset #1

open mutant_topology.pdb  
open mutant_aligned.xtc coordset #2

# Style cavity spheres
select :CAV
style sel sphere
color sel blue

select #2:CAV
style sel sphere  
color sel red

# Start trajectory animation
coordset slider
```

## Key Advantages

### Performance Benefits
- **Water-free trajectories**: Smaller file sizes and faster processing
- **Efficient alignment**: MDAnalysis optimized algorithms
- **Memory management**: Handles large trajectories efficiently

### Analysis Benefits  
- **Temporal cavity mapping**: Track cavity opening/closing events
- **Quantitative comparison**: Measure cavity size changes over time
- **Mutation impact**: Direct visualization of structural effects
- **Dynamic druggability**: Assess binding site accessibility variation

## Example Results

The analysis typically reveals:
- **Cavity breathing motions** throughout the simulation
- **Mutant-specific cavity changes** compared to wild-type
- **Transient binding sites** that appear temporarily  
- **Allosteric coupling** between distant cavities

## Technical Details

### Alignment Precision
- Uses MDAnalysis AlignTraj for robust superposition
- Supports various atom selections for alignment reference
- Maintains coordinate consistency across systems
- Preserves cavity-protein spatial relationships

### Cavity Validation
- Automatically detects CAV residues in both systems
- Reports cavity atom counts for quality control
- Validates coordinate ranges after alignment
- Warns about potential alignment issues

## File Formats Supported

### Input Formats
- **Topology**: PDB, PSF, PRMTOP, GRO
- **Trajectory**: XTC, DCD, TRR, NetCDF

### Output Formats  
- **XTC**: Compressed, widely compatible (recommended)
- **DCD**: CHARMM/NAMD format
- **PDB**: Multi-frame PDB for simple visualization

## Best Practices

### Alignment Selection
- **Metal centers**: Use `resname CU` for metalloenzymes
- **Backbone**: Use `protein and name CA` for standard proteins
- **Conserved regions**: Use specific residue ranges for membrane proteins
- **Active sites**: Use binding site residues for focused analysis

### Quality Control
- Check alignment RMSD values in output
- Verify cavity atom counts remain consistent  
- Inspect first/last frames visually for alignment quality
- Compare cavity coordinate ranges before/after alignment

## Troubleshooting

### Common Issues
- **No CAV residues found**: Ensure cavity detection was performed correctly
- **Different atom counts**: Check for sequence differences between systems
- **High alignment RMSD**: Consider different alignment selection
- **Missing trajectory frames**: Verify trajectory file integrity

### Performance Tips
- Use XTC format for large trajectories
- Consider trajectory striding for very long simulations
- Remove unnecessary atoms before analysis
- Use solid-state drive for improved I/O performance

## Applications in Drug Discovery

### Binding Site Analysis
- **Pocket druggability**: Assess cavity size and shape changes
- **Selectivity studies**: Compare cavity accessibility between variants
- **Allosteric sites**: Identify cryptic binding pockets
- **PROTAC applications**: Analyze degrader binding site dynamics

### Structure-Activity Relationships  
- **Mutation effects**: Quantify structural impacts on binding sites
- **Flexibility analysis**: Characterize pocket conformational dynamics
- **Pathway mapping**: Track cavity communication networks
- **Pharmacophore evolution**: Monitor key interaction points

## Integration with Other Tools

### Upstream (Cavity Detection)
- Import CAV residues from CASTp, fpocket, or POVME
- Compatible with custom cavity detection workflows
- Supports various cavity representation formats

### Downstream (Analysis)
- Output ready for cavity volume calculations
- Compatible with binding site analysis software
- Integrates with molecular visualization pipelines

## Citation

If you use this tool in your research, please cite:
[Your publication details]

## License  

[Add your preferred license]

## Support

For questions or issues:
- Check troubleshooting section above
- Verify input file formats and CAV residue presence
- Ensure MDAnalysis installation is current
- Review alignment selection for your specific system