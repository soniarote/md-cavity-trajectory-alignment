"""
Configuration file for MD Cavity Trajectory Alignment

This module provides configuration examples for different types of
molecular systems and alignment scenarios.
"""

# =============================================================================
# ALIGNMENT SELECTIONS
# =============================================================================

# Metal-centered alignment (for metalloenzymes)
METAL_ALIGNMENT = "resname CU"                    # Copper atoms
ZINC_ALIGNMENT = "resname ZN"                     # Zinc atoms  
GENERAL_METAL = "name CU or name ZN or name MG"   # Multiple metals

# Backbone alignment (standard proteins)
BACKBONE_ALIGNMENT = "protein and name CA"        # Alpha carbons
FULL_BACKBONE = "protein and backbone"            # All backbone atoms
HEAVY_BACKBONE = "protein and backbone and not name H*"  # Non-hydrogen backbone

# Region-specific alignment
ACTIVE_SITE = "resid 50:150 and name CA"         # Specific region
BINDING_SITE = "around 5.0 resname CU and protein and name CA"  # Around metal
CONSERVED_CORE = "resid 20:50 80:120 and name CA"  # Multiple conserved regions

# Special cases
TRANSMEMBRANE = "resid 100:200 and name CA"      # Stable TM region
BETA_SHEET = "resid 30:60 120:150 and name CA"   # Beta-sheet regions
ALPHA_HELIX = "resid 10:30 70:90 and name CA"    # Alpha-helical regions

# =============================================================================
# SYSTEM CONFIGURATIONS
# =============================================================================

# Your MCOA system configuration
MCOA_CONFIG = {
    'alignment_selection': 'resname CU',
    'output_format': 'xtc',
    'output_directory': 'mcoa_aligned_cavities',
    'description': 'MCOA metalloenzyme with copper-centered alignment'
}

# Standard protein system
STANDARD_PROTEIN_CONFIG = {
    'alignment_selection': 'protein and name CA',
    'output_format': 'xtc', 
    'output_directory': 'protein_aligned_cavities',
    'description': 'Standard protein with backbone alignment'
}

# Membrane protein system
MEMBRANE_PROTEIN_CONFIG = {
    'alignment_selection': 'resid 50:150 and name CA',  # Stable TM region
    'output_format': 'dcd',
    'output_directory': 'membrane_aligned_cavities',
    'description': 'Membrane protein aligned on transmembrane region'
}

# Multi-domain protein
MULTIDOMAIN_CONFIG = {
    'alignment_selection': 'resid 20:80 150:210 and name CA',  # Conserved domains
    'output_format': 'xtc',
    'output_directory': 'multidomain_aligned_cavities', 
    'description': 'Multi-domain protein aligned on stable domains'
}

# =============================================================================
# FILE PATH TEMPLATES
# =============================================================================

# Template for your specific files
YOUR_FILES_TEMPLATE = {
    'wt_topology': "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS_extended.pdb",
    'wt_trajectory': "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS.xtc",
    'mutant_topology': "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER_extended.pdb",
    'mutant_trajectory': "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER.xtc"
}

# Generic template
GENERIC_TEMPLATE = {
    'wt_topology': "/path/to/wt_topology.pdb",
    'wt_trajectory': "/path/to/wt_trajectory.xtc",
    'mutant_topology': "/path/to/mutant_topology.pdb", 
    'mutant_trajectory': "/path/to/mutant_trajectory.xtc"
}

# =============================================================================
# OUTPUT FORMAT CONFIGURATIONS
# =============================================================================

# XTC format (recommended - compressed and widely supported)
XTC_CONFIG = {
    'format': 'xtc',
    'advantages': ['Small file size', 'Fast I/O', 'VMD/ChimeraX compatible'],
    'use_case': 'General purpose, long trajectories'
}

# DCD format (CHARMM/NAMD standard)
DCD_CONFIG = {
    'format': 'dcd',
    'advantages': ['CHARMM/NAMD native', 'Good precision'],
    'use_case': 'CHARMM/NAMD workflows'
}

# PDB format (simple but large files)
PDB_CONFIG = {
    'format': 'pdb',
    'advantages': ['Human readable', 'Universal compatibility'],
    'use_case': 'Small trajectories, debugging'
}

# =============================================================================
# ANALYSIS CONFIGURATIONS
# =============================================================================

# Quick analysis (few frames)
QUICK_ANALYSIS = {
    'frames_to_analyze': [0, 50, 100],
    'alignment_selection': 'resname CU',
    'output_format': 'xtc',
    'description': 'Quick check of key frames'
}

# Comprehensive analysis (many frames)
COMPREHENSIVE_ANALYSIS = {
    'frames_to_analyze': list(range(0, 1000, 25)),  # Every 25th frame
    'alignment_selection': 'protein and name CA',
    'output_format': 'xtc',
    'description': 'Detailed trajectory analysis'
}

# Publication analysis (high quality)
PUBLICATION_ANALYSIS = {
    'frames_to_analyze': 'all',
    'alignment_selection': 'resname CU',
    'output_format': 'xtc',
    'quality_checks': True,
    'description': 'High-quality analysis for publication'
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_alignment_selection(system_type):
    """
    Get recommended alignment selection for different system types.
    
    Parameters
    ----------
    system_type : str
        Type of molecular system
        
    Returns
    -------
    str
        MDAnalysis selection string
    """
    selections = {
        'metalloenzyme': 'resname CU',
        'standard_protein': 'protein and name CA',
        'membrane_protein': 'resid 50:150 and name CA',
        'multidomain': 'resid 20:80 150:210 and name CA',
        'enzyme': 'around 5.0 resname CU and protein and name CA',
        'antibody': 'protein and name CA and not resid 20:30 50:60',  # Avoid loops
        'nucleic_acid': 'nucleic and name P',  # Phosphate backbone
        'rna': 'resname A U G C and name P',
        'dna': 'resname DA DT DG DC and name P'
    }
    
    return selections.get(system_type, 'protein and name CA')

def validate_file_paths(file_dict):
    """
    Validate that all required files exist.
    
    Parameters
    ----------
    file_dict : dict
        Dictionary with file paths
        
    Returns
    -------
    list
        List of missing files
    """
    import os
    
    missing_files = []
    required_keys = ['wt_topology', 'wt_trajectory', 'mutant_topology', 'mutant_trajectory']
    
    for key in required_keys:
        if key not in file_dict:
            missing_files.append(f"Missing key: {key}")
        elif not os.path.exists(file_dict[key]):
            missing_files.append(f"File not found: {file_dict[key]}")
    
    return missing_files

def create_analysis_script(file_paths, config):
    """
    Generate a complete analysis script with your specific configuration.
    
    Parameters
    ----------
    file_paths : dict
        Dictionary with file paths
    config : dict
        Analysis configuration
        
    Returns
    -------
    str
        Complete Python script as string
    """
    
    script_template = '''#!/usr/bin/env python3
"""
Generated analysis script for MD cavity trajectory alignment
"""

from md_cavity_alignment import CavityAnalyzer

def main():
    """Run cavity alignment analysis with custom configuration."""
    
    print("🚀 MD CAVITY TRAJECTORY ALIGNMENT")
    print("Configuration: {description}")
    print("="*50)
    
    # File paths
    file_paths = {file_paths}
    
    try:
        # Initialize analyzer
        analyzer = CavityAnalyzer(
            wt_topology=file_paths['wt_topology'],
            wt_trajectory=file_paths['wt_trajectory'],
            mutant_topology=file_paths['mutant_topology'],
            mutant_trajectory=file_paths['mutant_trajectory']
        )
        
        # Perform alignment
        analyzer.align_trajectories(selection="{alignment_selection}")
        
        # Analyze specific frames
        analyzer.get_cavity_info(frame_idx=0)
        mid_frame = len(analyzer.u_wt.trajectory) // 2
        analyzer.get_cavity_info(frame_idx=mid_frame)
        
        # Save results
        analyzer.save_aligned_trajectories(
            output_dir="{output_directory}",
            format="{output_format}"
        )
        
        print("\\n🎉 Analysis completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {{e}}")

if __name__ == "__main__":
    main()
'''
    
    return script_template.format(
        file_paths=repr(file_paths),
        description=config.get('description', 'Custom analysis'),
        alignment_selection=config.get('alignment_selection', 'protein and name CA'),
        output_directory=config.get('output_directory', 'aligned_cavities'),
        output_format=config.get('output_format', 'xtc')
    )

def print_system_recommendations():
    """Print alignment recommendations for different system types."""
    
    recommendations = {
        'Metalloenzymes': {
            'selection': 'resname CU (or ZN, MG, etc.)',
            'reason': 'Metal centers provide stable reference points',
            'example': 'Copper oxidases, zinc proteases'
        },
        'Standard Proteins': {
            'selection': 'protein and name CA',
            'reason': 'Backbone alignment captures overall fold',
            'example': 'Globular enzymes, antibodies'
        },
        'Membrane Proteins': {
            'selection': 'resid X:Y and name CA (transmembrane)',
            'reason': 'Transmembrane regions are most stable',
            'example': 'GPCRs, ion channels'
        },
        'Multi-domain Proteins': {
            'selection': 'resid X:Y Z:W and name CA (stable domains)',
            'reason': 'Avoid flexible linkers between domains',
            'example': 'Kinases, large enzymes'
        },
        'Intrinsically Disordered': {
            'selection': 'resid X:Y and name CA (structured region)',
            'reason': 'Focus on any structured regions',
            'example': 'Transcription factors with structured domains'
        }
    }
    
    print("ALIGNMENT RECOMMENDATIONS BY SYSTEM TYPE")
    print("="*50)
    
    for system_type, info in recommendations.items():
        print(f"\n{system_type}:")
        print(f"  Selection: {info['selection']}")  
        print(f"  Reason: {info['reason']}")
        print(f"  Example: {info['example']}")

# =============================================================================
# VISUALIZATION COMMAND TEMPLATES
# =============================================================================

VMD_COMMANDS = '''
# VMD commands for cavity visualization
mol new wt_topology.pdb type pdb
mol addfile wt_aligned.xtc type xtc waitfor all

mol new mutant_topology.pdb type pdb  
mol addfile mutant_aligned.xtc type xtc waitfor all

# Style cavities as spheres
mol selection "resname CAV"
mol representation VDW 1.4 8.0
mol color ColorID 0
mol addrep 0

mol selection "resname CAV"
mol representation VDW 1.4 8.0  
mol color ColorID 1
mol addrep 1

# Animate trajectory
animate forward
'''

CHIMERAX_COMMANDS = '''
# ChimeraX commands for cavity visualization
open wt_topology.pdb
open wt_aligned.xtc coordset #1

open mutant_topology.pdb
open mutant_aligned.xtc coordset #2

# Style cavities
select :CAV
style sel sphere
color sel blue

select #2:CAV
style sel sphere
color sel red

# Animation controls
coordset slider
'''

PYMOL_COMMANDS = '''
# PyMOL commands for cavity visualization
load wt_topology.pdb, wt_system
load_traj wt_aligned.xtc, wt_system

load mutant_topology.pdb, mut_system
load_traj mutant_aligned.xtc, mut_system

# Style cavities
select wt_cavities, wt_system and resn CAV
select mut_cavities, mut_system and resn CAV

show spheres, wt_cavities
show spheres, mut_cavities

color blue, wt_cavities
color red, mut_cavities

# Animation
mplay
'''

def save_visualization_commands(output_dir):
    """Save visualization command files."""
    import os
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save VMD commands
    with open(os.path.join(output_dir, 'load_in_vmd.tcl'), 'w') as f:
        f.write(VMD_COMMANDS)
    
    # Save ChimeraX commands
    with open(os.path.join(output_dir, 'load_in_chimerax.cxc'), 'w') as f:
        f.write(CHIMERAX_COMMANDS)
    
    # Save PyMOL commands  
    with open(os.path.join(output_dir, 'load_in_pymol.pml'), 'w') as f:
        f.write(PYMOL_COMMANDS)
    
    print(f"Visualization command files saved to {output_dir}/")

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    # Example 1: Check file paths
    file_paths = YOUR_FILES_TEMPLATE
    missing = validate_file_paths(file_paths)
    
    if missing:
        print("Missing files:")
        for item in missing:
            print(f"  - {item}")
    else:
        print("All files found!")
    
    # Example 2: Get alignment recommendation
    system_type = 'metalloenzyme'
    selection = get_alignment_selection(system_type)
    print(f"Recommended alignment for {system_type}: {selection}")
    
    # Example 3: Print all recommendations
    print_system_recommendations()
    
    # Example 4: Generate custom analysis script
    custom_config = MCOA_CONFIG
    script_content = create_analysis_script(file_paths, custom_config)
    
    # Save the script
    with open('custom_cavity_analysis.py', 'w') as f:
        f.write(script_content)
    print("Custom analysis script saved as 'custom_cavity_analysis.py'")
    
    # Example 5: Save visualization commands
    save_visualization_commands('visualization_commands')