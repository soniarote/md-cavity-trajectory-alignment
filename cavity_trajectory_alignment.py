#!/usr/bin/env python3
"""
MD Cavity Trajectory Alignment Tool

This module provides tools for aligning molecular dynamics trajectories 
containing cavity representations (CAV residues) for comparative analysis
between Wild-Type and Mutant protein systems.

Features:
- Robust trajectory alignment using customizable reference atoms
- CAV residue validation and tracking
- Multiple output formats for different visualization software
- Optimized for cavity dynamics analysis

Author: [Your Name]
Date: [Date]
Version: 1.0
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
import os
import warnings
from typing import Tuple, Optional

warnings.filterwarnings('ignore')


class CavityAnalyzer:
    """
    A tool for aligning MD trajectories with cavity representations for 
    comparative analysis between protein variants.
    
    This analyzer is specifically designed for trajectories containing
    CAV residues (cavity representations) and provides aligned outputs
    suitable for visualization and analysis.
    """
    
    def __init__(self, wt_topology: str, wt_trajectory: str, 
                 mutant_topology: str, mutant_trajectory: str):
        """
        Initialize the cavity analyzer with WT and mutant systems.
        
        Parameters
        ----------
        wt_topology : str
            Path to Wild-Type topology file (PDB, PSF, etc.)
        wt_trajectory : str
            Path to Wild-Type trajectory file (XTC, DCD, etc.)
        mutant_topology : str
            Path to Mutant topology file
        mutant_trajectory : str
            Path to Mutant trajectory file
            
        Note
        ----
        Both systems should contain CAV residues representing cavity spaces.
        Trajectories are typically water-free for optimal performance.
        """
        self.wt_topology = wt_topology
        self.wt_trajectory = wt_trajectory
        self.mutant_topology = mutant_topology
        self.mutant_trajectory = mutant_trajectory
        
        # Load molecular systems
        self.load_universes()
        
    def load_universes(self) -> None:
        """Load MDAnalysis Universe objects for both systems."""
        print("Loading molecular systems...")
        
        try:
            # Load Wild-Type system
            self.u_wt = mda.Universe(self.wt_topology, self.wt_trajectory)
            print(f"WT system loaded: {len(self.u_wt.trajectory)} frames")
            
            # Load Mutant system
            self.u_mutant = mda.Universe(self.mutant_topology, self.mutant_trajectory)
            print(f"Mutant system loaded: {len(self.u_mutant.trajectory)} frames")
            
            # Validate cavity content
            self.check_cavities()
            
        except Exception as e:
            raise RuntimeError(f"Error loading trajectories: {e}")
            
    def check_cavities(self) -> None:
        """Validate that both systems contain CAV residues."""
        # Check for CAV residues in both systems
        wt_cavities = self.u_wt.select_atoms("resname CAV")
        mutant_cavities = self.u_mutant.select_atoms("resname CAV")
        
        print(f"Cavity validation:")
        print(f"  WT cavities: {len(wt_cavities)} atoms")
        print(f"  Mutant cavities: {len(mutant_cavities)} atoms")
        
        if len(wt_cavities) == 0 or len(mutant_cavities) == 0:
            print("⚠️  Warning: No cavity spheres (CAV residues) found in one or both systems")
            print("    Ensure cavity detection has been performed and CAV residues are present")
            
    def align_trajectories(self, selection: str = "resname CU") -> None:
        """
        Align both trajectories using specified reference atoms.
        
        This method aligns all trajectory frames to maintain consistent
        spatial relationships for cavity comparison analysis.
        
        Parameters
        ----------
        selection : str, default="resname CU"
            MDAnalysis selection string for reference atoms.
            Common options:
            - "resname CU": Metal center alignment (default)
            - "protein and name CA": Backbone alignment
            - "resid 50:150": Specific residue range
            - "protein and backbone": Full backbone alignment
        """
        print(f"Aligning trajectories using selection: '{selection}'")
        
        # Select reference atoms for alignment
        ref_atoms_wt = self.u_wt.select_atoms(selection)
        ref_atoms_mutant = self.u_mutant.select_atoms(selection)
        
        print(f"  WT reference atoms: {len(ref_atoms_wt)}")
        print(f"  Mutant reference atoms: {len(ref_atoms_mutant)}")
        
        # Validate reference atoms
        if len(ref_atoms_wt) == 0 or len(ref_atoms_mutant) == 0:
            raise ValueError(f"No atoms found for selection: '{selection}'")
            
        if len(ref_atoms_wt) != len(ref_atoms_mutant):
            print("⚠️  Warning: Different number of reference atoms between systems")
            print("    This may indicate sequence differences or structural variations")
        
        try:
            # Move to first frame and store reference coordinates
            self.u_wt.trajectory[0]
            ref_coordinates = ref_atoms_wt.positions.copy()
            
            # Align WT trajectory to its first frame
            print("Aligning WT system...")
            align.AlignTraj(self.u_wt, self.u_wt, select=selection, 
                           in_memory=True, match_atoms=True).run()
            
            # Align mutant trajectory to WT first frame  
            print("Aligning mutant system...")
            align.AlignTraj(self.u_mutant, self.u_wt, select=selection,
                           in_memory=True, match_atoms=True).run()
            
            print("✅ Trajectory alignment completed successfully")
            
        except Exception as e:
            raise RuntimeError(f"Alignment failed: {e}")
        
    def get_cavity_info(self, frame_idx: int = 0) -> Tuple[mda.AtomGroup, mda.AtomGroup]:
        """
        Analyze cavity information for a specific trajectory frame.
        
        Parameters
        ----------
        frame_idx : int, default=0
            Frame index to analyze (0-based indexing)
            
        Returns
        -------
        tuple of AtomGroup
            (wt_cavities, mutant_cavities) - Cavity atoms for both systems
        """
        # Move to specified frame
        self.u_wt.trajectory[frame_idx]
        self.u_mutant.trajectory[frame_idx]
        
        # Select cavity atoms
        wt_cavities = self.u_wt.select_atoms("resname CAV")
        mutant_cavities = self.u_mutant.select_atoms("resname CAV")
        
        # Display frame analysis
        print(f"\nFrame {frame_idx} Analysis:")
        print(f"  WT cavity atoms: {len(wt_cavities)}")
        print(f"  Mutant cavity atoms: {len(mutant_cavities)}")
        
        if len(wt_cavities) > 0:
            x_range = wt_cavities.positions[:,0]
            print(f"  WT X-coordinate range: [{x_range.min():.2f}, {x_range.max():.2f}]")
            
        if len(mutant_cavities) > 0:
            x_range = mutant_cavities.positions[:,0]  
            print(f"  Mutant X-coordinate range: [{x_range.min():.2f}, {x_range.max():.2f}]")
            
        return wt_cavities, mutant_cavities
    
    def save_aligned_trajectories(self, output_dir: str = "aligned_trajectories", 
                                 format: str = "xtc") -> None:
        """
        Save aligned trajectories for visualization and analysis.
        
        This generates visualization-ready files that can be loaded in
        VMD, ChimeraX, PyMOL, or other molecular visualization software.
        
        Parameters
        ----------
        output_dir : str, default="aligned_trajectories"
            Output directory for aligned files
        format : str, default="xtc"
            Trajectory format ('xtc', 'dcd', 'pdb')
            
        Note
        ----
        XTC format is recommended for its compression and broad compatibility.
        All output files maintain the original topology information.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        print(f"Saving aligned trajectories in {format.upper()} format to {output_dir}/")
        
        # Select all atoms for output
        wt_all = self.u_wt.select_atoms("all")
        mutant_all = self.u_mutant.select_atoms("all")
        
        # Save trajectories in specified format
        self._write_trajectory_file(wt_all, f"{output_dir}/wt_aligned.{format}", format)
        self._write_trajectory_file(mutant_all, f"{output_dir}/mutant_aligned.{format}", format)
        
        # Save reference topology files (first frame as PDB)
        self.u_wt.trajectory[0]
        wt_all.write(f"{output_dir}/wt_topology.pdb")
        
        self.u_mutant.trajectory[0]
        mutant_all.write(f"{output_dir}/mutant_topology.pdb")
                
        print("✅ Aligned trajectories saved successfully")
        print(f"   📄 Trajectory files: wt_aligned.{format}, mutant_aligned.{format}")
        print(f"   📄 Topology files: wt_topology.pdb, mutant_topology.pdb")
        print(f"\n💡 Visualization tip:")
        print(f"   Load topology + trajectory in VMD/ChimeraX for cavity animation")
        
    def _write_trajectory_file(self, atoms: mda.AtomGroup, filename: str, format: str) -> None:
        """Write trajectory file in specified format."""
        if format.lower() == "xtc":
            with mda.Writer(filename, atoms.n_atoms) as writer:
                for ts in atoms.universe.trajectory:
                    writer.write(atoms)
        elif format.lower() == "dcd":
            with mda.Writer(filename, atoms.n_atoms) as writer:
                for ts in atoms.universe.trajectory:
                    writer.write(atoms)
        elif format.lower() == "pdb":
            # Multi-frame PDB
            with mda.Writer(filename, multiframe=True) as writer:
                for ts in atoms.universe.trajectory:
                    writer.write(atoms)
        else:
            raise ValueError(f"Unsupported format: {format}")


def main():
    """
    Example usage of the cavity analyzer.
    
    Update the file paths below to match your specific system.
    """
    print("🚀 MD CAVITY TRAJECTORY ALIGNMENT")
    print("="*50)
    
    # Configuration - UPDATE THESE PATHS TO YOUR FILES
    file_paths = {
        'wt_topology': "/path/to/WT_topology_with_cavities.pdb",
        'wt_trajectory': "/path/to/WT_trajectory.xtc", 
        'mutant_topology': "/path/to/mutant_topology_with_cavities.pdb",
        'mutant_trajectory': "/path/to/mutant_trajectory.xtc"
    }
    
    # Example with your specific files (commented out - uncomment and modify)
    # file_paths = {
    #     'wt_topology': "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS_extended.pdb",
    #     'wt_trajectory': "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS.xtc",
    #     'mutant_topology': "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER_extended.pdb",
    #     'mutant_trajectory': "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER.xtc"
    # }
    
    try:
        # Initialize analyzer
        analyzer = CavityAnalyzer(
            wt_topology=file_paths['wt_topology'],
            wt_trajectory=file_paths['wt_trajectory'],
            mutant_topology=file_paths['mutant_topology'],
            mutant_trajectory=file_paths['mutant_trajectory']
        )
        
        # Perform trajectory alignment
        # Options: "resname CU", "protein and name CA", "resid 50:150"
        analyzer.align_trajectories(selection="resname CU")
        
        # Analyze cavity information at key frames
        analyzer.get_cavity_info(frame_idx=0)  # First frame
        
        # Middle frame
        mid_frame = len(analyzer.u_wt.trajectory) // 2
        analyzer.get_cavity_info(frame_idx=mid_frame)
        
        # Save aligned trajectories for visualization
        analyzer.save_aligned_trajectories(format="xtc")
        
        print("\n🎉 Analysis completed successfully!")
        print("="*50)
        print("📋 Next steps:")
        print("   1. Load aligned files in VMD/ChimeraX for visualization")
        print("   2. Animate trajectory to observe cavity dynamics") 
        print("   3. Compare WT vs mutant cavity evolution")
        print("\n💡 Visualization commands saved in output directory")
        
    except FileNotFoundError as e:
        print(f"❌ File not found: {e}")
        print("Please update the file paths in the script to match your files")
    except Exception as e:
        print(f"❌ Analysis error: {e}")
        print("\n🔧 Troubleshooting:")
        print("   - Verify file paths are correct")
        print("   - Ensure CAV residues are present in topology files")
        print("   - Check that MDAnalysis is properly installed")
        print("   - Verify trajectory and topology file compatibility")


if __name__ == "__main__":
    main()
