#!/usr/bin/env python3
"""
SIMPLE RUN SCRIPT - MD Cavity Trajectory Alignment

The easiest way to run the analysis:
1. Update the file paths below (lines 20-27)
2. Run: python run_analysis.py

That's it!
"""

from md_cavity_alignment import CavityAnalyzer

def main():
    print("MD Cavity Trajectory Alignment - Simple Run")
    print("=" * 50)
    
    # =============================================================
    # EDIT THESE PATHS TO YOUR FILES:
    # =============================================================
    
    # Your file paths here:
    wt_topology = "/path/to/your/wt_topology.pdb"
    wt_trajectory = "/path/to/your/wt_trajectory.xtc"
    mutant_topology = "/path/to/your/mutant_topology.pdb"
    mutant_trajectory = "/path/to/your/mutant_trajectory.xtc"
    
    # Alignment selection (choose one):
    alignment = "resname CU"                    # For metalloenzymes
    # alignment = "protein and name CA"         # For standard proteins
    
    # =============================================================
    # RUN THE ANALYSIS:
    # =============================================================
    
    try:
        # Create analyzer
        analyzer = CavityAnalyzer(wt_topology, wt_trajectory, 
                                mutant_topology, mutant_trajectory)
        
        # Align trajectories
        analyzer.align_trajectories(selection=alignment)
        
        # Check a few frames
        analyzer.get_cavity_info(frame_idx=0)
        
        # Save results
        analyzer.save_aligned_trajectories(format="xtc")
        
        print("\nSUCCESS! Check 'aligned_trajectories' folder for results.")
        
    except Exception as e:
        print(f"ERROR: {e}")
        print("Make sure to update the file paths at the top of this script!")

if __name__ == "__main__":
    main()
