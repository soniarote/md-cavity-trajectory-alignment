#!/usr/bin/env python3
"""
Example analysis script for MD Cavity Trajectory Alignment

This script shows how to use the cavity analyzer with your specific files.
Simply update the file paths below and run this script.
"""

from md_cavity_alignment import CavityAnalyzer

def main():
    """
    Example analysis using the MD cavity trajectory alignment tool.
    
    STEP 1: Update the file paths below to point to your files
    STEP 2: Run this script with: python example_analysis.py
    """
    
    print("🚀 MD CAVITY TRAJECTORY ALIGNMENT - EXAMPLE")
    print("="*60)
    
    # =================================================================
    # STEP 1: UPDATE THESE PATHS TO YOUR FILES
    # =================================================================
    
    file_paths = {
        # Wild-Type system files
        'wt_topology': "/path/to/your/WT_topology_with_cavities.pdb",
        'wt_trajectory': "/path/to/your/WT_trajectory.xtc",
        
        # Mutant system files  
        'mutant_topology': "/path/to/your/mutant_topology_with_cavities.pdb",
        'mutant_trajectory': "/path/to/your/mutant_trajectory.xtc"
    }
    
    # Example for the MCOA system (uncomment and modify if needed):
    # file_paths = {
    #     'wt_topology': "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS_extended.pdb",
    #     'wt_trajectory': "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS.xtc",
    #     'mutant_topology': "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER_extended.pdb",
    #     'mutant_trajectory': "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER.xtc"
    # }
    
    # =================================================================
    # STEP 2: CONFIGURATION OPTIONS
    # =================================================================
    
    # Alignment options - choose the best for your system:
    alignment_selection = "resname CU"              # For metalloenzymes (default)
    # alignment_selection = "protein and name CA"   # For standard proteins
    # alignment_selection = "resid 50:150"          # For specific regions
    
    # Output options
    output_format = "xtc"                           # Recommended: compressed format
    output_directory = "aligned_cavity_results"    # Where to save results
    
    # Frames to analyze (optional - leave empty to skip detailed analysis)
    frames_to_check = [0, 50, 100]  # First, middle, and another frame
    
    # =================================================================
    # STEP 3: RUN THE ANALYSIS
    # =================================================================
    
    try:
        print(f"📁 Input files:")
        print(f"   WT topology: {file_paths['wt_topology']}")
        print(f"   WT trajectory: {file_paths['wt_trajectory']}")  
        print(f"   Mutant topology: {file_paths['mutant_topology']}")
        print(f"   Mutant trajectory: {file_paths['mutant_trajectory']}")
        print()
        
        # Initialize the analyzer
        print("🔧 Initializing cavity analyzer...")
        analyzer = CavityAnalyzer(
            wt_topology=file_paths['wt_topology'],
            wt_trajectory=file_paths['wt_trajectory'],
            mutant_topology=file_paths['mutant_topology'],
            mutant_trajectory=file_paths['mutant_trajectory']
        )
        print("✅ Analyzer initialized successfully!")
        print()
        
        # Perform trajectory alignment
        print(f"🔄 Aligning trajectories using: '{alignment_selection}'")
        analyzer.align_trajectories(selection=alignment_selection)
        print("✅ Trajectory alignment completed!")
        print()
        
        # Analyze cavity information at specific frames
        print("🔍 Analyzing cavity information...")
        for frame_idx in frames_to_check:
            try:
                analyzer.get_cavity_info(frame_idx=frame_idx)
            except:
                print(f"   ⚠️  Could not analyze frame {frame_idx} (may not exist)")
        print()
        
        # Save aligned trajectories for visualization
        print(f"💾 Saving results to: {output_directory}/")
        analyzer.save_aligned_trajectories(
            output_dir=output_directory,
            format=output_format
        )
        print("✅ Results saved successfully!")
        print()
        
        # =================================================================
        # STEP 4: NEXT STEPS
        # =================================================================
        
        print("🎉 ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*60)
        print("📋 NEXT STEPS:")
        print("   1. Check the output directory for aligned files")
        print("   2. Load files in molecular visualization software:")
        print("      • VMD: Load topology + trajectory")  
        print("      • ChimeraX: Open PDB then add trajectory")
        print("      • PyMOL: Load topology then load_traj")
        print("   3. Visualize CAV residues as spheres to see cavities")
        print("   4. Animate trajectory to observe cavity dynamics")
        print()
        print("🔬 VISUALIZATION TIP:")
        print("   Select 'resname CAV' and display as spheres")
        print("   Color WT and mutant cavities differently for comparison")
        print()
        
    except FileNotFoundError as e:
        print(f"❌ FILE ERROR: {e}")
        print("\n🔧 SOLUTION:")
        print("   Update the file paths at the top of this script")
        print("   Make sure all files exist and paths are correct")
        
    except Exception as e:
        print(f"❌ ANALYSIS ERROR: {e}")
        print("\n🔧 TROUBLESHOOTING:")
        print("   1. Check that CAV residues are present in topology files")
        print("   2. Verify trajectory and topology compatibility") 
        print("   3. Make sure MDAnalysis is installed: pip install MDAnalysis")
        print("   4. Try a different alignment selection if alignment fails")

def check_requirements():
    """Check if required packages are installed."""
    try:
        import MDAnalysis as mda
        import numpy as np
        print("✅ All required packages are installed")
        return True
    except ImportError as e:
        print(f"❌ Missing package: {e}")
        print("🔧 Install with: pip install MDAnalysis numpy")
        return False

if __name__ == "__main__":
    # Check requirements first
    if check_requirements():
        # Run the analysis
        main()
    else:
        print("Please install required packages first")
