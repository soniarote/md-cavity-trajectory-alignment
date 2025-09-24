#!/usr/bin/env python3
"""
Análisis de Cavidades en Trayectorias MD
Parte 1: Alineamiento de trayectorias WT y mutante
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
import os
import warnings
warnings.filterwarnings('ignore')

class CavityAnalyzer:
    def __init__(self, wt_topology, wt_trajectory, mutant_topology, mutant_trajectory):
        """
        Inicializa el analizador de cavidades
        
        Parameters:
        -----------
        wt_topology : str
            Ruta al archivo de topología del sistema WT
        wt_trajectory : str
            Ruta a la trayectoria del sistema WT
        mutant_topology : str
            Ruta al archivo de topología del sistema mutante
        mutant_trajectory : str
            Ruta a la trayectoria del sistema mutante
        """
        self.wt_topology = wt_topology
        self.wt_trajectory = wt_trajectory
        self.mutant_topology = mutant_topology
        self.mutant_trajectory = mutant_trajectory
        
        # Cargar universos
        self.load_universes()
        
    def load_universes(self):
        """Carga los universos MDAnalysis para ambos sistemas"""
        print("Cargando universos...")
        
        # Cargar sistema WT
        self.u_wt = mda.Universe(self.wt_topology, self.wt_trajectory)
        print(f"Sistema WT cargado: {len(self.u_wt.trajectory)} frames")
        
        # Cargar sistema mutante
        self.u_mutant = mda.Universe(self.mutant_topology, self.mutant_trajectory)
        print(f"Sistema mutante cargado: {len(self.u_mutant.trajectory)} frames")
        
        # Verificar que ambos sistemas tengan cavidades
        self.check_cavities()
        
    def check_cavities(self):
        """Verifica la presencia de residuos CAV en ambos sistemas"""
        # Cavidades en WT
        wt_cavities = self.u_wt.select_atoms("resname CAV")
        print(f"Cavidades en WT: {len(wt_cavities)} átomos")
        
        # Cavidades en mutante
        mutant_cavities = self.u_mutant.select_atoms("resname CAV")
        print(f"Cavidades en mutante: {len(mutant_cavities)} átomos")
        
        if len(wt_cavities) == 0 or len(mutant_cavities) == 0:
            print("⚠️  Advertencia: No se encontraron cavidades (resname CAV) en uno o ambos sistemas")
            
    def align_trajectories(self, selection="resname CU"):
        """
        Alinea las trayectorias usando la proteína como referencia
        
        Parameters:
        -----------
        selection : str
            Selección de átomos para el alineamiento (por defecto: alpha carbons de proteína)
        """
        print(f"Alineando trayectorias usando selección: {selection}")
        
        # Seleccionar átomos de referencia
        ref_atoms_wt = self.u_wt.select_atoms(selection)
        ref_atoms_mutant = self.u_mutant.select_atoms(selection)
        
        print(f"Átomos de referencia WT: {len(ref_atoms_wt)}")
        print(f"Átomos de referencia mutante: {len(ref_atoms_mutant)}")
        
        # Verificar que ambos sistemas tengan el mismo número de átomos de referencia
        if len(ref_atoms_wt) != len(ref_atoms_mutant):
            print("⚠️  Advertencia: Número diferente de átomos de referencia")
            print("    Esto puede indicar diferencias en la secuencia o estructura")
        
        # Alinear usando el primer frame del WT como referencia
        self.u_wt.trajectory[0]  # Ir al primer frame
        ref_coordinates = ref_atoms_wt.positions.copy()
        
        # Alinear todas las trayectorias
        print("Alineando sistema WT...")
        align.AlignTraj(self.u_wt, self.u_wt, select=selection, 
                       in_memory=True, match_atoms=True).run()
        
        print("Alineando sistema mutante...")
        align.AlignTraj(self.u_mutant, self.u_wt, select=selection,
                       in_memory=True, match_atoms=True).run()
        
        print("✅ Alineamiento completado")
        
    def get_cavity_info(self, frame_idx=0):
        """
        Obtiene información sobre las cavidades en un frame específico
        
        Parameters:
        -----------
        frame_idx : int
            Índice del frame a analizar
        """
        # Ir al frame especificado
        self.u_wt.trajectory[frame_idx]
        self.u_mutant.trajectory[frame_idx]
        
        # Obtener cavidades
        wt_cavities = self.u_wt.select_atoms("resname CAV")
        mutant_cavities = self.u_mutant.select_atoms("resname CAV")
        
        print(f"\nFrame {frame_idx}:")
        print(f"Cavidades WT: {len(wt_cavities)} átomos")
        print(f"Cavidades mutante: {len(mutant_cavities)} átomos")
        
        if len(wt_cavities) > 0:
            print(f"Rango coordenadas WT: X[{wt_cavities.positions[:,0].min():.2f}, {wt_cavities.positions[:,0].max():.2f}]")
            
        if len(mutant_cavities) > 0:
            print(f"Rango coordenadas mutante: X[{mutant_cavities.positions[:,0].min():.2f}, {mutant_cavities.positions[:,0].max():.2f}]")
            
        return wt_cavities, mutant_cavities
    
    def save_aligned_trajectories(self, output_dir="aligned_trajectories", format="xtc"):
        """
        Guarda las trayectorias alineadas
        
        Parameters:
        -----------
        output_dir : str
            Directorio donde guardar las trayectorias alineadas
        format : str
            Formato de salida ('xtc', 'dcd', 'pdb', etc.)
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        print(f"Guardando trayectorias alineadas en formato {format.upper()} en {output_dir}/")
        
        # Seleccionar todos los átomos
        wt_all = self.u_wt.select_atoms("all")
        mutant_all = self.u_mutant.select_atoms("all")
        
        # Guardar WT alineado
        if format.lower() == "xtc":
            with mda.Writer(f"{output_dir}/wt_aligned.xtc", wt_all.n_atoms) as W:
                for ts in self.u_wt.trajectory:
                    W.write(wt_all)
        elif format.lower() == "dcd":
            with mda.Writer(f"{output_dir}/wt_aligned.dcd", wt_all.n_atoms) as W:
                for ts in self.u_wt.trajectory:
                    W.write(wt_all)
        else:  # PDB o otros formatos
            with mda.Writer(f"{output_dir}/wt_aligned.{format}", multiframe=True) as W:
                for ts in self.u_wt.trajectory:
                    W.write(wt_all)
                
        # Guardar mutante alineado
        if format.lower() == "xtc":
            with mda.Writer(f"{output_dir}/mutant_aligned.xtc", mutant_all.n_atoms) as W:
                for ts in self.u_mutant.trajectory:
                    W.write(mutant_all)
        elif format.lower() == "dcd":
            with mda.Writer(f"{output_dir}/mutant_aligned.dcd", mutant_all.n_atoms) as W:
                for ts in self.u_mutant.trajectory:
                    W.write(mutant_all)
        else:  # PDB o otros formatos
            with mda.Writer(f"{output_dir}/mutant_aligned.{format}", multiframe=True) as W:
                for ts in self.u_mutant.trajectory:
                    W.write(mutant_all)
        
        # Guardar también las topologías de referencia (primer frame)
        self.u_wt.trajectory[0]
        wt_all.write(f"{output_dir}/wt_topology.pdb")
        
        self.u_mutant.trajectory[0]
        mutant_all.write(f"{output_dir}/mutant_topology.pdb")
                
        print("✅ Trayectorias alineadas guardadas")
        print(f"   - Trayectorias: wt_aligned.{format}, mutant_aligned.{format}")
        print(f"   - Topologías: wt_topology.pdb, mutant_topology.pdb")

# Función principal para probar el script
def main():
    """
    Función principal - ajusta las rutas según tus archivos
    """
    # AJUSTAR ESTAS RUTAS A TUS ARCHIVOS
    wt_topology = "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS_extended.pdb"  # o .psf, .prmtop, etc.
    wt_trajectory = "/HDD/sromero/MCOA/Cavity_Analizer/WT_all_replicas_allargat_CAVITY_ANALYSIS.xtc"  # o .dcd, .trr, etc.
    mutant_topology = "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER_extended.pdb"
    mutant_trajectory = "/HDD/sromero/MCOA/Cavity_Analizer/15e6_allargat_al_replicas_CAVITY_ANALYZER.xtc"
    
    try:
        # Crear analizador
        analyzer = CavityAnalyzer(wt_topology, wt_trajectory, 
                                mutant_topology, mutant_trajectory)
        
        # Alinear trayectorias
        analyzer.align_trajectories()
        
        # Verificar algunas cavidades
        analyzer.get_cavity_info(frame_idx=0)
        analyzer.get_cavity_info(frame_idx=len(analyzer.u_wt.trajectory)//2)
        
        # Guardar trayectorias alineadas en formato XTC
        analyzer.save_aligned_trajectories(format="xtc")
        
        print("\n🎉 Parte 1 completada exitosamente!")
        print("Próximo paso: Clustering de trayectorias")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("Verifica las rutas de los archivos y que MDAnalysis esté instalado")

if __name__ == "__main__":
    main()
