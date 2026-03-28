/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/bin/gnina.cuda12.8 \
    -r \
    /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/protein/receptor.pdb \
    -l \
    /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/input/ligand.sdf \
    --autobox_ligand \
    /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/reference/ref_ligand.sdf \
    --autobox_add \
    5 \
    --autobox_extend \
    1 \
    --flexres \
    A:180,A:181,A:215,A:235,A:240 \
    --num_modes \
    10 \
    --exhaustiveness \
    64 \
    --covalent_rec_atom \
    A:301:ZN \
    --covalent_lig_atom_pattern \
    [OX1;$([O]C=O)] \
    --covalent_lig_atom_position \
    6.739,10.721,31.893 \
    --covalent_fix_lig_atom_position \
    --covalent_optimize_lig \
    --cnn_scoring \
    none \
    --pose_sort_order \
    energy \
    --device \
    0 \
    --seed \
    42 \
    --atom_term_data \
    -o \
    /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/output/docked.sdf \
    --out_flex \
    /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/output/flex_residues.pdb \
    --log \
    /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/logs/gnina.log

    ==============================
*** Open Babel Warning  in Init
  Cannot initialize database 'space-groups.txt' which may cause further errors.


Usage error: Invalid solitary receptor atom A:301:ZN. Check bond lengths.
              _             
             (_)            
   __ _ _ __  _ _ __   __ _ 
  / _` | '_ \| | '_ \ / _` |
 | (_| | | | | | | | | (_| |
  \__, |_| |_|_|_| |_|\__,_|
   __/ |                    
  |___/                     

gnina v1.3.2 master:f23dd2b   Built Jul 29 2025.
gnina is based on smina and AutoDock Vina.
Please cite appropriately.
==============================
*** Open Babel Warning  in Init
  Cannot initialize database 'space-groups.txt' which may cause further errors.
==============================
*** Open Babel Warning  in PerceiveBondOrders
  Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders



Usage error: Invalid solitary receptor atom A:301:ZN. Check bond lengths.


Commandline: /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/bin/gnina.cuda12.8 -r /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/protein/receptor.pdb -l /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/input/ligand.sdf --autobox_ligand /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/reference/ref_ligand.sdf --autobox_add 5 --autobox_extend 1 --flexres A:180,A:181,A:215,A:235,A:240 --num_modes 10 --exhaustiveness 64 --covalent_rec_atom A:301:ZN --covalent_lig_atom_pattern [OX1;$([O]C=O)] --covalent_lig_atom_position 6.739,10.721,31.893 --covalent_fix_lig_atom_position --covalent_optimize_lig --cnn_scoring none --pose_sort_order energy --device 0 --seed 42 --atom_term_data -o /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/output/docked.sdf --out_flex /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/output/flex_residues.pdb --log /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/ligands/LIG_0001__N-((3-(4-(benzyloxy)phenyl)-1-phenyl-1H-pyrazol-4yl)methyl)-L-methionin/logs/gnina.log

