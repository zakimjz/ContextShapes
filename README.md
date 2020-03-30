# ContextShapes: Protein Docking and Partial Shape Matching

ContextShapes does rigid-body protein docking. It uses a novel contextshapes data structure to represent local surface regions/shapes on the protein. All critical points on both the receptor and ligand are represented via context shapes, and the best docking is found via pair-wise matching.

Zujun Shentu, Mohammad Al Hasan, Chris Bystroff, and Mohammed J. Zaki. Context Shapes: efficient complementary shape matching for protein-protein docking. Proteins: Structure, Function and Bioinformatics, 70(3):1056–1073, February 2008. doi:10.1002/prot.21600.

Dependencies:
=============

This code requires CGAL (https://www.cgal.org).

It also requires the boost libraries: archive and serialization (https://www.boost.org/).

For the molecular surface computation you also need to install the
binaries from http://mgltools.scripps.edu/packages/MSMS/

Finally, you need sqlite3 (https://www.sqlite.org/) for storing the docking results.

How to run:
===========

0) cd to example dir. Let's denote that $EXAMPLE

1) Copy the binary codes to $EXAMPLE. 
   The content of $EXAMPLE should be:

            $EXAMPLE/dir_database/
            $EXAMPLE/dir_output/
            $EXAMPLE/dir_PDB/
            $EXAMPLE/dir_table/
           
            params.txt
            smatch
            mtgen
            compgen
   
Under $EXAMPLE/dir_PDB/, there are the PDB files: 1BJ1_l_b.pdb and 1BJ1_r_b.pdb.
   
2) Download MSMS, the molecular surface calculator, from:
    http://mgltools.scripps.edu/packages/MSMS/
   
Before MSMS can be used to generate the molecular SES surface, the PDB files
should be cleaned. Usually the disconnected segments of HETATM should be deleted
from PDB files. The PDB files 1BJ1_l_b.pdb and 1BJ1_r_b.pdb have been cleaned.

Run pdb_to_xyzrn from the MSMS package to extract geometric info of the
proteins from PDB files:

        $> pdb_to_xyzrn 1BJ1_l_b.pdb > 1BJ1_l_b.xyzrn
        $> pdb_to_xyzrn 1BJ1_r_b.pdb > 1BJ1_r_b.xyzrn

Then run MSMS to generate SES:

        $> msms -probe_radius 1.4 -density 1 -if 1BJ1_l_b.xyzrn  -of 1BJ1_l_b

The generated SES is stored in two files: 1BJ1_l_b.vert and 1BJ1_l_b.face.

Do the same thing for 1BJ1_r_b:

        $> msms -probe_radius 1.4 -density 1 -if 1BJ1_r_b.xyzrn  -of 1BJ1_r_b

If you don't want to try out MSMS at this point, the SES files for this
example have been put in $EXAMPLE/dir_PDB/.
   
3) Generate Matching Tables:

        $> mtgen 15
   
The generated matching tables are in 5 files. Move these 5
ssave_matching_table_*.txt files into $EXAMPLE/dir_table/.
   
4) Compute context shapes for 1BJ1_r_b:

        $> smatch params.txt 410 1BJ1_r_b 1
    
   Compute context shapes for 1BJ1_l_b:

        $> smatch params.txt 410 1BJ1_l_b 2 
    
The context shapes are stored in $EXAMPLE/dir_database/.
    
5) Do the docking!

        $> smatch params.txt 500 1BJ1_r_b 1BJ1_l_b
    
The results will be stored in 
$EXAMPLE/dir_output/1BJ1_r_b_vs_1BJ1_l_b_SQL_INSERT.sql      
    
6) Analyze the results.


   Copy get_rotation_matrices.sql into $EXAMPLE/dir_output/. 
    
   cd $EXAMPLE/dir_output
         
   The 1BJ1_r_b_vs_1BJ1_l_b_SQL_INSERT.sql can be imported into any Relational 
   Database. If you don't have one yet, you can use sqlite3 for free:
   http://www.sqlite.org/.
   
   Create the db file:

       $> sqlite3 1BJ1.db < 1BJ1_r_b_vs_1BJ1_l_b_SQL_INSERT.sql

       (we have included the sql file in the dir for the example)
   
   Get the rotation matrices:

       $> sqlite3 -noheader -column 1BJ1.db < get_rotation_matrices.sql > 1BJ1_matrix.txt
   
       (we have included the txt file in the dir for the example)

   Generate the predictions:

       $> compgen ../dir_PDB/1BJ1_r_b.pdb ../dir_PDB/1BJ1_l_b.pdb 1BJ1_matrix.txt 2000 1BJ1
   
