# AMRsolver - General advection solver on AMR grids


 * Octree data structure 
 * Block-based AMR
 * Dynamic mesh generation


## Things to add

 * Add neighbors - done
 * boundary flags - done
 * Make neighbors vary in refinement by 2 at the max - done
 * Add linked list for every level - done
 * Gradient based refinement - done
 * Input file
 * Output file
 * Multigrid implementation
 * Parallel implementation

## Sample mesh output
<p align="left">
  <img src="images/3d.png" width="500"/>
  <figcaption>3D view.</figcaption>
</p>
<p align="right">
  <img src="images/2d.png" width="500"/>
  <figcaption>2D view.</figcaption>
</p>
<p align="center">
  <img src="images/gradientadapt.png"/>
  <figcaption>Gradient based refinement with block nesting.</figcaption>
</p>

