$general_keywords
  coordinate_system=cartesian number_of_spatial_dimensions=1
  coordinate_label_1=x sector_points sector_factors sector_polynomials
  weighting_function=legendre  reference_weight=none
  sector_matrices
$end

$cartesian_x
 reuse_space_data automate reuse_sector_information
 number_of_fixed_points=2 left_fixed_point right_fixed_point
 drop_left_function drop_right_function number_of_major_blocks=1
$end

$x_block_1
 number_of_subregions=2 left_boundary=0.d0
 right_boundary=1.0d0 default_order=20
$end

$potential
 potential_type_region_1=none 
$end 

$poisson_x
number_of_right_hand_sides=1 inhomogeneity=linear
$end










