:-style_check(-discontiguous).

%max_vars(2).
max_body(2).


head_pred(cancer,1).
%body_pred(atm,7).
body_pred(mean_radius,2).
body_pred(mean_texture,2).
body_pred(mean_perimeter,2).
body_pred(mean_area,2).
body_pred(mean_smoothness,2).
body_pred(mean_compactness,2).
body_pred(mean_concavity,2).
body_pred(mean_concave_points,2).
body_pred(mean_symmetry,2).
body_pred(mean_fractal_dimension,2).
body_pred(radius_error,2).
body_pred(texture_error,2).
body_pred(perimeter_error,2).
body_pred(area_error,2).
body_pred(smoothness_error,2).
body_pred(compactness_error,2).
body_pred(concavity_error,2).
body_pred(concave_points_error,2).
body_pred(symmetry_error,2).
body_pred(fractal_dimension_error,2).
body_pred(worst_radius,2).
body_pred(worst_texture,2).
body_pred(worst_perimeter,2).
body_pred(worst_area,2).
body_pred(worst_smoothness,2).
body_pred(worst_compactness,2).
body_pred(worst_concavity,2).
body_pred(worst_concave_points,2).
body_pred(worst_symmetry,2).
body_pred(worst_fractal_dimension,2).

type(cancer,(pat,)).


type(mean_radius,(pat,real)).
type(mean_texture,(pat,real)).
type(mean_perimeter,(pat,real)).
type(mean_area,(pat,real)).
type(mean_smoothness,(pat,real)).
type(mean_compactness,(pat,real)).
type(mean_concavity,(pat,real)).
type(mean_concave_points,(pat,real)).
type(mean_symmetry,(pat,real)).
type(mean_fractal_dimension,(pat,real)).
type(radius_error,(pat,real)).
type(texture_error,(pat,real)).
type(perimeter_error,(pat,real)).
type(area_error,(pat,real)).
type(smoothness_error,(pat,real)).
type(compactness_error,(pat,real)).
type(concavity_error,(pat,real)).
type(concave_points_error,(pat,real)).
type(symmetry_error,(pat,real)).
type(fractal_dimension_error,(pat,real)).
type(worst_radius,(pat,real)).
type(worst_texture,(pat,real)).
type(worst_perimeter,(pat,real)).
type(worst_area,(pat,real)).
type(worst_smoothness,(pat,real)).
type(worst_compactness,(pat,real)).
type(worst_concavity,(pat,real)).
type(worst_concave_points,(pat,real)).
type(worst_symmetry,(pat,real)).
type(worst_fractal_dimension,(pat,real)).

direction(cancer,(in,)).
%direction(atm,(in,out,out,out,out,out,out)).

direction(mean_radius,(in,out)).
direction(mean_texture,(in,out)).
direction(mean_perimeter,(in,out)).
direction(mean_area,(in,out)).
direction(mean_smoothness,(in,out)).
direction(mean_compactness,(in,out)).
direction(mean_concavity,(in,out)).
direction(mean_concave_points,(in,out)).
direction(mean_symmetry,(in,out)).
direction(mean_fractal_dimension,(in,out)).
direction(radius_error,(in,out)).
direction(texture_error,(in,out)).
direction(perimeter_error,(in,out)).
direction(area_error,(in,out)).
direction(smoothness_error,(in,out)).
direction(compactness_error,(in,out)).
direction(concavity_error,(in,out)).
direction(concave_points_error,(in,out)).
direction(symmetry_error,(in,out)).
direction(fractal_dimension_error,(in,out)).
direction(worst_radius,(in,out)).
direction(worst_texture,(in,out)).
direction(worst_perimeter,(in,out)).
direction(worst_area,(in,out)).
direction(worst_smoothness,(in,out)).
direction(worst_compactness,(in,out)).
direction(worst_concavity,(in,out)).
direction(worst_concave_points,(in,out)).
direction(worst_symmetry,(in,out)).
direction(worst_fractal_dimension,(in,out)).


magic_value_type(mean_radius).
magic_value_type(mean_texture).
magic_value_type(mean_perimeter).
magic_value_type(mean_area).
magic_value_type(mean_compactness).
magic_value_type(concavity_error).
magic_value_type(concavity_error).
magic_value_type(worst_radius).
magic_value_type(worst_smoothness).
magic_value_type(worst_concavity).
magic_value_type(worst_concave_points).
magic_value_type(worst_symmetry).
magic_value_type(worst_fractal_dimension).

numerical_pred(geq,2).
numerical_pred(leq,2).

type(geq,(real,real)).
type(leq,(real,real)).

direction(geq,(in, out)).
direction(leq,(in, out)).

numerical_pred(add,3).
type(add,(real, real, real)).
direction(add,(in,in,out)).

numerical_pred(mult,3).
type(mult,(real, real, real)).
direction(mult,(in,out,in)).

    
bounds(geq,1,(0,200)).
bounds(leq,1,(0,200)).
bounds(mult,1,(0,200)).
bounds(add,1,(0,200)).