:-style_check(-discontiguous).
:-use_module(numlogILP).
:-
    testAccuracy('examples/bk.pl','examples/test/examples.pl',(cancer(_25960698):-worst_fractal_dimension(_25960698,_25960710),inRange(_25960710,0.09048124999999861±0.035441249999998606),worst_area(_25960698,_25960752),leq(_25960752,1032.304477612246))).
