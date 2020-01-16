(TeX-add-style-hook
 "large_strain_mpm"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:large-stra-form"
    "sec:track-large-deform"
    "eq:def_gradient"
    "eq:update_defor_grad"
    "sec:large-stra-analys"
    "sec:eval-weight-funct"
    "sec:analyt-appr-weight"
    "sec:refin-mater-point"
    "sec:mater-point-splitt"
    "sec:mater-point-splitt-1"
    "sec:refin-comp-grid"))
 :latex)

