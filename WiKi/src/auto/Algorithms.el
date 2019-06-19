(TeX-add-style-hook
 "Algorithms"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:algor-mater-point"
    "sec:calc-natur-coord"
    "sec:local-search-gauss"
    "eq:NormaVectorElement"
    "fig:NormalVector"
    "eq:InOutSide"
    "fig:InOutPoint"
    "eq:VertexSeachDirection"))
 :latex)

