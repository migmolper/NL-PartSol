(TeX-add-style-hook
 "mpm"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("babel" "english") ("geometry" "left=2cm" "right=2cm" "top=2cm" "bottom=2cm")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "SECTIONS/TEX/INTRODUCTION/Introduction"
    "SECTIONS/TEX/GOVERNING_EQUATIONS/Governing_Equations"
    "SECTIONS/TEX/MATERIAL_POINT_METHOD/stresses_mpm"
    "SECTIONS/TEX/MATERIAL_POINT_METHOD/large_strain_mpm"
    "SECTIONS/TEX/ALGORITHMS/Algorithms"
    "article"
    "art10"
    "inputenc"
    "babel"
    "amsmath"
    "amsfonts"
    "amssymb"
    "makeidx"
    "graphicx"
    "fourier"
    "hyperref"
    "geometry")
   (LaTeX-add-bibliographies
    "/home/migmolper2/Documentos/PHD/BIBLIOGRAFIA/library"))
 :latex)

