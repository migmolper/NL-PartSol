(TeX-add-style-hook
 "Governing_Equations"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:deriv_eq"
    "eq:sf_conserv_momentum"
    "eq:sf_conserv_momentum_2"
    "eq:hilbert_space"
    "eq:def_psi"
    "eq:wf_conserv_momentum_1"
    "eq:wf_conserv_momentum_2"
    "eq:wf_conserv_momentum3"
    "sec:observations"))
 :latex)
