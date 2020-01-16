(TeX-add-style-hook
 "Material_Point_Method"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (LaTeX-add-labels
    "sec:state-art-material"
    "eq:mpm_discretization"
    "eq:mpm_density"
    "eq:specific_stress"
    "eq:specific_stress_intforces"
    "sec:orig-mater-point"
    "eq:part_uniti_chi_i"
    "eq:vpi_def"
    "eq:m_pi"
    "eq:p_pi"
    "sec:discr_sol_proc"
    "eq:mat_point_discretiz"
    "eq:wf_conserv_momentum_acc"
    "eq:wf_conserv_momentum_int_forces"
    "eq:wf_conserv_momentum_ext_forces"
    "eq:wf_conserv_momentum_GIMP"
    "eq:charact_volum"
    "eq:Partition_Unity"))
 :latex)

