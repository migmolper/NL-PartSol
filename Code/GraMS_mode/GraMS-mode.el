;;; GraMS-mode.el

;;; Copyright: (C) 2010 Stephane Popinet
;; 
;;     This program is free software; you can redistribute it and/or
;;     modify it under the terms of the GNU General Public License as
;;     published by the Free Software Foundation; either version 2 of
;;     the License, or (at your option) any later version.
;;     
;;     This program is distributed in the hope that it will be useful,
;;     but WITHOUT ANY WARRANTY; without even the implied warranty of
;;     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
;;     GNU General Public License for more details.
;;     
;;     You should have received a copy of the GNU General Public License
;;     along with GNU Emacs; if not, write to the Free Software
;;     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;     02110-1301 USA
;;
;; To use this package, you can save this file somewhere in your
;; load-path and put the following in your .emacs at a minimum:
;;
;; (require 'GraMS-mode)

(define-derived-mode GraMS-mode shell-script-mode "GraMS"
  "Major mode for editing GraMS simulation files."
  
  (require 'GraMS-keywords)

  ;; Define variables ;;

  (defvar GraMS-browse-base "http://gdf.sourceforge.net/wiki/index.php/"
    "First part of URL used to display documentation on the GraMS website.")

  (defvar GraMS-functions-regexp
    (eval-when-compile
      (concat "\\<" (regexp-opt GraMS-functions t) "\\>"))
    "Regular expression compiled GraMS keywords.")

  (defvar GraMS-formulations-regexp
    (eval-when-compile
      (concat "\\<" (regexp-opt GraMS-formulations t) "\\>"))
    "Regular expression compiled GraMS formulations.")

  (defvar GraMS-modules-regexp
    (eval-when-compile
      (concat "\\<" (regexp-opt GraMS-modules t) "\\>"))
    "Regular expression compiled GraMS modules.")
  
  (define-key GraMS-mode-map [mouse-2] 'GraMS-mode-mouse-2)
  (define-key GraMS-mode-map [follow-link] 'mouse-face)

  ;; Define functions ;;

  (defun GraMS-clickable-refs (limit)
    "Font-lock function which finds GraMS keywords and makes them clickable."
    (if	(re-search-forward (eval GraMS-functions-regexp) limit t)
	(progn
	  (add-text-properties (match-beginning 0) (match-end 0)
			       (list 'mouse-face 'highlight
				     'GraMS-keyword (match-string 0)
				     'help-echo "mouse-2: documentation"
				     'rear-nonsticky '(mouse-face GraMS-functions help-echo)))
	  t)))

  (defun GraMS-clickable-formulations (limit)
    "Font-lock function which finds GraMS formulations and makes them clickable."
    (if	(re-search-forward (eval GraMS-formulations-regexp) limit t)
	(progn
	  (add-text-properties (match-beginning 0) (match-end 0)
			       (list 'mouse-face 'highlight
				     'GraMS-formulations (match-string 0)
				     'help-echo "mouse-2: documentation"
				     'rear-nonsticky '(mouse-face GraMS-formulations help-echo)))
	  t)))

  (defun GraMS-clickable-modules (limit)
    "Font-lock function which finds GraMS modules and makes them clickable."
    (if	(re-search-forward (eval GraMS-modules-regexp) limit t)
	(progn
	  (add-text-properties (match-beginning 0) (match-end 0)
			       (list 'mouse-face 'highlight
				     'GraMS-modules (match-string 0)
				     'help-echo "mouse-2: documentation"
				     'rear-nonsticky '(mouse-face GraMS-modules help-echo)))
	  t)))

  (defun GraMS-comments (limit)
    "Font-lock function which finds GraMS comments."
    (re-search-forward "#.*$" limit t))

  (defconst GraMS-font-lock-keywords
    (list 
     '(GraMS-clickable-refs (0 'font-lock-function-name-face t))
     '(GraMS-clickable-formulations (0 'font-lock-type-face t))
     '(GraMS-clickable-modules (0 'font-lock-type-face t))
     '(GraMS-comments (0 'font-lock-comment-face t)))
    "Font-lock-keywords to be added when GraMS-mode is active.")

;;   (defun GraMS-url-create (ref-string module)
;;     "Returns REF-STRING without carriage returns and with spaces converted
;; to + signs, useful when creating a URL to lookup on the GraMS website."
;;     (with-temp-buffer
;;       (insert GraMS-browse-base)
;;       (if module
;; 	  (progn 
;; 	    (insert "Object_hierarchy#")
;; 	    (insert (capitalize ref-string)))
;; 	(progn 
;; 	  (unless (string= (substring ref-string 0 3) "GraMS")
;; 	    (insert "GraMS"))
;; 	  (insert ref-string)))
;;       (buffer-string)))

;;   (defun GraMS-browse-reference (reference &optional module)
;;     "Wrapper function to call standard Emacs browser function for REFERENCE."
;;     (message "Linking to GraMS website for %s..." reference)
;;     (browse-url (GraMS-url-create reference module)))
  
;;   (defun GraMS-mode-mouse-2 (event arg)
;;     "Fetch documentation for keyword under the mouse click."
;;     (interactive "e\nP")
;;     (let (my-keyword)
;;       (save-excursion
;; 	(set-buffer (window-buffer (posn-window (event-end event))))
;; 	(goto-char (posn-point (event-end event)))
;; 	(setq my-keyword (get-text-property (point) 'GraMS-keyword)))
;;       (if my-keyword
;; 	  (progn
;; 	    (select-window (posn-window (event-end event)))
;; 	    (GraMS-browse-reference my-keyword))
;; 	(progn
;; 	  (setq my-keyword (get-text-property (point) 'GraMS-module))
;; 	  (if my-keyword
;; 	      (progn
;; 		(select-window (posn-window (event-end event)))
;; 		(GraMS-browse-reference my-keyword t))
;; 	    (mouse-yank-at-click event arg)
;; 	    )))))

  (font-lock-add-keywords nil GraMS-font-lock-keywords)

  ;; load keywords for autocompletion with dabbrev ;;
  (find-file-noselect (locate-file "GraMS-keywords.el" load-path) t)
  (setq case-fold-search nil)

  (column-number-mode 1)
)

(add-to-list 'auto-mode-alist '("\\.gdf\\'" . GraMS-mode))

(provide 'GraMS-mode)
