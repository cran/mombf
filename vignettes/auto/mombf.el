(TeX-add-style-hook
 "mombf"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "Sweave"
    "amsmath"
    "amssymb"
    "bm"
    "graphicx"
    "verbatim"
    "color"
    "hyperref"
    "natbib"
    "epsf"
    "lscape")
   (TeX-add-symbols
    "qed")
   (LaTeX-add-labels
    "sec:priors"
    "def:nlp"
    "eq:ppmodel"
    "eq:nlp_from_lp"
    "eq:marglhood_nlp"
    "sec:prodvsadditive"
    "eq:additive_nlps"
    "eq:product_nlps"
    "fig:priorplot"
    "sec:varsellm"
    "sec:bma"
    "sec:block_diag"
    "fig:bmsorthopp"
    "fig:coolblock"
    "sec:bfglm")
   (LaTeX-add-environments
    "theorem"
    "lemma"
    "proposition"
    "corollary"
    "defn"
    "example")
   (LaTeX-add-bibliographies
    "references")))

