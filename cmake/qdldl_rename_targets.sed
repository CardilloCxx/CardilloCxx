# Renames qdldl's own CMakeLists.txt target names so a from-source build of
# the real upstream github.com/osqp/qdldl (what ConicXX wants) cannot collide
# with QOCO's separately vendored, API-incompatible copy of qdldl, which
# defines targets under these exact same unqualified names. \b matches a
# word boundary, so this leaves the already-renamed conicxx_qdldl* names
# alone if this patch is ever (harmlessly) re-applied.
s/\bqdldlobject\b/conicxx_qdldlobject/g
s/\bqdldlstatic\b/conicxx_qdldlstatic/g
