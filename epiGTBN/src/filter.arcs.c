#include "include/rcore.h"
#include "include/graph.h"
#include "include/matrix.h"
#include "include/bn.h"

/* remove duplicate arcs from the arc set. */
SEXP unique_arcs(SEXP arcs, SEXP nodes, SEXP warn) {

  return c_unique_arcs(arcs, nodes, isTRUE(warn));

}/*UNIQUE_ARCS*/

/* C-level interface to unique_arcs. */
SEXP c_unique_arcs(SEXP arcs, SEXP nodes, int warnlevel) {

int i = 0, j = 0, k = 0, nrow = 0, uniq_rows = 0, n = length(nodes);
int *checklist = NULL;
SEXP result, try, node, dup;

  if (isNull(arcs)) {

    /* use NULL as a special jolly value which returns all possible arcs
     * given the specified node ordering. */
    nrow = n * (n - 1)/2;

    /* allocate the return value. */
    PROTECT(result = allocMatrix(STRSXP, nrow, 2));

    /* fill in the nodes' labels. */
    for (i = 0; i < n; i++) {

      node = STRING_ELT(nodes, i);

      for (j = i + 1; j < n; j++) {

        SET_STRING_ELT(result, CMC(k, 0, nrow), node);
        SET_STRING_ELT(result, CMC(k, 1, nrow), STRING_ELT(nodes, j));
        k++;

      }/*FOR*/

    }/*FOR*/

  }/*THEN*/
  else if (length(arcs) == 0) {

    /* the arc set is empty, nothing to do. */
    return arcs;

  }/*THEN*/
  else {

    /* there really is a non-empty arc set, process it. */
    nrow = length(arcs)/2;

    /* match the node labels in the arc set. */
    PROTECT(try = arc_hash(arcs, nodes, FALSE, FALSE));
    /* check which are duplicated. */
    PROTECT(dup = duplicated(try, FALSE));
    checklist = INTEGER(dup);

    /* count how many are not. */
    for (i = 0; i < nrow; i++)
      if (checklist[i] == 0)
        uniq_rows++;

    /* if there is no duplicate arc simply return the original arc set. */
    if (uniq_rows == nrow) {

      UNPROTECT(2);
      return arcs;

    }/*THEN*/
    else {

      /* warn the user if told to do so. */
      if (warnlevel > 0)
        warning("removed %d duplicate arcs.", nrow - uniq_rows);

      /* allocate and initialize the return value. */
      PROTECT(result = allocMatrix(STRSXP, uniq_rows, 2));

      /* store the correct arcs in the return value. */
      for (i = 0, k = 0; i < nrow; i++) {

        if (checklist[i] == 0) {

          SET_STRING_ELT(result, k, STRING_ELT(arcs, i));
          SET_STRING_ELT(result, k + uniq_rows, STRING_ELT(arcs, i + nrow));
          k++;

        }/*THEN*/

      }/*FOR*/

    }/*ELSE*/

  }/*ELSE*/

  /* allocate, initialize and set the column names. */
  setDimNames(result, R_NilValue, mkStringVec(2, "from", "to"));

  if (uniq_rows == 0)
    UNPROTECT(1);
  else
    UNPROTECT(3);

  return result;

}/*C_UNIQUE_ARCS*/

/* determine which arcs are undirected. */
SEXP which_undirected(SEXP arcs, SEXP nodes) {

int i = 0, nrow = length(arcs)/2, nlvls = 0;
int *coords = NULL, *id = NULL;
SEXP result, labels, try, arc_id;

  /* get the node labels from the arcs, or use those passed down from R. */
  if (isNull(nodes))
    PROTECT(labels = unique(arcs));
  else
    labels = nodes;

  nlvls = length(labels);

  /* match the node labels in the arc set. */
  PROTECT(try = match(labels, arcs, 0));
  coords = INTEGER(try);

  /* initialize the checklist. */
  PROTECT(arc_id = allocVector(INTSXP, nrow));
  id = INTEGER(arc_id);

  /* fill the checklist with the UPTRI() coordinates, which uniquely
   * identify an arc modulo its direction. */
  for (i = 0; i < nrow; i++)
    id[i] = UPTRI(coords[i], coords[i + nrow], nlvls);

  PROTECT(result = dupe(arc_id));

  if (isNull(nodes))
    UNPROTECT(4);
  else
    UNPROTECT(3);

  return result;

}/*WHICH_UNDIRECTED*/

