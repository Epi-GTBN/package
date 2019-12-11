#include "include/rcore.h"
#include "include/globals.h"
#include "include/graph.h"
#include "include/scores.h"
#include "include/matrix.h"
#include "include/learning.h"
#include <math.h>
#include <time.h>

SEXP score_cache_fill(SEXP nodes, SEXP data, SEXP network, SEXP score,
                      SEXP extra, SEXP reference, SEXP equivalence, SEXP decomposability,
                      SEXP updated, SEXP amat, SEXP cache, SEXP blmat, SEXP debug)
{

    int *colsum = NULL, nnodes = length(nodes), lupd = length(updated);
    int *a = NULL, *upd = NULL, *b = NULL, debuglevel = isTRUE(debug);
    int i = 0, j = 0, k = 0;
    double *cache_value = NULL;
    SEXP arc, delta, op, temp;

    /* save a pointer to the adjacency matrix, the blacklist and the
   * updated nodes. */
    a = INTEGER(amat);
    b = INTEGER(blmat);
    upd = INTEGER(updated);

    /* if there are no nodes to update, return. */
    if (lupd == 0)
        return cache;

    /* set up row and column total to check for score equivalence;
   * zero means no parent nodes. */
    if (isTRUE(equivalence))
    {

        colsum = Calloc1D(nnodes, sizeof(int));

        for (i = 0; i < nnodes; i++)
            for (j = 0; j < nnodes; j++)
                colsum[j] += a[CMC(i, j, nnodes)];

    } /*THEN*/

    /* allocate and initialize the cache. */
    cache_value = REAL(cache); // 

    /* allocate a two-slot character vector. */
    PROTECT(arc = allocVector(STRSXP, 2));

    /* allocate and initialize the fake score delta. */
    PROTECT(delta = ScalarReal(0));

    /* allocate and initialize the score.delta() operator. */
    PROTECT(op = mkString("set"));

    for (i = 0; i < nnodes; i++)
    {

        for (j = 0; j < nnodes; j++)
        {

            /* incident nodes must be different from each other. */
            if (i == j)
                continue;

            /* if only one or two nodes' caches need updating, skip the rest. */
            for (k = 0; k < lupd; k++)
                if (upd[k] == j)
                    goto there;

            continue;

        there:

            /* no need to compute the score delta for blacklisted arcs. */
            if (b[CMC(i, j, nnodes)] == 1)
                continue;

            /* use score equivalence if possible to check only one orientation. */
            if (isTRUE(equivalence))
            {

                /* if the following conditions are met, look up the score delta of
          * the reverse of the current arc:
          *   1) that score delta has already been computed.
          *   2) both incident nodes have no parent, so the arc is really
          *      score equivalent (no v-structures).
          *   3) the reversed arc has not been blacklisted, as the score delta
          *      is not computed in this case. */
                if ((i > j) && (colsum[i] + colsum[j] == 0) && (b[CMC(j, i, nnodes)] == 0))
                {

                    cache_value[CMC(i, j, nnodes)] = cache_value[CMC(j, i, nnodes)];
                    continue;

                } /*THEN*/

            } /*THEN*/

            /* save the nodes incident on the arc. */
            SET_STRING_ELT(arc, 0, STRING_ELT(nodes, i));
            SET_STRING_ELT(arc, 1, STRING_ELT(nodes, j));

            /* if the arc is not present in the graph it should be added;
        * otherwise it should be removed. */
            if (a[CMC(i, j, nnodes)] == 0)
                SET_STRING_ELT(op, 0, mkChar("set"));
            else
            {
                SET_STRING_ELT(op, 0, mkChar("drop"));
            }
            /* checkpoint allocated memory. */
            /* evaluate the call to score.delta() for the arc. */
            PROTECT(temp = score_delta(arc, network, data, score, delta, reference,
                                       op, extra, decomposability));

            cache_value[CMC(i, j, nnodes)] = NUM(VECTOR_ELT(temp, 1));
            UNPROTECT(1);

            //  if (debuglevel > 0)
            //    Rprintf("* caching score delta for arc %s -> %s (%lf).\n",
            //      CHAR(STRING_ELT(nodes, i)), CHAR(STRING_ELT(nodes, j)),
            //       cache_value[CMC(i, j, nnodes)]);

        } /*FOR*/

    } /*FOR*/

    UNPROTECT(3);

    if (isTRUE(equivalence))
        Free1D(colsum);

    return cache;

} /*HC_CACHE_FILL*/

/* a single step of the optimized hill climbing (one arc addition/removal/reversal). */
SEXP hc_opt_step(SEXP amat, SEXP nodes, SEXP added, SEXP cache, SEXP reference,
                 SEXP wlmat, SEXP blmat, SEXP nparents, SEXP maxp, SEXP debug)
{

    int nnodes = length(nodes), i = 0, j = 0;
    int *am = NULL, *ad = NULL, *w = NULL, *b = NULL, debuglevel = isTRUE(debug);
    int counter = 0, update = 1, from = 0, to = 0, *path = NULL, *scratch = NULL;
    double *cache_value = NULL, temp = 0, max = 0, tol = MACHINE_TOL;
    double *mp = REAL(maxp), *np = REAL(nparents);
    SEXP bestop;

    /* allocate and initialize the return value (use FALSE as a canary value). */
    PROTECT(bestop = allocVector(VECSXP, 3));
    setAttrib(bestop, R_NamesSymbol, mkStringVec(3, "op", "from", "to"));

    /* allocate and initialize a dummy FALSE object. */
    SET_VECTOR_ELT(bestop, 0, ScalarLogical(FALSE));

    /* allocate buffers for c_has_path(). */
    path = Calloc1D(nnodes, sizeof(int));
    scratch = Calloc1D(nnodes, sizeof(int));

    /* save pointers to the numeric/integer matrices. */
    cache_value = REAL(cache);
    ad = INTEGER(added);
    am = INTEGER(amat);
    w = INTEGER(wlmat);
    b = INTEGER(blmat);

    if (debuglevel > 0)
    {

        /* count how may arcs are to be tested. */
        for (i = 0; i < nnodes * nnodes; i++)
            counter += ad[i];

        Rprintf("----------------------------------------------------------------\n");
        Rprintf("* trying to add one of %d arcs.\n", counter);

    } /*THEN*/

    for (i = 0; i < nnodes; i++)
    {

        for (j = 0; j < nnodes; j++)
        {

            /* nothing to see, move along. */
            if (ad[CMC(i, j, nnodes)] == 0)
                continue;

            /* retrieve the score delta from the cache. */
            temp = cache_value[CMC(i, j, nnodes)];

            if (debuglevel > 0)
            {

                Rprintf("  > trying to add %s -> %s.\n", NODE(i), NODE(j));
                Rprintf("    > delta between scores for nodes %s %s is %lf.\n",
                        NODE(i), NODE(j), temp);

            } /*THEN*/

            /* this score delta is the best one at the moment, so add the arc if it
       * does not introduce cycles in the graph. */
            if (temp - max > tol)
            {

                if (c_has_path(j, i, am, nnodes, nodes, FALSE, FALSE, path, scratch,
                               FALSE))
                {

                    if (debuglevel > 0)
                        Rprintf("    > not adding, introduce cycles in to the graph.\n");

                    continue;

                } /*THEN*/

                if (debuglevel > 0)
                    Rprintf("    @ adding %s -> %s.\n", NODE(i), NODE(j));

                /* update the return value. */
                bestop_update(bestop, "set", NODE(i), NODE(j));
                /* store the node indices to update the reference scores. */
                from = i;
                to = j;

                /* update the threshold score delta. */
                max = temp;

            } /*THEN*/

        } /*FOR*/

    } /*FOR*/

    if (debuglevel > 0)
    {

        /* count how may arcs are to be tested. */
        for (i = 0, counter = 0; i < nnodes * nnodes; i++)
            counter += am[i] * (1 - w[i]);

        Rprintf("----------------------------------------------------------------\n");
        Rprintf("* trying to remove one of %d arcs.\n", counter);

    } /*THEN*/

    for (i = 0; i < nnodes; i++)
    {

        for (j = 0; j < nnodes; j++)
        {

            /* nothing to see, move along. */
            if (am[CMC(i, j, nnodes)] == 0)
                continue;

            /* whitelisted arcs are not to be removed, ever. */
            if (w[CMC(i, j, nnodes)] == 1)
                continue;

            /* retrieve the score delta from the cache. */
            temp = cache_value[CMC(i, j, nnodes)];

            if (debuglevel > 0)
            {

                Rprintf("  > trying to remove %s -> %s.\n", NODE(i), NODE(j));
                Rprintf("    > delta between scores for nodes %s %s is %lf.\n",
                        NODE(i), NODE(j), temp);

            } /*THEN*/

            if (temp - max > tol)
            {

                if (debuglevel > 0)
                    Rprintf("    @ removing %s -> %s.\n", NODE(i), NODE(j));

                /* update the return value. */
                bestop_update(bestop, "drop", NODE(i), NODE(j));
                /* store the node indices to update the reference scores. */
                from = i;
                to = j;

                /* update the threshold score delta. */
                max = temp;

            } /*THEN*/

        } /*FOR*/

    } /*FOR*/

    if (debuglevel > 0)
    {

        /* count how may arcs are to be tested. */
        for (i = 0, counter = 0; i < nnodes; i++)
            for (j = 0; j < nnodes; j++)
                counter += am[CMC(i, j, nnodes)] * (1 - b[CMC(j, i, nnodes)]);

        Rprintf("----------------------------------------------------------------\n");
        Rprintf("* trying to reverse one of %d arcs.\n", counter);

    } /*THEN*/

    for (i = 0; i < nnodes; i++)
    {

        for (j = 0; j < nnodes; j++)
        {

            /* nothing to see, move along. */
            if (am[CMC(i, j, nnodes)] == 0)
                continue;

            /* don't reverse an arc if the one in the opposite direction is
       * blacklisted, ever. */
            if (b[CMC(j, i, nnodes)] == 1)
                continue;

            /* do not reverse an arc if that means violating the limit on the
       * maximum number of parents. */
            if (np[i] >= *mp)
                continue;

            /* retrieve the score delta from the cache. */
            temp = cache_value[CMC(i, j, nnodes)] + cache_value[CMC(j, i, nnodes)];
            /* nuke small values and negative zeroes. */
            if (fabs(temp) < tol)
                temp = 0;

            if (debuglevel > 0)
            {

                Rprintf("  > trying to reverse %s -> %s.\n", NODE(i), NODE(j));
                Rprintf("    > delta between scores for nodes %s %s is %lf.\n",
                        NODE(i), NODE(j), temp);

            } /*THEN*/

            if (temp - max > tol)
            {

                if (c_has_path(i, j, am, nnodes, nodes, FALSE, TRUE, path, scratch,
                               FALSE))
                {

                    if (debuglevel > 0)
                        Rprintf("    > not reversing, introduces cycles in the graph.\n");

                    continue;

                } /*THEN*/

                if (debuglevel > 0)
                    Rprintf("    @ reversing %s -> %s.\n", NODE(i), NODE(j));

                /* update the return value. */
                bestop_update(bestop, "reverse", NODE(i), NODE(j));
                /* store the node indices to update the reference scores. */
                from = i;
                to = j;
                /* both nodes' reference scores must be updated. */
                update = 2;

                /* update the threshold score delta. */
                max = temp;

            } /*THEN*/

        } /*FOR*/

    } /*FOR*/

    /* update the reference scores. */
    REAL(reference)
    [to] += cache_value[CMC(from, to, nnodes)];
    if (update == 2)
        REAL(reference)
        [from] += cache_value[CMC(to, from, nnodes)];

    Free1D(path);
    Free1D(scratch);

    UNPROTECT(1);

    return bestop;

} /*HC_OPT_STEP*/

void bestop_update(SEXP bestop, char *op, const char *from, const char *to)
{

    SET_VECTOR_ELT(bestop, 0, mkString(op));
    SET_VECTOR_ELT(bestop, 1, mkString(from));
    SET_VECTOR_ELT(bestop, 2, mkString(to));

} /*BESTOP_UPDATE*/

SEXP hc_to_be_added(SEXP arcs, SEXP blacklist, SEXP whitelist, SEXP nparents,
                    SEXP maxp, SEXP nodes, SEXP convert)
{

    int i = 0, j = 0, narcs = 0, dims = length(nodes);
    int *a = NULL, *coords = NULL;
    double *mp = REAL(maxp), *np = NULL;
    short int duplicated = 0;
    SEXP try
        , result = R_NilValue, result2;

    /* transform the arc set into an adjacency matrix, if it's not one already. */
    if (isInteger(arcs))
    {

        if ((duplicated = NAMED(arcs)) > 0)
            PROTECT(result = duplicate(arcs));

    } /*THEN*/
    else
    {

        PROTECT(result = arcs2amat(arcs, nodes));

    } /*ELSE*/

    /* dereference the adjacency matrix once and for all. */
    a = INTEGER(result);

    /* compute the number the parents of each node, unless provided. */
    if (nparents == R_NilValue)
    {

        np = Calloc1D(dims, sizeof(double));
        for (i = 0; i < dims; i++)
            for (j = 0; j < dims; j++)
                np[j] = a[CMC(i, j, dims)];

    } /*THEN*/
    else
    {

        np = REAL(nparents);

    } /*ELSE*/

    /* flip all the nondiagonal cells. */
    for (j = 0; j < dims; j++)
    {

        for (i = 0; i < dims; i++)
        {

            /* diagonal elements are always equal to zero, skip them. */
            if (i == j)
                continue;

            a[CMC(i, j, dims)] = 1 - a[CMC(i, j, dims)];

        } /*FOR*/

    } /*FOR*/

    /* if an arc is present in the graph in one direction, you cannot add it in
   * the other direction (it would be a reversal); flip both in the adjacency
   * matrix. */
    for (j = 0; j < dims; j++)
        for (i = j + 1; i < dims; i++)
            a[CMC(j, i, dims)] = a[CMC(i, j, dims)] = a[CMC(i, j, dims)] * a[CMC(j, i, dims)];

    /* if a node has already reached its maximum number parents, do not add
   * more arcs pointing to that node. */
    for (j = 0; j < dims; j++)
        if (np[j] >= *mp)
            memset(a + j * dims, '\0', dims * sizeof(int));

#define FLIP_FROM_LIST(list, value)                                         \
    if (!isNull(list))                                                      \
    {                                                                       \
        if (!isInteger(list))                                               \
        {                                                                   \
            PROTECT(try = match(nodes, list, 0));                           \
            coords = INTEGER(try);                                          \
            narcs = length(try) / 2;                                        \
            for (i = 0; i < narcs; i++)                                     \
                a[CMC(coords[i] - 1, coords[i + narcs] - 1, dims)] = value; \
            UNPROTECT(1);                                                   \
        } /*THEN*/                                                          \
        else                                                                \
        {                                                                   \
            coords = INTEGER(list);                                         \
            for (i = 0; i < dims * dims; i++)                               \
                if (coords[i] == 1)                                         \
                    a[i] = value;                                           \
        } /*ELSE*/                                                          \
    }     /*THEN*/

    /* now the blacklist gets involved. */
    FLIP_FROM_LIST(blacklist, 0);
    /* and, last but not least, the whitelist gets involved. */
    FLIP_FROM_LIST(whitelist, 1);

    if (nparents == R_NilValue)
        Free1D(np);

    /* return either the adjacency matrix or the arc set. */
    if (isTRUE(convert))
    {

        PROTECT(result2 = amat2arcs(result, nodes));

        if ((duplicated > 0) || !isInteger(arcs))
            UNPROTECT(2);
        else
            UNPROTECT(1);
        return result2;

    } /*THEN*/
    else
    {

        if ((duplicated > 0) || !isInteger(arcs))
            UNPROTECT(1);
        return result;

    } /*ELSE*/

} /*HC_TO_BE_ADDED*/

// gennerate one individual
void set_rand_seed(SEXP debug)
{
    int debuglevel = isTRUE(debug);
    if (debuglevel > 0)
        Rprintf("\n\n== in set_rand_seed: rand seed set! ==\n\n");
    srand((int)time(0)); // use srand instead of rand
}

void ga_generate_ind(SEXP amat, SEXP nodes, SEXP wlmat, SEXP blmat, SEXP debug)
{
    int nnodes = length(nodes);
    int *am = NULL, *w = NULL, *b = NULL, debuglevel = isTRUE(debug);
    int op = -1, from = 0, to = 0, *path = NULL, *scratch = NULL;
    int flag = 0;
    int max_iter = 8; /*if try max_iter times still not get different graph return the original graph*/
    if (debuglevel > 0)
        Rprintf("ga_generate_ind is called. Step into the ga_generate_ind\n");
    // SEXP result;

    // PROTECT(result = duplicate(amat));

    /* save pointers to the numeric/integer matrices. */
    am = INTEGER(amat);
    w = INTEGER(wlmat);
    b = INTEGER(blmat);

    /* allocate buffers for c_has_path(). */
    // c_has_path is in path.c
    path = Calloc1D(nnodes, sizeof(int));
    scratch = Calloc1D(nnodes, sizeof(int));

    while (max_iter--)
    {
        if (debuglevel > 0)
            Rprintf("step into the max_iter(8)-- cycle\n");

        from = rand() % nnodes;
        to = rand() % nnodes;
        /* 0:add, 1:drop, 2:reverse */
        op = rand() % 3;
        if (debuglevel > 0)
            Rprintf("args selected this time-->from: %d and to: %d and op: %d \n", from, to, op);

        if (from == to)
            continue;

        switch (op)
        {
        case 0:
            if (debuglevel > 0)
                Rprintf("in case 0\n");
            if (am[CMC(from, to, nnodes)])
            {
                flag = 0;
                if (debuglevel > 0)
                    Rprintf("in case 0 right now, arc already exist no need to add\n");
            }
            else if (b[CMC(from, to, nnodes)])
                flag = 0;
            else
            {
                if (c_has_path(to, from, am, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                {
                    flag = 0;
                    if (debuglevel > 0)
                        Rprintf("in case 0 right now, adding arc will result in c_has_path\n");
                }
                else
                {
                    am[CMC(from, to, nnodes)] = 1;
                    flag = 1;
                    if (debuglevel > 0)
                        Rprintf("    > generate a new individual by adding %s -> %s.\n", NODE(from), NODE(to));
                }
            }
            break;
        case 1:
            if (debuglevel > 0)
                Rprintf("in case 1\n");
            if (!am[CMC(from, to, nnodes)])
            {
                flag = 0;
                if (debuglevel > 0)
                    Rprintf("in case 1 right now, arc does not exist no need to drop\n");
            }
            else if (w[CMC(from, to, nnodes)])
                flag = 0;
            else
            {
                am[CMC(from, to, nnodes)] = 0;
                flag = 1;
                if (debuglevel > 0)
                    Rprintf("    > generate a new individual by dropping %s -> %s.\n", NODE(from), NODE(to));
            }
            break;
        case 2:
            if (debuglevel > 0)
                Rprintf("in case 2\n");
            if (!am[CMC(from, to, nnodes)] || am[CMC(to, from, nnodes)])
            {
                flag = 0;
                if (debuglevel > 0)
                    Rprintf("in case 2 right now, arc does not exist or exist-both-direction, no need to reverse\n");
            }
            else if (w[CMC(from, to, nnodes)] || b[CMC(to, from, nnodes)])
                flag = 0;
            else
            {
                if (c_has_path(from, to, am, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                {
                    flag = 0;
                    if (debuglevel > 0)
                        Rprintf("in case 2 right now, reversing arc will result in c_has_path\n");
                }
                else
                {
                    am[CMC(from, to, nnodes)] = 0;
                    am[CMC(to, from, nnodes)] = 1;
                    flag = 1;
                    if (debuglevel > 0)
                        Rprintf("    > generate a new individual by reversing %s -> %s.\n", NODE(from), NODE(to));
                }
            }
            break;
        }
        if (flag)
            break;
    }
    if ((debuglevel > 0) && (!flag))
        Rprintf("    > generate a new individual by copying the original one.\n");

    Free1D(path);
    Free1D(scratch);

} /*SEXP ga_generate_ind*/

void gnd_acyclic_network(SEXP amat1, SEXP amat2, SEXP nodes, SEXP wlmat, SEXP blmat, SEXP debug)
{
    int nnodes = length(nodes), i = 0, j = 0, debuglevel = isTRUE(debug);
    int *am1 = NULL, *am2 = NULL, *w = NULL, *b = NULL, *path = NULL, *scratch = NULL;

    /* save pointers to the numeric/integer matrices. */
    am1 = INTEGER(amat1);
    am2 = INTEGER(amat2);
    w = INTEGER(wlmat);
    b = INTEGER(blmat);

    /* allocate buffers for c_has_path(). */
    path = Calloc1D(nnodes, sizeof(int));
    scratch = Calloc1D(nnodes, sizeof(int));

    for (i = 0; i < nnodes; i++)
    {
        for (j = 0; j < nnodes; j++)
        {

            if (am2[CMC(i, j, nnodes)] == 1)
            {
                if (b[CMC(i, j, nnodes)] == 1)
                {
                    am2[CMC(i, j, nnodes)] = 0;
                    if (debuglevel > 0)
                    {
                        Rprintf(" > this arc is in the blacklist,drop it\n");
                    }
                    continue;
                }
            }
            else if (am2[CMC(i, j, nnodes)] == 0)
            {
                if (w[CMC(i, j, nnodes)] == 1)
                {
                    if (debuglevel > 0)
                    {
                        Rprintf(" > this arc is in the whitelist,add it\n");
                    }
                    am2[CMC(i, j, nnodes)] = 1;
                    if (c_has_path(j, i, am2, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    {
                        if (debuglevel > 0)
                        {
                            Rprintf("> not introducing cycles in to the graph, drop it\n");
                        }
                        am2[CMC(i, j, nnodes)] = 0;
                        continue;
                    }
                }
                else
                {
                    if (am1[CMC(i, j, nnodes)] == 1)
                    {
                        if (!c_has_path(j, i, am2, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                        {
                            if (debuglevel > 0)
                            {
                                Rprintf("> not introducing cycles in to the graph, add it\n");
                            }
                            am2[CMC(i, j, nnodes)] = 1;
                            continue;
                        }
                    }
                }
            }
        }
    }
    Free1D(path);
    Free1D(scratch);
}

// Exchange two columns of the matrix
// crossover which adds tabu list
void ga_cross(SEXP amat1, SEXP amat2, SEXP new_amat1, SEXP new_amat2, SEXP nodes, SEXP wlmat, SEXP blmat, SEXP tabu_list, SEXP current, SEXP debug)
{
    // Rprintf("begin ga cross\n ");
    int nnodes = length(nodes);
    int *am1 = NULL, *am2 = NULL, *new_am1 = NULL, *new_am2 = NULL, *w = NULL, *b = NULL, debuglevel = isTRUE(debug);
    int *path = NULL, *scratch = NULL;
    int flag; // a variable to control loop,continue recycle when flag is 1
    int *cur = NULL;
    int ntabu = length(tabu_list);

    // Rprintf("am1 am2 integer \n");
    am1 = INTEGER(amat1);
    am2 = INTEGER(amat2);
    cur = INTEGER(current);

    w = INTEGER(wlmat);
    b = INTEGER(blmat);

    // Rprintf("cur:%d\n ",*cur);

    // Rprintf(" begin compute narcs\n ");
    // calculate narcs first, which are the number of arcs in the two parent network
    int narcs1 = 0;
    int narcs2 = 0;
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < nnodes; j++)
        {
            if (am1[CMC(i, j, nnodes)] == 1)
                narcs1 = narcs1 + 1;
            if (am2[CMC(i, j, nnodes)] == 1)
                narcs2 = narcs2 + 1; // CMC(i,j,nnodes) = i+j*nrow
        }
    }

    do
    {
        // generate intersections randomly
        int f1 = rand() % nnodes, f2 = rand() % nnodes, s1 = rand() % nnodes, s2 = rand() % nnodes;

        // Rprintf("new am1 am2 integer \n");
        // new_am1 and new_am2 is matrix of two offspring, it needs to be initialized to the matrix of the parent network
        // Otherwise, changes to the child matrix may affect the parent matrix

        new_am1 = INTEGER(new_amat1);
        new_am2 = INTEGER(new_amat2);

        path = Calloc1D(nnodes, sizeof(int));
        scratch = Calloc1D(nnodes, sizeof(int));

        for (int i = 0; i < nnodes; i++)
        {
            // swap the first choosen col
            // if they are the same, they do not have to exchange
            if (am1[CMC(i, f1, nnodes)] == 0 && am2[CMC(i, s1, nnodes)] == 0)
                continue;
            else if (am1[CMC(i, f1, nnodes)] == 1 && am2[CMC(i, s1, nnodes)] == 1)
                continue;
            else if (am1[CMC(i, f1, nnodes)] == 0 && am2[CMC(i, s1, nnodes)] == 1)
            {
                // if the white or black list is violated after the exchange
                if (b[CMC(i, f1, nnodes)] || w[CMC(i, s1, nnodes)])
                    continue;
                // crossover can not produce rings
                if (!c_has_path(f1, i, new_am1, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am1[CMC(i, f1, nnodes)] = 1;
                if (!c_has_path(s1, i, new_am2, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am2[CMC(i, s1, nnodes)] = 0;
            }
            else if (am1[CMC(i, f1, nnodes)] == 1 && am2[CMC(i, s1, nnodes)] == 0)
            {
                if (w[CMC(i, f1, nnodes)] || b[CMC(i, s1, nnodes)])
                    continue;
                // crossover can not produce rings
                if (!c_has_path(f1, i, new_am1, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am1[CMC(i, f1, nnodes)] = 0;
                if (!c_has_path(s1, i, new_am2, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am2[CMC(i, s1, nnodes)] = 1;
            }

            // swap the second choosen col
            if (am1[CMC(i, f2, nnodes)] == 0 && am2[CMC(i, s2, nnodes)] == 0)
                continue;
            else if (am1[CMC(i, f2, nnodes)] == 1 && am2[CMC(i, s2, nnodes)] == 1)
                continue;
            else if (am1[CMC(i, f2, nnodes)] == 0 && am2[CMC(i, s2, nnodes)] == 1)
            {
                if (b[CMC(i, f2, nnodes)] || w[CMC(i, s2, nnodes)])
                    continue;
                // crossover can not produce rings
                if (!c_has_path(f2, i, new_am1, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am1[CMC(i, f2, nnodes)] = 1;
                if (!c_has_path(s2, i, new_am2, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am2[CMC(i, s2, nnodes)] = 0;
            }
            else if (am1[CMC(i, f2, nnodes)] == 1 && am2[CMC(i, s2, nnodes)] == 0)
            {
                if (w[CMC(i, f2, nnodes)] || b[CMC(i, s2, nnodes)])
                    continue;
                // crossover can not produce rings
                if (!c_has_path(f2, i, new_am1, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am1[CMC(i, f2, nnodes)] = 0;
                if (!c_has_path(s2, i, new_am2, nnodes, nodes, FALSE, FALSE, path, scratch, FALSE))
                    new_am2[CMC(i, s2, nnodes)] = 1;
            }
        }

        /* lookup in the tabu list.If it is not in the tabu list it will return 0 */
        int idx1 = tabu_match(tabu_list, cur, new_am1, &narcs1, &nnodes, debuglevel);
        if (idx1 > 0)
        {
            if (debuglevel > 0) Rprintf("   cross failure, network matches element %d in the tabu list.\n", idx1);
            flag = 1;
        }

        int idx2 = tabu_match(tabu_list, cur, new_am2, &narcs2, &nnodes, debuglevel);
        if (idx2 > 0)
        {
            if (debuglevel > 0) Rprintf("   cross failure, network matches element %d in the tabu list.\n", idx2);
            flag = 1;
        }

        if (idx1 == 0 && idx2 == 0)
        {
            flag = 0;
            // Update parents, assign new_amat to amat
            for (int i = 0; i < nnodes; i++)
            {
                for (int j = 0; j < nnodes; j++)
                {
                    am1[CMC(i, j, nnodes)] = new_am1[CMC(i, j, nnodes)];
                    am2[CMC(i, j, nnodes)] = new_am2[CMC(i, j, nnodes)];
                }
            }

            // After successful generation of children,
            // the amat of two descendants will be placed in the tabu list
            if (debuglevel > 0){
                Rprintf("   c call tabu_hash 1\n");
            }
            tabu_hash_crossover(new_amat1, nodes, tabu_list, *cur, debuglevel);
            *cur = (*cur + 1) % ntabu;
            if (debuglevel > 0){
                Rprintf("   c call tabu_hash 2\n");
            }
            tabu_hash_crossover(new_amat2, nodes, tabu_list, *cur, debuglevel);
            *cur = (*cur + 1) % ntabu;
        }

        Free1D(path);
        Free1D(scratch);

    } while (flag == 1);
}
