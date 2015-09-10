#include "postgres.h"

#include "access/spgist.h"
#include "access/stratnum.h"
#include "catalog/pg_type.h"
#include "utils/builtins.h"
#include "utils/datum.h"

/* SP-GiST API functions */
Datum		spg_box_quad_config(PG_FUNCTION_ARGS);
Datum		spg_box_quad_choose(PG_FUNCTION_ARGS);
Datum		spg_box_quad_picksplit(PG_FUNCTION_ARGS);
Datum		spg_box_quad_inner_consistent(PG_FUNCTION_ARGS);
Datum		spg_box_quad_leaf_consistent(PG_FUNCTION_ARGS);

/*
  #define DatumGetBoxP(X)    ((BOX *) DatumGetPointer(X)) :: Datum -> Box
  #define BoxPGetDatum(X)    PointerGetDatum(X)           :: Box -> Datum
*/


/*
 * SP-GiST 'config' interface function.
 */
Datum
spg_range_quad_config(PG_FUNCTION_ARGS)
{
	/* spgConfigIn *cfgin = (spgConfigIn *) PG_GETARG_POINTER(0); */
	spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);

	cfg->prefixType = BOXOID;
	cfg->labelType = VOIDOID;	/* we don't need node labels */
	cfg->canReturnData = true;
	cfg->longValuesOK = false;
	PG_RETURN_VOID();
}

/*
  quadrant is 8bits unsigned integer with bits: 
  [0,0,0,0,a,b,c,d] where a is one if inBox->low.x > centroid->low.x
                          b is one if inBox->high.x > centroid->high.x
                          c is one if inBox->low.y > centroid->low.y
                          d is one if inBox->high.y > centroid->high.y

*/

static uint8 getQuadrant(centroid, inBox){
	uint8 quadrant = 0;
	
	if (inBox->low.x > centroid->low.x)
		quadrant |= 1;

	quadrant <<= 1;
	
	if (inBox->high.x > centroid->high.x)
		quadrant |= 1;

	quadrant <<= 1;

	if (inBox->low.y > centroid->low.y)
		quadrant |= 1;

	quadrant <<= 1;
	
	if (inBox->high.y > centroid->high.y)
		quadrant |= 1;

	return quadrant;
}

Datum
spg_box_quad_choose(PG_FUNCTION_ARGS)
{
	spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
	spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
	BOX *inBox = DatumGetBoxP(in->datum), /* TODO what's difference between DatimGetBoxP and BoxPGetDatum */ 
		*centroid;

	uint8 quadrant;
	
	
	if (in->allTheSame)
	{
		out->resultType = spgMatchNode;
		/* nodeN will be set by core */
		out->result.matchNode.levelAdd = 0;
		out->result.matchNode.restDatum = BoxPGetDatum(inBox);
		PG_RETURN_VOID();
	}

	centroid = DatumGetBoxP(in->prefixDatum);
	quadrant = getQuadrant(centroid, inBox);
	
	out->resultType = spgMatchNode;
	out->result.matchNode.nodeN = quadrant;
	out->result.matchNode.levelAdd = 1;
	out->result.matchNode.restDatum = BoxPGetDatum(inBox);
	PG_RETURN_VOID();
}

static int double_cmp (const void *a, const void *b, const void *arg){
	if( *((double *) a) <= *((double *) b) ){
		return 1;
	} else {
		return -1;
	}
}

Datum
spg_range_quad_picksplit(PG_FUNCTION_ARGS)
{
	spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
	spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
	BOX  *centroid;
	int median;
	
	lowerBounds = palloc(sizeof(RangeBound) * in->nTuples);
	upperBounds = palloc(sizeof(RangeBound) * in->nTuples);
	
	lowXs  = palloc(sizeof(Double) * in->nTuples);
	highXs = palloc(sizeof(Double) * in->nTuples);
	lowYs  = palloc(sizeof(Double) * in->nTuples);
	highYs = palloc(sizeof(Double) * in->nTuples);
		
	for (int i = 0; i < in->nTuples; i++)
	{
		BOX *box = DatumGetBoxP(in->datums[i]);
		lowXs[i]  = box->low.x;
		highXs[i] = box->high.x;
		lowYs[i]  = box->low.y;
		highYs[i] = box->high.y;
		i++;
	}

	qsort_arg(lowXs, in->nTuples, sizeof(Double), double_cmp, NULL); // TODO зачем мне эти NULL и qsort_arg?
	qsort_arg(highXs, in->nTuples, sizeof(Double), double_cmp, NULL);
	qsort_arg(lowYs, in->nTuples, sizeof(Double), double_cmp, NULL);
	qsort_arg(highYs, in->nTuples, sizeof(Double), double_cmp, NULL);

	median = in->nTuples/2;
	
	centroid = palloc(sizeof(BOX));
	centroid->low.x = lowXs[median];
	centroid->high.x = highXs[median];
	centroid->low.y = lowYs[median];
	centroid->high.y = highYs[median];
	
	
	out->hasPrefix = true;
	out->prefixDatum = BoxPGetDatum(centroid);

	out->nNodes = 16;
	out->nodeLabels = NULL;		/* we don't need node labels */

	out->mapTuplesToNodes = palloc(sizeof(int) * in->nTuples);
	out->leafTupleDatums = palloc(sizeof(Datum) * in->nTuples);

	/*
	 * Assign ranges to corresponding nodes according to quadrants relative to
	 * "centroid" range.
	 */
	for (i = 0; i < in->nTuples; i++)
	{
		BOX *box = DatumGetBoxP(in->datums[i]);
		uint8 quadrant = getQuadrant(centroid, box);

		out->leafTupleDatums[i] = BoxPGetDatum(box);
		out->mapTuplesToNodes[i] = quadrant;
	}

	PG_RETURN_VOID();
}

Datum
spg_range_quad_inner_consistent(PG_FUNCTION_ARGS)
{
	spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
	spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
	int			which;
	int			i;

	/*
	 * For adjacent search we need also previous centroid (if any) to improve
	 * the precision of the consistent check. In this case needPrevious flag
	 * is set and centroid is passed into reconstructedValues. This is not the
	 * intended purpose of reconstructedValues (because we already have the
	 * full value available at the leaf), but it's a convenient place to store
	 * state while traversing the tree.
	 */
	bool		needPrevious = false;

	if (in->allTheSame)
	{
		/* Report that all nodes should be visited */
		out->nNodes = in->nNodes;
		out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
		for (i = 0; i < in->nNodes; i++)
			out->nodeNumbers[i] = i;
		PG_RETURN_VOID();
	}

	if (!in->hasPrefix)
	{
		/*
 		 * No centroid on this inner node. Such a node has two child nodes,
		 * the first for empty ranges, and the second for non-empty ones.
		 */
		Assert(in->nNodes == 2);

		/*
		 * Nth bit of which variable means that (N - 1)th node should be
		 * visited. Initially all bits are set. Bits of nodes which should be
		 * skipped will be unset.
		 */
		which = (1 << 1) | (1 << 2);
		for (i = 0; i < in->nkeys; i++)
		{
			StrategyNumber strategy = in->scankeys[i].sk_strategy;
			bool		empty;

			/*
			 * The only strategy when second argument of operator is not range
			 * is RANGESTRAT_CONTAINS_ELEM.
			 */
			if (strategy != RANGESTRAT_CONTAINS_ELEM)
				empty = RangeIsEmpty(
							 DatumGetRangeType(in->scankeys[i].sk_argument));
			else
				empty = false;

			switch (strategy)
			{
				case RANGESTRAT_BEFORE:
				case RANGESTRAT_OVERLEFT:
				case RANGESTRAT_OVERLAPS:
				case RANGESTRAT_OVERRIGHT:
				case RANGESTRAT_AFTER:
				case RANGESTRAT_ADJACENT:
					/* These strategies return false if any argument is empty */
					if (empty)
						which = 0;
					else
						which &= (1 << 2);
					break;

				case RANGESTRAT_CONTAINS:

					/*
					 * All ranges contain an empty range. Only non-empty
					 * ranges can contain a non-empty range.
					 */
					if (!empty)
						which &= (1 << 2);
					break;

				case RANGESTRAT_CONTAINED_BY:

					/*
					 * Only an empty range is contained by an empty range.
					 * Both empty and non-empty ranges can be contained by a
					 * non-empty range.
					 */
					if (empty)
						which &= (1 << 1);
					break;

				case RANGESTRAT_CONTAINS_ELEM:
					which &= (1 << 2);
					break;

				case RANGESTRAT_EQ:
					if (empty)
						which &= (1 << 1);
					else
						which &= (1 << 2);
					break;

				default:
					elog(ERROR, "unrecognized range strategy: %d", strategy);
					break;
			}
			if (which == 0)
				break;			/* no need to consider remaining conditions */
		}
	}
	else
	{
		RangeBound	centroidLower,
					centroidUpper;
		bool		centroidEmpty;
		TypeCacheEntry *typcache;
		RangeType  *centroid;

		/* This node has a centroid. Fetch it. */
		centroid = DatumGetRangeType(in->prefixDatum);
		typcache = range_get_typcache(fcinfo,
							   RangeTypeGetOid(DatumGetRangeType(centroid)));
		range_deserialize(typcache, centroid, &centroidLower, &centroidUpper,
						  &centroidEmpty);

		Assert(in->nNodes == 4 || in->nNodes == 5);

		/*
		 * Nth bit of which variable means that (N - 1)th node (Nth quadrant)
		 * should be visited. Initially all bits are set. Bits of nodes which
		 * can be skipped will be unset.
		 */
		which = (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5);

		for (i = 0; i < in->nkeys; i++)
		{
			StrategyNumber strategy;
			RangeBound	lower,
						upper;
			bool		empty;
			RangeType  *range = NULL;

			RangeType  *prevCentroid = NULL;
			RangeBound	prevLower,
						prevUpper;
			bool		prevEmpty;

			/* Restrictions on range bounds according to scan strategy */
			RangeBound *minLower = NULL,
					   *maxLower = NULL,
					   *minUpper = NULL,
					   *maxUpper = NULL;

			/* Are the restrictions on range bounds inclusive? */
			bool		inclusive = true;
			bool		strictEmpty = true;
			int			cmp,
						which1,
						which2;

			strategy = in->scankeys[i].sk_strategy;

			/*
			 * RANGESTRAT_CONTAINS_ELEM is just like RANGESTRAT_CONTAINS, but
			 * the argument is a single element. Expand the single element to
			 * a range containing only the element, and treat it like
			 * RANGESTRAT_CONTAINS.
			 */
			if (strategy == RANGESTRAT_CONTAINS_ELEM)
			{
				lower.inclusive = true;
				lower.infinite = false;
				lower.lower = true;
				lower.val = in->scankeys[i].sk_argument;

				upper.inclusive = true;
				upper.infinite = false;
				upper.lower = false;
				upper.val = in->scankeys[i].sk_argument;

				empty = false;

				strategy = RANGESTRAT_CONTAINS;
			}
			else
			{
				range = DatumGetRangeType(in->scankeys[i].sk_argument);
				range_deserialize(typcache, range, &lower, &upper, &empty);
			}

			/*
			 * Most strategies are handled by forming a bounding box from the
			 * search key, defined by a minLower, maxLower, minUpper,
			 * maxUpper. Some modify 'which' directly, to specify exactly
			 * which quadrants need to be visited.
			 *
			 * For most strategies, nothing matches an empty search key, and
			 * an empty range never matches a non-empty key. If a strategy
			 * does not behave like that wrt. empty ranges, set strictEmpty to
			 * false.
			 */
			switch (strategy)
			{
				case RANGESTRAT_BEFORE:

					/*
					 * Range A is before range B if upper bound of A is lower
					 * than lower bound of B.
					 */
					maxUpper = &lower;
					inclusive = false;
					break;

				case RANGESTRAT_OVERLEFT:

					/*
					 * Range A is overleft to range B if upper bound of A is
					 * less or equal to upper bound of B.
					 */
					maxUpper = &upper;
					break;

				case RANGESTRAT_OVERLAPS:

					/*
					 * Non-empty ranges overlap, if lower bound of each range
					 * is lower or equal to upper bound of the other range.
					 */
					maxLower = &upper;
					minUpper = &lower;
					break;

				case RANGESTRAT_OVERRIGHT:

					/*
					 * Range A is overright to range B if lower bound of A is
					 * greater or equal to lower bound of B.
					 */
					minLower = &lower;
					break;

				case RANGESTRAT_AFTER:

					/*
					 * Range A is after range B if lower bound of A is greater
					 * than upper bound of B.
					 */
					minLower = &upper;
					inclusive = false;
					break;

				case RANGESTRAT_ADJACENT:
					if (empty)
						break;	/* Skip to strictEmpty check. */

					/*
					 * Previously selected quadrant could exclude possibility
					 * for lower or upper bounds to be adjacent. Deserialize
					 * previous centroid range if present for checking this.
					 */
					if (in->reconstructedValue != (Datum) 0)
					{
						prevCentroid = DatumGetRangeType(in->reconstructedValue);
						range_deserialize(typcache, prevCentroid,
										  &prevLower, &prevUpper, &prevEmpty);
					}

					/*
					 * For a range's upper bound to be adjacent to the
					 * argument's lower bound, it will be found along the line
					 * adjacent to (and just below) Y=lower. Therefore, if the
					 * argument's lower bound is less than the centroid's
					 * upper bound, the line falls in quadrants 2 and 3; if
					 * greater, the line falls in quadrants 1 and 4. (see
					 * adjacent_cmp_bounds for description of edge cases).
					 */
					cmp = adjacent_inner_consistent(typcache, &lower,
													&centroidUpper,
										   prevCentroid ? &prevUpper : NULL);
					if (cmp > 0)
						which1 = (1 << 1) | (1 << 4);
					else if (cmp < 0)
						which1 = (1 << 2) | (1 << 3);
					else
						which1 = 0;

					/*
					 * Also search for ranges's adjacent to argument's upper
					 * bound. They will be found along the line adjacent to
					 * (and just right of) X=upper, which falls in quadrants 3
					 * and 4, or 1 and 2.
					 */
					cmp = adjacent_inner_consistent(typcache, &upper,
													&centroidLower,
										   prevCentroid ? &prevLower : NULL);
					if (cmp > 0)
						which2 = (1 << 1) | (1 << 2);
					else if (cmp < 0)
						which2 = (1 << 3) | (1 << 4);
					else
						which2 = 0;

					/* We must chase down ranges adjacent to either bound. */
					which &= which1 | which2;

					needPrevious = true;
					break;

				case RANGESTRAT_CONTAINS:

					/*
					 * Non-empty range A contains non-empty range B if lower
					 * bound of A is lower or equal to lower bound of range B
					 * and upper bound of range A is greater or equal to upper
					 * bound of range A.
					 *
					 * All non-empty ranges contain an empty range.
					 */
					strictEmpty = false;
					if (!empty)
					{
						which &= (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4);
						maxLower = &lower;
						minUpper = &upper;
					}
					break;

				case RANGESTRAT_CONTAINED_BY:
					/* The opposite of contains. */
					strictEmpty = false;
					if (empty)
					{
						/* An empty range is only contained by an empty range */
						which &= (1 << 5);
					}
					else
					{
						minLower = &lower;
						maxUpper = &upper;
					}
					break;

				case RANGESTRAT_EQ:

					/*
					 * Equal range can be only in the same quadrant where
					 * argument would be placed to.
					 */
					strictEmpty = false;
					which &= (1 << getQuadrant(typcache, centroid, range));
					break;

				default:
					elog(ERROR, "unrecognized range strategy: %d", strategy);
					break;
			}

			if (strictEmpty)
			{
				if (empty)
				{
					/* Scan key is empty, no branches are satisfying */
					which = 0;
					break;
				}
				else
				{
					/* Shouldn't visit tree branch with empty ranges */
					which &= (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4);
				}
			}

			/*
			 * Using the bounding box, see which quadrants we have to descend
			 * into.
			 */
			if (minLower)
			{
				/*
				 * If the centroid's lower bound is less than or equal to the
				 * minimum lower bound, anything in the 3rd and 4th quadrants
				 * will have an even smaller lower bound, and thus can't
				 * match.
				 */
				if (range_cmp_bounds(typcache, &centroidLower, minLower) <= 0)
					which &= (1 << 1) | (1 << 2) | (1 << 5);
			}
			if (maxLower)
			{
				/*
				 * If the centroid's lower bound is greater than the maximum
				 * lower bound, anything in the 1st and 2nd quadrants will
				 * also have a greater than or equal lower bound, and thus
				 * can't match. If the centroid's lower bound is equal to the
				 * maximum lower bound, we can still exclude the 1st and 2nd
				 * quadrants if we're looking for a value strictly greater
				 * than the maximum.
				 */
				int			cmp;

				cmp = range_cmp_bounds(typcache, &centroidLower, maxLower);
				if (cmp > 0 || (!inclusive && cmp == 0))
					which &= (1 << 3) | (1 << 4) | (1 << 5);
			}
			if (minUpper)
			{
				/*
				 * If the centroid's upper bound is less than or equal to the
				 * minimum upper bound, anything in the 2nd and 3rd quadrants
				 * will have an even smaller upper bound, and thus can't
				 * match.
				 */
				if (range_cmp_bounds(typcache, &centroidUpper, minUpper) <= 0)
					which &= (1 << 1) | (1 << 4) | (1 << 5);
			}
			if (maxUpper)
			{
				/*
				 * If the centroid's upper bound is greater than the maximum
				 * upper bound, anything in the 1st and 4th quadrants will
				 * also have a greater than or equal upper bound, and thus
				 * can't match. If the centroid's upper bound is equal to the
				 * maximum upper bound, we can still exclude the 1st and 4th
				 * quadrants if we're looking for a value strictly greater
				 * than the maximum.
				 */
				int			cmp;

				cmp = range_cmp_bounds(typcache, &centroidUpper, maxUpper);
				if (cmp > 0 || (!inclusive && cmp == 0))
					which &= (1 << 2) | (1 << 3) | (1 << 5);
			}

			if (which == 0)
				break;			/* no need to consider remaining conditions */
		}
	}

	/* We must descend into the quadrant(s) identified by 'which' */
	out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
	if (needPrevious)
		out->reconstructedValues = (Datum *) palloc(sizeof(Datum) * in->nNodes);
	out->nNodes = 0;
	for (i = 1; i <= in->nNodes; i++)
	{
		if (which & (1 << i))
		{
			/* Save previous prefix if needed */
			if (needPrevious)
				out->reconstructedValues[out->nNodes] = in->prefixDatum;
			out->nodeNumbers[out->nNodes++] = i - 1;
		}
	}

	PG_RETURN_VOID();
}
