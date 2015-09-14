#include "postgres.h"

#include "access/spgist.h"
#include "access/stratnum.h"
#include "catalog/pg_type.h"
#include "utils/builtins.h"
#include "utils/datum.h"

// typedef struct Box4D  


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
	
	if (in->allTheSame)
	{
		/* Report that all nodes should be visited */
		out->nNodes = in->nNodes;
		out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
		for (i = 0; i < in->nNodes; i++)
			out->nodeNumbers[i] = i;
		PG_RETURN_VOID();
	}
	
	
	centroid = DatumGetBoxP(in->prefixDatum);
	
	for (int i = 0; i < in->nkeys; i++)
	{
		StrategyNumber strategy;
		// RTOverlapStrategyNumber
		strategy = in->scankeys[i].sk_strategy;
		if (strategy == RTOverlapStrategyNumber) // TODO правильно ли я выбрал определение стратегии? 
		{
			// тут пишем основной сейчас код
			

			
			PG_RETURN_VOID();
		}
	}

	out->nNodes = 0;
	out->nodeNumbers = NULL;
	PG_RETURN_VOID();
}

Datum
spg_range_quad_leaf_consistent(PG_FUNCTION_ARGS)
{
	spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
	spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
	BOX *leafBox = DatumGetBoxP(in->leafDatum);

	/* all tests are exact */
	out->recheck = false;

	/* leafDatum is what it is... */
	out->leafValue = in->leafDatum;

	/* Perform the required comparison(s) */
	res = false;
	for (i = 0; i < in->nkeys; i++)
	{
		Datum  keyDatum = in->scankeys[i].sk_argument;
		if(in->scankeys[i].sk_strategy == RTOverlapStrategyNumber){
			PG_RETURN_BOOL(box_ov(leafBox, DatumGetBoxP(keyDatum)));
		}
	}

	PG_RETURN_BOOL(false);
}
