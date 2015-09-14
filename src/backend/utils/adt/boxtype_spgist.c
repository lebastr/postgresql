#include "postgres.h"

#include "access/spgist.h"
#include "access/stratnum.h"
#include "catalog/pg_type.h"
#include "utils/builtins.h"
#include "utils/datum.h"
#include "utils/geo_decls.h"

const int NegInf = -1;
const int PosInf = 1;
const int LT = 1;
const int GE = -1;

typedef struct {
	int infFlag;
	double val;
} InfR;

typedef struct {
	InfR low;
	InfR high;
} Interval;

typedef struct {
	Interval lowX;
	Interval highX;
	Interval lowY;
	Interval highY;
} Box4D;	

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
spg_box_quad_config(PG_FUNCTION_ARGS)
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

static uint8 getQuadrant(BOX *centroid, BOX *inBox){
	uint8 quadrant = 0;
	
	if (inBox->low.x > centroid->low.x)
		quadrant |= 0x8;

	if (inBox->high.x > centroid->high.x)
		quadrant |= 0x4;

	if (inBox->low.y > centroid->low.y)
		quadrant |= 0x2;

	if (inBox->high.y > centroid->high.y)
		quadrant |= 0x1;

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
spg_box_quad_picksplit(PG_FUNCTION_ARGS)
{
	spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
	spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
	BOX  *centroid;
	int median, i;
	
	double *lowXs  = palloc(sizeof(double) * in->nTuples);
	double *highXs = palloc(sizeof(double) * in->nTuples);
	double *lowYs  = palloc(sizeof(double) * in->nTuples);
	double *highYs = palloc(sizeof(double) * in->nTuples);
		
	for (i = 0; i < in->nTuples; i++)
	{
		BOX *box = DatumGetBoxP(in->datums[i]);
		lowXs[i]  = box->low.x;
		highXs[i] = box->high.x;
		lowYs[i]  = box->low.y;
		highYs[i] = box->high.y;
		i++;
	}

	qsort_arg(lowXs, in->nTuples, sizeof(double), double_cmp, NULL); // TODO зачем мне эти NULL и qsort_arg?
	qsort_arg(highXs, in->nTuples, sizeof(double), double_cmp, NULL);
	qsort_arg(lowYs, in->nTuples, sizeof(double), double_cmp, NULL);
	qsort_arg(highYs, in->nTuples, sizeof(double), double_cmp, NULL);

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

inline static InfR toInfR(double v){
	InfR r;

	r.infFlag = 0;
	r.val = v;
	return r;
}

inline static InfR negInf (void){
	InfR r;
	r.infFlag = NegInf;
	return r;
}

inline static InfR posInf (void){
	InfR r;
	r.infFlag = PosInf;
	return r;
}

// TODO правильно ли аллоцироать память внутри служебных функций?
static Box4D * splitBox(Box4D *box4d, BOX *centroid, int node){

	Box4D *new_box4d;
	int lowXFlag, highXFlag, lowYFlag, highYFlag;
	
	lowXFlag = node & 0x8;
	highXFlag = node & 0x4;
	lowYFlag = node & 0x2;
	highYFlag = node & 0x1;
	
	if (lowXFlag == 0){
		new_box4d->lowX.high = toInfR(centroid->low.x);
	} else {
		new_box4d->lowX.low = toInfR(centroid->low.x);
	}

	if (highXFlag == 0){
		new_box4d->highX.high = toInfR(centroid->high.x);
	} else {
		new_box4d->highX.low = toInfR(centroid->high.x);
	}

	if (lowYFlag == 0){
		new_box4d->lowY.high = toInfR(centroid->low.y);
	} else {
		new_box4d->lowY.low = toInfR(centroid->low.y);
	}

	if (highYFlag == 0){
		new_box4d->highY.high = toInfR(centroid->high.y);
	} else {
		new_box4d->highY.low = toInfR(centroid->high.y);
	}

	return new_box4d;
}

static Box4D * allBoundBox4D (void)
{
	Box4D *box = (Box4D *)palloc(sizeof(Box4D));
	box->lowX.low = negInf();
	box->lowX.high = posInf();

	box->highX.low = negInf();
	box->highX.high = posInf();

	box->lowY.low = negInf();
	box->lowY.high = posInf();

	box->highY.low = negInf();
	box->highY.high = posInf();

	return box;
}

static int cmp_InfR_r(InfR inf_r, double r){
	if (inf_r.infFlag == 0){
		if (inf_r.val < r){
			return LT;
		} else {
			return GE;
		}
	} else if(inf_r.infFlag == PosInf){
		return GE;
	} else {
		return LT;
	}
}

static int intersectBox2D(double low, double high, Interval lowInterval, Interval highInterval){
	double a = low;
	double b = high;

	InfR x0 = lowInterval.low;
	InfR x1 = lowInterval.high;
	InfR y0 = highInterval.low;
	InfR y1 = highInterval.high;
	
	if (cmp_InfR_r(x1, b) == LT && cmp_InfR_r(y0, a) == GE){
		return 1;
	} else {
		if (cmp_InfR_r(y1, a) == LT && cmp_InfR_r(x0, b) == GE){
			return 0;
		} else {
			return 1;
		}
	}
}


Datum
spg_box_quad_inner_consistent(PG_FUNCTION_ARGS)
{
	spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
	spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
	int i, quadrantNumber;

	BOX *centroid, *query_box;
	Box4D *box4d;
	MemoryContext oldCtx;
	int p1,p2;
	
	if (in->allTheSame)
	{
		/* Report that all nodes should be visited */
		out->nNodes = in->nNodes;
		out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
		for (i = 0; i < in->nNodes; i++)
			out->nodeNumbers[i] = i;
		PG_RETURN_VOID();
	}
	
	for (i = 0; i < in->nkeys; i++)
	{
		StrategyNumber strategy;
		// RTOverlapStrategyNumber
		strategy = in->scankeys[i].sk_strategy;
		if (strategy == RTOverlapStrategyNumber) // TODO правильно ли я выбрал определение стратегии? 
		{
			// тут пишем основной сейчас код

			centroid = DatumGetBoxP(in->prefixDatum);

			out->traversalValues = (double **) palloc(sizeof(void) * in->nNodes);
			out->nNodes = 0;
			
			if(in->traversalValue){
				box4d = in->traversalValue;
			} else {
				box4d = allBoundBox4D();
			}

			
			query_box = DatumGetBoxP(in->scankeys[i].sk_argument);

			// Переключаем контекст для аллокации под traversalValue в traversalMemoryContext
			oldCtx = MemoryContextSwitchTo(in->traversalMemoryContext);

			for (quadrantNumber = 0; quadrantNumber < in->nNodes; quadrantNumber++){
				out->traversalValues[quadrantNumber] = splitBox(box4d, centroid, i);
			}

			
			MemoryContextSwitchTo(oldCtx);
			
			for (quadrantNumber = 0; quadrantNumber < in->nNodes; quadrantNumber++)
			{
				p1 = intersectBox2D(query_box->low.x, query_box->high.x, box4d->lowX, box4d->highX);
				p2 = intersectBox2D(query_box->low.y, query_box->high.y, box4d->lowY, box4d->highY);
				
				if(p1 == 1 && p2 == 1)
				{
					out->nodeNumbers[out->nNodes++] = quadrantNumber;
				}
			}
			
			PG_RETURN_VOID();
		}
	}

	out->nNodes = 0;
	out->nodeNumbers = NULL;
	PG_RETURN_VOID();
}

Datum
spg_box_quad_leaf_consistent(PG_FUNCTION_ARGS)
{
	spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
	spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
	BOX *leafBox = DatumGetBoxP(in->leafDatum);
	int res, i;
	
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
