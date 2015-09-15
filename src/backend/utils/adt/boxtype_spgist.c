#include "postgres.h"

#include "access/spgist.h"
#include "access/stratnum.h"
#include "catalog/pg_type.h"
#include "utils/builtins.h"
#include "utils/datum.h"
#include "utils/geo_decls.h"

const int NegInf = -1;
const int PosInf = 1;
const int NotInf = 0;
const int LT = -1;
const int GT = 1;
const int EQ = 0;

// InfR объединяет в себя числа из R и +- бесконечность

typedef struct {
	int infFlag;
	double val;
} InfR;

// Интервал, возможно с концами в бесконечности
typedef struct {
	InfR low;
	InfR high;
} Interval;

// 2-х мерный куб - есть Interval x Interval
typedef struct {
	Interval low_interval;     // интервал изменения левого конца отрезка
	Interval high_interval;    // интервал изменения правого конца отрезка
} BoundBox2D;

// 4-х мерный гиперкуб - есть BoundBox2D x BoundBox2D
typedef struct {
	BoundBox2D x_bound_box2d;
	BoundBox2D y_bound_box2d;
} BoundBox4D;	

// Отрезки представляем точками в 2-х мерном пространстве
typedef struct {
	double low;
	double high;
} PInterval;

// Прямоугольник - декартово произведение двух отрезков
typedef struct {
	PInterval x_pinterval;
	PInterval y_pinterval;
} PRectangle;

inline static InfR toInfR(double v){
	InfR r;
	r.infFlag = NotInf;
	r.val = v;
	return r;
}

inline static InfR negInf (void){
	InfR r;
	r.infFlag = NegInf;
	r.val = 0;
	return r;
}

inline static InfR posInf (void){
	InfR r;
	r.infFlag = PosInf;
	r.val = 0;
	return r;
}

/*
  quadrant is 8bits unsigned integer with bits: 
  [0,0,0,0,a,b,c,d] where a is one if inBox->low.x > centroid->low.x
                          b is one if inBox->high.x > centroid->high.x
                          c is one if inBox->low.y > centroid->low.y
                          d is one if inBox->high.y > centroid->high.y

*/

static uint8 getQuadrant(const BOX *centroid, const BOX *inBox){
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

static int double_cmp(const void *a, const void *b){
  double x = *(double*)a;
  double y = *(double*)b;
  if (x < y) return LT;
  if (x > y) return GT;
  return EQ;
}

/* Отрезок можно представить точкой в 2-х мерном постранстве (PInterval), которая может разделить 
   2-x мерный гиперкуб (BoundBox2D) на 4-е квадранта (BoundBox2D). 
   Эта функция вычисляет один из этих квадрантов,
   который кодируется значением переменных half1, half2.
   Если half1 == 0, то мы рассматриваем ту часть пространства, в которой первая координата отрицательна.
   Если half2 == 0, то мы рассматриваем ту часть пространства, в которой вторая координата отрицательна.
*/

static BoundBox2D evalBoundBox2D(const BoundBox2D bound_box2d, const PInterval p_interval, const int half1, const int half2){
	BoundBox2D new_bound_box2d;

	if (half1 == 0)
		new_bound_box2d.low_interval.high = toInfR(p_interval.low);
	else 
		new_bound_box2d.low_interval.low = toInfR(p_interval.low);


	if (half2 == 0)
		new_bound_box2d.high_interval.high = toInfR(p_interval.high);
	else 
		new_bound_box2d.high_interval.low = toInfR(p_interval.high);

	return new_bound_box2d;
}

/* Прямоугольник можно представить точкой в 4-х мерном постранстве (PRectangle), которая может разделить 
   4-x мерный гиперкуб (BoundBox4D) на 16-ть квадрантов (BoundBox4D). 
   Эта функция вычисляет один из этих квадрантов,
   который кодируется значением битов переменной quadrant. Смотрите функцию getQuadrant
*/

static BoundBox4D evalBoundBox4D(const BoundBox4D bound_box4d, const PRectangle p_rectangle, const uint8 quadrant){
	BoundBox4D new_bound_box4d;

	const int half1 = quadrant & 0x8;
	const int half2 = quadrant & 0x4;
	const int half3 = quadrant & 0x2;
	const int half4 = quadrant & 0x1;
   
	const BoundBox2D xb = evalBoundBox2D(bound_box4d.x_bound_box2d, p_rectangle.x_pinterval, half1, half2);
	const BoundBox2D yb = evalBoundBox2D(bound_box4d.y_bound_box2d, p_rectangle.y_pinterval, half3, half4);

	new_bound_box4d.x_bound_box2d = xb;
	new_bound_box4d.y_bound_box2d = yb;

	return new_bound_box4d;
}

// Инициализация 1-мерного гиперкуба (Interval), который покроет все пространство
inline static Interval allBoundInterval(void){
	Interval interval;
	interval.low = negInf();
	interval.high = posInf();
	return interval;
}

// Инициализация 2-х мерного гиперкуба, который покроет все пространство.
inline static BoundBox2D allBoundBox2D(void){
	BoundBox2D bound_box2d;
	bound_box2d.low_interval = allBoundInterval();
	bound_box2d.high_interval = allBoundInterval();
	return bound_box2d;
}

// Инициализация 4-х мерного гиперкуба, который покроет все пространство.
inline static BoundBox4D allBoundBox4D(void)
{
	BoundBox4D bound_box4d;
	bound_box4d.x_bound_box2d = allBoundBox2D();
	bound_box4d.y_bound_box2d = allBoundBox2D();
	return bound_box4d;
}

static int cmp_InfR_r(const InfR infVal, const double val){
	if (infVal.infFlag == PosInf)
		return GT;

	else if(infVal.infFlag == NegInf)
		return LT;

	else {
		double val0 = infVal.val;
		if(val0 < val) return LT;
		if(val0 > val) return GT;
	}

	return EQ;
}

static int intersect2D(const PInterval p_interval, const BoundBox2D bound_box2d){
	const InfR x0 = bound_box2d.low_interval.low;
//	const InfR x1 = bound_box2d.low_interval.high;

//	const InfR y0 = bound_box2d.high_interval.low;
	const InfR y1 = bound_box2d.high_interval.high;

	const double a = p_interval.low;
	const double b = p_interval.high;

	if ((cmp_InfR_r(y1, a) == LT) || (cmp_InfR_r(x0, b) == GT))
		return 0;
	else
		return 1;
}

static int intersect4D(const PRectangle p_rectangle, const BoundBox4D bound_box4d){
	const int px = intersect2D(p_rectangle.x_pinterval, bound_box4d.x_bound_box2d);
	const int py = intersect2D(p_rectangle.y_pinterval, bound_box4d.y_bound_box2d);
	return (px && py);
}

inline static PRectangle boxPointerToPRectangle(BOX *box){
	PRectangle p_rectangle;

	p_rectangle.x_pinterval.low = box->low.x;
	p_rectangle.x_pinterval.high = box->high.x;

	p_rectangle.y_pinterval.low = box->low.y;
	p_rectangle.y_pinterval.high = box->high.y;

	return p_rectangle;
}


/* SP-GiST API functions */
Datum		spg_box_quad_config(PG_FUNCTION_ARGS);
Datum		spg_box_quad_choose(PG_FUNCTION_ARGS);
Datum		spg_box_quad_picksplit(PG_FUNCTION_ARGS);
Datum		spg_box_quad_inner_consistent(PG_FUNCTION_ARGS);
Datum		spg_box_quad_leaf_consistent(PG_FUNCTION_ARGS);

/*
  TODO: подсказка
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


Datum
spg_box_quad_choose(PG_FUNCTION_ARGS)
{
	const spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
	spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);

	const BOX *inBox = DatumGetBoxP(in->datum);
	const BOX *centroid = DatumGetBoxP(in->prefixDatum);

	uint8 quadrant;
	
	if (in->allTheSame)
	{
		out->resultType = spgMatchNode;
		/* nodeN will be set by core */
		out->result.matchNode.levelAdd = 0;
		out->result.matchNode.restDatum = BoxPGetDatum(inBox);
		PG_RETURN_VOID();
	}

	quadrant = getQuadrant(centroid, inBox);
	
	out->resultType = spgMatchNode;
	out->result.matchNode.nodeN = quadrant;
	out->result.matchNode.levelAdd = 1;
	out->result.matchNode.restDatum = BoxPGetDatum(inBox);
	PG_RETURN_VOID();
}


Datum
spg_box_quad_picksplit(PG_FUNCTION_ARGS)
{
	const spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
	spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

	BOX  *centroid;
	int median, i;
	
	double *lowXs  = palloc(sizeof(double) * in->nTuples);
	double *highXs = palloc(sizeof(double) * in->nTuples);
	double *lowYs  = palloc(sizeof(double) * in->nTuples);
	double *highYs = palloc(sizeof(double) * in->nTuples);
		
	for (i = 0; i < in->nTuples; i++)
	{
		const BOX *box = DatumGetBoxP(in->datums[i]);
		lowXs[i]  = box->low.x;
		highXs[i] = box->high.x;
		lowYs[i]  = box->low.y;
		highYs[i] = box->high.y;
		i++;
	}

	qsort(lowXs, in->nTuples, sizeof(double), double_cmp);
	qsort(highXs, in->nTuples, sizeof(double), double_cmp);
	qsort(lowYs, in->nTuples, sizeof(double), double_cmp);
	qsort(highYs, in->nTuples, sizeof(double), double_cmp);

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
		const BOX *box = DatumGetBoxP(in->datums[i]);
		const uint8 quadrant = getQuadrant(centroid, box);

		out->leafTupleDatums[i] = BoxPGetDatum(box);
		out->mapTuplesToNodes[i] = quadrant;
	}

	PG_RETURN_VOID();
}



Datum
spg_box_quad_inner_consistent(PG_FUNCTION_ARGS)
{
	spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
	spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
	int i;

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

			MemoryContext oldCtx;
	
			uint8 quadrant;
			const PRectangle p_rectangle_centroid = boxPointerToPRectangle(DatumGetBoxP(in->prefixDatum));
			const PRectangle p_query_rect = boxPointerToPRectangle(DatumGetBoxP(in->scankeys[i].sk_argument));

			BoundBox4D *bound_box4d;
			
			out->traversalValues = (double **) palloc(sizeof(void) * in->nNodes);
			out->nNodes = 0;
			
			if(in->traversalValue){
				bound_box4d = in->traversalValue;
			} else {
				bound_box4d = (BoundBox4D *)palloc(sizeof(BoundBox4D));
				*bound_box4d = allBoundBox4D();
			}


			// Переключаем контекст для аллокации под traversalValue в traversalMemoryContext
			oldCtx = MemoryContextSwitchTo(in->traversalMemoryContext);

			for (quadrant = 0; quadrant < in->nNodes; quadrant++){
				BoundBox4D *new_bound_box4d;

				new_bound_box4d = (BoundBox4D *)palloc(sizeof(BoundBox4D));
				*new_bound_box4d = evalBoundBox4D(*bound_box4d, p_rectangle_centroid, quadrant);

				out->traversalValues[quadrant] = new_bound_box4d;
			}

			
			MemoryContextSwitchTo(oldCtx);
			
			for (quadrant = 0; quadrant < in->nNodes; quadrant++)
			{
				if(intersect4D(p_query_rect, *bound_box4d))
					out->nodeNumbers[out->nNodes++] = quadrant;
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
	int res, i, retval;
	
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
			retval = DatumGetBool(DirectFunctionCall2(box_overlap,
													  PointerGetDatum(leafBox),
													  keyDatum));
			/* TODO: в чем отличие между PointerGetDatum и BoxPGetDatum???
			   Какое лучше использовать тут ?
			*/
			
			PG_RETURN_BOOL(retval);
		}
	}

	PG_RETURN_BOOL(false);
}
