/*
 * Copyright (c) 2004-2005 Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * MIT grants permission to use, copy, modify, and distribute this software and
 * its documentation for NON-COMMERCIAL purposes and without fee, provided that
 * this copyright notice appears in all copies.
 *
 * MIT provides this software "as is," without representations or warranties of
 * any kind, either expressed or implied, including but not limited to the
 * implied warranties of merchantability, fitness for a particular purpose, and
 * noninfringement.  MIT shall not be liable for any damages arising from any
 * use of this software.
 *
 * Author: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)
 */

/*
  The main functionality of the LSH scheme is in this file (all except
  the hashing of the buckets). This file includes all the functions
  for processing a PRNearNeighborStructT data structure, which is the
  main R-NN data structure based on LSH scheme. The particular
  functions are: initializing a DS, adding new points to the DS, and
  responding to queries on the DS.
 */

#include "headers.h"

void printRNNParameters(FILE *output, RNNParametersT parameters){
  ASSERT(output != NULL);
  fprintf(output, "R\n");
  fprintf(output, "%0.9lf\n", parameters.parameterR);
  fprintf(output, "Success probability\n");
  fprintf(output, "%0.9lf\n", parameters.successProbability);
  fprintf(output, "Dimension\n");
  fprintf(output, "%d\n", parameters.dimension);
  fprintf(output, "R^2\n");
  fprintf(output, "%0.9lf\n", parameters.parameterR2);
  fprintf(output, "Use <u> functions\n");
  fprintf(output, "%d\n", parameters.useUfunctions);
  fprintf(output, "k\n");
  fprintf(output, "%d\n", parameters.parameterK);
  fprintf(output, "m [# independent tuples of LSH functions]\n");
  fprintf(output, "%d\n", parameters.parameterM);
  fprintf(output, "L\n");
  fprintf(output, "%d\n", parameters.parameterL);
  fprintf(output, "W\n");
  fprintf(output, "%0.9lf\n", parameters.parameterW);
  fprintf(output, "T\n");
  fprintf(output, "%d\n", parameters.parameterT);
  fprintf(output, "typeHT\n");
  fprintf(output, "%d\n", parameters.typeHT);
}

RNNParametersT readRNNParameters(FILE *input){
  ASSERT(input != NULL);
  RNNParametersT parameters;
  char s[1000];// TODO: possible buffer overflow

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  FSCANF_REAL(input, &parameters.parameterR);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  FSCANF_REAL(input, &parameters.successProbability);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.dimension);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  FSCANF_REAL(input, &parameters.parameterR2);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.useUfunctions);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.parameterK);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.parameterM);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.parameterL);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  FSCANF_REAL(input, &parameters.parameterW);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.parameterT);

  fscanf(input, "\n");fscanf(input, "%[^\n]\n", s);
  fscanf(input, "%d", &parameters.typeHT);

  return parameters;
}

// Creates the LSH hash functions for the R-near neighbor structure
// <nnStruct>. The functions fills in the corresponding field of
// <nnStruct>.
void initHashFunctions(PRNearNeighborStructT nnStruct){
  ASSERT(nnStruct != NULL);
  LSHFunctionT **lshFunctions;
  // allocate memory for the functions
  FAILIF(NULL == (lshFunctions = (LSHFunctionT**)MALLOC(nnStruct->nHFTuples * sizeof(LSHFunctionT*)))); //根据使用哈希函数的个数L分配内存
  for(IntT i = 0; i < nnStruct->nHFTuples; i++){
    FAILIF(NULL == (lshFunctions[i] = (LSHFunctionT*)MALLOC(nnStruct->hfTuplesLength * sizeof(LSHFunctionT)))); //根据k来分配内存
    for(IntT j = 0; j < nnStruct->hfTuplesLength; j++){
      FAILIF(NULL == (lshFunctions[i][j].a = (RealT*)MALLOC(nnStruct->dimension * sizeof(RealT))));  //根据维度来分配内存，a为首地址
    }
  }

  // initialize the LSH functions
  for(IntT i = 0; i < nnStruct->nHFTuples; i++){   //将产生L个hash族函数，每个族函数由k个hash函数组成，这k个hash函数对每一维进行处理
    for(IntT j = 0; j < nnStruct->hfTuplesLength; j++){
      // vector a
      for(IntT d = 0; d < nnStruct->dimension; d++){
#ifdef USE_L1_DISTANCE
	lshFunctions[i][j].a[d] = genCauchyRandom();    //1-stable //对于L1 norm的没有实现
#else
	lshFunctions[i][j].a[d] = genGaussianRandom();   //2-stable //返回一个正态分布随机因子
#endif
      }
      // b
      lshFunctions[i][j].b = genUniformRandom(0, nnStruct->parameterW);   //产生一个随机数
    }
  }

  nnStruct->lshFunctions = lshFunctions;
}

// Initializes the fields of a R-near neighbors data structure except
// the hash tables for storing the buckets.    //对RNN数据结果哈希表存储桶的初始化
PRNearNeighborStructT initializePRNearNeighborFields(RNNParametersT algParameters, Int32T nPointsEstimate){
  PRNearNeighborStructT nnStruct;
  FAILIF(NULL == (nnStruct = (PRNearNeighborStructT)MALLOC(sizeof(RNearNeighborStructT))));
  nnStruct->parameterR = algParameters.parameterR;
  nnStruct->parameterR2 = algParameters.parameterR2;
  nnStruct->useUfunctions = algParameters.useUfunctions;
  nnStruct->parameterK = algParameters.parameterK;
  if (!algParameters.useUfunctions) {           //如果u函数为0，使用普通g函数
    // Use normal <g> functions.
    nnStruct->parameterL = algParameters.parameterL;
    nnStruct->nHFTuples = algParameters.parameterL;          //使用L个元组，即使用L个哈希函数
    nnStruct->hfTuplesLength = algParameters.parameterK;     //每个哈希函数中元组长度为k个
  }else{
    // Use <u> hash functions; a <g> function is a pair of 2 <u> functions.
    nnStruct->parameterL = algParameters.parameterL;
    nnStruct->nHFTuples = algParameters.parameterM;             //使用M个元组，hashing函数
    nnStruct->hfTuplesLength = algParameters.parameterK / 2;    //每个哈希函数中元组长度为k/2个
  }
  nnStruct->parameterT = algParameters.parameterT;              //数据集的数量
  nnStruct->dimension = algParameters.dimension;
  nnStruct->parameterW = algParameters.parameterW;

  nnStruct->nPoints = 0;
  nnStruct->pointsArraySize = nPointsEstimate;                  //根据点的个数，确定点数组的大小

  FAILIF(NULL == (nnStruct->points = (PPointT*)MALLOC(nnStruct->pointsArraySize * sizeof(PPointT))));   //内存分配失败

  // create the hash functions
  initHashFunctions(nnStruct);

  // init fields that are used only in operations ("temporary" variables for operations).

  // init the vector <pointULSHVectors> and the vector
  // <precomputedHashesOfULSHs>
  //<pointULSHVectors>是用来存储函数u<hfTuplesLength>-tuple of LSH fuctions,即k,用来连接g函数。
  //<precomputedHashesOfULSHs>Precomputed hashes of each of the <nHFTuples> of <u> functions (to be used by the bucket hashing module).
  FAILIF(NULL == (nnStruct->pointULSHVectors = (Uns32T**)MALLOC(nnStruct->nHFTuples * sizeof(Uns32T*))));
  for(IntT i = 0; i < nnStruct->nHFTuples; i++){
    FAILIF(NULL == (nnStruct->pointULSHVectors[i] = (Uns32T*)MALLOC(nnStruct->hfTuplesLength * sizeof(Uns32T))));  //使用u函数分配？
  }
  FAILIF(NULL == (nnStruct->precomputedHashesOfULSHs = (Uns32T**)MALLOC(nnStruct->nHFTuples * sizeof(Uns32T*))));    //使用普通函数分配？
  for(IntT i = 0; i < nnStruct->nHFTuples; i++){
    FAILIF(NULL == (nnStruct->precomputedHashesOfULSHs[i] = (Uns32T*)MALLOC(N_PRECOMPUTED_HASHES_NEEDED * sizeof(Uns32T))));
  }
  // init the vector <reducedPoint>
  FAILIF(NULL == (nnStruct->reducedPoint = (RealT*)MALLOC(nnStruct->dimension * sizeof(RealT))));
  // init the vector <nearPoints>
  nnStruct->sizeMarkedPoints = nPointsEstimate;
  FAILIF(NULL == (nnStruct->markedPoints = (BooleanT*)MALLOC(nnStruct->sizeMarkedPoints * sizeof(BooleanT))));
  for(IntT i = 0; i < nnStruct->sizeMarkedPoints; i++){
    nnStruct->markedPoints[i] = FALSE;
  }
  // init the vector <nearPointsIndeces>
  FAILIF(NULL == (nnStruct->markedPointsIndeces = (Int32T*)MALLOC(nnStruct->sizeMarkedPoints * sizeof(Int32T))));

  nnStruct->reportingResult = TRUE;

  return nnStruct;
}

// Constructs a new empty R-near-neighbor data structure.
PRNearNeighborStructT initLSH(RNNParametersT algParameters, Int32T nPointsEstimate){
  ASSERT(algParameters.typeHT == HT_LINKED_LIST || algParameters.typeHT == HT_STATISTICS);
  PRNearNeighborStructT nnStruct = initializePRNearNeighborFields(algParameters, nPointsEstimate);

  // initialize second level hashing (bucket hashing)
  FAILIF(NULL == (nnStruct->hashedBuckets = (PUHashStructureT*)MALLOC(nnStruct->parameterL * sizeof(PUHashStructureT))));
  Uns32T *mainHashA = NULL, *controlHash1 = NULL;
  BooleanT uhashesComputedAlready = FALSE;
  for(IntT i = 0; i < nnStruct->parameterL; i++){
    nnStruct->hashedBuckets[i] = newUHashStructure(algParameters.typeHT, nPointsEstimate, nnStruct->parameterK, uhashesComputedAlready, mainHashA, controlHash1, NULL);
    uhashesComputedAlready = TRUE;
  }

  return nnStruct;
}

void preparePointAdding(PRNearNeighborStructT nnStruct, PUHashStructureT uhash, PPointT point);


// Construct PRNearNeighborStructT given the data set <dataSet> (all
// the points <dataSet> will be contained in the resulting DS).
// Currenly only type HT_HYBRID_CHAINS is supported for this
// operation.
PRNearNeighborStructT initLSH_WithDataSet(RNNParametersT algParameters, Int32T nPoints, PPointT *dataSet){   //至298行
  ASSERT(algParameters.typeHT == HT_HYBRID_CHAINS);
  ASSERT(dataSet != NULL);
  ASSERT(USE_SAME_UHASH_FUNCTIONS);

  PRNearNeighborStructT nnStruct = initializePRNearNeighborFields(algParameters, nPoints);  //参数集，点的个数

  // Set the fields <nPoints> and <points>.
  nnStruct->nPoints = nPoints;
  for(Int32T i = 0; i < nPoints; i++){
    nnStruct->points[i] = dataSet[i];
  }
  
  // initialize second level hashing (bucket hashing)
  FAILIF(NULL == (nnStruct->hashedBuckets = (PUHashStructureT*)MALLOC(nnStruct->parameterL * sizeof(PUHashStructureT))));
  Uns32T *mainHashA = NULL, *controlHash1 = NULL;
  //BucketHashing.cpp  Line 70. 使用第1种拉链存储，不使用u函数，生成h1和h2,搭出hash Table的框架。
  PUHashStructureT modelHT = newUHashStructure(HT_LINKED_LIST, nPoints, nnStruct->parameterK, FALSE, mainHashA, controlHash1, NULL);
  
  Uns32T **(precomputedHashesOfULSHs[nnStruct->nHFTuples]);   //预处理hashing,L ，一个两维指针数组
  for(IntT l = 0; l < nnStruct->nHFTuples; l++){              //给L的每一个分配n个指针大小的内存
    FAILIF(NULL == (precomputedHashesOfULSHs[l] = (Uns32T**)MALLOC(nPoints * sizeof(Uns32T*)))); //每次hash需要这么多内存，指向类型点的指针
    for(IntT i = 0; i < nPoints; i++){         //对L中的每个点分配4个哈希值的空间
      FAILIF(NULL == (precomputedHashesOfULSHs[l][i] = (Uns32T*)MALLOC(N_PRECOMPUTED_HASHES_NEEDED * sizeof(Uns32T))));
    }
  }

  for(IntT i = 0; i < nPoints; i++){
    preparePointAdding(nnStruct, modelHT, dataSet[i]);   //对每个数据点进行处理，将d维的数据映射成整数。对点的一级映射、二级映射、物理块的映射
    for(IntT l = 0; l < nnStruct->nHFTuples; l++){   //L桶的值
      for(IntT h = 0; h < N_PRECOMPUTED_HASHES_NEEDED; h++){   //N_Precomputed_hashes_needed是Result的值，是4
	    precomputedHashesOfULSHs[l][i][h] = nnStruct->precomputedHashesOfULSHs[l][h];  //将每个映射的点映射到对应的hashing桶中
      }
    }
  }

  //DPRINTF("Allocated memory(modelHT and precomputedHashesOfULSHs just a.): %lld\n", totalAllocatedMemory);

  // Initialize the counters for defining the pair of <u> functions used for <g> functions. //初始化计数器
  IntT firstUComp = 0;
  IntT secondUComp = 1;
  for(IntT i = 0; i < nnStruct->parameterL; i++){
    // build the model HT.
    for(IntT p = 0; p < nPoints; p++){
      // Add point <dataSet[p]> to modelHT.
      if (!nnStruct->useUfunctions) {
	// Use usual <g> functions (truly independent; <g>s are precisly
	// <u>s).
	addBucketEntry(modelHT, 1, precomputedHashesOfULSHs[i][p], NULL, p);  //把点加到桶中，k  //i<=L,所以是L个独立的g函数
      } else {          // k/2,因为g函数不独立，重复计算两部分的值？？g函数是怎样不独立的？？
	// Use <u> functions (<g>s are pairs of <u> functions).
	//该函数产生第三个参数，并写入参数文件。即为-c产生，作为-p的第十个参数  //firstUComp 和 secondUComp
	addBucketEntry(modelHT, 2, precomputedHashesOfULSHs[firstUComp][p], precomputedHashesOfULSHs[secondUComp][p], p);
      }
    }

    //ASSERT(nAllocatedGBuckets <= nPoints);
    //ASSERT(nAllocatedBEntries <= nPoints);

    // compute what is the next pair of <u> functions.
    secondUComp++;  //先计算行，在计算列，所以是secondUComp++，共L行k列
    if (secondUComp == nnStruct->nHFTuples) {
      firstUComp++;
      secondUComp = firstUComp + 1;
    }

    // copy the model HT into the actual (packed) HT. copy the uhash function too.
	//之前的模型用的是HT_Linked_List,但其中的uhash 函数可以复用。现在将真正要使用的哈希表类型加进去，为混合链表类型
    nnStruct->hashedBuckets[i] = newUHashStructure(algParameters.typeHT, nPoints, nnStruct->parameterK, TRUE, mainHashA, controlHash1, modelHT);

    // clear the model HT for the next iteration.
    clearUHashStructure(modelHT);
  }

  freeUHashStructure(modelHT, FALSE); // do not free the uhash functions since they are used by nnStruct->hashedBuckets[i]

  // freeing precomputedHashesOfULSHs
  for(IntT l = 0; l < nnStruct->nHFTuples; l++){
    for(IntT i = 0; i < nPoints; i++){
      FREE(precomputedHashesOfULSHs[l][i]);
    }
    FREE(precomputedHashesOfULSHs[l]);
  }

  return nnStruct;
}



// // Packed version (static).
// PRNearNeighborStructT buildPackedLSH(RealT R, BooleanT useUfunctions, IntT k, IntT LorM, RealT successProbability, IntT dim, IntT T, Int32T nPoints, PPointT *points){
//   ASSERT(points != NULL);
//   PRNearNeighborStructT nnStruct = initializePRNearNeighborFields(R, useUfunctions, k, LorM, successProbability, dim, T, nPoints);

//   // initialize second level hashing (bucket hashing)
//   FAILIF(NULL == (nnStruct->hashedBuckets = (PUHashStructureT*)MALLOC(nnStruct->parameterL * sizeof(PUHashStructureT))));
//   Uns32T *mainHashA = NULL, *controlHash1 = NULL;
//   PUHashStructureT modelHT = newUHashStructure(HT_STATISTICS, nPoints, nnStruct->parameterK, FALSE, mainHashA, controlHash1, NULL);
//   for(IntT i = 0; i < nnStruct->parameterL; i++){
//     // build the model HT.
//     for(IntT p = 0; p < nPoints; p++){
//       // addBucketEntry(modelHT, );
//     }



//     // copy the model HT into the actual (packed) HT.
//     nnStruct->hashedBuckets[i] = newUHashStructure(HT_PACKED, nPointsEstimate, nnStruct->parameterK, TRUE, mainHashA, controlHash1, modelHT);

//     // clear the model HT for the next iteration.
//     clearUHashStructure(modelHT);
//   }

//   return nnStruct;
// }


// Optimizes the nnStruct (non-agressively, i.e., without changing the
// parameters).
void optimizeLSH(PRNearNeighborStructT nnStruct){
  ASSERT(nnStruct != NULL);

  PointsListEntryT *auxList = NULL;
  for(IntT i = 0; i < nnStruct->parameterL; i++){
    optimizeUHashStructure(nnStruct->hashedBuckets[i], auxList);
  }
  FREE(auxList);
}

// Frees completely all the memory occupied by the <nnStruct>
// structure.
void freePRNearNeighborStruct(PRNearNeighborStructT nnStruct){
  if (nnStruct == NULL){
    return;
  }

  if (nnStruct->points != NULL) {
    free(nnStruct->points);
  }
  
  if (nnStruct->lshFunctions != NULL) {
    for(IntT i = 0; i < nnStruct->nHFTuples; i++){
      for(IntT j = 0; j < nnStruct->hfTuplesLength; j++){
	free(nnStruct->lshFunctions[i][j].a);
      }
      free(nnStruct->lshFunctions[i]);
    }
    free(nnStruct->lshFunctions);
  }
  
  if (nnStruct->precomputedHashesOfULSHs != NULL) {
    for(IntT i = 0; i < nnStruct->nHFTuples; i++){
      free(nnStruct->precomputedHashesOfULSHs[i]);
    }
    free(nnStruct->precomputedHashesOfULSHs);
  }

  freeUHashStructure(nnStruct->hashedBuckets[0], TRUE);
  for(IntT i = 1; i < nnStruct->parameterL; i++){
    freeUHashStructure(nnStruct->hashedBuckets[i], FALSE);
  }
  free(nnStruct->hashedBuckets);

  if (nnStruct->pointULSHVectors != NULL){
    for(IntT i = 0; i < nnStruct->nHFTuples; i++){
      free(nnStruct->pointULSHVectors[i]);
    }
    free(nnStruct->pointULSHVectors);
  }

  if (nnStruct->reducedPoint != NULL){
    free(nnStruct->reducedPoint);
  }

  if (nnStruct->markedPoints != NULL){
    free(nnStruct->markedPoints);
  }

  if (nnStruct->markedPointsIndeces != NULL){
    free(nnStruct->markedPointsIndeces);
  }
}

// If <reportingResult> == FALSe, no points are reported back in a
// <get> function. In particular any point that is found in the bucket
// is considered to be outside the R-ball of the query point.  If
// <reportingResult> == TRUE, then the structure behaves normally.
void setResultReporting(PRNearNeighborStructT nnStruct, BooleanT reportingResult){
  ASSERT(nnStruct != NULL);
  nnStruct->reportingResult = reportingResult;
}

// Compute the value of a hash function u=lshFunctions[gNumber] (a
// vector of <hfTuplesLength> LSH functions) in the point <point>. The
// result is stored in the vector <vectorValue>. <vectorValue> must be
// already allocated (and have space for <hfTuplesLength> Uns32T-words).  //一个g函数相当k个哈希函数
//对每个L中的g函数中的k个哈希函数进行处理.*vectorValue中存储的是k个值的首地址
inline void computeULSH(PRNearNeighborStructT nnStruct, IntT gNumber, RealT *point, Uns32T *vectorValue){
  CR_ASSERT(nnStruct != NULL);
  CR_ASSERT(point != NULL);
  CR_ASSERT(vectorValue != NULL);

  for(IntT i = 0; i < nnStruct->hfTuplesLength; i++){  //一个表中有k个哈希函数。
    RealT value = 0;
    for(IntT d = 0; d < nnStruct->dimension; d++){   //用的是reducePoint。
      value += point[d] * nnStruct->lshFunctions[gNumber][i].a[d];  //对一个点做内积：每个维度的值对位内积。a[d]里面存储的是一个服从p-stable分布的随机数，b里面存储的也是一个随机数[0,w]之间
    }
  
    vectorValue[i] = (Uns32T)(FLOOR_INT32((value + nnStruct->lshFunctions[gNumber][i].b) / nnStruct->parameterW) /* - MIN_INT32T*/);  //公式(a+b)/w
  }
}
//line 243
inline void preparePointAdding(PRNearNeighborStructT nnStruct, PUHashStructureT uhash, PPointT point){  //uhash为hash通的结构，point为点的具体值
  ASSERT(nnStruct != NULL);
  ASSERT(uhash != NULL);
  ASSERT(point != NULL);

  TIMEV_START(timeComputeULSH);   //开始计时
  for(IntT d = 0; d < nnStruct->dimension; d++){
    nnStruct->reducedPoint[d] = point->coordinates[d] / nnStruct->parameterR;  //将点的维度除以半径，用于放大或者缩小点的值
  }                                                                            //对所有的点做了规范化

  // Compute all ULSH functions.
  for(IntT i = 0; i < nnStruct->nHFTuples; i++){    //进行L次，将每个d维点，映射成一个数 //nHFTuples是L，Tupleslength是k
      //生成gi的值，g_i=(h_1^((i) ),…,h_k^((i) ))
	  computeULSH(nnStruct, i, nnStruct->reducedPoint, nnStruct->pointULSHVectors[i]);  //pointULSHVectors是第i个Table的具体的哈希函数
	                                                                            //nnStruct->pointULSHVectors[i]为k个hash函数的指针头地址
 }

  // Compute data for <precomputedHashesOfULSHs>.
  if (USE_SAME_UHASH_FUNCTIONS) {
    for(IntT i = 0; i < nnStruct->nHFTuples; i++){   //重复L次，hash表及其桶的结构准备好，并将点存入桶中 /
		//下列函数运算的结果是根据哈希公式，生成k个大数存入nnStruct->precomputedHashesOfULSHs[i]中，循环L次。
      precomputeUHFsForULSH(uhash, nnStruct->pointULSHVectors[i], nnStruct->hfTuplesLength, nnStruct->precomputedHashesOfULSHs[i]);
	  //哈希表的结构，LSH表的值，k,L组哈希函数的头指针
    }
  }

  TIMEV_END(timeComputeULSH);
}

inline void batchAddRequest(PRNearNeighborStructT nnStruct, IntT i, IntT &firstIndex, IntT &secondIndex, PPointT point){
//   Uns32T *(gVector[4]);
//   if (!nnStruct->useUfunctions) {
//     // Use usual <g> functions (truly independent).
//     gVector[0] = nnStruct->pointULSHVectors[i];
//     gVector[1] = nnStruct->precomputedHashesOfULSHs[i];
//     addBucketEntry(nnStruct->hashedBuckets[firstIndex], gVector, 1, point);
//   } else {
//     // Use <u> functions (<g>s are pairs of <u> functions).
//     gVector[0] = nnStruct->pointULSHVectors[firstIndex];
//     gVector[1] = nnStruct->pointULSHVectors[secondIndex];
//     gVector[2] = nnStruct->precomputedHashesOfULSHs[firstIndex];
//     gVector[3] = nnStruct->precomputedHashesOfULSHs[secondIndex];
    
//     // compute what is the next pair of <u> functions.
//     secondIndex++;
//     if (secondIndex == nnStruct->nHFTuples) {
//       firstIndex++;
//       secondIndex = firstIndex + 1;
//     }
    
//     addBucketEntry(nnStruct->hashedBuckets[i], gVector, 2, point);
//   }
  ASSERT(1 == 0);
}

// Adds a new point to the LSH data structure, that is for each
// i=0..parameterL-1, the point is added to the bucket defined by
// function g_i=lshFunctions[i].  //将一个新的点加入到数据结构中
void addNewPointToPRNearNeighborStruct(PRNearNeighborStructT nnStruct, PPointT point){
  ASSERT(nnStruct != NULL);
  ASSERT(point != NULL);
  ASSERT(nnStruct->reducedPoint != NULL);
  ASSERT(!nnStruct->useUfunctions || nnStruct->pointULSHVectors != NULL);
  ASSERT(nnStruct->hashedBuckets[0]->typeHT == HT_LINKED_LIST || nnStruct->hashedBuckets[0]->typeHT == HT_STATISTICS);

  nnStruct->points[nnStruct->nPoints] = point;
  nnStruct->nPoints++;

  preparePointAdding(nnStruct, nnStruct->hashedBuckets[0], point);

  // Initialize the counters for defining the pair of <u> functions used for <g> functions.
  IntT firstUComp = 0;
  IntT secondUComp = 1;

  TIMEV_START(timeBucketIntoUH);
  for(IntT i = 0; i < nnStruct->parameterL; i++){
    if (!nnStruct->useUfunctions) {
      // Use usual <g> functions (truly independent; <g>s are precisly
      // <u>s).
      addBucketEntry(nnStruct->hashedBuckets[i], 1, nnStruct->precomputedHashesOfULSHs[i], NULL, nnStruct->nPoints - 1);
    } else {
      // Use <u> functions (<g>s are pairs of <u> functions).
      addBucketEntry(nnStruct->hashedBuckets[i], 2, nnStruct->precomputedHashesOfULSHs[firstUComp], nnStruct->precomputedHashesOfULSHs[secondUComp], nnStruct->nPoints - 1);

      // compute what is the next pair of <u> functions.
      secondUComp++;
      if (secondUComp == nnStruct->nHFTuples) {
	firstUComp++;
	secondUComp = firstUComp + 1;
      }
    }
    //batchAddRequest(nnStruct, i, firstUComp, secondUComp, point);
  }
  TIMEV_END(timeBucketIntoUH);

  // Check whether the vectors <nearPoints> & <nearPointsIndeces> is still big enough.
  if (nnStruct->nPoints > nnStruct->sizeMarkedPoints) {  //将标记点的数量的大小设为比点的数量大。
    nnStruct->sizeMarkedPoints = 2 * nnStruct->nPoints;
    FAILIF(NULL == (nnStruct->markedPoints = (BooleanT*)REALLOC(nnStruct->markedPoints, nnStruct->sizeMarkedPoints * sizeof(BooleanT))));
    for(IntT i = 0; i < nnStruct->sizeMarkedPoints; i++){
      nnStruct->markedPoints[i] = FALSE;
    }
    FAILIF(NULL == (nnStruct->markedPointsIndeces = (Int32T*)REALLOC(nnStruct->markedPointsIndeces, nnStruct->sizeMarkedPoints * sizeof(Int32T))));
  }
}

// Returns TRUE iff |p1-p2|_2^2 <= threshold
inline BooleanT isDistanceSqrLeq(IntT dimension, PPointT p1, PPointT p2, RealT threshold){
  RealT result = 0;
  nOfDistComps++;

  TIMEV_START(timeDistanceComputation);
  for (IntT i = 0; i < dimension; i++){
    RealT temp = p1->coordinates[i] - p2->coordinates[i];   //计算每一维度的差值
#ifdef USE_L1_DISTANCE
    result += ABS(temp);
#else
    result += SQR(temp);
#endif
    if (result > threshold){
      // TIMEV_END(timeDistanceComputation);
      return 0;
    }
  }
  TIMEV_END(timeDistanceComputation);

  //return result <= threshold;
  return 1;
}

// // Returns TRUE iff |p1-p2|_2^2 <= threshold
// inline BooleanT isDistanceSqrLeq(IntT dimension, PPointT p1, PPointT p2, RealT threshold){
//   RealT result = 0;
//   nOfDistComps++;

//   //TIMEV_START(timeDistanceComputation);
//   for (IntT i = 0; i < dimension; i++){
//     result += p1->coordinates[i] * p2->coordinates[i];
//   }
//   //TIMEV_END(timeDistanceComputation);

//   return p1->sqrLength + p2->sqrLength - 2 * result <= threshold;
// }

// Returns the list of near neighbors of the point <point> (with a
// certain success probability). Near neighbor is defined as being a
// point within distance <parameterR>. Each near neighbor from the
// data set is returned is returned with a certain probability,
// dependent on <parameterK>, <parameterL>, and <parameterT>. The
// returned points are kept in the array <result>. If result is not
// allocated, it will be allocated to at least some minimum size
// (RESULT_INIT_SIZE). If number of returned points is bigger than the
// size of <result>, then the <result> is resized (to up to twice the
// number of returned points). The return value is the number of
// points found.
Int32T getNearNeighborsFromPRNearNeighborStruct(PRNearNeighborStructT nnStruct, PPointT query, PPointT *(&result), Int32T &resultSize){
  ASSERT(nnStruct != NULL);
  ASSERT(query != NULL);
  ASSERT(nnStruct->reducedPoint != NULL);
  ASSERT(!nnStruct->useUfunctions || nnStruct->pointULSHVectors != NULL);  //本代码中U函数为Ture，所以pointULSHVector不会为空

  PPointT point = query;

  if (result == NULL){
    resultSize = RESULT_INIT_SIZE;
    FAILIF(NULL == (result = (PPointT*)MALLOC(resultSize * sizeof(PPointT))));
  }
  
  preparePointAdding(nnStruct, nnStruct->hashedBuckets[0], point);  //将查询点加入到已经建立好的数据结构中

  Uns32T precomputedHashesOfULSHs[nnStruct->nHFTuples][N_PRECOMPUTED_HASHES_NEEDED];
  for(IntT i = 0; i < nnStruct->nHFTuples; i++){
    for(IntT j = 0; j < N_PRECOMPUTED_HASHES_NEEDED; j++){   //点的数的两倍
      precomputedHashesOfULSHs[i][j] = nnStruct->precomputedHashesOfULSHs[i][j];
    }
  }
  TIMEV_START(timeTotalBuckets);   //生成桶的时间？

  BooleanT oldTimingOn = timingOn;
  if (noExpensiveTiming) {
    timingOn = FALSE;
  }
  
  // Initialize the counters for defining the pair of <u> functions used for <g> functions.
  IntT firstUComp = 0;
  IntT secondUComp = 1;

  Int32T nNeighbors = 0;// the number of near neighbors found so far.
  Int32T nMarkedPoints = 0;// the number of marked points
  for(IntT i = 0; i < nnStruct->parameterL; i++){ 
    TIMEV_START(timeGetBucket);   //取得桶的时间？？
    GeneralizedPGBucket gbucket;
    if (!nnStruct->useUfunctions) {   //使用独立的g函数
      // Use usual <g> functions (truly independent; <g>s are precisly
      // <u>s).
      gbucket = getGBucket(nnStruct->hashedBuckets[i], 1, precomputedHashesOfULSHs[i], NULL);  //取得一个桶
    } else {
      // Use <u> functions (<g>s are pairs of <u> functions).
      gbucket = getGBucket(nnStruct->hashedBuckets[i], 2, precomputedHashesOfULSHs[firstUComp], precomputedHashesOfULSHs[secondUComp]);

      // compute what is the next pair of <u> functions.
      secondUComp++;
      if (secondUComp == nnStruct->nHFTuples) {
	firstUComp++;
	secondUComp = firstUComp + 1;
      }
    }
    TIMEV_END(timeGetBucket);  //取得一个hash表中桶的时间

    PGBucketT bucket;

    TIMEV_START(timeCycleBucket);
    switch (nnStruct->hashedBuckets[i]->typeHT){
    case HT_LINKED_LIST:
      bucket = gbucket.llGBucket;
      if (bucket != NULL){
	// circle through the bucket and add to <result> the points that are near.
	PBucketEntryT bucketEntry = &(bucket->firstEntry);
	//TIMEV_START(timeCycleProc);
	while (bucketEntry != NULL){
	  //TIMEV_END(timeCycleProc);
	  //ASSERT(bucketEntry->point != NULL);
	  //TIMEV_START(timeDistanceComputation);
	  Int32T candidatePIndex = bucketEntry->pointIndex;
	  PPointT candidatePoint = nnStruct->points[candidatePIndex];
	  if (isDistanceSqrLeq(nnStruct->dimension, point, candidatePoint, nnStruct->parameterR2) && nnStruct->reportingResult){
	    //TIMEV_END(timeDistanceComputation);
	    if (nnStruct->markedPoints[candidatePIndex] == FALSE) {
	      //TIMEV_START(timeResultStoring);
	      // a new R-NN point was found (not yet in <result>).
	      if (nNeighbors >= resultSize){    //扩大结果集的空间
		// run out of space => resize the <result> array.
		resultSize = 2 * resultSize;
		result = (PPointT*)REALLOC(result, resultSize * sizeof(PPointT));
	      }
	      result[nNeighbors] = candidatePoint;   //将选择的点加入结果集合
	      nNeighbors++;
	      nnStruct->markedPointsIndeces[nMarkedPoints] = candidatePIndex;
	      nnStruct->markedPoints[candidatePIndex] = TRUE; // do not include more points with the same index
	      nMarkedPoints++;
	      //TIMEV_END(timeResultStoring);
	    }
	  }else{
	    //TIMEV_END(timeDistanceComputation);
	  }
	  //TIMEV_START(timeCycleProc);
	  bucketEntry = bucketEntry->nextEntry;  //获取下一个桶
	}
	//TIMEV_END(timeCycleProc);
      }
      break;
    case HT_STATISTICS:
      ASSERT(FALSE); // HT_STATISTICS not supported anymore

//       if (gbucket.linkGBucket != NULL && gbucket.linkGBucket->indexStart != INDEX_START_EMPTY){
// 	Int32T position;
// 	PointsListEntryT *pointsList = nnStruct->hashedBuckets[i]->bucketPoints.pointsList;
// 	position = gbucket.linkGBucket->indexStart;
// 	// circle through the bucket and add to <result> the points that are near.
// 	while (position != INDEX_START_EMPTY){
// 	  PPointT candidatePoint = pointsList[position].point;
// 	  if (isDistanceSqrLeq(nnStruct->dimension, point, candidatePoint, nnStruct->parameterR2) && nnStruct->reportingResult){
// 	    if (nnStruct->nearPoints[candidatePoint->index] == FALSE) {
// 	      // a new R-NN point was found (not yet in <result>).
// 	      if (nNeighbors >= resultSize){
// 		// run out of space => resize the <result> array.
// 		resultSize = 2 * resultSize;
// 		result = (PPointT*)REALLOC(result, resultSize * sizeof(PPointT));
// 	      }
// 	      result[nNeighbors] = candidatePoint;
// 	      nNeighbors++;
// 	      nnStruct->nearPoints[candidatePoint->index] = TRUE; // do not include more points with the same index
// 	    }
// 	  }
// 	  // Int32T oldP = position;
// 	  position = pointsList[position].nextPoint;
// 	  // ASSERT(position == INDEX_START_EMPTY || position == oldP + 1);
// 	}
//       }
      break;
    case HT_HYBRID_CHAINS:
      if (gbucket.hybridGBucket != NULL){
	PHybridChainEntryT hybridPoint = gbucket.hybridGBucket;
	Uns32T offset = 0;
	if (hybridPoint->point.bucketLength == 0){
	  // there are overflow points in this bucket.
	  offset = 0;
	  for(IntT j = 0; j < N_FIELDS_PER_INDEX_OF_OVERFLOW; j++){
	    offset += ((Uns32T)((hybridPoint + 1 + j)->point.bucketLength) << (j * N_BITS_FOR_BUCKET_LENGTH));  //具体移动的位数还要弄懂？
	  }
	}
	Uns32T index = 0;
	BooleanT done = FALSE;
	while(!done){
	  if (index == MAX_NONOVERFLOW_POINTS_PER_BUCKET){
	    //CR_ASSERT(hybridPoint->point.bucketLength == 0);
	    index = index + offset;
	  }
	  Int32T candidatePIndex = (hybridPoint + index)->point.pointIndex;
	  CR_ASSERT(candidatePIndex >= 0 && candidatePIndex < nnStruct->nPoints);
	  done = (hybridPoint + index)->point.isLastPoint == 1 ? TRUE : FALSE;
	  index++;
	  if (nnStruct->markedPoints[candidatePIndex] == FALSE){
	    // mark the point first.
	    nnStruct->markedPointsIndeces[nMarkedPoints] = candidatePIndex;
	    nnStruct->markedPoints[candidatePIndex] = TRUE; // do not include more points with the same index
	    nMarkedPoints++;

	    PPointT candidatePoint = nnStruct->points[candidatePIndex];
	    if (isDistanceSqrLeq(nnStruct->dimension, point, candidatePoint, nnStruct->parameterR2) && nnStruct->reportingResult){
	      //if (nnStruct->markedPoints[candidatePIndex] == FALSE) {
	      // a new R-NN point was found (not yet in <result>).
	      //TIMEV_START(timeResultStoring);
	      if (nNeighbors >= resultSize){
		// run out of space => resize the <result> array.
		resultSize = 2 * resultSize;
		result = (PPointT*)REALLOC(result, resultSize * sizeof(PPointT));
	      }
	      result[nNeighbors] = candidatePoint;
	      nNeighbors++;
	      //TIMEV_END(timeResultStoring);
	      //nnStruct->markedPointsIndeces[nMarkedPoints] = candidatePIndex;
	      //nnStruct->markedPoints[candidatePIndex] = TRUE; // do not include more points with the same index
	      //nMarkedPoints++;
	      //}
	    }
	  }else{
	    // the point was already marked (& examined)
	  }
	}
      }
      break;
    default:
      ASSERT(FALSE);
    }
    TIMEV_END(timeCycleBucket);
    
  }

  timingOn = oldTimingOn;
  TIMEV_END(timeTotalBuckets);

  // we need to clear the array nnStruct->nearPoints for the next query.
  for(Int32T i = 0; i < nMarkedPoints; i++){
    ASSERT(nnStruct->markedPoints[nnStruct->markedPointsIndeces[i]] == TRUE);
    nnStruct->markedPoints[nnStruct->markedPointsIndeces[i]] = FALSE;
  }
  DPRINTF("nMarkedPoints: %d\n", nMarkedPoints);

  return nNeighbors;
}
