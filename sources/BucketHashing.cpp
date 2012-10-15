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

#include "headers.h"


// Creates a new bucket with specified fields. The new bucket contains
// only a single entry -- bucketEntry. bucketEntry->nextEntry is
// expected to be NULL.
inline PGBucketT newGBucket(PUHashStructureT uhash, Uns32T control1, /*PPointT point, */ Int32T pointIndex, PGBucketT nextGBucket){
  PGBucketT bucket;
  if (uhash != NULL && uhash->unusedPGBuckets != NULL){
    bucket = uhash->unusedPGBuckets;
    uhash->unusedPGBuckets = uhash->unusedPGBuckets->nextGBucketInChain;
  } else {
    FAILIF(NULL == (bucket = (PGBucketT)MALLOC(sizeof(GBucketT))));
    nAllocatedGBuckets++;
  }
  ASSERT(bucket != NULL);
  bucket->controlValue1 = control1;
  bucket->firstEntry.pointIndex = pointIndex;
  bucket->firstEntry.nextEntry = NULL;
  bucket->nextGBucketInChain = nextGBucket;

  nGBuckets++;
  return bucket;
}

// Adds the entry <bucketEntry> to the bucket <bucket>.
inline void addPointToGBucket(PUHashStructureT uhash, PGBucketT bucket/*, PPointT point*/ , Int32T pointIndex){
  ASSERT(bucket != NULL);
  ASSERT(uhash != NULL);

  // create a new bucket entry for the point
  TIMEV_START(timeBucketCreation);
  PBucketEntryT bucketEntry;
  if (uhash->unusedPBucketEntrys != NULL){
    bucketEntry = uhash->unusedPBucketEntrys;
    uhash->unusedPBucketEntrys = uhash->unusedPBucketEntrys->nextEntry;
  }else{
    FAILIF(NULL == (bucketEntry = (PBucketEntryT)MALLOC(sizeof(BucketEntryT))));
    nAllocatedBEntries++;
  }
  ASSERT(bucketEntry != NULL);
  bucketEntry->pointIndex = pointIndex;
  TIMEV_END(timeBucketCreation);
  
  bucketEntry->nextEntry = bucket->firstEntry.nextEntry;
  bucket->firstEntry.nextEntry = bucketEntry;
}

// Creates a new UH structure (initializes the hash table and the hash
// functions used). If <typeHT>==HT_PACKED or HT_HYBRID_CHAINS, then
// <modelHT> gives the sizes of all the static arrays that are
// used. Otherwise parameter <modelHT> is not used.
//newUHashStructure(HT_LINKED_LIST, nPoints, nnStruct->parameterK, FALSE, mainHashA, controlHash1, NULL)
PUHashStructureT newUHashStructure(IntT typeHT, Int32T hashTableSize, IntT bucketVectorLength, BooleanT useExternalUHFs, Uns32T *(&mainHashA), Uns32T *(&controlHash1), PUHashStructureT modelHT){
  PUHashStructureT uhash;
  FAILIF(NULL == (uhash = (PUHashStructureT)MALLOC(sizeof(UHashStructureT))));
  uhash->typeHT = typeHT;
  uhash->hashTableSize = hashTableSize;  //哈希表的大小就是点的个数。
  uhash->nHashedBuckets = 0;
  uhash->nHashedPoints = 0;
  uhash->unusedPGBuckets = NULL;
  uhash->unusedPBucketEntrys = NULL;

  uhash->prime = UH_PRIME_DEFAULT;
  uhash->hashedDataLength = bucketVectorLength;   //这边指k。

  uhash->chainSizes = NULL;
  uhash->bucketPoints.pointsArray = NULL;
  uhash->hybridChainsStorage = NULL;

  Int32T totalN = 0;
  Int32T indexInStorage = 0;
  Int32T lastIndexInSt = 0;
  switch (typeHT) {
  case HT_LINKED_LIST:
    FAILIF(NULL == (uhash->hashTable.llHashTable = (PGBucketT*)MALLOC(hashTableSize * sizeof(PGBucketT))));
    uhash->chainSizes = NULL;
    for(Int32T i = 0; i < hashTableSize; i++){  //每个点都会给一个gi的值
      uhash->hashTable.llHashTable[i] = NULL;
    }
    break;
  case HT_STATISTICS:
    ASSERT(FALSE); // Not supported
    FAILIF(NULL == (uhash->hashTable.linkHashTable = (LinkPackedGBucketT**)MALLOC(hashTableSize * sizeof(LinkPackedGBucketT*))));
    FAILIF(NULL == (uhash->chainSizes = (IntT*)(MALLOC(hashTableSize * sizeof(IntT)))));
    FAILIF(NULL == (uhash->bucketPoints.pointsList = (PointsListEntryT*)MALLOC(hashTableSize * sizeof(PointsListEntryT))));
    for(Int32T i = 0; i < hashTableSize; i++){
      uhash->chainSizes[i] = CHAIN_INIT_SIZE;
      FAILIF(NULL == (uhash->hashTable.linkHashTable[i] = (LinkPackedGBucketT*)MALLOC(uhash->chainSizes[i] * sizeof(LinkPackedGBucketT))) && (uhash->chainSizes[i] > 0));
      if (uhash->chainSizes[i] > 0) {
	// the first bucket is empty
	uhash->hashTable.linkHashTable[i][0].indexStart = INDEX_START_EMPTY; 
      }
    }
    break;
  case HT_PACKED:
    ASSERT(FALSE); // Not supported
//     ASSERT(modelHT != NULL);
//     ASSERT(modelHT->typeHT == HT_STATISTICS);
//     ASSERT(modelHT->nHashedPoints == hashTableSize); // TODO
//     FAILIF(NULL == (uhash->hashTable.packedHashTable = (PackedGBucketT**)MALLOC(hashTableSize * sizeof(PackedGBucketT*))));
//     FAILIF(NULL == (uhash->chainSizes = (IntT*)(MALLOC(hashTableSize * sizeof(IntT)))));
//     FAILIF(NULL == (uhash->bucketPoints.pointsArray = (PPointT*)MALLOC(hashTableSize * sizeof(PPointT))));
//     totalN = 0; // total number of points hashed so far.
//     for(Int32T i = 0; i < hashTableSize; i++){

//       // TODO: NOT TESTED AT ALL

//       // count first size of the chain:
//       IntT j;
//       for(j = 0; j < modelHT->chainSizes[i] && modelHT->hashTable.linkHashTable[i][j].indexStart != INDEX_START_EMPTY; j++)
// 	;
//       uhash->chainSizes[i] = j;

//       if (j == 0){
// 	uhash->hashTable.packedHashTable[i] = NULL;
// 	continue;
//       }

//       // allocate memory for the chain
//       FAILIF(NULL == (uhash->hashTable.packedHashTable[i] = (PackedGBucketT*)MALLOC(uhash->chainSizes[i] * sizeof(PackedGBucketT))));

//       // copy each bucket in the chain:
//       for(j = 0; j < uhash->chainSizes[i]; j++){
// 	// "general" info for the bucket
// 	uhash->hashTable.packedHashTable[i][j].controlValue1 = modelHT->hashTable.linkHashTable[i][j].controlValue1;
// 	uhash->hashTable.packedHashTable[i][j].indexStart = totalN;
	
// 	Int32T p = modelHT->hashTable.linkHashTable[i][j].indexStart;
// 	ASSERT(p != INDEX_START_EMPTY);
// 	IntT count = 0;
// 	while(p != INDEX_START_EMPTY){
// 	  uhash->bucketPoints.pointsArray[totalN] = modelHT->bucketPoints.pointsList[p].point;
// 	  totalN++;
// 	  count++;
// 	  p = modelHT->bucketPoints.pointsList[p].nextPoint;
// 	}
// 	uhash->hashTable.packedHashTable[i][j].nPointsInBucket = count;
//       }
//     }
    break;
  case HT_HYBRID_CHAINS:
    ASSERT(modelHT != NULL);
    ASSERT(modelHT->typeHT == HT_LINKED_LIST);
	//hash表项中只存储指针，不存储具体的点，存放真正的h1和h2的值和指针
    FAILIF(NULL == (uhash->hashTable.hybridHashTable = (PHybridChainEntryT*)MALLOC(hashTableSize * sizeof(PHybridChainEntryT))));
	//有多少个表项和点，直接拉链，存放指针。
	//在modeHT开始时 还不知道桶的个数即nHashedBuckets设为0.
	//之后，真正使用混合链表时，已经知道需要多少个桶
    FAILIF(NULL == (uhash->hybridChainsStorage = (HybridChainEntryT*)MALLOC((modelHT->nHashedPoints + modelHT->nHashedBuckets) * sizeof(HybridChainEntryT))));
    
    // the index of the first unoccupied entry in <uhash->hybridChainsStorage>.
    indexInStorage = 0; 
    
    // the index of the last unoccupied entry in <uhash->hybridChainsStorage>
    // (the entries of <uhash->hybridChainsStorage> are filled from start and from end:
    // at the beginning we fill the normal buckets with their points;
    // at the end we fill the "overflow" points of buckets (additional points of buckets that have
    // more than MAX_NONOVERFLOW_POINTS_PER_BUCKET points).
    lastIndexInSt = modelHT->nHashedPoints + modelHT->nHashedBuckets - 1;   //最后一个桶要用来存储超过数量的点
	               
    for(Int32T i = 0; i < hashTableSize; i++){
      PGBucketT bucket = modelHT->hashTable.llHashTable[i];  //使用链形成的每个桶进行处理
      if (bucket != NULL){
	uhash->hashTable.hybridHashTable[i] = uhash->hybridChainsStorage + indexInStorage; // the position where the bucket starts
      }else{
	uhash->hashTable.hybridHashTable[i] = NULL;
      }
      while(bucket != NULL){
	// Compute number of points in the current bucket.
	Int32T nPointsInBucket = 1;
	PBucketEntryT bucketEntry = bucket->firstEntry.nextEntry;
	while(bucketEntry != NULL){
	  nPointsInBucket++;  //计算桶中点的数量
	  bucketEntry = bucketEntry->nextEntry;
	}

	//混合链表中的一个结点放h1的值，接下来一个结点存放是否是最后一个桶，桶长度，是否是最后一个点，点的索引这些值。
	// Copy the points from the bucket to the new HT.
	ASSERT(nPointsInBucket > 0);
	
	uhash->hybridChainsStorage[indexInStorage].controlValue1 = bucket->controlValue1;
	indexInStorage++;
	uhash->hybridChainsStorage[indexInStorage].point.isLastBucket = (bucket->nextGBucketInChain == NULL ? 1 : 0);
	uhash->hybridChainsStorage[indexInStorage].point.bucketLength = (nPointsInBucket <= MAX_NONOVERFLOW_POINTS_PER_BUCKET ?
									 nPointsInBucket : 
									 0); // 0 means there are "overflow" points
	uhash->hybridChainsStorage[indexInStorage].point.isLastPoint = (nPointsInBucket == 1 ? 1 : 0);
	uhash->hybridChainsStorage[indexInStorage].point.pointIndex = bucket->firstEntry.pointIndex;
	indexInStorage++;

	// Store all other points in the storage
	Uns32T currentIndex = indexInStorage; // index where the "current" point will be stored.
	Uns32T nOverflow = 0;
	Uns32T overflowStart = lastIndexInSt;
	if (nPointsInBucket <= MAX_NONOVERFLOW_POINTS_PER_BUCKET){
	  indexInStorage = indexInStorage + nPointsInBucket - 1;
	}else{
	  // bucket too large.
	  // store the overflow points at the end of the array <uhash->hybridChainsStorage>.
	  nOverflow = nPointsInBucket - MAX_NONOVERFLOW_POINTS_PER_BUCKET;  //溢出的点的数量
	  overflowStart = lastIndexInSt - nOverflow + 1;   //溢出点的起始位置
	  lastIndexInSt = overflowStart - 1;    //重置lastIndexInSt。本来lastIndexInSt设为总内存空间中倒数第二个点的位置。

	  // specify the offset of the start of overflow points in the
	  // fields <bucketLength> of points 2, 3, ... of the space
	  // immediately after the bucket. //指定在桶里面溢出点的偏移量 //？？还有有点不清楚
	  Uns32T value = overflowStart - (currentIndex - 1 + MAX_NONOVERFLOW_POINTS_PER_BUCKET);
	  //j<4，用4个节点来表示溢出点存储的偏移量，在接下来的第2、3、4...个点中
	  for(IntT j = 0; j < N_FIELDS_PER_INDEX_OF_OVERFLOW; j++){
	    uhash->hybridChainsStorage[currentIndex + j].point.bucketLength = value & ((1U << N_BITS_FOR_BUCKET_LENGTH) - 1);
	    value = value >> N_BITS_FOR_BUCKET_LENGTH;
	  }

	  // update <indexInStorage>
	  indexInStorage = indexInStorage + MAX_NONOVERFLOW_POINTS_PER_BUCKET - 1;
	  ASSERT(indexInStorage <= lastIndexInSt + 1);

	  //FAILIFWR(nPointsInBucket > MAX_NONOVERFLOW_POINTS_PER_BUCKET, "Too many points in a bucket -- feature not implemented yet. Try to lower N_BITS_PER_POINT_INDEX as much as possible.");// TODO: not implemented yet
	}

	bucketEntry = bucket->firstEntry.nextEntry;
	while(bucketEntry != NULL){
	  uhash->hybridChainsStorage[currentIndex].point.pointIndex = bucketEntry->pointIndex;
	  uhash->hybridChainsStorage[currentIndex].point.isLastPoint = 0;
	  bucketEntry = bucketEntry->nextEntry;

	  currentIndex++;
	  if (currentIndex == indexInStorage && nPointsInBucket > MAX_NONOVERFLOW_POINTS_PER_BUCKET){
	    // finished the normal alloted space -> going to the space reserved at end of the table.
	    currentIndex = overflowStart;
	  }
	}

	// set the <isLastBucket> field of the last point = 1.
	uhash->hybridChainsStorage[currentIndex - 1].point.isLastPoint = 1;
	
	bucket = bucket->nextGBucketInChain;
	//ASSERT((uhash->hashTable.hybridHashTable[i] + 1)->point.bucketLength > 0);
      }
    }
    ASSERT(indexInStorage == lastIndexInSt + 1);
    uhash->nHashedPoints = modelHT->nHashedPoints;
    uhash->nHashedBuckets = modelHT->nHashedBuckets;
    break;
  default:
    ASSERT(FALSE);
  }

  // Initializing the main hash function.
  if (!useExternalUHFs){
    FAILIF(NULL == (uhash->mainHashA = (Uns32T*)MALLOC(uhash->hashedDataLength * sizeof(Uns32T))));
    for(IntT i = 0; i < uhash->hashedDataLength; i++){    //hashedDataLength为k
      uhash->mainHashA[i] = genRandomUns32(1, MAX_HASH_RND);  //
    }
    mainHashA = uhash->mainHashA;
  } else {
    uhash->mainHashA = mainHashA;//是一个指针，指向一个有k个Unsigned int型整数的空间
  }

  // Initializing the control hash functions.
  if (!useExternalUHFs){
    FAILIF(NULL == (uhash->controlHash1 = (Uns32T*)MALLOC(uhash->hashedDataLength * sizeof(Uns32T))));
    for(IntT i = 0; i < uhash->hashedDataLength; i++){
      uhash->controlHash1[i] = genRandomUns32(1, MAX_HASH_RND);   //取一个随机数
    }
    controlHash1 = uhash->controlHash1;
  } else {
    uhash->controlHash1 = controlHash1;
  }

  return uhash;
}

// Removes all the buckets/points from the hash table. Used only for
// HT_LINKED_LIST.
void clearUHashStructure(PUHashStructureT uhash){
  ASSERT(uhash != NULL);
  switch (uhash->typeHT) {
  case HT_LINKED_LIST:
    for(Int32T i = 0; i < uhash->hashTableSize; i++){
      PGBucketT bucket = uhash->hashTable.llHashTable[i];
      while(bucket != NULL){
	PGBucketT tempBucket = bucket;
	bucket = bucket->nextGBucketInChain;
	tempBucket->nextGBucketInChain = uhash->unusedPGBuckets;
	uhash->unusedPGBuckets = tempBucket;
	
	PBucketEntryT bucketEntry = tempBucket->firstEntry.nextEntry;
	while(bucketEntry != NULL){
	  PBucketEntryT tempEntry = bucketEntry;
	  bucketEntry = bucketEntry->nextEntry;
	  tempEntry->nextEntry = uhash->unusedPBucketEntrys;
	  uhash->unusedPBucketEntrys = tempEntry;
	}
      }
      uhash->hashTable.llHashTable[i] = NULL;
    }
    break;
  case HT_STATISTICS:
    ASSERT(FALSE); // Not supported.
//     for(IntT i = 0; i < uhash->hashTableSize; i++){
//       uhash->hashTable.linkHashTable[i][0].indexStart = INDEX_START_EMPTY;
//     }
    break;
  default:
    ASSERT(FALSE);
  }
  uhash->nHashedPoints = 0;
  uhash->nHashedBuckets = 0;
}

// Reorders uhash (of type HT_STATISTICS) to optimize for cache
// behavior. If <auxPtsList> is NULL, it is allocated, but not freed at
// the end (so, the caller has to free <auxPtsList>).
void optimizeUHashStructure(PUHashStructureT uhash, PointsListEntryT *(&auxPtsList)){
  ASSERT(FALSE); // HT_STATISTICS not supported
  ASSERT(uhash->typeHT == HT_STATISTICS);
  
  if (auxPtsList == NULL){
    FAILIF(NULL == (auxPtsList = (PointsListEntryT*)MALLOC(uhash->hashTableSize * sizeof(PointsListEntryT))));
  }
  
  Int32T newP = 0;
  for(Int32T i = 0; i < uhash->hashTableSize; i++){
    for(IntT j = 0; j < uhash->chainSizes[i] && uhash->hashTable.linkHashTable[i][j].indexStart != INDEX_START_EMPTY; j++){
      Int32T p = uhash->hashTable.linkHashTable[i][j].indexStart;
      uhash->hashTable.linkHashTable[i][j].indexStart = newP;
      while (p != INDEX_START_EMPTY) {
	auxPtsList[newP].point = uhash->bucketPoints.pointsList[p].point;
	auxPtsList[newP].nextPoint = newP + 1;
	newP++;
	p = uhash->bucketPoints.pointsList[p].nextPoint;
      }
      auxPtsList[newP - 1].nextPoint = INDEX_START_EMPTY;
    }
  }
  PointsListEntryT *tempList = auxPtsList;
  auxPtsList = uhash->bucketPoints.pointsList;
  uhash->bucketPoints.pointsList = tempList;
}

// Frees the <uhash> structure. If <freeHashFunctions>==FALSE, then
// the hash functions are not freed (because they might be reused, and
// therefore shared by several PUHashStructureT structures).
void freeUHashStructure(PUHashStructureT uhash, BooleanT freeHashFunctions){
  if (uhash == NULL){
    return;
  }

  switch (uhash->typeHT) {
  case HT_LINKED_LIST:
    for(IntT i = 0; i < uhash->hashTableSize; i++){
      PGBucketT bucket = uhash->hashTable.llHashTable[i];
      while (bucket != NULL){
	PGBucketT tempBucket = bucket;
	bucket = bucket->nextGBucketInChain;
	PBucketEntryT bucketEntry = tempBucket->firstEntry.nextEntry;
	while (bucketEntry != NULL){
	  PBucketEntryT tempEntry = bucketEntry;
	  bucketEntry = bucketEntry->nextEntry;
	  free(tempEntry);
	}
	free(tempBucket);
      }
    }
    free(uhash->hashTable.llHashTable);
    if (uhash->unusedPGBuckets != NULL){
      PGBucketT bucket = uhash->unusedPGBuckets;
      while (bucket != NULL){
	PGBucketT tempBucket = bucket;
	bucket = bucket->nextGBucketInChain;
	free(tempBucket);
      }
    }
    if (uhash->unusedPBucketEntrys != NULL){
      PBucketEntryT bucketEntry = uhash->unusedPBucketEntrys;
      while (bucketEntry != NULL){
	PBucketEntryT tempEntry = bucketEntry;
	bucketEntry = bucketEntry->nextEntry;
	free(tempEntry);
      }
    }
    ASSERT(uhash->chainSizes == NULL);
    break;
  case HT_HYBRID_CHAINS:
    free(uhash->hashTable.hybridHashTable);
    free(uhash->hybridChainsStorage);
    ASSERT(uhash->chainSizes == NULL);
    break;
  default:
    ASSERT(FALSE);
  }

  if (freeHashFunctions){
    free(uhash->mainHashA);
    free(uhash->controlHash1);
  }
  free(uhash);
}

// Computes (a.b)mod UH_PRIME_DEFAULT. b is coded as <nBPieces> blocks
// of size totaling <size>. <a> is of length <size>.
//size是一个桶的大小，nBPieces是磁盘块的大小，b是向量坐标点的具体指，a是hash桶大下的一个随机值
//computeBlockProductModDefaultPrime(rndVector, data, nDataPieces, size)
//mainHashA,
inline Uns32T computeBlockProductModDefaultPrime(Uns32T *a, Uns32T *(b[]), IntT nBPieces, IntT size){
  LongUns64T h = 0;
  IntT j = 0;
  IntT bPiece = 0;
  IntT bPieceSize = size / nBPieces;           //一个桶需要多少个磁盘块
  IntT i = 0;
  for(IntT bPiece = 0; bPiece < nBPieces; bPiece++){      //如果一个磁盘块没有装满
    for(IntT j = 0; j < bPieceSize; j++, i++){            //一个桶还没有装满
      h = h + (LongUns64T)a[i] * (LongUns64T)b[bPiece][j];
      h = (h & TWO_TO_32_MINUS_1) + 5 * (h >> 32);
      if (h >= UH_PRIME_DEFAULT) {
	        h = h - UH_PRIME_DEFAULT;                      //h2(a1)=h2(a1)-prime
      }
      CR_ASSERT(h < UH_PRIME_DEFAULT);
    }
  }
  return h;
}

// Computes (a.b)mod UH_PRIME_DEFAULT.  公式h2
//h2(a1)=(〖r^''〗_1 a_1 )mod(2^32-5)=(low[〖r^''〗_1 a_1 ]+5high[〖r^''〗_1 a_1 ] )mod(2^32-5)
//哈希函数，每个哈希函数存储的首地址，k维
//a是一个随机值，b是K维哈希函数的首地址也可以视为一个随机值，size为k   //a是一个分量  //这边只是计算了H2的值
//
inline Uns32T computeProductModDefaultPrime(Uns32T *a, Uns32T *b, IntT size){
  LongUns64T h = 0;
  IntT i = 0;
  for(IntT i = 0; i < size; i++){
    h = h + (LongUns64T)a[i] * (LongUns64T)b[i];  //(∑_(i=1)^k▒〖〖r'〗_i a_i 〗)
    h = (h & TWO_TO_32_MINUS_1) + 5 * (h >> 32);
    if (h >= UH_PRIME_DEFAULT) {
      h = h - UH_PRIME_DEFAULT;
    }
    CR_ASSERT(h < UH_PRIME_DEFAULT);
  }
  return h;
}

// Compute fuction ((rndVector . data)mod prime)mod hashTableSize
// Vectors <rndVector> and <data> are assumed to have length <size>.
//size:一个桶的大小，nDataPieces:数据块的大小
//computeUHashFunction(uhash->mainHashA, tempVector, nBucketVectorPieces, uhash->hashedDataLength, uhash->prime, uhash->hashTableSize);
//h1的哈希值，第一个桶与第二个桶指针的首地址，第几个数据块，长度，素数，哈希表大小
inline Uns32T computeUHashFunction(Uns32T *rndVector, Uns32T *(data[]), IntT nDataPieces, IntT size, Uns32T prime, Int32T hashTableSize){
  ASSERT(prime == UH_PRIME_DEFAULT);
  ASSERT(rndVector != NULL);
  ASSERT(data != NULL);

  Uns32T h = computeBlockProductModDefaultPrime(rndVector, data, nDataPieces, size) % hashTableSize;

  ASSERT(h >= 0 && h < hashTableSize);

  return h;
}
//combinePrecomputedHashes(firstBucketVector, secondBucketVector, nBucketVectorPieces, UHF_MAIN_INDEX)
//nBucketVectorPieces的值是2，UHF_MAIN_INDEX的值是0.
inline Uns32T combinePrecomputedHashes(Uns32T *firstBucketVector, Uns32T *secondBucketVector, IntT nBucketVectorPieces, IntT uhfIndex){
  // CR_ASSERT(bucketVector != NULL);
//   if (nBucketVectorPieces == 1) {
//     // using normal <g> functions.
//     CR_ASSERT(bucketVector[1] != NULL);
//     return (bucketVector[1][uhfIndex] % UH_PRIME_DEFAULT);
//   } else {
//     CR_ASSERT(nBucketVectorPieces == 2); // each of the <g> functions is a pair of 2 <u> functions.
//     //CR_ASSERT(bucketVector[2] != NULL);
//     //CR_ASSERT(bucketVector[3] != NULL);
//     LongUns64T r = (LongUns64T)(bucketVector[2][uhfIndex]) + (LongUns64T)(bucketVector[3][uhfIndex + UHF_NUMBER_OF_HASHES]);
//     if (r >= UH_PRIME_DEFAULT) {
//       r -= UH_PRIME_DEFAULT;
//     }
//     CR_ASSERT(r < UH_PRIME_DEFAULT);
//     return (Uns32T)r;
//   }
  if (nBucketVectorPieces == 1) {  //对正常g函数使用
    // using normal <g> functions.
    Uns32T h = firstBucketVector[uhfIndex];
    if (h > UH_PRIME_DEFAULT){
      h = h - UH_PRIME_DEFAULT;
    }
    return h;
  } else {     //每个g函数有一对的 2<u>函数组成
    CR_ASSERT(nBucketVectorPieces == 2); // each of the <g> functions is a pair of 2 <u> functions.
    //CR_ASSERT(bucketVector[2] != NULL);
    //CR_ASSERT(bucketVector[3] != NULL);
    LongUns64T r = (LongUns64T)(firstBucketVector[uhfIndex]) + (LongUns64T)(secondBucketVector[uhfIndex + UHF_NUMBER_OF_HASHES]);  //2
    if (r >= UH_PRIME_DEFAULT) {
      r -= UH_PRIME_DEFAULT;
    }
    CR_ASSERT(r < UH_PRIME_DEFAULT);
    return (Uns32T)r;
  }
}

// Adds the bucket entry (a point <point>) to the bucket defined by
// bucketVector in the uh structure with number uhsNumber. If no such
// bucket exists, then it is first created.
//LocalitySensitiveHashing.cpp Line 506
//addBucketEntry(nnStruct->hashedBuckets[i], 2, nnStruct->precomputedHashesOfULSHs[firstUComp], nnStruct->precomputedHashesOfULSHs[secondUComp], nnStruct->nPoints - 1);
//哈希桶的结构，桶向量的块数，第一个桶向量，第二个桶向量，点的索引号。
void addBucketEntry(PUHashStructureT uhash, IntT nBucketVectorPieces, Uns32T firstBucketVector[], Uns32T secondBucketVector[]/*, PPointT point*/ , Int32T pointIndex){
  CR_ASSERT(uhash != NULL);
  // CR_ASSERT(bucketVector != NULL);

  Uns32T hIndex;    //为了求h1，h2(a1,a2,..,av)，整个h2的值
	Uns32T control1;   //为了求h2，h(a1)的值

  if (!USE_PRECOMPUTED_HASHES){
    // if not using the same hash functions across multiple
    // UHashStructureT, then we need to compute explicitly the hases.
    Uns32T *tempVector[2];
    tempVector[0] = firstBucketVector;
    tempVector[1] = secondBucketVector;
	//计算h1的值和 h2的值
    hIndex = computeUHashFunction(uhash->mainHashA, tempVector, nBucketVectorPieces, uhash->hashedDataLength, uhash->prime, uhash->hashTableSize);
    control1 = computeBlockProductModDefaultPrime(uhash->controlHash1, tempVector, nBucketVectorPieces, uhash->hashedDataLength);
  } else {
    // if using the same hash functions across multiple
    // UHashStructureT, then we can use the (possibly partially)
    // precomputed hash values.
    CR_ASSERT(uhash->prime == UH_PRIME_DEFAULT);        //计算hj;hj=（left+right）mod A mod b ,对一个点进行计算
    hIndex = combinePrecomputedHashes(firstBucketVector, secondBucketVector, nBucketVectorPieces, UHF_MAIN_INDEX) % uhash->hashTableSize;
    control1 = combinePrecomputedHashes(firstBucketVector, secondBucketVector, nBucketVectorPieces, UHF_CONTROL1_INDEX);
  }

  PGBucketT p;
  BooleanT found;
  Int32T j;
  Int32T temp;
  switch (uhash->typeHT) {    //桶结构的整理
  case HT_LINKED_LIST:    //先找下是否链中已有这个h1的值，如果没有,新建一个bucket，且bucket中的值为这个点的h1的值
	                     //如果有，就将这个点加到与其h1值相等的点中。
    p = uhash->hashTable.llHashTable[hIndex];
    while(p != NULL && 
	  (p->controlValue1 != control1)) {
      p = p->nextGBucketInChain;
    }
    if (p == NULL) {
      // new bucket to add to the hash table
      uhash->nHashedBuckets++;
      uhash->hashTable.llHashTable[hIndex] = newGBucket(uhash,
							control1, 
							pointIndex, 
							uhash->hashTable.llHashTable[hIndex]);
    } else {
      // add this bucket entry to the existing bucket
      addPointToGBucket(uhash, p, pointIndex);
    }
    break;
  case HT_PACKED:
//     // The bucket should already exist.
//     IntT i;
//     for(i = 0; i < uhash->chainSizes[hIndex] && uhash->hashTable.packedHashTable[hIndex][i].nPointsInBucket > 0; i++){
//       if (uhash->hashTable.packedHashTable[hIndex][i].controlValue1 == control1){
// 	break;
//       }
//     }
//     uhash->hashTable.packedHashTable[hIndex][i].nPointsInBucket++;
    break;
  case HT_STATISTICS:
    ASSERT(FALSE);
//     uhash->bucketPoints.pointsList[uhash->nHashedPoints].point = point;
//     found = FALSE;
//     for(j = 0; j < uhash->chainSizes[hIndex] && uhash->hashTable.linkHashTable[hIndex][j].indexStart != INDEX_START_EMPTY; j++){
//       if (uhash->hashTable.linkHashTable[hIndex][j].controlValue1 == control1){
// 	found = TRUE;
// 	break;
//       }
//     }
//     if (!found) {
//       // new bucket
//       if (j >= uhash->chainSizes[hIndex]) {
// 	// dont have enough space in pREALLOCated memory.
// 	if (uhash->chainSizes[hIndex] > 0) {
// 	  uhash->chainSizes[hIndex] = CEIL(uhash->chainSizes[hIndex] * CHAIN_RESIZE_RATIO);
// 	}else{
// 	  uhash->chainSizes[hIndex] = 1;
// 	}
// 	uhash->hashTable.linkHashTable[hIndex] = (LinkPackedGBucketT*)REALLOC(uhash->hashTable.linkHashTable[hIndex], uhash->chainSizes[hIndex] * sizeof(LinkPackedGBucketT));
// 	uhash->hashTable.linkHashTable[hIndex][j].controlValue1 = control1;
//       }
//       uhash->hashTable.linkHashTable[hIndex][j].controlValue1 = control1;
//       uhash->hashTable.linkHashTable[hIndex][j].indexStart = INDEX_START_EMPTY;
//       uhash->nHashedBuckets++;
//       if (j + 1 < uhash->chainSizes[hIndex]){
// 	uhash->hashTable.linkHashTable[hIndex][j + 1].indexStart = INDEX_START_EMPTY;
//       }
//     }
    
//     temp = uhash->hashTable.linkHashTable[hIndex][j].indexStart;
//     uhash->hashTable.linkHashTable[hIndex][j].indexStart = uhash->nHashedPoints;
//     uhash->bucketPoints.pointsList[uhash->nHashedPoints].nextPoint = temp;

    break;
  default:
    ASSERT(FALSE);
  }
  uhash->nHashedPoints++;
}

// Returns the bucket defined by the vector <bucketVector> in the UH
// structure number <uhsNumber>.
GeneralizedPGBucket getGBucket(PUHashStructureT uhash, IntT nBucketVectorPieces, Uns32T firstBucketVector[], Uns32T secondBucketVector[]){
  Uns32T hIndex;
  Uns32T control1;

  //TIMEV_START(timeGBHash);
  if (!USE_PRECOMPUTED_HASHES){
    // if not using the same hash functions across multiple
    // UHashStructureT, then we need to compute explicitly the hases.
	//使用不同的hash函数
    Uns32T *tempVector[2];
    tempVector[0] = firstBucketVector;
    tempVector[1] = secondBucketVector;
    hIndex = computeUHashFunction(uhash->mainHashA, tempVector, nBucketVectorPieces, uhash->hashedDataLength, uhash->prime, uhash->hashTableSize);
    control1 = computeBlockProductModDefaultPrime(uhash->controlHash1, tempVector, nBucketVectorPieces, uhash->hashedDataLength);
  } else {
  //使用相同的hash函数
    // if using the same hash functions across multiple
    // UHashStructureT, then we can use the (possibly partially)
    // precomputed hash values.
    CR_ASSERT(uhash->prime == UH_PRIME_DEFAULT);
    hIndex = combinePrecomputedHashes(firstBucketVector, secondBucketVector, nBucketVectorPieces, UHF_MAIN_INDEX) % uhash->hashTableSize;
    control1 = combinePrecomputedHashes(firstBucketVector, secondBucketVector, nBucketVectorPieces, UHF_CONTROL1_INDEX);
  }
  //TIMEV_END(timeGBHash);

  GeneralizedPGBucket result;
  PGBucketT p;
  PHybridChainEntryT indexHybrid = NULL;
  //TIMEV_START(timeChainTraversal);
  switch(uhash->typeHT) {
  case HT_LINKED_LIST:
    p = uhash->hashTable.llHashTable[hIndex];
    while(p != NULL && 
	  (p->controlValue1 != control1)) {
      p = p->nextGBucketInChain;
      nBucketsInChains++;
    }
    result.llGBucket = p;
    return result;
  case HT_STATISTICS:
    ASSERT(FALSE); // HT_STATISTICS not supported anymore
    for(IntT j = 0; j < uhash->chainSizes[hIndex] && uhash->hashTable.linkHashTable[hIndex][j].indexStart != INDEX_START_EMPTY; j++){
      if (uhash->hashTable.linkHashTable[hIndex][j].controlValue1 == control1){
	result.linkGBucket = &(uhash->hashTable.linkHashTable[hIndex][j]);
	return result;
      }
    }
    result.linkGBucket = NULL;
    return result;
  case HT_PACKED:
    ASSERT(FALSE); // HT_PACKED not supported anymore
    for(IntT j = 0; j < uhash->chainSizes[hIndex]; j++){
      if (uhash->hashTable.packedHashTable[hIndex][j].controlValue1 == control1){
	result.packedGBucket = &(uhash->hashTable.packedHashTable[hIndex][j]);
	return result;
      }
    }
    result.packedGBucket = NULL;
    return result;
  case HT_HYBRID_CHAINS:
    indexHybrid = uhash->hashTable.hybridHashTable[hIndex];
    while (indexHybrid != NULL){ 
      if (indexHybrid->controlValue1 == control1){
	result.hybridGBucket = indexHybrid + 1;
	return result;
      }else{
	indexHybrid = indexHybrid + 1;
	if (indexHybrid->point.isLastBucket != 0){
	  result.hybridGBucket = NULL;
	  return result;
	}
	indexHybrid = indexHybrid + indexHybrid->point.bucketLength;
      }
    }
    result.hybridGBucket = NULL;
    return result;
    break;
  default:
    ASSERT(FALSE);
  }
  //TIMEV_END(timeChainTraversal);

}
//precomputeUHFsForULSH(uhash, nnStruct->pointULSHVectors[i], nnStruct->hfTuplesLength, nnStruct->precomputedHashesOfULSHs[i]);
//LocalitySensitiveHashing.cpp Line 444 ,共有L个g函数
//哈希桶的结构，g函数的头指针地址g(v)=(h1(v),h2(v),...,hk(v))，k的长度，具体g函数的值
void precomputeUHFsForULSH(PUHashStructureT uhash, Uns32T *uVector, IntT length, Uns32T *result){  //result可能是一个真正的hash表
  if (length == uhash->hashedDataLength){  //传进来的k值==哈希表真实的k
    //建立物理索引
    result[UHF_MAIN_INDEX] = computeProductModDefaultPrime(uhash->mainHashA, uVector, length);
    result[UHF_CONTROL1_INDEX] = computeProductModDefaultPrime(uhash->controlHash1, uVector, length);
  } else {     
    ASSERT(2 * length == uhash->hashedDataLength); // the length is 1/2 of the bucket length
	   //k/2 ，??怎么建立两个hash表 
	//UHF_MAIN_INDEX 主索引  //并没有存实际的k维的哈希值，最存了实际的值
	//  一个随机值line265，每个哈希函数存储的首地址，k维
    result[UHF_MAIN_INDEX] = computeProductModDefaultPrime(uhash->mainHashA, uVector, length); //uVector是g函数中的k维的向量
    result[UHF_CONTROL1_INDEX] = computeProductModDefaultPrime(uhash->controlHash1, uVector, length);
    result[UHF_MAIN_INDEX + UHF_NUMBER_OF_HASHES] = computeProductModDefaultPrime(uhash->mainHashA + length, uVector, length);
    result[UHF_CONTROL1_INDEX + UHF_NUMBER_OF_HASHES] = computeProductModDefaultPrime(uhash->controlHash1 + length, uVector, length);
  }
}

